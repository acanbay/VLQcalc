'''
Author: Ali Can Canbay - acanbay
Date  : 21 May 2024

This module contains the VLQ object. It debugs commands and communicates with other modules.
'''
from VLQcalc import decays
from VLQcalc import modelConverter
from VLQcalc import functions
import math, numbers, warnings
############
def isNumber(value):
    return isinstance(value, numbers.Real) and not isinstance(value, bool)

def toNumber(value):
    # numpy scalars are converted to plain int/float so that cards and files look the same
    return value.item() if hasattr(value, 'item') else value

def checkBRs(BRs, where, text, name='branchings'):
    # Three non-negative fractions whose sum is 1 (within rounding)
    if type(BRs) not in [list,tuple] or len(BRs) != 3 or not all(isNumber(BR) for BR in BRs):
        raise ValueError(f'{where}: {text}')
    BRs = [toNumber(BR) for BR in BRs]
    if any(BR < 0 for BR in BRs) or not math.isclose(sum(BRs), 1, abs_tol=1e-9):
        raise ValueError(f'{where}: the sum of all {name} must be 1 and each of them must be non-negative')
    return BRs

def checkZeta(zeta, where):
    # Generation fractions [zeta1, zeta2, zeta3]; the default couples to the 3rd generation only
    if zeta is None:
        return [0, 0, 1]
    return checkBRs(zeta, where, 'zeta must have 3 generation fractions -> [zeta1, zeta2, zeta3]', 'generation fractions (zeta)')

class VLQ:
    '''
    A vector-like quark: VLQ(VLQtype, FNS=4, LR=False)

    VLQtype: 'T', 'B', 'X' or 'Y' (also 'tp', 'bp', 'VLT', ...)
    FNS    : flavour scheme, 4 or 5 (m_b = 0 in 5FNS)
    LR     : False for left-handed couplings only, True for left- and right-handed couplings (kappa_R = kappa_L)
    '''
    def __init__(self,VLQtype,FNS=4,LR=False):
        # Checking the VLQ type
        if type(VLQtype)!=str:
            raise ValueError('VLQ: VLQ type must be string')
        if VLQtype.lower() not in ['t', 'vlt', 'tp', 'b', 'vlb', 'bp', 'x', 'vlx', 'y', 'vly']:
            raise ValueError('VLQ: VLQ type must be one of the following: '
                             'tp (charge: 2/3), bp (charge: -1/3), x (charge: 5/3), y (charge: -4/3)')
        if VLQtype.lower() in ['t', 'vlt', 'tp']:
            self.type = 'VLT'
            self.BRtext = 'BR(Ht),BR(Wb),BR(Zt)'
            self.NewKappatext = 'kH,kW,kZ'
            self.OldKappatext = 'kT'
        elif VLQtype.lower() in ['b', 'vlb', 'bp']:
            self.type = 'VLB'
            self.BRtext = 'BR(Hb),BR(Wt),BR(Zb)'
            self.NewKappatext = 'kH,kW,kZ'
            self.OldKappatext = 'kB'
        elif VLQtype.lower() in ['x', 'vlx']:
            self.type = 'VLX'
            self.BRtext = 'BR(Wt)'
            self.NewKappatext = 'kW'
            self.OldKappatext = 'kX'
        elif VLQtype.lower() in ['y', 'vly']:
            self.type = 'VLY'
            self.BRtext = 'BR(Wb)'
            self.NewKappatext = 'kW'
            self.OldKappatext = 'kY'

        # Checking the FNS value
        if FNS not in [4,5] or isinstance(FNS, bool):
            raise ValueError('VLQ: Flavor value must be 4 or 5')

        if FNS == 4:
            if self.type == 'VLT':
                self.calcDecay = decays.decayTP4
            elif self.type == 'VLB':
                self.calcDecay = decays.decayBP4
            elif self.type == 'VLX':
                self.calcDecay = decays.decayX4
            elif self.type == 'VLY':
                self.calcDecay = decays.decayY4
        elif FNS == 5:
            if self.type == 'VLT':
                self.calcDecay = decays.decayTP5
            elif self.type == 'VLB':
                self.calcDecay = decays.decayBP5
            elif self.type == 'VLX':
                self.calcDecay = decays.decayX5
            elif self.type == 'VLY':
                self.calcDecay = decays.decayY5

        # Checking the LR value
        if type(LR)!=bool:
            raise ValueError('VLQ: LR must be True (both left and right handed) or False (only left handed)')
        self.LR = LR

        self.Mass = None
        self.Branchings = None
        self.KappaOld = None
        self.KappaH = None
        self.KappaW = None
        self.KappaZ = None
        self.FNS = FNS
        self.BRs = None
        self.Ratio = None
        self.zeta = [0, 0, 1]
        # Couplings to the 1st and 2nd generation quarks: {(boson, generation): [kappa per mass]}
        self.KappaLight = {}

    @property
    def bosons(self):
        return list(decays.CHANNELS[self.type])

    @property
    def Kappas(self):
        '''All couplings named as in the VLQ_v4/VLQ_v5 models, e.g. {'KTLh1': [...], ..., 'KTRz3': [...]}'''
        if self.Mass == None:
            return {}
        zeros = [0]*len(self.Mass)
        third = {'h': self.KappaH, 'w': self.KappaW, 'z': self.KappaZ}
        letter = self.type[-1]
        Kappas = {}
        for boson in self.bosons:
            for gen in [1,2,3]:
                if gen == 3:
                    values = third[boson] if third[boson] != None else zeros
                else:
                    values = self.KappaLight.get((boson,gen), zeros)
                name = boson if self.type in ['VLT','VLB'] else ''
                Kappas[f'K{letter}L{name}{gen}'] = list(values)
                Kappas[f'K{letter}R{name}{gen}'] = list(values) if self.LR else list(zeros)
        return Kappas

    def setKappas(self, Kappas):
        '''Stores couplings given as {(boson, generation): [kappa per mass]}'''
        if 'h' in self.bosons:
            self.KappaH = Kappas.get(('h',3), [0]*len(self.Mass))
            self.KappaZ = Kappas.get(('z',3), [0]*len(self.Mass))
        self.KappaW = Kappas.get(('w',3), [0]*len(self.Mass))
        self.KappaLight = {channel: values for channel, values in Kappas.items() if channel[1] != 3}

    def calcWidths(self, index):
        '''Partial widths {(boson, generation): width} in GeV for Mass[index] with the stored couplings'''
        Kappas = {}
        third = {'h': self.KappaH, 'w': self.KappaW, 'z': self.KappaZ}
        for boson in self.bosons:
            for gen in [1,2,3]:
                values = third[boson] if gen == 3 else self.KappaLight.get((boson,gen))
                if values != None:
                    Kappas[(boson,gen)] = values[index] if len(values) > index else values[0]
        return decays.partialWidths(self.type, self.FNS, self.Mass[index], Kappas, self.LR)

    def setMass(self, Masses):
        '''Sets the mass in GeV: a number, or a list, tuple, range or numpy array of numbers'''
        if isNumber(Masses):
            Masses = [Masses]
        try:
            Masses = [toNumber(Mass) for Mass in Masses]
        except TypeError:
            raise ValueError('VLQ.setMass: Mass must be number or number list') from None
        if len(Masses) == 0 or not all(isNumber(Mass) and Mass > 0 for Mass in Masses):
            raise ValueError('VLQ.setMass: Mass must be a positive number or a list of positive numbers')
        self.Mass = Masses

    def convertModel(self,Kappas,BRs=1,reverse=False,zeta=None):
        '''
        Converts couplings between VLQ_UFO and VLQ_v4/VLQ_v5.

        reverse=False: Kappas is kappa_Q of VLQ_UFO, BRs = [xiH, xiW, xiZ] (T, B), zeta = [zeta1, zeta2, zeta3];
                       the results are stored in KappaH, KappaW, KappaZ and KappaLight.
        reverse=True : Kappas = [kH, kW, kZ] (T, B) or kW (X, Y), one value per mass;
                       the result is stored in KappaOld.
        '''
        if self.Mass == None:
            raise ValueError('VLQ.convertModel: It was not specified for which mass the new kappas would be calculated. '
                             'Please define mass to VLQ object with setMass method')
        if type(reverse)!=bool:
            raise ValueError('VLQ.convertModel: reverse must be False (VLQ_UFO -> VLQ_v4) or True (VLQ_v4 -> VLQ_UFO)')
        ##
        if not reverse:
            if self.type in ['VLX','VLY']:
                BRs = 1
            else:
                BRs = checkBRs(BRs, 'VLQ.convertModel',
                               f'{self.type} must have 3 float branching values (xi) in old model -> [xiH, xiW, xiZ]')
            self.zeta = checkZeta(zeta, 'VLQ.convertModel')
            if type(Kappas) in [list,tuple] and len(Kappas)==1:
                Kappas = Kappas[0]
            if not isNumber(Kappas):
                raise ValueError(f'VLQ.convertModel: kappa value to be converted must be a number -> {self.OldKappatext}')
            if self.KappaH!=None or self.KappaW!=None or self.KappaZ!=None:
                warnings.warn(f'VLQ.convertModel: deleting old kappa values for {self.NewKappatext}', stacklevel=2)
            self.KappaOld = [toNumber(Kappas)]
            self.BRs = BRs
            self.setKappas(modelConverter.convert(self,toNumber(Kappas),BRs,reverse))
        else:
            if zeta is not None:
                raise ValueError('VLQ.convertModel: zeta is only used for VLQ_UFO -> VLQ_v4; '
                                 'for VLQ_v4 -> VLQ_UFO the generation couplings come from the VLQ object')
            n = len(self.Mass)
            def asList(values):
                # one value per mass
                values = [values] if isNumber(values) else values
                if type(values) not in [list,tuple] or len(values) != n or not all(isNumber(k) for k in values):
                    raise ValueError(f'VLQ.convertModel: kappa values to be converted must be numbers for {self.type} -> '
                                     f'{self.NewKappatext}, one value per mass ({n})')
                return [toNumber(k) for k in values]
            if self.type in ['VLX','VLY']:
                if type(Kappas) in [list,tuple] and len(Kappas) == 1 and n != 1:
                    Kappas = Kappas[0]
                self.KappaW = asList(Kappas)
            else:
                if type(Kappas) not in [list,tuple] or len(Kappas) != 3:
                    raise ValueError(f'VLQ.convertModel: {self.type} needs 3 kappa values (or lists) -> [{self.NewKappatext}]')
                self.KappaH, self.KappaW, self.KappaZ = [asList(ks) for ks in Kappas]
            if self.KappaOld!=None:
                warnings.warn(f'VLQ.convertModel: deleting old kappa values for {self.OldKappatext}', stacklevel=2)
            self.KappaOld = modelConverter.convert(self,None,reverse=reverse)

    def calcRatioKappas(self,BRs=None,Ratio=None,zeta=None):
        '''
        Computes the couplings for each mass from the width-to-mass ratio and the branching ratios.

        BRs  : [BR(H), BR(W), BR(Z)] for T and B (sum 1); not used for X and Y
        Ratio: Gamma/M
        zeta : generation fractions [zeta1, zeta2, zeta3] (sum 1); default [0, 0, 1]
        '''
        if self.Mass == None:
            raise ValueError('VLQ.calcRatioKappas: Please define mass to VLQ object with setMass method')
        if not isNumber(Ratio) or Ratio <= 0:
            raise ValueError(f'VLQ.calcRatioKappas: Gamma/M({self.type}) ratio must be a positive number (int or float)')
        if self.type in ['VLX','VLY']:
            BRs = 1
        else:
            BRs = checkBRs(BRs, 'VLQ.calcRatioKappas',
                           f'{self.type} must have 3 float branching values for {self.type} -> [{self.BRtext}]')
        self.Ratio = toNumber(Ratio)
        self.BRs = BRs
        self.zeta = checkZeta(zeta, 'VLQ.calcRatioKappas')
        functions.calcRatioKappas(self)

    def table(self, style='md'):
        '''Returns the results as a Markdown (style='md') or aligned text (style='txt') table'''
        from VLQcalc import output
        return output.table(self, style)

    def save(self, filename):
        '''Writes the results to a file; the format follows the extension: .md, .txt, .csv or .json'''
        from VLQcalc import output
        output.save(self, filename)

    def plot(self, quantity='kappa', save=None, ax=None, logy=None):
        '''
        Plots the couplings ('kappa'), branching ratios ('BR') or partial widths ('width') against the mass.

        save: file name of the figure; ax: existing axes (single panel only);
        logy: True or False for the y scale, chosen from the values if None. Needs matplotlib.
        '''
        from VLQcalc import output
        return output.plot(self, quantity, save, ax, logy)
