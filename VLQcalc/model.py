'''
Author : Ali Can Canbay - acanbay
Created: 21 May 2024
Updated: 30 October 2026

This module contains the VLQ object. It debugs commands and communicates with other modules.
'''
from VLQcalc import decays
from VLQcalc import modelConverter
from VLQcalc import functions
from VLQcalc import constants
import math, numbers, warnings
############
def isNumber(value):
    return isinstance(value, numbers.Real) and not isinstance(value, bool)

def toNumber(value):
    # numpy scalars are converted to plain int/float so that cards and files look the same
    return value.item() if hasattr(value, 'item') else value

def isFinite(value):
    return isNumber(value) and math.isfinite(value)

def perMass(values, n, where):
    # A list with one value per mass; a single value is used for every mass
    if isNumber(values):
        values = [values]
    values = list(values)
    if len(values) == 1:
        return values*n
    if len(values) != n:
        raise ValueError(f'{where}: {len(values)} values for {n} masses; give one value per mass')
    return values

def checkBRs(BRs, where, text, name='branchings'):
    # Three non-negative fractions whose sum is 1 (within rounding)
    try:
        BRs = [toNumber(BR) for BR in BRs] if not isinstance(BRs, (str, dict, set)) else None
    except TypeError:
        BRs = None
    if BRs is None or len(BRs) != 3 or not all(isFinite(BR) for BR in BRs):
        raise ValueError(f'{where}: {text}')
    # the tolerance accepts fractions rounded to 6 significant digits, as VLQcalc writes them
    if any(BR < 0 for BR in BRs) or not math.isclose(sum(BRs), 1, abs_tol=1e-5):
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
        if not isNumber(FNS) or FNS not in [4,5]:
            raise ValueError('VLQ: Flavor value must be 4 or 5')
        FNS = int(FNS)

        if FNS == 4:
            if self.type == 'VLT':
                self.decay = decays.decayTP4
            elif self.type == 'VLB':
                self.decay = decays.decayBP4
            elif self.type == 'VLX':
                self.decay = decays.decayX4
            elif self.type == 'VLY':
                self.decay = decays.decayY4
        elif FNS == 5:
            if self.type == 'VLT':
                self.decay = decays.decayTP5
            elif self.type == 'VLB':
                self.decay = decays.decayBP5
            elif self.type == 'VLX':
                self.decay = decays.decayX5
            elif self.type == 'VLY':
                self.decay = decays.decayY5

        # Checking the LR value
        if type(LR)!=bool:
            raise ValueError('VLQ: LR must be True (both left and right handed) or False (only left handed)')
        self.LR = LR

        self.Mass = None
        self.Branchings = None
        self.KappaOld = None
        # KappaOld computed by convertModel(..., reverse=True); it belongs to the masses it was computed for
        self.KappaOldConverted = None
        self.KappaH = None
        self.KappaW = None
        self.KappaZ = None
        self.FNS = FNS
        self.BRs = None
        self.Ratio = None
        self.zeta = [0, 0, 1]
        # Couplings to the 1st and 2nd generation quarks: {(boson, generation): [kappa per mass]}
        self.KappaLight = {}
        # SM inputs the couplings were computed with (constants.setSM may change them later)
        self.SM = None

    def calcDecay(self, M, Kappas, LR=None):
        '''
        Widths of the 3rd generation channels for the mass M and Kappas = [kH, kW, kZ] (T, B) or kW (X, Y):
        (H, W, Z, total) for T and B, the W width for X and Y. LR is the one of the object unless given.
        '''
        return self.decay(M, Kappas, self.LR if LR is None else LR)

    @property
    def bosons(self):
        return list(decays.CHANNELS[self.type])

    @property
    def Kappas(self):
        '''All couplings named as in the VLQ_v4/VLQ_v5 models, e.g. {'KTLh1': [...], ..., 'KTRz3': [...]}'''
        if self.Mass is None:
            return {}
        zeros = [0]*len(self.Mass)
        third = {'h': self.KappaH, 'w': self.KappaW, 'z': self.KappaZ}
        letter = self.type[-1]
        Kappas = {}
        for boson in self.bosons:
            for gen in [1,2,3]:
                if gen == 3:
                    values = perMass(third[boson], len(self.Mass), 'VLQ.Kappas') if third[boson] is not None else zeros
                else:
                    values = perMass(self.KappaLight.get((boson,gen), zeros), len(self.Mass), 'VLQ.Kappas')
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
                if values is not None:
                    Kappas[(boson,gen)] = perMass(values, len(self.Mass), 'VLQ.calcWidths')[index]
        with constants.using(self.SM):
            return decays.partialWidths(self.type, self.FNS, self.Mass[index], Kappas, self.LR)

    def setMass(self, Masses):
        '''
        Sets the mass in GeV: a number, or a list, tuple, range or numpy array of numbers.
        Results computed for other masses are deleted.
        '''
        if getattr(Masses, 'ndim', None) == 0:
            Masses = Masses.item()
        if isNumber(Masses):
            Masses = [Masses]
        if isinstance(Masses, (str, dict, set, frozenset)):
            raise ValueError('VLQ.setMass: Mass must be number or number list')
        try:
            Masses = [toNumber(Mass) for Mass in Masses]
        except (TypeError, ValueError):
            raise ValueError('VLQ.setMass: Mass must be number or number list') from None
        if len(Masses) == 0 or not all(isFinite(Mass) and Mass > 0 for Mass in Masses):
            raise ValueError('VLQ.setMass: Mass must be a positive number or a list of positive numbers')
        if self.Mass is not None and Masses != self.Mass:
            self.clearResults('VLQ.setMass')
        self.Mass = Masses

    def clearResults(self, where):
        # Couplings and branchings computed for the previous masses do not belong to the new ones.
        # A single VLQ_UFO coupling and a single set of branchings (e.g. set by hand) apply to every mass.
        perMassResults = [self.KappaH, self.KappaW, self.KappaZ] + list(self.KappaLight.values())
        if any(values is not None for values in perMassResults):
            warnings.warn(f'{where}: the mass changed; deleting the couplings computed for the previous masses', stacklevel=3)
        self.KappaH = self.KappaW = self.KappaZ = None
        self.KappaLight = {}
        self.Ratio = None
        if self.KappaOld is not None and (self.KappaOld is self.KappaOldConverted
                                          or (not isNumber(self.KappaOld) and len(self.KappaOld) != 1)):
            self.KappaOld = None
        self.KappaOldConverted = None
        self.SM = None
        if type(self.BRs) == list and any(type(BR) == list for BR in self.BRs):
            self.BRs = None

    def convertModel(self,Kappas,BRs=1,reverse=False,zeta=None):
        '''
        Converts couplings between VLQ_UFO and VLQ_v4/VLQ_v5.

        reverse=False: Kappas is kappa_Q of VLQ_UFO, BRs = [xiH, xiW, xiZ] (T, B), zeta = [zeta1, zeta2, zeta3];
                       the results are stored in KappaH, KappaW, KappaZ and KappaLight.
        reverse=True : Kappas = [kH, kW, kZ] (T, B) or kW (X, Y), one value per mass;
                       kappa_Q is stored in KappaOld, the xi of VLQ_UFO in BRs (T, B) and its zeta in zeta.
        '''
        if self.Mass is None:
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
            zeta = checkZeta(zeta, 'VLQ.convertModel')
            if hasattr(Kappas, '__len__') and not isinstance(Kappas, (str, dict)) and len(Kappas)==1:
                Kappas = Kappas[0]
            if not isFinite(Kappas):
                raise ValueError(f'VLQ.convertModel: kappa value to be converted must be a number -> {self.OldKappatext}')
            Kappas = toNumber(Kappas)
            converted = modelConverter.convert(self,Kappas,BRs,reverse,zeta)
            if self.KappaH is not None or self.KappaW is not None or self.KappaZ is not None:
                warnings.warn(f'VLQ.convertModel: deleting old kappa values for {self.NewKappatext}', stacklevel=2)
            self.zeta = zeta
            self.KappaOld = [Kappas]
            self.KappaOldConverted = None
            self.BRs = BRs
            self.Ratio = None
            self.setKappas(converted)
            self.SM = constants.inputs()
            # VLQ_UFO neglects the quark masses, so a channel can be open in VLQ_UFO but closed in VLQ_v4/VLQ_v5
            closed = set()
            for i, Mass in enumerate(self.Mass):
                unit = decays.partialWidths(self.type, self.FNS, Mass, {channel: 1 for channel in converted}, self.LR)
                closed |= {(channel, Mass) for channel, values in converted.items() if values[i] != 0 and unit[channel] <= 0}
            if closed:
                text = ', '.join(f'{self.type} -> {boson.upper()} + generation {gen} quark at M = {Mass} GeV'
                                 for (boson, gen), Mass in sorted(closed))
                warnings.warn(f'VLQ.convertModel: channels kinematically closed in VLQ_v4/VLQ_v5, whose couplings give '
                              f'no width there: {text}', stacklevel=2)
        else:
            if zeta is not None:
                raise ValueError('VLQ.convertModel: zeta is only used for VLQ_UFO -> VLQ_v4; '
                                 'for VLQ_v4 -> VLQ_UFO the generation couplings come from the VLQ object')
            n = len(self.Mass)
            def asList(values):
                # one value per mass
                values = [values] if isNumber(values) else values
                values = list(values) if hasattr(values, '__len__') and not isinstance(values, (str, dict)) else values
                if type(values) != list or len(values) != n or not all(isFinite(k) for k in values):
                    raise ValueError(f'VLQ.convertModel: kappa values to be converted must be numbers for {self.type} -> '
                                     f'{self.NewKappatext}, one value per mass ({n})')
                return [toNumber(k) for k in values]
            if self.type in ['VLX','VLY']:
                self.KappaW = asList(Kappas)
            else:
                if hasattr(Kappas, '__len__') and not isinstance(Kappas, (str, dict)):
                    Kappas = list(Kappas)
                if type(Kappas) != list or len(Kappas) != 3:
                    raise ValueError(f'VLQ.convertModel: {self.type} needs 3 kappa values (or lists) -> [{self.NewKappatext}]')
                self.KappaH, self.KappaW, self.KappaZ = [asList(ks) for ks in Kappas]
            if self.KappaOld is not None:
                warnings.warn(f'VLQ.convertModel: deleting old kappa values for {self.OldKappatext}', stacklevel=2)
            self.KappaOld = modelConverter.convert(self,None,reverse=reverse)
            self.KappaOldConverted = self.KappaOld
            self.SM = constants.inputs()
            # the couplings are given directly, so a Gamma/M from an earlier calcRatioKappas does not belong to them
            self.Ratio = None

    def calcRatioKappas(self,BRs=None,Ratio=None,zeta=None):
        '''
        Computes the couplings for each mass from the width-to-mass ratio and the branching ratios.

        BRs  : [BR(H), BR(W), BR(Z)] for T and B (sum 1); not used for X and Y
        Ratio: Gamma/M
        zeta : generation fractions [zeta1, zeta2, zeta3] (sum 1); default [0, 0, 1]
        '''
        if self.Mass is None:
            raise ValueError('VLQ.calcRatioKappas: Please define mass to VLQ object with setMass method')
        if not isFinite(Ratio) or Ratio <= 0:
            raise ValueError(f'VLQ.calcRatioKappas: Gamma/M({self.type}) ratio must be a positive number (int or float)')
        if self.type in ['VLX','VLY']:
            BRs = 1
        else:
            BRs = checkBRs(BRs, 'VLQ.calcRatioKappas',
                           f'{self.type} must have 3 float branching values for {self.type} -> [{self.BRtext}]')
        zeta = checkZeta(zeta, 'VLQ.calcRatioKappas')
        # the object is changed only if the calculation succeeds
        Kappas = functions.calcRatioKappas(self, BRs, toNumber(Ratio), zeta)
        self.Ratio = toNumber(Ratio)
        self.BRs = BRs
        self.zeta = zeta
        self.setKappas(Kappas)
        self.SM = constants.inputs()
        # a VLQ_UFO coupling from an earlier conversion does not belong to the new couplings
        self.KappaOld = None
        self.KappaOldConverted = None

    def table(self, style='md'):
        '''Returns the results as a Markdown (style='md') or aligned text (style='txt') table'''
        from VLQcalc import output
        return output.table(self, style)

    def save(self, filename):
        '''Writes the results to a file; the format follows the extension: .md, .txt, .csv or .json'''
        from VLQcalc import output
        output.save(self, filename)

    def plot(self, quantity='kappa', save=None, ax=None, logy=None, combined=False):
        '''
        Plots the couplings ('kappa'), branching ratios ('BR') or partial widths ('width') against the mass.

        save: file name of the figure; ax: existing axes (single panel only);
        logy: True or False for the y scale, chosen from the values if None;
        combined: with couplings to lighter quarks, True draws all channels in one panel
        (one colour per boson) instead of one panel per boson. Needs matplotlib.
        '''
        from VLQcalc import output
        return output.plot(self, quantity, save, ax, logy, combined)
