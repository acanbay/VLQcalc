'''
Author : Ali Can Canbay - acanbay
Created: 21 May 2024
Updated: 25 October 2026

This module is used to create MG5 object. MG5 object contains methods to create input file for 
MadGraph5_aMC@NLO [1].

====REFERENCES====
[1] Alwall, J., Herquet, M., Maltoni, F., Mattelaer, O., & Stelzer, T. (2011). 
    MadGraph 5: going beyond. 
    Journal of High Energy Physics, 2011(6), 128.
    https://launchpad.net/mg5amcnlo
'''
import re, os, warnings

# Multiparticle aliases written to the MG5 card
ALIASES = {
    'WW': ['w+', 'w-'], 'tt': ['t', 't~'], 'bb': ['b', 'b~'],
    'VLT': ['tp', 'tp~'], 'VLB': ['bp', 'bp~'], 'VLX': ['x', 'x~'], 'VLY': ['y', 'y~'],
    'ee': ['e-', 'e+'], 'mumu': ['mu-', 'mu+'], 'll': ['l-', 'l+'],
    'veve': ['ve', 've~'], 'vmvm': ['vm', 'vm~'], 'vtvt': ['vt', 'vt~'], 'vlvl': ['vl', 'vl~'],
}
LABELS = {p: alias for alias, ps in ALIASES.items() for p in ps}
# A label matches only as a whole particle name, so "tp~" is never read as "t"
PATTERN = re.compile(
    r'(?<![\w~+-])('
    + '|'.join(re.escape(p) for p in sorted(LABELS, key=len, reverse=True))
    + r')(?![\w~+-])')

# Parameter names of the VLQ_UFO model
OLD_MODEL = {'VLT': {'mass':'MTP', 'kappa':'KT', 'zeta':['zetaTu','zetaTc','zetaTt'], 'xi':['xitph','xitpw','xitpz']},
             'VLB': {'mass':'MBP', 'kappa':'KB', 'zeta':['zetaBd','zetaBs','zetaBb'], 'xi':['xibph','xibpw','xibpz']},
             'VLX': {'mass':'MX',  'kappa':'KX', 'zeta':['zetaXu','zetaXc','zetaXt'], 'xi':[]},
             'VLY': {'mass':'MY',  'kappa':'KY', 'zeta':['zetaYd','zetaYs','zetaYb'], 'xi':[]}}

class MG5:
    '''
    MadGraph5 input card for a VLQ object: MG5(VLQ, model='VLQ_v4_UFO')

    model: 'VLQ_UFO', 'VLQ_v4_UFO' or 'VLQ_v5_UFO'
    Run settings: shower, detector, analysis, madspin, reweight (default 'OFF')
    '''
    def __init__(self,VLQ,model='VLQ_v4_UFO'):
        if getattr(VLQ, 'type', None) not in OLD_MODEL:
            raise ValueError('MG5: the first argument must be a VLQ object (VLQcalc.model.VLQ)')
        self.VLQ = VLQ

        # Checking the model
        if type(model) != str or model.upper() not in ['VLQ_UFO', 'VLQ_V4_UFO', 'VLQ_V5_UFO']:
            raise ValueError('MG5: model must be one of the following: VLQ_UFO, VLQ_v4_UFO, VLQ_v5_UFO')
        if model.upper() == 'VLQ_UFO':
            self.model = 'VLQ_UFO'
            self.text = 'import model VLQ_UFO\n'
        elif model.upper() == 'VLQ_V4_UFO':
            self.model = 'VLQ_v4'
            # The full model is imported (5FNS has no restriction cards, and couplings to all generations
            # can be used); every coupling of the VLQ is written to the card in createMG5Input
            self.text = f'import model {self.model}_{VLQ.FNS}FNS_UFO\n'
        elif model.upper() == 'VLQ_V5_UFO':
            # VLQ_v5 has the same couplings to SM quarks as VLQ_v4 (and the same decay formulas). Its couplings
            # between two VLQs (KWTBL, KZTTL, ...) are zero by default and are not changed here.
            self.model = 'VLQ_v5'
            self.text = f'import model {self.model}_{VLQ.FNS}FNS_NLO_UFO\n'

        self.Identication = ''
        self.Process = []
        self.shower = 'OFF' #OFF, pythia8
        self.detector = 'OFF' #OFF, Delphes
        self.analysis = 'OFF' #OFF, ExRoot, MadAnalysis
        self.madspin = 'OFF' #OFF, ON, onshell, full
        self.reweight = 'OFF' #OFF, ON
        self.runInputs = []

    def correctText(self,text):
        '''Replaces every particle name with its particle+antiparticle label and defines the label once'''
        def replace(match):
            alias = LABELS[match.group(1)]
            if f'define {alias} ' not in self.text:
                self.text += f'define {alias} = {" ".join(ALIASES[alias])}\n'
            return alias
        return PATTERN.sub(replace, text)

    def splitProcess(self, process):
        # Kept for backward compatibility; labels are matched as whole names, so no splitting is needed
        return self.correctText(process)

    def setProcess(self,process):
        '''Defines the main process, e.g. 'p p > tp tp~' '''
        if type(process) != str:
            raise ValueError(f'MG5.setProcess: process must be a string -> {process}')
        process = process.strip()
        if not process.startswith('generate '):
            process = 'generate '+process
        if '>' not in process:
            raise ValueError(f'MG5.setProcess: Wrong process -> {process}')
        process = self.splitProcess(process)
        self.Process.append( process )

    def addProcess(self,process):
        '''Adds a process after the main one'''
        if type(process) != str:
            raise ValueError(f'MG5.addProcess: process must be a string -> {process}')
        process = process.strip()
        if not process.startswith('add process '):
            process = 'add process '+process
        if '>' not in process:
            raise ValueError(f'MG5.addProcess: Wrong process -> {process}')
        if len(self.Process) == 0:
            raise ValueError('MG5.addProcess: addProcess cannot be used without defining the main process')
        process = self.splitProcess(process)
        self.Process.append( process )

    def addInput(self, input):
        '''Adds a setting for the run cards, e.g. 'nevents 10000' '''
        if type(input) != str:
            raise ValueError(f'MG5.addInput: input must be a string -> {input}')
        input = input.strip()
        if not input.startswith('set '):
            input = 'set '+input
        self.runInputs.append(input)

    def setLine(self, name, values):
        # A single value, or one value per mass written as a scan; numbers with 6 significant digits
        from VLQcalc.model import isNumber, toNumber
        def number(value):
            value = toNumber(value)
            return str(value) if isinstance(value, int) else f'{value:.6g}'
        values = [values] if isNumber(values) else list(values)
        if all(value == 0 for value in values):
            return f'set {name} 0\n'
        if len(values) == 1 or all(value == values[0] for value in values):
            return f'set {name} {number(values[0])}\n'
        if self.VLQ.Mass is None or len(values) != len(self.VLQ.Mass):
            raise ValueError(f'MG5.createMG5Input: {name} has {len(values)} values; give one value per mass')
        return f'set {name} scan1:[{", ".join(number(value) for value in values)}]\n'

    def couplingNames(self):
        # Names of all couplings of the VLQ to SM quarks in VLQ_v4/VLQ_v5, e.g. KTLh1, ..., KTRz3
        letter = self.VLQ.type[-1]
        names = []
        for boson in self.VLQ.bosons:
            for gen in [1,2,3]:
                name = boson if self.VLQ.type in ['VLT','VLB'] else ''
                names += [f'K{letter}L{name}{gen}', f'K{letter}R{name}{gen}']
        return names

    def createMG5Input(self, name='mg5_input'):
        '''Writes the card to name.dat'''
        if self.Process == []:
            raise ValueError('MG5.createMG5Input: Process not defined')
        name = os.fspath(name)
        if name.endswith('.dat'):
            name = name[:-4]
        VLQ = self.VLQ
        # The card is built in a local string and written only if it is complete
        text = self.text
        for process in self.Process:
            text += process+'\n'

        if name!='mg5_input':
            text += f'output {name}\n'
            text += f'launch {name}\n'
        else:
            text += 'launch\n'

        ###
        text += f'shower={self.shower}\n'
        text += f'detector={self.detector}\n'
        text += f'analysis={self.analysis}\n'
        text += f'madspin={self.madspin}\n'
        text += f'reweight={self.reweight}\n'
        text += 'done\n'

        ###############################
        names = OLD_MODEL[VLQ.type]
        if self.model in ['VLQ_v4', 'VLQ_v5']:
            if VLQ.Mass is not None:
                text += self.setLine(names['mass'], list(VLQ.Mass))
            # Every coupling of the VLQ is written, since the defaults of the full model are not zero.
            # The 3rd generation couplings are skipped only if they were never set (model defaults are kept).
            third = {'h': VLQ.KappaH, 'w': VLQ.KappaW, 'z': VLQ.KappaZ}
            Kappas = VLQ.Kappas
            for parameter in self.couplingNames():
                boson = parameter[3] if VLQ.type in ['VLT','VLB'] else 'w'
                if parameter.endswith('3') and third[boson] is None:
                    continue
                text += self.setLine(parameter, Kappas[parameter]) if parameter in Kappas else f'set {parameter} 0\n'

        else:
            # zeta: generation fractions; with LR they are shared equally between both chiralities
            share = 0.5 if VLQ.LR else 1.
            for zetaName, zeta in zip(names['zeta'], VLQ.zeta):
                text += self.setLine(f'{zetaName}L', zeta*share)
                text += self.setLine(f'{zetaName}R', zeta*share if VLQ.LR else 0)
            if names['xi']:
                if VLQ.BRs is None:
                    raise ValueError(f'MG5.createMG5Input: branchings (xi) of {VLQ.type} are not defined for VLQ_UFO')
                for xi, BR in zip(names['xi'], VLQ.BRs):
                    text += self.setLine(xi, BR)
            if VLQ.Mass is not None:
                text += self.setLine(names['mass'], list(VLQ.Mass))
            if VLQ.KappaOld is not None:
                text += self.setLine(names['kappa'], VLQ.KappaOld)
            else:
                warnings.warn(f'MG5.createMG5Input: no VLQ_UFO coupling ({names["kappa"]}) is set, so the default of the '
                              'model is used; set VLQ.KappaOld or use convertModel(..., reverse=True)', stacklevel=2)

        # Settings given with addInput come last, so that they can also change a mass or a coupling
        for input in self.runInputs:
            text += input+'\n'

        text += 'done\n'

        with open(f'{name}.dat', 'w') as f:
            f.write(text)

        print(f'\n{name}.dat is created!\n')
