'''
Author: Ali Can Canbay - acanbay
Date  : 21 May 2024

This module is used to create MG5 object. MG5 object contains methods to create input file for 
MadGraph5_aMC@NLO [1].

====REFERENCES====
[1] Alwall, J., Herquet, M., Maltoni, F., Mattelaer, O., & Stelzer, T. (2011). 
    MadGraph 5: going beyond. 
    Journal of High Energy Physics, 2011(6), 1-40.
    https://launchpad.net/mg5amcnlo
'''
import re

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
        if 'generate' not in process:
            process = 'generate '+process
        if '>' not in process:
            raise ValueError(f'MG5.setProcess: Wrong process -> {process}')
        process = self.splitProcess(process)
        self.Process.append( process )

    def addProcess(self,process):
        '''Adds a process after the main one'''
        if 'add process' not in process:
            process = 'add process '+process
        if '>' not in process:
            raise ValueError(f'MG5.addProcess: Wrong process -> {process}')
        if len(self.Process) == 0:
            raise ValueError('MG5.addProcess: addProcess cannot be used without defining the main process')
        process = self.splitProcess(process)
        self.Process.append( process )

    def addInput(self, input):
        '''Adds a setting for the run cards, e.g. 'nevents 10000' '''
        if 'set ' not in input:
            input = 'set '+input
        self.runInputs.append(input)

    def setLine(self, name, values):
        # One value per mass is written as a scan, otherwise the first value is used
        if type(values) != list:
            return f'set {name} {values}\n'
        if all(value == 0 for value in values):
            return f'set {name} 0\n'
        if len(values) == len(self.VLQ.Mass):
            return f'set {name} scan1:{values}\n'
        return f'set {name} {values[0]}\n'

    def createMG5Input(self, name='mg5_input'):
        '''Writes the card to name.dat'''
        if self.Process == []:
            raise ValueError('MG5.createMG5Input: Process not defined')
        # Header (model import and defines) is restored at the end, so repeated calls do not duplicate the card
        header = self.text

        for process in self.Process:
            self.text += process+'\n'

        if name!='mg5_input':
            self.text += f'output {name}\n'
            self.text += f'launch {name}\n'
        else:
            self.text += 'launch\n'

        ###
        self.text += f'shower={self.shower}\n'
        self.text += f'detector={self.detector}\n'
        self.text += f'analysis={self.analysis}\n'
        self.text += f'madspin={self.madspin}\n'
        self.text += f'reweight={self.reweight}\n'
        self.text += 'done\n'

        ###Cards
        for input in self.runInputs:
            self.text += input+'\n'
        
        ###############################
        names = OLD_MODEL[self.VLQ.type]
        if self.model in ['VLQ_v4', 'VLQ_v5']:
            if self.VLQ.Mass != None:
                self.text += self.setLine(names['mass'], list(self.VLQ.Mass))
            # Every coupling of the VLQ is written, since the defaults of the full model are not zero.
            # The 3rd generation couplings are skipped only if they were never set (model defaults are kept).
            third = {'h': self.VLQ.KappaH, 'w': self.VLQ.KappaW, 'z': self.VLQ.KappaZ}
            if self.VLQ.Mass != None:
                for parameter, values in self.VLQ.Kappas.items():
                    boson = parameter[3] if self.VLQ.type in ['VLT','VLB'] else 'w'
                    if parameter.endswith('3') and third[boson] == None:
                        continue
                    self.text += self.setLine(parameter, values)
            else:
                for parameter in [p for p in self.VLQ.Kappas if not p.endswith('3')]:
                    self.text += f'set {parameter} 0\n'
        
        else:
            # zeta: generation fractions; with LR they are shared equally between both chiralities
            share = 0.5 if self.VLQ.LR else 1.
            for zetaName, zeta in zip(names['zeta'], self.VLQ.zeta):
                self.text += f'set {zetaName}L {zeta*share:.6e}\n'
                self.text += f'set {zetaName}R {zeta*share if self.VLQ.LR else 0:.6e}\n'
            if names['xi']:
                if self.VLQ.BRs == None:
                    raise ValueError(f'MG5.createMG5Input: branchings (xi) of {self.VLQ.type} are not defined for VLQ_UFO')
                for xi, BR in zip(names['xi'], self.VLQ.BRs):
                    self.text += self.setLine(xi, BR)
            if self.VLQ.Mass != None:
                self.text += self.setLine(names['mass'], list(self.VLQ.Mass))
            if self.VLQ.KappaOld != None:
                self.text += self.setLine(names['kappa'], self.VLQ.KappaOld)
            
        self.text += 'done\n'

        with open(f'{name}.dat', 'w') as f:
            f.write(self.text)
        self.text = header
        
        print(f'\n{name}.dat is created!\n')
