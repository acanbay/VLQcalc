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
class MG5:
    def __init__(self,VLQ,model='VLQ_v4_UFO'):
        self.VLQ = VLQ

        # Checking the model
        if model.upper() not in ['VLQ_UFO', 'VLQ_V4']:
            print('\n!! Error: VLQ !!\nmodel must be one of the following:\n')
            print('  VLQ_UFO')
            print('  VLQ_v4_UFO')
            exit()
        if model.upper() == 'VLQ_UFO':
            self.model = 'VLQ_UFO'
            self.text = 'import model VLQ_UFO\n'
        elif model.upper() == 'VLQ_V4':
            self.model = 'VLQ_v4_UFO'
            self.text = f'import model {self.model}_{VLQ.FNS}FNS_UFO-3rd'
            if VLQ.LR:
                self.text += '\n'
            else:
                self.text += 'L\n'

        self.Identication = ''
        self.Process = []
        self.shower = 'OFF' #OFF, pythia8
        self.detector = 'OFF' #OFF, Delphes
        self.analysis = 'OFF' #OFF, ExRoot, MadAnalysis
        self.madspin = 'OFF' #OFF, ON, onshell, full
        self.reweight = 'OFF' #OFF, ON
        self.runInputs = []

    def correctText(self,text):
        #W
        if 'w+' in text or 'w-' in text:
            if 'define WW' not in self.text:
                self.text += 'define WW = w+ w-\n'
            if 'w-' in text:
                text = text.replace('w-', 'WW')
            if 'w+' in text:
                text = text.replace('w+', 'WW')
        #t
        if 't ' in text or 't,' in text or 't)' in text or 't>' in text or 't~' in text:
            if 'define tt' not in self.text:
                self.text += 'define tt = t t~\n'
            if 't ' in text:
                text = text.replace('t ', 'tt ')
            if 't,' in text:
                text = text.replace('t,', 'tt,')
            if 't)' in text:
                text = text.replace('t)', 'tt)')
            if 't>' in text:
                text = text.replace('t>', 'tt>')
            if 't~' in text:
                text = text.replace('t~', 'tt')
        #b
        if 'b ' in text or 'b,' in text or 'b)' in text or 'b>' in text or 'b~' in text:
            if 'define bb' not in self.text:
                self.text += 'define bb = b b~\n'
            if 'b ' in text:
                text = text.replace('b ', 'bb ')
            if 'b,' in text:
                text = text.replace('b,', 'bb,')
            if 'b)' in text:
                text = text.replace('b)', 'bb)')
            if 'b>' in text:
                text = text.replace('b>', 'bb >')
            if 'b~' in text:
                text = text.replace('b~', 'bb')
        #tp
        if 'tp' in text or 'tp~' in text:
            if 'define VLT' not in self.text:
                self.text += 'define VLT = tp tp~\n'
            if 'tp' in text:
                text = text.replace('tp', 'VLT')
            if 'tp~' in text:
                text = text.replace('tp~', 'VLT')
        #bp
        if 'bp' in text or 'bp~' in text:
            if 'define VLB' not in self.text:
                self.text += 'define VLB = bp bp~\n'
            if 'bp' in text:
                text = text.replace('bp', 'VLB')
            if 'bp~' in text:
                text = text.replace('bp~', 'VLB')
        #x
        if 'x' in text or 'x~' in text:
            if 'define VLX' not in self.text:
                self.text += 'define VLX = x x~\n'
            if 'x' in text:
                text = text.replace('x', 'VLX')
            if 'x~' in text:
                text = text.replace('x~', 'VLX')
        #y
        if 'y' in text or 'y~' in text:
            if 'define VLY' not in self.text:
                self.text += 'define VLY = y y~\n'
            if 'y' in text:
                text = text.replace('y', 'VLY')
            if 'y~' in text:
                text = text.replace('y~', 'VLY')

        #e
        if 'e-' in text or 'e+' in text:
            if 'define ee' not in self.text:
                self.text += 'define ee = e- e+\n'
            if 'e-' in text:
                text = text.replace('e-', 'ee')
            if 'e+' in text:
                text = text.replace('e+', 'ee')

        #mu
        if 'mu-' in text or 'mu+' in text:
            if 'define mumu' not in self.text:
                self.text += 'define mumu = mu- mu+\n'
            if 'mu-' in text:
                text = text.replace('mu-', 'mumu')
            if 'mu+' in text:
                text = text.replace('mu+', 'mumu')

        #l
        if 'l-' in text or 'l+' in text:
            if 'define ll' not in self.text:
                self.text += 'define ll = l- l+\n'
            if 'l-' in text:
                text = text.replace('l-', 'll')
            if 'l+' in text:
                text = text.replace('l+', 'll')

        #ve
        if 've' in text or 've~' in text:
            if 'define veve' not in self.text:
                self.text += 'define veve = ve ve~\n'
            if 've' in text:
                text = text.replace('ve ', 'veve')
            if 've~' in text:
                text = text.replace('ve~', 'veve')

        #vm
        if 'vm' in text or 'vm~' in text:
            if 'define vmvm' not in self.text:
                self.text += 'define vmvm = vm vm~\n'
            if 'vm' in text:
                text = text.replace('vm ', 'vmvm')
            if 'vm~' in text:
                text = text.replace('vm~', 'vmvm')

        #vt
        if 'vt' in text or 'vt~' in text:
            if 'define vtvt' not in self.text:
                self.text += 'define vtvt = vt vt~\n'
            if 'vt' in text:
                text = text.replace('vt ', 'vtvt')
            if 'vt~' in text:
                text = text.replace('vt~', 'vtvt')

        #vl
        if 'vl' in text or 'vl~' in text:
            if 'define vlvl' not in self.text:
                self.text += 'define vlvl = vl vl~\n'
            if 'vl' in text:
                text = text.replace('vl ', 'vlvl')
            if 'vl~' in text:
                text = text.replace('vl~', 'vlvl')
        
        return text

    def splitProcess(self, process):
        proc= ''
        veto = ''
        decay = ''
        if '/' in process:
            idx_veto = process.find('/')
            proc = self.correctText( process[:idx_veto])

            if '(' in process[idx_veto:]:
                idx_decay = process.find('(')
                veto = self.correctText( process[idx_veto:idx_decay] )
                decay = self.correctText( process[idx_decay:] )
        
        elif '(' in process:
            idx_decay = process.find('(')
            proc = self.correctText( process[:idx_decay] )
            decay = self.correctText( process[idx_decay:] )

        else:
            proc = self.correctText( process )

        return proc+veto+decay

    def setProcess(self,process):
        if 'generate' not in process:
            process = 'generate '+process
        if '>' in process:
            check = True
        else:
            check = False
        #
        if not check:
            print('\n!! Error: setProcess !!\nWrong process\n')
            exit()
        #
        process = self.splitProcess(process)
        self.Process.append( process )

    def addProcess(self,process):
        if 'add process' not in process:
            process = 'add process '+process
        #
        if '>' in process:
            check = True
        else:
            check = False
        if not check:
            print('\n!! Error: addProcess !!\nWrong process\n')
            exit()
        #
        if len(self.Process) > 0:
            process = self.splitProcess(process)
            self.Process.append( process )
        else:
            print('\n!! Error: addProcess !!\naddProcess cannot be used without defining the main process\n')
            exit()
            
    def addInput(self, input):
        if 'set ' not in input:
            input = 'set '+input
        self.runInputs.append(input)

    def createMG5Input(self, name='mg5_input'):
        if self.Process == []:
            print('\n!! Error: createMG5Input !!\nProcess not defined\n')
            exit()
        
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
        if self.model[:-1] == 'VLQ_v4_UFO':
            if self.VLQ.type == 'VLT':
                if self.VLQ.Mass != None:
                    self.text += f'set MTP scan1:{list(self.VLQ.Mass)}\n'
                if self.VLQ.KappaH != None and self.VLQ.KappaW != None and self.VLQ.KappaZ != None:
                    if len(self.VLQ.KappaH)+len(self.VLQ.KappaW)+len(self.VLQ.KappaZ) != 3*len(self.VLQ.Mass):
                        self.text += f'set KTLh3 {self.VLQ.KappaH[0]}\n'
                        if self.VLQ.LR == True:
                            self.text += f'set KTRh3 {self.VLQ.KappaH[0]}\n'
                        #
                        self.text += f'set KTLw3 {self.VLQ.KappaW[0]}\n'
                        if self.VLQ.LR == True:
                            self.text += f'set KTRw3 {self.VLQ.KappaW[0]}\n'
                        #
                        self.text += f'set KTLz3 {self.VLQ.KappaZ[0]}\n'
                        if self.VLQ.LR == True:
                            self.text += f'set KTRz3 {self.VLQ.KappaZ[0]}\n'
                    else:
                        self.text += f'set KTLh3 scan1:{self.VLQ.KappaH}\n'
                        if self.VLQ.LR == True:
                            self.text += f'set KTRh3 scan1:{self.VLQ.KappaH}\n'
                        #
                        self.text += f'set KTLw3 scan1:{self.VLQ.KappaW}\n'
                        if self.VLQ.LR == True:
                            self.text += f'set KTRw3 scan1:{self.VLQ.KappaW}\n'
                        #
                        self.text += f'set KTLz3 scan1:{self.VLQ.KappaZ}\n'
                        if self.VLQ.LR == True:
                            self.text += f'set KTRz3 scan1:{self.VLQ.KappaZ}\n'

            elif self.VLQ.type == 'VLB':
                if self.VLQ.Mass != None:
                    self.text += f'set MBP scan1:{list(self.VLQ.Mass)}\n'
                if self.VLQ.KappaH != None and self.VLQ.KappaW != None and self.VLQ.KappaZ != None:
                    if len(self.VLQ.KappaH)+len(self.VLQ.KappaW)+len(self.VLQ.KappaZ) != 3*len(self.VLQ.Mass):
                        self.text += f'set KBLh3 {self.VLQ.KappaH[0]}\n'
                        if self.VLQ.LR == True:
                            self.text += f'set KBRh3 {self.VLQ.KappaH[0]}\n'
                        #
                        self.text += f'set KBLw3 {self.VLQ.KappaW[0]}\n'
                        if self.VLQ.LR == True:
                            self.text += f'set KBRw3 {self.VLQ.KappaW[0]}\n'
                        #
                        self.text += f'set KBLz3 {self.VLQ.KappaZ[0]}\n'
                        if self.VLQ.LR == True:
                            self.text += f'set KBRz3 {self.VLQ.KappaZ[0]}\n'
                    else:
                        self.text += f'set KBLh3 scan1:{self.VLQ.KappaH}\n'
                        if self.VLQ.LR == True:
                            self.text += f'set KBRh3 scan1:{self.VLQ.KappaH}\n'
                        #
                        self.text += f'set KBLw3 scan1:{self.VLQ.KappaW}\n'
                        if self.VLQ.LR == True:
                            self.text += f'set KBRw3 scan1:{self.VLQ.KappaW}\n'
                        #
                        self.text += f'set KBLz3 scan1:{self.VLQ.KappaZ}\n'
                        if self.VLQ.LR == True:
                            self.text += f'set KBRz3 scan1:{self.VLQ.KappaZ}\n'

            elif self.VLQ.type == 'VLX':
                if self.VLQ.Mass != None:
                    self.text += f'set MX scan1:{list(self.VLQ.Mass)}\n'
                if self.VLQ.KappaW != None:
                    if len(self.VLQ.KappaW) != len(self.VLQ.Mass):
                        self.text += f'set KXL3 {self.VLQ.KappaW[0]}\n'
                        if self.VLQ.LR == True:
                            self.text += f'set KXR3 {self.VLQ.KappaW[0]}\n'
                    else:
                        self.text += f'set KXL3 scan1:{self.VLQ.KappaW}\n'
                        if self.VLQ.LR == True:
                            self.text += f'set KXR3 scan1:{self.VLQ.KappaW}\n'

            elif self.VLQ.type == 'VLY':
                if self.VLQ.Mass != None:
                    self.text += f'set MY scan1:{list(self.VLQ.Mass)}\n'
                if self.VLQ.KappaW != None:
                    if len(self.VLQ.KappaW) != len(self.VLQ.Mass):
                        self.text += f'set KYL3 {self.VLQ.KappaW[0]}\n'
                        if self.VLQ.LR == True:
                            self.text += f'set KYR3 {self.VLQ.KappaW[0]}\n'
                    else:
                        self.text += f'set KYL3 scan1:{self.VLQ.KappaW}\n'
                        if self.VLQ.LR == True:
                            self.text += f'set KYR3 scan1:{self.VLQ.KappaW}\n'
        
        else:
            if self.VLQ.type == 'VLT':
                self.text += 'set zetaTuL 0.000000e-01\n'
                self.text += 'set zetaTuR 0.000000e-01\n'
                self.text += 'set zetaTcL 0.000000e-01\n'
                self.text += 'set zetaTcR 0.000000e-01\n'
                self.text += 'set zetaTtL 1.000000e+00\n'
                if self.VLQ.LR:
                    self.text += 'set zetaTtR 1.000000e-01\n'
                else:
                    self.text += 'set zetaTtR 0.000000e-01\n'
                self.text += f'set xitph {self.VLQ.BRs[0]}\n'
                self.text += f'set xitpw {self.VLQ.BRs[1]}\n'
                self.text += f'set xitpz {self.VLQ.BRs[2]}\n'
                if self.VLQ.Mass != None:
                    self.text += f'set MTP scan1:{list(self.VLQ.Mass)}\n'
                if self.VLQ.KappaOld != None:
                    if len(self.VLQ.KappaOld) != len(self.VLQ.Mass): 
                        self.text += f'set KT {self.VLQ.KappaOld[0]}\n'
                    else:
                        self.text += f'set KT scan1:{self.VLQ.KappaOld}\n'

            elif self.VLQ.type == 'VLB':
                self.text += 'set zetaBdL 0.000000e-01\n'
                self.text += 'set zetaBdR 0.000000e-01\n'
                self.text += 'set zetaBsL 0.000000e-01\n'
                self.text += 'set zetaBsR 0.000000e-01\n'
                self.text += 'set zetaBbL 1.000000e+00\n'
                if self.VLQ.LR:
                    self.text += 'set zetaBbR 1.000000e-01\n'
                else:
                    self.text += 'set zetaBbR 0.000000e-01\n'
                self.text += f'set xibph {self.VLQ.BRs[0]}\n'
                self.text += f'set xibpw {self.VLQ.BRs[1]}\n'
                self.text += f'set xibpz {self.VLQ.BRs[2]}\n'
                if self.VLQ.Mass != None:
                    self.text += f'set MBP scan1:{list(self.VLQ.Mass)}\n'
                if self.VLQ.KappaOld != None:
                    if len(self.VLQ.KappaOld) != len(self.VLQ.Mass): 
                        self.text += f'set KB {self.VLQ.KappaOld[0]}\n'
                    else:
                        self.text += f'set KB scan1:{self.VLQ.KappaOld}\n'
                
            elif self.VLQ.type == 'VLX':
                self.text += 'set zetaXuL 0.000000e-01\n'
                self.text += 'set zetaXuR 0.000000e-01\n'
                self.text += 'set zetaXcL 0.000000e-01\n'
                self.text += 'set zetaXcR 0.000000e-01\n'
                self.text += 'set zetaXtL 1.000000e+00\n'
                if self.VLQ.LR:
                    self.text += 'set zetaXtR 1.000000e-01\n'
                else:
                    self.text += 'set zetaXtR 0.000000e-01\n'
                if self.VLQ.Mass != None:
                    self.text += f'set MX scan1:{list(self.VLQ.Mass)}\n'
                if self.VLQ.KappaOld != None:
                    if len(self.VLQ.KappaOld) != len(self.VLQ.Mass): 
                        self.text += f'set KX {self.VLQ.KappaOld[0]}\n'
                    else:
                        self.text += f'set KX scan1:{self.VLQ.KappaOld}\n'

            elif self.VLQ.type == 'VLY':
                self.text += 'set zetaYdL 0.000000e-01\n'
                self.text += 'set zetaYdR 0.000000e-01\n'
                self.text += 'set zetaYsL 0.000000e-01\n'
                self.text += 'set zetaYsR 0.000000e-01\n'
                self.text += 'set zetaYbL 1.000000e+00\n'
                if self.VLQ.LR:
                    self.text += 'set zetaYbR 1.000000e-01\n'
                else:
                    self.text += 'set zetaYbR 0.000000e-01\n'
                if self.VLQ.Mass != None:
                    self.text += f'set MY scan1:{list(self.VLQ.Mass)}\n'
                if self.VLQ.KappaOld != None:
                    if len(self.VLQ.KappaOld) != len(self.VLQ.Mass): 
                        self.text += f'set KY {self.VLQ.KappaOld[0]}\n'
                    else:
                        self.text += f'set KY scan1:{self.VLQ.KappaOld}\n'
            
        self.text += 'done\n'

        with open(f'{name}.dat', 'w') as f:
            f.write(self.text)