"""
Author: Ali Can Canbay - acanbay
Date  : 21 May 2024

This module contains the VLQ object. It debugs commands and communicates with other modules.
"""
from VLQcalc import decays
from VLQcalc import modelConverter
from VLQcalc import functions
############
class VLQ:
    def __init__( self, VLQtype, FNS=4, LR=False ):
        # Checking the VLQ type
        if type(VLQtype)!=str:
            print("\n!! Error: VLQ !!\nVLQ type must be string\n")
            exit()
        if VLQtype.lower() not in ["t", "vlt", "tp", "b", "vlb", "bp", "x", "vlx", "y", "vly"]:
            print("\n!! Error: VLQ !!\nVLQ type must be one of the following:\n")
            print("  tp (charge:  2/3)")
            print("  bp (charge: -1/3)")
            print("  x  (charge:  5/3)")
            print("  y  (charge: -4/3)\n")
            exit()
        if VLQtype.lower() in ["t", "vlt", "tp"]:
            self.type = "VLT"
            self.BRtext = "BR(Ht),BR(Wb),BR(Zt)"
            self.NewKappatext = "kH,kW,kZ"
            self.OldKappatext = "kT"
        elif VLQtype.lower() in ["b", "vlb", "bp"]:
            self.type = "VLB"
            self.BRtext = "BR(Hb),BR(Wt),BR(Zb)"
            self.NewKappatext = "kH,kW,kZ"
            self.OldKappatext = "kB"
        elif VLQtype.lower() in ["x", "vlx"]:
            self.type = "VLX"
            self.BRtext = "BR(Wb)"
            self.NewKappatext = "kW"
            self.OldKappatext = "kX"
        elif VLQtype.lower() in ["y", "vly"]:
            self.type = "VLY"
            self.BRtext = "BR(Wt)"
            self.NewKappatext = "kW"
            self.OldKappatext = "kY"

        # Checking the FNS value
        if FNS not in [4,5]:
            print("\n!! Error: VLQ !!\nFlavor value must be 4 or 5\n")
            exit()

        if FNS == 4:
            if self.type == "VLT":
                self.calcDecay = decays.decayTP4
            elif self.type == "VLB":
                self.calcDecay = decays.decayBP4
            elif self.type == "VLX":
                self.calcDecay = decays.decayX4
            elif self.type == "VLY":
                self.calcDecay = decays.decayY4
        elif FNS == 5:
            if self.type == "VLT":
                self.calcDecay = decays.decayTP5
            elif self.type == "VLB":
                self.calcDecay = decays.decayBP5
            elif self.type == "VLX":
                self.calcDecay = decays.decayX5
            elif self.type == "VLY":
                self.calcDecay = decays.decayY5

        # Checking the LR value
        if type(LR)!=bool:
            print("\n!! Error: VLQ !!\nLR must be True (both left and right handed) or False (only left handed)\n")
            exit()
        else:
            self.LR = LR

        self.Mass = None
        self.Branchings = None
        self.KappaOld = None
        self.KappaH = None
        self.KappaW = None
        self.KappaZ = None
        self.FNS = FNS
        self.BRs = None

    def setMass(self, Masses):
        # Setting Mass
        check = True
        if type(Masses) in [int,float]:
            Masses = [Masses]
        elif type(Masses) in [list,range]:
            for Mass in Masses:
                if type(Mass) not in [int,float]:
                    check = False
                    break
        else:
            check = False
        if not check:
            print("\n!! Error: VLQ.setMass !!\nMass must be number or number list\n")
            exit()
        self.Mass = Masses

    def convertModel(self,Kappas,BRs=1,reverse=False):
        # Converts kappas VLQ_UFO <-> VLQ_v4
        if self.Mass == None:
            print(f"\n!! Error: VLQ.convertModel !!\nIt was not specified for which mass the new kappas would be calculated.")
            print(f"Please define mass to VLQ object with setMass method\n")
            exit()
        ##
        if type(reverse)!=bool:
            print("\n!! Error: VLQ.convertModel !!\nSource must be False (old->new) or True (new->False)'\n")
            exit()
        ##
        check=True
        if not reverse:
            if self.type in ["VLX","VLY"]:
                BRs = 1
            else:
                if type(BRs)==list:
                    if len(BRs)==3:
                        for BR in BRs:  
                            if type(BR) not in [int,float]:
                                check = False
                                break
                    else:
                        check=False
                else:
                    check=False
            if not check:
                print(f"\n!! Error: VLQ.convertModel !!\n{self.type} must have 3 float branching values (xi) in old model -> [xiH, xiW, xiZ]\n")
                exit()
            if type(BRs)==list:
                if sum(BRs)!=1:
                    print(f"\n!! Error: VLQ.convertModel !!\nFor {self.type} , the sum of all branchings (xiH+xiW+xiZ) must be 1\n")
                    exit()
            if type(Kappas) == list:
                if len(Kappas)==1:
                    Kappas=Kappas[0]
            if type(Kappas) not in [int,float]:
                print(f"\n!! Error: VLQ.convertModel !!\nkappa value to be converted must be a number -> {self.OldKappatext}\n")
                exit()
            self.KappaOld = Kappas
            if self.type in ["VLX","VLY"]:
                if self.KappaH!=None or self.KappaW!=None or self.KappaZ!=None:
                    print(f"\nWarning !!\nDeleting old kappa values for kW\n")
                self.KappaW = modelConverter.convert(self,Kappas,reverse=reverse)
            else:
                if self.KappaH!=None or self.KappaW!=None or self.KappaZ!=None:
                    print(f"\nWarning !!\nDeleting old kappa values for kH, kW and kZ\n")
                KappaHs, KappaWs, kappaZs = modelConverter.convert(self,Kappas,BRs,reverse)
                self.KappaH = KappaHs
                self.KappaW = KappaWs
                self.KappaZ = kappaZs
        else:
            if len(self.Mass)==1:
                if self.type in ["VLX","VLY"]:
                    if type(Kappas) in [int,float]:
                        Kappas=[Kappas]
                    elif type(Kappas)==list:
                        if len(Kappas)!=len(self.Mass):
                            check=False
                        elif type(Kappas[0]) not in [int,float]:
                            check=False
                    else:
                        check=False
                else:
                    if type(Kappas) in [int,float]:
                        check=False
                    elif type(Kappas)==list:
                        if len(Kappas)!=3:
                            check=False
                        else:
                            for k in Kappas:
                                if type(k) not in [int,float]:
                                    check=False
                                    break
                    else:
                        check=False
                    Kappas = [ [Kappas[0]],[Kappas[1]],[Kappas[2]] ]
                if not check:
                    print(f"\n!! Error: VLQ.convertModel !!\nkappa value to be converted must be a number for {self.type} -> {self.NewKappatext}\n")
                    exit()
            else:
                if self.type in ["VLX","VLY"]:
                    if type(Kappas) in [int,float]:
                        check=False
                    elif type(Kappas)==list:
                        if len(Kappas)!=len(self.Mass):
                            check=False
                        else:
                            for k in Kappas:
                                if type(k) not in [int,float]:
                                    check=False
                                    break
                    else:
                        check=False
                else:
                    if type(Kappas) in [int,float]:
                        check=False
                    elif type(Kappas)==list:
                        if len(Kappas)!=3:
                            check=False
                        else:
                            for ks in Kappas:
                                if len(ks)!=len(self.Mass):
                                    check=False
                                    break
                                for k in ks:
                                    if type(k) not in [int,float]:
                                        check=False
                                        break
                    else:
                        check=False
                if not check:
                    print(f"\n!! Error: VLQ.convertModel !!\nThe number of {self.NewKappatext} should be equal to the number of mass to be scanned\n")
                    exit()
            if self.KappaOld!=None:
                print(f"\nWarning !\nDeleting old kappa values for {self.OldKappatext}\n")
            self.KappaOld = modelConverter.convert(self,Kappas,reverse=reverse)

    def calcRatioKappas(self,BRs=None,Ratio=None):
        # Calculates kappas based on Gamma/M ratios
        check=True
        if type(Ratio) in [int,float]:
                self.Ratio=Ratio
        else:
            print(f"\n!! Error: VLQ.calcRatioKappas !!\nGamma/M({self.type}) ratio must be number (int or float)\n")
            exit()
        #
        check=True
        if self.type in ["VLX","VLY"]:
            BRs = 1
        else:
            if type(BRs)==list:
                if len(BRs)==3:
                    for BR in BRs:  
                        if type(BR) not in [int,float]:
                            check = False
                            break
                else:
                    check=False
            else:
                check=False
        if not check:
            print(f"\n!! Error: VLQ.calcRatioKappas !!\n{self.type} must have 3 float branching values for {self.type} -> [{self.BRtext}]\n")
            exit()
        if type(BRs)==list:
            if sum(BRs)!=1:
                print(f"\n!! Error: VLQ.calcRatioKappas !!\nFor {self.type} , the sum of all branchings ( {self.BRtext} ) must be 1\n")
                exit()
        self.BRs = BRs
        functions.calcRatioKappas(self)