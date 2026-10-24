'''
Author: Ali Can Canbay - acanbay
Date  : 21 May 2024

Calculates decay widths for VLQ_v4_4FNS_UFO and VLQ_v4_5FNS_UFO models [1,2].
All formulas are taken directly from the models; SM inputs are read from constants.py.

====REFERENCES====
[1] Alloul, A., Christensen, N. D., Degrande, C., Duhr, C., & Fuks, B. (2014). 
    FeynRules 2.0—A complete toolbox for tree-level phenomenology. 
    Computer Physics Communications, 185(8), 2250-2300.
    https://feynrules.irmp.ucl.ac.be/
[2] Fuks, B., & Shao, H. S. (2017). 
    QCD next-to-leading-order predictions matched to parton showers for vector-like quark models. 
    The European Physical Journal C, 77(2), 1-21.
    https://feynrules.irmp.ucl.ac.be/wiki/NLOModels
'''
import cmath
from VLQcalc import constants as c
#######
# 4FNS
#######
#VLT
def decayTP4(MtP,Kappas,LR=False):
    KTLh3 = Kappas[0]
    KTLw3 = Kappas[1]
    KTLz3 = Kappas[2]
    if LR == True:
        KTRh3 = Kappas[0]
        KTRw3 = Kappas[1]
        KTRz3 = Kappas[2]
    else:
        KTRh3 = 0
        KTRw3 = 0
        KTRz3 = 0
    decayTP_Ht = ((-3*KTLh3**2*c.MH**2 - 3*KTRh3**2*c.MH**2 + 3*KTLh3**2*c.Mt**2 + 3*KTRh3**2*c.Mt**2 + 12*KTLh3*KTRh3*c.Mt*MtP + 3*KTLh3**2*MtP**2 + 3*KTRh3**2*MtP**2)*cmath.sqrt(c.MH**4 - 2*c.MH**2*c.Mt**2 + c.Mt**4 - 2*c.MH**2*MtP**2 - 2*c.Mt**2*MtP**2 + MtP**4))/(96.*cmath.pi*abs(MtP)**3)
    decayTP_Wb = (((3*c.ee**2*KTLw3**2*c.Mb**2)/(2.*c.sw**2) + (3*c.ee**2*KTRw3**2*c.Mb**2)/(2.*c.sw**2) - (18*c.ee**2*KTLw3*KTRw3*c.Mb*MtP)/c.sw**2 + (3*c.ee**2*KTLw3**2*MtP**2)/(2.*c.sw**2) + (3*c.ee**2*KTRw3**2*MtP**2)/(2.*c.sw**2) + (3*c.ee**2*KTLw3**2*c.Mb**4)/(2.*c.MW**2*c.sw**2) + (3*c.ee**2*KTRw3**2*c.Mb**4)/(2.*c.MW**2*c.sw**2) - (3*c.ee**2*KTLw3**2*c.Mb**2*MtP**2)/(c.MW**2*c.sw**2) - (3*c.ee**2*KTRw3**2*c.Mb**2*MtP**2)/(c.MW**2*c.sw**2) + (3*c.ee**2*KTLw3**2*MtP**4)/(2.*c.MW**2*c.sw**2) + (3*c.ee**2*KTRw3**2*MtP**4)/(2.*c.MW**2*c.sw**2) - (3*c.ee**2*KTLw3**2*c.MW**2)/c.sw**2 - (3*c.ee**2*KTRw3**2*c.MW**2)/c.sw**2)*cmath.sqrt(c.Mb**4 - 2*c.Mb**2*MtP**2 + MtP**4 - 2*c.Mb**2*c.MW**2 - 2*MtP**2*c.MW**2 + c.MW**4))/(96.*cmath.pi*abs(MtP)**3)
    decayTP_Zt = (((3*c.ee**2*KTLz3**2*c.Mt**2)/(4.*c.cw**2*c.sw**2) + (3*c.ee**2*KTRz3**2*c.Mt**2)/(4.*c.cw**2*c.sw**2) - (9*c.ee**2*KTLz3*KTRz3*c.Mt*MtP)/(c.cw**2*c.sw**2) + (3*c.ee**2*KTLz3**2*MtP**2)/(4.*c.cw**2*c.sw**2) + (3*c.ee**2*KTRz3**2*MtP**2)/(4.*c.cw**2*c.sw**2) + (3*c.ee**2*KTLz3**2*c.Mt**4)/(4.*c.cw**2*c.MZ**2*c.sw**2) + (3*c.ee**2*KTRz3**2*c.Mt**4)/(4.*c.cw**2*c.MZ**2*c.sw**2) - (3*c.ee**2*KTLz3**2*c.Mt**2*MtP**2)/(2.*c.cw**2*c.MZ**2*c.sw**2) - (3*c.ee**2*KTRz3**2*c.Mt**2*MtP**2)/(2.*c.cw**2*c.MZ**2*c.sw**2) + (3*c.ee**2*KTLz3**2*MtP**4)/(4.*c.cw**2*c.MZ**2*c.sw**2) + (3*c.ee**2*KTRz3**2*MtP**4)/(4.*c.cw**2*c.MZ**2*c.sw**2) - (3*c.ee**2*KTLz3**2*c.MZ**2)/(2.*c.cw**2*c.sw**2) - (3*c.ee**2*KTRz3**2*c.MZ**2)/(2.*c.cw**2*c.sw**2))*cmath.sqrt(c.Mt**4 - 2*c.Mt**2*MtP**2 + MtP**4 - 2*c.Mt**2*c.MZ**2 - 2*MtP**2*c.MZ**2 + c.MZ**4))/(96.*cmath.pi*abs(MtP)**3)
    Gamma = decayTP_Ht + decayTP_Wb + decayTP_Zt
    return decayTP_Ht.real, decayTP_Wb.real, decayTP_Zt.real, Gamma.real 
#VLB
def decayBP4(MbP,Kappas,LR=False):
    KBLh3 = Kappas[0]
    KBLw3 = Kappas[1]
    KBLz3 = Kappas[2]
    if LR == True:
        KBRh3 = Kappas[0]
        KBRw3 = Kappas[1]
        KBRz3 = Kappas[2]
    else:
        KBRh3 = 0
        KBRw3 = 0
        KBRz3 = 0
    decayBP_Hb = ((3*KBLh3**2*c.Mb**2 + 3*KBRh3**2*c.Mb**2 + 12*KBLh3*KBRh3*c.Mb*MbP + 3*KBLh3**2*MbP**2 + 3*KBRh3**2*MbP**2 - 3*KBLh3**2*c.MH**2 - 3*KBRh3**2*c.MH**2)*cmath.sqrt(c.Mb**4 - 2*c.Mb**2*MbP**2 + MbP**4 - 2*c.Mb**2*c.MH**2 - 2*MbP**2*c.MH**2 + c.MH**4))/(96.*cmath.pi*abs(MbP)**3)
    decayBP_Wt = (((3*c.ee**2*KBLw3**2*MbP**2)/(2.*c.sw**2) + (3*c.ee**2*KBRw3**2*MbP**2)/(2.*c.sw**2) - (18*c.ee**2*KBLw3*KBRw3*MbP*c.Mt)/c.sw**2 + (3*c.ee**2*KBLw3**2*c.Mt**2)/(2.*c.sw**2) + (3*c.ee**2*KBRw3**2*c.Mt**2)/(2.*c.sw**2) + (3*c.ee**2*KBLw3**2*MbP**4)/(2.*c.MW**2*c.sw**2) + (3*c.ee**2*KBRw3**2*MbP**4)/(2.*c.MW**2*c.sw**2) - (3*c.ee**2*KBLw3**2*MbP**2*c.Mt**2)/(c.MW**2*c.sw**2) - (3*c.ee**2*KBRw3**2*MbP**2*c.Mt**2)/(c.MW**2*c.sw**2) + (3*c.ee**2*KBLw3**2*c.Mt**4)/(2.*c.MW**2*c.sw**2) + (3*c.ee**2*KBRw3**2*c.Mt**4)/(2.*c.MW**2*c.sw**2) - (3*c.ee**2*KBLw3**2*c.MW**2)/c.sw**2 - (3*c.ee**2*KBRw3**2*c.MW**2)/c.sw**2)*cmath.sqrt(MbP**4 - 2*MbP**2*c.Mt**2 + c.Mt**4 - 2*MbP**2*c.MW**2 - 2*c.Mt**2*c.MW**2 + c.MW**4))/(96.*cmath.pi*abs(MbP)**3)
    decayBP_Zb = (((3*c.ee**2*KBLz3**2*c.Mb**2)/(4.*c.cw**2*c.sw**2) + (3*c.ee**2*KBRz3**2*c.Mb**2)/(4.*c.cw**2*c.sw**2) - (9*c.ee**2*KBLz3*KBRz3*c.Mb*MbP)/(c.cw**2*c.sw**2) + (3*c.ee**2*KBLz3**2*MbP**2)/(4.*c.cw**2*c.sw**2) + (3*c.ee**2*KBRz3**2*MbP**2)/(4.*c.cw**2*c.sw**2) + (3*c.ee**2*KBLz3**2*c.Mb**4)/(4.*c.cw**2*c.MZ**2*c.sw**2) + (3*c.ee**2*KBRz3**2*c.Mb**4)/(4.*c.cw**2*c.MZ**2*c.sw**2) - (3*c.ee**2*KBLz3**2*c.Mb**2*MbP**2)/(2.*c.cw**2*c.MZ**2*c.sw**2) - (3*c.ee**2*KBRz3**2*c.Mb**2*MbP**2)/(2.*c.cw**2*c.MZ**2*c.sw**2) + (3*c.ee**2*KBLz3**2*MbP**4)/(4.*c.cw**2*c.MZ**2*c.sw**2) + (3*c.ee**2*KBRz3**2*MbP**4)/(4.*c.cw**2*c.MZ**2*c.sw**2) - (3*c.ee**2*KBLz3**2*c.MZ**2)/(2.*c.cw**2*c.sw**2) - (3*c.ee**2*KBRz3**2*c.MZ**2)/(2.*c.cw**2*c.sw**2))*cmath.sqrt(c.Mb**4 - 2*c.Mb**2*MbP**2 + MbP**4 - 2*c.Mb**2*c.MZ**2 - 2*MbP**2*c.MZ**2 + c.MZ**4))/(96.*cmath.pi*abs(MbP)**3)
    Gamma = decayBP_Hb + decayBP_Wt + decayBP_Zb
    return decayBP_Hb.real, decayBP_Wt.real, decayBP_Zb.real, Gamma.real
#VLX
def decayX4(MX,kappa,LR=False):
    KXL3 = kappa
    if LR == True:
        KXR3 = kappa
    else:
        KXR3 = 0
    return ((((3*c.ee**2*KXL3**2*c.Mt**2)/(2.*c.sw**2) + (3*c.ee**2*KXR3**2*c.Mt**2)/(2.*c.sw**2) + (3*c.ee**2*KXL3**2*c.Mt**4)/(2.*c.MW**2*c.sw**2) + (3*c.ee**2*KXR3**2*c.Mt**4)/(2.*c.MW**2*c.sw**2) - (3*c.ee**2*KXL3**2*c.MW**2)/c.sw**2 - (3*c.ee**2*KXR3**2*c.MW**2)/c.sw**2 - (18*c.ee**2*KXL3*KXR3*c.Mt*MX)/c.sw**2 + (3*c.ee**2*KXL3**2*MX**2)/(2.*c.sw**2) + (3*c.ee**2*KXR3**2*MX**2)/(2.*c.sw**2) - (3*c.ee**2*KXL3**2*c.Mt**2*MX**2)/(c.MW**2*c.sw**2) - (3*c.ee**2*KXR3**2*c.Mt**2*MX**2)/(c.MW**2*c.sw**2) + (3*c.ee**2*KXL3**2*MX**4)/(2.*c.MW**2*c.sw**2) + (3*c.ee**2*KXR3**2*MX**4)/(2.*c.MW**2*c.sw**2))*cmath.sqrt(c.Mt**4 - 2*c.Mt**2*c.MW**2 + c.MW**4 - 2*c.Mt**2*MX**2 - 2*c.MW**2*MX**2 + MX**4))/(96.*cmath.pi*abs(MX)**3)).real
#VLY
def decayY4(MY,kappa,LR=False):
    KYL3 = kappa
    if LR == True:
        KYR3 = kappa
    else:
        KYR3 = 0
    return ((((3*c.ee**2*KYL3**2*c.Mb**2)/(2.*c.sw**2) + (3*c.ee**2*KYR3**2*c.Mb**2)/(2.*c.sw**2) + (3*c.ee**2*KYL3**2*c.Mb**4)/(2.*c.MW**2*c.sw**2) + (3*c.ee**2*KYR3**2*c.Mb**4)/(2.*c.MW**2*c.sw**2) - (3*c.ee**2*KYL3**2*c.MW**2)/c.sw**2 - (3*c.ee**2*KYR3**2*c.MW**2)/c.sw**2 - (18*c.ee**2*KYL3*KYR3*c.Mb*MY)/c.sw**2 + (3*c.ee**2*KYL3**2*MY**2)/(2.*c.sw**2) + (3*c.ee**2*KYR3**2*MY**2)/(2.*c.sw**2) - (3*c.ee**2*KYL3**2*c.Mb**2*MY**2)/(c.MW**2*c.sw**2) - (3*c.ee**2*KYR3**2*c.Mb**2*MY**2)/(c.MW**2*c.sw**2) + (3*c.ee**2*KYL3**2*MY**4)/(2.*c.MW**2*c.sw**2) + (3*c.ee**2*KYR3**2*MY**4)/(2.*c.MW**2*c.sw**2))*cmath.sqrt(c.Mb**4 - 2*c.Mb**2*c.MW**2 + c.MW**4 - 2*c.Mb**2*MY**2 - 2*c.MW**2*MY**2 + MY**4))/(96.*cmath.pi*abs(MY)**3)).real
#######
# 5FNS
#######
#VLT
def decayTP5(MtP,Kappas,LR=False):
    KTLh3 = Kappas[0]
    KTLw3 = Kappas[1]
    KTLz3 = Kappas[2]
    if LR == True:
        KTRh3 = Kappas[0]
        KTRw3 = Kappas[1]
        KTRz3 = Kappas[2]
    else:
        KTRh3 = 0
        KTRw3 = 0
        KTRz3 = 0
    decayTP_Ht = ((-3*KTLh3**2*c.MH**2 - 3*KTRh3**2*c.MH**2 + 3*KTLh3**2*c.Mt**2 + 3*KTRh3**2*c.Mt**2 + 12*KTLh3*KTRh3*c.Mt*MtP + 3*KTLh3**2*MtP**2 + 3*KTRh3**2*MtP**2)*cmath.sqrt(c.MH**4 - 2*c.MH**2*c.Mt**2 + c.Mt**4 - 2*c.MH**2*MtP**2 - 2*c.Mt**2*MtP**2 + MtP**4))/(96.*cmath.pi*abs(MtP)**3)
    decayTP_Wb = ((MtP**2 - c.MW**2)*((3*c.ee**2*KTLw3**2*MtP**2)/(2.*c.sw**2) + (3*c.ee**2*KTRw3**2*MtP**2)/(2.*c.sw**2) + (3*c.ee**2*KTLw3**2*MtP**4)/(2.*c.MW**2*c.sw**2) + (3*c.ee**2*KTRw3**2*MtP**4)/(2.*c.MW**2*c.sw**2) - (3*c.ee**2*KTLw3**2*c.MW**2)/c.sw**2 - (3*c.ee**2*KTRw3**2*c.MW**2)/c.sw**2))/(96.*cmath.pi*abs(MtP)**3)      
    decayTP_Zt = (((3*c.ee**2*KTLz3**2*c.Mt**2)/(4.*c.cw**2*c.sw**2) + (3*c.ee**2*KTRz3**2*c.Mt**2)/(4.*c.cw**2*c.sw**2) - (9*c.ee**2*KTLz3*KTRz3*c.Mt*MtP)/(c.cw**2*c.sw**2) + (3*c.ee**2*KTLz3**2*MtP**2)/(4.*c.cw**2*c.sw**2) + (3*c.ee**2*KTRz3**2*MtP**2)/(4.*c.cw**2*c.sw**2) + (3*c.ee**2*KTLz3**2*c.Mt**4)/(4.*c.cw**2*c.MZ**2*c.sw**2) + (3*c.ee**2*KTRz3**2*c.Mt**4)/(4.*c.cw**2*c.MZ**2*c.sw**2) - (3*c.ee**2*KTLz3**2*c.Mt**2*MtP**2)/(2.*c.cw**2*c.MZ**2*c.sw**2) - (3*c.ee**2*KTRz3**2*c.Mt**2*MtP**2)/(2.*c.cw**2*c.MZ**2*c.sw**2) + (3*c.ee**2*KTLz3**2*MtP**4)/(4.*c.cw**2*c.MZ**2*c.sw**2) + (3*c.ee**2*KTRz3**2*MtP**4)/(4.*c.cw**2*c.MZ**2*c.sw**2) - (3*c.ee**2*KTLz3**2*c.MZ**2)/(2.*c.cw**2*c.sw**2) - (3*c.ee**2*KTRz3**2*c.MZ**2)/(2.*c.cw**2*c.sw**2))*cmath.sqrt(c.Mt**4 - 2*c.Mt**2*MtP**2 + MtP**4 - 2*c.Mt**2*c.MZ**2 - 2*MtP**2*c.MZ**2 + c.MZ**4))/(96.*cmath.pi*abs(MtP)**3)
    Gamma = decayTP_Ht + decayTP_Wb + decayTP_Zt
    return decayTP_Ht.real, decayTP_Wb.real, decayTP_Zt.real, Gamma.real 
#VLB
def decayBP5(MbP,Kappas,LR=False):
    KBLh3 = Kappas[0]
    KBLw3 = Kappas[1]
    KBLz3 = Kappas[2]
    if LR == True:
        KBRh3 = Kappas[0]
        KBRw3 = Kappas[1]
        KBRz3 = Kappas[2]
    else:
        KBRh3 = 0
        KBRw3 = 0
        KBRz3 = 0
    decayBP_Hb = ((MbP**2 - c.MH**2)*(3*KBLh3**2*MbP**2 + 3*KBRh3**2*MbP**2 - 3*KBLh3**2*c.MH**2 - 3*KBRh3**2*c.MH**2))/(96.*cmath.pi*abs(MbP)**3)
    decayBP_Wt = (((3*c.ee**2*KBLw3**2*MbP**2)/(2.*c.sw**2) + (3*c.ee**2*KBRw3**2*MbP**2)/(2.*c.sw**2) - (18*c.ee**2*KBLw3*KBRw3*MbP*c.Mt)/c.sw**2 + (3*c.ee**2*KBLw3**2*c.Mt**2)/(2.*c.sw**2) + (3*c.ee**2*KBRw3**2*c.Mt**2)/(2.*c.sw**2) + (3*c.ee**2*KBLw3**2*MbP**4)/(2.*c.MW**2*c.sw**2) + (3*c.ee**2*KBRw3**2*MbP**4)/(2.*c.MW**2*c.sw**2) - (3*c.ee**2*KBLw3**2*MbP**2*c.Mt**2)/(c.MW**2*c.sw**2) - (3*c.ee**2*KBRw3**2*MbP**2*c.Mt**2)/(c.MW**2*c.sw**2) + (3*c.ee**2*KBLw3**2*c.Mt**4)/(2.*c.MW**2*c.sw**2) + (3*c.ee**2*KBRw3**2*c.Mt**4)/(2.*c.MW**2*c.sw**2) - (3*c.ee**2*KBLw3**2*c.MW**2)/c.sw**2 - (3*c.ee**2*KBRw3**2*c.MW**2)/c.sw**2)*cmath.sqrt(MbP**4 - 2*MbP**2*c.Mt**2 + c.Mt**4 - 2*MbP**2*c.MW**2 - 2*c.Mt**2*c.MW**2 + c.MW**4))/(96.*cmath.pi*abs(MbP)**3)
    decayBP_Zb = ((MbP**2 - c.MZ**2)*((3*c.ee**2*KBLz3**2*MbP**2)/(4.*c.cw**2*c.sw**2) + (3*c.ee**2*KBRz3**2*MbP**2)/(4.*c.cw**2*c.sw**2) + (3*c.ee**2*KBLz3**2*MbP**4)/(4.*c.cw**2*c.MZ**2*c.sw**2) + (3*c.ee**2*KBRz3**2*MbP**4)/(4.*c.cw**2*c.MZ**2*c.sw**2) - (3*c.ee**2*KBLz3**2*c.MZ**2)/(2.*c.cw**2*c.sw**2) - (3*c.ee**2*KBRz3**2*c.MZ**2)/(2.*c.cw**2*c.sw**2)))/(96.*cmath.pi*abs(MbP)**3)
    Gamma = decayBP_Hb + decayBP_Wt + decayBP_Zb
    return decayBP_Hb.real, decayBP_Wt.real, decayBP_Zb.real, Gamma.real
#VLX
def decayX5(MX,Kappa,LR=False):
    KXL3 = Kappa
    if LR == True:
        KXR3 = Kappa
    else:
        KXR3 = 0
    return ((((3*c.ee**2*KXL3**2*c.Mt**2)/(2.*c.sw**2) + (3*c.ee**2*KXR3**2*c.Mt**2)/(2.*c.sw**2) + (3*c.ee**2*KXL3**2*c.Mt**4)/(2.*c.MW**2*c.sw**2) + (3*c.ee**2*KXR3**2*c.Mt**4)/(2.*c.MW**2*c.sw**2) - (3*c.ee**2*KXL3**2*c.MW**2)/c.sw**2 - (3*c.ee**2*KXR3**2*c.MW**2)/c.sw**2 - (18*c.ee**2*KXL3*KXR3*c.Mt*MX)/c.sw**2 + (3*c.ee**2*KXL3**2*MX**2)/(2.*c.sw**2) + (3*c.ee**2*KXR3**2*MX**2)/(2.*c.sw**2) - (3*c.ee**2*KXL3**2*c.Mt**2*MX**2)/(c.MW**2*c.sw**2) - (3*c.ee**2*KXR3**2*c.Mt**2*MX**2)/(c.MW**2*c.sw**2) + (3*c.ee**2*KXL3**2*MX**4)/(2.*c.MW**2*c.sw**2) + (3*c.ee**2*KXR3**2*MX**4)/(2.*c.MW**2*c.sw**2))*cmath.sqrt(c.Mt**4 - 2*c.Mt**2*c.MW**2 + c.MW**4 - 2*c.Mt**2*MX**2 - 2*c.MW**2*MX**2 + MX**4))/(96.*cmath.pi*abs(MX)**3)).real
#VLY
def decayY5(MY,Kappa,LR=False):
    KYL3 = Kappa
    if LR == True:
        KYR3 = Kappa
    else:
        KYR3 = 0
    return (((-c.MW**2 + MY**2)*((-3*c.ee**2*KYL3**2*c.MW**2)/c.sw**2 - (3*c.ee**2*KYR3**2*c.MW**2)/c.sw**2 + (3*c.ee**2*KYL3**2*MY**2)/(2.*c.sw**2) + (3*c.ee**2*KYR3**2*MY**2)/(2.*c.sw**2) + (3*c.ee**2*KYL3**2*MY**4)/(2.*c.MW**2*c.sw**2) + (3*c.ee**2*KYR3**2*MY**4)/(2.*c.MW**2*c.sw**2)))/(96.*cmath.pi*abs(MY)**3)).real
#######
# All generations
#######
# Decay channels: boson -> SM quarks of generation 1, 2, 3
UP, DOWN = ['u','c','t'], ['d','s','b']
CHANNELS = {'VLT': {'h':UP, 'w':DOWN, 'z':UP},
            'VLB': {'h':DOWN, 'w':UP, 'z':DOWN},
            'VLX': {'w':UP},
            'VLY': {'w':DOWN}}

def quarkMass(quark,FNS):
    if quark == 't':
        return c.Mt
    if quark == 'b' and FNS == 4:
        return c.Mb
    return 0.

def sqrtLambda(M,m1,m2):
    return cmath.sqrt(M**4 + m1**4 + m2**4 - 2*M**2*m1**2 - 2*M**2*m2**2 - 2*m1**2*m2**2)

def widthH(M,mq,kL,kR):
    # Q -> H q, same expression as the 3rd generation formulas above with a general quark mass
    return ((3*(kL**2 + kR**2)*(mq**2 + M**2 - c.MH**2) + 12*kL*kR*mq*M)*sqrtLambda(M,mq,c.MH)/(96.*cmath.pi*abs(M)**3)).real

def widthV(M,mq,MV,g2,kL,kR):
    # Q -> V q (V = W, Z); g2 is the squared gauge factor of the vertex
    return ((3*g2*(kL**2 + kR**2)*(mq**2 + M**2 + (M**2 - mq**2)**2/MV**2 - 2*MV**2) - 36*g2*kL*kR*mq*M)*sqrtLambda(M,mq,MV)/(96.*cmath.pi*abs(M)**3)).real

def partialWidths(VLQtype,FNS,M,Kappas,LR=False):
    '''Partial widths {(boson, generation): width} in GeV for couplings {(boson, generation): kappa}'''
    widths = {}
    for boson, quarks in CHANNELS[VLQtype].items():
        for gen, quark in enumerate(quarks, start=1):
            kL = Kappas.get((boson,gen), 0)
            kR = kL if LR else 0
            mq = quarkMass(quark,FNS)
            if boson == 'h':
                widths[(boson,gen)] = widthH(M,mq,kL,kR)
            elif boson == 'w':
                widths[(boson,gen)] = widthV(M,mq,c.MW,(c.ee**2/(2.*c.sw**2)).real,kL,kR)
            else:
                widths[(boson,gen)] = widthV(M,mq,c.MZ,(c.ee**2/(4.*c.cw**2*c.sw**2)).real,kL,kR)
    return widths
