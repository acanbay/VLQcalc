'''
Author: Ali Can Canbay - acanbay
Date  : 21 May 2024

Calculates decay widths for VLQ_v4_4FNS_UFO and VLQ_v4_5FNS_UFO models [1,2].
All constants and formulas are taken directly from the models.

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
############
MH=125.10 
MW=80.379 
MZ=91.1876 
Mt=172.9 
Mb=4.18
aEW=1./127.9
ee=2*cmath.sqrt(aEW)*cmath.sqrt(cmath.pi)
sw2=1 - MW**2/MZ**2
sw=cmath.sqrt(sw2)
cw=cmath.sqrt(1 - sw2)
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
    decayTP_Ht = ((-3*KTLh3**2*MH**2 - 3*KTRh3**2*MH**2 + 3*KTLh3**2*Mt**2 + 3*KTRh3**2*Mt**2 + 12*KTLh3*KTRh3*Mt*MtP + 3*KTLh3**2*MtP**2 + 3*KTRh3**2*MtP**2)*cmath.sqrt(MH**4 - 2*MH**2*Mt**2 + Mt**4 - 2*MH**2*MtP**2 - 2*Mt**2*MtP**2 + MtP**4))/(96.*cmath.pi*abs(MtP)**3)
    decayTP_Wb = (((3*ee**2*KTLw3**2*Mb**2)/(2.*sw**2) + (3*ee**2*KTRw3**2*Mb**2)/(2.*sw**2) - (18*ee**2*KTLw3*KTRw3*Mb*MtP)/sw**2 + (3*ee**2*KTLw3**2*MtP**2)/(2.*sw**2) + (3*ee**2*KTRw3**2*MtP**2)/(2.*sw**2) + (3*ee**2*KTLw3**2*Mb**4)/(2.*MW**2*sw**2) + (3*ee**2*KTRw3**2*Mb**4)/(2.*MW**2*sw**2) - (3*ee**2*KTLw3**2*Mb**2*MtP**2)/(MW**2*sw**2) - (3*ee**2*KTRw3**2*Mb**2*MtP**2)/(MW**2*sw**2) + (3*ee**2*KTLw3**2*MtP**4)/(2.*MW**2*sw**2) + (3*ee**2*KTRw3**2*MtP**4)/(2.*MW**2*sw**2) - (3*ee**2*KTLw3**2*MW**2)/sw**2 - (3*ee**2*KTRw3**2*MW**2)/sw**2)*cmath.sqrt(Mb**4 - 2*Mb**2*MtP**2 + MtP**4 - 2*Mb**2*MW**2 - 2*MtP**2*MW**2 + MW**4))/(96.*cmath.pi*abs(MtP)**3)
    decayTP_Zt = (((3*ee**2*KTLz3**2*Mt**2)/(4.*cw**2*sw**2) + (3*ee**2*KTRz3**2*Mt**2)/(4.*cw**2*sw**2) - (9*ee**2*KTLz3*KTRz3*Mt*MtP)/(cw**2*sw**2) + (3*ee**2*KTLz3**2*MtP**2)/(4.*cw**2*sw**2) + (3*ee**2*KTRz3**2*MtP**2)/(4.*cw**2*sw**2) + (3*ee**2*KTLz3**2*Mt**4)/(4.*cw**2*MZ**2*sw**2) + (3*ee**2*KTRz3**2*Mt**4)/(4.*cw**2*MZ**2*sw**2) - (3*ee**2*KTLz3**2*Mt**2*MtP**2)/(2.*cw**2*MZ**2*sw**2) - (3*ee**2*KTRz3**2*Mt**2*MtP**2)/(2.*cw**2*MZ**2*sw**2) + (3*ee**2*KTLz3**2*MtP**4)/(4.*cw**2*MZ**2*sw**2) + (3*ee**2*KTRz3**2*MtP**4)/(4.*cw**2*MZ**2*sw**2) - (3*ee**2*KTLz3**2*MZ**2)/(2.*cw**2*sw**2) - (3*ee**2*KTRz3**2*MZ**2)/(2.*cw**2*sw**2))*cmath.sqrt(Mt**4 - 2*Mt**2*MtP**2 + MtP**4 - 2*Mt**2*MZ**2 - 2*MtP**2*MZ**2 + MZ**4))/(96.*cmath.pi*abs(MtP)**3)
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
    decayBP_Hb = ((3*KBLh3**2*Mb**2 + 3*KBRh3**2*Mb**2 + 12*KBLh3*KBRh3*Mb*MbP + 3*KBLh3**2*MbP**2 + 3*KBRh3**2*MbP**2 - 3*KBLh3**2*MH**2 - 3*KBRh3**2*MH**2)*cmath.sqrt(Mb**4 - 2*Mb**2*MbP**2 + MbP**4 - 2*Mb**2*MH**2 - 2*MbP**2*MH**2 + MH**4))/(96.*cmath.pi*abs(MbP)**3)
    decayBP_Wt = (((3*ee**2*KBLw3**2*MbP**2)/(2.*sw**2) + (3*ee**2*KBRw3**2*MbP**2)/(2.*sw**2) - (18*ee**2*KBLw3*KBRw3*MbP*Mt)/sw**2 + (3*ee**2*KBLw3**2*Mt**2)/(2.*sw**2) + (3*ee**2*KBRw3**2*Mt**2)/(2.*sw**2) + (3*ee**2*KBLw3**2*MbP**4)/(2.*MW**2*sw**2) + (3*ee**2*KBRw3**2*MbP**4)/(2.*MW**2*sw**2) - (3*ee**2*KBLw3**2*MbP**2*Mt**2)/(MW**2*sw**2) - (3*ee**2*KBRw3**2*MbP**2*Mt**2)/(MW**2*sw**2) + (3*ee**2*KBLw3**2*Mt**4)/(2.*MW**2*sw**2) + (3*ee**2*KBRw3**2*Mt**4)/(2.*MW**2*sw**2) - (3*ee**2*KBLw3**2*MW**2)/sw**2 - (3*ee**2*KBRw3**2*MW**2)/sw**2)*cmath.sqrt(MbP**4 - 2*MbP**2*Mt**2 + Mt**4 - 2*MbP**2*MW**2 - 2*Mt**2*MW**2 + MW**4))/(96.*cmath.pi*abs(MbP)**3)
    decayBP_Zb = (((3*ee**2*KBLz3**2*Mb**2)/(4.*cw**2*sw**2) + (3*ee**2*KBRz3**2*Mb**2)/(4.*cw**2*sw**2) - (9*ee**2*KBLz3*KBRz3*Mb*MbP)/(cw**2*sw**2) + (3*ee**2*KBLz3**2*MbP**2)/(4.*cw**2*sw**2) + (3*ee**2*KBRz3**2*MbP**2)/(4.*cw**2*sw**2) + (3*ee**2*KBLz3**2*Mb**4)/(4.*cw**2*MZ**2*sw**2) + (3*ee**2*KBRz3**2*Mb**4)/(4.*cw**2*MZ**2*sw**2) - (3*ee**2*KBLz3**2*Mb**2*MbP**2)/(2.*cw**2*MZ**2*sw**2) - (3*ee**2*KBRz3**2*Mb**2*MbP**2)/(2.*cw**2*MZ**2*sw**2) + (3*ee**2*KBLz3**2*MbP**4)/(4.*cw**2*MZ**2*sw**2) + (3*ee**2*KBRz3**2*MbP**4)/(4.*cw**2*MZ**2*sw**2) - (3*ee**2*KBLz3**2*MZ**2)/(2.*cw**2*sw**2) - (3*ee**2*KBRz3**2*MZ**2)/(2.*cw**2*sw**2))*cmath.sqrt(Mb**4 - 2*Mb**2*MbP**2 + MbP**4 - 2*Mb**2*MZ**2 - 2*MbP**2*MZ**2 + MZ**4))/(96.*cmath.pi*abs(MbP)**3)
    Gamma = decayBP_Hb + decayBP_Wt + decayBP_Zb
    return decayBP_Hb.real, decayBP_Wt.real, decayBP_Zb.real, Gamma.real
#VLX
def decayX4(MX,kappa,LR=False):
    KXL3 = kappa
    if LR == True:
        KXR3 = kappa
    else:
        KXR3 = 0
    return ((((3*ee**2*KXL3**2*Mt**2)/(2.*sw**2) + (3*ee**2*KXR3**2*Mt**2)/(2.*sw**2) + (3*ee**2*KXL3**2*Mt**4)/(2.*MW**2*sw**2) + (3*ee**2*KXR3**2*Mt**4)/(2.*MW**2*sw**2) - (3*ee**2*KXL3**2*MW**2)/sw**2 - (3*ee**2*KXR3**2*MW**2)/sw**2 - (18*ee**2*KXL3*KXR3*Mt*MX)/sw**2 + (3*ee**2*KXL3**2*MX**2)/(2.*sw**2) + (3*ee**2*KXR3**2*MX**2)/(2.*sw**2) - (3*ee**2*KXL3**2*Mt**2*MX**2)/(MW**2*sw**2) - (3*ee**2*KXR3**2*Mt**2*MX**2)/(MW**2*sw**2) + (3*ee**2*KXL3**2*MX**4)/(2.*MW**2*sw**2) + (3*ee**2*KXR3**2*MX**4)/(2.*MW**2*sw**2))*cmath.sqrt(Mt**4 - 2*Mt**2*MW**2 + MW**4 - 2*Mt**2*MX**2 - 2*MW**2*MX**2 + MX**4))/(96.*cmath.pi*abs(MX)**3)).real
#VLY
def decayY4(MY,kappa,LR=False):
    KYL3 = kappa
    if LR == True:
        KYR3 = kappa
    else:
        KYR3 = 0
    return ((((3*ee**2*KYL3**2*Mb**2)/(2.*sw**2) + (3*ee**2*KYR3**2*Mb**2)/(2.*sw**2) + (3*ee**2*KYL3**2*Mb**4)/(2.*MW**2*sw**2) + (3*ee**2*KYR3**2*Mb**4)/(2.*MW**2*sw**2) - (3*ee**2*KYL3**2*MW**2)/sw**2 - (3*ee**2*KYR3**2*MW**2)/sw**2 - (18*ee**2*KYL3*KYR3*Mb*MY)/sw**2 + (3*ee**2*KYL3**2*MY**2)/(2.*sw**2) + (3*ee**2*KYR3**2*MY**2)/(2.*sw**2) - (3*ee**2*KYL3**2*Mb**2*MY**2)/(MW**2*sw**2) - (3*ee**2*KYR3**2*Mb**2*MY**2)/(MW**2*sw**2) + (3*ee**2*KYL3**2*MY**4)/(2.*MW**2*sw**2) + (3*ee**2*KYR3**2*MY**4)/(2.*MW**2*sw**2))*cmath.sqrt(Mb**4 - 2*Mb**2*MW**2 + MW**4 - 2*Mb**2*MY**2 - 2*MW**2*MY**2 + MY**4))/(96.*cmath.pi*abs(MY)**3)).real
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
    decayTP_Ht = ((-3*KTLh3**2*MH**2 - 3*KTRh3**2*MH**2 + 3*KTLh3**2*Mt**2 + 3*KTRh3**2*Mt**2 + 12*KTLh3*KTRh3*Mt*MtP + 3*KTLh3**2*MtP**2 + 3*KTRh3**2*MtP**2)*cmath.sqrt(MH**4 - 2*MH**2*Mt**2 + Mt**4 - 2*MH**2*MtP**2 - 2*Mt**2*MtP**2 + MtP**4))/(96.*cmath.pi*abs(MtP)**3)
    decayTP_Wb = ((MtP**2 - MW**2)*((3*ee**2*KTLw3**2*MtP**2)/(2.*sw**2) + (3*ee**2*KTRw3**2*MtP**2)/(2.*sw**2) + (3*ee**2*KTLw3**2*MtP**4)/(2.*MW**2*sw**2) + (3*ee**2*KTRw3**2*MtP**4)/(2.*MW**2*sw**2) - (3*ee**2*KTLw3**2*MW**2)/sw**2 - (3*ee**2*KTRw3**2*MW**2)/sw**2))/(96.*cmath.pi*abs(MtP)**3)      
    decayTP_Zt = (((3*ee**2*KTLz3**2*Mt**2)/(4.*cw**2*sw**2) + (3*ee**2*KTRz3**2*Mt**2)/(4.*cw**2*sw**2) - (9*ee**2*KTLz3*KTRz3*Mt*MtP)/(cw**2*sw**2) + (3*ee**2*KTLz3**2*MtP**2)/(4.*cw**2*sw**2) + (3*ee**2*KTRz3**2*MtP**2)/(4.*cw**2*sw**2) + (3*ee**2*KTLz3**2*Mt**4)/(4.*cw**2*MZ**2*sw**2) + (3*ee**2*KTRz3**2*Mt**4)/(4.*cw**2*MZ**2*sw**2) - (3*ee**2*KTLz3**2*Mt**2*MtP**2)/(2.*cw**2*MZ**2*sw**2) - (3*ee**2*KTRz3**2*Mt**2*MtP**2)/(2.*cw**2*MZ**2*sw**2) + (3*ee**2*KTLz3**2*MtP**4)/(4.*cw**2*MZ**2*sw**2) + (3*ee**2*KTRz3**2*MtP**4)/(4.*cw**2*MZ**2*sw**2) - (3*ee**2*KTLz3**2*MZ**2)/(2.*cw**2*sw**2) - (3*ee**2*KTRz3**2*MZ**2)/(2.*cw**2*sw**2))*cmath.sqrt(Mt**4 - 2*Mt**2*MtP**2 + MtP**4 - 2*Mt**2*MZ**2 - 2*MtP**2*MZ**2 + MZ**4))/(96.*cmath.pi*abs(MtP)**3)
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
    decayBP_Hb = ((MbP**2 - MH**2)*(3*KBLh3**2*MbP**2 + 3*KBRh3**2*MbP**2 - 3*KBLh3**2*MH**2 - 3*KBRh3**2*MH**2))/(96.*cmath.pi*abs(MbP)**3)
    decayBP_Wt = (((3*ee**2*KBLw3**2*MbP**2)/(2.*sw**2) + (3*ee**2*KBRw3**2*MbP**2)/(2.*sw**2) - (18*ee**2*KBLw3*KBRw3*MbP*Mt)/sw**2 + (3*ee**2*KBLw3**2*Mt**2)/(2.*sw**2) + (3*ee**2*KBRw3**2*Mt**2)/(2.*sw**2) + (3*ee**2*KBLw3**2*MbP**4)/(2.*MW**2*sw**2) + (3*ee**2*KBRw3**2*MbP**4)/(2.*MW**2*sw**2) - (3*ee**2*KBLw3**2*MbP**2*Mt**2)/(MW**2*sw**2) - (3*ee**2*KBRw3**2*MbP**2*Mt**2)/(MW**2*sw**2) + (3*ee**2*KBLw3**2*Mt**4)/(2.*MW**2*sw**2) + (3*ee**2*KBRw3**2*Mt**4)/(2.*MW**2*sw**2) - (3*ee**2*KBLw3**2*MW**2)/sw**2 - (3*ee**2*KBRw3**2*MW**2)/sw**2)*cmath.sqrt(MbP**4 - 2*MbP**2*Mt**2 + Mt**4 - 2*MbP**2*MW**2 - 2*Mt**2*MW**2 + MW**4))/(96.*cmath.pi*abs(MbP)**3)
    decayBP_Zb = ((MbP**2 - MZ**2)*((3*ee**2*KBLz3**2*MbP**2)/(4.*cw**2*sw**2) + (3*ee**2*KBRz3**2*MbP**2)/(4.*cw**2*sw**2) + (3*ee**2*KBLz3**2*MbP**4)/(4.*cw**2*MZ**2*sw**2) + (3*ee**2*KBRz3**2*MbP**4)/(4.*cw**2*MZ**2*sw**2) - (3*ee**2*KBLz3**2*MZ**2)/(2.*cw**2*sw**2) - (3*ee**2*KBRz3**2*MZ**2)/(2.*cw**2*sw**2)))/(96.*cmath.pi*abs(MbP)**3)
    Gamma = decayBP_Hb + decayBP_Wt + decayBP_Zb
    return decayBP_Hb.real, decayBP_Wt.real, decayBP_Zb.real, Gamma.real
#VLX
def decayX5(MX,Kappa,LR=False):
    KXL3 = Kappa
    if LR == True:
        KXR3 = Kappa
    else:
        KXR3 = 0
    return ((((3*ee**2*KXL3**2*Mt**2)/(2.*sw**2) + (3*ee**2*KXR3**2*Mt**2)/(2.*sw**2) + (3*ee**2*KXL3**2*Mt**4)/(2.*MW**2*sw**2) + (3*ee**2*KXR3**2*Mt**4)/(2.*MW**2*sw**2) - (3*ee**2*KXL3**2*MW**2)/sw**2 - (3*ee**2*KXR3**2*MW**2)/sw**2 - (18*ee**2*KXL3*KXR3*Mt*MX)/sw**2 + (3*ee**2*KXL3**2*MX**2)/(2.*sw**2) + (3*ee**2*KXR3**2*MX**2)/(2.*sw**2) - (3*ee**2*KXL3**2*Mt**2*MX**2)/(MW**2*sw**2) - (3*ee**2*KXR3**2*Mt**2*MX**2)/(MW**2*sw**2) + (3*ee**2*KXL3**2*MX**4)/(2.*MW**2*sw**2) + (3*ee**2*KXR3**2*MX**4)/(2.*MW**2*sw**2))*cmath.sqrt(Mt**4 - 2*Mt**2*MW**2 + MW**4 - 2*Mt**2*MX**2 - 2*MW**2*MX**2 + MX**4))/(96.*cmath.pi*abs(MX)**3)).real
#VLY
def decayY5(MY,Kappa,LR=False):
    KYL3 = Kappa
    if LR == True:
        KYR3 = Kappa
    else:
        KYR3 = 0
    return (((-MW**2 + MY**2)*((-3*ee**2*KYL3**2*MW**2)/sw**2 - (3*ee**2*KYR3**2*MW**2)/sw**2 + (3*ee**2*KYL3**2*MY**2)/(2.*sw**2) + (3*ee**2*KYR3**2*MY**2)/(2.*sw**2) + (3*ee**2*KYL3**2*MY**4)/(2.*MW**2*sw**2) + (3*ee**2*KYR3**2*MY**4)/(2.*MW**2*sw**2)))/(96.*cmath.pi*abs(MY)**3)).real
