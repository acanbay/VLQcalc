"""
Author: Ali Can Canbay - acanbay
Date  : 21 May 2024

Calculates decay widths for VLQ_v4_4FNS_UFO and VLQ_v4_4FNS_UFO models [1,2].
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
"""
import cmath
############
MB=4.18
MH=125.10 
MW=80.379 
MZ=91.1876 
MT=172.9 
aEW=1./127.9
ee=2*cmath.sqrt(aEW)*cmath.sqrt(cmath.pi)
sw2=1 - MW**2/MZ**2
sw=cmath.sqrt(sw2)
cw=cmath.sqrt(1 - sw2)
#######
# 4FNS
#######
#VLT
def decayTP4(MTP,kappas,LR=None):
    KTLh3 = kappas[0]
    KTLw3 = kappas[1]
    KTLz3 = kappas[2]
    if LR == True:
        KTRh3 = kappas[0]
        KTRw3 = kappas[1]
        KTRz3 = kappas[2]
    else:
        KTRh3 = 0
        KTRw3 = 0
        KTRz3 = 0
    decayTP_Ht = ((-3*KTLh3**2*MH**2 - 3*KTRh3**2*MH**2 + 3*KTLh3**2*MT**2 + 3*KTRh3**2*MT**2 + 12*KTLh3*KTRh3*MT*MTP + 3*KTLh3**2*MTP**2 + 3*KTRh3**2*MTP**2)*cmath.sqrt(MH**4 - 2*MH**2*MT**2 + MT**4 - 2*MH**2*MTP**2 - 2*MT**2*MTP**2 + MTP**4))/(96.*cmath.pi*abs(MTP)**3)
    decayTP_Wb = (((3*ee**2*KTLw3**2*MB**2)/(2.*sw**2) + (3*ee**2*KTRw3**2*MB**2)/(2.*sw**2) - (18*ee**2*KTLw3*KTRw3*MB*MTP)/sw**2 + (3*ee**2*KTLw3**2*MTP**2)/(2.*sw**2) + (3*ee**2*KTRw3**2*MTP**2)/(2.*sw**2) + (3*ee**2*KTLw3**2*MB**4)/(2.*MW**2*sw**2) + (3*ee**2*KTRw3**2*MB**4)/(2.*MW**2*sw**2) - (3*ee**2*KTLw3**2*MB**2*MTP**2)/(MW**2*sw**2) - (3*ee**2*KTRw3**2*MB**2*MTP**2)/(MW**2*sw**2) + (3*ee**2*KTLw3**2*MTP**4)/(2.*MW**2*sw**2) + (3*ee**2*KTRw3**2*MTP**4)/(2.*MW**2*sw**2) - (3*ee**2*KTLw3**2*MW**2)/sw**2 - (3*ee**2*KTRw3**2*MW**2)/sw**2)*cmath.sqrt(MB**4 - 2*MB**2*MTP**2 + MTP**4 - 2*MB**2*MW**2 - 2*MTP**2*MW**2 + MW**4))/(96.*cmath.pi*abs(MTP)**3)
    decayTP_Zt = (((3*ee**2*KTLz3**2*MT**2)/(4.*cw**2*sw**2) + (3*ee**2*KTRz3**2*MT**2)/(4.*cw**2*sw**2) - (9*ee**2*KTLz3*KTRz3*MT*MTP)/(cw**2*sw**2) + (3*ee**2*KTLz3**2*MTP**2)/(4.*cw**2*sw**2) + (3*ee**2*KTRz3**2*MTP**2)/(4.*cw**2*sw**2) + (3*ee**2*KTLz3**2*MT**4)/(4.*cw**2*MZ**2*sw**2) + (3*ee**2*KTRz3**2*MT**4)/(4.*cw**2*MZ**2*sw**2) - (3*ee**2*KTLz3**2*MT**2*MTP**2)/(2.*cw**2*MZ**2*sw**2) - (3*ee**2*KTRz3**2*MT**2*MTP**2)/(2.*cw**2*MZ**2*sw**2) + (3*ee**2*KTLz3**2*MTP**4)/(4.*cw**2*MZ**2*sw**2) + (3*ee**2*KTRz3**2*MTP**4)/(4.*cw**2*MZ**2*sw**2) - (3*ee**2*KTLz3**2*MZ**2)/(2.*cw**2*sw**2) - (3*ee**2*KTRz3**2*MZ**2)/(2.*cw**2*sw**2))*cmath.sqrt(MT**4 - 2*MT**2*MTP**2 + MTP**4 - 2*MT**2*MZ**2 - 2*MTP**2*MZ**2 + MZ**4))/(96.*cmath.pi*abs(MTP)**3)
    Gamma = decayTP_Ht + decayTP_Wb + decayTP_Zt
    return decayTP_Ht.real, decayTP_Wb.real, decayTP_Zt.real, Gamma.real 
#VLB
def decayBP4(MBP,kappas,LR=None):
    KBLh3 = kappas[0]
    KBLw3 = kappas[1]
    KBLz3 = kappas[2]
    if LR == True:
        KBRh3 = kappas[0]
        KBRw3 = kappas[1]
        KBRz3 = kappas[2]
    else:
        KBRh3 = 0
        KBRw3 = 0
        KBRz3 = 0
    decayBP_Hb = ((3*KBLh3**2*MB**2 + 3*KBRh3**2*MB**2 + 12*KBLh3*KBRh3*MB*MBP + 3*KBLh3**2*MBP**2 + 3*KBRh3**2*MBP**2 - 3*KBLh3**2*MH**2 - 3*KBRh3**2*MH**2)*cmath.sqrt(MB**4 - 2*MB**2*MBP**2 + MBP**4 - 2*MB**2*MH**2 - 2*MBP**2*MH**2 + MH**4))/(96.*cmath.pi*abs(MBP)**3)
    decayBP_Wt = (((3*ee**2*KBLw3**2*MBP**2)/(2.*sw**2) + (3*ee**2*KBRw3**2*MBP**2)/(2.*sw**2) - (18*ee**2*KBLw3*KBRw3*MBP*MT)/sw**2 + (3*ee**2*KBLw3**2*MT**2)/(2.*sw**2) + (3*ee**2*KBRw3**2*MT**2)/(2.*sw**2) + (3*ee**2*KBLw3**2*MBP**4)/(2.*MW**2*sw**2) + (3*ee**2*KBRw3**2*MBP**4)/(2.*MW**2*sw**2) - (3*ee**2*KBLw3**2*MBP**2*MT**2)/(MW**2*sw**2) - (3*ee**2*KBRw3**2*MBP**2*MT**2)/(MW**2*sw**2) + (3*ee**2*KBLw3**2*MT**4)/(2.*MW**2*sw**2) + (3*ee**2*KBRw3**2*MT**4)/(2.*MW**2*sw**2) - (3*ee**2*KBLw3**2*MW**2)/sw**2 - (3*ee**2*KBRw3**2*MW**2)/sw**2)*cmath.sqrt(MBP**4 - 2*MBP**2*MT**2 + MT**4 - 2*MBP**2*MW**2 - 2*MT**2*MW**2 + MW**4))/(96.*cmath.pi*abs(MBP)**3)
    decayBP_Zb = (((3*ee**2*KBLz3**2*MB**2)/(4.*cw**2*sw**2) + (3*ee**2*KBRz3**2*MB**2)/(4.*cw**2*sw**2) - (9*ee**2*KBLz3*KBRz3*MB*MBP)/(cw**2*sw**2) + (3*ee**2*KBLz3**2*MBP**2)/(4.*cw**2*sw**2) + (3*ee**2*KBRz3**2*MBP**2)/(4.*cw**2*sw**2) + (3*ee**2*KBLz3**2*MB**4)/(4.*cw**2*MZ**2*sw**2) + (3*ee**2*KBRz3**2*MB**4)/(4.*cw**2*MZ**2*sw**2) - (3*ee**2*KBLz3**2*MB**2*MBP**2)/(2.*cw**2*MZ**2*sw**2) - (3*ee**2*KBRz3**2*MB**2*MBP**2)/(2.*cw**2*MZ**2*sw**2) + (3*ee**2*KBLz3**2*MBP**4)/(4.*cw**2*MZ**2*sw**2) + (3*ee**2*KBRz3**2*MBP**4)/(4.*cw**2*MZ**2*sw**2) - (3*ee**2*KBLz3**2*MZ**2)/(2.*cw**2*sw**2) - (3*ee**2*KBRz3**2*MZ**2)/(2.*cw**2*sw**2))*cmath.sqrt(MB**4 - 2*MB**2*MBP**2 + MBP**4 - 2*MB**2*MZ**2 - 2*MBP**2*MZ**2 + MZ**4))/(96.*cmath.pi*abs(MBP)**3)
    Gamma = decayBP_Hb + decayBP_Wt + decayBP_Zb
    return decayBP_Hb.real, decayBP_Wt.real, decayBP_Zb.real, Gamma.real
#VLX
def decayX4(MX,kappa,LR=None):
    KXL3 = kappa
    if LR == True:
        KXR3 = kappa
    else:
        KXR3 = 0
    return ((((3*ee**2*KXL3**2*MT**2)/(2.*sw**2) + (3*ee**2*KXR3**2*MT**2)/(2.*sw**2) + (3*ee**2*KXL3**2*MT**4)/(2.*MW**2*sw**2) + (3*ee**2*KXR3**2*MT**4)/(2.*MW**2*sw**2) - (3*ee**2*KXL3**2*MW**2)/sw**2 - (3*ee**2*KXR3**2*MW**2)/sw**2 - (18*ee**2*KXL3*KXR3*MT*MX)/sw**2 + (3*ee**2*KXL3**2*MX**2)/(2.*sw**2) + (3*ee**2*KXR3**2*MX**2)/(2.*sw**2) - (3*ee**2*KXL3**2*MT**2*MX**2)/(MW**2*sw**2) - (3*ee**2*KXR3**2*MT**2*MX**2)/(MW**2*sw**2) + (3*ee**2*KXL3**2*MX**4)/(2.*MW**2*sw**2) + (3*ee**2*KXR3**2*MX**4)/(2.*MW**2*sw**2))*cmath.sqrt(MT**4 - 2*MT**2*MW**2 + MW**4 - 2*MT**2*MX**2 - 2*MW**2*MX**2 + MX**4))/(96.*cmath.pi*abs(MX)**3)).real
#VLY
def decayY4(MY,kappa,LR=None):
    KYL3 = kappa
    if LR == True:
        KYR3 = kappa
    else:
        KYR3 = 0
    return ((((3*ee**2*KYL3**2*MB**2)/(2.*sw**2) + (3*ee**2*KYR3**2*MB**2)/(2.*sw**2) + (3*ee**2*KYL3**2*MB**4)/(2.*MW**2*sw**2) + (3*ee**2*KYR3**2*MB**4)/(2.*MW**2*sw**2) - (3*ee**2*KYL3**2*MW**2)/sw**2 - (3*ee**2*KYR3**2*MW**2)/sw**2 - (18*ee**2*KYL3*KYR3*MB*MY)/sw**2 + (3*ee**2*KYL3**2*MY**2)/(2.*sw**2) + (3*ee**2*KYR3**2*MY**2)/(2.*sw**2) - (3*ee**2*KYL3**2*MB**2*MY**2)/(MW**2*sw**2) - (3*ee**2*KYR3**2*MB**2*MY**2)/(MW**2*sw**2) + (3*ee**2*KYL3**2*MY**4)/(2.*MW**2*sw**2) + (3*ee**2*KYR3**2*MY**4)/(2.*MW**2*sw**2))*cmath.sqrt(MB**4 - 2*MB**2*MW**2 + MW**4 - 2*MB**2*MY**2 - 2*MW**2*MY**2 + MY**4))/(96.*cmath.pi*abs(MY)**3)).real
#######
# 5FNS
#######
#VLT
def decayTP5(MTP,kappas,LR=None):
    KTLh3 = kappas[0]
    KTLw3 = kappas[1]
    KTLz3 = kappas[2]
    if LR == True:
        KTRh3 = kappas[0]
        KTRw3 = kappas[1]
        KTRz3 = kappas[2]
    else:
        KTRh3 = 0
        KTRw3 = 0
        KTRz3 = 0
    decayTP_Ht = ((-3*KTLh3**2*MH**2 - 3*KTRh3**2*MH**2 + 3*KTLh3**2*MT**2 + 3*KTRh3**2*MT**2 + 12*KTLh3*KTRh3*MT*MTP + 3*KTLh3**2*MTP**2 + 3*KTRh3**2*MTP**2)*cmath.sqrt(MH**4 - 2*MH**2*MT**2 + MT**4 - 2*MH**2*MTP**2 - 2*MT**2*MTP**2 + MTP**4))/(96.*cmath.pi*abs(MTP)**3)
    decayTP_Wb = ((MTP**2 - MW**2)*((3*ee**2*KTLw3**2*MTP**2)/(2.*sw**2) + (3*ee**2*KTRw3**2*MTP**2)/(2.*sw**2) + (3*ee**2*KTLw3**2*MTP**4)/(2.*MW**2*sw**2) + (3*ee**2*KTRw3**2*MTP**4)/(2.*MW**2*sw**2) - (3*ee**2*KTLw3**2*MW**2)/sw**2 - (3*ee**2*KTRw3**2*MW**2)/sw**2))/(96.*cmath.pi*abs(MTP)**3)      
    decayTP_Zt = (((3*ee**2*KTLz3**2*MT**2)/(4.*cw**2*sw**2) + (3*ee**2*KTRz3**2*MT**2)/(4.*cw**2*sw**2) - (9*ee**2*KTLz3*KTRz3*MT*MTP)/(cw**2*sw**2) + (3*ee**2*KTLz3**2*MTP**2)/(4.*cw**2*sw**2) + (3*ee**2*KTRz3**2*MTP**2)/(4.*cw**2*sw**2) + (3*ee**2*KTLz3**2*MT**4)/(4.*cw**2*MZ**2*sw**2) + (3*ee**2*KTRz3**2*MT**4)/(4.*cw**2*MZ**2*sw**2) - (3*ee**2*KTLz3**2*MT**2*MTP**2)/(2.*cw**2*MZ**2*sw**2) - (3*ee**2*KTRz3**2*MT**2*MTP**2)/(2.*cw**2*MZ**2*sw**2) + (3*ee**2*KTLz3**2*MTP**4)/(4.*cw**2*MZ**2*sw**2) + (3*ee**2*KTRz3**2*MTP**4)/(4.*cw**2*MZ**2*sw**2) - (3*ee**2*KTLz3**2*MZ**2)/(2.*cw**2*sw**2) - (3*ee**2*KTRz3**2*MZ**2)/(2.*cw**2*sw**2))*cmath.sqrt(MT**4 - 2*MT**2*MTP**2 + MTP**4 - 2*MT**2*MZ**2 - 2*MTP**2*MZ**2 + MZ**4))/(96.*cmath.pi*abs(MTP)**3)
    Gamma = decayTP_Ht + decayTP_Wb + decayTP_Zt
    return decayTP_Ht.real, decayTP_Wb.real, decayTP_Zt.real, Gamma.real 
#VLB
def decayBP5(MBP,kappas,LR=None):
    KBLh3 = kappas[0]
    KBLw3 = kappas[1]
    KBLz3 = kappas[2]
    if LR == True:
        KBRh3 = kappas[0]
        KBRw3 = kappas[1]
        KBRz3 = kappas[2]
    else:
        KBRh3 = 0
        KBRw3 = 0
        KBRz3 = 0
    decayBP_Hb = ((MBP**2 - MH**2)*(3*KBLh3**2*MBP**2 + 3*KBRh3**2*MBP**2 - 3*KBLh3**2*MH**2 - 3*KBRh3**2*MH**2))/(96.*cmath.pi*abs(MBP)**3)
    decayBP_Wt = (((3*ee**2*KBLw3**2*MBP**2)/(2.*sw**2) + (3*ee**2*KBRw3**2*MBP**2)/(2.*sw**2) - (18*ee**2*KBLw3*KBRw3*MBP*MT)/sw**2 + (3*ee**2*KBLw3**2*MT**2)/(2.*sw**2) + (3*ee**2*KBRw3**2*MT**2)/(2.*sw**2) + (3*ee**2*KBLw3**2*MBP**4)/(2.*MW**2*sw**2) + (3*ee**2*KBRw3**2*MBP**4)/(2.*MW**2*sw**2) - (3*ee**2*KBLw3**2*MBP**2*MT**2)/(MW**2*sw**2) - (3*ee**2*KBRw3**2*MBP**2*MT**2)/(MW**2*sw**2) + (3*ee**2*KBLw3**2*MT**4)/(2.*MW**2*sw**2) + (3*ee**2*KBRw3**2*MT**4)/(2.*MW**2*sw**2) - (3*ee**2*KBLw3**2*MW**2)/sw**2 - (3*ee**2*KBRw3**2*MW**2)/sw**2)*cmath.sqrt(MBP**4 - 2*MBP**2*MT**2 + MT**4 - 2*MBP**2*MW**2 - 2*MT**2*MW**2 + MW**4))/(96.*cmath.pi*abs(MBP)**3)
    decayBP_Zb = ((MBP**2 - MZ**2)*((3*ee**2*KBLz3**2*MBP**2)/(4.*cw**2*sw**2) + (3*ee**2*KBRz3**2*MBP**2)/(4.*cw**2*sw**2) + (3*ee**2*KBLz3**2*MBP**4)/(4.*cw**2*MZ**2*sw**2) + (3*ee**2*KBRz3**2*MBP**4)/(4.*cw**2*MZ**2*sw**2) - (3*ee**2*KBLz3**2*MZ**2)/(2.*cw**2*sw**2) - (3*ee**2*KBRz3**2*MZ**2)/(2.*cw**2*sw**2)))/(96.*cmath.pi*abs(MBP)**3)
    Gamma = decayBP_Hb + decayBP_Wt + decayBP_Zb
    return decayBP_Hb.real, decayBP_Wt.real, decayBP_Zb.real, Gamma.real
#VLX
def decayX5(MX,kappa,LR=None):
    KXL3 = kappa
    if LR == True:
        KXR3 = kappa
    else:
        KXR3 = 0
    return ((((3*ee**2*KXL3**2*MT**2)/(2.*sw**2) + (3*ee**2*KXR3**2*MT**2)/(2.*sw**2) + (3*ee**2*KXL3**2*MT**4)/(2.*MW**2*sw**2) + (3*ee**2*KXR3**2*MT**4)/(2.*MW**2*sw**2) - (3*ee**2*KXL3**2*MW**2)/sw**2 - (3*ee**2*KXR3**2*MW**2)/sw**2 - (18*ee**2*KXL3*KXR3*MT*MX)/sw**2 + (3*ee**2*KXL3**2*MX**2)/(2.*sw**2) + (3*ee**2*KXR3**2*MX**2)/(2.*sw**2) - (3*ee**2*KXL3**2*MT**2*MX**2)/(MW**2*sw**2) - (3*ee**2*KXR3**2*MT**2*MX**2)/(MW**2*sw**2) + (3*ee**2*KXL3**2*MX**4)/(2.*MW**2*sw**2) + (3*ee**2*KXR3**2*MX**4)/(2.*MW**2*sw**2))*cmath.sqrt(MT**4 - 2*MT**2*MW**2 + MW**4 - 2*MT**2*MX**2 - 2*MW**2*MX**2 + MX**4))/(96.*cmath.pi*abs(MX)**3)).real
#VLY
def decayY5(MY,kappa,LR=None):
    KYL3 = kappa
    if LR == True:
        KYR3 = kappa
    else:
        KYR3 = 0
    return (((-MW**2 + MY**2)*((-3*ee**2*KYL3**2*MW**2)/sw**2 - (3*ee**2*KYR3**2*MW**2)/sw**2 + (3*ee**2*KYL3**2*MY**2)/(2.*sw**2) + (3*ee**2*KYR3**2*MY**2)/(2.*sw**2) + (3*ee**2*KYL3**2*MY**4)/(2.*MW**2*sw**2) + (3*ee**2*KYR3**2*MY**4)/(2.*MW**2*sw**2)))/(96.*cmath.pi*abs(MY)**3)).real