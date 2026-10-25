'''
Author : Ali Can Canbay - acanbay
Created: 21 May 2024
Updated: 25 October 2026

Calculates decay widths for VLQ_v4_4FNS_UFO and VLQ_v4_5FNS_UFO models [1,2].
The formulas are those of the decays.py files of the models, written for a general quark mass;
SM inputs are read from constants.py.

====REFERENCES====
[1] Alloul, A., Christensen, N. D., Degrande, C., Duhr, C., & Fuks, B. (2014). 
    FeynRules 2.0—A complete toolbox for tree-level phenomenology. 
    Computer Physics Communications, 185(8), 2250-2300.
    https://feynrules.irmp.ucl.ac.be/
[2] Fuks, B., & Shao, H. S. (2017). 
    QCD next-to-leading-order predictions matched to parton showers for vector-like quark models. 
    The European Physical Journal C, 77(2), 135.
    https://feynrules.irmp.ucl.ac.be/wiki/NLOModels
'''
import cmath
from VLQcalc import constants as c
#######
# 3rd generation (VLQ.calcDecay)
#######
def thirdGeneration(VLQtype,FNS,M,Kappas,LR):
    # Partial widths of the 3rd generation channels for Kappas = [kH, kW, kZ] (T, B) or kW (X, Y)
    bosons = list(CHANNELS[VLQtype])
    Kappas = list(Kappas) if len(bosons) > 1 else [Kappas]
    widths = partialWidths(VLQtype, FNS, M, {(boson,3): kappa for boson, kappa in zip(bosons, Kappas)}, LR)
    return [widths[(boson,3)] for boson in bosons]

def decayT(FNS):
    def decay(MtP,Kappas,LR=False):
        widthH, widthW, widthZ = thirdGeneration('VLT', FNS, MtP, Kappas, LR)
        return widthH, widthW, widthZ, widthH + widthW + widthZ
    return decay

def decayB(FNS):
    def decay(MbP,Kappas,LR=False):
        widthH, widthW, widthZ = thirdGeneration('VLB', FNS, MbP, Kappas, LR)
        return widthH, widthW, widthZ, widthH + widthW + widthZ
    return decay

def decayXY(VLQtype,FNS):
    def decay(M,kappa,LR=False):
        return thirdGeneration(VLQtype, FNS, M, kappa, LR)[0]
    return decay

# (H, W, Z, total) for T and B, the W width for X and Y
decayTP4, decayBP4, decayX4, decayY4 = decayT(4), decayB(4), decayXY('VLX',4), decayXY('VLY',4)
decayTP5, decayBP5, decayX5, decayY5 = decayT(5), decayB(5), decayXY('VLX',5), decayXY('VLY',5)
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
    # Q -> H q, as in decays.py of the models with a general quark mass; zero below the threshold
    if M <= mq + c.MH:
        return 0.
    return ((3*(kL**2 + kR**2)*(mq**2 + M**2 - c.MH**2) + 12*kL*kR*mq*M)*sqrtLambda(M,mq,c.MH)/(96.*cmath.pi*abs(M)**3)).real

def widthV(M,mq,MV,g2,kL,kR):
    # Q -> V q (V = W, Z); g2 is the squared gauge factor of the vertex; zero below the threshold
    if M <= mq + MV:
        return 0.
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
