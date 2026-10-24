'''
Author: Ali Can Canbay - acanbay
Date  : 21 May 2024

Converts parameters between VLQ_UFO model [1,2] and VLQ_v4_#FNS_UFO (#: 4 or 5) models [1,3].
    reverse=False: VLQ_UFO --> VLQ_v4_#FNS_UFO
    reverse=True : VLQ_v4_#FNS_UFO --> VLQ_UFO

====REFERENCES====
[1] Alloul, A., Christensen, N. D., Degrande, C., Duhr, C., & Fuks, B. (2014). 
    FeynRules 2.0—A complete toolbox for tree-level phenomenology. 
    Computer Physics Communications, 185(8), 2250-2300.
    https://feynrules.irmp.ucl.ac.be/
[2] Buchkremer, M., Cacciapaglia, G., Deandrea, A., & Panizzi, L. (2013). 
    Model-independent framework for searches of top partners. 
    Nuclear Physics B, 876(2), 376-417.
    https://feynrules.irmp.ucl.ac.be/wiki/VLQ
[3] Fuks, B., & Shao, H. S. (2017). 
    QCD next-to-leading-order predictions matched to parton showers for vector-like quark models. 
    The European Physical Journal C, 77(2), 1-21.
    https://feynrules.irmp.ucl.ac.be/wiki/NLOModels

'''
import math
from VLQcalc import constants
############
# SM inputs of the VLQ_UFO param_card [2]; MW, gWeak and v are derived as in that model
MH=125.
MZ=91.1876
aEW=1./127.9
Gf=1.166e-5
MW=math.sqrt(MZ**2/2. + math.sqrt(MZ**4/4. - (aEW*math.pi*MZ**2)/(Gf*math.sqrt(2))))
sw=math.sqrt(1 - MW**2/MZ**2)
cw=math.sqrt(1 - sw**2)
gWeak = 2*math.sqrt(aEW*math.pi)/sw
v = 2*MW/gWeak
def gaugeRatios():
    # Ratios of the W and Z gauge factors of VLQ_UFO to those of VLQ_v4 (constants.py)
    gV4 = (constants.ee/constants.sw).real
    return gWeak/gV4, (gWeak/cw)/(gV4/constants.cw.real)
############
def Lambda(a,b,c):
    return a**2 + b**2 + c**2 - 2*a*b - 2*a*c - 2*b*c
def GammaH(MVLQ,mq=0):
    a = 1
    b = mq**2/MVLQ**2
    c = MH**2/MVLQ**2
    return math.sqrt(Lambda(a,b,c))*( 1 + b - c )/2
def GammaW(MVLQ,mq=0):
    a = 1
    b = mq**2/MVLQ**2
    c = MW**2/MVLQ**2
    return math.sqrt(Lambda(a,b,c))*( (1-b)**2 + c - 2*c**2 + b*c )
def GammaZ(MVLQ,mq=0):
    a = 1
    b = mq**2/MVLQ**2
    c = MZ**2/MVLQ**2
    return math.sqrt(Lambda(a,b,c))*( (1-b)**2 + c - 2*c**2 + b*c )/2
############
def zetaShare(VLQ):
    # VLQ_UFO normalises the zeta parameters to 1; with LR the coupling is shared equally (zetaL = zetaR = 0.5)
    return 0.5 if VLQ.LR else 1.

def vertexFactor(boson,MVLQ):
    # kappa(VLQ_v4) = kappa(VLQ_UFO) * sqrt(zeta_i * xi_V) * vertexFactor, from matching the vertices of both models
    if boson == 'h':
        return (MVLQ/v)*math.sqrt(1/GammaH(MVLQ))
    if boson == 'w':
        return math.sqrt(1/GammaW(MVLQ))*gaugeRatios()[0]
    return math.sqrt(1/GammaZ(MVLQ))*gaugeRatios()[1]

def convert(VLQ,Kappas,BRs=1,reverse=False):
    # reverse=False: returns {(boson, generation): [kappa per mass]} for VLQ_v4 from the VLQ_UFO kappa (Kappas)
    # reverse=True : returns the VLQ_UFO kappa per mass from the VLQ_v4 couplings stored in VLQ
    z = zetaShare(VLQ)
    bosons = VLQ.bosons
    xis = BRs if type(BRs) == list else [BRs]
    if not reverse:
        kappas = {}
        for boson, xi in zip(bosons, xis):
            for gen, zeta in enumerate(VLQ.zeta, start=1):
                kappas[(boson,gen)] = [Kappas*math.sqrt(z*zeta*xi)*vertexFactor(boson,MVLQ) for MVLQ in VLQ.Mass]
        return kappas
    ######
    KappaOld=[]
    xisPerMass=[[] for boson in bosons]
    for i, MVLQ in enumerate(VLQ.Mass):
        widths = VLQ.calcWidths(i)
        Gamma = sum(widths.values())
        couplings = {}
        third = {'h': VLQ.KappaH, 'w': VLQ.KappaW, 'z': VLQ.KappaZ}
        for (boson,gen) in widths:
            values = third[boson] if gen == 3 else VLQ.KappaLight.get((boson,gen), [0]*len(VLQ.Mass))
            couplings[(boson,gen)] = values[i]
        kappa = []
        for channel, width in widths.items():
            BR = width/Gamma if Gamma > 0 else 0
            if BR != 0:
                kappa.append( couplings[channel]/(math.sqrt(z*BR)*vertexFactor(channel[0],MVLQ)) )
        for xi, boson in zip(xisPerMass, bosons):
            xi.append( round(sum(w for (b,g), w in widths.items() if b == boson)/Gamma, 6) if Gamma > 0 else 0 )
        if len(kappa)!=0:
            KappaOld.append( sum(kappa)/len(kappa) )
        else:
            KappaOld.append(0)
    if VLQ.BRs == None and VLQ.type in ['VLT','VLB']:
        # Branchings per mass, needed as xi parameters for the VLQ_UFO card
        VLQ.BRs = xisPerMass
    return KappaOld
