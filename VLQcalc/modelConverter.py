'''
Author : Ali Can Canbay - acanbay
Created: 21 May 2024
Updated: 30 October 2026

Converts parameters between VLQ_UFO model [1,2] and VLQ_v4_#FNS_UFO (#: 4 or 5) models [1,3].
    reverse=False: VLQ_UFO --> VLQ_v4_#FNS_UFO
    reverse=True : VLQ_v4_#FNS_UFO --> VLQ_UFO

====REFERENCES====
[1] Alloul, A., Christensen, N. D., Degrande, C., Duhr, C., & Fuks, B. (2014). 
    FeynRules 2.0 - A complete toolbox for tree-level phenomenology. 
    Computer Physics Communications, 185(8), 2250-2300.
    https://feynrules.irmp.ucl.ac.be/
[2] Buchkremer, M., Cacciapaglia, G., Deandrea, A., & Panizzi, L. (2013). 
    Model-independent framework for searches of top partners. 
    Nuclear Physics B, 876(2), 376-417.
    https://feynrules.irmp.ucl.ac.be/wiki/VLQ
[3] Fuks, B., & Shao, H. S. (2017). 
    QCD next-to-leading-order predictions matched to parton showers for vector-like quark models. 
    The European Physical Journal C, 77(2), 135.
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
    if MVLQ <= {'h': MH, 'w': MW, 'z': MZ}[boson]:
        raise ValueError(f'VLQ.convertModel: the decay into {boson.upper()} is kinematically closed in VLQ_UFO '
                         f'for M = {MVLQ} GeV; set its branching ratio (xi) to zero or use a larger mass')
    if boson == 'h':
        return (MVLQ/v)*math.sqrt(1/GammaH(MVLQ))
    if boson == 'w':
        return math.sqrt(1/GammaW(MVLQ))*gaugeRatios()[0]
    return math.sqrt(1/GammaZ(MVLQ))*gaugeRatios()[1]

def convert(VLQ,Kappas,BRs=1,reverse=False,zeta=None):
    # reverse=False: returns {(boson, generation): [kappa per mass]} for VLQ_v4 from the VLQ_UFO kappa (Kappas)
    #                and the generation fractions zeta (default: those stored in VLQ)
    # reverse=True : returns the VLQ_UFO kappa per mass from the VLQ_v4 couplings stored in VLQ;
    #                the xi parameters of VLQ_UFO are stored in VLQ.BRs (T, B) and the zeta in VLQ.zeta
    z = zetaShare(VLQ)
    bosons = VLQ.bosons
    xis = BRs if type(BRs) == list else [BRs]
    if not reverse:
        zeta = VLQ.zeta if zeta is None else zeta
        kappas = {}
        for boson, xi in zip(bosons, xis):
            for gen, zg in enumerate(zeta, start=1):
                kappas[(boson,gen)] = [Kappas*math.sqrt(z*zg*xi)*vertexFactor(boson,MVLQ) if zg*xi != 0 else 0.
                                       for MVLQ in VLQ.Mass]
        return kappas
    ######
    # Exact inverse of the vertex matching: (kappa_Vi/F_V)^2 / s = kappa_Q^2 * zeta_i * xi_V, and since
    # sum(zeta_i) = sum(xi_V) = 1, the sums over channels give kappa_Q, xi_V and zeta_i directly.
    from VLQcalc.model import perMass
    import warnings
    third = {'h': VLQ.KappaH, 'w': VLQ.KappaW, 'z': VLQ.KappaZ}
    n = len(VLQ.Mass)
    couplings = {}
    for boson in bosons:
        for gen in [1,2,3]:
            values = third[boson] if gen == 3 else VLQ.KappaLight.get((boson,gen))
            couplings[(boson,gen)] = perMass(values, n, 'VLQ.convertModel') if values is not None else [0]*n
    KappaOld = []
    xisPerMass = [[] for boson in bosons]
    zetas = []
    exact = True
    for i, MVLQ in enumerate(VLQ.Mass):
        q = {channel: (values[i]/vertexFactor(channel[0],MVLQ))**2/z if values[i] != 0 else 0.
             for channel, values in couplings.items()}
        total = sum(q.values())
        if total == 0:
            KappaOld.append(0)
            for xi in xisPerMass:
                xi.append(0)
            zetas.append(None)
            continue
        xi = {boson: sum(v for (b,g), v in q.items() if b == boson)/total for boson in bosons}
        zeta = [sum(v for (b,g), v in q.items() if g == gen)/total for gen in [1,2,3]]
        # VLQ_UFO describes only couplings of the form kappa_Q^2 * zeta_i * xi_V
        if any(not math.isclose(v, total*xi[b]*zeta[g-1], rel_tol=1e-4, abs_tol=1e-12*total) for (b,g), v in q.items()):
            exact = False
        KappaOld.append( float(f'{math.sqrt(total):.6g}') )
        for values, boson in zip(xisPerMass, bosons):
            values.append( float(f'{xi[boson]:.6g}') )
        zetas.append([float(f'{zg:.6g}') for zg in zeta])
    known = [zg for zg in zetas if zg is not None]
    sameZeta = all(all(math.isclose(a, b, rel_tol=1e-4, abs_tol=1e-9) for a, b in zip(zg, known[0])) for zg in known)
    if exact and known and sameZeta:
        # generation fractions written to the VLQ_UFO card
        VLQ.zeta = known[0]
    elif known:
        warnings.warn('VLQ.convertModel: VLQ_UFO cannot reproduce every coupling exactly (its couplings have the form '
                      'kappa_Q*sqrt(zeta_i*xi_V); with couplings to lighter quarks from calcRatioKappas the quark masses '
                      'break this form). kappa_Q and xi reproduce the total coupling strength and its share per boson; '
                      'the zeta of the VLQ object are kept', stacklevel=3)
    if VLQ.type in ['VLT','VLB']:
        # xi parameters of VLQ_UFO per mass, written to the VLQ_UFO card
        VLQ.BRs = xisPerMass
    return KappaOld
