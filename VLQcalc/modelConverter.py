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
############
MH=125.10 
MW=80.379 
MZ=91.1876
gWeak =  0.65
v = 2*MW/gWeak
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
def convert(VLQ,Kappas,BRs=1,reverse=False):
    if not reverse:
        if VLQ.type in ['VLX','VLY']:
            kappaWs=[]
            for MVLQ in VLQ.Mass:
                kappaWs.append( Kappas*(math.sqrt(1/GammaW(MVLQ))) )
            return kappaWs
        ###
        else:
            kappaHs=[]
            kappaWs=[]
            kappaZs=[]
            for MVLQ in VLQ.Mass:
                kappaHs.append( Kappas*((MVLQ/v)*math.sqrt(BRs[0]/GammaH(MVLQ))) )
                kappaWs.append( Kappas*(math.sqrt(BRs[1]/GammaW(MVLQ))) )
                kappaZs.append( Kappas*(math.sqrt(BRs[2]/GammaZ(MVLQ))) )
            return kappaHs, kappaWs, kappaZs
    ######
    else:
        if VLQ.type in ['VLX','VLY']:
            KappaOld=[]
            for i in range(len(VLQ.Mass)):
                KappaOld.append( Kappas[i]/(math.sqrt(1/GammaW(VLQ.Mass[i]))) )
            return KappaOld
        ###
        else:
            KappaOld=[]
            for i in range(len(VLQ.Mass)):
                DH, DW, DZ, Gamma = VLQ.calcDecay(VLQ.Mass[i],[ Kappas[0][i],Kappas[1][i],Kappas[2][i] ],VLQ.LR)
                BRH = DH/Gamma
                BRW = DW/Gamma
                BRZ = DZ/Gamma
                kappa = []
                if BRH!=0:
                    kappa.append( Kappas[0][i]/((VLQ.Mass[i]/v)*math.sqrt(BRH/GammaH(VLQ.Mass[i]))) )
                if BRW!=0:
                    kappa.append( Kappas[1][i]/(math.sqrt(BRW/GammaW(VLQ.Mass[i]))) )
                if BRZ!=0:
                    kappa.append( Kappas[2][i]/(math.sqrt(BRZ/GammaZ(VLQ.Mass[i]))) )
                if len(kappa)!=0:
                    KappaOld.append( sum(kappa)/len(kappa) )
                else:
                    KappaOld.append(0)
            return KappaOld