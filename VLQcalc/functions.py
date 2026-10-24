'''
Author: Ali Can Canbay - acanbay
Date  : 21 May 2024

This module contains the kappa calculation methods.
calcRatioKappas calculates kappas based on Gamma/Mass ratios.
Every partial width is proportional to kappa^2 (also with LR, where kappa_L = kappa_R), so for a target
branching ratio BR_c of channel c and a ratio Gamma/M the coupling is obtained directly:
    kappa_c = sqrt( BR_c * (Gamma/M) * M / Gamma_c(kappa=1) )
Versions up to 0.1.2 used a modified genetic algorithm [1] for the same task.

====REFERENCES====
[1] Canbay, A. C. modifiedGA (Version 0.1.0) [Computer software].
    https://doi.org/10.5281/zenodo.12569505
    https://github.com/acanbay/modifiedGA

'''
import math
from VLQcalc import decays
###########################
def channelBRs(VLQ):
    # Target branching ratio of every channel: BR(V q_i) = zeta_i * BR(V)
    bosons = list(decays.CHANNELS[VLQ.type])
    BRs = VLQ.BRs if type(VLQ.BRs) == list else [VLQ.BRs]
    return {(boson,gen): BR*zeta for boson, BR in zip(bosons,BRs) for gen, zeta in enumerate(VLQ.zeta, start=1)}

def calcRatioKappas(VLQ):
    targets = channelBRs(VLQ)
    Kappas = {channel: [] for channel in targets}
    for Mass in VLQ.Mass:
        unit = decays.partialWidths(VLQ.type, VLQ.FNS, Mass, {channel: 1 for channel in targets}, VLQ.LR)
        for channel, BR in targets.items():
            if BR == 0:
                kappa = 0
            elif unit[channel] <= 0:
                boson, gen = channel
                raise ValueError(f'{VLQ.type} -> {boson.upper()} + generation {gen} quark is kinematically closed for M = {Mass} GeV, '
                                 f'but its branching ratio is {BR}')
            else:
                # 6 significant digits, so that small couplings to light quarks keep their precision
                kappa = float(f'{math.sqrt(BR*VLQ.Ratio*Mass/unit[channel]):.6g}')
            Kappas[channel].append(kappa)
    VLQ.setKappas(Kappas)
