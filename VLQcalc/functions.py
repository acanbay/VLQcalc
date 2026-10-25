'''
Author : Ali Can Canbay - acanbay
Created: 21 May 2024
Updated: 25 October 2026

This module contains the kappa calculation methods.
calcRatioKappas calculates kappas based on Gamma/Mass ratios.
Every partial width is proportional to kappa^2 (also with LR, where kappa_L = kappa_R), so for a target
branching ratio BR_c of channel c and a ratio Gamma/M the coupling is obtained directly:
    kappa_c = sqrt( BR_c * (Gamma/M) * M / Gamma_c(kappa=1) )
'''
import math
from VLQcalc import decays
###########################
def channelBRs(VLQ, BRs=None, zeta=None):
    # Target branching ratio of every channel: BR(V q_i) = zeta_i * BR(V)
    bosons = list(decays.CHANNELS[VLQ.type])
    BRs = VLQ.BRs if BRs is None else BRs
    zeta = VLQ.zeta if zeta is None else zeta
    BRs = BRs if type(BRs) == list else [BRs]
    return {(boson,gen): BR*z for boson, BR in zip(bosons,BRs) for gen, z in enumerate(zeta, start=1)}

def calcRatioKappas(VLQ, BRs, Ratio, zeta):
    # Returns {(boson, generation): [kappa per mass]}; the VLQ object is not changed
    targets = channelBRs(VLQ, BRs, zeta)
    Kappas = {channel: [] for channel in targets}
    for Mass in VLQ.Mass:
        unit = decays.partialWidths(VLQ.type, VLQ.FNS, Mass, {channel: 1 for channel in targets}, VLQ.LR)
        for channel, BR in targets.items():
            if BR == 0:
                kappa = 0
            elif unit[channel] <= 0:
                boson, gen = channel
                raise ValueError(f'VLQ.calcRatioKappas: {VLQ.type} -> {boson.upper()} + generation {gen} quark is kinematically '
                                 f'closed for M = {Mass} GeV, but its branching ratio is {BR}')
            else:
                # 6 significant digits, so that small couplings to light quarks keep their precision
                kappa = float(f'{math.sqrt(BR*Ratio*Mass/unit[channel]):.6g}')
            Kappas[channel].append(kappa)
    return Kappas
