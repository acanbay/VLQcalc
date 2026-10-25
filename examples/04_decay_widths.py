# Total width and branching ratios of a T quark for given couplings,
# with left-handed couplings only (LR=False) and with kappa_R = kappa_L (LR=True).
import VLQcalc.model as model

kappas = [1.0, 0.3, 0.3]                     # [kH, kW, kZ] of the 3rd generation

for LR in [False, True]:
    vlq = model.VLQ('T', LR=LR)
    print(f'\nLR = {LR}')
    print('M [GeV]   Gamma [GeV]   BR(Ht)   BR(Wb)   BR(Zt)')
    for M in range(1000, 2001, 250):
        widthH, widthW, widthZ, Gamma = vlq.calcDecay(M, kappas, LR)
        print(f'{M:7d}   {Gamma:11.4f}   {widthH/Gamma:6.4f}   {widthW/Gamma:6.4f}   {widthZ/Gamma:6.4f}')
