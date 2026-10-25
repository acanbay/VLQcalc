# The same as example 01 with left- and right-handed couplings (kappa_R = kappa_L).
# Both chiralities contribute to the width, so each coupling is smaller.
import VLQcalc.model as model

vlq = model.VLQ('T', LR=True)
vlq.setMass(range(1000, 2001, 250))
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05) # [BR(Ht), BR(Wb), BR(Zt)], Gamma/M

print(vlq.table('txt'))
