# Couplings of a T quark for a chosen width-to-mass ratio and chosen branching ratios.
# Only left-handed couplings to 3rd generation quarks (default settings).
import VLQcalc.model as model

vlq = model.VLQ('T')                         # FNS=4, LR=False by default
vlq.setMass(range(1000, 2001, 250))          # masses in GeV

branchings = [0.5, 0.25, 0.25]               # [BR(Ht), BR(Wb), BR(Zt)], sum 1
vlq.calcRatioKappas(branchings, 0.05)        # Gamma/M = 0.05

print(vlq.table('txt'))
print('\nkH:', vlq.KappaH)                   # one value per mass
print('kW:', vlq.KappaW)
print('kZ:', vlq.KappaZ)
