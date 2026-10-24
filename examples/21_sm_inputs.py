import VLQcalc.model as model
import VLQcalc.constants as constants

print('default SM inputs:', constants.getSM())

vlq = model.VLQ('T')
vlq.setMass( [1000, 1500, 2000] )
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05)
print('kH, kZ with default inputs:', vlq.KappaH, vlq.KappaZ)

# if the param_card of MadGraph5 is changed, give the same values to VLQcalc
constants.setSM(Mt=173.0, MH=125.1)
print('\nchanged SM inputs:', constants.getSM())

vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05)
print('kH, kZ with changed inputs:', vlq.KappaH, vlq.KappaZ) # T -> Ht and T -> Zt depend on Mt and MH

constants.setSM() # back to the default inputs
