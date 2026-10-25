# SM inputs: the widths use the inputs of the VLQ_v4 and VLQ_v5 restriction cards.
# If the param_card of MadGraph5 is changed, the same values have to be given to VLQcalc.
import VLQcalc.model as model
import VLQcalc.constants as constants

print('default SM inputs:', constants.getSM())

vlq = model.VLQ('T')
vlq.setMass([1000, 1500, 2000])
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05)  # [BR(Ht), BR(Wb), BR(Zt)], Gamma/M
print('kH, kZ with the default inputs:', vlq.KappaH, vlq.KappaZ)

constants.setSM(Mt=173.0, MH=125.1)          # only the given inputs change
print('\nchanged SM inputs:', constants.getSM())
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05)
print('kH, kZ with the changed inputs:', vlq.KappaH, vlq.KappaZ)  # T -> Ht and T -> Zt depend on Mt and MH

constants.setSM()                            # back to the default inputs
