# Conversion of VLQ_v4_UFO couplings to the coupling kappa_T of VLQ_UFO (reverse=True).
import VLQcalc.model as model

vlq = model.VLQ('T')
vlq.setMass([1000, 1500, 2000])
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.1)  # [BR(Ht), BR(Wb), BR(Zt)], Gamma/M

kappas = [vlq.KappaH, vlq.KappaW, vlq.KappaZ]
vlq.convertModel(kappas, reverse=True)

print(vlq.table('txt'))
print('\nkT of VLQ_UFO:', vlq.KappaOld)      # one value per mass
