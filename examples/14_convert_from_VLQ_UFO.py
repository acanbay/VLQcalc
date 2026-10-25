# Conversion of the coupling kappa_B of VLQ_UFO to the couplings of VLQ_v4_UFO and VLQ_v5_UFO.
import VLQcalc.model as model

vlq = model.VLQ('B')
vlq.setMass([1000, 1500, 2000])

kB = 0.5
xi = [0.25, 0.5, 0.25]                       # branching ratios of VLQ_UFO: [xi(H), xi(W), xi(Z)]
vlq.convertModel(kB, xi)

print(vlq.table('txt'))
