# Conversion from VLQ_UFO with generation fractions zeta: the couplings to lighter quarks
# of VLQ_v4_UFO and VLQ_v5_UFO are computed as well.
import VLQcalc.model as model

vlq = model.VLQ('B')
vlq.setMass([1000, 1500, 2000])
vlq.convertModel(0.5, [0.25, 0.5, 0.25], zeta=[0, 0.3, 0.7])

print(vlq.table('txt'))
