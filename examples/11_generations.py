# Couplings to 1st and 2nd generation quarks with the generation fractions zeta:
# BR(Q -> V q_i) = zeta_i * BR(V). Without zeta, all decays go to 3rd generation quarks.
import VLQcalc.model as model

vlq = model.VLQ('T')
vlq.setMass([1000, 1500, 2000])

# 10% of the decays into 1st, 20% into 2nd and 70% into 3rd generation quarks
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05, zeta=[0.1, 0.2, 0.7])

print(vlq.table('txt'))

# every coupling, named as in the VLQ_v4 and VLQ_v5 models
print('\nKTLw2 (T-W-s, left-handed):', vlq.Kappas['KTLw2'])
print('KTLw3 (T-W-b, left-handed):', vlq.Kappas['KTLw3'])
print('\ncouplings to lighter quarks:', vlq.KappaLight)  # {(boson, generation): values}
