import VLQcalc.model as model

vlq = model.VLQ('T')
vlq.setMass( range(1000, 2001, 250) )

# 10% of the decays into 1st, 20% into 2nd and 70% into 3rd generation quarks
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05, zeta=[0.1, 0.2, 0.7])

print(vlq.table())

# every coupling with the name used in the VLQ_v4 and VLQ_v5 models
for name, values in vlq.Kappas.items():
    print(name, values)
