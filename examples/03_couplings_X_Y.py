# X and Y quarks decay only through W (X -> W t, Y -> W b), so only Gamma/M is given.
import VLQcalc.model as model

for VLQ_type in ['X', 'Y']:
    vlq = model.VLQ(VLQ_type)
    vlq.setMass(range(1000, 2001, 250))
    vlq.calcRatioKappas(Ratio=0.05)
    print(vlq.table('txt'), '\n')
