# Conversion for X and Y quarks, which decay only through W, so no branching ratios are given.
import VLQcalc.model as model

for VLQ_type in ['X', 'Y']:
    vlq = model.VLQ(VLQ_type, FNS=5, LR=True)
    vlq.setMass([1000, 1500, 2000])
    vlq.convertModel(0.5)                    # kappa_X or kappa_Y of VLQ_UFO
    print(vlq.table('txt'), '\n')
