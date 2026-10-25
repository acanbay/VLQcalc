# MadGraph5 input card for X pair production with VLQ_v5_UFO.
import VLQcalc.model as model
import VLQcalc.madgraph as madgraph

vlq = model.VLQ('X', FNS=5)
vlq.setMass(range(1000, 2001, 200))
vlq.calcRatioKappas(Ratio=0.1)

mg5 = madgraph.MG5(vlq, 'VLQ_v5_UFO')
mg5.setProcess('p p > x x~, x > w+ t')

mg5.addInput('nevents 1000')
mg5.addInput('WX auto')                      # MadGraph5 computes the width of X from the couplings

mg5.createMG5Input('x_pair')
