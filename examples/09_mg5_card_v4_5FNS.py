# The same card as example 08 in the 5-flavour scheme (massless b quark, b in the proton).
import VLQcalc.model as model
import VLQcalc.madgraph as madgraph

vlq = model.VLQ('B', FNS=5, LR=True)
vlq.setMass(range(1000, 2001, 200))
vlq.calcRatioKappas([0.25, 0.5, 0.25], 0.1)  # [BR(Hb), BR(Wt), BR(Zb)], Gamma/M

mg5 = madgraph.MG5(vlq, 'VLQ_v4_UFO')
mg5.setProcess('p p > bp j, (bp > w- t, t > w+ b)')
mg5.addProcess('p p > bp j b, (bp > w- t, t > w+ b)')

mg5.addInput('nevents 1000')
mg5.addInput('WBP auto')

mg5.createMG5Input('bp_single_5FNS')
