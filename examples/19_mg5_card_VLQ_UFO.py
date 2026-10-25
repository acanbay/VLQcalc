# MadGraph5 card for VLQ_UFO with a coupling kappa_B and branching ratios xi set by hand.
import VLQcalc.model as model
import VLQcalc.madgraph as madgraph

vlq = model.VLQ('B')
vlq.setMass(range(1000, 2001, 200))
vlq.BRs = [0.25, 0.5, 0.25]                  # [xi(H), xi(W), xi(Z)]
vlq.KappaOld = [0.5]                         # kappa_B, the same for every mass

mg5 = madgraph.MG5(vlq, 'VLQ_UFO')
mg5.setProcess('p p > bp j, (bp > w- t, t > w+ b)')
mg5.addInput('nevents 1000')

mg5.createMG5Input('bp_single_VLQ_UFO')
