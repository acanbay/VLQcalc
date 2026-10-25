# VLQ_UFO card for a chosen width-to-mass ratio: the couplings are computed for VLQ_v4_UFO
# and converted to VLQ_UFO before the card is written.
import VLQcalc.model as model
import VLQcalc.madgraph as madgraph

vlq = model.VLQ('B')
vlq.setMass(range(1000, 2001, 200))
vlq.calcRatioKappas([0.25, 0.5, 0.25], 0.1)  # [BR(Hb), BR(Wt), BR(Zb)], Gamma/M
vlq.convertModel([vlq.KappaH, vlq.KappaW, vlq.KappaZ], reverse=True)

mg5 = madgraph.MG5(vlq, 'VLQ_UFO')
mg5.setProcess('p p > bp j, (bp > w- t, t > w+ b)')
mg5.addInput('nevents 1000')

mg5.createMG5Input('bp_single_VLQ_UFO_converted')
