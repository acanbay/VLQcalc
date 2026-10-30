# MadGraph5 input card for single B production with VLQ_v4_UFO in the 4-flavour scheme.
# The card is run with: ./bin/mg5_aMC bp_single_4FNS.dat
import VLQcalc.model as model
import VLQcalc.madgraph as madgraph

vlq = model.VLQ('B', FNS=4, LR=True)
vlq.setMass(range(1000, 2001, 200))          # one MadGraph5 run per mass
vlq.calcRatioKappas([0.25, 0.5, 0.25], 0.1)  # [BR(Hb), BR(Wt), BR(Zb)], Gamma/M

mg5 = madgraph.MG5(vlq, 'VLQ_v4_UFO')
mg5.setProcess('p p > bp j, (bp > w- t, t > w+ b)')
mg5.addProcess('p p > bp j b, (bp > w- t, t > w+ b)')

mg5.shower = 'Pythia8'
mg5.detector = 'Delphes'

mg5.addInput('nevents 1000')                 # settings for the run cards
# The card sets WBP auto: MadGraph5 computes the width of B from the couplings

mg5.createMG5Input('bp_single_4FNS')         # writes bp_single_4FNS.dat
