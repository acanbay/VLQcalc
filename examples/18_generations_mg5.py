import VLQcalc.model as model
import VLQcalc.madgraph as madgraph

vlq = model.VLQ('T')
vlq.setMass( range(1000, 2001, 250) )

# 20% of the decays into 2nd and 80% into 3rd generation quarks
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05, zeta=[0, 0.2, 0.8])

mg5 = madgraph.MG5(vlq, 'VLQ_v4_UFO')

mg5.setProcess( 'p p > tp tp~' )

mg5.addInput('WTP auto')     # MadGraph5 computes the width of T from the couplings

mg5.createMG5Input('tp_generations') # every coupling of T (1st, 2nd and 3rd generation) is written
