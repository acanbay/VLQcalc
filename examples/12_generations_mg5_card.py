# MadGraph5 card with couplings to lighter quarks. Every coupling of the VLQ
# (1st, 2nd and 3rd generation, both chiralities) is written to the card.
import VLQcalc.model as model
import VLQcalc.madgraph as madgraph

vlq = model.VLQ('T')
vlq.setMass(range(1000, 2001, 250))
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05, zeta=[0, 0.2, 0.8])  # [BR(Ht), BR(Wb), BR(Zt)]

mg5 = madgraph.MG5(vlq, 'VLQ_v4_UFO')
mg5.setProcess('p p > tp tp~')

mg5.createMG5Input('tp_pair_generations')
