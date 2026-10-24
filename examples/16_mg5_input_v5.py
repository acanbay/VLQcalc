import VLQcalc.model as model
import VLQcalc.madgraph as madgraph

vlq = model.VLQ('X', FNS=5, LR=True)
vlq.setMass( range(1000, 2001, 200) )
vlq.calcRatioKappas(Ratio=0.1)

mg5 = madgraph.MG5(vlq, 'VLQ_v5_UFO')

mg5.setProcess( 'p p > x x~, x > w+ t' )

mg5.addInput('nevents 1000') # inputs for cards
mg5.addInput('WX auto')      # MadGraph5 computes the width of X from the couplings

mg5.createMG5Input('mg5_input') # name of input file
