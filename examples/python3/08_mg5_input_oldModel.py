import VLQcalc.model as model
import VLQcalc.madgraph as madgraph

vlq = model.VLQ('b')
vlq.setMass( range(1000,2001,200) )
vlq.BRs = [0.25,0.5,0.25]
vlq.KappaOld = [0.5]

mg5 = madgraph.MG5(vlq,'VLQ_UFO')

mg5.setProcess( 'p p > bp j, (bp > w- t, t > w+ b)' )

mg5.shower = 'Pythia8'
mg5.detector = 'Delphes'

mg5.addInput('Nevents 1000') # inputs for cards

mg5.createMG5Input('mg5_input') # name of input file