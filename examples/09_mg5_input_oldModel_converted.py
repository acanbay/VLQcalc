import VLQcalc.model as model
import VLQcalc.madgraph as madgraph

Mass = range(1000, 2001, 200)
branchings = [0.25, 0.5, 0.25] # [BR(Hb), BR(Wt), BR(Zb)]
GM = 0.1 # Gamma/Mass ratio

vlq = model.VLQ('B')
vlq.setMass(Mass)
vlq.calcRatioKappas(branchings, GM)

newKappas = [vlq.KappaH, vlq.KappaW, vlq.KappaZ]
vlq.convertModel(newKappas, branchings, reverse=True)  # converts kappas VLQ_v4_UFO to VLQ_UFO

mg5 = madgraph.MG5(vlq, 'VLQ_UFO')

mg5.setProcess( 'p p > bp j, (bp > w- t, t > w+ b)' )

mg5.shower = 'Pythia8'
mg5.detector = 'Delphes'

mg5.addInput('Nevents 1000') # inputs for cards

mg5.createMG5Input('mg5_input') # name of input file