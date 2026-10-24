import VLQcalc.model as model

vlq = model.VLQ('B')
vlq.setMass( [1000, 1500, 2000] )

###################################
kB = 0.5
branchings = [0.25, 0.5, 0.25] # [xiH, xiW, xiZ] of VLQ_UFO
zeta = [0, 0.3, 0.7]           # generation fractions of VLQ_UFO
vlq.convertModel(kB, branchings, zeta=zeta)  # converts kappa from VLQ_UFO to VLQ_v4_UFO

print(vlq.table())
