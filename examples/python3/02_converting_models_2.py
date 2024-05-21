import VLQcalc.model as model

Mass = [1000,1500,2000]
vlq = model.VLQ("y",FNS=5,LR=True)
vlq.setMass( Mass )

###################################
kY = 0.5
vlq.convertModel( kY )  # converts kappa VLQ_UFO to VLQ_v4

print('\nVLQ_v4:')
print('MY [GeV]\tkW')
for i in range(len(Mass)):
    M = Mass[i]
    kW = round(vlq.KappaW[i],2)
    print(f' {M}\t\t{kW}')

###################################
newKappas = vlq.KappaW
vlq.convertModel(newKappas,reverse=True)  # converts kappas VLQ_v4 to VLQ_UFO

print('VLQ_UFO:')
print('MY [GeV]\tkY')
for i in range(len(Mass)):
    M = Mass[i]
    kB = round(vlq.KappaOld[i],2)
    print(f' {M}\t\t{kB}')
print()