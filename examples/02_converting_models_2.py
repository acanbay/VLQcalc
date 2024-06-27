import VLQcalc.model as model

Mass = [1000, 1500, 2000]
vlq = model.VLQ('Y', FNS=5, LR=True)
vlq.setMass(Mass)

###################################
kY = 0.5
vlq.convertModel(kY)  # converts kappa from VLQ_UFO to VLQ_v4_UFO

print('\nVLQ_v4_UFO:')
print('mY [GeV]\tkW')
for i in range(len(Mass)):
    M = Mass[i]
    kW = round(vlq.KappaW[i], 2)
    print(f'{M}\t\t{kW}')

###################################
newKappas = vlq.KappaW
vlq.convertModel(newKappas, reverse=True)  # converts kappas from VLQ_v4_UFO to VLQ_UFO

print('VLQ_UFO:')
print('mY [GeV]\tkY')
for i in range(len(Mass)):
    M = Mass[i]
    kY = round(vlq.KappaOld[i], 2)
    print(f'{M}\t\t{kY}')
print()