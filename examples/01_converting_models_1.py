import VLQcalc.model as model

Mass = [1000, 1500, 2000]
branchings = [0.25, 0.5, 0.25] # [BR(Hq), BR(Wq), BR(Zq)]
vlq = model.VLQ('B')
vlq.setMass( Mass )

###################################
kB = 0.5
vlq.convertModel(kB, branchings)  # converts kappa from VLQ_UFO to VLQ_v4_UFO

print('\nVLQ_v4_UFO:')
print('mB [GeV]\tkH\tkW\tkZ')
for i in range(len(Mass)):
    M = Mass[i]
    kH = round(vlq.KappaH[i], 2)
    kW = round(vlq.KappaW[i], 2)
    kZ = round(vlq.KappaZ[i], 2)
    print(f'{M}\t\t{kH}\t{kW}\t{kZ}')

###################################
newKappas = [vlq.KappaH, vlq.KappaW, vlq.KappaZ]
vlq.convertModel(newKappas, branchings, True)  # converts kappas from VLQ_v4_UFO to VLQ_UFO

print('VLQ_UFO:')
print('mB [GeV]\tkB')
for i in range(len(Mass)):
    M = Mass[i]
    kB = round(vlq.KappaOld[i], 2)
    print(f'{M}\t\t{kB}')
print()