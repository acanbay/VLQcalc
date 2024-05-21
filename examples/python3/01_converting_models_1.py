import VLQcalc.model as model

Mass = [1000,1500,2000]
branchings = [0.25,0.5,0.25] # [BR(Hb),BR(Wt),BR(Zb)]
vlq = model.VLQ('bp')
vlq.setMass( Mass )

###################################
kB = 0.5
vlq.convertModel( kB, branchings )  # converts kappa VLQ_UFO to VLQ_v4

print('\nVLQ_v4:')
print('MBP [GeV]\tkH\tkW\tkZ')
for i in range(len(Mass)):
    M = Mass[i]
    kH = round(vlq.KappaH[i],2)
    kW = round(vlq.KappaW[i],2)
    kZ = round(vlq.KappaZ[i],2)
    print(f'  {M}\t\t{kH}\t{kW}\t{kZ}')

###################################
newKappas = [vlq.KappaH,vlq.KappaW,vlq.KappaZ]
vlq.convertModel(newKappas,branchings,reverse=True)  # converts kappas VLQ_v4 to VLQ_UFO

print('VLQ_UFO:')
print('MBP [GeV]\tkB')
for i in range(len(Mass)):
    M = Mass[i]
    kB = round(vlq.KappaOld[i],2)
    print(f'  {M}\t\t{kB}')
print()