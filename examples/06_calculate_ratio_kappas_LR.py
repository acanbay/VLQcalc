import VLQcalc.model as model

Mass = range(1000, 2001, 200)
branchings = [0.25, 0.5, 0.25] # [BR(Hb), BR(Wt), BR(Zb)]
GM = 0.1 # Gamma/Mass ratio

vlq = model.VLQ('T', LR=True)
vlq.setMass(Mass)
vlq.calcRatioKappas(branchings, GM)

print('\nmT [GeV]\tkH\tkW\tkZ')
for i in range(len(Mass)):
    M = Mass[i]
    kH = round(vlq.KappaH[i], 3)
    kW = round(vlq.KappaW[i], 3)
    kZ = round(vlq.KappaZ[i], 3)
    print(f'  {M}\t\t{kH}\t{kW}\t{kZ}')
print()