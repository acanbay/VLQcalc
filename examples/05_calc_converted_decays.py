
import VLQcalc.model as model

Mass = range(1000, 2001, 200)
branchings = [0.25, 0.5, 0.25] # [BR(Hb), BR(Wt), BR(Zb)]
vlq = model.VLQ('T')
vlq.setMass( Mass )

kT = 0.5
vlq.convertModel(kT, branchings)  # converts kappa from VLQ_UFO to VLQ_v4_UFO

kappas = []
for i in range(len(Mass)):
    kappas.append( [vlq.KappaH[i], vlq.KappaW[i], vlq.KappaZ[i]] )

print('\nmT [GeV]\tGamma [GeV]\tBR(Ht)\tBR(Wb)\tBR(Zt)')
for i in range(len(Mass)):
    decayH, decayW, decayZ, Gamma = vlq.calcDecay(Mass[i], kappas[i])
    BRH = round(decayH/Gamma, 3)
    BRW = round(decayW/Gamma, 3)
    BRZ = round(decayZ/Gamma, 3)
    Gamma = round(Gamma, 3)
    print(f'{Mass[i]}\t\t{Gamma}\t\t{BRH}\t{BRW}\t{BRZ}')
print()
