import VLQcalc.model as model

vlq = model.VLQ('T', LR=True)
kappas = [0.984, 0.272, 0.285]

print('\nmT [GeV]\tGamma [GeV]\tBR(Ht)\tBR(Wb)\tBR(Zt)')
for mass in range(1000,2001,200):
    decayH, decayW, decayZ, Gamma = vlq.calcDecay(mass, kappas, vlq.LR)
    BRH = round(decayH/Gamma, 3)
    BRW = round(decayW/Gamma, 3)
    BRZ = round(decayZ/Gamma, 3)
    Gamma = round(Gamma, 3)
    print(f'  {mass}\t\t{Gamma}\t\t{BRH}\t{BRW}\t{BRZ}')
print()