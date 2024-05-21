import VLQcalc.model as model

vlq = model.VLQ('tp')
kappas = [1.612, 0.384, 0.402]

print('\nMTP [GeV]\tGamma [GeV]\tBR(Ht)\tBR(Wb)\tBR(Zt)')
for mass in range(1000,2001,200):
    decayH, decayW, decayZ, Gamma = vlq.calcDecay(mass,kappas)
    BRH = round(decayH/Gamma,3)
    BRW = round(decayW/Gamma,3)
    BRZ = round(decayZ/Gamma,3)
    Gamma = round(Gamma,3)
    print(f'{mass}\t\t{Gamma}\t\t{BRH}\t{BRW}\t{BRZ}')
print()