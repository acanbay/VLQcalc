import VLQcalc.model as model

Mass = range(1000,2001,200)
branchings = [0.25,0.5,0.25] # [BR(Hb),BR(Wt),BR(Zb)]
GM = 0.1 # Gamma/Mass ratio

vlq = model.VLQ('tp')
vlq.setMass( Mass )
vlq.calcRatioKappas(branchings,GM)

newKappas = [vlq.KappaH,vlq.KappaW,vlq.KappaZ]
vlq.convertModel(newKappas,branchings,reverse=True)  # converts kappas VLQ_v4 to VLQ_UFO

print('MBP [GeV]\tkB')
for i in range(len(Mass)):
    M = Mass[i]
    kB = round(vlq.KappaOld[i],2)
    print(f'  {M}\t\t{kB}')
print()