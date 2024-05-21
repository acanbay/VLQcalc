import VLQcalc.model as model

vlq = model.VLQ('bp')
vlq.setMass( range(1000,2001,100) )
vlq.calcRatioKappas([0.25,0.5,0.25],0.1)

output_file = open(f'./save_kappas.txt', 'w')
output_file.writelines( [f'VLQ type: BP\n',
                         f'Branchings:\n',
                         f'  BR(Hb)={vlq.BRs[0]}\n',
                         f'  BR(Wt)={vlq.BRs[1]}\n',
                         f'  BR(Zb)={vlq.BRs[2]}\n',
                         '_____________\n\n'])

for i in range(len(vlq.Mass)):
    output_file.writelines( [f'{vlq.Mass[i]} GeV\n',
                             f'\tkH = {round(vlq.KappaH[i],4)}\n',
                             f'\tkW = {round(vlq.KappaW[i],4)}\n',
                             f'\tkZ = {round(vlq.KappaZ[i],4)}\n\n'])
output_file.close()