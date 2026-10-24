import VLQcalc.model as model

vlq = model.VLQ('T')
vlq.setMass( [1000, 1500, 2000] )
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05, zeta=[0, 0.2, 0.8])

# partial widths of all channels for each mass, with the couplings stored in the object
for i, M in enumerate(vlq.Mass):
    widths = vlq.calcWidths(i)
    total = sum(widths.values())
    print(f'\nM = {M} GeV, total width = {total:.4f} GeV')
    print('channel\t\twidth [GeV]\tBR')
    for (boson, generation), width in widths.items():
        if width > 0:
            print(f'{boson.upper()} + gen. {generation}\t{width:.4f}\t\t{width/total:.4f}')
