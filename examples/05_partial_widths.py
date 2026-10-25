# Partial width and branching ratio of every channel, with the couplings stored in the object.
import VLQcalc.model as model

vlq = model.VLQ('T')
vlq.setMass([1000, 1500, 2000])
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05)

for i, M in enumerate(vlq.Mass):
    widths = vlq.calcWidths(i)               # {(boson, generation): width in GeV}
    total = sum(widths.values())
    print(f'\nM = {M} GeV, total width = {total:.4f} GeV')
    print('channel        width [GeV]   BR')
    for (boson, generation), width in widths.items():
        if width > 0:
            print(f'{boson.upper()}, gen. {generation}      {width:10.4f}   {width/total:.4f}')
