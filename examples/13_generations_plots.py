# Plots with couplings to lighter quarks: one panel per boson (default),
# or all channels in one panel with combined=True (the colour gives the boson).
import VLQcalc.model as model

vlq = model.VLQ('B', LR=True)
vlq.setMass(range(1000, 3001, 100))
vlq.calcRatioKappas([0.25, 0.5, 0.25], 0.1, zeta=[0, 0.3, 0.7])  # [BR(Hb), BR(Wt), BR(Zb)]

vlq.plot(save='kappas_panels.pdf')
vlq.plot(combined=True, save='kappas_combined.pdf')
vlq.plot('BR', combined=True, save='branchings_combined.pdf')
vlq.plot('width', save='widths_panels.pdf')
vlq.save('kappas_generations.json')
