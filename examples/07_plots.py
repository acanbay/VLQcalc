# Plots of the couplings, branching ratios and partial widths against the mass.
# Plots need matplotlib: pip install "VLQcalc[plot]"
import VLQcalc.model as model

vlq = model.VLQ('B')
vlq.setMass(range(1000, 2001, 50))
vlq.calcRatioKappas([0.25, 0.5, 0.25], 0.1)  # [BR(Hb), BR(Wt), BR(Zb)], Gamma/M

vlq.plot(save='kappas.pdf')
vlq.plot('BR', save='branchings.pdf')
vlq.plot('width', save='widths.pdf')
