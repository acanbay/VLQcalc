# Width-to-mass ratio and branching ratios that a VLQ_UFO coupling corresponds to.
import VLQcalc.model as model

vlq = model.VLQ('T')
vlq.setMass([1000, 1500, 2000])
vlq.convertModel(0.5, [0.5, 0.25, 0.25])     # kappa_T and [xi(H), xi(W), xi(Z)] of VLQ_UFO

print('M [GeV]   Gamma/M   BR(Ht)   BR(Wb)   BR(Zt)')
for i, M in enumerate(vlq.Mass):
    widths = vlq.calcWidths(i)
    total = sum(widths.values())
    BRs = [widths[(boson, 3)] / total for boson in ['h', 'w', 'z']]
    print(f'{M:7d}   {total/M:7.4f}   {BRs[0]:6.4f}   {BRs[1]:6.4f}   {BRs[2]:6.4f}')
