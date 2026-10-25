# Tables and files of the results. Every table and file starts with the settings
# needed to reproduce it (VLQ, FNS, chirality, Gamma/M, branching ratios, zeta, SM inputs).
import VLQcalc.model as model

vlq = model.VLQ('B')
vlq.setMass(range(1000, 2001, 250))
vlq.calcRatioKappas([0.25, 0.5, 0.25], 0.1)  # [BR(Hb), BR(Wt), BR(Zb)], Gamma/M

print(vlq.table())                           # Markdown table
print(vlq.table('txt'))                      # aligned text table

vlq.save('kappas.md')                        # Markdown
vlq.save('kappas.txt')                       # text, settings as comment lines
vlq.save('kappas.csv')                       # CSV, settings as comment lines
vlq.save('kappas.json')                      # settings, couplings with model names, widths
