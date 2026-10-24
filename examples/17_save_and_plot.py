import VLQcalc.model as model

vlq = model.VLQ('B', LR=True)
vlq.setMass( range(1000, 3001, 100) )
vlq.calcRatioKappas([0.25, 0.5, 0.25], 0.1, zeta=[0, 0.3, 0.7])

vlq.save('kappas.md')    # Markdown table
vlq.save('kappas.csv')   # CSV, settings as comment lines
vlq.save('kappas.json')  # settings, couplings with model names, widths

# plots need matplotlib: pip install "VLQcalc[plot]"
vlq.plot(save='kappas.pdf')
vlq.plot('BR', save='branchings.pdf')
vlq.plot('width', save='widths.pdf')
