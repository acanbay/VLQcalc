import VLQcalc.model as model
import matplotlib.pyplot as plt

vlq = model.VLQ('B')
vlq.setMass( range(1000, 2001, 100) )
vlq.calcRatioKappas([0.25, 0.5, 0.25], 0.1)

plt.plot(vlq.Mass, vlq.KappaW, '-b', label='$\kappa_{W}$')
plt.plot(vlq.Mass, vlq.KappaH, '-.g', label='$\kappa_{H}$')
plt.plot(vlq.Mass, vlq.KappaZ, '--r', label='$\kappa_{Z}$')
plt.legend()
plt.xlabel('$\mathbf{m_{B}}$ [GeV]', horizontalalignment='right', x=1.0, weight='bold')
plt.ylabel('$\mathbf{\kappa}$', horizontalalignment='right', y=1.0, weight='bold')
plt.xlim([min(vlq.Mass),max(vlq.Mass)])
plt.grid()
# plt.savefig('kappa_plot.png',dpi=720)
plt.show()