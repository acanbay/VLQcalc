# VLQcalc documentation

| Page | Content |
|-|-|
| [API reference](api.md) | Classes, methods, arguments and attributes |
| [Physics conventions](physics.md) | Couplings in the Lagrangian, decay widths, coupling calculation, SM inputs, choice of settings |
| [Models](models.md) | Supported UFO models, their parameters, conversion between models |
| [MadGraph5 cards](madgraph.md) | Writing input cards, what is written to them |
| [Limitations and common questions](faq.md) | What VLQcalc does not do, and answers to frequent problems |

For a first reading, start with the [README](../README.md) for the installation and a short example, then read [Models](models.md) to choose a model and [Physics conventions](physics.md) for the meaning of the inputs. The [tutorial](../binder/tutorial.ipynb) runs every step in a notebook, and the [examples](../examples) folder has one script per task.

## Typical workflows

### Couplings for a given Γ/M and branching ratios

```python
import VLQcalc.model as model
import VLQcalc.madgraph as madgraph

vlq = model.VLQ('T', FNS=4, LR=False)
vlq.setMass(range(1000, 2001, 250))
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05)   # BR(Ht), BR(Wb), BR(Zt); Gamma/M
vlq.save('kappas.csv')                          # couplings and settings

mg5 = madgraph.MG5(vlq, 'VLQ_v4_UFO')
mg5.setProcess('p p > tp tp~')
mg5.addInput('WTP auto')                        # MadGraph5 computes the width from the couplings
mg5.createMG5Input('tp_pair')
```

The card `tp_pair.dat` runs one point per mass. Examples: `06_calculate_ratio_kappas.py`, `11_mg5_input_v4_4FNS.py`.

### Converting a coupling of VLQ_UFO

```python
import VLQcalc.model as model

vlq = model.VLQ('B')
vlq.setMass([1000, 1500, 2000])
vlq.convertModel(0.5, [0.25, 0.5, 0.25])        # kappa_B and the branching ratios xi of VLQ_UFO
print(vlq.KappaH, vlq.KappaW, vlq.KappaZ)      # couplings of VLQ_v4/VLQ_v5
```

`convertModel(..., reverse=True)` converts in the other direction. Examples: `01_converting_models_1.py`, `02_converting_models_2.py`, `08_convert_ratio_kappas.py`.

### Couplings to 1st and 2nd generation quarks

```python
import VLQcalc.model as model

vlq = model.VLQ('T')
vlq.setMass([1000, 1500])
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05, zeta=[0.1, 0.2, 0.7])
vlq.Kappas['KTLw2']                             # coupling to the s quark
```

The generation fractions ζ<sub>i</sub> give the share of the decays into each generation. Example: `15_generations.py`.

### SM inputs different from the defaults

If the param_card of MadGraph5 is changed, give the same values to VLQcalc before computing the couplings:

```python
import VLQcalc.constants as constants
constants.setSM(Mt=173.0)
```

See [SM inputs](physics.md#sm-inputs).
