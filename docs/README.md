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

### Couplings for a given width-to-mass ratio and branching ratios

```python
import VLQcalc.model as model
import VLQcalc.madgraph as madgraph

vlq = model.VLQ('T', FNS=4, LR=False)
vlq.setMass(range(1000, 2001, 250))
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05)   # BR(Ht), BR(Wb), BR(Zt); Gamma/M
vlq.save('kappas.csv')                          # couplings and settings

mg5 = madgraph.MG5(vlq, 'VLQ_v4_UFO')
mg5.setProcess('p p > tp tp~')
mg5.createMG5Input('tp_pair')
```

The card `tp_pair.dat` runs one point per mass. See examples 01-03 and 08-10.

### Converting a coupling of VLQ_UFO

```python
import VLQcalc.model as model

vlq = model.VLQ('B')
vlq.setMass([1000, 1500, 2000])
vlq.convertModel(0.5, [0.25, 0.5, 0.25])        # kappa_B and the branching ratios xi of VLQ_UFO
print(vlq.KappaH, vlq.KappaW, vlq.KappaZ)      # couplings of VLQ_v4/VLQ_v5
```

`convertModel(..., reverse=True)` converts in the other direction. See examples 14-18.

### Couplings to 1st and 2nd generation quarks

```python
import VLQcalc.model as model

vlq = model.VLQ('T')
vlq.setMass([1000, 1500])
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05, zeta=[0.1, 0.2, 0.7])
vlq.Kappas['KTLw2']                             # coupling to the s quark
```

The generation fractions $\zeta_i$ give the share of the decays into each generation. `vlq.calcWidths(i)` gives the partial width of every channel for the i-th mass (example 05). See examples 11-13 and 18.

### SM inputs different from the defaults

If the param_card of MadGraph5 is changed, give the same values to VLQcalc before computing the couplings:

```python
import VLQcalc.constants as constants
constants.setSM(Mt=173.0)
```

See [SM inputs](physics.md#sm-inputs) and example 21.

## Examples

The [examples](../examples) folder has one script per task, ordered from the basic calculation to the conversion between models:

| Script | Content |
|-|-|
| `01_couplings_from_width.py` | Couplings of T for a chosen $\Gamma/M$ and branching ratios |
| `02_couplings_from_width_LR.py` | The same with left- and right-handed couplings |
| `03_couplings_X_Y.py` | Couplings of X and Y, which decay only through W |
| `04_decay_widths.py` | Width and branching ratios for given couplings |
| `05_partial_widths.py` | Partial width of every channel |
| `06_tables_and_files.py` | Tables and files (Markdown, text, CSV, JSON) |
| `07_plots.py` | Plots of couplings, branching ratios and widths |
| `08_mg5_card_v4_4FNS.py` | MadGraph5 card, VLQ_v4_UFO, 4-flavour scheme |
| `09_mg5_card_v4_5FNS.py` | MadGraph5 card, VLQ_v4_UFO, 5-flavour scheme |
| `10_mg5_card_v5.py` | MadGraph5 card, VLQ_v5_UFO |
| `11_generations.py` | Couplings to 1st and 2nd generation quarks |
| `12_generations_mg5_card.py` | MadGraph5 card with couplings to lighter quarks |
| `13_generations_plots.py` | Plots with couplings to lighter quarks, one panel per boson or combined |
| `14_convert_from_VLQ_UFO.py` | Conversion of a VLQ_UFO coupling (T, B) |
| `15_convert_from_VLQ_UFO_X_Y.py` | Conversion of a VLQ_UFO coupling (X, Y) |
| `16_convert_to_VLQ_UFO.py` | Conversion to VLQ_UFO |
| `17_converted_decay_widths.py` | $\Gamma/M$ and branching ratios of a VLQ_UFO coupling |
| `18_convert_generations.py` | Conversion with generation fractions |
| `19_mg5_card_VLQ_UFO.py` | MadGraph5 card, VLQ_UFO, coupling set by hand |
| `20_mg5_card_VLQ_UFO_converted.py` | MadGraph5 card, VLQ_UFO, coupling for a chosen $\Gamma/M$ |
| `21_sm_inputs.py` | SM inputs different from the defaults |
