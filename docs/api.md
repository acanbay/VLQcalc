# API reference

Masses and widths are in GeV; couplings and branching ratios are dimensionless. Wrong arguments raise `ValueError` with a message that names the method and the expected input.

## `VLQcalc.model.VLQ`

```python
import VLQcalc.model as model
vlq = model.VLQ(VLQ_type, FNS=4, LR=False)
```

| Argument | Type | Values | Default |
|-|-|-|-|
| `VLQ_type` | str | `'T'`, `'B'`, `'X'`, `'Y'` (also `'tp'`, `'bp'`, `'VLT'`, ...; case-insensitive) | – |
| `FNS` | int | `4`, `5` (in 5FNS $m_b = 0$) | `4` |
| `LR` | bool | `False`: left-handed couplings only; `True`: left- and right-handed couplings with $\kappa_R = \kappa_L$ | `False` |

### `setMass(Masses)`

Sets the mass in GeV: a number, or a list, tuple, range or numpy array of numbers. All results are given for each mass.

```python
vlq.setMass(1000)
vlq.setMass(range(1000, 2001, 100))
```

### `calcRatioKappas(BRs, Ratio, zeta=None)`

Computes the couplings that give the width-to-mass ratio `Ratio` and the branching ratios `BRs` for each mass. The formula is given in [Physics conventions](physics.md#couplings-from-γm-and-branching-ratios).

| Argument | T, B | X, Y |
|-|-|-|
| `BRs` | T: `[BR(Ht), BR(Wb), BR(Zt)]`<br>B: `[BR(Hb), BR(Wt), BR(Zb)]`<br>sum 1 | not used (`None`) |
| `Ratio` | $\Gamma/M$ | $\Gamma/M$ |
| `zeta` | generation fractions `[ζ1, ζ2, ζ3]`, sum 1; default `[0, 0, 1]` | same |

```python
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05)                        # T or B
vlq.calcRatioKappas(Ratio=0.05)                                      # X or Y
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05, zeta=[0.1, 0.2, 0.7])  # with lighter quarks
```

### `convertModel(Kappas, BRs=1, reverse=False, zeta=None)`

Converts couplings between VLQ_UFO and VLQ_v4/VLQ_v5; see [Models](models.md#conversion-between-vlq_ufo-and-vlq_v4vlq_v5).

- `reverse=False` (VLQ_UFO → VLQ_v4/VLQ_v5): `Kappas` is $\kappa_Q$ (a number), `BRs` is `[ξH, ξW, ξZ]` for T and B, and `zeta` gives the generation fractions. The results are stored in `KappaH`, `KappaW`, `KappaZ` and `KappaLight`; `KappaOld` is set to `[Kappas]`.
- `reverse=True` (VLQ_v4/VLQ_v5 → VLQ_UFO): `Kappas` is `[kH, kW, kZ]` for T and B, or `kW` for X and Y, with one value per mass (a number is accepted for a single mass). Couplings to lighter quarks stored in the object are included. The result is stored in `KappaOld`.

### `calcDecay(Mass, Kappas, LR)`

Partial widths for the 3rd generation:

```python
decayH, decayW, decayZ, Gamma = vlq.calcDecay(1500, [kH, kW, kZ], False)   # T, B
decayW = vlq.calcDecay(1500, kW, False)                                     # X, Y
```

### `calcWidths(index)`

Partial widths of all channels for the mass `vlq.Mass[index]`, with the couplings stored in the object, as a dictionary `{(boson, generation): width}`:

```python
widths = vlq.calcWidths(0)
widths[('w', 3)]          # T -> W b
sum(widths.values())      # total width
```

### `table(style='md')`

Returns the results as a Markdown table (`style='md'`) or an aligned text table (`style='txt'`). The table starts with the settings needed to reproduce it (VLQ, FNS, chirality, $\Gamma/M$, branching ratios, $\zeta$, SM inputs) and has one row per mass with the couplings, the total width and, if present, $\kappa_Q$ of VLQ_UFO.

### `save(filename)`

Writes the results to a file; the format follows the extension:

| Extension | Content |
|-|-|
| `.md` | Markdown table, as `table()` |
| `.txt` | aligned text table, settings as `#` comment lines |
| `.csv` | comma-separated table, settings as `#` comment lines |
| `.json` | settings, masses, all non-zero couplings with the model names, widths, `KappaOld` |

A JSON file for a T quark with two masses (the widths are rounded here; the file has full precision):

```json
{
  "VLQ": "T (charge 2/3)",
  "FNS": 4,
  "couplings": "left-handed",
  "Gamma/M": 0.05,
  "branchings": {"BR(Ht)": 0.5, "BR(Wb)": 0.25, "BR(Zt)": 0.25},
  "zeta": [0, 0, 1],
  "SM inputs": {"MZ": 91.1876, "aEWM1": 127.9, "Gf": 1.16637e-05, "Mt": 172.0, "Mb": 4.7, "MH": 125.0, "MW": 79.8244},
  "mass": [1000, 1500],
  "kappas": {"KTLh3": [1.61204, 1.59671], "KTLw3": [0.195189, 0.130117], "KTLz3": [0.288764, 0.187702]},
  "width": [50.000005, 74.999937]
}
```

`kappa VLQ_UFO` is added after a conversion. The widths are total widths in GeV, computed from the rounded couplings, so they agree with Γ/M × M to about 6 digits.

### `plot(quantity='kappa', save=None, ax=None, logy=None)`

Plots the couplings (`'kappa'`), branching ratios (`'BR'`) or partial widths (`'width'`) against the mass and returns the matplotlib axes. With couplings to lighter quarks there is one panel per boson, and the generations are told apart by the line style. `save` writes the figure (for example `'kappas.pdf'`); `ax` draws into existing axes (single panel only); `logy` forces a logarithmic (`True`) or linear (`False`) axis, otherwise the scale is chosen from the range of the values. Needs matplotlib (`pip install "VLQcalc[plot]"`).

### Attributes

| Attribute | Content |
|-|-|
| `type` | `'VLT'`, `'VLB'`, `'VLX'` or `'VLY'` |
| `FNS`, `LR` | as given |
| `Mass` | list of masses |
| `BRs`, `Ratio`, `zeta` | inputs of the last calculation |
| `KappaH`, `KappaW`, `KappaZ` | couplings to 3rd generation quarks, one value per mass (`KappaH`, `KappaZ` only for T and B) |
| `KappaLight` | couplings to 1st and 2nd generation quarks, `{(boson, generation): [values]}` |
| `Kappas` | all couplings with the names of VLQ_v4/VLQ_v5 (read only), for example `vlq.Kappas['KTLw2']` |
| `KappaOld` | $\kappa_Q$ of VLQ_UFO |
| `bosons` | bosons of the decay channels, `['h', 'w', 'z']` or `['w']` |

The couplings can also be set by hand, for example `vlq.KappaOld = [0.5]` before writing a VLQ_UFO card.

## `VLQcalc.madgraph.MG5`

```python
import VLQcalc.madgraph as madgraph
mg5 = madgraph.MG5(vlq, model='VLQ_v4_UFO')
```

| Method or attribute | Content |
|-|-|
| `model` argument | `'VLQ_UFO'`, `'VLQ_v4_UFO'` or `'VLQ_v5_UFO'` |
| `setProcess(process)` | main process, for example `'p p > tp tp~'` |
| `addProcess(process)` | additional process (after `setProcess`) |
| `addInput(text)` | a setting for the run cards, for example `'nevents 10000'` |
| `shower`, `detector`, `analysis`, `madspin`, `reweight` | run settings (default `'OFF'`) |
| `createMG5Input(name='mg5_input')` | writes `name.dat` |

Details are in [MadGraph5 cards](madgraph.md).

## `VLQcalc.constants`

| Function | Content |
|-|-|
| `setSM(**inputs)` | changes SM inputs (`MZ`, `aEWM1`, `Gf`, `Mt`, `Mb`, `MH`); without arguments restores the defaults |
| `getSM()` | current inputs and the derived `MW` |

## `VLQcalc.decays`

| Function | Content |
|-|-|
| `partialWidths(VLQtype, FNS, M, Kappas, LR=False)` | partial widths `{(boson, generation): width}` for couplings `{(boson, generation): kappa}` |
| `decayTP4`, `decayBP4`, `decayX4`, `decayY4`, `decayTP5`, ... | 3rd generation widths, as `VLQ.calcDecay` |
