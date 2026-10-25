# MadGraph5 cards

The `MG5` class writes a MadGraph5 input card that imports the model, generates the processes, launches the run and sets the masses and couplings of the VLQ object.

```python
import VLQcalc.model as model
import VLQcalc.madgraph as madgraph

vlq = model.VLQ('B', FNS=4, LR=False)
vlq.setMass([1000, 1500])
vlq.calcRatioKappas([0.25, 0.5, 0.25], 0.1)

mg5 = madgraph.MG5(vlq, 'VLQ_v4_UFO')
mg5.setProcess('p p > bp j, (bp > w- t, t > w+ b)')
mg5.shower = 'Pythia8'
mg5.addInput('nevents 10000')
mg5.createMG5Input('bp_single')
```

The card `bp_single.dat` is run with `./bin/mg5_aMC bp_single.dat`:

```
import model VLQ_v4_4FNS_UFO
define VLB = bp bp~
define WW = w+ w-
define tt = t t~
define bb = b b~
generate p p > VLB j, (VLB > WW tt, tt > WW bb)
output bp_single
launch bp_single
shower=Pythia8
detector=OFF
analysis=OFF
madspin=OFF
reweight=OFF
done
set MBP scan1:[1000, 1500]
set KBLh1 0
set KBRh1 0
set KBLh2 0
set KBRh2 0
set KBLh3 scan1:[1.6105, 1.59642]
set KBRh3 0
set KBLw1 0
set KBRw1 0
set KBLw2 0
set KBRw2 0
set KBLw3 scan1:[0.408354, 0.265449]
set KBRw3 0
set KBLz1 0
set KBRz1 0
set KBLz2 0
set KBRz2 0
set KBLz3 scan1:[0.390394, 0.260237]
set KBRz3 0
set nevents 10000
done
```

With more than one mass, every mass-dependent parameter is written as a `scan1` list, and MadGraph5 runs one point per mass. Numbers are written with 6 significant digits.

## Processes

`setProcess` defines the main process and `addProcess` additional ones. The particle names are replaced by multiparticle labels that contain the particle and the antiparticle, and the labels are defined once at the top of the card:

| Name in the process | Label | Name in the process | Label |
|-|-|-|-|
| `tp`, `tp~` | `VLT` | `w+`, `w-` | `WW` |
| `bp`, `bp~` | `VLB` | `e-`, `e+` | `ee` |
| `x`, `x~` | `VLX` | `mu-`, `mu+` | `mumu` |
| `y`, `y~` | `VLY` | `l-`, `l+` | `ll` |
| `t`, `t~` | `tt` | `ve`, `ve~` (also `vm`, `vt`, `vl`) | `veve` (`vmvm`, `vtvt`, `vlvl`) |
| `b`, `b~` | `bb` | | |

Only whole particle names are replaced, so `p p > x~ j / x` becomes `p p > VLX j / VLX`. Other parts of the process (`j`, `$$`, `QED=2`, ...) are written unchanged.

Because every name stands for the particle and the antiparticle, a process always includes its charge conjugate: `p p > tp b~ j` also generates $\bar{T} b j$, and `tp > w+ b` lets $\bar{T}$ decay as well. A process for one charge only has to be edited in the card.

## Run settings

| Attribute | Values | Default |
|-|-|-|
| `shower` | `OFF`, `Pythia8`, ... | `OFF` |
| `detector` | `OFF`, `Delphes`, ... | `OFF` |
| `analysis` | `OFF`, `ExRoot`, `MadAnalysis`, ... | `OFF` |
| `madspin` | `OFF`, `ON`, `onshell`, `full` | `OFF` |
| `reweight` | `OFF`, `ON` | `OFF` |

`addInput(text)` adds a line to the card settings (`set ` is added if missing), for example `mg5.addInput('ebeam1 7000')`. These lines are written after the masses and couplings, so they can also change a value written by VLQcalc, for example `mg5.addInput('KTRw3 0.2')`.

`createMG5Input(name)` writes `name.dat` (a name ending in `.dat` is used as it is). With a name other than the default `mg5_input`, the output folder of MadGraph5 gets the same name. The method can be called several times; each call writes the full card again. If the card cannot be completed (for example a VLQ_UFO card for T without branching ratios), a `ValueError` is raised and no file is written.

## What is written for each model

### VLQ_v4_UFO and VLQ_v5_UFO

- `import model VLQ_v4_4FNS_UFO` (or `5FNS`; `VLQ_v5_4FNS_NLO_UFO` for VLQ_v5).
- The mass (`MTP`, `MBP`, `MX` or `MY`).
- Every coupling of the VLQ to SM quarks (`vlq.Kappas`), including the ones that are zero.

The full model is imported, without a restriction card such as `-3rd`. The VLQ_v4 5FNS model has no restriction cards, and the restriction cards of the 4FNS model remove the right-handed couplings and the couplings to the 1st and 2nd generations, which VLQcalc can use. In the full VLQ_v4 model the default couplings to the 1st and 2nd generations are not zero (for example `KTLh1 = 0.479`), so VLQcalc writes all couplings explicitly. Writing them gives the same couplings as the restriction card would; MadGraph5 may take somewhat longer to generate processes with the full model, because it also generates the diagrams whose couplings are set to zero later.

If no couplings were computed, the 3rd generation couplings are not written and the defaults of the model are used; the couplings to lighter quarks are still set to zero. If no mass was set, the mass is not written either.

### VLQ_UFO

- `import model VLQ_UFO`.
- The generation fractions `zeta...L` and `zeta...R` from `vlq.zeta`. With `LR=True` each is $\zeta_i/2$; see [Models](models.md#both-chiralities).
- For T and B, the parameters `xi...` from `vlq.BRs`. A VLQ_UFO card for T or B needs them; `convertModel(..., reverse=True)` sets them together with `KappaOld` (see [Models](models.md#from-vlq_v4vlq_v5-to-vlq_ufo)).
- The mass and the overall coupling (`KT`, `KB`, `KX`, `KY`) from `vlq.KappaOld`.

`vlq.KappaOld` is set by `convertModel` or by hand. For a chosen $\Gamma/M$, compute the couplings with `calcRatioKappas` and convert them with `vlq.convertModel([vlq.KappaH, vlq.KappaW, vlq.KappaZ], reverse=True)` (`vlq.KappaW` for X and Y) before writing the card. Without `KappaOld`, no coupling line is written, the default of the model is used, and a warning is issued.

## Widths in MadGraph5

The card does not set the VLQ width. To compute it in MadGraph5 with the same couplings, add `mg5.addInput('WTP auto')` (`WBP`, `WX`, `WY` for the other VLQs). The width computed by MadGraph5 is the one VLQcalc computes; see [Physics conventions](physics.md#validation).
