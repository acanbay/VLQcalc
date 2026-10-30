# Changelog

## 0.2.1 (2026-10-30)

Scripts written for 0.1.x and 0.2.0 keep working. MadGraph5 cards now set the VLQ width.

### Changed

- **The VLQ width in MadGraph5 cards.** Cards now set the width of the VLQ to `auto` (`set WTP auto`, `WBP`, `WX` or `WY`), so MadGraph5 computes it from the couplings written to the card; with VLQ_v4 and VLQ_v5 this is the width VLQcalc computes. Before, no width was written unless `mg5.addInput('WTP auto')` was added, and MadGraph5 used the fixed width of the param_card (5 GeV in VLQ_v4, 11 GeV in VLQ_v5), which does not belong to the computed couplings and changes the branching ratios in decay chains. `mg5.width` sets a width in GeV (one value per mass for a scan), or `None` to keep the value of the param_card. Scripts that already add `WTP auto` write it twice, which gives the same result.
- **SM inputs.** The VLQ object keeps the SM inputs its couplings were computed with. Tables, files and `calcWidths` use them, also after a later `constants.setSM(...)`. Before, a later `setSM` changed the widths and the SM inputs shown next to couplings computed with other values.
- `vlq.calcDecay(M, Kappas)` uses the `LR` setting of the object when `LR` is not given; before, it always used `LR=False`.

### Fixed

- After `convertModel(..., reverse=True)` for one mass, changing the mass kept the VLQ_UFO coupling `KappaOld` of the old mass, and it was written to VLQ_UFO cards. It is now deleted with the other couplings of the old mass. A `KappaOld` set by hand is still kept, also as a single number (before, a number was deleted and a list with one value was kept).
- After `convertModel(..., reverse=True)`, tables and saved files still showed the $\Gamma/M$ of an earlier `calcRatioKappas`, although the converted couplings give a different width.
- The $\xi_V$ and $\zeta_i$ written by the reverse conversion are rounded to 6 significant digits, so their sum can differ from 1 by about $10^{-7}$; given back to `convertModel` or `calcRatioKappas`, they were rejected. Sums within $10^{-5}$ of 1 are now accepted.
- `convertModel` did not accept numpy arrays for the couplings of T and B (`reverse=True`) or for a single $\kappa_Q$ (`np.array([0.3])`); they are now accepted like lists.
- With couplings to lighter quarks (`zeta`), tables and saved files labelled the branching ratios per boson with the names of 3rd generation channels (`BR(Ht)`, `BR(Wt)`); they are now labelled `BR(H)`, `BR(W)`, `BR(Z)`.

### Added

- `mg5.width = 'vlqcalc'` writes the widths VLQcalc computes (VLQ_v4/VLQ_v5 cards), one value per mass in the order of the masses, so that each point of a mass scan gets its own width. A width given as a number is checked: it must be a non-negative number, or one per mass.
- Warnings:
  - `convertModel` from VLQ_UFO warns when a channel that is open in VLQ_UFO, which neglects the quark masses, is kinematically closed in VLQ_v4/VLQ_v5 (for example $X \to Wt$ at 200 GeV); the coupling of this channel gives no width there.
  - VLQ_UFO cards with `FNS=5` warn that VLQ_UFO exists only in the 4-flavour scheme (massive b quark), so MadGraph5 runs in the 4-flavour scheme; the card itself does not depend on `FNS`.
  - VLQ_v4/VLQ_v5 cards warn when no couplings to the 3rd generation are set, so that the non-zero defaults of the model are used.
- Documentation: the widths of VLQ_UFO and VLQ_v4 for converted couplings differ slightly because of the quark masses of their param_cards; the increase of the T width with `LR=True` in the conversion depends on $\xi_H$; $\zeta_i$ of VLQcalc and of VLQ_UFO agree for massless quarks; the reverse conversion does not keep the signs of the couplings.

## 0.2.0 (2026-10-25)

Scripts written for 0.1.x keep working. Some numerical results change, for the reasons listed under *Results that change*, and a few behaviours change, as listed under *Changed*.

### Results that change

- **SM inputs.** The widths now use the default SM inputs of the VLQ_v4/VLQ_v5 models, and $m_W$ is derived from $G_F$ as in the models ($m_W = 79.824$ GeV instead of 80.379 GeV, $m_t = 172$ GeV instead of 172.9 GeV, $m_b = 4.7$ GeV instead of 4.18 GeV, $m_H = 125$ GeV instead of 125.1 GeV). The widths computed by VLQcalc are now identical to the ones MadGraph5 computes with these models. With the old inputs, the W and Z widths were 3.3% larger, so a coupling computed for a given $\Gamma/M$ gave a different $\Gamma/M$ in MadGraph5. For a given $\Gamma/M$ and branching ratios, $\kappa_W$ and $\kappa_Z$ are 1.6% smaller, and $\kappa_H$ changes by less than 0.1% above about 500 GeV (more close to the $Ht$ threshold).
- **VLQ_UFO conversion.** The Higgs vacuum expectation value and $m_W$ are now derived from the VLQ_UFO inputs ($v = 246.26$ GeV instead of 247.3 GeV). This changes the converted $\kappa_H$ by 0.4%; the converted couplings now reproduce the VLQ_UFO vertices exactly.
- **Couplings from $\Gamma/M$.** The genetic algorithm is replaced by an exact solution. Every partial width is proportional to $\kappa^2$, so the coupling of each channel follows directly from its branching ratio and $\Gamma/M$. The results agree with the targets to the stated precision (6 significant digits) and no longer change from run to run.
- **Conversion from VLQ_v4 to VLQ_UFO.** The reverse conversion is now the exact inverse of the vertex matching, so VLQ_UFO $\to$ VLQ_v4 $\to$ VLQ_UFO returns the same $\kappa_Q$, $\xi_V$ and $\zeta_i$. Before, every channel was inverted with the branching ratios of VLQ_v4, which include the quark masses, and the results were averaged; for T with `LR=True` the round trip differed by 11% at 500 GeV and 3.4% at 1 TeV. The $\xi_V$ and $\zeta_i$ written to VLQ_UFO cards are now the vertex parameters that belong to the written $\kappa_Q$.
- **VLQ_UFO cards with `LR=True`.** Both chiralities now get half of the generation fraction ($\zeta_L = \zeta_R = 0.5$ for the 3rd generation), so that the sum of all $\zeta$ stays 1 as in the model. Before, $\zeta_R$ was set to 0.1 without being used in the conversion.

### Fixed

- MadGraph5 cards for VLQ_v4_UFO: the model name was written as `VLQ_v4__4FNS_UFO` (two underscores), the parameters of VLQ_UFO were written instead of the couplings of VLQ_v4, and the couplings were never written.
- MadGraph5 process text: particle names are now matched as whole names. Before, `tp~` became `VLT~`, a `b` at the end of a process was not replaced, any `x` or `y` inside other names was replaced, and a veto (`/ x`) without a decay chain was dropped.
- Calling `createMG5Input` twice wrote the card twice into the file.
- VLQ_v4 5FNS: the official model has no `restrict_3rd` card, so `import model VLQ_v4_5FNS_UFO-3rd` failed. The full model is now imported for both flavour schemes, and every coupling of the VLQ is written to the card (see [MadGraph5 cards](docs/madgraph.md)). This also makes `LR=True` work: the `-3rd` restriction removed the right-handed couplings.
- `convertModel(..., reverse=True)` followed by a VLQ_UFO card failed, because the branchings were not defined; they are now set by the conversion.
- A VLQ_UFO card after a forward conversion failed, because `KappaOld` was a number instead of a list.
- `convertModel(..., reverse=True)` crashed with `TypeError` for some wrong inputs instead of reporting the error.
- Branchings such as `[0.1, 0.2, 0.7]` were rejected because their floating-point sum is not exactly 1.
- `setMass` now accepts tuples and numpy arrays and rejects non-positive masses.
- The labels of the X and Y branchings were swapped (X decays to $Wt$, Y to $Wb$).
- An example script and the tutorial used invalid escape sequences in plot labels.

### Added

- Couplings to 1st and 2nd generation quarks with the generation fractions `zeta=[zeta1, zeta2, zeta3]` in `calcRatioKappas` and `convertModel`; the default `[0, 0, 1]` keeps the 3rd generation only.
- `vlq.Kappas`: every coupling of the VLQ, named as in the VLQ_v4/VLQ_v5 models (`KTLh1`, ..., `KTRz3`).
- `vlq.calcWidths(i)`: partial widths of all channels for the i-th mass.
- Support for VLQ_v5_UFO in `MG5`.
- `VLQcalc.constants` with `setSM(...)` and `getSM()` to change the SM inputs, for example when the param_card of MadGraph5 is changed.
- `vlq.table()`, `vlq.save(filename)` and `vlq.plot(...)` for tables, files (Markdown, CSV, JSON, text) and black-and-white plots; `vlq.plot(..., combined=True)` draws all generations in one panel.
- Installation with pip (`pyproject.toml`) or `python setup.py install`, and releases on PyPI.
- A test suite and documentation in `docs/`.
- The examples are rewritten and ordered from the basic calculation to the conversion between models; new examples cover X and Y, partial widths, couplings to lighter quarks and SM inputs.

### Changed

- The license is MIT instead of GPL-3.0. Versions up to 0.1.2 remain under GPL-3.0.
- Wrong inputs raise `ValueError` with a message that names the method, instead of printing the message and stopping Python with `sys.exit`. Scripts stop in the same way, and notebooks keep running.
- After a conversion from VLQ_UFO, `vlq.KappaOld` is a list with one value (`[0.5]` instead of `0.5`), like after the reverse conversion; a number set by hand is still accepted.
- `setMass` with new masses deletes the couplings computed for the previous masses, with a warning. Before, the old couplings were kept and written for the new masses.
- Channels below their threshold have zero width; a branching ratio for such a channel raises `ValueError`.
- Settings given with `MG5.addInput` are written after the masses and couplings, so that they can change them.
- Warnings about replaced couplings are issued with the `warnings` module.
- The Binder environment uses Python 3.12.
