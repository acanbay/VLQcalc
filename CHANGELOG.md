# Changelog

## 0.2.0 (2026-10-25)

Scripts written for 0.1.x keep working without changes. Some numerical results change; the reasons are listed under *Results that change*.

### Results that change

- **SM inputs.** The widths now use the SM inputs of the VLQ_v4/VLQ_v5 restriction cards, and $m_W$ is derived from $G_F$ as in the models ($m_W = 79.824$ GeV instead of 80.379 GeV; $m_t = 172$ GeV, $m_b = 4.7$ GeV, $m_H = 125$ GeV). The widths computed by VLQcalc are now identical to the ones MadGraph5 computes with these models. With the old inputs, the W and Z widths were 3.3% larger, so a coupling computed for a given $\Gamma/M$ gave a different $\Gamma/M$ in MadGraph5. For a given $\Gamma/M$ and branching ratios, $\kappa_W$ and $\kappa_Z$ change by 1.7% and $\kappa_H$ by less than 0.1%.
- **VLQ_UFO conversion.** The Higgs vacuum expectation value and $m_W$ are now derived from the VLQ_UFO inputs ($v = 246.26$ GeV instead of 247.3 GeV). This changes the converted $\kappa_H$ by 0.4%; the converted couplings now reproduce the VLQ_UFO vertices exactly.
- **Couplings from $\Gamma/M$.** The genetic algorithm is replaced by an exact solution. Every partial width is proportional to $\kappa^2$, so the coupling of each channel follows directly from its branching ratio and $\Gamma/M$. The results agree with the targets to the stated precision (6 significant digits) and no longer change from run to run.
- **VLQ_UFO cards with `LR=True`.** Both chiralities now get half of the generation fraction ($\zeta_L = \zeta_R = 0.5$ for the 3rd generation), so that the sum of all $\zeta$ stays 1 as in the model. Before, $\zeta_R$ was set to 0.1 without being used in the conversion.

### Fixed

- MadGraph5 cards for VLQ_v4_UFO: the model name was written as `VLQ_v4__4FNS_UFO` (two underscores), the parameters of VLQ_UFO were written instead of the couplings of VLQ_v4, and the couplings were never written.
- MadGraph5 process text: particle names are now matched as whole names. Before, `tp~` became `VLT~`, a `b` at the end of a process was not replaced, any `x` or `y` inside other names was replaced, and a veto (`/ x`) without a decay chain was dropped.
- Calling `createMG5Input` twice wrote the card twice into the file.
- VLQ_v4 5FNS: the official model has no `restrict_3rd` card, so `import model VLQ_v4_5FNS_UFO-3rd` failed. The full model is now imported for both flavour schemes, and every coupling of the VLQ is written to the card (see [MadGraph5 cards](docs/madgraph.md)). This also makes `LR=True` work: the `-3rd` restriction removed the right-handed couplings.
- `convertModel(..., reverse=True)` followed by a VLQ_UFO card failed, because the branchings were not defined; they are now taken from the widths.
- A VLQ_UFO card after a forward conversion failed, because `KappaOld` was a number instead of a list.
- `convertModel(..., reverse=True)` crashed with `TypeError` for some wrong inputs instead of reporting the error.
- Branchings such as `[0.1, 0.2, 0.7]` were rejected because their floating-point sum is not exactly 1.
- `setMass` now accepts tuples and numpy arrays and rejects non-positive masses.
- The labels of the X and Y branchings were swapped (X decays to $Wt$, Y to $Wb$).
- `examples/14_plot_kappas.py` and the tutorial used invalid escape sequences in plot labels.

### Added

- Couplings to 1st and 2nd generation quarks with the generation fractions `zeta=[ζ1, ζ2, ζ3]` in `calcRatioKappas` and `convertModel`; the default `[0, 0, 1]` keeps the 3rd generation only.
- `vlq.Kappas`: every coupling of the VLQ, named as in the VLQ_v4/VLQ_v5 models (`KTLh1`, ..., `KTRz3`).
- `vlq.calcWidths(i)`: partial widths of all channels for the i-th mass.
- Support for VLQ_v5_UFO in `MG5`.
- `VLQcalc.constants` with `setSM(...)` and `getSM()` to change the SM inputs, for example when the param_card of MadGraph5 is changed.
- `vlq.table()`, `vlq.save(filename)` and `vlq.plot(...)` for tables, files (Markdown, CSV, JSON, text) and plots.
- Installation with pip (`pyproject.toml`) or `python setup.py install`, and releases on PyPI.
- A test suite and documentation in `docs/`.

### Changed

- Wrong inputs raise `ValueError` with the same messages as before, instead of printing the message and stopping Python with `sys.exit`. Scripts stop in the same way, and notebooks keep running.
- Warnings about replaced couplings are issued with the `warnings` module.
- The Binder environment uses Python 3.12.
