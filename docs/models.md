# Models

VLQcalc works with three FeynRules models of vector-like quarks. The models are not part of VLQcalc; they have to be installed in MadGraph5 separately.

| Model in VLQcalc | UFO folder imported in MadGraph5 | Source |
|-|-|-|
| `VLQ_UFO` | `VLQ_UFO` | [FeynRules: VLQ](https://feynrules.irmp.ucl.ac.be/wiki/VLQ) [1] |
| `VLQ_v4_UFO` | `VLQ_v4_4FNS_UFO`, `VLQ_v4_5FNS_UFO` | [FeynRules: NLOModels](https://feynrules.irmp.ucl.ac.be/wiki/NLOModels) [2] |
| `VLQ_v5_UFO` | `VLQ_v5_4FNS_NLO_UFO`, `VLQ_v5_5FNS_NLO_UFO` | [FeynRules: NLOModels](https://feynrules.irmp.ucl.ac.be/wiki/NLOModels) [2] |

The flavour scheme of the folder is chosen with the `FNS` argument of the VLQ object.

## VLQ_v4_UFO and VLQ_v5_UFO

Every coupling of a VLQ to a SM quark is a separate parameter:

| VLQ | Parameters | Example |
|-|-|-|
| T | `KT` + chirality (`L`, `R`) + boson (`h`, `w`, `z`) + generation (1, 2, 3) | `KTLw3`: $T$–$W$–$b$, left-handed |
| B | `KB` + chirality + boson + generation | `KBRh1`: $B$–$H$–$d$, right-handed |
| X | `KX` + chirality + generation | `KXL3`: $X$–$W$–$t$, left-handed |
| Y | `KY` + chirality + generation | `KYL2`: $Y$–$W$–$s$, left-handed |

The masses are `MTP`, `MBP`, `MX` and `MY`. The decay formulas and SM inputs are the same in both models. VLQ_v5 adds the full CKM matrix and couplings between two VLQs (`KWTBL`, `KZTTL`, ...); these couplings are zero by default and VLQcalc does not change them.

In VLQcalc, `vlq.Kappas` holds all couplings of the VLQ with these names, one value per mass:

```python
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05, zeta=[0, 0.2, 0.8])
vlq.Kappas['KTLw3']   # same as vlq.KappaW
vlq.Kappas['KTLw2']   # coupling to the s quark
vlq.Kappas['KTRw3']   # zeros, because LR=False
```

`vlq.KappaH`, `vlq.KappaW` and `vlq.KappaZ` hold the couplings to 3rd generation quarks, as in version 0.1.x. The couplings to lighter quarks are in `vlq.KappaLight`, a dictionary with keys `(boson, generation)`.

## VLQ_UFO

VLQ_UFO describes a VLQ with an overall coupling $\kappa_Q$ (`KT`, `KB`, `KX`, `KY`), the branching ratios of the bosons $\xi_V$ (`xitph`, `xitpw`, `xitpz` for T; `xibph`, `xibpw`, `xibpz` for B) and the generation fractions $\zeta_i$ (`zetaTuL`, `zetaTcL`, `zetaTtL` and the right-handed `zeta...R` for T, and so on).

## Conversion between VLQ_UFO and VLQ_v4/VLQ_v5

`convertModel` matches the vertices of both models. From VLQ_UFO to VLQ_v4/VLQ_v5:

$$\kappa_{V,i} = \kappa_Q \sqrt{s\,\zeta_i\,\xi_V}\; F_V(M),$$

with

$$F_H = \frac{M}{v}\sqrt{\frac{1}{\Gamma^0_H}}, \qquad F_W = \sqrt{\frac{1}{\Gamma^0_W}}\;\frac{g}{g'}, \qquad F_Z = \sqrt{\frac{1}{\Gamma^0_Z}}\;\frac{g/c_W}{g'/c'_W},$$

$$\Gamma^0_H = \frac{(1 - x_H)^2}{2}, \qquad \Gamma^0_W = (1 - x_W)(1 + x_W - 2x_W^2), \qquad \Gamma^0_Z = \frac{(1 - x_Z)(1 + x_Z - 2x_Z^2)}{2}, \qquad x_V = \frac{m_V^2}{M^2}.$$

$\Gamma^0_V$ are the kinematic factors of VLQ_UFO. $v$, $g = e/s_W$ and $c_W$ are computed from the SM inputs of VLQ_UFO ($G_F = 1.166 \times 10^{-5}$ GeV⁻², $m_Z = 91.1876$ GeV, $1/\alpha_{EW} = 127.9$, $m_H = 125$ GeV), which give $v = 246.26$ GeV; $g'$ and $c'_W$ are those of VLQ_v4/VLQ_v5. For X and Y, $\xi_W = 1$.

### Both chiralities

The VLQ_UFO model is built for one chirality and normalises the generation fractions to $\sum_i \zeta_i = 1$ [1]. With `LR=True` VLQcalc shares this normalisation equally between both chiralities: $\zeta_{iL} = \zeta_{iR} = \zeta_i/2$, which gives the factor $s = 1/2$ above ($s = 1$ with `LR=False`). With this choice, a given $\kappa_Q$ means about the same total width with and without `LR`. The model paper does not define the case with both chiralities, so this is a convention of VLQcalc.

### From VLQ_v4/VLQ_v5 to VLQ_UFO

The reverse conversion (`reverse=True`) uses the couplings stored in the VLQ object. For each mass it computes the branching ratio of every channel with the VLQ_v4 formulas and inverts the formula above for every channel with a non-zero branching ratio; $\kappa_Q$ is the mean of these values. If the object has no branching ratios yet, the branching ratios of the bosons are stored in `vlq.BRs` (one value per mass) and written to the VLQ_UFO card as $\xi_V$.

For X and Y the round trip VLQ_UFO → VLQ_v4 → VLQ_UFO is exact. For T and B it differs by a few per mille at 1 TeV, because the branching ratios of VLQ_v4 include the quark masses and those of VLQ_UFO do not.

## References

The numbers refer to the list in the [README](https://github.com/acanbay/VLQcalc#references).
