# Models

VLQcalc works with three FeynRules models of vector-like quarks. The models are not part of VLQcalc; they have to be installed in MadGraph5 separately.

| Model in VLQcalc | UFO folder imported in MadGraph5 | Source |
|-|-|-|
| `VLQ_UFO` | `VLQ_UFO` (4FNS only) | [FeynRules: VLQ](https://feynrules.irmp.ucl.ac.be/wiki/VLQ) [1] |
| `VLQ_v4_UFO` | `VLQ_v4_4FNS_UFO`, `VLQ_v4_5FNS_UFO` | [FeynRules: NLOModels](https://feynrules.irmp.ucl.ac.be/wiki/NLOModels) [2] |
| `VLQ_v5_UFO` | `VLQ_v5_4FNS_NLO_UFO`, `VLQ_v5_5FNS_NLO_UFO` | [FeynRules: NLOModels](https://feynrules.irmp.ucl.ac.be/wiki/NLOModels) [2] |

The flavour scheme of the folder is chosen with the `FNS` argument of the VLQ object.

## VLQ_v4_UFO and VLQ_v5_UFO

Every coupling of a VLQ to a SM quark is a separate parameter:

| VLQ | Parameters | Example |
|-|-|-|
| T | `KT` + chirality (`L`, `R`) + boson (`h`, `w`, `z`) + generation (1, 2, 3) | `KTLw3`: T-W-b, left-handed |
| B | `KB` + chirality + boson + generation | `KBRh1`: B-H-d, right-handed |
| X | `KX` + chirality + generation | `KXL3`: X-W-t, left-handed |
| Y | `KY` + chirality + generation | `KYL2`: Y-W-s, left-handed |

The masses are `MTP`, `MBP`, `MX` and `MY`. The decay formulas and SM inputs are the same in both models. VLQ_v5 adds the full CKM matrix and couplings between two VLQs (`KWTBL`, `KZTTL`, ...); these couplings are zero by default and VLQcalc does not change them.

In VLQcalc, `vlq.Kappas` holds all couplings of the VLQ with these names, one value per mass:

```python
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05, zeta=[0, 0.2, 0.8])
vlq.Kappas['KTLw3']   # same as vlq.KappaW
vlq.Kappas['KTLw2']   # coupling to the s quark
vlq.Kappas['KTRw3']   # zeros, because LR=False
```

`vlq.KappaH`, `vlq.KappaW` and `vlq.KappaZ` hold the couplings to 3rd generation quarks. The couplings to lighter quarks are in `vlq.KappaLight`, a dictionary with keys `(boson, generation)`.

## VLQ_UFO

VLQ_UFO describes a VLQ with an overall coupling $\kappa_Q$ (`KT`, `KB`, `KX`, `KY`), the branching ratios of the bosons $\xi_V$ (`xitph`, `xitpw`, `xitpz` for T; `xibph`, `xibpw`, `xibpz` for B) and the generation fractions $\zeta_i$ (`zetaTuL`, `zetaTcL`, `zetaTtL` and the right-handed `zeta...R` for T, and so on).

## Conversion between VLQ_UFO and VLQ_v4/VLQ_v5

`convertModel` matches the vertices of both models. From VLQ_UFO to VLQ_v4/VLQ_v5:

$$\kappa_{V,i} = \kappa_Q \sqrt{s\thinspace \zeta_i\thinspace \xi_V}\thickspace  F_V(M),$$

with

$$F_H = \frac{M}{v}\sqrt{\frac{1}{\Gamma^0_H}}, \qquad F_W = \sqrt{\frac{1}{\Gamma^0_W}}\thickspace \frac{g}{g'}, \qquad F_Z = \sqrt{\frac{1}{\Gamma^0_Z}}\thickspace \frac{g/c_W}{g'/c'_W},$$

$$\Gamma^0_H = \frac{(1 - x_H)^2}{2}, \qquad \Gamma^0_W = (1 - x_W)(1 + x_W - 2x_W^2), \qquad \Gamma^0_Z = \frac{(1 - x_Z)(1 + x_Z - 2x_Z^2)}{2}, \qquad x_V = \frac{m_V^2}{M^2}.$$

$\Gamma^0_V$ are the kinematic factors of VLQ_UFO. $v$, $g = e/s_W$ and $c_W$ are computed from the SM inputs of VLQ_UFO ($G_F = 1.166 \times 10^{-5}\ \mathrm{GeV}^{-2}$, $m_Z = 91.1876$ GeV, $1/\alpha_{EW} = 127.9$, $m_H = 125$ GeV), which give $v = 246.26$ GeV; $g'$ and $c'_W$ are those of VLQ_v4/VLQ_v5. For X and Y, $\xi_W = 1$.

### Both chiralities

The VLQ_UFO model is built for one chirality and normalises the generation fractions to $\sum_i \zeta_i = 1$ [1]. With `LR=True` VLQcalc shares this normalisation equally between both chiralities: $\zeta_{iL} = \zeta_{iR} = \zeta_i/2$, which gives the factor $s = 1/2$ above ($s = 1$ with `LR=False`). With this choice, a given $\kappa_Q$ gives about the same total width with and without `LR` for B and Y (within 2% above 500 GeV) and for X above 1 TeV. For T the interference term of $T \to Ht$ ($\propto \kappa_L \kappa_R m_t M$), which VLQ_UFO does not have, makes the width with `LR=True` larger. The increase grows with $\xi_H$: for $\xi_H = 0.5$ it is 33% at 500 GeV, 17% at 1 TeV and 9% at 2 TeV, for $\xi_H = 0.25$ it is 15.5%, 8% and 4%. The model paper does not define the case with both chiralities, so this is a convention of VLQcalc.

### From VLQ_v4/VLQ_v5 to VLQ_UFO

The reverse conversion (`reverse=True`) uses the couplings stored in the VLQ object and inverts the formula above exactly. With

$$q_{V,i} = \frac{1}{s}\left(\frac{\kappa_{V,i}}{F_V(M)}\right)^2 = \kappa_Q^2\thinspace \zeta_i\thinspace \xi_V$$

and $\sum_i \zeta_i = \sum_V \xi_V = 1$, the sums over the channels give

$$\kappa_Q = \sqrt{\sum_{V,i} q_{V,i}}, \qquad \xi_V = \frac{1}{\kappa_Q^2}\sum_i q_{V,i}, \qquad \zeta_i = \frac{1}{\kappa_Q^2}\sum_V q_{V,i}.$$

$\kappa_Q$ is stored in `vlq.KappaOld`, $\xi_V$ in `vlq.BRs` (one value per mass, for T and B) and $\zeta_i$ in `vlq.zeta`, and all of them are written to the VLQ_UFO card. The round trip VLQ_UFO $\to$ VLQ_v4 $\to$ VLQ_UFO returns the same $\kappa_Q$, $\xi_V$ and $\zeta_i$. The reverse conversion uses the squares of the couplings, so their signs are not kept: VLQ_UFO has positive couplings only. The widths do not depend on the signs; they matter only where amplitudes with different VLQ vertices interfere.

$\xi_V$ are parameters of the VLQ_UFO vertices. They differ from the branching ratios, which include the quark masses: for example, the couplings of T for $\mathrm{BR}(Ht) = 0.5$ at 1 TeV give $\xi_H = 0.489$. The vertices of both models agree for the converted couplings. The param_card of VLQ_UFO has $m_t = 174.3$ GeV and $m_b = 4.2$ GeV, against 172 and 4.7 GeV in VLQ_v4, so the widths MadGraph5 computes with the two models differ by up to 2% at 400 GeV and 0.25% at 1 TeV.

VLQ_UFO describes only couplings of the form $\kappa_Q\sqrt{\zeta_i\thinspace \xi_V}\thinspace F_V$. Couplings to lighter quarks computed with `calcRatioKappas` do not have exactly this form, because the quark masses enter the widths of the 3rd generation only. In this case a warning is issued; $\kappa_Q$ and $\xi_V$ reproduce the total coupling and its share per boson, and `vlq.zeta` is kept.

### Chirality of the Higgs coupling

For leading left-handed mixing, the Higgs coupling in the Lagrangian of VLQ_UFO is $\bar{T}_R H q_L$ (Eqs. 2.1 and 3.2 of [1]), the structure that follows from the Yukawa origin of the mixing. The left-handed Higgs couplings of VLQ_v4 and VLQ_v5 (`KTLh`, `KBLh`) have the same structure, $\bar{T} P_L q$. VLQcalc therefore maps the left-handed couplings of VLQ_UFO to the left-handed couplings of VLQ_v4/VLQ_v5 for all bosons. The FeynRules implementation of VLQ_UFO (VLQ.fr, version 1.2.5) attaches its left-handed Higgs coupling to $\bar{T} P_R t$ instead. The widths do not depend on this, so the couplings, widths and cards of VLQcalc are not affected; only the helicity of the quark in $Q \to Hq$ differs between events generated with the two UFO models. Since the chirality of this vertex is fixed in the VLQ_UFO files and not by a parameter, it cannot be changed through the card. For studies that depend on the helicity of the quark in $Q \to Hq$, for example angular distributions or the top-quark polarisation, VLQ_v4_UFO or VLQ_v5_UFO should be used.

## References

The numbers refer to the list in the [README](https://github.com/acanbay/VLQcalc#references).
