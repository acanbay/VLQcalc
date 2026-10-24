# Physics conventions

This page lists the formulas and conventions that VLQcalc uses. The widths and couplings follow the VLQ_v4_UFO and VLQ_v5_UFO models [2], so the numbers are the same as the ones MadGraph5 computes with these models.

## Decay channels

A vector-like quark $Q$ decays into a boson $V = H, W, Z$ and a quark $q_i$ of generation $i$:

| VLQ | Charge | $H$ | $W$ | $Z$ |
|-|:-:|:-:|:-:|:-:|
| T | 2/3 | $u, c, t$ | $d, s, b$ | $u, c, t$ |
| B | −1/3 | $d, s, b$ | $u, c, t$ | $d, s, b$ |
| X | 5/3 | – | $u, c, t$ | – |
| Y | −4/3 | – | $d, s, b$ | – |

The 1st and 2nd generation quarks are massless. The top quark has $m_t$ in both flavour schemes; the bottom quark has $m_b$ in the 4-flavour scheme (`FNS=4`) and is massless in the 5-flavour scheme (`FNS=5`).

Decays into another VLQ (for example $T \to W B$ in VLQ_v5) are not included.

## Partial widths

Each channel has a left-handed coupling $\kappa_L$ and a right-handed coupling $\kappa_R$. With `LR=False` only $\kappa_L$ is used ($\kappa_R = 0$); with `LR=True` both are used and $\kappa_R = \kappa_L$.

With $M$ the VLQ mass, $m_q$ the quark mass and

$$\lambda^{1/2}(M, m_q, m_V) = \sqrt{M^4 + m_q^4 + m_V^4 - 2M^2m_q^2 - 2M^2m_V^2 - 2m_q^2m_V^2},$$

the partial widths are

$$\Gamma(Q \to H q) = \frac{3\left[(\kappa_L^2 + \kappa_R^2)(m_q^2 + M^2 - m_H^2) + 4\kappa_L\kappa_R\, m_q M\right]}{96\pi M^3}\,\lambda^{1/2}(M, m_q, m_H),$$

$$\Gamma(Q \to V q) = \frac{3g_V^2\left[(\kappa_L^2 + \kappa_R^2)\left(m_q^2 + M^2 + \frac{(M^2 - m_q^2)^2}{m_V^2} - 2m_V^2\right) - 12\kappa_L\kappa_R\, m_q M\right]}{96\pi M^3}\,\lambda^{1/2}(M, m_q, m_V)$$

for $V = W, Z$, with the gauge factors

$$g_W^2 = \frac{e^2}{2 s_W^2}, \qquad g_Z^2 = \frac{e^2}{4 c_W^2 s_W^2}.$$

These are the expressions in the `decays.py` files of the VLQ_v4 models, written for a general quark mass. A channel below its threshold has zero width.

The total width is the sum over all channels, $\Gamma = \sum_{V,i} \Gamma(Q \to V q_i)$.

## Couplings from Γ/M and branching ratios

`calcRatioKappas(BRs, Ratio, zeta)` takes

- the branching ratios of the bosons, $\mathrm{BR}_V$ (for X and Y only $W$, so $\mathrm{BR}_W = 1$),
- the width-to-mass ratio $\Gamma/M$,
- the generation fractions $\zeta = (\zeta_1, \zeta_2, \zeta_3)$, with $\sum_i \zeta_i = 1$ (default $(0, 0, 1)$).

The branching ratio of each channel is

$$\mathrm{BR}(Q \to V q_i) = \zeta_i\, \mathrm{BR}_V,$$

which is the same factorisation as in VLQ_UFO [1].

Every partial width is proportional to $\kappa^2$; this also holds with `LR=True`, where $\kappa_L = \kappa_R = \kappa$ and the interference term is proportional to $\kappa^2$ as well. The coupling of each channel is therefore

$$\kappa_{V,i} = \sqrt{\frac{\zeta_i\, \mathrm{BR}_V\, (\Gamma/M)\, M}{\Gamma(Q \to V q_i)\big|_{\kappa = 1}}}.$$

The couplings are rounded to 6 significant digits. A channel with a non-zero branching ratio that is kinematically closed (for example $T \to Ht$ for $M < m_H + m_t$) raises a `ValueError`.

Up to version 0.1.2 the couplings were found with a modified genetic algorithm [6]. With the same SM inputs, the couplings of version 0.1.2 differ from the exact solution by less than $10^{-5}$ (relative), so the change of method does not change the results; the SM inputs do (see the [changelog](../CHANGELOG.md)).

## SM inputs

The default inputs are those of the restriction cards of VLQ_v4 and VLQ_v5:

| Input | Value | Name in `setSM` |
|-|-|-|
| $m_Z$ | 91.1876 GeV | `MZ` |
| $1/\alpha_{EW}$ | 127.9 | `aEWM1` |
| $G_F$ | 1.16637 × 10⁻⁵ GeV⁻² | `Gf` |
| $m_t$ | 172 GeV | `Mt` |
| $m_b$ | 4.7 GeV (4FNS only) | `Mb` |
| $m_H$ | 125 GeV | `MH` |

The other parameters are derived as in the models:

$$m_W = \sqrt{\frac{m_Z^2}{2} + \sqrt{\frac{m_Z^4}{4} - \frac{\pi\alpha_{EW} m_Z^2}{\sqrt{2}\,G_F}}} = 79.824\ \mathrm{GeV}, \qquad e = \sqrt{4\pi\alpha_{EW}}, \qquad s_W^2 = 1 - \frac{m_W^2}{m_Z^2}.$$

If the param_card of MadGraph5 is changed, the same values should be given to VLQcalc, otherwise the widths differ:

```python
import VLQcalc.constants as constants

constants.setSM(Mt=173.0, MH=125.1)  # only the given inputs change
constants.getSM()                    # current inputs and the derived MW
constants.setSM()                    # back to the defaults
```

$m_W$ cannot be set directly, because it is not an independent parameter of the models.

## Validation

- Every partial width (all VLQs, both flavour schemes, `LR` on and off, all generations) equals the expression in the `decays.py` file of VLQ_v4_4FNS_UFO or VLQ_v4_5FNS_UFO to machine precision. The test suite checks this with reference values.
- For couplings computed by VLQcalc, the widths computed by MadGraph5 (`compute_widths`, version 3.7.2) with VLQ_v4 and VLQ_v5 give the requested $\Gamma/M$ to 6 digits, and every branching ratio within $5 \times 10^{-8}$.

## References

The numbers refer to the list in the [README](https://github.com/acanbay/VLQcalc#references).
