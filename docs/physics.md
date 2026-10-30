# Physics conventions

This page lists the formulas and conventions that VLQcalc uses. The widths and couplings follow the VLQ_v4_UFO and VLQ_v5_UFO models [2], so the numbers are the same as the ones MadGraph5 computes with these models.

## Decay channels

A vector-like quark $Q$ decays into a boson $V = H, W, Z$ and a quark $q_i$ of generation $i$:

| VLQ | Charge | $H$ | $W$ | $Z$ |
|-|:-:|:-:|:-:|:-:|
| T | 2/3 | $u, c, t$ | $d, s, b$ | $u, c, t$ |
| B | $-1/3$ | $d, s, b$ | $u, c, t$ | $d, s, b$ |
| X | 5/3 | - | $u, c, t$ | - |
| Y | $-4/3$ | - | $d, s, b$ | - |

The 1st and 2nd generation quarks are massless. The top quark has $m_t$ in both flavour schemes; the bottom quark has $m_b$ in the 4-flavour scheme (`FNS=4`) and is massless in the 5-flavour scheme (`FNS=5`).

Decays into another VLQ (for example $T \to W B$ in VLQ_v5) are not included.

## Couplings in the Lagrangian

In VLQ_v4_UFO and VLQ_v5_UFO [2] the interactions of a T quark with the SM quarks are

$$\mathcal{L}_T = \frac{g}{\sqrt{2}}\thinspace  W^+_\mu\thinspace  \bar{T} \gamma^\mu \left(\kappa^W_{L,i} P_L + \kappa^W_{R,i} P_R\right) d_i + \frac{g}{2c_W}\thinspace  Z_\mu\thinspace  \bar{T} \gamma^\mu \left(\kappa^Z_{L,i} P_L + \kappa^Z_{R,i} P_R\right) u_i + H\thinspace  \bar{T} \left(\kappa^H_{L,i} P_L + \kappa^H_{R,i} P_R\right) u_i + \mathrm{h.c.},$$

summed over the generations $i = 1, 2, 3$, with $g = e/s_W$ and $P_{L,R} = (1 \mp \gamma^5)/2$. For B the up- and down-type quarks are exchanged ($W^-$ with $u_i$; $Z$ and $H$ with $d_i$); X couples through $W^+$ to $u_i$ and Y through $W^-$ to $d_i$.

Each $\kappa$ is one parameter of the model, for example $\kappa^W_{L,3}$ of T is `KTLw3` (see [Models](models.md#vlq_v4_ufo-and-vlq_v5_ufo)). The Higgs couplings $\kappa^H$ multiply the Yukawa-type vertex directly, without a factor $M/v$. These are the couplings that VLQcalc computes.

## Partial widths

Each channel has a left-handed coupling $\kappa_L$ and a right-handed coupling $\kappa_R$. With `LR=False` only $\kappa_L$ is used ($\kappa_R = 0$); with `LR=True` both are used and $\kappa_R = \kappa_L$.

With $M$ the VLQ mass, $m_q$ the quark mass and

$$\lambda^{1/2}(M, m_q, m_V) = \sqrt{M^4 + m_q^4 + m_V^4 - 2M^2m_q^2 - 2M^2m_V^2 - 2m_q^2m_V^2},$$

the partial widths are

$$\Gamma(Q \to H q) = \frac{3\left[(\kappa_L^2 + \kappa_R^2)(m_q^2 + M^2 - m_H^2) + 4\kappa_L\kappa_R\thinspace  m_q M\right]}{96\pi M^3}\thinspace \lambda^{1/2}(M, m_q, m_H),$$

$$\Gamma(Q \to V q) = \frac{3g_V^2\left[(\kappa_L^2 + \kappa_R^2)\left(m_q^2 + M^2 + \frac{(M^2 - m_q^2)^2}{m_V^2} - 2m_V^2\right) - 12\kappa_L\kappa_R\thinspace  m_q M\right]}{96\pi M^3}\thinspace \lambda^{1/2}(M, m_q, m_V)$$

for $V = W, Z$, with the gauge factors

$$g_W^2 = \frac{e^2}{2 s_W^2}, \qquad g_Z^2 = \frac{e^2}{4 c_W^2 s_W^2}.$$

These are the expressions in the `decays.py` files of the VLQ_v4 models, written for a general quark mass. A channel below its threshold has zero width.

The total width is the sum over all channels, $\Gamma = \sum_{V,i} \Gamma(Q \to V q_i)$.

## Couplings from the width-to-mass ratio and branching ratios

`calcRatioKappas(BRs, Ratio, zeta)` takes

- the branching ratios of the bosons, $\mathrm{BR}_V$ (for X and Y only $W$, so $\mathrm{BR}_W = 1$),
- the width-to-mass ratio $\Gamma/M$,
- the generation fractions $\zeta = (\zeta_1, \zeta_2, \zeta_3)$, with $\sum_i \zeta_i = 1$ (default `[0, 0, 1]`).

The branching ratio of each channel is

$$\mathrm{BR}(Q \to V q_i) = \zeta_i\thinspace  \mathrm{BR}_V,$$

which is the same factorisation as in VLQ_UFO [1] for massless quarks. In VLQ_UFO, $\zeta_i$ is a fraction of the coupling squared and neglects the quark masses, so with a massive top or bottom quark the two definitions differ slightly; see [Models](models.md#conversion-between-vlq_ufo-and-vlq_v4vlq_v5).

Every partial width is proportional to $\kappa^2$; this also holds with `LR=True`, where $\kappa_L = \kappa_R = \kappa$ and the interference term is proportional to $\kappa^2$ as well. The coupling of each channel is therefore

$$\kappa_{V,i} = \sqrt{\frac{\zeta_i\thinspace  \mathrm{BR}_V\thinspace  (\Gamma/M)\thinspace  M}{\Gamma(Q \to V q_i)\big|_{\kappa = 1}}}.$$

The couplings are rounded to 6 significant digits. A channel with a non-zero branching ratio that is kinematically closed (for example $T \to Ht$ for $M \lt  m_H + m_t$) raises a `ValueError`.

## SM inputs

The default inputs are the default parameters of VLQ_v4 and VLQ_v5 (the same in all their restriction cards):

| Input | Value | Name in `setSM` |
|-|-|-|
| $m_Z$ | 91.1876 GeV | `MZ` |
| $1/\alpha_{EW}$ | 127.9 | `aEWM1` |
| $G_F$ | $1.16637 \times 10^{-5}\ \mathrm{GeV}^{-2}$ | `Gf` |
| $m_t$ | 172 GeV | `Mt` |
| $m_b$ | 4.7 GeV (4FNS only) | `Mb` |
| $m_H$ | 125 GeV | `MH` |

The other parameters are derived as in the models:

$$m_W = \sqrt{\frac{m_Z^2}{2} + \sqrt{\frac{m_Z^4}{4} - \frac{\pi\alpha_{EW} m_Z^2}{\sqrt{2}\thinspace G_F}}} = 79.824\ \mathrm{GeV}, \qquad e = \sqrt{4\pi\alpha_{EW}}, \qquad s_W^2 = 1 - \frac{m_W^2}{m_Z^2}.$$

If the param_card of MadGraph5 is changed, the same values should be given to VLQcalc, otherwise the widths differ:

```python
import VLQcalc.constants as constants

constants.setSM(Mt=173.0, MH=125.1)  # only the given inputs change
constants.getSM()                    # current inputs and the derived MW
constants.setSM()                    # back to the defaults
```

$m_W$ cannot be set directly, because it is not an independent parameter of the models.

The couplings are computed with the SM inputs set at that time, and the object keeps them: tables, files and widths of the object use them, also after a later `setSM`. MadGraph5 cards do not change SM inputs.

## Choosing the settings

**Flavour scheme (`FNS`).** In the 4-flavour scheme the b quark is massive and is not a parton of the proton; in the 5-flavour scheme it is massless and has a parton density. The scheme should be the same as the one of the model folder and the PDF set used for the event generation. In VLQcalc it changes only $m_b$ in the widths: for $\Gamma/M = 0.05$ with `LR=False` the couplings of the two schemes differ by less than 0.04% at $M = 300$ GeV and less than 0.004% at $M = 1$ TeV. With `LR=True` the interference term is linear in $m_b$, and the couplings of B differ by up to 1.8% at 300 GeV and 0.5% at 1 TeV (T and Y: 0.4% and 0.01%).

**Chirality (`LR`).** In the minimal models with one VLQ multiplet, the mixing with the SM quarks is mainly left-handed for the singlets (T, B) and mainly right-handed for the doublets ((X T), (T B), (B Y)); the other chirality is suppressed by $m_q/M$ [6]. With `LR=False` VLQcalc computes one coupling per channel and writes it to the left-handed parameters (`KTLw3`, ...). Because the widths depend on $\kappa_L^2 + \kappa_R^2$ and on the product $\kappa_L \kappa_R$, a purely right-handed coupling gives the same widths as a purely left-handed coupling of the same value; the angular distributions of the decay products differ. For a right-handed scenario, the computed couplings can therefore be used for the right-handed parameters: settings given with `mg5.addInput` are written after the couplings, so `mg5.addInput('KTLw3 0')` and `mg5.addInput('KTRw3 0.195')` (or a `scan1:[...]` list, one value per mass) replace the written values. `LR=True` sets $\kappa_R = \kappa_L$.

**Width (`Ratio`).** The couplings are computed from the total width at the given mass. For small $\Gamma/M$ the production and the decay of the VLQ factorise (narrow width approximation). For larger widths the off-shell contributions and the interference with SM processes change the cross sections and the distributions, and the full process should be simulated [7].

## Validation

- Every partial width (all VLQs, both flavour schemes, `LR` on and off, all generations) equals the expression in the `decays.py` file of VLQ_v4_4FNS_UFO or VLQ_v4_5FNS_UFO to machine precision. The test suite checks this with reference values.
- For couplings computed by VLQcalc, the widths computed by MadGraph5 (`compute_widths`, version 3.7.2) with VLQ_v4 and VLQ_v5 give the requested $\Gamma/M$ to 6 digits, and every branching ratio within $5 \times 10^{-8}$.

## References

The numbers refer to the list in the [README](https://github.com/acanbay/VLQcalc#references).
