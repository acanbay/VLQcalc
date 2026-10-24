# Limitations and common questions

## Limitations

- **Decays.** Only the decays into a SM boson ($H$, $W$, $Z$) and a SM quark are included. Decays into another VLQ (for example $T \to W B$ in VLQ_v5) and couplings between two VLQs are not computed; VLQcalc leaves these parameters at their default values (zero in VLQ_v5).
- **Chirality.** A coupling is either left-handed only (`LR=False`) or has $\kappa_R = \kappa_L$ (`LR=True`). Purely right-handed couplings are not written to the card automatically; see [Choosing the settings](physics.md#choosing-the-settings).
- **Models.** VLQ_UFO, VLQ_v4_UFO and VLQ_v5_UFO are supported. Other models, such as models with additional scalars, are not.
- **Width.** The couplings are computed at a fixed mass from $\Gamma/M$. Whether the narrow width approximation holds for the chosen width has to be checked separately [8].
- **Mixing.** The couplings of the different channels are independent inputs. VLQcalc does not relate them to the mixing angles of a specific multiplet.

## Common questions

**The width computed by MadGraph5 is not the requested Γ/M.**
The card does not set the VLQ width, so MadGraph5 uses the value of the param_card unless `mg5.addInput('WTP auto')` (`WBP`, `WX`, `WY`) is added. If the SM inputs of the param_card were changed, the same values have to be given with `constants.setSM(...)` before the couplings are computed.

**`ValueError: ... is kinematically closed`.**
A channel with a non-zero branching ratio is below its threshold, for example $T \to H t$ for $M < m_H + m_t = 297$ GeV. Set the branching ratio of this channel to zero or use a larger mass.

**MadGraph5 generates the processes slowly.**
The full model is imported, without a restriction card, so MadGraph5 also generates the diagrams with couplings to 1st and 2nd generation quarks, even when they are set to zero. See [MadGraph5 cards](madgraph.md#vlq_v4_ufo-and-vlq_v5_ufo).

**How can a result be reproduced later?**
`vlq.save(...)` and `vlq.table()` write the settings (VLQ, flavour scheme, chirality, Γ/M, branching ratios, ζ and SM inputs) together with the couplings.

**Can VLQcalc be used without pip?**
Yes. It has no dependencies; copy the `VLQcalc` folder next to the script or add the main folder of the repository to `PYTHONPATH`.

## References

The numbers refer to the list in the [README](https://github.com/acanbay/VLQcalc#references).
