<h1 align="center">VLQcalc</h1>

<p align="center"><b>a Python Module for Calculating Vector-like Quark Couplings</b></p>

<p align="center">
  <a href="https://pypi.org/project/VLQcalc/"><img src="https://img.shields.io/pypi/v/VLQcalc.svg" alt="PyPI"></a>
  <a href="https://www.python.org"><img src="https://img.shields.io/badge/python-3.9%2B-blue.svg" alt="Python 3.9+"></a>
  <a href="https://github.com/acanbay/VLQcalc/blob/main/LICENSE"><img src="https://img.shields.io/badge/License-GPL--3.0-blue.svg" alt="License: GPL-3.0"></a>
  <a href="https://doi.org/10.5281/zenodo.12594620"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.12594620.svg" alt="DOI"></a>
  <a href="https://mybinder.org/v2/gh/acanbay/VLQcalc/HEAD?labpath=binder%2Ftutorial.ipynb"><img src="https://mybinder.org/badge_logo.svg" alt="Binder"></a>
</p>

<p align="center">
  Ali Can Canbay <a href="https://orcid.org/0000-0003-4602-473X"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="ORCID" width="15"></a>
  &emsp;&middot;&emsp;
  Orhan Cakir <a href="https://orcid.org/0000-0002-9016-138X"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="ORCID" width="15"></a>
</p>

<!---
[![arXiv](https://img.shields.io/badge/arXiv-1234.56789-b31b1b.svg?style=flat&logo=arxiv&logoColor=red)](https://arxiv.org/abs/1234.56789)
-->

<br>

VLQcalc works with the vector-like quark (VLQ) models VLQ_UFO [1], VLQ_v4_UFO and VLQ_v5_UFO [2], which were prepared with FeynRules [3]. It can:
* compute the couplings that give a chosen width-to-mass ratio (Γ/M) and chosen branching ratios,
* compute partial and total decay widths, with the same formulas and SM inputs as the models,
* include couplings to the 1st and 2nd generation quarks,
* convert couplings between VLQ_UFO and VLQ_v4_UFO/VLQ_v5_UFO,
* write MadGraph5 [4] input cards,
* write the results as tables (Markdown, CSV, JSON, text) and plots.

The [tutorial](https://github.com/acanbay/VLQcalc/blob/main/binder/tutorial.ipynb) can be run online with Binder [5] by clicking the *launch binder* button above.

> Version 0.2.0 changes some numerical results with respect to 0.1.x, because the SM inputs now follow the models. See the [changelog](https://github.com/acanbay/VLQcalc/blob/main/CHANGELOG.md).

<br>

## Installation

```console
pip install VLQcalc            # core package, no dependencies
pip install "VLQcalc[plot]"    # with matplotlib for plots
```

Python 3.9 or newer is needed. The package name is not case-sensitive for pip (`pip install vlqcalc` works as well); in Python it is imported as `VLQcalc`.

From a copy of this repository, run one of the following in its main folder:

```console
pip install .
python setup.py install        # older method; setuptools prints a deprecation warning
```

VLQcalc has no dependencies, so it can also be used without installation, for example on clusters where pip is not available: copy the `VLQcalc` folder next to your script, or add the main folder of the repository to `PYTHONPATH`.

<br>

## Quick start

```python
import VLQcalc.model as model
import VLQcalc.madgraph as madgraph

vlq = model.VLQ('T', FNS=4, LR=False)        # T quark, 4-flavour scheme, left-handed couplings
vlq.setMass(range(1000, 2001, 100))          # masses in GeV
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05) # [BR(Ht), BR(Wb), BR(Zt)] and Gamma/M

print(vlq.KappaH, vlq.KappaW, vlq.KappaZ)    # couplings to 3rd generation quarks, one value per mass
print(vlq.table())                           # Markdown table of the results
vlq.save('kappas.csv')                       # .md, .txt, .csv or .json
vlq.plot(save='kappas.pdf')                  # needs matplotlib

mg5 = madgraph.MG5(vlq, 'VLQ_v4_UFO')
mg5.setProcess('p p > tp tp~, tp > w+ b')
mg5.createMG5Input('tp_pair')                # writes tp_pair.dat
```

Couplings to lighter quarks are set with the generation fractions ζ<sub>i</sub>, for example 20% of the decays into 2nd generation quarks:

```python
vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05, zeta=[0, 0.2, 0.8])
vlq.Kappas['KTLw2']                          # every coupling, named as in the model
```

More examples are in the [examples](https://github.com/acanbay/VLQcalc/tree/main/examples) folder.

<br>

## Documentation

| Page | Content |
|-|-|
| [API reference](https://github.com/acanbay/VLQcalc/blob/main/docs/api.md) | Classes, methods, arguments and attributes |
| [Physics conventions](https://github.com/acanbay/VLQcalc/blob/main/docs/physics.md) | Couplings in the Lagrangian, decay widths, coupling calculation, SM inputs, choice of settings |
| [Models](https://github.com/acanbay/VLQcalc/blob/main/docs/models.md) | Supported UFO models, their parameters, conversion between models |
| [MadGraph5 cards](https://github.com/acanbay/VLQcalc/blob/main/docs/madgraph.md) | Writing input cards, what is written to them |
| [Limitations and common questions](https://github.com/acanbay/VLQcalc/blob/main/docs/faq.md) | What VLQcalc does not do, and answers to frequent problems |

The [documentation overview](https://github.com/acanbay/VLQcalc/tree/main/docs) shows typical workflows step by step.

<br>

## Testing

The tests are not part of the installed package. Run them from a copy of this repository:

```console
git clone https://github.com/acanbay/VLQcalc.git
cd VLQcalc
pip install ".[test]"
pytest
```

The tests compare the widths with the formulas shipped with the UFO models and check the cards, the conversions and the output files.

**AI assistance.** Version 0.2.0 was developed with the help of AI assistants (Claude, Anthropic). All output was reviewed, tested and validated by the authors.

<br>

---

### References

1. M. Buchkremer, G. Cacciapaglia, A. Deandrea, and L. Panizzi. Model-independent framework for searches of top partners. *Nuclear Physics B*, 876(2):376–417, 2013.
2. B. Fuks and H. S. Shao. QCD next-to-leading-order predictions matched to parton showers for vector-like quark models. *The European Physical Journal C*, 77(2):1–21, 2017.
3. N. D. Christensen and C. Duhr. Feynrules–Feynman rules made easy. *Computer Physics Communications*,180(9):1614–1641, 2009.
4. J. Alwall, M. Herquet, F. Maltoni, O. Mattelaer, and T. Stelzer. Madgraph 5: going beyond. *Journal of High Energy Physics*, 2011(6):1–40, 2011.
5. Jupyter et al. Binder 2.0 - Reproducible, Interactive, Sharable Environments for Science at Scale. *Proceedings of the 17th Python in Science Conference*. 2018.
6. A. C. Canbay. modifiedGA (Version 0.1.0) [Computer software]. [DOI: 10.5281/zenodo.12569505](https://doi.org/10.5281/zenodo.12569505), [GitHub: acanbay/modifiedGA](https://github.com/acanbay/modifiedGA). Used for the couplings up to version 0.1.2.
7. J. A. Aguilar-Saavedra, R. Benbrik, S. Heinemeyer, and M. Pérez-Victoria. Handbook of vectorlike quarks: Mixing and single production. *Physical Review D*, 88(9):094010, 2013.
8. S. Moretti, D. O'Brien, L. Panizzi, and H. Prager. Production of extra quarks at the Large Hadron Collider beyond the narrow width approximation. *Physical Review D*, 96(7):075035, 2017.
