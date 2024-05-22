# VLQcalc [![arXiv](https://img.shields.io/badge/arXiv-1234.56789-b31b1b.svg?&logo=arxiv&logoColor=red)](https://arxiv.org/abs/1234.56789) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/acanbay/VLQcalc/HEAD)



>**a Python Module for Calculating Vector-like Quark Parameters**

It is capable of:
* Converting coupling constants between models VLQ_UFO $^{[1]}$ and VLQ_v4_UFO (or VLQ_v5_5FNS_UFO) $^{[2]}$ generated with FeynRules $^{[3]}$,
* Calculating VLQ decay widths,
* Computing Narrow Width Approximation coupling constants,
* Generating MadGraph5 $^{[4]}$ input card.
---

### VLQ Class

It is defined in the submodule named *model*.
```python
import VLQcalc.model as model
```
<br>

Creating a VLQ object:

```python
vlq = model.VLQ( VLQ_type, FNS, LR )
```

| Parameter | Format | Values | Default Value |
|-|-|:-:|:-:|
|VLQ_type|string|X, T, B, Y|-|
|FNS|integer|4, 5|4|
|LR|bool|True, False|False|

* $m_{b}=0$ when `FNS=4`
* ```LR=False``` allows calculations with only left-handed couplings, while ```LR=True``` allows calculations with both left and right-handed couplings.
<br>

#### `setMass` method:
The mass is passed to the method as a single value or a list of values in integer or float data type in GeV unit.

```python
vlq.setMass( 1000 )
vlq.setMAss( [1000,1500,2000] )
```
<br>

#### `convertModel` method:

```python
vlq.convertModel( Kappas, BRs, reverse )
```

| Parameter | Format | Values | Default Value |
|-|-|:-:|:-:|
|Kappas|float/integer or list|any|-|
|BRs|float/integer or list|any|1|
|reverse|bool|True, False|False|

* For T and B, the parameter named Kappas is given in the list as $[\kappa_H, \kappa_W, \kappa_Z]$. Only $\kappa_W$ is written for the X and Y.
* For T and B, the parameter named BRs is given in the list as [BR(Q→Hq), BR(Q→Wq), BR(Q→Zq)] where Q=X,T,B,Y, and q represents the 3rd family quarks of the Standard Model. Only BR(Q→Wq) is specified for X and Y.
* When the `reverse` value is `False`, the transformation is from VLQ_UFO to VLQ_v4_UFO; when `True`, it is from VLQ_v4_UFO to VLQ_UFO.

<br>

#### `calcDecay` method:

For T and B:
```python
decayH, decayW, decayZ, Gamma = vlq.calcDecay( mass, Kappas, LR )
```
For X and Y:
```python
decayW = vlq.calcDecay( Mass, Kappa, LR )
```
| Parameter | Format | Values | Default Value |
|-|-|:-:|:-:|
|Mass|float/integer or list|any|-|
|Kappas|list|any|-|
|Kappa|float/integer|any|-|
|LR|bool|True, False|False|

<br>

#### `calcRatioKappas` method:

```python
vlq.calcRatioKappas( BRs, Ratio )
```

| Parameter | Format | Values | Default Value |
|-|-|:-:|:-:|
|BRs|float/integer or list|any|-|
|Ratio|float or integer|any|-|

* For T and B, the parameter named BRs is given in the list as [BR(Q→Hq), BR(Q→Wq), BR(Q→Zq)] where Q=X,T,B,Y, and q represents the 3rd family quarks of the Standard Model. Only BR(Q→Wq) is specified for X and Y.
* Ratio is $\Gamma_Q/m_Q$ value.

### MG5 Class

It is defined in the submodule named *madgraph*.
```python
import VLQcalc.madgraph as madgraph
```

Creating a MG5 object:

```python
mg5 = madgraph.MG5(VLQ, model)
```

| Parameter | Format | Values | Default Value |
|-|-|:-:|:-:|
|VLQ|VLQ object|-|-|
|model|string|VLQ_UFO, VLQ_v4, VLQ_v5|VLQ_v4|

### References
1. M. Buchkremer, G. Cacciapaglia, A. Deandrea, and L. Panizzi. Model-independent framework for searches of top partners. *Nuclear Physics B*, 876(2):376–417, 2013.
2. B. Fuks and H. S. Shao. QCD next-to-leading-order predictions matched to parton showers for vector-like quark models. *The European Physical Journal C*, 77(2):1–21, 2017.
3. N. D. Christensen and C. Duhr. Feynrules–Feynman rules made easy. *Computer Physics Communications*,180(9):1614–1641, 2009.
4. J. Alwall, M. Herquet, F. Maltoni, O. Mattelaer, and T. Stelzer. Madgraph 5: going beyond. *Journal of High Energy Physics*, 2011(6):1–40, 2011.
