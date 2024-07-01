# VLQcalc [![DOI](https://zenodo.org/badge/803805406.svg)](https://zenodo.org/doi/10.5281/zenodo.12594620) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/acanbay/VLQcalc/HEAD?labpath=binder%2Ftutorial.ipynb)

<!--- 
[![arXiv](https://img.shields.io/badge/arXiv-1234.56789-b31b1b.svg?style=flat&logo=arxiv&logoColor=red)](https://arxiv.org/abs/1234.56789)
-->

>**a Python Module for Calculating Vector-like Quark Couplings**

**Contributors:**
Ali Can Canbay<a href="https://orcid.org/0000-0003-4602-473X"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="Ali Can Canbay ORCID" width="15"></a>
&
Orhan Cakir<a href="https://orcid.org/0000-0002-9016-138X"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="Orhan Cakir ORCID" width="15"></a> 

<br>

It is capable of:
* Converting coupling constants between models VLQ_UFO $^{[1]}$ and VLQ_v4_UFO $^{[2]}$ prepared within the FeynRules $^{[3]}$ framework,
* Calculating VLQ decay widths,
* Computing coupling constants for width-to-mass ratio ($\Gamma/m$),
* Generating MadGraph5 $^{[4]}$ input card.

<br>

By clicking on the '*launch binder*' button above, the module can be run online with Binder $^{[5]}$.

[Click here](https://github.com/acanbay/VLQcalc/blob/main/binder/tutorial.ipynb) for the tutorial.

<br>

**Installation:**

Download the [latest release](https://github.com/acanbay/VLQcalc/releases), extract it, enter the extracted file, and run the following command via the console.
```console
python setup.py install
```

---

<br>

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

* $m_{b}=0$ when `FNS=5`
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
It converts the couplings in VLQ_UFO and VLQ_v4_UFO models to each other.

```python
vlq.convertModel( Kappas, BRs, reverse )
```

| Parameter | Format | Values | Default Value |
|-|-|:-:|:-:|
|Kappas|float/integer or list|any|-|
|BRs|float/integer or list|any|1|
|reverse|bool|True, False|False|

* `Kappas` is
    * $\kappa_Q$ where Q is X, T, B or Y for converting couplings from VLQ_UFO to VLQ_v4_UFO.
    * given in the list as $[\kappa_H, \kappa_W, \kappa_Z]$ for T and B, while it is only $\kappa_W$ for X and Y when converting couplings from VLQ_v4_UFO to VLQ_UFO.
* For T and B, the parameter named `BRs` is given in the list as [BR(Q→Hq), BR(Q→Wq), BR(Q→Zq)] where q represents the 3rd family quarks of the Standard Model. Only BR(Q→Wq) is specified for X and Y.
* When the `reverse` value is `False`, the conversion is from VLQ_UFO to VLQ_v4_UFO; when `True`, it is from VLQ_v4_UFO to VLQ_UFO.

<br>

#### `calcDecay` method:
*It only works with the **VLQ_v4_UFO** model. Converted couplings can be used for VLQ_UFO model*
<br><br>

For T and B:
```python
decayH, decayW, decayZ, Gamma = vlq.calcDecay( Mass, Kappas, LR )
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
It is used to calculate the couplings according to the $\Gamma_Q/m_Q$ ratio.<br>
*It only works with the **VLQ_v4_UFO** model. Couplings can be converted to VLQ_UFO model after calculations are made according to this model.*

The **modified Genetic Algorithm (mGA)** $^{[5]}$ is used to calculate couplings according to the $\Gamma_Q/m_Q$ ratio.

```python
vlq.calcRatioKappas( BRs, Ratio )
```

* Ratio is $\Gamma_Q/m_Q$ value.

<br>

### MG5 Class

It is defined in the submodule named *madgraph*.
```python
import VLQcalc.madgraph as madgraph
```

Creating a MG5 object:
```python
mg5 = madgraph.MG5(VLQ, model)
```
* `model` parameter can be VLQ_UFO or VLQ_v4_UFO in string format.

<br>

In the MG5 object, there are two different methods for entering a process: `setProcess` and `addProcess` methods, which take their parameters as strings, are used respectively to define the main process and additional processes.
```python
mg5.setProcess( process )
mg5.addProcess( process )
```

<br>
MG5 object also have the properties in string format listed in the table below.

| Property | Values | Default Value |
|-|-|:-:|
|shower|OFF, pythia8, ...|OFF|
|detector|OFF, Delphes, ...|OFF|
|analysis|OFF, ExRoot, MadAnalysis, ...|OFF|
|madspin|OFF, ON, onshell, full|OFF|
|reweight|OFF, ON|OFF|

<br>

`addInput` method, which takes a string parameter, is used to define inputs that allow modifications to be made on simulation cards. After all definitions have been made, an MG5 input card can be generated using the `createMG5Input` method, which takes the output file name (in string format) as a parameter.
```python
mg5.addInput( input )
mg5.createMG5Input( file_name )
```

<br>

When creating an input card, the processes defined in `setProcess` and `addProcess` are expanded to include both particles and antiparticles.

<br>

---

### References

1. M. Buchkremer, G. Cacciapaglia, A. Deandrea, and L. Panizzi. Model-independent framework for searches of top partners. *Nuclear Physics B*, 876(2):376–417, 2013.
2. B. Fuks and H. S. Shao. QCD next-to-leading-order predictions matched to parton showers for vector-like quark models. *The European Physical Journal C*, 77(2):1–21, 2017.
3. N. D. Christensen and C. Duhr. Feynrules–Feynman rules made easy. *Computer Physics Communications*,180(9):1614–1641, 2009.
4. J. Alwall, M. Herquet, F. Maltoni, O. Mattelaer, and T. Stelzer. Madgraph 5: going beyond. *Journal of High Energy Physics*, 2011(6):1–40, 2011.
5. Jupyter et al. Binder 2.0 - Reproducible, Interactive, Sharable Environments for Science at Scale. *Proceedings of the 17th Python in Science Conference*. 2018.
6. A. C. Canbay. modifiedGA (Version 0.1.0) [Computer software]. [DOI: 10.5281/zenodo.12569505](https://doi.org/10.5281/zenodo.12569505), [GitHub: acanbay/modifiedGA](https://github.com/acanbay/modifiedGA)

