'''
Author : Ali Can Canbay - acanbay
Created: 25 October 2026

SM inputs used in the decay width calculations.
Default values are the SM inputs of the VLQ_v4_UFO restrict cards [1,2]; MW, ee, sw and cw are derived
from them exactly as in the UFO model, so the widths match the ones computed by MadGraph5.
If the param_card used in MadGraph5 is changed, the same values should be set here with setSM.

====REFERENCES====
[1] Alloul, A., Christensen, N. D., Degrande, C., Duhr, C., & Fuks, B. (2014).
    FeynRules 2.0—A complete toolbox for tree-level phenomenology.
    Computer Physics Communications, 185(8), 2250-2300.
    https://feynrules.irmp.ucl.ac.be/
[2] Fuks, B., & Shao, H. S. (2017).
    QCD next-to-leading-order predictions matched to parton showers for vector-like quark models.
    The European Physical Journal C, 77(2), 135.
    https://feynrules.irmp.ucl.ac.be/wiki/NLOModels
'''
import cmath, math, numbers
############
# Independent inputs
DEFAULTS = {'MZ': 91.1876, 'aEWM1': 127.9, 'Gf': 1.16637e-5, 'Mt': 172.0, 'Mb': 4.7, 'MH': 125.0}

def _derive():
    # Derived parameters, as in parameters.py of the UFO model
    global aEW, MW, ee, sw2, sw, cw
    aEW = 1./aEWM1
    MW = cmath.sqrt(MZ**2/2. + cmath.sqrt(MZ**4/4. - (aEW*cmath.pi*MZ**2)/(Gf*cmath.sqrt(2)))).real
    ee = 2*cmath.sqrt(aEW)*cmath.sqrt(cmath.pi)
    sw2 = 1 - MW**2/MZ**2
    sw = cmath.sqrt(sw2)
    cw = cmath.sqrt(1 - sw2)

def setSM(**inputs):
    '''Changes SM inputs (MZ, aEWM1, Gf, Mt, Mb, MH), e.g. setSM(Mt=173.0); setSM() restores the defaults'''
    unknown = [name for name in inputs if name not in DEFAULTS]
    if unknown:
        raise ValueError(f'Unknown SM input(s) {unknown}; allowed inputs are {list(DEFAULTS)}')
    for name, value in inputs.items():
        if not isinstance(value, numbers.Real) or isinstance(value, bool) or not math.isfinite(value) or value <= 0:
            raise ValueError(f'SM input {name} must be a positive number')
    inputs = {name: float(value) for name, value in inputs.items()}
    values = dict(DEFAULTS) if not inputs else {name: globals()[name] for name in DEFAULTS}
    values.update(inputs)
    globals().update(values)
    _derive()

def getSM():
    '''Returns the current SM inputs together with the derived MW'''
    values = {name: globals()[name] for name in DEFAULTS}
    values['MW'] = MW
    return values

setSM()
