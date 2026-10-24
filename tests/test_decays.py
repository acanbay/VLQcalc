import math
import pytest
from VLQcalc import decays, constants
from reference_widths import MODEL_WIDTHS

KAPPA = {1: 0.3, 2: 0.2, 3: 0.5}

@pytest.fixture(autouse=True)
def defaultSM():
    constants.setSM()
    yield
    constants.setSM()

@pytest.mark.parametrize('VLQtype, FNS, LR, boson, gen, width', MODEL_WIDTHS)
def test_widths_match_model_formulas(VLQtype, FNS, LR, boson, gen, width):
    Kappas = {(b, g): KAPPA[g] for b in decays.CHANNELS[VLQtype] for g in (1, 2, 3)}
    widths = decays.partialWidths(VLQtype, FNS, 1500., Kappas, LR)
    assert widths[(boson, gen)] == pytest.approx(width, rel=1e-12)

@pytest.mark.parametrize('FNS', [4, 5])
@pytest.mark.parametrize('LR', [False, True])
@pytest.mark.parametrize('M', [300., 700., 1500., 4000.])
def test_general_formulas_equal_3rd_generation_functions(FNS, LR, M):
    old = {'VLT': decays.decayTP4, 'VLB': decays.decayBP4, 'VLX': decays.decayX4, 'VLY': decays.decayY4} if FNS == 4 else \
          {'VLT': decays.decayTP5, 'VLB': decays.decayBP5, 'VLX': decays.decayX5, 'VLY': decays.decayY5}
    for VLQtype, function in old.items():
        if VLQtype in ['VLT', 'VLB']:
            kappas = [0.3, 0.7, 0.2]
            reference = function(M, kappas, LR)[:3]
            widths = decays.partialWidths(VLQtype, FNS, M, {('h', 3): 0.3, ('w', 3): 0.7, ('z', 3): 0.2}, LR)
            values = [widths[('h', 3)], widths[('w', 3)], widths[('z', 3)]]
        else:
            reference = [function(M, 0.4, LR)]
            values = [decays.partialWidths(VLQtype, FNS, M, {('w', 3): 0.4}, LR)[('w', 3)]]
        assert values == pytest.approx(reference, rel=1e-12, abs=1e-15)

def test_madgraph_reference_value():
    # Total width of X computed by MadGraph5 with VLQ_v4_4FNS_UFO for M = 1000 GeV, KXL3 = 1
    assert decays.decayX4(1000, 1.0) == pytest.approx(299.8443, abs=1e-4)

def test_closed_channel_has_zero_width():
    widths = decays.partialWidths('VLT', 4, 200., {('h', 3): 1, ('w', 3): 1, ('z', 3): 1})
    assert widths[('h', 3)] == 0 and widths[('z', 3)] == 0 and widths[('w', 3)] > 0

def test_setSM_changes_and_restores_inputs():
    reference = decays.decayX4(1000, 1.0)
    constants.setSM(Mt=173.0)
    assert constants.Mt == 173.0 and constants.MH == 125.0
    assert decays.decayX4(1000, 1.0) != pytest.approx(reference)
    constants.setSM(MH=125.1)
    assert constants.Mt == 173.0 and constants.MH == 125.1
    constants.setSM()
    assert decays.decayX4(1000, 1.0) == pytest.approx(reference, rel=1e-14)
    assert constants.getSM()['MW'] == pytest.approx(79.82436, abs=1e-5)

@pytest.mark.parametrize('inputs', [{'MW': 80.}, {'Mt': -1.}, {'Mt': 'x'}, {'MH': True}])
def test_setSM_rejects_wrong_inputs(inputs):
    with pytest.raises(ValueError):
        constants.setSM(**inputs)
