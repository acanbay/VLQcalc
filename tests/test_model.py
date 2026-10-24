import math
import warnings
import pytest
import VLQcalc.model as model
from VLQcalc import modelConverter

warnings.simplefilter('ignore')

CASES = [('T', 4, False, [0.5, 0.25, 0.25]), ('B', 5, True, [0.25, 0.5, 0.25]),
         ('X', 4, True, None), ('Y', 5, False, None)]

def targetsHit(vlq, BRs, ratio):
    worst = 0
    for i, M in enumerate(vlq.Mass):
        widths = vlq.calcWidths(i)
        total = sum(widths.values())
        worst = max(worst, abs(total/M/ratio - 1))
        for (boson, gen), width in widths.items():
            target = (BRs or [1])[vlq.bosons.index(boson)] * vlq.zeta[gen-1]
            worst = max(worst, abs(width/total - target))
    return worst

@pytest.mark.parametrize('VLQtype, FNS, LR, BRs', CASES)
def test_ratio_kappas_hit_targets(VLQtype, FNS, LR, BRs):
    vlq = model.VLQ(VLQtype, FNS=FNS, LR=LR)
    vlq.setMass([800, 1500, 3000])
    vlq.calcRatioKappas(BRs, 0.08)
    assert targetsHit(vlq, BRs, 0.08) < 2e-5

@pytest.mark.parametrize('VLQtype, FNS, LR, BRs', CASES)
def test_ratio_kappas_with_generations(VLQtype, FNS, LR, BRs):
    vlq = model.VLQ(VLQtype, FNS=FNS, LR=LR)
    vlq.setMass([800, 1500, 3000])
    vlq.calcRatioKappas(BRs, 0.05, zeta=[0.001, 0.2, 0.799])
    assert targetsHit(vlq, BRs, 0.05) < 2e-5
    assert len(vlq.Kappas) == (18 if VLQtype in 'TB' else 6)

def test_default_zeta_gives_3rd_generation_only():
    vlq = model.VLQ('T')
    vlq.setMass(1000)
    vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05)
    assert vlq.zeta == [0, 0, 1]
    assert all(values == [0] for name, values in vlq.Kappas.items() if not name.endswith('3') or 'R' in name)
    assert vlq.Kappas['KTLh3'] == vlq.KappaH

def test_closed_channel_raises():
    vlq = model.VLQ('T')
    vlq.setMass(200)
    with pytest.raises(ValueError, match='kinematically closed'):
        vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05)

def test_set_mass_inputs():
    vlq = model.VLQ('T')
    vlq.setMass((1000, 1500))
    assert vlq.Mass == [1000, 1500]
    vlq.setMass(range(1000, 1201, 100))
    assert vlq.Mass == [1000, 1100, 1200]
    np = pytest.importorskip('numpy')
    vlq.setMass(np.arange(1000, 1600, 250))
    assert vlq.Mass == [1000, 1250, 1500] and type(vlq.Mass[0]) == int

def test_branching_sum_tolerance():
    vlq = model.VLQ('T')
    vlq.setMass(1000)
    vlq.calcRatioKappas([0.1, 0.2, 0.7], 0.05)

@pytest.mark.parametrize('call', [
    lambda: model.VLQ('Q'), lambda: model.VLQ(5), lambda: model.VLQ('T', FNS=3), lambda: model.VLQ('T', LR='yes'),
    lambda: model.VLQ('T').setMass(-5), lambda: model.VLQ('T').setMass('a'), lambda: model.VLQ('T').setMass([]),
    lambda: model.VLQ('T').calcRatioKappas([0.5, 0.25, 0.25], 0.1),
])
def test_wrong_inputs_raise_value_error(call):
    with pytest.raises(ValueError):
        call()

@pytest.mark.parametrize('BRs, zeta', [([0.5, 0.5, 0.5], None), ([0.5, 0.5], None), ([0.5, 0.25, 0.25], [1, 1, 0]),
                                       ([0.5, 0.25, 0.25], [0.5, 0.5]), ([1.2, -0.2, 0], None)])
def test_wrong_fractions_raise_value_error(BRs, zeta):
    vlq = model.VLQ('T')
    vlq.setMass(1000)
    with pytest.raises(ValueError):
        vlq.calcRatioKappas(BRs, 0.1, zeta=zeta)

def test_convert_requires_mass_and_valid_kappas():
    vlq = model.VLQ('T')
    with pytest.raises(ValueError):
        vlq.convertModel(0.5, [0.5, 0.25, 0.25])
    vlq.setMass([1000, 1500])
    with pytest.raises(ValueError):
        vlq.convertModel(5, reverse=True)
    with pytest.raises(ValueError):
        vlq.convertModel([[0.1], [0.1], [0.1]], reverse=True)

def oldModelVertex(KT, xi, zeta, boson, M, LR):
    # VLQ_UFO vertex (per chirality) divided by the VLQ_v4 gauge factor, from the VLQ_UFO parameters.py
    share = 0.5 if LR else 1.
    return KT * math.sqrt(share*zeta*xi) * modelConverter.vertexFactor(boson, M)

@pytest.mark.parametrize('LR', [False, True])
def test_conversion_matches_vertices(LR):
    vlq = model.VLQ('T', LR=LR)
    vlq.setMass([1000, 2000])
    vlq.convertModel(0.1, [0.25, 0.5, 0.25], zeta=[0.1, 0.2, 0.7])
    for i, M in enumerate(vlq.Mass):
        for boson, xi, values in [('h', 0.25, vlq.KappaH), ('w', 0.5, vlq.KappaW), ('z', 0.25, vlq.KappaZ)]:
            assert values[i] == pytest.approx(oldModelVertex(0.1, xi, 0.7, boson, M, LR), rel=1e-12)
            assert vlq.KappaLight[(boson, 1)][i] == pytest.approx(oldModelVertex(0.1, xi, 0.1, boson, M, LR), rel=1e-12)

def test_conversion_vertex_factors():
    # H: M/v * sqrt(1/Gamma0_H) with v of VLQ_UFO; W, Z: sqrt(1/Gamma0) times the ratio of the gauge couplings
    assert modelConverter.v == pytest.approx(246.2596, abs=1e-4)
    M = 1500.
    c = modelConverter.MH**2/M**2
    assert modelConverter.vertexFactor('h', M) == pytest.approx(M/modelConverter.v*math.sqrt(2/(1-c)**2), rel=1e-12)

@pytest.mark.parametrize('VLQtype, BRs', [('T', [0.25, 0.5, 0.25]), ('B', [0.5, 0.25, 0.25]), ('X', None), ('Y', None)])
def test_conversion_round_trip(VLQtype, BRs):
    vlq = model.VLQ(VLQtype)
    vlq.setMass([1000, 1500, 2000])
    vlq.convertModel(0.2, BRs) if BRs else vlq.convertModel(0.2)
    kappas = [vlq.KappaH, vlq.KappaW, vlq.KappaZ] if BRs else vlq.KappaW
    vlq.convertModel(kappas, reverse=True)
    # exact for X and Y; for T and B the branchings of VLQ_v4 include the quark masses
    assert vlq.KappaOld == pytest.approx([0.2]*3, rel=5e-3 if BRs else 1e-12)

def test_reverse_conversion_sets_branchings():
    vlq = model.VLQ('B')
    vlq.setMass([1000, 1500])
    vlq.convertModel([[0.3, 0.4], [0.2, 0.25], [0.2, 0.25]], reverse=True)
    assert len(vlq.BRs) == 3 and all(len(xi) == 2 for xi in vlq.BRs)
    assert [sum(x) for x in zip(*vlq.BRs)] == pytest.approx([1, 1], abs=1e-5)
    assert len(vlq.KappaOld) == 2

def test_forward_conversion_stores_list():
    vlq = model.VLQ('X')
    vlq.setMass(1000)
    vlq.convertModel(0.2)
    assert vlq.KappaOld == [0.2]
