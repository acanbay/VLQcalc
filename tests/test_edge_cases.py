import math
import warnings
import numpy as np
import pytest

import VLQcalc.model as model
import VLQcalc.madgraph as madgraph
import VLQcalc.constants as constants
from VLQcalc import decays


def computed(VLQtype='T', masses=(1000, 1500), **kwargs):
    vlq = model.VLQ(VLQtype, **kwargs)
    vlq.setMass(list(masses))
    if VLQtype in ['T', 'B']:
        vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.1)
    else:
        vlq.calcRatioKappas(Ratio=0.1)
    return vlq


def card(tmp_path, vlq, model_name='VLQ_v4_UFO', process='p p > tp tp~', inputs=()):
    mg5 = madgraph.MG5(vlq, model_name)
    mg5.setProcess(process)
    for line in inputs:
        mg5.addInput(line)
    path = tmp_path / 'card'
    mg5.createMG5Input(str(path))
    return (tmp_path / 'card.dat').read_text()


# Results that belong to other masses

def test_new_masses_delete_old_couplings():
    vlq = computed()
    with pytest.warns(UserWarning, match='mass changed'):
        vlq.setMass([1000, 1500, 2000])
    assert vlq.KappaH is None and vlq.KappaW is None and vlq.Ratio is None
    assert 'k_H' not in vlq.table('txt')


def test_same_masses_keep_couplings():
    vlq = computed()
    kW = list(vlq.KappaW)
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        vlq.setMass([1000, 1500])
    assert vlq.KappaW == kW


def test_coupling_set_by_hand_survives_new_masses():
    vlq = model.VLQ('B')
    vlq.setMass(1000)
    vlq.BRs = [0.25, 0.5, 0.25]
    vlq.KappaOld = [0.5]
    vlq.setMass([1000, 2000])
    assert vlq.KappaOld == [0.5] and vlq.BRs == [0.25, 0.5, 0.25]


def test_wrong_number_of_values_raises(tmp_path):
    vlq = computed()
    vlq.KappaW = [0.1, 0.2, 0.3]
    with pytest.raises(ValueError, match='one value per mass'):
        vlq.calcWidths(0)
    with pytest.raises(ValueError, match='one value per mass'):
        card(tmp_path, vlq)


def test_calcRatioKappas_deletes_converted_coupling():
    vlq = computed()
    vlq.convertModel([vlq.KappaH, vlq.KappaW, vlq.KappaZ], reverse=True)
    vlq.calcRatioKappas([0, 1, 0], 0.3)
    assert vlq.KappaOld is None


def test_forward_conversion_deletes_ratio():
    vlq = computed()
    with pytest.warns(UserWarning):
        vlq.convertModel(0.1, [0.5, 0.25, 0.25])
    assert vlq.Ratio is None
    assert 'Gamma/M' not in vlq.table()


def test_failed_calculation_keeps_previous_settings():
    vlq = model.VLQ('T')
    vlq.setMass(280)
    vlq.calcRatioKappas([0, 1, 0], 0.1)
    with pytest.raises(ValueError, match='kinematically closed'):
        vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.2)
    assert vlq.Ratio == 0.1 and vlq.BRs == [0, 1, 0]


# Reverse conversion and VLQ_UFO cards

def test_reverse_conversion_updates_xi(tmp_path):
    vlq = model.VLQ('T')
    vlq.setMass(1000)
    vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.1)
    vlq.convertModel([[0], [0.3], [0]], reverse=True)
    assert vlq.BRs == [[0], [1.0], [0]]
    text = card(tmp_path, vlq, 'VLQ_UFO')
    assert 'set xitph 0\n' in text and 'set xitpw 1\n' in text


def test_VLQ_UFO_card_without_coupling_warns(tmp_path):
    vlq = computed()
    with pytest.warns(UserWarning, match='no VLQ_UFO coupling'):
        text = card(tmp_path, vlq, 'VLQ_UFO')
    assert 'set KT ' not in text


def test_coupling_given_as_number_or_tuple(tmp_path):
    vlq = model.VLQ('T')
    vlq.setMass([1000, 1500])
    vlq.BRs = [0.5, 0.25, 0.25]
    vlq.KappaOld = 0.1
    assert 'k_T (VLQ_UFO)' in vlq.table('txt')
    vlq.KappaOld = (0.1, 0.2)
    assert 'set KT scan1:[0.1, 0.2]\n' in card(tmp_path, vlq, 'VLQ_UFO')


def test_forward_conversion_closed_channel():
    vlq = model.VLQ('T')
    vlq.setMass(100)
    with pytest.raises(ValueError, match='kinematically closed'):
        vlq.convertModel(0.1, [0.5, 0.5, 0])
    vlq.setMass(100)
    vlq.convertModel(0.1, [0, 1, 0])
    assert vlq.KappaH == [0] and vlq.KappaW[0] > 0


# MG5 cards

def test_failed_card_does_not_leave_a_half_card(tmp_path):
    vlq = model.VLQ('T')
    vlq.setMass(1000)
    mg5 = madgraph.MG5(vlq, 'VLQ_UFO')
    mg5.setProcess('p p > tp tp~')
    with pytest.raises(ValueError):
        mg5.createMG5Input(str(tmp_path / 'card'))
    vlq.convertModel(0.1, [0.5, 0.25, 0.25])
    mg5.createMG5Input(str(tmp_path / 'card'))
    text = (tmp_path / 'card.dat').read_text()
    assert text.count('generate') == 1 and text.count('import model') == 1


def test_card_without_mass_sets_light_couplings_to_zero(tmp_path):
    vlq = model.VLQ('T')
    text = card(tmp_path, vlq)
    assert 'set KTLh1 0\n' in text and 'set KTRz2 0\n' in text
    assert 'KTLh3' not in text and 'MTP' not in text


def test_run_inputs_can_change_couplings(tmp_path):
    vlq = computed()
    text = card(tmp_path, vlq, inputs=['KTRw3 0.2'])
    assert text.index('set KTRw3 0\n') < text.index('set KTRw3 0.2\n')


def test_card_numbers_and_names(tmp_path):
    vlq = model.VLQ('T', FNS=4.0)
    vlq.setMass([1000, 1500])
    vlq.convertModel(0.123456789, [0.5, 0.25, 0.25])
    text = card(tmp_path, vlq, 'VLQ_UFO')
    assert 'set KT 0.123457\n' in text and 'set zetaTtL 1\n' in text and 'set zetaTuL 0\n' in text
    assert vlq.FNS == 4 and 'VLQ_v4_4FNS_UFO' in card(tmp_path, vlq)


def test_card_file_name_with_extension(tmp_path):
    vlq = computed()
    mg5 = madgraph.MG5(vlq)
    mg5.setProcess('p p > tp tp~')
    mg5.createMG5Input(str(tmp_path / 'x.dat'))
    assert (tmp_path / 'x.dat').exists() and not (tmp_path / 'x.dat.dat').exists()


def test_prefixes_are_matched_at_the_start(tmp_path):
    vlq = computed()
    mg5 = madgraph.MG5(vlq)
    mg5.setProcess('g g > tp tp~')
    mg5.addInput('offset 5')
    assert mg5.Process == ['generate g g > VLT VLT']
    assert mg5.runInputs == ['set offset 5']


@pytest.mark.parametrize('call', [lambda m: m.setProcess(5), lambda m: m.addProcess(5), lambda m: m.addInput(5)])
def test_non_string_inputs(call):
    mg5 = madgraph.MG5(computed())
    mg5.setProcess('p p > tp tp~')
    with pytest.raises(ValueError):
        call(mg5)


def test_MG5_needs_a_VLQ_object():
    with pytest.raises(ValueError):
        madgraph.MG5(None)


# Inputs

def test_numpy_arrays_as_attributes(tmp_path):
    vlq = model.VLQ('T')
    vlq.Mass = np.array([1000., 1500.])
    vlq.calcRatioKappas(np.array([0.5, 0.25, 0.25]), np.float64(0.1), zeta=np.array([0, 0, 1]))
    vlq.KappaOld = np.array([0.1, 0.2])
    vlq.KappaH = np.array(vlq.KappaH)
    assert len(vlq.Kappas['KTLh3']) == 2 and 'k_T (VLQ_UFO)' in vlq.table('txt')
    card(tmp_path, vlq, 'VLQ_UFO')


@pytest.mark.parametrize('masses', [float('inf'), [1000, float('nan')], {1000, 2000}, np.array([[1000, 2000]]), 'abc'])
def test_bad_masses(masses):
    with pytest.raises(ValueError):
        model.VLQ('T').setMass(masses)


def test_zero_dimensional_array_mass():
    vlq = model.VLQ('T')
    vlq.setMass(np.array(1000))
    assert vlq.Mass == [1000]


@pytest.mark.parametrize('ratio', [float('nan'), float('inf'), 0, -0.1])
def test_bad_ratio(ratio):
    vlq = model.VLQ('X')
    vlq.setMass(1000)
    with pytest.raises(ValueError):
        vlq.calcRatioKappas(Ratio=ratio)


@pytest.mark.parametrize('FNS', [3, 4.5, '4', True, None])
def test_bad_FNS(FNS):
    with pytest.raises(ValueError):
        model.VLQ('T', FNS=FNS)


def test_setSM_accepts_numpy_and_rejects_non_finite():
    try:
        constants.setSM(Mt=np.float64(173.0))
        assert constants.getSM()['Mt'] == 173.0
        with pytest.raises(ValueError):
            constants.setSM(Mt=float('nan'))
    finally:
        constants.setSM()


def test_unknown_table_style():
    with pytest.raises(ValueError):
        computed().table('xyz')


def test_save_accepts_path_objects(tmp_path):
    computed().save(tmp_path / 'kappas.csv')
    assert (tmp_path / 'kappas.csv').exists()


def test_xi_are_labelled_after_forward_conversion():
    vlq = model.VLQ('B')
    vlq.setMass(1000)
    vlq.convertModel(0.5, [0.25, 0.5, 0.25])
    text = vlq.table('txt')
    assert 'xi (VLQ_UFO): xi(H) = 0.25' in text and 'BR(Hb)' not in text


# Widths below the threshold

@pytest.mark.parametrize('VLQtype, M', [('VLB', 85), ('VLT', 70), ('VLT', 100), ('VLX', 200)])
def test_closed_channels_have_zero_width(VLQtype, M):
    widths = decays.partialWidths(VLQtype, 4, M, {(b, g): 1 for b in decays.CHANNELS[VLQtype] for g in [1, 2, 3]})
    for (boson, gen), width in widths.items():
        quark = decays.CHANNELS[VLQtype][boson][gen-1]
        threshold = decays.quarkMass(quark, 4) + {'h': constants.MH, 'w': constants.MW, 'z': constants.MZ}[boson]
        if M <= threshold:
            assert width == 0
        else:
            assert width > 0


def test_closed_channel_in_calcRatioKappas():
    vlq = model.VLQ('B')
    vlq.setMass(85)
    with pytest.raises(ValueError, match='kinematically closed'):
        vlq.calcRatioKappas([0, 1, 0], 0.05)


def test_BR_plot_with_zero_width():
    pytest.importorskip('matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    vlq = model.VLQ('T')
    vlq.setMass([1000, 1500])
    vlq.convertModel([[0, 0], [0, 0], [0, 0]], reverse=True)
    vlq.plot('BR')
