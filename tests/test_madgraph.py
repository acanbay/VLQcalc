import re
import warnings
import pytest
import VLQcalc.model as model
import VLQcalc.madgraph as madgraph

warnings.simplefilter('ignore')

def card(tmp_path, vlq, modelName, process='p p > tp tp~'):
    mg5 = madgraph.MG5(vlq, modelName)
    mg5.setProcess(process)
    mg5.createMG5Input(str(tmp_path / 'card'))
    return (tmp_path / 'card.dat').read_text()

def settings(text):
    return dict(re.findall(r'^set (\w+) (.*)$', text, re.M))

@pytest.mark.parametrize('text, expected', [
    ('p p > x~ j / x', 'p p > VLX j / VLX'),
    ('p p > tp~ b j', 'p p > VLT bb j'),
    ('p p > t t~ ve ve~', 'p p > tt tt veve veve'),
    ('p p > x~ j / x, (x > w+ t, t > b l+ vl)', 'p p > VLX j / VLX, (VLX > WW tt, tt > bb ll vlvl)'),
    ('p p > y y~ QED=2 $$ w+', 'p p > VLY VLY QED=2 $$ WW'),
])
def test_particle_labels(text, expected):
    vlq = model.VLQ('X')
    assert madgraph.MG5(vlq).correctText(text) == expected

def test_defines_written_once():
    vlq = model.VLQ('T')
    mg5 = madgraph.MG5(vlq)
    mg5.setProcess('p p > tp tp~, tp > w+ b')
    mg5.addProcess('p p > tp b j, tp > w+ b')
    assert mg5.text.count('define bb') == 1 and mg5.text.count('define VLT') == 1

@pytest.mark.parametrize('FNS', [4, 5])
@pytest.mark.parametrize('modelName, importName', [('VLQ_v4_UFO', 'VLQ_v4_{}FNS_UFO'), ('VLQ_v5_UFO', 'VLQ_v5_{}FNS_NLO_UFO')])
def test_new_model_cards(tmp_path, FNS, modelName, importName):
    vlq = model.VLQ('B', FNS=FNS)
    vlq.setMass([1000, 1500])
    vlq.calcRatioKappas([0.25, 0.5, 0.25], 0.1)
    text = card(tmp_path, vlq, modelName)
    assert text.splitlines()[0] == f'import model {importName.format(FNS)}'
    values = settings(text)
    assert values['MBP'] == 'scan1:[1000, 1500]'
    assert values['KBLw3'] == f'scan1:{vlq.KappaW}'
    # every other coupling of B is set to zero explicitly (the defaults of the full model are not zero)
    for name in ['KBLh1', 'KBLw2', 'KBLz1', 'KBRh3', 'KBRw3', 'KBRz3', 'KBRw1']:
        assert values[name] == '0'
    assert not any(name.startswith('zeta') or name.startswith('xi') for name in values)

def test_new_model_card_with_LR_and_generations(tmp_path):
    vlq = model.VLQ('T', LR=True)
    vlq.setMass([1000, 1500])
    vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05, zeta=[0, 0.3, 0.7])
    values = settings(card(tmp_path, vlq, 'VLQ_v4_UFO'))
    assert values['KTRh3'] == values['KTLh3'] == f'scan1:{vlq.KappaH}'
    assert values['KTLw2'] == values['KTRw2'] == f'scan1:{vlq.KappaLight[("w", 2)]}'
    assert values['KTLz1'] == values['KTRz1'] == '0'

def test_old_model_card(tmp_path):
    vlq = model.VLQ('B', LR=True)
    vlq.setMass([1000, 1500])
    vlq.convertModel(0.3, [0.5, 0.25, 0.25], zeta=[0.1, 0.1, 0.8])
    text = card(tmp_path, vlq, 'VLQ_UFO')
    assert text.splitlines()[0] == 'import model VLQ_UFO'
    values = settings(text)
    assert float(values['zetaBbL']) == float(values['zetaBbR']) == pytest.approx(0.4)
    assert float(values['zetaBdL']) == pytest.approx(0.05)
    assert values['xibph'] == '0.5' and values['KB'] == '0.3'

def test_old_model_card_after_reverse_conversion(tmp_path):
    vlq = model.VLQ('B')
    vlq.setMass([1000, 1500])
    vlq.convertModel([[0.3, 0.4], [0.2, 0.25], [0.2, 0.25]], reverse=True)
    values = settings(card(tmp_path, vlq, 'VLQ_UFO'))
    assert values['xibpw'].startswith('scan1:[') and values['KB'].startswith('scan1:[')
    assert float(values['zetaBbL']) == 1 and float(values['zetaBbR']) == 0

def test_repeated_card_creation_does_not_duplicate(tmp_path):
    vlq = model.VLQ('X')
    vlq.setMass(1000)
    vlq.calcRatioKappas(Ratio=0.1)
    mg5 = madgraph.MG5(vlq)
    mg5.setProcess('p p > x x~')
    mg5.createMG5Input(str(tmp_path / 'a'))
    mg5.createMG5Input(str(tmp_path / 'a'))
    assert (tmp_path / 'a.dat').read_text().count('generate') == 1

@pytest.mark.parametrize('call', [lambda m: m.setProcess('p p'), lambda m: m.addProcess('p p > x x~'), lambda m: m.createMG5Input()])
def test_wrong_usage_raises_value_error(call):
    with pytest.raises(ValueError):
        call(madgraph.MG5(model.VLQ('X')))

def test_unknown_model_raises_value_error():
    with pytest.raises(ValueError):
        madgraph.MG5(model.VLQ('X'), 'VLQ_v9_UFO')
