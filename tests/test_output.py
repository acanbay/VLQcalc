import csv
import json
import warnings
import pytest
import VLQcalc.model as model

warnings.simplefilter('ignore')

@pytest.fixture
def vlq():
    vlq = model.VLQ('T', LR=True)
    vlq.setMass([1000, 1500])
    vlq.calcRatioKappas([0.5, 0.25, 0.25], 0.05, zeta=[0, 0.2, 0.8])
    return vlq

def test_markdown_table(vlq):
    text = vlq.table()
    assert '- Gamma/M: 0.05' in text and '- zeta: (0, 0.2, 0.8)' in text
    header = [line for line in text.splitlines() if line.startswith('| $M$')][0]
    # 3rd and 2nd generation channels, no 1st generation (zeta1 = 0), total width
    assert header.count(r'\kappa') == 6 and r'\kappa_{H,u}' not in header and r'$\Gamma$ [GeV]' in header
    assert len([line for line in text.splitlines() if line.startswith('| 1')]) == 2

def test_text_table(vlq):
    lines = vlq.table('txt').splitlines()
    assert lines[0].startswith('# VLQ: T') and lines[-1].split()[0] == '1500'

@pytest.mark.parametrize('extension', ['md', 'txt', 'csv', 'json'])
def test_save_formats(tmp_path, vlq, extension):
    path = tmp_path / f'kappas.{extension}'
    vlq.save(str(path))
    assert path.read_text().strip()

def test_csv_is_readable(tmp_path, vlq):
    path = tmp_path / 'kappas.csv'
    vlq.save(str(path))
    rows = list(csv.reader(line for line in path.read_text().splitlines() if not line.startswith('#')))
    assert rows[0][0] == 'M [GeV]' and len(rows) == 3 and all(len(row) == len(rows[0]) for row in rows)

def test_json_is_readable(tmp_path, vlq):
    path = tmp_path / 'kappas.json'
    vlq.save(str(path))
    data = json.loads(path.read_text())
    assert data['mass'] == [1000, 1500] and data['zeta'] == [0, 0.2, 0.8]
    assert data['kappas']['KTLw2'] == vlq.KappaLight[('w', 2)] and 'KTLw1' not in data['kappas']
    assert data['width'] == pytest.approx([50, 75], rel=1e-4)

def test_unknown_extension(tmp_path, vlq):
    with pytest.raises(ValueError):
        vlq.save(str(tmp_path / 'kappas.xlsx'))

@pytest.mark.parametrize('quantity', ['kappa', 'BR', 'width'])
def test_plots(tmp_path, vlq, quantity):
    pytest.importorskip('matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    axes = vlq.plot(quantity, save=str(tmp_path / f'{quantity}.png'))
    assert len(axes) == 3 and (tmp_path / f'{quantity}.png').exists()

def test_single_panel_plot(tmp_path):
    pytest.importorskip('matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    vlq = model.VLQ('X')
    vlq.setMass(range(1000, 2001, 100))
    vlq.calcRatioKappas(Ratio=0.1)
    ax = vlq.plot()
    assert len(ax.get_lines()) == 1 and ax.get_lines()[0].get_linestyle() == '-.'
    # one channel: the total equals it and is not drawn
    assert [line.get_label() for line in vlq.plot('width').get_lines()] == [r'$Wt$']

@pytest.mark.parametrize('quantity', ['kappa', 'BR', 'width'])
def test_combined_plot(tmp_path, vlq, quantity):
    pytest.importorskip('matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    ax = vlq.plot(quantity, combined=True, save=str(tmp_path / f'{quantity}.png'))
    assert (tmp_path / f'{quantity}.png').exists()
    # every channel is in the legend, alone or in a joint label; colours mark single bosons
    labels = ' '.join(line.get_label() for line in ax.get_lines())
    for name in (['H', 'W', 'Z'] if quantity != 'kappa' else ['{H,', '{W,', '{Z,']):
        assert labels.count(name) == 2
    assert {line.get_color() for line in ax.get_lines()} <= {'#0072B2', '#D55E00', '#009E73', 'black'}

def test_combined_plot_into_given_axes(vlq):
    plt = pytest.importorskip('matplotlib.pyplot')
    fig, ax = plt.subplots()
    assert vlq.plot(combined=True, ax=ax) is ax
    with pytest.raises(ValueError):
        vlq.plot(ax=ax)
    with pytest.raises(ValueError):
        vlq.plot(combined='yes')

def test_line_styles():
    pytest.importorskip('matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    vlq = model.VLQ('T')
    vlq.setMass(range(1000, 2001, 100))
    vlq.calcRatioKappas([0.6, 0.3, 0.1], 0.1)
    lines = vlq.plot('width').get_lines()
    styles = {line.get_label(): (line.get_linestyle(), line.get_color()) for line in lines}
    # a solid black line is the total width only; the channels are dashed, dash-dotted and dotted
    assert styles['total'] == ('-', 'black')
    assert [styles[name][0] for name in [r'$Ht$', r'$Wb$', r'$Zt$']] == ['--', '-.', ':']
    # all lines have the same width
    assert len({line.get_linewidth() for line in lines}) == 1


def test_equal_curves_are_drawn_once():
    pytest.importorskip('matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    vlq = model.VLQ('B')
    vlq.setMass(range(1000, 2001, 100))
    vlq.calcRatioKappas([0.25, 0.5, 0.25], 0.1)
    assert sorted(line.get_label() for line in vlq.plot('BR').get_lines()) == [r'$Hb$ = $Zb$', r'$Wt$']
    # curves that only cannot be told apart in the plot are joined with $\approx$
    labels = [line.get_label() for line in vlq.plot().get_lines()]
    assert r'$\kappa_{W,t}$ $\approx$ $\kappa_{Z,b}$' in labels


def test_joint_line_of_different_bosons_is_black():
    pytest.importorskip('matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    vlq = model.VLQ('B', LR=True)
    vlq.setMass(range(1000, 2001, 100))
    vlq.calcRatioKappas([0.25, 0.5, 0.25], 0.1, zeta=[0.1, 0.3, 0.6])
    colors = {line.get_label(): line.get_color() for line in vlq.plot('BR', combined=True).get_lines()}
    assert colors[r'$Hb$ = $Wc$ = $Zb$'] == 'black' and colors[r'$Wt$'] == '#D55E00'


def test_log_axis_compares_ratios():
    pytest.importorskip('matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    vlq = model.VLQ('B', LR=True)
    vlq.setMass(range(1000, 3001, 100))
    vlq.calcRatioKappas([0.25, 0.5, 0.25], 0.1, zeta=[0.1, 0.3, 0.6])
    axes = vlq.plot('width')
    assert axes[0].get_yscale() == 'log'
    # widths that differ by a factor 3 are separate lines on a log axis
    assert [line.get_label() for line in axes[0].get_lines()] == [r'$Hd$', r'$Hs$', r'$Hb$', 'total']


def test_combined_plot_keeps_line_styles(vlq):
    pytest.importorskip('matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    panels = [line.get_linestyle() for ax in vlq.plot('width') for line in ax.get_lines()]
    combined = [line.get_linestyle() for line in vlq.plot('width', combined=True).get_lines()]
    assert sorted(set(panels)) == sorted(set(combined))


def test_constant_values_start_at_zero():
    pytest.importorskip('matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    vlq = model.VLQ('B')
    vlq.setMass(range(1000, 2001, 100))
    vlq.calcRatioKappas([1/3, 1/3, 1/3], 0.1)
    ax = vlq.plot('BR')
    assert ax.get_ylim()[0] == 0 and ax.get_ylim()[1] > 0.3

def test_branchings_are_labelled_per_boson_with_lighter_quarks(vlq):
    text = vlq.table()
    assert '- branchings: BR(H) = 0.5, BR(W) = 0.25, BR(Z) = 0.25' in text

@pytest.mark.parametrize('zeta, label', [(None, 'BR(Wt) = 1'), ([0.2, 0.3, 0.5], 'BR(W) = 1')])
def test_X_branching_label(zeta, label):
    vlq = model.VLQ('X')
    vlq.setMass(1000)
    vlq.calcRatioKappas(Ratio=0.05, zeta=zeta)
    assert f'- branchings: {label}' in vlq.table()
