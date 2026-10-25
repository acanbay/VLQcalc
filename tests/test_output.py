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
    assert len(ax.get_lines()) == 1 and ax.get_lines()[0].get_linestyle() == '-'

@pytest.mark.parametrize('quantity', ['kappa', 'BR', 'width'])
def test_combined_plot(tmp_path, vlq, quantity):
    pytest.importorskip('matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    ax = vlq.plot(quantity, combined=True, save=str(tmp_path / f'{quantity}.png'))
    channels = len([line for line in ax.get_lines() if line.get_label() != 'total'])
    assert channels == 6 and (tmp_path / f'{quantity}.png').exists()
    colors = {line.get_color() for line in ax.get_lines() if line.get_label() != 'total'}
    assert len(colors) == 3

def test_combined_plot_into_given_axes(vlq):
    plt = pytest.importorskip('matplotlib.pyplot')
    fig, ax = plt.subplots()
    assert vlq.plot(combined=True, ax=ax) is ax
    with pytest.raises(ValueError):
        vlq.plot(ax=ax)
    with pytest.raises(ValueError):
        vlq.plot(combined='yes')

def test_monochrome_plot_draws_solid_line_widest():
    pytest.importorskip('matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    vlq = model.VLQ('B')
    vlq.setMass(range(1000, 2001, 100))
    vlq.calcRatioKappas([0.25, 0.5, 0.25], 0.1)
    lines = vlq.plot('BR').get_lines()
    assert {line.get_color() for line in lines} == {'black'}
    solid = [line for line in lines if line.get_linestyle() == '-']
    others = [line for line in lines if line.get_linestyle() != '-']
    assert all(s.get_linewidth() > o.get_linewidth() and s.get_zorder() < o.get_zorder()
               for s in solid for o in others)
