'''
Author : Ali Can Canbay - acanbay
Created: 25 October 2026

Tables, files and plots of the results stored in a VLQ object.
Plots need matplotlib, which is an optional dependency (pip install VLQcalc[plot]).
'''
import json, os
from VLQcalc import constants
from VLQcalc.model import perMass
############
CHARGE = {'VLT': '2/3', 'VLB': '-1/3', 'VLX': '5/3', 'VLY': '-4/3'}
NAME = {'VLT': 'T', 'VLB': 'B', 'VLX': 'X', 'VLY': 'Y'}

def fmt(value):
    return f'{value:.6g}'

def channels(VLQ):
    # Channels with a coupling: 3rd generation always, lighter generations only if their zeta is not zero
    from VLQcalc import decays
    out = []
    for boson in VLQ.bosons:
        for gen in [1,2,3]:
            if gen == 3 or VLQ.zeta[gen-1] != 0:
                out.append((boson, gen, decays.CHANNELS[VLQ.type][boson][gen-1]))
    return out

def kappaValues(VLQ, boson, gen):
    values = {'h': VLQ.KappaH, 'w': VLQ.KappaW, 'z': VLQ.KappaZ}[boson] if gen == 3 else VLQ.KappaLight.get((boson,gen))
    if values is None:
        return None
    return perMass(values, len(VLQ.Mass), 'VLQ')

def columns(VLQ):
    # [(plain name, LaTeX name, values)] for every column of the table
    if VLQ.Mass is None:
        raise ValueError('VLQ: Mass is not defined; use setMass first')
    cols = [('M [GeV]', r'$M$ [GeV]', list(VLQ.Mass))]
    hasKappas = False
    for boson, gen, quark in channels(VLQ):
        values = kappaValues(VLQ, boson, gen)
        if values is not None:
            hasKappas = True
            cols.append((f'k_{boson.upper()}({quark})', rf'$\kappa_{{{boson.upper()},{quark}}}$', values))
    if hasKappas:
        widths = [sum(VLQ.calcWidths(i).values()) for i in range(len(VLQ.Mass))]
        cols.append(('Gamma [GeV]', r'$\Gamma$ [GeV]', widths))
    if VLQ.KappaOld is not None:
        old = perMass(VLQ.KappaOld, len(VLQ.Mass), 'VLQ.KappaOld')
        cols.append((f'k_{NAME[VLQ.type]} (VLQ_UFO)', rf'$\kappa_{{{NAME[VLQ.type]}}}$ (VLQ_UFO)', old))
    return cols

def metadata(VLQ):
    # Settings needed to reproduce the numbers
    meta = {'VLQ': f'{NAME[VLQ.type]} (charge {CHARGE[VLQ.type]})',
            'FNS': VLQ.FNS,
            'couplings': 'left- and right-handed (kappa_L = kappa_R)' if VLQ.LR else 'left-handed',
            'Gamma/M': VLQ.Ratio,
            'branchings': None,
            'zeta': VLQ.zeta,
            'SM inputs': {name: float(fmt(value)) for name, value in constants.getSM().items()}}
    if type(VLQ.BRs) == list and all(type(BR) != list for BR in VLQ.BRs):
        if VLQ.Ratio is None and VLQ.KappaOld is not None:
            # after convertModel the inputs are the branching ratios xi of VLQ_UFO, which neglect the quark masses
            meta['xi (VLQ_UFO)'] = dict(zip(['xi(H)', 'xi(W)', 'xi(Z)'], VLQ.BRs))
        else:
            meta['branchings'] = dict(zip(VLQ.BRtext.split(','), VLQ.BRs))
    elif VLQ.type in ['VLX','VLY']:
        meta['branchings'] = {VLQ.BRtext: 1}
    return {key: value for key, value in meta.items() if value is not None}

def metadataLines(VLQ):
    lines = []
    for key, value in metadata(VLQ).items():
        if type(value) == dict:
            value = ', '.join(f'{k} = {fmt(v)}' for k, v in value.items())
        elif type(value) == list:
            value = '(' + ', '.join(fmt(v) for v in value) + ')'
        elif type(value) == float:
            value = fmt(value)
        lines.append(f'{key}: {value}')
    return lines

def table(VLQ, style='md'):
    # Markdown (style='md') or aligned plain text (style='txt') table
    if style not in ['md', 'txt']:
        raise ValueError("VLQ.table: style must be 'md' (Markdown) or 'txt' (text)")
    cols = columns(VLQ)
    rows = [[fmt(v) for v in values] for values in zip(*[c[2] for c in cols])]
    if style == 'md':
        text = ''.join(f'- {line}\n' for line in metadataLines(VLQ)) + '\n'
        text += '| ' + ' | '.join(c[1] for c in cols) + ' |\n'
        text += '|' + '|'.join('---:' for c in cols) + '|\n'
        text += ''.join('| ' + ' | '.join(row) + ' |\n' for row in rows)
        return text
    widths = [max(len(c[0]), *(len(row[i]) for row in rows)) for i, c in enumerate(cols)]
    text = ''.join(f'# {line}\n' for line in metadataLines(VLQ)) + '\n'
    text += '  '.join(c[0].rjust(w) for c, w in zip(cols, widths)) + '\n'
    text += ''.join('  '.join(v.rjust(w) for v, w in zip(row, widths)) + '\n' for row in rows)
    return text

def save(VLQ, filename):
    # The format is chosen from the extension: .md, .txt, .csv or .json
    filename = os.fspath(filename)
    extension = filename.rsplit('.', 1)[-1].lower() if '.' in filename else ''
    if extension in ['md', 'txt']:
        text = table(VLQ, extension)
    elif extension == 'csv':
        cols = columns(VLQ)
        text = ''.join(f'# {line}\n' for line in metadataLines(VLQ))
        text += ','.join(c[0] for c in cols) + '\n'
        text += ''.join(','.join(fmt(v) for v in row) + '\n' for row in zip(*[c[2] for c in cols]))
    elif extension == 'json':
        data = metadata(VLQ)
        data['mass'] = list(VLQ.Mass)
        data['kappas'] = {name: values for name, values in VLQ.Kappas.items() if any(v != 0 for v in values)}
        if VLQ.KappaOld is not None:
            data['kappa VLQ_UFO'] = [float(fmt(k)) for k in perMass(VLQ.KappaOld, len(VLQ.Mass), 'VLQ.KappaOld')]
        if VLQ.KappaW is not None:
            data['width'] = [float(fmt(sum(VLQ.calcWidths(i).values()))) for i in range(len(VLQ.Mass))]
        text = json.dumps(data, indent=2) + '\n'
    else:
        raise ValueError(f'VLQ.save: unknown file format "{filename}"; use .md, .txt, .csv or .json')
    with open(filename, 'w') as f:
        f.write(text)

# Line styles: the solid line is drawn first and wider, so that a line drawn on top of it stays visible
# when two curves coincide. Colours are used only with combined=True, where they mark the boson.
STYLE = {'h': ('-', 2.0), 'w': ('--', 1.2), 'z': ('-.', 1.2)}
GENSTYLE = {3: ('-', 2.0), 2: ('--', 1.2), 1: (':', 1.6)}
COLOR = {'h': '#0072B2', 'w': '#D55E00', 'z': '#009E73'}
WIDTHSCALE = {'h': 1.6, 'w': 1.1, 'z': 0.7}

def plot(VLQ, quantity='kappa', save=None, ax=None, logy=None, combined=False):
    # quantity: 'kappa', 'BR' or 'width'. With couplings to lighter generations there is one panel per boson,
    # or a single panel with one colour per boson if combined=True.
    # Returns the matplotlib axes (an array of axes for several panels).
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patheffects as effects
    except ImportError:
        raise ImportError('VLQ.plot needs matplotlib: pip install VLQcalc[plot]') from None
    if quantity not in ['kappa', 'BR', 'width']:
        raise ValueError("VLQ.plot: quantity must be 'kappa', 'BR' or 'width'")
    if type(combined) != bool:
        raise ValueError('VLQ.plot: combined must be True (one panel) or False (one panel per boson)')
    if VLQ.Mass is None or VLQ.KappaW is None:
        raise ValueError('VLQ.plot: no couplings to plot; use calcRatioKappas or convertModel first')
    widths = [VLQ.calcWidths(i) for i in range(len(VLQ.Mass))]
    light = any(VLQ.zeta[g] != 0 for g in [0,1])
    panels = VLQ.bosons if light and not combined else [None]
    if ax is None:
        fig, axes = plt.subplots(1, len(panels), figsize=(6 if len(panels) == 1 else 4*len(panels), 4.2),
                                 sharey=(quantity != 'kappa'), squeeze=False)
        axes = list(axes[0])
    else:
        if len(panels) > 1:
            raise ValueError('VLQ.plot: ax can only be given for a single panel; use combined=True '
                             'for couplings to lighter quarks')
        axes = [ax]
    single = len(VLQ.Mass) == 1
    for axis, panel in zip(axes, panels):
        for boson, gen, quark in channels(VLQ):
            if panel is not None and boson != panel:
                continue
            if quantity == 'kappa':
                values = kappaValues(VLQ, boson, gen)
                label = rf'$\kappa_{{{boson.upper()},{quark}}}$'
            else:
                if quantity == 'BR':
                    values = [w[(boson,gen)]/sum(w.values()) if sum(w.values()) > 0 else 0 for w in widths]
                else:
                    values = [w[(boson,gen)] for w in widths]
                label = rf'${boson.upper()}{quark}$'
            if light:
                style, width = GENSTYLE[gen]
            elif len(VLQ.bosons) > 1:
                style, width = STYLE[boson]
            else:
                style, width = '-', 1.4
            color = 'black'
            if combined and len(VLQ.bosons) > 1:
                # wider lines below narrower ones, so coinciding curves of different bosons stay visible
                color, width = COLOR[boson], width * WIDTHSCALE[boson]
            # a white edge keeps dashed and dotted lines visible on top of a solid line
            edge = [effects.withStroke(linewidth=width + 2, foreground='white')] if style != '-' else None
            axis.plot(VLQ.Mass, values, style, color=color, linewidth=width, zorder=4 - width,
                      marker='o' if single else None, markersize=4, label=label, path_effects=edge)
        if quantity == 'width':
            axis.plot(VLQ.Mass, [sum(w.values()) for w in widths], '-', color='0.7', linewidth=3,
                      zorder=1, marker='o' if single else None, markersize=4, label='total')
        positive = [v for line in axis.get_lines() for v in line.get_ydata() if v > 0]
        useLog = logy if logy is not None else (len(positive) > 0 and max(positive)/min(positive) > 100)
        if useLog:
            axis.set_yscale('log')
        axis.set_xlabel(rf'$m_{{{NAME[VLQ.type]}}}$ [GeV]', loc='right')
        if not single:
            axis.set_xlim(min(VLQ.Mass), max(VLQ.Mass))
        axis.grid(True, color='0.88', linewidth=0.6, zorder=0)
        axis.tick_params(direction='in', top=True, right=True)
        if combined and light:
            axis.legend(frameon=False, fontsize=9, loc='upper left', bbox_to_anchor=(1.02, 1))
        else:
            axis.legend(frameon=False, fontsize=9, loc='best')
    axes[0].set_ylabel({'kappa': r'$\kappa$', 'BR': 'BR', 'width': r'$\Gamma$ [GeV]'}[quantity], loc='top')
    info = f'{NAME[VLQ.type]} ({VLQ.FNS}FNS, {"LR" if VLQ.LR else "L"})'
    if VLQ.Ratio is not None:
        info += rf', $\Gamma/M = {fmt(VLQ.Ratio)}$'
    if light:
        info += ', $\\zeta = ({})$'.format(', '.join(fmt(z) for z in VLQ.zeta))
    axes[0].set_title(info, loc='left', fontsize=10)
    figure = axes[0].figure
    if len(axes) > 1:
        figure.tight_layout()
    if save:
        figure.savefig(save, dpi=300, bbox_inches='tight')
    return axes[0] if len(axes) == 1 else axes
