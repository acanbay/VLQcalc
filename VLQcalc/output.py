'''
Author : Ali Can Canbay - acanbay
Created: 25 October 2026
Updated: 30 October 2026

Tables, files and plots of the results stored in a VLQ object.
Plots need matplotlib, which is an optional dependency (pip install VLQcalc[plot]).
'''
import json, math, os
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

def SMinputs(VLQ):
    # The SM inputs the couplings were computed with, together with the derived MW
    with constants.using(VLQ.SM):
        return {name: float(fmt(value)) for name, value in constants.getSM().items()}

def metadata(VLQ):
    # Settings needed to reproduce the numbers
    meta = {'VLQ': f'{NAME[VLQ.type]} (charge {CHARGE[VLQ.type]})',
            'FNS': VLQ.FNS,
            'couplings': 'left- and right-handed (kappa_L = kappa_R)' if VLQ.LR else 'left-handed',
            'Gamma/M': VLQ.Ratio,
            'branchings': None,
            'zeta': VLQ.zeta,
            'SM inputs': SMinputs(VLQ)}
    # With couplings to lighter quarks the branchings are the totals per boson, summed over the generations
    labels = VLQ.BRtext.split(',') if list(VLQ.zeta) == [0, 0, 1] else [f'BR({boson.upper()})' for boson in VLQ.bosons]
    if type(VLQ.BRs) == list and all(type(BR) != list for BR in VLQ.BRs):
        if VLQ.Ratio is None and VLQ.KappaOld is not None:
            # after convertModel the inputs are the branching ratios xi of VLQ_UFO, which neglect the quark masses
            meta['xi (VLQ_UFO)'] = dict(zip(['xi(H)', 'xi(W)', 'xi(Z)'], VLQ.BRs))
        else:
            meta['branchings'] = dict(zip(labels, VLQ.BRs))
    elif VLQ.type in ['VLX','VLY']:
        meta['branchings'] = {labels[0]: 1}
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

# Line styles: a solid black line is always the total width; the channels are dashed, dash-dotted and
# dotted. All lines have the same width. Colours are used only with combined=True, where they mark the
# boson; the line styles do not change. Curves that cannot be told apart are drawn once (see plot).
STYLE = {'h': '--', 'w': '-.', 'z': ':'}
GENSTYLE = {3: '--', 2: '-.', 1: ':'}
COLOR = {'h': '#0072B2', 'w': '#D55E00', 'z': '#009E73'}
LINEWIDTH = 1.5

def plot(VLQ, quantity='kappa', save=None, ax=None, logy=None, combined=False):
    # quantity: 'kappa', 'BR' or 'width'. With couplings to lighter generations there is one panel per boson,
    # or a single panel with one colour per boson if combined=True.
    # Returns the matplotlib axes (an array of axes for several panels).
    try:
        import matplotlib.pyplot as plt
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
    sharedPositive = [w[channel]/(sum(w.values()) if quantity == 'BR' else 1) for w in widths
                      for channel in [(b,g) for b, g, q in channels(VLQ)] if w[channel] > 0 and sum(w.values()) > 0]
    if quantity == 'width' and len(channels(VLQ)) > 1:
        sharedPositive += [sum(w.values()) for w in widths if sum(w.values()) > 0]
    for axis, panel in zip(axes, panels):
        curves = []
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
            style = GENSTYLE[gen] if light else STYLE[boson]
            color = COLOR[boson] if combined and len(VLQ.bosons) > 1 else 'black'
            curves.append({'values': list(values), 'labels': [label], 'style': style, 'color': color, 'exact': True})
        # The y scale is chosen before the curves are compared: logarithmic for values over more than two decades
        total = [sum(w.values()) for w in widths] if quantity == 'width' and len(channels(VLQ)) > 1 else []
        allValues = [v for c in curves for v in c['values']] + total
        positive = [v for v in allValues if v > 0]
        # panels of BR and widths share the y axis, so one decision is taken for all of them
        decide = positive if quantity == 'kappa' else sharedPositive
        useLog = logy if logy is not None else (len(decide) > 0 and max(decide)/min(decide) > 100)
        # Curves that cannot be told apart in the plot are drawn once, with a joint label: '=' if the values
        # agree to the shown precision, '$\approx$' if they differ by less than 2% of the axis (in decades on a log axis)
        if useLog and positive:
            position = lambda v: math.log10(v) if v > 0 else math.log10(min(positive)) - 1
            span = math.log10(max(positive)) - math.log10(min(positive)) + 1
        else:
            position = lambda v: v
            span = max([abs(v) for v in allValues] + [0])
        span = max(span, 1e-30)
        groups = []
        for curve in curves:
            for group in groups:
                difference = max(abs(position(a) - position(b)) for a, b in zip(curve['values'], group['values']))
                if difference <= 1e-4*span:
                    group['labels'] += curve['labels']
                    group['colors'].add(curve['color'])
                    break
                if difference <= 0.02*span:
                    group['labels'] += curve['labels']
                    group['colors'].add(curve['color'])
                    group['exact'] = False
                    break
            else:
                curve['colors'] = {curve['color']}
                groups.append(curve)
        for group in groups:
            label = (' = ' if group['exact'] else r' $\approx$ ').join(group['labels'])
            # a joint line of different bosons is black, since the colour marks one boson
            color = group['color'] if len(group['colors']) == 1 else 'black'
            axis.plot(VLQ.Mass, group['values'], group['style'], color=color, linewidth=LINEWIDTH, zorder=3,
                      marker='o' if single else None, markersize=4, label=label)
        if total:
            # the total equals the only channel if there is one channel, so it is drawn only for several
            axis.plot(VLQ.Mass, total, '-', color='black', linewidth=LINEWIDTH,
                      zorder=2, marker='o' if single else None, markersize=4, label='total')
        if useLog:
            axis.set_yscale('log')
            if positive:
                # half a decade of space above the largest value, for the legend
                axis.set_ylim(top=3*max(sharedPositive if quantity != 'kappa' else positive))
        else:
            # linear axes start at zero, with some space above the largest value; constant values are
            # drawn as flat lines, not as rounding noise
            axis.ticklabel_format(axis='y', useOffset=False)
            largest = max([v for line in axis.get_lines() for v in line.get_ydata()], default=0)
            axis.set_ylim(0, 1.08*largest if largest > 0 else 1)
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
