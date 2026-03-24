import matplotlib.pyplot as plt
import matplotlib as mpl

# Configurações globais

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
mpl.rcParams['font.size'] = 12


COLORS = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D', '#44BBA4', '#E94F37', '#3B1F2B', '#393E41']


FABRICANTE_COLORS = {
    'Airbus Helicopters': '#2E86AB',
    'Airbus':             '#2E86AB',
    'Robinson':           '#A23B72',
    'Bell':               '#F18F01',
    'Bell Helicopter':    '#F18F01',
    'Leonardo':           '#C73E1D',
    'Sikorsky':           '#44BBA4',
    'Helibras':           '#E94F37',
    'MD Helicopters':     '#3B1F2B',
    'Enstrom':            '#393E41',
    'Schweizer':          '#5C4D7D',
    'Kopter':             '#7FB069',
    'Outros':             '#BBBBBB',
}

_EXTRA = ['#D4A373', '#6B4226', '#1B998B', '#FF6B6B', '#5C4D7D', '#7FB069', '#8338EC', '#FB5607']


def get_fabricante_colors(labels):
    used = {v for k, v in FABRICANTE_COLORS.items() if k in labels}
    extra_it = (c for c in _EXTRA if c not in used)
    colors = []
    for label in labels:
        if label in FABRICANTE_COLORS:
            colors.append(FABRICANTE_COLORS[label])
        else:
            colors.append(next(extra_it, '#888888'))
    return colors



def pizza(data, filename, callout=None):
    """
    callout: lista de labels que devem ser anotados externamente
             com linha em L para a esquerda (ex: ['Outros', 'MD Helicopters']).
    """
    import numpy as np

    labels = list(data.keys())
    values = list(data.values())
    total = sum(values)
    colors = get_fabricante_colors(labels)

    callout_set = set(callout or [])
    callout_indices = {i for i, lab in enumerate(labels) if lab in callout_set}

    fig, ax = plt.subplots(figsize=(9, 7))

    explode = [0.02] * len(labels)
    explode[0] = 0.03

    def _autopct(pct, idx_holder=[0]):
        i = idx_holder[0]
        idx_holder[0] += 1
        if i in callout_indices:
            return ''
        return f'{pct:.1f}%\n({int(round(pct / 100 * total))})'

    wedges, texts, autotexts = ax.pie(
        values,
        autopct=_autopct,
        colors=colors,
        explode=explode, startangle=90,
        pctdistance=0.65,
        wedgeprops=dict(edgecolor='white', linewidth=2),
    )

    for at in autotexts:
        at.set_fontsize(12)
        at.set_fontweight('bold')
        at.set_color('white')

    # ── Callouts em L para a esquerda 
    if callout_indices:
        edge_r = 1.05
        knee_r = 1.20
        text_x = -1.3

        items = []
        for idx in callout_indices:
            wedge = wedges[idx]
            # angulo do centro do setor
            ang = np.deg2rad((wedge.theta2 + wedge.theta1) / 2)
            ky = np.sin(ang) * knee_r
            items.append((idx, ang, ky))

        items.sort(key=lambda t: -t[2])

        gap = 0.1
        for i, (idx, ang, ky) in enumerate(items):
            color = '#888888' if labels[idx] == 'Outros' else colors[idx]
            knee_r_i = knee_r - i * gap

            ex = np.cos(ang) * edge_r
            ey = np.sin(ang) * edge_r
            kx = np.cos(ang) * knee_r_i
            ky_raw = np.sin(ang) * knee_r_i

            ty = np.sin(ang) * knee_r - i * gap

            ax.plot([ex, kx], [ey, ky_raw], '-', color=color, lw=1.5, clip_on=False)
            ax.plot([kx, kx], [ky_raw, ty], '-', color=color, lw=1.5, clip_on=False)
            ax.plot([kx, text_x], [ty, ty], '-', color=color, lw=1.5, clip_on=False)

            pct = values[idx] / total * 100
            txt = f'{labels[idx]}  {pct:.1f}% ({values[idx]})'
            ax.text(text_x - 0.05, ty, txt, fontsize=12, fontweight='bold',
                    color=color, ha='right', va='center')

    ax.legend(wedges, labels, loc='center left', bbox_to_anchor=(-0.3, 0.5),
              fontsize=12, labelcolor='#333',
              edgecolor='#ccc', fancybox=True, shadow=False)

    plt.tight_layout()
    plt.savefig(filename, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)



def barras(data, xlabel, ylabel, filename):
    cats = list(data.keys())
    values = list(data.values())

    fig, ax = plt.subplots(figsize=(10, 6))

    cmap = plt.cm.get_cmap('Blues')
    bar_colors = [cmap(0.35 + 0.55 * i / max(len(cats) - 1, 1)) for i in range(len(cats))]

    bars = ax.bar(
        cats, values,
        color=bar_colors,
        edgecolor='white', linewidth=1.5, width=0.6,
    )

    for bar in bars:
        val = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2, val + max(values) * 0.015,
            f'{int(val)}', ha='center', fontweight='bold', fontsize=11,
            color='#333',
        )

    ax.set_xlabel(xlabel, fontweight='bold', fontsize=13, color='#333')
    ax.set_ylabel(ylabel, fontweight='bold', fontsize=13, color='#333')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', labelsize=12, colors='#555')
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)

    plt.tight_layout()
    plt.savefig(filename, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
