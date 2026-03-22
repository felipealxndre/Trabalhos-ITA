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



def pizza(data, filename):
    labels = list(data.keys())
    values = list(data.values())
    total = sum(values)
    colors = get_fabricante_colors(labels)

    fig, ax = plt.subplots(figsize=(9, 7))

    explode = [0.02] * len(labels)
    explode[0] = 0.03

    wedges, texts, autotexts = ax.pie(
        values, labels=labels,
        autopct=lambda pct: f'{pct:.1f}%\n({int(round(pct / 100 * total))})',
        colors=colors,
        explode=explode, startangle=90,
        pctdistance=0.65, labeldistance=1.10,
        wedgeprops=dict(edgecolor='white', linewidth=2),
    )

    for t in texts:
        t.set_fontsize(14)
        t.set_fontweight('bold')
        t.set_color('#333333')
    for at in autotexts:
        at.set_fontsize(12)
        at.set_fontweight('bold')
        at.set_color('white')

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
