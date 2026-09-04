"""Publication-ready matplotlib style for Syn2bANI paper figures.

Nature Methods figure guidelines:
- Single column: 8.5 cm (3.35 in)
- Double column: 17.8 cm (7 in)
- Max height: 23.7 cm
- Vector output preferred (PDF/SVG); raster at 600 dpi for final
- Fonts: sans-serif (Arial/Helvetica), 7–9 pt labels
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from cycler import cycler

# Physical sizes in inches
SINGLE_COL = 3.35
DOUBLE_COL = 7.0
GOLDEN = 0.618

# Color palette: Okabe-Ito colorblind-safe
COLORS = {
    'orange': '#E69F00',
    'sky_blue': '#56B4E9',
    'bluish_green': '#009E73',
    'yellow': '#F0E442',
    'blue': '#0072B2',
    'vermillion': '#D55E00',
    'reddish_purple': '#CC79A7',
    'black': '#000000',
    'grey': '#999999',
    'light_grey': '#CCCCCC',
}

COLOR_CYCLE = [
    COLORS['orange'], COLORS['sky_blue'], COLORS['bluish_green'],
    COLORS['yellow'], COLORS['blue'], COLORS['vermillion'],
    COLORS['reddish_purple'], COLORS['black'],
]


def set_publication_style():
    """Apply Nature Methods-compatible defaults."""
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8,
        'axes.labelsize': 9,
        'axes.titlesize': 9,
        'axes.linewidth': 0.8,
        'axes.edgecolor': '#333333',
        'axes.labelcolor': '#333333',
        'axes.prop_cycle': cycler('color', COLOR_CYCLE),
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.minor.width': 0.5,
        'ytick.minor.width': 0.5,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.size': 3.5,
        'ytick.major.size': 3.5,
        'xtick.minor.size': 2.0,
        'ytick.minor.size': 2.0,
        'lines.linewidth': 1.2,
        'lines.markersize': 4,
        'legend.fontsize': 7,
        'legend.frameon': False,
        'legend.borderpad': 0.3,
        'legend.labelspacing': 0.3,
        'figure.dpi': 300,
        'savefig.dpi': 600,
        'savefig.format': 'png',
        'savefig.facecolor': 'white',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })
    plt.rcParams.update(mpl.rcParams)


def cm_to_inches(cm):
    return cm / 2.54


def figure_size(width_cm, aspect=GOLDEN):
    """Return (width, height) in inches for a given width in cm."""
    width = cm_to_inches(width_cm)
    return (width, width * aspect)


def label_panel(ax, label, x=-0.12, y=1.05, fontweight='bold', fontsize=10):
    """Add a bold panel label (a, b, c...) to an axes."""
    ax.text(x, y, label, transform=ax.transAxes, fontsize=fontsize,
            fontweight=fontweight, va='bottom', ha='right')


def save_figure(fig, path, formats=('png', 'pdf')):
    """Save figure in multiple formats at publication resolution."""
    for fmt in formats:
        out = f"{path}.{fmt}"
        fig.savefig(out, bbox_inches='tight', pad_inches=0.02)
        print(f"Saved {out}")


if __name__ == '__main__':
    set_publication_style()
    print("Publication style loaded.")
