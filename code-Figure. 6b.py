import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import to_rgba

FIGSIZE = (8.0, 3.8)
FONT = "Arial"
FONT_SIZE = 11
LABEL_SIZE = 12
BAR_WIDTH = 0.80
BAR_EDGE_COLOR = "#242424"
BAR_EDGE_WIDTH = 0.65
AXIS_WIDTH = 1.05
Y_LIM = (-0.24, 0.265)
Y_TICKS = [-0.2, -0.1, 0, 0.1, 0.2]
Y_LABEL = "Standardized abundance"
COHORT_NAME_SIZE = 12
SHOW_COHORT_NAMES = True
SHOW_PANEL_LETTER = True

COLORS = {
    "NGT":  "#5EB3A8",
    "IFG":  "#B0B0B0",
    "IGT":  "#6BB56A",
    "CGI":  "#D97A5A",
    "T2DM": "#E8B85E",
}
BAR_ALPHA = 0.95

DATA = {'IGT cohort': {'NGT': -0.076, 'IFG': -0.056, 'IGT': 0.109, 'CGI': -0.002, 'T2DM': 0.191},
 'SCAPIS cohort': {'NGT': -0.201, 'IGT': 0.122, 'CGI': 0.23, 'T2DM': 0.096}}

def plot_b():
    plt.rcParams.update({
        "font.family": FONT, "font.size": FONT_SIZE,
        "axes.labelsize": LABEL_SIZE, "axes.linewidth": AXIS_WIDTH,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })
    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE, gridspec_kw={"width_ratios": [5.3, 4.3]})
    fig.subplots_adjust(left=0.11, right=0.98, bottom=0.18, top=0.83, wspace=0.25)
    for ax, (cohort, values) in zip(axes, DATA.items()):
        groups = list(values)
        x = np.arange(len(groups))
        ax.bar(x, list(values.values()), width=BAR_WIDTH,
               color=[to_rgba(COLORS[g], BAR_ALPHA) for g in groups], edgecolor=BAR_EDGE_COLOR,
               linewidth=BAR_EDGE_WIDTH, zorder=3)
        ax.axhline(0, color=BAR_EDGE_COLOR, linewidth=0.8, zorder=2)
        ax.set_xticks(x, groups)
        ax.set_xlim(-0.65, len(groups)-0.35)
        ax.set_ylim(*Y_LIM)
        ax.set_yticks(Y_TICKS)
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        ax.tick_params(axis="x", length=0, pad=8)
        ax.tick_params(axis="y", width=AXIS_WIDTH, length=4, pad=5)
        ax.spines[["top", "right"]].set_visible(False)
        if SHOW_COHORT_NAMES:
            ax.set_title(cohort, fontsize=COHORT_NAME_SIZE, pad=13)
    axes[0].set_ylabel(Y_LABEL, labelpad=10)
    if SHOW_PANEL_LETTER:
        fig.text(0.018, 0.97, "b", fontsize=20, fontweight="bold", va="top")
    return fig, axes

if __name__ == "__main__":
    fig, axes = plot_b()
    plt.show()
