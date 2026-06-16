"""
Small utilities for loading per-(block, subject) result grids and plotting
summary curves.
"""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def plot_grid_curves(
    subjects,
    blocks,
    summary_grid,
    plot_bbs,
    out_fig,
    *,
    title,
    xlabel,
    ylabel,
    horizontal_line=None,
    ylim=None,
    legend_title="Block Size",
    colors=None,
    markers=None,
    linestyles=None,
):
    """
    Plot selected block-condition curves from a (block x subject) summary grid.
    """
    if colors is None:
        colors = ['#2a3b8b', '#8e44ad', '#66c2a5']
    if markers is None:
        markers = ['o', 's', '^']
    if linestyles is None:
        linestyles = ['-', '--', ':']

    plt.rcParams.update({
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'axes.labelsize': 14
    })

    fig, ax = plt.subplots(figsize=(8, 5))

    for i, bb in enumerate(plot_bbs):
        ax.plot(
            subjects,
            summary_grid[bb, :],
            label=f"{int(blocks[bb])} blocks",
            color=colors[i % len(colors)],
            marker=markers[i % len(markers)],
            linestyle=linestyles[i % len(linestyles)],
            linewidth=2,
            markersize=6
        )

    if horizontal_line is not None:
        ax.axhline(horizontal_line, color='grey', linestyle='--', linewidth=1.5, alpha=0.75)

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if ylim is not None:
        ax.set_ylim(*ylim)
    else:
        ax.set_ylim([0, 1])


    ax.legend(
        title=legend_title,
        loc="best",
        fontsize=10,
        title_fontsize=12,
        frameon=True
    )

    plt.tight_layout()

    out_fig = Path(out_fig)
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_fig, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved plot to {out_fig}")
