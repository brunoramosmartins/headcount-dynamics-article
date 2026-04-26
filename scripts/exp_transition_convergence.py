"""Experiment: visualise convergence of P^n to the rank-1 stationary matrix.

Generates ``figures/transition_convergence_heatmap.png`` showing $P^n$ for
$n \\in \\{1, 5, 10, 50, 100\\}$ alongside the limiting matrix $\\mathbf{1}\\pi$.
Each panel uses the same colour scale so geometric convergence to the
stationary distribution is visually obvious.

Run with:

    python scripts/exp_transition_convergence.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from src.dtmc import MarkovChain

SEED = 2026
N_STEPS = (1, 5, 10, 50, 100)
LABELS = ("Junior", "Mid", "Senior", "Exit")
FIG_PATH = Path("figures/transition_convergence_heatmap.png")
DPI = 300


def headcount_chain() -> MarkovChain:
    P = np.array(
        [
            [0.93, 0.03, 0.00, 0.04],
            [0.00, 0.96, 0.02, 0.02],
            [0.00, 0.00, 0.99, 0.01],
            [0.50, 0.00, 0.00, 0.50],
        ]
    )
    return MarkovChain(P, state_labels=list(LABELS))


def main() -> None:
    np.random.default_rng(SEED)  # for any future stochastic step
    mc = headcount_chain()
    pi = mc.stationary()
    limiting = np.tile(pi, (mc.n_states, 1))

    panels = [(n, mc.n_step(n)) for n in N_STEPS] + [("∞", limiting)]

    sns.set_style("white")
    fig, axes = plt.subplots(1, len(panels), figsize=(3.0 * len(panels), 3.2))

    for ax, (label, M) in zip(axes, panels, strict=True):
        sns.heatmap(
            M,
            ax=ax,
            cmap="rocket_r",
            vmin=0.0,
            vmax=1.0,
            cbar=False,
            annot=True,
            fmt=".2f",
            annot_kws={"fontsize": 7},
            xticklabels=LABELS,
            yticklabels=LABELS,
            square=True,
            linewidths=0.4,
            linecolor="white",
        )
        title = f"$P^{{{label}}}$" if isinstance(label, int) else r"$\mathbf{1}\pi^\top$"
        ax.set_title(title, fontsize=11)
        ax.tick_params(axis="x", rotation=45)
        ax.tick_params(axis="y", rotation=0)

    fig.suptitle(
        "Convergence of $P^n$ to the rank-1 stationary matrix",
        fontsize=12,
        y=1.02,
    )
    fig.tight_layout()
    FIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_PATH, dpi=DPI, bbox_inches="tight")
    print(f"wrote {FIG_PATH}")


if __name__ == "__main__":
    main()
