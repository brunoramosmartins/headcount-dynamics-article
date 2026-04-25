"""Experiment: visualise expected time to absorption per starting level
in the headcount chain (Exit absorbing).

Generates ``figures/absorption_times.png`` with two panels:

1. Expected absorption time per starting level, with $\\pm 1\\sigma$ bars
   computed from the fundamental-matrix variance formula.
2. Probability of reaching Senior before Exit per starting level
   (computed from the modified chain in which Senior is also made absorbing).

Run with:

    python scripts/exp_absorption.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from src.absorption import (
    absorption_times,
    absorption_variances,
    hitting_probability,
)

LABELS = ("Junior", "Mid", "Senior")
FIG_PATH = Path("figures/absorption_times.png")
DPI = 300

# Headcount chain with Exit absorbing.
P = np.array(
    [
        [0.93, 0.03, 0.00, 0.04],
        [0.00, 0.96, 0.02, 0.02],
        [0.00, 0.00, 0.99, 0.01],
        [0.00, 0.00, 0.00, 1.00],
    ]
)


def main() -> None:
    sns.set_style("white")

    t = absorption_times(P)
    sigma = np.sqrt(absorption_variances(P))
    p_reach_senior = np.array(
        [hitting_probability(P, source=i, target=2) for i in range(3)]
    )

    fig, (ax_t, ax_p) = plt.subplots(1, 2, figsize=(10, 4.2))

    # Left panel — expected absorption time with +/- 1 sigma bars.
    bars = ax_t.bar(
        LABELS,
        t,
        yerr=sigma,
        capsize=8,
        color=["#1D76DB", "#5319E7", "#0E8A16"],
        edgecolor="black",
    )
    for b, v in zip(bars, t):
        ax_t.text(
            b.get_x() + b.get_width() / 2,
            v + 4,
            f"{v:.1f} mo",
            ha="center",
            fontsize=10,
        )
    ax_t.set_ylim(0, max(t + sigma) * 1.15)
    ax_t.set_ylabel("Expected months until exit")
    ax_t.set_title(
        "Expected absorption time by starting level\n"
        "(error bars: ±1 standard deviation)"
    )

    # Right panel — probability of reaching Senior before Exit.
    bars2 = ax_p.bar(
        LABELS,
        p_reach_senior,
        color=["#1D76DB", "#5319E7", "#0E8A16"],
        edgecolor="black",
    )
    for b, v in zip(bars2, p_reach_senior):
        ax_p.text(
            b.get_x() + b.get_width() / 2,
            v + 0.02,
            f"{v:.1%}",
            ha="center",
            fontsize=10,
        )
    ax_p.set_ylim(0, 1.1)
    ax_p.set_ylabel("Probability")
    ax_p.set_title("P(reach Senior before Exit) by starting level")

    fig.suptitle(
        "Headcount chain with Exit absorbing — absorption diagnostics",
        fontsize=12,
        y=1.02,
    )
    fig.tight_layout()
    FIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_PATH, dpi=DPI, bbox_inches="tight")
    print(f"wrote {FIG_PATH}")


if __name__ == "__main__":
    main()
