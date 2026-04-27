"""Shared figure-style constants for all experiment scripts.

Every experiment in ``scripts/`` calls :func:`apply` once at module top to
pin a uniform DPI, palette, font size, and seaborn theme. This guarantees
that the whole figure set in ``figures/`` shares a visual vocabulary —
required for publication-grade output.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import seaborn as sns

DPI = 300
SEED = 2026  # default seed for any stochastic experiment

# Project palette. Five named colours so every figure speaks the same language.
PALETTE = {
    "primary": "#1D76DB",   # blues — analytical / mean
    "accent": "#B60205",    # reds  — targets / lambda_2 / "danger"
    "support": "#0E8A16",   # green — successes, lambda_1 = 1
    "tertiary": "#5319E7",  # purple — secondary categories (Mid)
    "neutral": "#374151",   # near-black for axis annotations
}

LEVEL_COLOURS = {
    "Junior": "#1D76DB",
    "Mid": "#5319E7",
    "Senior": "#0E8A16",
    "Exit": "#B60205",
}


def apply() -> None:
    """Apply the shared style. Call once at the top of any experiment script."""
    sns.set_style("white")
    plt.rcParams.update(
        {
            "savefig.dpi": DPI,
            "figure.dpi": 100,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "font.family": "DejaVu Sans",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
