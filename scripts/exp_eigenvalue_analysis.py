"""Experiment: visualise the spectrum of the headcount transition matrix.

Generates ``figures/eigenvalue_spectrum.png`` showing the eigenvalues of $P$
in the complex plane against the unit circle, with the spectral gap and
mixing-time annotation.

Run with:

    python scripts/exp_eigenvalue_analysis.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from src.dtmc import MarkovChain

LABELS = ("Junior", "Mid", "Senior", "Exit")
FIG_PATH = Path("figures/eigenvalue_spectrum.png")
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
    mc = headcount_chain()
    eigs = mc.eigenvalues()
    gap = mc.spectral_gap()
    t_mix = mc.mixing_time()

    fig, (ax_disc, ax_bar) = plt.subplots(1, 2, figsize=(9, 4))

    # Left panel: complex plane with unit circle and eigenvalues
    theta = np.linspace(0, 2 * np.pi, 360)
    ax_disc.plot(np.cos(theta), np.sin(theta), color="black", lw=1, alpha=0.6)
    ax_disc.axhline(0, color="grey", lw=0.5, alpha=0.6)
    ax_disc.axvline(0, color="grey", lw=0.5, alpha=0.6)
    ax_disc.scatter(eigs.real, eigs.imag, s=80, color="#B60205", zorder=5)
    for lam in eigs:
        ax_disc.annotate(
            f"{lam.real:.3f}" + (f"{lam.imag:+.3f}i" if abs(lam.imag) > 1e-10 else ""),
            xy=(lam.real, lam.imag),
            xytext=(8, 6),
            textcoords="offset points",
            fontsize=8,
        )
    ax_disc.set_xlim(-1.2, 1.2)
    ax_disc.set_ylim(-1.2, 1.2)
    ax_disc.set_aspect("equal")
    ax_disc.set_title("Eigenvalues of $P$ in the complex plane")
    ax_disc.set_xlabel("Re($\\lambda$)")
    ax_disc.set_ylabel("Im($\\lambda$)")

    # Right panel: bar chart of |lambda_k|
    moduli = np.abs(eigs)
    indices = np.arange(1, len(moduli) + 1)
    bars = ax_bar.bar(indices, moduli, color="#1D76DB", edgecolor="black")
    bars[0].set_color("#0E8A16")  # highlight lambda_1 = 1
    if len(bars) > 1:
        bars[1].set_color("#D93F0B")  # highlight lambda_2 (sets the gap)
    ax_bar.axhline(1.0, color="grey", lw=0.5, ls="--")
    ax_bar.set_ylim(0, 1.15)
    ax_bar.set_xticks(indices)
    ax_bar.set_xticklabels([f"$\\lambda_{i}$" for i in indices])
    ax_bar.set_ylabel(r"$|\lambda_k|$")
    ax_bar.set_title(
        f"Spectral gap $\\gamma = 1 - |\\lambda_2| \\approx {gap:.3f}$  |  "
        f"$t_{{\\mathrm{{mix}}}}(1/4) = {t_mix}$ months"
    )
    for i, m in enumerate(moduli):
        ax_bar.text(i + 1, m + 0.02, f"{m:.3f}", ha="center", fontsize=9)

    fig.suptitle("Spectrum of the headcount transition matrix", fontsize=12, y=1.02)
    fig.tight_layout()
    FIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_PATH, dpi=DPI, bbox_inches="tight")
    print(f"wrote {FIG_PATH}")


if __name__ == "__main__":
    main()
