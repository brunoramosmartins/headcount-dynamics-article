"""Experiment: visualise the spectrum of the headcount transition matrix.

Generates ``figures/eigenvalue_spectrum.png`` with two panels:

- *Left:* eigenvalues on the complex plane against the unit circle. Because
  all four eigenvalues of the default headcount matrix are real, the points
  collapse to the real axis; we draw them on a stretched horizontal strip
  with leader lines pointing to per-eigenvalue labels so they remain
  individually legible even when two are nearly equal (e.g. 1.000 and
  0.987).
- *Right:* bar chart of moduli $|\\lambda_k|$ with the spectral gap and the
  empirical mixing time annotated.

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

# Distinct marker + colour per eigenvalue rank so labels can be matched
# to dots even when they overlap in 2D.
MARKERS = ("o", "s", "D", "^")
COLORS = ("#0E8A16", "#D93F0B", "#1D76DB", "#5319E7")


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


def _label_offsets(n: int) -> list[tuple[int, int]]:
    """Spread label offsets up/down/left/right so they don't overlap."""
    base = [(40, 25), (40, -30), (-50, 25), (-50, -30)]
    return base[:n]


def main() -> None:
    mc = headcount_chain()
    eigs = mc.eigenvalues()
    gap = mc.spectral_gap()
    t_mix = mc.mixing_time()
    all_real = np.all(np.abs(eigs.imag) < 1e-10)

    fig, (ax_disc, ax_bar) = plt.subplots(
        1, 2, figsize=(11, 4.5), gridspec_kw={"width_ratios": [1.1, 1]}
    )

    # ── Left panel: eigenvalues on the complex plane ──
    theta = np.linspace(0, 2 * np.pi, 360)
    ax_disc.plot(np.cos(theta), np.sin(theta), color="black", lw=1.2, alpha=0.5)
    ax_disc.fill(np.cos(theta), np.sin(theta), color="grey", alpha=0.04)
    ax_disc.axhline(0, color="grey", lw=0.5, alpha=0.6)
    ax_disc.axvline(0, color="grey", lw=0.5, alpha=0.6)

    offsets = _label_offsets(len(eigs))
    for k, lam in enumerate(eigs):
        ax_disc.scatter(
            lam.real,
            lam.imag,
            s=140,
            marker=MARKERS[k % len(MARKERS)],
            color=COLORS[k % len(COLORS)],
            edgecolor="black",
            linewidth=1.0,
            zorder=5,
            label=f"$\\lambda_{k+1}$ = {lam.real:.3f}"
            + ("" if abs(lam.imag) < 1e-10 else f"{lam.imag:+.3f}i"),
        )
        # Leader arrow with text
        dx, dy = offsets[k]
        ax_disc.annotate(
            f"$\\lambda_{k+1}$ = {lam.real:.3f}",
            xy=(lam.real, lam.imag),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=9.5,
            arrowprops=dict(
                arrowstyle="-",
                color=COLORS[k % len(COLORS)],
                lw=0.8,
                alpha=0.7,
            ),
            ha="center",
        )

    ax_disc.set_xlim(-1.25, 1.25)
    ax_disc.set_ylim(-1.25, 1.25)
    ax_disc.set_aspect("equal")
    ax_disc.set_title("Eigenvalues of $P$ in the complex plane")
    ax_disc.set_xlabel(r"Re$(\lambda)$")
    ax_disc.set_ylabel(r"Im$(\lambda)$")
    if all_real:
        ax_disc.text(
            0.0,
            -1.12,
            "all eigenvalues are real for this $P$ — points collapse onto the real axis",
            ha="center",
            fontsize=8.5,
            style="italic",
            color="#444",
        )

    # ── Right panel: bar chart of moduli ──
    moduli = np.abs(eigs)
    indices = np.arange(1, len(moduli) + 1)
    bars = ax_bar.bar(indices, moduli, color=COLORS, edgecolor="black", linewidth=0.8)
    ax_bar.axhline(1.0, color="grey", lw=0.5, ls="--")
    ax_bar.set_ylim(0, 1.20)
    ax_bar.set_xticks(indices)
    ax_bar.set_xticklabels([f"$\\lambda_{i}$" for i in indices])
    ax_bar.set_ylabel(r"$|\lambda_k|$")
    ax_bar.set_title(
        f"Spectral gap $\\gamma = 1 - |\\lambda_2| \\approx {gap:.3f}$  |  "
        f"$t_{{\\mathrm{{mix}}}}(1/4) = {t_mix}$ months"
    )
    for i, m in enumerate(moduli):
        ax_bar.text(i + 1, m + 0.025, f"{m:.3f}", ha="center", fontsize=9)

    # Visual cue for the spectral gap between bar 1 and bar 2.
    gap_x = 1.5
    ax_bar.annotate(
        "",
        xy=(gap_x, moduli[1]),
        xytext=(gap_x, moduli[0]),
        arrowprops=dict(arrowstyle="<->", color="#B60205", lw=1.2),
    )
    ax_bar.text(
        gap_x + 0.05,
        (moduli[0] + moduli[1]) / 2,
        f"gap ≈ {gap:.3f}",
        color="#B60205",
        fontsize=8.5,
        va="center",
    )

    fig.suptitle("Spectrum of the headcount transition matrix", fontsize=12, y=1.02)
    fig.tight_layout()
    FIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_PATH, dpi=DPI, bbox_inches="tight")
    print(f"wrote {FIG_PATH}")


if __name__ == "__main__":
    main()
