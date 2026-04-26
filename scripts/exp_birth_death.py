"""Experiment: stationary distribution and Kolmogorov trajectory for the
M/M/inf headcount model.

Generates two figures:

- ``figures/birth_death_stationary.png`` — bar chart of the stationary
  distribution computed via detailed balance, overlaid with the analytical
  Poisson(rho) PMF.
- ``figures/kolmogorov_trajectory.png`` — expected team size E[n(t)] from
  the closed-form ODE solution alongside an empirical band derived from the
  exact transient distribution at each time (5th-95th percentiles).

Run with:

    python scripts/exp_birth_death.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import poisson

from src.ctmc import BirthDeathProcess

LAMBDA = 2.0  # hires per month
MU = 0.025    # per-capita attrition per month
N0 = 45       # initial team size
N_MAX = 200   # generator truncation
DPI = 300

OUT_STATIONARY = Path("figures/birth_death_stationary.png")
OUT_TRAJECTORY = Path("figures/kolmogorov_trajectory.png")


def make_chain() -> BirthDeathProcess:
    return BirthDeathProcess.from_constant_hiring(
        lambda_rate=LAMBDA, mu_rate=MU, n_max=N_MAX
    )


# ──────────────────────────────────────────────────────────────────────────
# Figure 1 — stationary distribution
# ──────────────────────────────────────────────────────────────────────────


def plot_stationary(bd: BirthDeathProcess) -> None:
    pi = bd.stationary_birth_death()
    rho = LAMBDA / MU
    n_plot = np.arange(0, 130)
    pi_plot = pi[n_plot]
    poisson_pmf = poisson.pmf(n_plot, mu=rho)

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.bar(
        n_plot,
        pi_plot,
        color="#1D76DB",
        alpha=0.65,
        edgecolor="white",
        label="Detailed-balance $\\pi_n$",
    )
    ax.plot(
        n_plot,
        poisson_pmf,
        color="#B60205",
        lw=1.8,
        label=f"Poisson($\\rho$ = {rho:.0f})",
    )
    ax.axvline(rho, color="black", ls="--", lw=1, alpha=0.6)
    ax.text(
        rho + 1,
        ax.get_ylim()[1] * 0.95,
        f"$\\rho$ = {rho:.0f}",
        fontsize=10,
    )
    ax.set_xlabel("Team size $n$")
    ax.set_ylabel("Stationary probability $\\pi_n$")
    ax.set_title(
        f"Birth-death stationary distribution\n"
        f"$\\lambda = {LAMBDA}$/month, $\\mu_n = n \\cdot {MU}$ per month — "
        f"matches Poisson($\\rho$ = {rho:.0f})"
    )
    ax.legend()
    fig.tight_layout()
    OUT_STATIONARY.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_STATIONARY, dpi=DPI, bbox_inches="tight")
    print(f"wrote {OUT_STATIONARY}")
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────────
# Figure 2 — Kolmogorov trajectory
# ──────────────────────────────────────────────────────────────────────────


def plot_trajectory(bd: BirthDeathProcess) -> None:
    times = np.linspace(0, 60, 61)
    rho = LAMBDA / MU

    # Exact transient distribution at each time t: P(t)[n0, :].
    states = np.arange(N_MAX + 1)
    means = np.empty_like(times)
    p05 = np.empty_like(times)
    p50 = np.empty_like(times)
    p95 = np.empty_like(times)

    for k, t in enumerate(times):
        if t == 0:
            row = np.zeros(N_MAX + 1)
            row[N0] = 1.0
        else:
            row = bd.transition_matrix(t)[N0]
        means[k] = float(np.sum(states * row))
        cum = np.cumsum(row)
        # Use generalized inverse of CDF to find quantiles
        p05[k] = float(states[np.searchsorted(cum, 0.05)])
        p50[k] = float(states[np.searchsorted(cum, 0.50)])
        p95[k] = float(states[np.searchsorted(cum, 0.95)])

    closed_form = bd.expected_trajectory(n0=N0, times=times)

    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.fill_between(
        times,
        p05,
        p95,
        color="#1D76DB",
        alpha=0.18,
        label="5–95% band (exact transient)",
    )
    ax.plot(
        times,
        means,
        color="#1D76DB",
        lw=1.8,
        label="$\\mathbb{E}[n(t)]$ from $e^{Qt}$",
    )
    ax.plot(
        times,
        closed_form,
        color="#B60205",
        lw=1.4,
        ls="--",
        label="$\\rho + (n_0 - \\rho) e^{-\\mu t}$ (ODE)",
    )
    ax.axhline(
        rho, color="black", ls=":", lw=1, alpha=0.6,
        label=f"Steady-state mean $\\rho$ = {rho:.0f}",
    )
    ax.scatter([0], [N0], color="black", zorder=5, label=f"$n_0$ = {N0}")

    # Annotation: half-life.
    half = np.log(2) / MU
    ax.axvline(half, color="grey", ls="--", lw=0.7, alpha=0.5)
    ax.text(
        half + 0.5,
        N0 + 1,
        f"half-life $\\ln 2 / \\mu \\approx {half:.1f}$ months",
        fontsize=9,
        color="grey",
    )

    ax.set_xlim(0, times[-1])
    ax.set_xlabel("$t$ (months)")
    ax.set_ylabel("Team size $n$")
    ax.set_title(
        "Kolmogorov trajectory — expected team size and 5–95% band\n"
        f"$\\lambda = {LAMBDA}$/month, $\\mu_n = n \\cdot {MU}$, $n_0 = {N0}$"
    )
    ax.legend(loc="lower right", framealpha=0.95)
    fig.tight_layout()
    OUT_TRAJECTORY.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_TRAJECTORY, dpi=DPI, bbox_inches="tight")
    print(f"wrote {OUT_TRAJECTORY}")
    plt.close(fig)


def main() -> None:
    sns.set_style("white")
    bd = make_chain()
    plot_stationary(bd)
    plot_trajectory(bd)


if __name__ == "__main__":
    main()
