"""Experiment D — simulation-based steady-state distribution.

Simulates many long trajectories of the M/M/inf headcount chain via
Gillespie steps, discards the burn-in transient, and pools the remaining
team-size samples. The empirical histogram is overlaid with the analytical
Poisson(rho) PMF and the detailed-balance distribution from
:meth:`BirthDeathProcess.stationary_birth_death` for cross-validation.

Generates ``figures/birth_death_steady_simulation.png``.

Run with:

    python scripts/exp_birth_death_simulation.py
"""

from __future__ import annotations

import sys
from pathlib import Path

# Make the project root importable when this script is invoked directly
# (no PYTHONPATH needed).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.stats import poisson  # noqa: E402

from scripts._style import DPI, PALETTE, SEED, apply  # noqa: E402
from src.headcount import BirthDeathParams, HeadcountModel  # noqa: E402

OUT = Path("figures/birth_death_steady_simulation.png")

# Run many independent paths. A burn-in of 5 / mu = 200 months reaches
# stationarity within machine precision; collect samples after that point.
N_PATHS = 200
BURN_IN_MONTHS = 200
SAMPLE_MONTHS = 200  # take one sample per month after burn-in
LAMBDA = 2.0
MU = 0.025


def main() -> None:
    apply()
    rho = LAMBDA / MU
    m = HeadcountModel(
        n0=int(round(rho)),  # start at the asymptote to minimise burn-in
        birth_death_params=BirthDeathParams(
            lambda_rate=LAMBDA, mu_rate=MU, n_max=200
        ),
    )
    paths = m.simulate_trajectories(
        n_paths=N_PATHS,
        months=BURN_IN_MONTHS + SAMPLE_MONTHS,
        seed=SEED,
    )
    samples = paths[:, BURN_IN_MONTHS:].ravel()

    # Analytical references.
    n_grid = np.arange(0, 130)
    poisson_pmf = poisson.pmf(n_grid, mu=rho)
    detailed_balance_pi = m.birth_death.stationary_birth_death()[: n_grid.size]

    fig, ax = plt.subplots(figsize=(9.5, 4.6))
    bins = np.arange(int(samples.min()), int(samples.max()) + 2) - 0.5
    ax.hist(
        samples,
        bins=bins,
        density=True,
        color=PALETTE["primary"],
        alpha=0.55,
        edgecolor="white",
        label=f"Simulated  ($N$ = {N_PATHS} paths × {SAMPLE_MONTHS} months)",
    )
    ax.plot(
        n_grid,
        detailed_balance_pi,
        color=PALETTE["tertiary"],
        lw=1.4,
        ls="--",
        label="Detailed-balance $\\pi_n$",
    )
    ax.plot(
        n_grid,
        poisson_pmf,
        color=PALETTE["accent"],
        lw=2.0,
        label=f"Poisson($\\rho$ = {rho:.0f})",
    )
    ax.axvline(rho, color=PALETTE["neutral"], ls=":", lw=1, alpha=0.6)
    ax.set_xlim(40, 130)
    ax.set_xlabel("Team size $n$")
    ax.set_ylabel("Density / probability")
    ax.set_title(
        "Experiment D — empirical steady-state team size matches "
        f"Poisson($\\rho$ = {rho:.0f})"
    )
    ax.legend(loc="upper right")

    # Print summary statistics for the commit log / notebook reference.
    print(f"empirical mean:     {samples.mean():.4f}  (target {rho:.0f})")
    print(f"empirical variance: {samples.var():.4f}  (target {rho:.0f})")
    ks = np.max(
        np.abs(np.cumsum(np.bincount(samples, minlength=n_grid.size)[: n_grid.size]
                         / samples.size) - poisson.cdf(n_grid, mu=rho))
    )
    print(f"sup |F_emp - F_Poisson|: {ks:.4f}")

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=DPI, bbox_inches="tight")
    print(f"wrote {OUT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
