"""Generate one figure per applied headcount scenario (S1–S5).

Each figure shows the analytical mean trajectory, an exact 5–95% band
derived from the truncated CTMC's transient distribution, and 30 simulated
sample paths. All scripts are deterministic given the seed.

Run with:

    python scripts/exp_headcount_scenarios.py
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from src.ctmc import BirthDeathProcess
from src.headcount import BirthDeathParams, DTMCParams, HeadcountModel

SEED = 2026
DPI = 300
N_SIM_PATHS = 30
FIG_DIR = Path("figures")


# ──────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────


def transient_band(
    bd: BirthDeathProcess, n0: int, times: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Exact mean and 5/95 quantiles from the CTMC's transient distribution."""
    states = np.arange(bd.n_states)
    means = np.empty_like(times, dtype=float)
    p05 = np.empty_like(times, dtype=float)
    p95 = np.empty_like(times, dtype=float)
    for k, t in enumerate(times):
        if t == 0:
            row = np.zeros(bd.n_states)
            row[n0] = 1.0
        else:
            row = bd.transition_matrix(float(t))[n0]
        means[k] = float((states * row).sum())
        cum = np.cumsum(row)
        p05[k] = float(states[np.searchsorted(cum, 0.05)])
        p95[k] = float(states[np.searchsorted(cum, 0.95)])
    return means, p05, p95


def plot_scenario(
    title: str,
    subtitle: str,
    out_path: Path,
    times: np.ndarray,
    paths: np.ndarray,
    mean: np.ndarray,
    p05: np.ndarray,
    p95: np.ndarray,
    annotate: Callable[[plt.Axes], None] | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.fill_between(
        times, p05, p95, color="#1D76DB", alpha=0.18,
        label="5–95% (exact transient)",
    )
    for path in paths:
        ax.plot(times, path, color="grey", alpha=0.25, lw=0.7)
    ax.plot(times, mean, color="#1D76DB", lw=2.0, label="$\\mathbb{E}[n(t)]$")
    ax.scatter([times[0]], [paths[0, 0]], color="black", zorder=5,
               label=f"$n_0$ = {int(paths[0, 0])}")

    if annotate is not None:
        annotate(ax)

    ax.set_xlim(times[0], times[-1])
    ax.set_xlabel("$t$ (months)")
    ax.set_ylabel("Team size $n(t)$")
    ax.set_title(f"{title}\n{subtitle}", fontsize=11)
    ax.legend(loc="best", framealpha=0.95)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    print(f"wrote {out_path}")
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────────
# Scenarios
# ──────────────────────────────────────────────────────────────────────────


def scenario_s1_steady_growth() -> None:
    """S1: target 60 by month 12. Compare base vs boosted lambda."""
    months = 24
    times = np.arange(months + 1)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    for ax, (title, lam) in zip(
        axes, [("Base ($\\lambda$ = 2.0)", 2.0), ("Boosted ($\\lambda$ = 2.31)", 2.31)],
        strict=True,
    ):
        m = HeadcountModel(
            n0=45,
            birth_death_params=BirthDeathParams(lambda_rate=lam, mu_rate=0.025),
        )
        paths = m.simulate_trajectories(N_SIM_PATHS, months, seed=SEED)
        mean, p05, p95 = transient_band(m.birth_death, 45, times)
        ax.fill_between(times, p05, p95, color="#1D76DB", alpha=0.18)
        for path in paths:
            ax.plot(times, path, color="grey", alpha=0.25, lw=0.7)
        ax.plot(times, mean, color="#1D76DB", lw=2.0, label="$\\mathbb{E}[n(t)]$")
        ax.axhline(60, color="#B60205", ls="--", lw=1, label="target = 60")
        ax.axvline(12, color="grey", ls=":", lw=0.8)
        prob = m.prob_reach_target(60, 12)
        ax.set_title(f"{title}  |  $P(n(12) \\geq 60)$ = {prob:.2f}")
        ax.set_xlabel("$t$ (months)")
        ax.set_ylabel("Team size")
        ax.set_xlim(0, months)
        ax.legend(loc="lower right")
    fig.suptitle("S1 — Steady growth: hitting 60 by month 12", fontsize=12, y=1.02)
    fig.tight_layout()
    out = FIG_DIR / "scenario_s1_steady_growth.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    print(f"wrote {out}")
    plt.close(fig)


def scenario_s2_hiring_freeze() -> None:
    """S2: lambda = 0; pure exponential decay."""
    months = 36
    times = np.arange(months + 1)
    m = HeadcountModel(
        n0=45,
        birth_death_params=BirthDeathParams(lambda_rate=0.0, mu_rate=0.025),
    )
    paths = m.simulate_trajectories(N_SIM_PATHS, months, seed=SEED)
    mean, p05, p95 = transient_band(m.birth_death, 45, times)

    def annotate(ax: plt.Axes) -> None:
        # Threshold crossing at 40
        t_cross = -np.log(40.0 / 45.0) / 0.025
        ax.axhline(40, color="#B60205", ls="--", lw=1, label="threshold = 40")
        ax.axvline(t_cross, color="#B60205", ls=":", lw=0.8)
        ax.text(t_cross + 0.3, 41, f"$\\mathbb{{E}}[n] = 40$ at $t \\approx {t_cross:.1f}$",
                fontsize=9, color="#B60205")

    plot_scenario(
        title="S2 — Hiring freeze",
        subtitle="$\\lambda = 0$, $\\mu = 0.025$ — exponential decay $\\mathbb{E}[n(t)] = 45 e^{-0.025 t}$",
        out_path=FIG_DIR / "scenario_s2_hiring_freeze.png",
        times=times,
        paths=paths,
        mean=mean,
        p05=p05,
        p95=p95,
        annotate=annotate,
    )


def scenario_s3_layoff_recovery() -> None:
    """S3: start at n0 = 35, recover toward steady state."""
    months = 36
    times = np.arange(months + 1)
    m = HeadcountModel(n0=35)
    paths = m.simulate_trajectories(N_SIM_PATHS, months, seed=SEED)
    mean, p05, p95 = transient_band(m.birth_death, 35, times)

    t_recover = m.expected_time_to_target(45)

    def annotate(ax: plt.Axes) -> None:
        ax.axhline(45, color="#B60205", ls="--", lw=1, label="pre-layoff = 45")
        ax.axhline(80, color="black", ls=":", lw=0.8, label="$\\rho$ = 80")
        ax.axvline(t_recover, color="#B60205", ls=":", lw=0.8)
        ax.text(
            t_recover + 0.3,
            46,
            f"$\\mathbb{{E}}[n] = 45$ at $t \\approx {t_recover:.1f}$",
            fontsize=9,
            color="#B60205",
        )

    plot_scenario(
        title="S3 — Layoff recovery",
        subtitle="$n_0$ = 35 after a cut; $\\lambda$, $\\mu$ unchanged — pulled back toward $\\rho = 80$",
        out_path=FIG_DIR / "scenario_s3_layoff_recovery.png",
        times=times,
        paths=paths,
        mean=mean,
        p05=p05,
        p95=p95,
        annotate=annotate,
    )


def scenario_s4_composition_shift() -> None:
    """S4: vary Senior attrition; show effect on stationary composition."""
    p34_grid = np.linspace(0.005, 0.02, 9)
    senior_share = []
    junior_share = []
    mid_share = []
    for p34 in p34_grid:
        m = HeadcountModel(dtmc_params=DTMCParams(p_se=float(p34)))
        comp = m.steady_state_composition()
        senior_share.append(comp["Senior"])
        junior_share.append(comp["Junior"])
        mid_share.append(comp["Mid"])

    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.plot(p34_grid, senior_share, marker="o", color="#0E8A16", lw=2,
            label="Senior $\\pi_S$")
    ax.plot(p34_grid, junior_share, marker="s", color="#1D76DB", lw=1.5,
            label="Junior $\\pi_J$")
    ax.plot(p34_grid, mid_share, marker="^", color="#5319E7", lw=1.5,
            label="Mid $\\pi_M$")
    ax.axvline(0.01, color="grey", ls=":", lw=0.8)
    ax.text(0.01 + 0.0003, 0.05, "base $p_{34}$ = 0.01", fontsize=9, color="grey")
    ax.set_xlabel("Senior attrition rate $p_{34}$")
    ax.set_ylabel("Stationary fraction $\\pi_i$")
    ax.set_title("S4 — Composition shift via Senior attrition")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = FIG_DIR / "scenario_s4_composition_shift.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    print(f"wrote {out}")
    plt.close(fig)


def scenario_s5_seasonal_hiring() -> None:
    """S5: piecewise lambda(t). Compose CTMC step-wise."""
    months = 24
    n0 = 45
    n_max = 200
    mu = 0.025

    high, low = 3.0, 1.0
    quarters = []
    for q in range(months):
        # Quarter boundaries: months 0-2 = Q1 (high), 3-5 Q2 (low),
        # 6-8 Q3 (high), 9-11 Q4 (low), repeat next year.
        season = (q // 3) % 4
        lam = high if season in (0, 2) else low
        quarters.append(lam)

    # Compose piecewise CTMC matrices to get distribution at each month.
    states = np.arange(n_max + 1)
    e_initial = np.zeros(n_max + 1)
    e_initial[n0] = 1.0
    pi_t = e_initial.copy()
    means = [n0]
    p05 = [n0]
    p95 = [n0]
    for q in range(months):
        bd_q = BirthDeathProcess.from_constant_hiring(quarters[q], mu, n_max=n_max)
        pi_t = pi_t @ bd_q.transition_matrix(1.0)
        means.append(float((states * pi_t).sum()))
        cum = np.cumsum(pi_t)
        p05.append(float(states[np.searchsorted(cum, 0.05)]))
        p95.append(float(states[np.searchsorted(cum, 0.95)]))
    times = np.arange(months + 1)
    means = np.array(means)
    p05 = np.array(p05)
    p95 = np.array(p95)

    # Compare with constant-lambda baseline.
    bd_const = BirthDeathProcess.from_constant_hiring(2.0, mu, n_max=n_max)
    means_const, p05_const, p95_const = transient_band(bd_const, n0, times)

    fig, ax = plt.subplots(figsize=(9, 4.8))
    # Background bands for quarters
    for q in range(months):
        season = (q // 3) % 4
        if season in (0, 2):
            ax.axvspan(q, q + 1, color="#0E8A16", alpha=0.06)
        else:
            ax.axvspan(q, q + 1, color="#B60205", alpha=0.06)

    ax.fill_between(times, p05, p95, color="#1D76DB", alpha=0.18,
                    label="5–95% seasonal")
    ax.plot(times, means, color="#1D76DB", lw=2,
            label="$\\mathbb{E}[n(t)]$ seasonal $\\bar\\lambda = 2$")
    ax.plot(times, means_const, color="black", lw=1.4, ls="--",
            label="constant $\\lambda = 2$")
    ax.set_xlim(0, months)
    ax.set_xlabel("$t$ (months)")
    ax.set_ylabel("Team size $n(t)$")
    ax.set_title("S5 — Seasonal hiring  ($\\lambda$ = 3 in Q1/Q3, 1 in Q2/Q4)")
    ax.legend(loc="lower right")
    fig.tight_layout()
    out = FIG_DIR / "scenario_s5_seasonal_hiring.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    print(f"wrote {out}")
    plt.close(fig)


def main() -> None:
    sns.set_style("white")
    scenario_s1_steady_growth()
    scenario_s2_hiring_freeze()
    scenario_s3_layoff_recovery()
    scenario_s4_composition_shift()
    scenario_s5_seasonal_hiring()


if __name__ == "__main__":
    main()
