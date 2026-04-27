"""Experiment G — animated team-size evolution.

Simulates 100 trajectories of the M/M/inf headcount chain over 36 months and
renders an animated GIF that reveals each successive month frame by frame:

- 100 grey sample paths grow forward in time
- the analytical mean line ``E[n(t)] = rho + (n0 - rho) e^{-mu t}`` overlays
  in blue
- the 5–95% band from the exact transient distribution shades the
  background

The GIF uses a coarse 12 fps × 37 frames so file size stays comfortably
under 5 MB on the project palette and figure size.

Run with:

    python scripts/exp_trajectory_animation.py
"""

from __future__ import annotations

import sys
from pathlib import Path

# Make the project root importable when this script is invoked directly
# (no PYTHONPATH needed).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib.animation as animation  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from scripts._style import PALETTE, SEED, apply  # noqa: E402
from src.headcount import BirthDeathParams, HeadcountModel  # noqa: E402

OUT = Path("figures/trajectory_animation.gif")
N_PATHS = 100
MONTHS = 36
N0 = 45
LAMBDA = 2.0
MU = 0.025


def transient_band(model: HeadcountModel, n0: int, times: np.ndarray):
    bd = model.birth_death
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


def main() -> None:
    apply()
    model = HeadcountModel(
        n0=N0,
        birth_death_params=BirthDeathParams(
            lambda_rate=LAMBDA, mu_rate=MU, n_max=200
        ),
    )
    times = np.arange(MONTHS + 1)
    paths = model.simulate_trajectories(
        n_paths=N_PATHS, months=MONTHS, seed=SEED
    )
    means, p05, p95 = transient_band(model, N0, times)

    rho = LAMBDA / MU

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.set_xlim(0, MONTHS)
    ax.set_ylim(min(p05.min(), paths.min()) - 2, max(p95.max(), paths.max()) + 2)
    ax.axhline(rho, color=PALETTE["neutral"], ls=":", lw=0.8, alpha=0.6,
               label=f"$\\rho$ = {rho:.0f}")
    ax.set_xlabel("$t$ (months)")
    ax.set_ylabel("Team size $n(t)$")
    ax.set_title(
        f"Headcount trajectory animation — {N_PATHS} paths, "
        f"$\\lambda$ = {LAMBDA}, $\\mu_n = n \\cdot {MU}$, $n_0 = {N0}$"
    )

    # Painted up front so the legend has a stable entry; the actual band is
    # redrawn in :func:`update` once we have real data for each frame.
    ax.fill_between([], [], [], color=PALETTE["primary"], alpha=0.18, label="5–95%")
    mean_line, = ax.plot([], [], color=PALETTE["primary"], lw=2.0,
                         label="$\\mathbb{E}[n(t)]$")
    path_lines = [
        ax.plot([], [], color="grey", alpha=0.30, lw=0.7)[0]
        for _ in range(N_PATHS)
    ]
    annotation = ax.text(
        0.02, 0.95, "", transform=ax.transAxes, fontsize=10,
        color=PALETTE["neutral"], va="top",
    )
    ax.legend(loc="lower right", framealpha=0.95)

    def init():  # type: ignore[no-untyped-def]
        for line in path_lines:
            line.set_data([], [])
        mean_line.set_data([], [])
        annotation.set_text("")
        return [mean_line, annotation, *path_lines]

    def update(frame: int):  # type: ignore[no-untyped-def]
        sub_t = times[: frame + 1]
        # Refresh band by removing the previous PolyCollection and redrawing.
        for coll in list(ax.collections):
            coll.remove()
        ax.fill_between(
            sub_t, p05[: frame + 1], p95[: frame + 1],
            color=PALETTE["primary"], alpha=0.18,
        )
        for k, line in enumerate(path_lines):
            line.set_data(sub_t, paths[k, : frame + 1])
        mean_line.set_data(sub_t, means[: frame + 1])
        annotation.set_text(
            f"$t$ = {frame:>2} months\n"
            f"$\\mathbb{{E}}[n]$ = {means[frame]:.1f}\n"
            f"5–95% = [{p05[frame]:.0f}, {p95[frame]:.0f}]"
        )
        return [mean_line, annotation, *path_lines]

    anim = animation.FuncAnimation(
        fig, update, init_func=init,
        frames=MONTHS + 1, interval=200, blit=False, repeat=False,
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    # PillowWriter keeps the file size small for the chosen colour count.
    writer = animation.PillowWriter(fps=6)
    anim.save(OUT, writer=writer, dpi=120)
    plt.close(fig)
    size_kb = OUT.stat().st_size / 1024
    print(f"wrote {OUT}  ({size_kb:.1f} KB)")
    if size_kb > 5_000:
        raise RuntimeError(
            f"GIF exceeds 5 MB target ({size_kb / 1024:.2f} MB) — "
            "lower fps, frames, dpi, or path count."
        )


if __name__ == "__main__":
    main()
