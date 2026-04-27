"""Run every experiment script in order and report figure sizes.

Each experiment is invoked as a Python module so a single seed and style
module are shared across all of them. Run from the project root:

    python scripts/run_all_experiments.py

The script halts on the first failure and prints a summary table at the end
(experiment ID, figure path, file size in KB).
"""

from __future__ import annotations

import importlib
import sys
import time
from collections.abc import Callable
from pathlib import Path

# Make the project root importable when this script is invoked as
# `python scripts/run_all_experiments.py` (no PYTHONPATH needed).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

EXPERIMENTS: list[tuple[str, str, str, list[Path]]] = [
    (
        "A",
        "Transition convergence",
        "scripts.exp_transition_convergence",
        [Path("figures/transition_convergence_heatmap.png")],
    ),
    (
        "B",
        "Eigenvalue spectrum",
        "scripts.exp_eigenvalue_analysis",
        [Path("figures/eigenvalue_spectrum.png")],
    ),
    (
        "C",
        "Absorption times",
        "scripts.exp_absorption",
        [Path("figures/absorption_times.png")],
    ),
    (
        "D-analytical",
        "Birth-death stationary (analytical)",
        "scripts.exp_birth_death",
        [
            Path("figures/birth_death_stationary.png"),
            Path("figures/kolmogorov_trajectory.png"),
        ],
    ),
    (
        "D-simulated",
        "Birth-death steady state (simulation histogram)",
        "scripts.exp_birth_death_simulation",
        [Path("figures/birth_death_steady_simulation.png")],
    ),
    (
        "E",
        "Headcount scenarios",
        "scripts.exp_headcount_scenarios",
        [
            Path("figures/scenario_s1_steady_growth.png"),
            Path("figures/scenario_s2_hiring_freeze.png"),
            Path("figures/scenario_s3_layoff_recovery.png"),
            Path("figures/scenario_s4_composition_shift.png"),
            Path("figures/scenario_s5_seasonal_hiring.png"),
        ],
    ),
    (
        "F",
        "Sensitivity tornado",
        "scripts.exp_sensitivity",
        [Path("figures/sensitivity_tornado.png")],
    ),
    (
        "G",
        "Trajectory animation",
        "scripts.exp_trajectory_animation",
        [Path("figures/trajectory_animation.gif")],
    ),
]


def _resolve_main(module_name: str) -> Callable[[], None]:
    """Import a script module and return its ``main`` callable."""
    module = importlib.import_module(module_name)
    if not hasattr(module, "main"):
        raise AttributeError(f"{module_name} has no `main()` function.")
    return module.main


def main() -> None:
    rows: list[tuple[str, str, str, str]] = []
    print("Running all Phase 6 experiments\n" + "=" * 50)

    for exp_id, description, module_name, outputs in EXPERIMENTS:
        print(f"\n[{exp_id}] {description}")
        print(f"    module: {module_name}")
        t0 = time.perf_counter()
        run = _resolve_main(module_name)
        run()
        elapsed = time.perf_counter() - t0
        print(f"    elapsed: {elapsed:.2f}s")
        for fig in outputs:
            if not fig.exists():
                raise FileNotFoundError(
                    f"expected output {fig} from experiment {exp_id} was not produced"
                )
            size_kb = fig.stat().st_size / 1024
            rows.append((exp_id, description, str(fig), f"{size_kb:,.1f} KB"))

    print("\n" + "=" * 50)
    print("Summary")
    print("=" * 50)
    width_id = max(len(r[0]) for r in rows)
    width_fig = max(len(r[2]) for r in rows)
    for exp_id, _description, fig_path, size in rows:
        print(f"  [{exp_id:<{width_id}}]  {fig_path:<{width_fig}}  {size}")


if __name__ == "__main__":
    main()
