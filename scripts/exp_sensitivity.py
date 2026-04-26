"""Sensitivity tornado chart for Senior steady-state share.

Vary each rate in the headcount model by ±50% from its base value, recompute
the Senior fraction (DTMC) and the steady-state mean team size (birth-death),
and rank parameters by the magnitude of the impact.

Generates ``figures/sensitivity_tornado.png``.

Run with:

    python scripts/exp_sensitivity.py
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from src.headcount import BirthDeathParams, DTMCParams, HeadcountModel

OUT = Path("figures/sensitivity_tornado.png")
DPI = 300

# Base case
BASE_DTMC = DTMCParams()
BASE_BD = BirthDeathParams()

# Parameters to vary, with display labels.
DTMC_PARAMS: list[tuple[str, str]] = [
    ("p_jm", "$p_{12}$ (J->M promotion)"),
    ("p_je", "$p_{14}$ (J attrition)"),
    ("p_ms", "$p_{23}$ (M->S promotion)"),
    ("p_me", "$p_{24}$ (M attrition)"),
    ("p_se", "$p_{34}$ (S attrition)"),
]

BD_PARAMS: list[tuple[str, str]] = [
    ("lambda_rate", "$\\lambda$ (hiring rate)"),
    ("mu_rate", "$\\mu$ (per-capita attrition)"),
]


def vary_dtmc(field: str, factor: float) -> float:
    """Return Senior steady-state fraction with ``field`` scaled by ``factor``."""
    base_value = float(getattr(BASE_DTMC, field))
    new_value = base_value * factor
    new = replace(BASE_DTMC, **{field: new_value})
    m = HeadcountModel(dtmc_params=new)
    return m.steady_state_composition()["Senior"]


def vary_bd(field: str, factor: float) -> float:
    """Return steady-state mean team size with ``field`` scaled by ``factor``."""
    base_value = float(getattr(BASE_BD, field))
    new_value = base_value * factor
    new = replace(BASE_BD, **{field: new_value})
    return new.lambda_rate / new.mu_rate


def main() -> None:
    sns.set_style("white")
    base_senior = HeadcountModel().steady_state_composition()["Senior"]
    base_rho = BASE_BD.lambda_rate / BASE_BD.mu_rate

    fig, (ax_senior, ax_rho) = plt.subplots(1, 2, figsize=(13, 4.8))

    # ── Tornado for Senior share ──
    senior_rows = []
    for field, label in DTMC_PARAMS:
        low = vary_dtmc(field, 0.5)
        high = vary_dtmc(field, 1.5)
        senior_rows.append((label, low - base_senior, high - base_senior))
    senior_rows.sort(key=lambda r: -max(abs(r[1]), abs(r[2])))

    y = np.arange(len(senior_rows))
    for i, (label, dlow, dhigh) in enumerate(senior_rows):
        # Plot the bar from min to max delta around the base.
        ax_senior.barh(
            i, dlow, left=0, color="#B60205", alpha=0.75, edgecolor="black",
            label="−50%" if i == 0 else None,
        )
        ax_senior.barh(
            i, dhigh, left=0, color="#0E8A16", alpha=0.75, edgecolor="black",
            label="+50%" if i == 0 else None,
        )
        ax_senior.text(
            dlow - 0.005,
            i,
            f"{base_senior + dlow:.2f}",
            va="center",
            ha="right",
            fontsize=9,
            color="#B60205",
        )
        ax_senior.text(
            dhigh + 0.005,
            i,
            f"{base_senior + dhigh:.2f}",
            va="center",
            ha="left",
            fontsize=9,
            color="#0E8A16",
        )

    ax_senior.axvline(0, color="black", lw=0.8)
    ax_senior.set_yticks(y)
    ax_senior.set_yticklabels([row[0] for row in senior_rows])
    ax_senior.invert_yaxis()
    ax_senior.set_xlabel("Δ Senior $\\pi_S$  (base = " + f"{base_senior:.3f})")
    ax_senior.set_title("Senior fraction sensitivity (DTMC rates ±50%)")
    ax_senior.legend(loc="lower right", framealpha=0.95)
    ax_senior.grid(axis="x", alpha=0.3)

    # ── Tornado for steady-state mean team size ρ = λ/μ ──
    rho_rows = []
    for field, label in BD_PARAMS:
        low = vary_bd(field, 0.5)
        high = vary_bd(field, 1.5)
        rho_rows.append((label, low - base_rho, high - base_rho))
    rho_rows.sort(key=lambda r: -max(abs(r[1]), abs(r[2])))

    y = np.arange(len(rho_rows))
    for i, (label, dlow, dhigh) in enumerate(rho_rows):
        ax_rho.barh(
            i, dlow, left=0, color="#B60205", alpha=0.75, edgecolor="black",
            label="−50%" if i == 0 else None,
        )
        ax_rho.barh(
            i, dhigh, left=0, color="#0E8A16", alpha=0.75, edgecolor="black",
            label="+50%" if i == 0 else None,
        )
        ax_rho.text(dlow - 1.5, i, f"{base_rho + dlow:.0f}", va="center",
                    ha="right", fontsize=9, color="#B60205")
        ax_rho.text(dhigh + 1.5, i, f"{base_rho + dhigh:.0f}", va="center",
                    ha="left", fontsize=9, color="#0E8A16")
    ax_rho.axvline(0, color="black", lw=0.8)
    ax_rho.set_yticks(y)
    ax_rho.set_yticklabels([row[0] for row in rho_rows])
    ax_rho.invert_yaxis()
    ax_rho.set_xlabel(f"Δ steady-state mean team size  (base = {base_rho:.0f})")
    ax_rho.set_title("Team-size sensitivity (birth-death rates ±50%)")
    ax_rho.legend(loc="lower right", framealpha=0.95)
    ax_rho.grid(axis="x", alpha=0.3)

    fig.suptitle(
        "Sensitivity analysis — which rates move the steady state?",
        fontsize=12,
        y=1.02,
    )
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=DPI, bbox_inches="tight")
    print(f"wrote {OUT}")
    print()
    print(f"base Senior fraction = {base_senior:.4f}")
    for label, dlow, dhigh in senior_rows:
        print(f"  {label:<32} +/- {max(abs(dlow), abs(dhigh)):.4f}")
    print(f"\nbase rho = {base_rho:.0f}")
    for label, dlow, dhigh in rho_rows:
        print(f"  {label:<32} range [{base_rho + dlow:.1f}, {base_rho + dhigh:.1f}]")


if __name__ == "__main__":
    main()
