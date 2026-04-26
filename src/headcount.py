"""Applied headcount model combining DTMC career progression with a
birth-death team-size process.

The model separates two questions:

- **Composition.** What fraction of the team sits at each career level? This
  is the discrete-time Markov chain on {Junior, Mid, Senior, Exit} from
  Phases 1–3.
- **Total size.** How many people are on the team at time t? This is the
  birth-death CTMC from Phase 4.

A `HeadcountModel` combines both. Methods answer the planning questions
listed in the roadmap:

- ``forecast(months)`` — distribution of team size at a future month
- ``prob_reach_target(target, deadline)`` — probability of n(t) >= target
- ``expected_time_to_target(target)`` — first time E[n(t)] crosses target
- ``steady_state_composition()`` — long-run % per career level
- ``simulate_trajectories(n_paths, months, seed)`` — sample paths
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from src.ctmc import BirthDeathProcess
from src.dtmc import MarkovChain

ArrayF = NDArray[np.float64]
ArrayI = NDArray[np.int64]


# ──────────────────────────────────────────────────────────────────────────
# Default parameters from the model design (docs/model-design.md)
# ──────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class DTMCParams:
    """Career-level transition rates (per month).

    The rate ``p_jm`` is the monthly probability a Junior is promoted to Mid;
    ``p_je`` is the monthly probability a Junior leaves the company; etc.
    The Exit row defines replacement-hiring behaviour: ``p_ej`` > 0 routes
    Exit back to Junior at that rate (recycling). For the absorbing variant
    set ``p_ej = 0`` and ``p_ee = 1``.
    """

    p_jm: float = 0.03  # Junior -> Mid promotion
    p_je: float = 0.04  # Junior attrition
    p_ms: float = 0.02  # Mid -> Senior promotion
    p_me: float = 0.02  # Mid attrition
    p_se: float = 0.01  # Senior attrition
    p_ej: float = 0.50  # Exit -> Junior (recycling). 0 for absorbing variant.

    def transition_matrix(self) -> ArrayF:
        """Build the 4x4 row-stochastic matrix on (J, M, S, E)."""
        p_ee = 1.0 - self.p_ej
        return np.array(
            [
                [1 - self.p_jm - self.p_je, self.p_jm, 0.0, self.p_je],
                [0.0, 1 - self.p_ms - self.p_me, self.p_ms, self.p_me],
                [0.0, 0.0, 1 - self.p_se, self.p_se],
                [self.p_ej, 0.0, 0.0, p_ee],
            ]
        )


@dataclass(frozen=True)
class BirthDeathParams:
    """Team-size birth-death parameters.

    Default M/M/inf: lambda hires per month and per-capita attrition mu.
    """

    lambda_rate: float = 2.0       # hires per month
    mu_rate: float = 0.025         # per-capita attrition per month
    n_max: int = 200               # generator truncation


# ──────────────────────────────────────────────────────────────────────────
# HeadcountModel
# ──────────────────────────────────────────────────────────────────────────


@dataclass
class HeadcountModel:
    """Applied workforce model wrapping a DTMC and a birth-death CTMC.

    The two processes are independent in this formulation: composition
    (career-level fractions) is governed by ``dtmc``; total team size is
    governed by ``birth_death``. The expected number of people at level k at
    time t is therefore ``E[n(t)] * fraction_k(t)`` where the second factor
    comes from the DTMC and the first from the birth-death.

    Attributes:
        dtmc_params: Career-level transition rates.
        birth_death_params: Hiring/attrition rates and truncation.
        n0: Initial team size (used for trajectory queries).
        labels: Career-level labels in the order (J, M, S, E).
    """

    dtmc_params: DTMCParams = field(default_factory=DTMCParams)
    birth_death_params: BirthDeathParams = field(default_factory=BirthDeathParams)
    n0: float = 45.0
    labels: tuple[str, ...] = ("Junior", "Mid", "Senior", "Exit")

    def __post_init__(self) -> None:
        self._dtmc = MarkovChain(
            self.dtmc_params.transition_matrix(),
            state_labels=list(self.labels),
        )
        self._birth_death = BirthDeathProcess.from_constant_hiring(
            lambda_rate=self.birth_death_params.lambda_rate,
            mu_rate=self.birth_death_params.mu_rate,
            n_max=self.birth_death_params.n_max,
        )

    # ──────────────────────────────────────────────────────────────────
    # Accessors
    # ──────────────────────────────────────────────────────────────────

    @property
    def dtmc(self) -> MarkovChain:
        return self._dtmc

    @property
    def birth_death(self) -> BirthDeathProcess:
        return self._birth_death

    # ──────────────────────────────────────────────────────────────────
    # Team size — birth-death view
    # ──────────────────────────────────────────────────────────────────

    def forecast(self, months: float) -> ArrayF:
        """Distribution over team sizes at month ``months`` given size ``n0``.

        Returns:
            1D array of length ``n_max + 1`` summing to 1.
        """
        if months < 0:
            raise ValueError("months must be non-negative.")
        n_max = self.birth_death_params.n_max
        n0_int = int(round(self.n0))
        if not (0 <= n0_int <= n_max):
            raise ValueError(f"n0 must be in [0, {n_max}].")
        if months == 0:
            row = np.zeros(n_max + 1)
            row[n0_int] = 1.0
            return row
        return self._birth_death.transition_matrix(months)[n0_int]

    def expected_team_size(self, months: float) -> float:
        """Closed-form expected team size at time ``months``."""
        lam = self.birth_death_params.lambda_rate
        mu = self.birth_death_params.mu_rate
        rho = lam / mu if mu > 0 else float("inf")
        if not np.isfinite(rho):
            return float(self.n0 + lam * months)
        return float(rho + (self.n0 - rho) * np.exp(-mu * months))

    def prob_reach_target(self, target: int, deadline: float) -> float:
        """Probability the team size is at least ``target`` at month ``deadline``."""
        if target < 0:
            raise ValueError("target must be non-negative.")
        dist = self.forecast(deadline)
        return float(dist[target:].sum())

    def expected_time_to_target(
        self, target: float, max_months: float = 1_000.0
    ) -> float:
        """Smallest t such that ``E[n(t)] >= target`` (or ``<= target`` if
        the target is below ``n0``).

        For the M/M/inf model this has a closed form. Returns +inf when the
        steady-state mean ``rho = lambda / mu`` cannot reach the target from
        ``n0``.
        """
        lam = self.birth_death_params.lambda_rate
        mu = self.birth_death_params.mu_rate
        if mu <= 0:
            # Pure-birth: n(t) = n0 + lambda * t
            if target < self.n0 and lam == 0:
                return 0.0
            if lam <= 0:
                return float("inf")
            return float(max(0.0, (target - self.n0) / lam))
        rho = lam / mu
        if self.n0 == target:
            return 0.0
        if (target - rho) * (self.n0 - rho) <= 0 or abs(target - rho) > abs(self.n0 - rho):
            # target lies between n0 and rho (reachable monotonically) OR
            # target is strictly outside the asymptote (not reachable).
            if abs(target - rho) >= abs(self.n0 - rho):
                # Trying to grow past the asymptote, or shrink past it.
                return float("inf")
        # E[n(t)] = rho + (n0 - rho) * exp(-mu t) -> solve for t.
        ratio = (target - rho) / (self.n0 - rho)
        if ratio <= 0:
            return float("inf")
        t = -np.log(ratio) / mu
        if t < 0 or t > max_months:
            return float("inf")
        return float(t)

    # ──────────────────────────────────────────────────────────────────
    # Composition — DTMC view
    # ──────────────────────────────────────────────────────────────────

    def composition_forecast(
        self, months: int, initial_composition: Sequence[float]
    ) -> ArrayF:
        """Career-level composition (J, M, S, E fractions) at month ``months``.

        Args:
            months: Non-negative integer number of months.
            initial_composition: Length-4 vector summing to 1.

        Returns:
            Length-4 vector summing to 1.
        """
        return self._dtmc.evolve(initial_composition, months)

    def steady_state_composition(self) -> dict[str, float]:
        """Long-run fraction of headcount-time at each career level."""
        pi = self._dtmc.stationary()
        return dict(zip(self.labels, pi.tolist(), strict=True))

    # ──────────────────────────────────────────────────────────────────
    # Simulation
    # ──────────────────────────────────────────────────────────────────

    def simulate_trajectories(
        self,
        n_paths: int,
        months: int,
        seed: int | None = None,
    ) -> ArrayI:
        """Simulate ``n_paths`` trajectories of total team size sampled at
        each integer month from 0 to ``months``.

        Uses Gillespie steps for the birth-death CTMC: at each integer month
        boundary the state is recorded.

        Args:
            n_paths: Number of independent paths.
            months: Horizon in integer months.
            seed: Optional RNG seed.

        Returns:
            Array of shape ``(n_paths, months + 1)`` with team size at each
            month boundary.
        """
        if n_paths <= 0:
            raise ValueError("n_paths must be positive.")
        if months < 0:
            raise ValueError("months must be non-negative.")
        rng = np.random.default_rng(seed)
        lam = self.birth_death_params.lambda_rate
        mu = self.birth_death_params.mu_rate
        n_max = self.birth_death_params.n_max
        n0 = int(round(self.n0))

        out = np.empty((n_paths, months + 1), dtype=np.int64)

        for k in range(n_paths):
            n = n0
            t = 0.0
            out[k, 0] = n
            recorded = 1
            while recorded <= months:
                rate_birth = lam if n < n_max else 0.0
                rate_death = mu * n
                total_rate = rate_birth + rate_death
                if total_rate <= 0:
                    # Absorbed at boundary; fill remainder with current state.
                    out[k, recorded:] = n
                    break
                dt = rng.exponential(1.0 / total_rate)
                t += dt
                # Record state at every integer month crossed before t.
                while recorded <= months and recorded <= t:
                    out[k, recorded] = n
                    recorded += 1
                if recorded > months:
                    break
                # Apply jump.
                if rng.random() < rate_birth / total_rate:
                    n = min(n + 1, n_max)
                else:
                    n = max(n - 1, 0)
        return out

    # ──────────────────────────────────────────────────────────────────
    # Budget bridge (link to Article 1)
    # ──────────────────────────────────────────────────────────────────

    def expected_total_salary(
        self,
        months: int,
        salary_mean: float,
        salary_var: float,
    ) -> tuple[float, float]:
        """Expected value and variance of the annual salary cost.

        Uses the law of total expectation/variance with N = headcount at the
        given month and i.i.d. salaries S_i with mean ``salary_mean`` and
        variance ``salary_var``:

            E[sum S_i over months M] = E[N] * E[S] * months
            Var[sum]                  = months^2 * (E[N]^2 Var[S]
                                                    + E[S]^2 Var[N]
                                                    + Var[N] Var[S])

        This connects Article 3 (headcount process) to Article 1
        (Monte Carlo budget simulation). The result is a closed-form analog
        of what Article 1 computes by simulation.

        Args:
            months: Number of months over which salaries accrue.
            salary_mean: ``E[S]`` for one salary.
            salary_var: ``Var[S]`` for one salary.

        Returns:
            Tuple ``(E[total], Var[total])``.
        """
        if months <= 0:
            raise ValueError("months must be positive.")
        if salary_mean < 0 or salary_var < 0:
            raise ValueError("salary moments must be non-negative.")
        dist = self.forecast(months)
        n_states = np.arange(dist.size)
        e_n = float((n_states * dist).sum())
        var_n = float(((n_states - e_n) ** 2 * dist).sum())
        e_total = e_n * salary_mean * months
        var_total = (months ** 2) * (
            e_n ** 2 * salary_var
            + salary_mean ** 2 * var_n
            + var_n * salary_var
        )
        return e_total, var_total
