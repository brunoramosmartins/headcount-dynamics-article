"""Tests for src/headcount.py — HeadcountModel."""

from __future__ import annotations

import numpy as np
import pytest

from src.headcount import BirthDeathParams, DTMCParams, HeadcountModel


# Construction & defaults ────────────────────────────────────────────────


class TestConstruction:
    def test_defaults(self) -> None:
        m = HeadcountModel()
        assert m.n0 == 45.0
        assert m.dtmc.n_states == 4
        assert m.birth_death.n_states == 201

    def test_dtmc_is_row_stochastic(self) -> None:
        m = HeadcountModel()
        np.testing.assert_allclose(m.dtmc.P.sum(axis=1), np.ones(4), atol=1e-12)

    def test_dtmc_recycling_default(self) -> None:
        # Default p_ej = 0.5, so the chain is irreducible.
        m = HeadcountModel()
        assert m.dtmc.is_irreducible() is True

    def test_dtmc_absorbing_when_p_ej_zero(self) -> None:
        m = HeadcountModel(dtmc_params=DTMCParams(p_ej=0.0))
        assert m.dtmc.is_irreducible() is False
        assert m.dtmc.is_absorbing_state(3) is True


# Forecasts ──────────────────────────────────────────────────────────────


class TestForecast:
    def test_forecast_zero_is_concentrated(self) -> None:
        m = HeadcountModel(n0=45)
        dist = m.forecast(0)
        assert dist[45] == 1.0
        assert dist.sum() == 1.0

    def test_forecast_sums_to_one(self) -> None:
        m = HeadcountModel()
        dist = m.forecast(12)
        assert np.isclose(dist.sum(), 1.0)
        assert (dist >= 0).all()

    def test_expected_team_size_closed_form(self) -> None:
        m = HeadcountModel(n0=45)
        # E[n(12)] = 80 + (45 - 80) * exp(-0.025*12)
        expected = 80.0 + (45.0 - 80.0) * np.exp(-0.3)
        assert np.isclose(m.expected_team_size(12), expected)

    def test_expected_team_size_at_zero(self) -> None:
        m = HeadcountModel(n0=45)
        assert m.expected_team_size(0) == 45

    def test_pure_birth_extreme(self) -> None:
        # mu = 0 -> pure-birth: E[n(t)] = n0 + lambda * t, but the chain is
        # not positive recurrent. The model still computes E[n(t)] correctly.
        m = HeadcountModel(
            birth_death_params=BirthDeathParams(lambda_rate=2.0, mu_rate=1e-12, n_max=200),
            n0=45,
        )
        # With essentially zero attrition, mean should drift roughly linearly
        # over short horizons.
        size_at_5 = m.expected_team_size(5)
        assert 50 < size_at_5 < 60


class TestProbReachTarget:
    def test_probabilities_in_zero_one(self) -> None:
        m = HeadcountModel(n0=45)
        for target in [40, 50, 60, 80]:
            p = m.prob_reach_target(target, deadline=12)
            assert 0.0 <= p <= 1.0

    def test_motivating_question(self) -> None:
        # P(X_12 >= 50 | X_0 = 45) for default params; matches Phase 4 notes.
        m = HeadcountModel(n0=45)
        p = m.prob_reach_target(50, deadline=12)
        assert 0.78 < p < 0.82

    def test_very_high_target_is_unlikely(self) -> None:
        m = HeadcountModel(n0=45)
        assert m.prob_reach_target(199, deadline=12) < 1e-6


class TestExpectedTimeToTarget:
    def test_target_at_n0_is_zero(self) -> None:
        m = HeadcountModel(n0=45)
        assert m.expected_time_to_target(45) == 0.0

    def test_target_at_rho_is_infinite(self) -> None:
        # n0 = 45, rho = 80. Mean trajectory approaches rho asymptotically
        # but never reaches it.
        m = HeadcountModel(n0=45)
        assert np.isinf(m.expected_time_to_target(80))

    def test_target_above_rho_is_infinite(self) -> None:
        m = HeadcountModel(n0=45)
        assert np.isinf(m.expected_time_to_target(85))

    def test_target_70(self) -> None:
        m = HeadcountModel(n0=45)
        # 80 - 35 * exp(-0.025 t) = 70 -> t = ln(35/10) / 0.025
        expected = np.log(35.0 / 10.0) / 0.025
        assert np.isclose(m.expected_time_to_target(70), expected)


# Steady-state composition ───────────────────────────────────────────────


class TestSteadyStateComposition:
    def test_keys_match_labels(self) -> None:
        m = HeadcountModel()
        comp = m.steady_state_composition()
        assert set(comp.keys()) == {"Junior", "Mid", "Senior", "Exit"}

    def test_sums_to_one(self) -> None:
        m = HeadcountModel()
        comp = m.steady_state_composition()
        assert np.isclose(sum(comp.values()), 1.0)

    def test_matches_phase2_hand_solution(self) -> None:
        # Phase 2: pi_J = 1/3.39, pi_M = 0.75 pi_J, pi_S = 1.5 pi_J,
        # pi_E = 0.14 pi_J.
        m = HeadcountModel()
        comp = m.steady_state_composition()
        pi_J = 1.0 / 3.39
        assert np.isclose(comp["Junior"], pi_J, atol=1e-3)
        assert np.isclose(comp["Senior"], 1.5 * pi_J, atol=1e-3)


class TestCompositionForecast:
    def test_evolves_correctly(self) -> None:
        m = HeadcountModel()
        pi0 = np.array([1.0, 0.0, 0.0, 0.0])
        pi1 = m.composition_forecast(1, pi0)
        np.testing.assert_allclose(pi1, np.array([0.93, 0.03, 0.0, 0.04]), atol=1e-12)

    def test_long_run_converges_to_stationary(self) -> None:
        m = HeadcountModel()
        pi0 = np.array([0.0, 0.0, 0.0, 1.0])  # all in Exit
        pi_long = m.composition_forecast(2_000, pi0)
        stationary = np.array(list(m.steady_state_composition().values()))
        np.testing.assert_allclose(pi_long, stationary, atol=1e-6)


# Simulation ─────────────────────────────────────────────────────────────


class TestSimulation:
    def test_shape(self) -> None:
        m = HeadcountModel()
        paths = m.simulate_trajectories(n_paths=20, months=12, seed=2026)
        assert paths.shape == (20, 13)

    def test_starts_at_n0(self) -> None:
        m = HeadcountModel(n0=45)
        paths = m.simulate_trajectories(n_paths=10, months=6, seed=0)
        assert (paths[:, 0] == 45).all()

    def test_states_in_valid_range(self) -> None:
        m = HeadcountModel()
        paths = m.simulate_trajectories(n_paths=50, months=24, seed=42)
        assert paths.min() >= 0
        assert paths.max() <= m.birth_death_params.n_max

    def test_seed_reproducibility(self) -> None:
        m = HeadcountModel()
        a = m.simulate_trajectories(n_paths=5, months=12, seed=7)
        b = m.simulate_trajectories(n_paths=5, months=12, seed=7)
        np.testing.assert_array_equal(a, b)

    def test_empirical_mean_close_to_analytical(self) -> None:
        m = HeadcountModel(n0=45)
        paths = m.simulate_trajectories(n_paths=2_000, months=12, seed=2026)
        empirical_means = paths.mean(axis=0)
        analytical_means = np.array(
            [m.expected_team_size(t) for t in range(13)]
        )
        # Standard error scales like sqrt(rho/n_paths) ~ 0.2; allow 0.5.
        np.testing.assert_allclose(empirical_means, analytical_means, atol=0.6)

    def test_invalid_n_paths(self) -> None:
        m = HeadcountModel()
        with pytest.raises(ValueError, match="n_paths"):
            m.simulate_trajectories(n_paths=0, months=1, seed=0)


# Scenarios ──────────────────────────────────────────────────────────────


class TestScenarios:
    def test_s1_growth_target(self) -> None:
        """Increasing lambda should push prob(reach 60 by month 12) up."""
        base = HeadcountModel(n0=45)
        boosted = HeadcountModel(
            n0=45,
            birth_death_params=BirthDeathParams(lambda_rate=2.4, mu_rate=0.025),
        )
        assert boosted.prob_reach_target(60, 12) > base.prob_reach_target(60, 12)

    def test_s2_freeze_decays(self) -> None:
        """With lambda = 0, expected size decays exponentially."""
        m = HeadcountModel(
            n0=45,
            birth_death_params=BirthDeathParams(lambda_rate=0.0, mu_rate=0.025),
        )
        assert m.expected_team_size(12) < 45
        # n0 * exp(-mu * t)
        assert np.isclose(m.expected_team_size(12), 45 * np.exp(-0.3), atol=1e-9)

    def test_s3_layoff_recovers(self) -> None:
        """With n0 = 35 < rho, mean grows back up."""
        m = HeadcountModel(n0=35)
        assert m.expected_team_size(12) > 35
        # E[time to 45] solvable in closed form
        t = m.expected_time_to_target(45)
        # ln(45/35) / 0.025 = 10.04 ish
        assert 9.5 < t < 10.5


# Budget bridge ──────────────────────────────────────────────────────────


class TestBudgetBridge:
    def test_basic_moments(self) -> None:
        m = HeadcountModel(n0=45)
        e_total, var_total = m.expected_total_salary(
            months=12,
            salary_mean=10_000.0,
            salary_var=4_000_000.0,  # std = 2000
        )
        # E[N] ~ 54, so E[total] ~ 54 * 10000 * 12 = 6.48M
        assert 6.0e6 < e_total < 7.0e6
        assert var_total > 0

    def test_zero_salary_var(self) -> None:
        # When salaries are deterministic, only headcount variance matters.
        m = HeadcountModel(n0=45)
        e_total, var_total = m.expected_total_salary(
            months=12, salary_mean=10_000.0, salary_var=0.0
        )
        # Var[total] = 12^2 * E[S]^2 * Var[N]
        dist = m.forecast(12)
        states = np.arange(dist.size)
        e_n = float((states * dist).sum())
        var_n = float(((states - e_n) ** 2 * dist).sum())
        expected_var = 144 * (10_000 ** 2) * var_n
        assert np.isclose(var_total, expected_var, rtol=1e-9)

    def test_invalid_months(self) -> None:
        m = HeadcountModel()
        with pytest.raises(ValueError, match="months"):
            m.expected_total_salary(months=0, salary_mean=1.0, salary_var=0.0)
