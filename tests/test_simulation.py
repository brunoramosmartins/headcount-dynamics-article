"""Tests for src/simulation.py — path simulation."""

from __future__ import annotations

import numpy as np
import pytest

from src.dtmc import MarkovChain
from src.simulation import simulate_path, simulate_paths, state_distribution_at


@pytest.fixture
def chain() -> MarkovChain:
    P = np.array(
        [
            [0.93, 0.03, 0.00, 0.04],
            [0.00, 0.96, 0.02, 0.02],
            [0.00, 0.00, 0.99, 0.01],
            [0.50, 0.00, 0.00, 0.50],
        ]
    )
    return MarkovChain(P, state_labels=["Junior", "Mid", "Senior", "Exit"])


# simulate_path ──────────────────────────────────────────────────────────


class TestSimulatePath:
    def test_path_shape(self, chain: MarkovChain) -> None:
        path = simulate_path(chain, initial=0, n_steps=10, seed=42)
        assert path.shape == (11,)

    def test_starts_at_initial(self, chain: MarkovChain) -> None:
        for s in range(chain.n_states):
            path = simulate_path(chain, initial=s, n_steps=5, seed=0)
            assert path[0] == s

    def test_visits_only_valid_states(self, chain: MarkovChain) -> None:
        path = simulate_path(chain, initial=0, n_steps=200, seed=7)
        assert path.min() >= 0
        assert path.max() < chain.n_states

    def test_respects_zero_transitions(self) -> None:
        # Junior can never directly become Senior — that arc is forbidden.
        P = np.array(
            [
                [0.5, 0.5, 0.0],
                [0.0, 0.5, 0.5],
                [0.0, 0.0, 1.0],
            ]
        )
        mc = MarkovChain(P)
        path = simulate_path(mc, initial=0, n_steps=500, seed=1)
        # No 0 → 2 jump in a single step.
        for t in range(len(path) - 1):
            if path[t] == 0:
                assert path[t + 1] in (0, 1)

    def test_seeded_paths_are_reproducible(self, chain: MarkovChain) -> None:
        a = simulate_path(chain, initial=0, n_steps=50, seed=123)
        b = simulate_path(chain, initial=0, n_steps=50, seed=123)
        np.testing.assert_array_equal(a, b)

    def test_different_seeds_differ(self, chain: MarkovChain) -> None:
        a = simulate_path(chain, initial=0, n_steps=200, seed=1)
        b = simulate_path(chain, initial=0, n_steps=200, seed=2)
        # Two long paths from different seeds should disagree somewhere.
        assert not np.array_equal(a, b)

    def test_invalid_initial(self, chain: MarkovChain) -> None:
        with pytest.raises(ValueError, match="initial"):
            simulate_path(chain, initial=99, n_steps=1, seed=0)

    def test_n_steps_zero(self, chain: MarkovChain) -> None:
        path = simulate_path(chain, initial=2, n_steps=0, seed=0)
        assert path.shape == (1,)
        assert path[0] == 2

    def test_negative_n_steps(self, chain: MarkovChain) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            simulate_path(chain, initial=0, n_steps=-1, seed=0)


# simulate_paths ─────────────────────────────────────────────────────────


class TestSimulatePaths:
    def test_paths_shape(self, chain: MarkovChain) -> None:
        paths = simulate_paths(chain, initial=0, n_steps=10, n_paths=50, seed=0)
        assert paths.shape == (50, 11)

    def test_all_paths_start_at_initial(self, chain: MarkovChain) -> None:
        paths = simulate_paths(chain, initial=1, n_steps=5, n_paths=20, seed=0)
        assert (paths[:, 0] == 1).all()

    def test_all_visited_states_valid(self, chain: MarkovChain) -> None:
        paths = simulate_paths(chain, initial=0, n_steps=20, n_paths=100, seed=0)
        assert paths.min() >= 0
        assert paths.max() < chain.n_states

    def test_invalid_n_paths(self, chain: MarkovChain) -> None:
        with pytest.raises(ValueError, match="n_paths"):
            simulate_paths(chain, initial=0, n_steps=1, n_paths=0, seed=0)


# Empirical → analytical convergence ─────────────────────────────────────


class TestEmpiricalConvergence:
    def test_distribution_converges_to_pn_e_initial(
        self, chain: MarkovChain
    ) -> None:
        """Empirical distribution at step n approaches P^n applied to e_initial."""
        n_steps = 5
        n_paths = 5000
        initial = 0
        paths = simulate_paths(
            chain, initial=initial, n_steps=n_steps, n_paths=n_paths, seed=2026
        )
        empirical = state_distribution_at(paths, step=n_steps, n_states=chain.n_states)
        e_initial = np.zeros(chain.n_states)
        e_initial[initial] = 1.0
        analytical = e_initial @ chain.n_step(n_steps)
        # 5000 paths gives Monte Carlo error ~ 1/sqrt(5000) ~ 0.014.
        np.testing.assert_allclose(empirical, analytical, atol=0.03)


class TestStateDistributionAt:
    def test_distribution_sums_to_one(self, chain: MarkovChain) -> None:
        paths = simulate_paths(chain, initial=0, n_steps=3, n_paths=100, seed=0)
        dist = state_distribution_at(paths, step=2, n_states=chain.n_states)
        assert np.isclose(dist.sum(), 1.0)
        assert (dist >= 0).all()
        assert dist.shape == (chain.n_states,)

    def test_step_at_zero_is_concentrated(self, chain: MarkovChain) -> None:
        paths = simulate_paths(chain, initial=2, n_steps=1, n_paths=10, seed=0)
        dist = state_distribution_at(paths, step=0, n_states=chain.n_states)
        expected = np.array([0.0, 0.0, 1.0, 0.0])
        np.testing.assert_allclose(dist, expected)

    def test_invalid_step(self, chain: MarkovChain) -> None:
        paths = simulate_paths(chain, initial=0, n_steps=2, n_paths=5, seed=0)
        with pytest.raises(ValueError, match="step"):
            state_distribution_at(paths, step=99, n_states=chain.n_states)
