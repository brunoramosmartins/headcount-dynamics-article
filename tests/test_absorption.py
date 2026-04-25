"""Tests for src/absorption.py — fundamental matrix and hitting times."""

from __future__ import annotations

import numpy as np
import pytest

from src.absorption import (
    absorption_probabilities,
    absorption_times,
    absorption_variances,
    canonical_form,
    first_passage_time,
    fundamental_matrix,
    hitting_probability,
    is_absorbing_chain,
    mean_return_time,
)
from src.dtmc import MarkovChain


# Fixtures ────────────────────────────────────────────────────────────────


@pytest.fixture
def three_state_wpe() -> np.ndarray:
    """Working / Promoted / Exit — the worked example in the notes."""
    return np.array(
        [
            [0.90, 0.05, 0.05],
            [0.00, 1.00, 0.00],
            [0.00, 0.00, 1.00],
        ]
    )


@pytest.fixture
def headcount_absorbing() -> np.ndarray:
    """Full 4-state headcount with Exit absorbing (no replacement)."""
    return np.array(
        [
            [0.93, 0.03, 0.00, 0.04],
            [0.00, 0.96, 0.02, 0.02],
            [0.00, 0.00, 0.99, 0.01],
            [0.00, 0.00, 0.00, 1.00],
        ]
    )


@pytest.fixture
def headcount_recycling() -> MarkovChain:
    """Headcount with replacement hiring (no absorbing state)."""
    P = np.array(
        [
            [0.93, 0.03, 0.00, 0.04],
            [0.00, 0.96, 0.02, 0.02],
            [0.00, 0.00, 0.99, 0.01],
            [0.50, 0.00, 0.00, 0.50],
        ]
    )
    return MarkovChain(P, state_labels=["Junior", "Mid", "Senior", "Exit"])


# Canonical form ──────────────────────────────────────────────────────────


class TestCanonicalForm:
    def test_no_absorbing_raises(self, headcount_recycling: MarkovChain) -> None:
        with pytest.raises(ValueError, match="no absorbing state"):
            canonical_form(headcount_recycling)

    def test_three_state_layout(self, three_state_wpe: np.ndarray) -> None:
        cf = canonical_form(three_state_wpe)
        np.testing.assert_array_equal(cf.transient, [0])
        np.testing.assert_array_equal(cf.absorbing, [1, 2])
        np.testing.assert_allclose(cf.Q, [[0.90]])
        np.testing.assert_allclose(cf.R, [[0.05, 0.05]])

    def test_full_headcount_layout(self, headcount_absorbing: np.ndarray) -> None:
        cf = canonical_form(headcount_absorbing)
        np.testing.assert_array_equal(cf.transient, [0, 1, 2])
        np.testing.assert_array_equal(cf.absorbing, [3])

    def test_accepts_markov_chain(self, three_state_wpe: np.ndarray) -> None:
        mc = MarkovChain(three_state_wpe)
        cf = canonical_form(mc)
        assert cf.transient.size == 1


class TestIsAbsorbingChain:
    def test_recycling_is_not_absorbing(
        self, headcount_recycling: MarkovChain
    ) -> None:
        assert is_absorbing_chain(headcount_recycling) is False

    def test_absorbing_chain(self, headcount_absorbing: np.ndarray) -> None:
        assert is_absorbing_chain(headcount_absorbing) is True


# Fundamental matrix and absorption time ─────────────────────────────────


class TestFundamentalMatrix:
    def test_three_state(self, three_state_wpe: np.ndarray) -> None:
        N = fundamental_matrix(three_state_wpe)
        np.testing.assert_allclose(N, np.array([[10.0]]))

    def test_full_headcount(self, headcount_absorbing: np.ndarray) -> None:
        N = fundamental_matrix(headcount_absorbing)
        expected = np.array(
            [
                [1 / 0.07, 0.03 / (0.07 * 0.04), 0.0006 / (0.07 * 0.04 * 0.01)],
                [0.0, 1 / 0.04, 0.02 / (0.04 * 0.01)],
                [0.0, 0.0, 1 / 0.01],
            ]
        )
        np.testing.assert_allclose(N, expected, atol=1e-9)

    def test_diag_entries_at_least_one(
        self, headcount_absorbing: np.ndarray
    ) -> None:
        N = fundamental_matrix(headcount_absorbing)
        assert (np.diag(N) >= 1.0 - 1e-9).all()


class TestAbsorptionTimes:
    def test_three_state(self, three_state_wpe: np.ndarray) -> None:
        t = absorption_times(three_state_wpe)
        np.testing.assert_allclose(t, [10.0])

    def test_full_headcount(self, headcount_absorbing: np.ndarray) -> None:
        t = absorption_times(headcount_absorbing)
        np.testing.assert_allclose(
            t,
            [46.428571428, 75.0, 100.0],
            atol=1e-6,
        )

    def test_senior_matches_one_over_attrition(
        self, headcount_absorbing: np.ndarray
    ) -> None:
        # Senior has self-loop 0.99 and exit 0.01, so expected time to exit
        # from Senior is 1/0.01 = 100.
        t = absorption_times(headcount_absorbing)
        assert np.isclose(t[2], 1.0 / 0.01)


class TestAbsorptionVariances:
    def test_three_state_geometric(self, three_state_wpe: np.ndarray) -> None:
        # T ~ Geometric(0.1) on {1, 2, ...} so Var = 0.9 / 0.01 = 90.
        v = absorption_variances(three_state_wpe)
        np.testing.assert_allclose(v, [90.0], atol=1e-6)

    def test_non_negative(self, headcount_absorbing: np.ndarray) -> None:
        v = absorption_variances(headcount_absorbing)
        assert (v >= -1e-9).all()


# Absorption probabilities ───────────────────────────────────────────────


class TestAbsorptionProbabilities:
    def test_three_state(self, three_state_wpe: np.ndarray) -> None:
        B = absorption_probabilities(three_state_wpe)
        np.testing.assert_allclose(B, [[0.5, 0.5]], atol=1e-9)

    def test_rows_sum_to_one(self, headcount_absorbing: np.ndarray) -> None:
        B = absorption_probabilities(headcount_absorbing)
        np.testing.assert_allclose(B.sum(axis=1), np.ones(B.shape[0]), atol=1e-9)

    def test_full_headcount_only_exit(
        self, headcount_absorbing: np.ndarray
    ) -> None:
        # Exit is the only absorbing state, so B is the all-ones column.
        B = absorption_probabilities(headcount_absorbing)
        np.testing.assert_allclose(B, np.ones((3, 1)), atol=1e-9)


# First passage and hitting probability ──────────────────────────────────


class TestFirstPassageTime:
    def test_recycling_junior_to_senior_finite(
        self, headcount_recycling: MarkovChain
    ) -> None:
        t = first_passage_time(headcount_recycling, source=0, target=2)
        assert np.isfinite(t)
        # Sanity: pretty large due to small spectral gap, but should be in tens.
        assert 10.0 < t < 1_000.0

    def test_absorbing_variant_junior_to_senior_infinite(
        self, headcount_absorbing: np.ndarray
    ) -> None:
        # In the absorbing variant Exit is pre-existing absorbing and
        # reachable from Junior, so unconditional E[T_S] is infinite.
        t = first_passage_time(headcount_absorbing, source=0, target=2)
        assert np.isinf(t)

    def test_self_returns_zero(self, headcount_recycling: MarkovChain) -> None:
        assert first_passage_time(headcount_recycling, source=2, target=2) == 0.0

    def test_invalid_source(self, three_state_wpe: np.ndarray) -> None:
        with pytest.raises(ValueError, match="state indices"):
            first_passage_time(three_state_wpe, source=99, target=0)


class TestHittingProbability:
    def test_self_is_one(self, headcount_recycling: MarkovChain) -> None:
        assert hitting_probability(headcount_recycling, source=0, target=0) == 1.0

    def test_recycling_always_hits_senior(
        self, headcount_recycling: MarkovChain
    ) -> None:
        # Recycling chain is irreducible — every state visits every other
        # eventually, so hitting prob from Junior to Senior is 1.
        p = hitting_probability(headcount_recycling, source=0, target=2)
        assert np.isclose(p, 1.0, atol=1e-9)

    def test_absorbing_variant_junior_to_senior(
        self, headcount_absorbing: np.ndarray
    ) -> None:
        # With Exit absorbing, only ~21% of Juniors ever reach Senior.
        p = hitting_probability(headcount_absorbing, source=0, target=2)
        assert np.isclose(p, 0.21428571, atol=1e-6)

    def test_absorbing_variant_mid_to_senior(
        self, headcount_absorbing: np.ndarray
    ) -> None:
        # Mid: 50% reach Senior, 50% exit first.
        p = hitting_probability(headcount_absorbing, source=1, target=2)
        assert np.isclose(p, 0.5, atol=1e-6)

    def test_columns_consistent_with_b(
        self, headcount_absorbing: np.ndarray
    ) -> None:
        # Hitting probability to the absorbing state matches B.
        B = absorption_probabilities(headcount_absorbing)
        for source in range(3):
            assert np.isclose(
                hitting_probability(headcount_absorbing, source=source, target=3),
                B[source, 0],
            )


# Mean return time bridges to Phase 2 ─────────────────────────────────────


class TestMeanReturnTime:
    def test_matches_one_over_pi(
        self, headcount_recycling: MarkovChain
    ) -> None:
        pi = headcount_recycling.stationary()
        for i in range(headcount_recycling.n_states):
            assert np.isclose(
                mean_return_time(headcount_recycling, i),
                1.0 / pi[i],
                rtol=1e-9,
            )

    def test_invalid_state(self, headcount_recycling: MarkovChain) -> None:
        with pytest.raises(ValueError, match="state index"):
            mean_return_time(headcount_recycling, 99)
