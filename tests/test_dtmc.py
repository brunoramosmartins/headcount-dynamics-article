"""Tests for src/dtmc.py — MarkovChain class."""

from __future__ import annotations

import numpy as np
import pytest

from src.dtmc import MarkovChain

# Shared fixtures ──────────────────────────────────────────────────────────


@pytest.fixture
def headcount_recycling() -> MarkovChain:
    """Headcount model with replacement hiring (irreducible)."""
    P = np.array(
        [
            [0.93, 0.03, 0.00, 0.04],
            [0.00, 0.96, 0.02, 0.02],
            [0.00, 0.00, 0.99, 0.01],
            [0.50, 0.00, 0.00, 0.50],
        ]
    )
    return MarkovChain(P, state_labels=["Junior", "Mid", "Senior", "Exit"])


@pytest.fixture
def headcount_absorbing() -> MarkovChain:
    """Headcount model with Exit absorbing (no replacement hiring)."""
    P = np.array(
        [
            [0.93, 0.03, 0.00, 0.04],
            [0.00, 0.96, 0.02, 0.02],
            [0.00, 0.00, 0.99, 0.01],
            [0.00, 0.00, 0.00, 1.00],
        ]
    )
    return MarkovChain(P, state_labels=["Junior", "Mid", "Senior", "Exit"])


@pytest.fixture
def three_state_upper() -> MarkovChain:
    """3-state upper-triangular chain from worked example in notes."""
    P = np.array(
        [
            [0.7, 0.3, 0.0],
            [0.0, 0.8, 0.2],
            [0.0, 0.0, 1.0],
        ]
    )
    return MarkovChain(P)


# Construction & validation ────────────────────────────────────────────────


class TestConstruction:
    def test_rejects_non_square(self) -> None:
        with pytest.raises(ValueError, match="square"):
            MarkovChain(np.array([[0.5, 0.5]]))

    def test_rejects_negative_entries(self) -> None:
        P = np.array([[0.5, 0.6, -0.1], [1, 0, 0], [0, 0, 1]])
        with pytest.raises(ValueError, match="negative"):
            MarkovChain(P)

    def test_rejects_non_stochastic_rows(self) -> None:
        P = np.array([[0.5, 0.4], [0.0, 1.0]])
        with pytest.raises(ValueError, match="sum to 1"):
            MarkovChain(P)

    def test_default_labels(self, three_state_upper: MarkovChain) -> None:
        assert three_state_upper.state_labels == ["0", "1", "2"]

    def test_label_length_mismatch(self) -> None:
        P = np.eye(3)
        with pytest.raises(ValueError, match="length"):
            MarkovChain(P, state_labels=["a", "b"])

    def test_n_states(self, headcount_recycling: MarkovChain) -> None:
        assert headcount_recycling.n_states == 4


# Theorem: P^1 = P, P^n is row-stochastic, distribution evolution ─────────


class TestNStep:
    def test_p1_equals_p(self, three_state_upper: MarkovChain) -> None:
        np.testing.assert_allclose(three_state_upper.n_step(1), three_state_upper.P)

    def test_p0_is_identity(self, three_state_upper: MarkovChain) -> None:
        np.testing.assert_allclose(
            three_state_upper.n_step(0), np.eye(three_state_upper.n_states)
        )

    @pytest.mark.parametrize("n", [1, 2, 3, 5, 10, 25, 100])
    def test_pn_rows_sum_to_one(
        self, headcount_recycling: MarkovChain, n: int
    ) -> None:
        Pn = headcount_recycling.n_step(n)
        np.testing.assert_allclose(Pn.sum(axis=1), np.ones(Pn.shape[0]), atol=1e-10)

    @pytest.mark.parametrize("n", [1, 2, 5, 10])
    def test_pn_non_negative(
        self, headcount_recycling: MarkovChain, n: int
    ) -> None:
        Pn = headcount_recycling.n_step(n)
        assert (Pn >= -1e-12).all()

    def test_negative_n_rejected(self, three_state_upper: MarkovChain) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            three_state_upper.n_step(-1)

    def test_p2_matches_worked_example(self, three_state_upper: MarkovChain) -> None:
        expected = np.array(
            [
                [0.49, 0.45, 0.06],
                [0.00, 0.64, 0.36],
                [0.00, 0.00, 1.00],
            ]
        )
        np.testing.assert_allclose(three_state_upper.n_step(2), expected, atol=1e-12)

    def test_p3_matches_worked_example(self, three_state_upper: MarkovChain) -> None:
        expected = np.array(
            [
                [0.343, 0.507, 0.150],
                [0.000, 0.512, 0.488],
                [0.000, 0.000, 1.000],
            ]
        )
        np.testing.assert_allclose(three_state_upper.n_step(3), expected, atol=1e-12)

    def test_chapman_kolmogorov(self, headcount_recycling: MarkovChain) -> None:
        """P^{m+n} = P^m P^n for several (m, n)."""
        for m, n in [(1, 1), (2, 3), (5, 7), (10, 4)]:
            Pmn = headcount_recycling.n_step(m + n)
            Pm_Pn = headcount_recycling.n_step(m) @ headcount_recycling.n_step(n)
            np.testing.assert_allclose(Pmn, Pm_Pn, atol=1e-10)


class TestEvolve:
    def test_evolve_one_step(self, three_state_upper: MarkovChain) -> None:
        pi0 = np.array([1.0, 0.0, 0.0])
        pi1 = three_state_upper.evolve(pi0, 1)
        np.testing.assert_allclose(pi1, np.array([0.7, 0.3, 0.0]))

    def test_evolve_matches_pn(self, headcount_recycling: MarkovChain) -> None:
        pi0 = np.array([1.0, 0.0, 0.0, 0.0])
        for n in [1, 5, 12]:
            pi_n_evolve = headcount_recycling.evolve(pi0, n)
            pi_n_matrix = pi0 @ headcount_recycling.n_step(n)
            np.testing.assert_allclose(pi_n_evolve, pi_n_matrix, atol=1e-12)

    def test_invalid_pi0(self, three_state_upper: MarkovChain) -> None:
        with pytest.raises(ValueError, match="probability"):
            three_state_upper.evolve(np.array([0.5, 0.4, 0.0]), 1)


# Stationary distribution (placeholder) ───────────────────────────────────


class TestStationaryPlaceholder:
    def test_stationary_satisfies_pi_p_equals_pi(
        self, headcount_recycling: MarkovChain
    ) -> None:
        pi = headcount_recycling.stationary()
        np.testing.assert_allclose(pi @ headcount_recycling.P, pi, atol=1e-8)

    def test_stationary_is_a_distribution(
        self, headcount_recycling: MarkovChain
    ) -> None:
        pi = headcount_recycling.stationary()
        assert (pi >= 0).all()
        assert np.isclose(pi.sum(), 1.0)


# Structural properties ──────────────────────────────────────────────────


class TestStructure:
    def test_recycling_is_irreducible(
        self, headcount_recycling: MarkovChain
    ) -> None:
        assert headcount_recycling.is_irreducible() is True

    def test_absorbing_is_not_irreducible(
        self, headcount_recycling: MarkovChain, headcount_absorbing: MarkovChain
    ) -> None:
        assert headcount_absorbing.is_irreducible() is False

    def test_recycling_has_one_communication_class(
        self, headcount_recycling: MarkovChain
    ) -> None:
        classes = headcount_recycling.communication_classes()
        assert len(classes) == 1
        assert classes[0] == {0, 1, 2, 3}

    def test_absorbing_communication_classes(
        self, headcount_absorbing: MarkovChain
    ) -> None:
        classes = headcount_absorbing.communication_classes()
        # Each transient state and the absorbing state form their own class.
        as_frozen = {frozenset(c) for c in classes}
        assert as_frozen == {
            frozenset({0}),
            frozenset({1}),
            frozenset({2}),
            frozenset({3}),
        }

    def test_is_absorbing_state(
        self, headcount_absorbing: MarkovChain, headcount_recycling: MarkovChain
    ) -> None:
        assert headcount_absorbing.is_absorbing_state(3) is True
        for i in range(3):
            assert headcount_absorbing.is_absorbing_state(i) is False
        for i in range(headcount_recycling.n_states):
            assert headcount_recycling.is_absorbing_state(i) is False

    def test_is_absorbing_chain(
        self, headcount_absorbing: MarkovChain, headcount_recycling: MarkovChain
    ) -> None:
        assert headcount_absorbing.is_absorbing() is True
        assert headcount_recycling.is_absorbing() is False

    def test_label_lookup(self, headcount_recycling: MarkovChain) -> None:
        assert headcount_recycling.label(2) == "Senior"


# Phase 2 — Spectral analysis ─────────────────────────────────────────────


class TestEigenvalues:
    def test_one_is_eigenvalue(self, headcount_recycling: MarkovChain) -> None:
        eigs = headcount_recycling.eigenvalues()
        assert np.isclose(np.abs(eigs[0]), 1.0)

    def test_all_eigenvalues_in_unit_disk(
        self, headcount_recycling: MarkovChain
    ) -> None:
        eigs = headcount_recycling.eigenvalues()
        assert np.all(np.abs(eigs) <= 1.0 + 1e-10)

    def test_eigenvalues_sorted_by_modulus(
        self, headcount_recycling: MarkovChain
    ) -> None:
        eigs = headcount_recycling.eigenvalues()
        moduli = np.abs(eigs)
        assert np.all(moduli[:-1] >= moduli[1:] - 1e-12)


class TestStationary:
    def test_linear_solver_satisfies_pi_p(
        self, headcount_recycling: MarkovChain
    ) -> None:
        pi = headcount_recycling.stationary(method="linear")
        np.testing.assert_allclose(
            pi @ headcount_recycling.P, pi, atol=1e-9
        )

    def test_eigen_solver_satisfies_pi_p(
        self, headcount_recycling: MarkovChain
    ) -> None:
        pi = headcount_recycling.stationary(method="eigen")
        np.testing.assert_allclose(
            pi @ headcount_recycling.P, pi, atol=1e-9
        )

    def test_methods_agree(self, headcount_recycling: MarkovChain) -> None:
        a = headcount_recycling.stationary(method="linear")
        b = headcount_recycling.stationary(method="eigen")
        np.testing.assert_allclose(a, b, atol=1e-8)

    def test_unknown_method(self, headcount_recycling: MarkovChain) -> None:
        with pytest.raises(ValueError, match="Unknown method"):
            headcount_recycling.stationary(method="bogus")

    def test_stationary_matches_hand_solution(
        self, headcount_recycling: MarkovChain
    ) -> None:
        """Compare against the hand calculation in notes/phase2-dtmc-analysis.md."""
        pi = headcount_recycling.stationary()
        # Junior fraction is 1 / (1 + 0.75 + 1.5 + 0.14).
        expected_junior = 1.0 / 3.39
        assert np.isclose(pi[0], expected_junior, atol=1e-3)


class TestSpectralGap:
    def test_gap_in_unit_interval(self, headcount_recycling: MarkovChain) -> None:
        gap = headcount_recycling.spectral_gap()
        assert 0.0 <= gap <= 1.0

    def test_iid_chain_has_full_gap(self) -> None:
        """If every row of P equals pi, |lambda_2| = 0 and gap = 1."""
        P = np.tile(np.array([0.2, 0.3, 0.5]), (3, 1))
        mc = MarkovChain(P)
        assert np.isclose(mc.spectral_gap(), 1.0, atol=1e-10)

    def test_high_diagonal_has_small_gap(
        self, headcount_recycling: MarkovChain
    ) -> None:
        # Headcount chain has p_ii close to 1 for transient states, so the
        # second eigenvalue is also close to 1 and the gap is small.
        assert headcount_recycling.spectral_gap() < 0.05


class TestMixingTime:
    def test_iid_chain_mixes_in_one_step(self) -> None:
        P = np.tile(np.array([0.2, 0.3, 0.5]), (3, 1))
        mc = MarkovChain(P)
        assert mc.mixing_time() == 1

    def test_invalid_eps(self, headcount_recycling: MarkovChain) -> None:
        with pytest.raises(ValueError, match="eps"):
            headcount_recycling.mixing_time(eps=1.5)

    def test_smaller_eps_takes_longer(
        self, headcount_recycling: MarkovChain
    ) -> None:
        t_loose = headcount_recycling.mixing_time(eps=0.25)
        t_tight = headcount_recycling.mixing_time(eps=0.01)
        assert t_tight >= t_loose


# Periodicity ─────────────────────────────────────────────────────────────


class TestPeriod:
    def test_self_loop_period_is_one(self, headcount_recycling: MarkovChain) -> None:
        for i in range(headcount_recycling.n_states):
            assert headcount_recycling.period(i) == 1

    def test_period_two_chain(self) -> None:
        # Two-state alternating chain: P = [[0,1],[1,0]] has period 2.
        P = np.array([[0.0, 1.0], [1.0, 0.0]])
        mc = MarkovChain(P)
        assert mc.period(0) == 2
        assert mc.period(1) == 2
        assert mc.is_aperiodic() is False

    def test_aperiodic(self, headcount_recycling: MarkovChain) -> None:
        assert headcount_recycling.is_aperiodic() is True
