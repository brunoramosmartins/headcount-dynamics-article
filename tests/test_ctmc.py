"""Tests for src/ctmc.py — CTMC and BirthDeathProcess."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.linalg import expm
from scipy.stats import poisson

from src.ctmc import CTMC, BirthDeathProcess


# Two-state generator for hand-checks ─────────────────────────────────────


@pytest.fixture
def two_state_generator() -> np.ndarray:
    """Two-state CTMC: state 0 -> 1 at rate alpha, 1 -> 0 at rate beta."""
    alpha, beta = 1.0, 2.0
    return np.array([[-alpha, alpha], [beta, -beta]])


@pytest.fixture
def two_state_ctmc(two_state_generator: np.ndarray) -> CTMC:
    return CTMC(two_state_generator, state_labels=["A", "B"])


# Construction & validation ──────────────────────────────────────────────


class TestCTMCConstruction:
    def test_rejects_non_square(self) -> None:
        with pytest.raises(ValueError, match="square"):
            CTMC(np.array([[-1.0, 1.0]]))

    def test_rejects_negative_offdiagonals(self) -> None:
        Q = np.array([[1.0, -1.0], [1.0, -1.0]])
        with pytest.raises(ValueError, match="non-negative"):
            CTMC(Q)

    def test_rejects_non_zero_row_sums(self) -> None:
        Q = np.array([[-1.0, 0.5], [1.0, -1.0]])
        with pytest.raises(ValueError, match="sum to 0"):
            CTMC(Q)

    def test_label_length_mismatch(self) -> None:
        Q = -np.eye(3) + np.eye(3, k=1)
        Q[-1, 0] = 1.0
        Q[-1, -1] = -1.0
        with pytest.raises(ValueError, match="length"):
            CTMC(Q, state_labels=["a"])

    def test_n_states(self, two_state_ctmc: CTMC) -> None:
        assert two_state_ctmc.n_states == 2


# Transition matrix ──────────────────────────────────────────────────────


class TestTransitionMatrix:
    def test_at_t_zero_is_identity(self, two_state_ctmc: CTMC) -> None:
        np.testing.assert_allclose(
            two_state_ctmc.transition_matrix(0), np.eye(2)
        )

    def test_matches_scipy_expm(
        self, two_state_ctmc: CTMC, two_state_generator: np.ndarray
    ) -> None:
        for t in [0.5, 1.0, 5.0, 10.0]:
            np.testing.assert_allclose(
                two_state_ctmc.transition_matrix(t),
                expm(two_state_generator * t),
                atol=1e-12,
            )

    def test_rows_stochastic_at_all_times(self, two_state_ctmc: CTMC) -> None:
        for t in [0.1, 1.0, 100.0]:
            P = two_state_ctmc.transition_matrix(t)
            np.testing.assert_allclose(P.sum(axis=1), np.ones(P.shape[0]), atol=1e-10)
            assert (P >= -1e-12).all()

    def test_rejects_negative_t(self, two_state_ctmc: CTMC) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            two_state_ctmc.transition_matrix(-1.0)

    def test_chapman_kolmogorov_continuous(self, two_state_ctmc: CTMC) -> None:
        """P(s + t) = P(s) P(t)."""
        for s, t in [(0.5, 0.5), (1.0, 2.0), (3.0, 4.0)]:
            lhs = two_state_ctmc.transition_matrix(s + t)
            rhs = two_state_ctmc.transition_matrix(s) @ two_state_ctmc.transition_matrix(t)
            np.testing.assert_allclose(lhs, rhs, atol=1e-10)


# Stationary ─────────────────────────────────────────────────────────────


class TestStationary:
    def test_two_state_closed_form(
        self, two_state_ctmc: CTMC, two_state_generator: np.ndarray
    ) -> None:
        """For Q = [[-a,a],[b,-b]], stationary is (b/(a+b), a/(a+b))."""
        a = -two_state_generator[0, 0]
        b = -two_state_generator[1, 1]
        expected = np.array([b / (a + b), a / (a + b)])
        np.testing.assert_allclose(two_state_ctmc.stationary(), expected, atol=1e-10)

    def test_pi_q_equals_zero(self, two_state_ctmc: CTMC) -> None:
        pi = two_state_ctmc.stationary()
        np.testing.assert_allclose(pi @ two_state_ctmc.Q, np.zeros(2), atol=1e-10)


# Birth-death process ────────────────────────────────────────────────────


class TestBirthDeathConstruction:
    def test_invalid_lengths(self) -> None:
        with pytest.raises(ValueError, match="same length"):
            BirthDeathProcess(lambdas=[1.0, 1.0], mus=[0.0])

    def test_negative_rate(self) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            BirthDeathProcess(lambdas=[1.0, -1.0], mus=[0.0, 1.0])

    def test_tridiagonal_structure(self) -> None:
        bd = BirthDeathProcess(
            lambdas=[2.0, 2.0, 0.0],
            mus=[0.0, 1.0, 2.0],
        )
        # Q[i, j] = 0 for |i - j| > 1
        Q = bd.Q
        for i in range(Q.shape[0]):
            for j in range(Q.shape[1]):
                if abs(i - j) > 1:
                    assert Q[i, j] == 0.0

    def test_constant_hiring_helper(self) -> None:
        bd = BirthDeathProcess.from_constant_hiring(
            lambda_rate=2.0, mu_rate=0.025, n_max=10
        )
        # lambda_n = 2 except at the boundary where it is 0.
        assert np.allclose(bd.lambdas[:-1], 2.0)
        assert bd.lambdas[-1] == 0.0
        # mu_n = n * 0.025
        np.testing.assert_allclose(
            bd.mus, np.arange(11) * 0.025, atol=1e-12
        )


class TestBirthDeathStationary:
    def test_poisson_steady_state(self) -> None:
        """M/M/inf model: stationary equals Poisson(rho)."""
        lam, mu = 2.0, 0.025
        n_max = 200
        bd = BirthDeathProcess.from_constant_hiring(
            lambda_rate=lam, mu_rate=mu, n_max=n_max
        )
        pi = bd.stationary_birth_death()
        rho = lam / mu
        expected = poisson.pmf(np.arange(n_max + 1), mu=rho)
        np.testing.assert_allclose(pi, expected, atol=1e-10)

    def test_mean_matches_rho(self) -> None:
        bd = BirthDeathProcess.from_constant_hiring(2.0, 0.025, n_max=200)
        pi = bd.stationary_birth_death()
        mean = np.sum(np.arange(201) * pi)
        assert np.isclose(mean, 80.0, atol=1e-3)

    def test_variance_matches_rho(self) -> None:
        bd = BirthDeathProcess.from_constant_hiring(2.0, 0.025, n_max=200)
        pi = bd.stationary_birth_death()
        mean = np.sum(np.arange(201) * pi)
        var = np.sum((np.arange(201) - mean) ** 2 * pi)
        assert np.isclose(var, 80.0, atol=1e-3)

    def test_solver_methods_agree(self) -> None:
        """Detailed-balance solver agrees with the generic linear solver
        for the M/M/inf chain."""
        bd = BirthDeathProcess.from_constant_hiring(2.0, 0.025, n_max=200)
        pi_db = bd.stationary_birth_death()
        pi_lin = bd.stationary()
        np.testing.assert_allclose(pi_db, pi_lin, atol=1e-8)

    def test_state_dependent_rates(self) -> None:
        """Sanity check on a non-Poisson birth-death (logistic-like).

        Rates: lambda_n = lambda * (1 - n/N), mu_n = n * mu.
        Detailed balance gives a closed form that we verify by inserting
        the stationary back into pi Q = 0.
        """
        N, lam, mu = 5, 1.0, 0.5
        lambdas = np.array([lam * (1 - n / N) for n in range(N + 1)])
        mus = np.array([n * mu for n in range(N + 1)])
        bd = BirthDeathProcess(lambdas=lambdas, mus=mus)
        pi = bd.stationary_birth_death()
        np.testing.assert_allclose(
            pi @ bd.Q, np.zeros(N + 1), atol=1e-9
        )


class TestExpectedTrajectory:
    def test_starts_at_n0(self) -> None:
        bd = BirthDeathProcess.from_constant_hiring(2.0, 0.025, n_max=200)
        traj = bd.expected_trajectory(n0=45, times=[0])
        assert np.isclose(traj[0], 45.0)

    def test_long_run_is_rho(self) -> None:
        bd = BirthDeathProcess.from_constant_hiring(2.0, 0.025, n_max=200)
        traj = bd.expected_trajectory(n0=45, times=[10_000])
        assert np.isclose(traj[0], 80.0, atol=1e-6)

    def test_closed_form(self) -> None:
        """Match the analytical formula rho + (n0 - rho) * exp(-mu t)."""
        lam, mu = 2.0, 0.025
        bd = BirthDeathProcess.from_constant_hiring(lam, mu, n_max=200)
        rho = lam / mu
        n0 = 45.0
        for t in [1.0, 6.0, 12.0, 50.0]:
            expected = rho + (n0 - rho) * np.exp(-mu * t)
            got = bd.expected_trajectory(n0=n0, times=[t])[0]
            assert np.isclose(got, expected, atol=1e-12)


class TestTransientHeadcount:
    def test_p_x12_at_least_50(self) -> None:
        """The motivating headcount question: P(X_12 >= 50 | X_0 = 45)."""
        bd = BirthDeathProcess.from_constant_hiring(2.0, 0.025, n_max=200)
        P12 = bd.transition_matrix(12.0)
        prob = P12[45, 50:].sum()
        # Notes claim ~0.80; allow a bit of slack for the truncation choice.
        assert 0.78 < prob < 0.82

    def test_initial_distribution_is_concentrated(self) -> None:
        bd = BirthDeathProcess.from_constant_hiring(2.0, 0.025, n_max=50)
        e0 = np.zeros(51)
        e0[45] = 1.0
        np.testing.assert_allclose(bd.evolve(e0, t=0.0), e0)


class TestUniformization:
    def test_consistent_with_dtmc_step(self) -> None:
        """exp(Q t) should match the embedded DTMC view for an integer
        number of expected jumps in some limit. This is a soft check that
        the matrix exponential is consistent with the Q passed in."""
        bd = BirthDeathProcess(
            lambdas=[1.0, 1.0, 0.0],
            mus=[0.0, 0.5, 1.0],
        )
        P_small = bd.transition_matrix(1e-4)
        np.testing.assert_allclose(
            P_small, np.eye(3) + bd.Q * 1e-4, atol=1e-7
        )
