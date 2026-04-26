"""Continuous-time Markov chain (CTMC) module.

Wraps a generator matrix Q and provides:

- the transition matrix P(t) = exp(Q t)
- the stationary distribution pi as the left null vector of Q
- a :class:`BirthDeathProcess` subclass with closed-form detailed-balance
  stationary and tridiagonal-generator construction

Conventions:

- States are indexed 0, 1, ..., K-1.
- Distributions are row vectors and evolve by right multiplication:
  pi(t) = pi(0) exp(Q t).
- Off-diagonal entries Q[i, j] (i != j) are non-negative transition rates;
  diagonals Q[i, i] = -sum_{j != i} Q[i, j] so each row sums to 0.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import expm

ArrayF = NDArray[np.float64]


# ──────────────────────────────────────────────────────────────────────────
# Base CTMC
# ──────────────────────────────────────────────────────────────────────────


@dataclass
class CTMC:
    """A finite-state, time-homogeneous continuous-time Markov chain.

    Attributes:
        Q: Generator matrix of shape (K, K). Off-diagonals must be
            non-negative; rows must sum to 0 (within ``atol``).
        state_labels: Optional human-readable labels of length K.
        atol: Numerical tolerance for the row-sum and stationary checks.
    """

    Q: ArrayF
    state_labels: list[str] = field(default_factory=list)
    atol: float = 1e-9

    def __post_init__(self) -> None:
        Q = np.asarray(self.Q, dtype=np.float64)
        if Q.ndim != 2 or Q.shape[0] != Q.shape[1]:
            raise ValueError(f"Q must be a square 2D array; got shape {Q.shape}.")
        K = Q.shape[0]
        # Off-diagonals must be non-negative.
        offdiag = Q.copy()
        np.fill_diagonal(offdiag, 0.0)
        if np.any(offdiag < -self.atol):
            raise ValueError("Off-diagonal entries of Q must be non-negative.")
        # Rows must sum to 0.
        row_sums = Q.sum(axis=1)
        if not np.allclose(row_sums, 0.0, atol=self.atol):
            raise ValueError(
                f"Rows of Q must sum to 0 (within atol={self.atol}); got {row_sums}."
            )
        # Clip any small negative noise on off-diagonals.
        np.fill_diagonal(offdiag, 0.0)
        offdiag = np.clip(offdiag, 0.0, None)
        diag = -offdiag.sum(axis=1)
        self.Q = offdiag + np.diag(diag)

        if not self.state_labels:
            self.state_labels = [str(i) for i in range(K)]
        elif len(self.state_labels) != K:
            raise ValueError(
                f"state_labels has length {len(self.state_labels)}, expected {K}."
            )

    # ──────────────────────────────────────────────────────────────────
    # Basic structure
    # ──────────────────────────────────────────────────────────────────

    @property
    def n_states(self) -> int:
        return int(self.Q.shape[0])

    def transition_matrix(self, t: float) -> ArrayF:
        """Return P(t) = exp(Q t).

        Args:
            t: Non-negative time horizon.

        Returns:
            Array of shape (K, K) which is row-stochastic and equal to the
            CTMC's transition matrix at time t.
        """
        if t < 0:
            raise ValueError("t must be non-negative.")
        return expm(self.Q * t)

    def evolve(self, pi0: Sequence[float], t: float) -> ArrayF:
        """Evolve a distribution forward by time t: pi(t) = pi0 exp(Q t)."""
        pi = np.asarray(pi0, dtype=np.float64).reshape(-1)
        if pi.size != self.n_states:
            raise ValueError(f"pi0 must have length {self.n_states}.")
        if not np.isclose(pi.sum(), 1.0, atol=self.atol) or np.any(pi < -self.atol):
            raise ValueError("pi0 must be a probability vector.")
        return pi @ self.transition_matrix(t)

    # ──────────────────────────────────────────────────────────────────
    # Stationary distribution
    # ──────────────────────────────────────────────────────────────────

    def stationary(self) -> ArrayF:
        """Return a stationary distribution pi solving pi Q = 0, sum(pi) = 1.

        Solves the augmented linear system [(Q^T)^stack(1...1)] pi = e_K
        via least squares. For irreducible chains the solution is unique;
        otherwise one valid stationary is returned.
        """
        K = self.n_states
        A = np.vstack([self.Q.T, np.ones(K)])
        b = np.zeros(K + 1)
        b[-1] = 1.0
        pi, *_ = np.linalg.lstsq(A, b, rcond=None)
        pi = np.clip(pi, 0.0, None)
        s = pi.sum()
        if s == 0.0:
            raise RuntimeError("Stationary solver produced the zero vector.")
        return pi / s

    # ──────────────────────────────────────────────────────────────────
    # Convenience
    # ──────────────────────────────────────────────────────────────────

    def label(self, i: int) -> str:
        return self.state_labels[i]

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"CTMC(n_states={self.n_states}, labels={self.state_labels})"


# ──────────────────────────────────────────────────────────────────────────
# Birth-Death specialisation
# ──────────────────────────────────────────────────────────────────────────


class BirthDeathProcess(CTMC):
    """Birth-death process on states {0, 1, ..., n_max}.

    Constructs the tridiagonal generator from per-state birth rates
    ``lambdas[n]`` (rate n -> n+1) and death rates ``mus[n]`` (rate n -> n-1).
    By convention ``mus[0] = 0`` and ``lambdas[n_max] = 0`` (boundary).

    For state-dependent rates the user supplies the full arrays; for
    M/M/inf-style hiring with constant lambda and per-capita mu use
    :meth:`from_constant_hiring`.
    """

    def __init__(
        self,
        lambdas: Sequence[float],
        mus: Sequence[float],
        state_labels: list[str] | None = None,
    ) -> None:
        lam = np.asarray(lambdas, dtype=np.float64)
        mu = np.asarray(mus, dtype=np.float64)
        if lam.shape != mu.shape:
            raise ValueError("lambdas and mus must have the same length.")
        if lam.ndim != 1:
            raise ValueError("lambdas and mus must be 1D arrays.")
        if np.any(lam < 0) or np.any(mu < 0):
            raise ValueError("Rates must be non-negative.")
        n_max = lam.size - 1  # states 0..n_max
        Q = np.zeros((n_max + 1, n_max + 1))
        for n in range(n_max + 1):
            if n + 1 <= n_max:
                Q[n, n + 1] = lam[n]
            if n - 1 >= 0:
                Q[n, n - 1] = mu[n]
            Q[n, n] = -(Q[n].sum())
        labels = state_labels or [str(n) for n in range(n_max + 1)]
        super().__init__(Q=Q, state_labels=labels)
        self.lambdas = lam
        self.mus = mu

    @classmethod
    def from_constant_hiring(
        cls,
        lambda_rate: float,
        mu_rate: float,
        n_max: int,
    ) -> BirthDeathProcess:
        """Convenience constructor: lambda_n = lambda, mu_n = n * mu.

        Models a team where hires arrive at constant rate and each employee
        leaves independently at rate mu per unit time. The truncation at
        ``n_max`` should be large enough that the stationary distribution
        has negligible mass beyond n_max (e.g. n_max >> lambda / mu + 5 sigma).
        """
        if lambda_rate < 0 or mu_rate < 0:
            raise ValueError("rates must be non-negative.")
        if n_max < 1:
            raise ValueError("n_max must be at least 1.")
        lambdas = np.full(n_max + 1, lambda_rate, dtype=np.float64)
        lambdas[-1] = 0.0  # absorbing boundary at n_max
        mus = np.array([n * mu_rate for n in range(n_max + 1)], dtype=np.float64)
        return cls(lambdas=lambdas, mus=mus)

    # ──────────────────────────────────────────────────────────────────
    # Stationary via detailed balance
    # ──────────────────────────────────────────────────────────────────

    def stationary_birth_death(self) -> ArrayF:
        """Stationary distribution from the detailed-balance recurrence.

        Computes pi_n = pi_0 * prod_{k<n} lambda_k / mu_{k+1} and normalises.
        Numerically more stable than the eigenvalue solver for very small rates.
        """
        n_max = self.n_states - 1
        pi = np.empty(n_max + 1)
        pi[0] = 1.0
        for n in range(n_max):
            mu_next = self.mus[n + 1]
            if mu_next <= 0:
                raise ValueError(
                    f"mu[{n + 1}] is zero — detailed-balance recursion blows up."
                )
            pi[n + 1] = pi[n] * self.lambdas[n] / mu_next
        s = pi.sum()
        if s == 0.0 or not np.isfinite(s):
            raise RuntimeError(
                "Detailed-balance product diverged; increase n_max or check rates."
            )
        return pi / s

    def expected_trajectory(
        self,
        n0: float,
        times: Sequence[float],
        lambda_rate: float | None = None,
        mu_rate: float | None = None,
    ) -> ArrayF:
        """Closed-form expected trajectory under the M/M/inf model.

        For lambda_n = lambda (constant) and mu_n = n * mu, the mean satisfies
        dn/dt = lambda - mu * n with solution
        ``E[n(t)] = rho + (n0 - rho) * exp(-mu * t)`` where rho = lambda/mu.

        Args:
            n0: Initial team size (a real number; expectation has no integer
                constraint).
            times: 1D array of non-negative times.
            lambda_rate: Override hiring rate; defaults to ``self.lambdas[0]``.
            mu_rate: Override per-capita attrition rate; defaults to
                ``self.mus[1]`` (which equals ``mu`` under M/M/inf).

        Returns:
            Array of shape ``(len(times),)`` with the expected team size at
            each time.
        """
        lam = float(self.lambdas[0]) if lambda_rate is None else float(lambda_rate)
        mu = float(self.mus[1]) if mu_rate is None else float(mu_rate)
        if mu <= 0:
            raise ValueError("mu_rate must be positive for the M/M/inf trajectory.")
        rho = lam / mu
        t = np.asarray(times, dtype=np.float64)
        return rho + (n0 - rho) * np.exp(-mu * t)
