"""Discrete-time Markov chain (DTMC) module.

Wraps a row-stochastic transition matrix P and provides:

- n-step transitions P^n
- structural queries: irreducibility, absorption, communication classes,
  periodicity
- spectral analysis: stationary distribution, eigenvalues, spectral gap,
  mixing-time bound

Conventions:

- States are indexed 0, 1, ..., K-1.
- Distributions are row vectors. They evolve by right-multiplication: pi_{n+1} = pi_n P.
- Each row of P is a conditional distribution and sums to 1.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from math import gcd

import numpy as np
from numpy.typing import NDArray

ArrayF = NDArray[np.float64]
ArrayC = NDArray[np.complex128]


@dataclass
class MarkovChain:
    """A finite-state, time-homogeneous discrete-time Markov chain.

    Attributes:
        P: Transition matrix of shape (K, K). Must be row-stochastic.
        state_labels: Optional human-readable labels, length K.
        atol: Numerical tolerance used in row-sum and equality checks.
    """

    P: ArrayF
    state_labels: list[str] = field(default_factory=list)
    atol: float = 1e-9

    def __post_init__(self) -> None:
        P = np.asarray(self.P, dtype=np.float64)
        if P.ndim != 2 or P.shape[0] != P.shape[1]:
            raise ValueError(f"P must be a square 2D array; got shape {P.shape}.")
        if np.any(P < -self.atol):
            raise ValueError("P contains negative entries.")
        row_sums = P.sum(axis=1)
        if not np.allclose(row_sums, 1.0, atol=self.atol):
            raise ValueError(
                f"Rows of P must sum to 1 (within atol={self.atol}); got {row_sums}."
            )
        # Clip tiny negative noise to zero for safety.
        self.P = np.clip(P, 0.0, None)

        K = self.P.shape[0]
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
        """Number of states K."""
        return int(self.P.shape[0])

    def n_step(self, n: int) -> ArrayF:
        """Return the n-step transition matrix P^n.

        Args:
            n: Non-negative integer.

        Returns:
            Array of shape (K, K) equal to P^n. P^0 is the identity.
        """
        if n < 0:
            raise ValueError("n must be non-negative.")
        return np.linalg.matrix_power(self.P, n)

    def evolve(self, pi0: Sequence[float], n: int) -> ArrayF:
        """Evolve a distribution forward by n steps: pi_n = pi_0 P^n.

        Args:
            pi0: Initial distribution, length K, non-negative, sums to 1.
            n: Non-negative integer.

        Returns:
            Row vector pi_n.
        """
        pi = np.asarray(pi0, dtype=np.float64).reshape(-1)
        if pi.size != self.n_states:
            raise ValueError(f"pi0 must have length {self.n_states}.")
        if not np.isclose(pi.sum(), 1.0, atol=self.atol) or np.any(pi < -self.atol):
            raise ValueError("pi0 must be a probability vector.")
        return pi @ self.n_step(n)

    # ──────────────────────────────────────────────────────────────────
    # Spectral analysis
    # ──────────────────────────────────────────────────────────────────

    def eigenvalues(self) -> ArrayC:
        """Return the eigenvalues of P, sorted by descending modulus.

        For a row-stochastic matrix the spectrum lies in the closed unit
        disk and contains 1. Eigenvalues are returned as complex numbers;
        callers can take ``.real`` when periodicity is not a concern.
        """
        eigs = np.linalg.eigvals(self.P)
        order = np.argsort(-np.abs(eigs), kind="stable")
        return eigs[order]

    def stationary(self, method: str = "linear") -> ArrayF:
        """Return a stationary distribution pi such that pi P = pi.

        Args:
            method: One of ``"linear"`` (default — solve the augmented linear
                system pi (P - I) = 0 with sum(pi) = 1) or ``"eigen"`` (take
                the left eigenvector of P with eigenvalue closest to 1).

        Returns:
            Row vector of length K, non-negative, summing to 1.

        Notes:
            For an irreducible aperiodic chain the stationary distribution
            is unique and either method recovers it. For reducible chains
            the linear solver returns one stationary distribution but
            uniqueness is not guaranteed.
        """
        if method == "linear":
            return self._stationary_linear()
        if method == "eigen":
            return self._stationary_eigen()
        raise ValueError(f"Unknown method: {method!r}. Use 'linear' or 'eigen'.")

    def _stationary_linear(self) -> ArrayF:
        K = self.n_states
        A = np.vstack([(self.P.T - np.eye(K)), np.ones(K)])
        b = np.zeros(K + 1)
        b[-1] = 1.0
        pi, *_ = np.linalg.lstsq(A, b, rcond=None)
        pi = np.clip(pi, 0.0, None)
        s = pi.sum()
        if s == 0.0:
            raise RuntimeError("Stationary solver produced the zero vector.")
        return pi / s

    def _stationary_eigen(self) -> ArrayF:
        # Left eigenvectors of P are right eigenvectors of P.T.
        eigvals, eigvecs = np.linalg.eig(self.P.T)
        idx = int(np.argmin(np.abs(eigvals - 1.0)))
        v = eigvecs[:, idx]
        # The eigenvector should be (numerically) real for an irreducible
        # aperiodic chain; drop a tiny imaginary part.
        v = np.real_if_close(v, tol=1e6)
        v = np.asarray(v, dtype=np.float64)
        # Sign and scale so the result is a probability vector.
        if v.sum() < 0:
            v = -v
        v = np.clip(v, 0.0, None)
        s = v.sum()
        if s == 0.0:
            raise RuntimeError("Eigen stationary solver produced the zero vector.")
        return v / s

    def spectral_gap(self) -> float:
        """Return the spectral gap 1 - |lambda_2|.

        ``lambda_2`` is the eigenvalue with the second-largest modulus. For an
        irreducible aperiodic finite chain this is in (0, 1] and controls
        the geometric rate of convergence to the stationary distribution.
        """
        eigs = self.eigenvalues()
        if eigs.size < 2:
            return 1.0
        return float(1.0 - np.abs(eigs[1]))

    def mixing_time(self, eps: float = 0.25) -> int:
        """Smallest n such that ||P^n[i, :] - pi||_TV <= eps for all i.

        Uses total variation distance, defined as half the L1 distance between
        rows of P^n and the stationary distribution. The bound is computed
        empirically by powering P; for chains with very small spectral gap
        the result may exceed ``max_iter`` (in which case the upper bound is
        returned).

        Args:
            eps: Target total variation tolerance (default 1/4, the standard
                convention from Levin–Peres–Wilmer).

        Returns:
            Integer n.
        """
        if not (0.0 < eps < 1.0):
            raise ValueError("eps must lie in (0, 1).")
        pi = self.stationary()
        K = self.n_states
        max_iter = 10_000
        Pn = np.eye(K)
        for n in range(1, max_iter + 1):
            Pn = Pn @ self.P
            tv = 0.5 * np.abs(Pn - pi).sum(axis=1).max()
            if tv <= eps:
                return n
        return max_iter

    # ──────────────────────────────────────────────────────────────────
    # Structural properties
    # ──────────────────────────────────────────────────────────────────

    def _reachability(self) -> NDArray[np.bool_]:
        """Boolean matrix R where R[i, j] = True iff j is accessible from i."""
        K = self.n_states
        R = self.P > 0.0
        # Diagonal: i is accessible from i in 0 steps.
        np.fill_diagonal(R, True)
        # Transitive closure (Warshall).
        for k in range(K):
            R = R | (R[:, [k]] & R[[k], :])
        return R

    def is_irreducible(self) -> bool:
        """True iff every state is accessible from every other state."""
        return bool(self._reachability().all())

    def is_absorbing_state(self, i: int) -> bool:
        """True iff state i is absorbing (p_ii == 1)."""
        return bool(np.isclose(self.P[i, i], 1.0, atol=self.atol))

    def is_absorbing(self) -> bool:
        """True iff the chain has at least one absorbing state and every
        state can eventually reach some absorbing state.
        """
        absorbing = [i for i in range(self.n_states) if self.is_absorbing_state(i)]
        if not absorbing:
            return False
        R = self._reachability()
        for i in range(self.n_states):
            if not any(R[i, a] for a in absorbing):
                return False
        return True

    def communication_classes(self) -> list[set[int]]:
        """Return the equivalence classes of the 'communicates with' relation.

        Two states i, j are in the same class iff i is accessible from j and
        j is accessible from i.
        """
        R = self._reachability()
        # i communicates with j iff R[i,j] and R[j,i].
        comm = R & R.T
        K = self.n_states
        seen: set[int] = set()
        classes: list[set[int]] = []
        for i in range(K):
            if i in seen:
                continue
            cls = {j for j in range(K) if comm[i, j]}
            classes.append(cls)
            seen |= cls
        return classes

    # ──────────────────────────────────────────────────────────────────
    # Periodicity
    # ──────────────────────────────────────────────────────────────────

    def period(self, i: int, max_n: int | None = None) -> int:
        """Return the period d(i) = gcd{n >= 1 : p_ii^(n) > 0}.

        Args:
            i: State index.
            max_n: Largest power of P to inspect. Defaults to ``2 * K``,
                which is sufficient for finite chains because cycle lengths
                are bounded by the number of states in any communication
                class.

        Returns:
            Positive integer (the period). Returns 0 if no return path
            exists (transient state with no self-recurrence).
        """
        K = self.n_states
        if max_n is None:
            max_n = max(2 * K, 4)
        Pn = np.eye(K)
        d = 0
        for n in range(1, max_n + 1):
            Pn = Pn @ self.P
            if Pn[i, i] > self.atol:
                d = gcd(d, n)
                if d == 1:
                    return 1
        return d

    def is_aperiodic(self) -> bool:
        """True iff every state has period 1."""
        return all(self.period(i) == 1 for i in range(self.n_states))

    # ──────────────────────────────────────────────────────────────────
    # Convenience
    # ──────────────────────────────────────────────────────────────────

    def label(self, i: int) -> str:
        """Return the human-readable label of state i."""
        return self.state_labels[i]

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (
            f"MarkovChain(n_states={self.n_states}, "
            f"labels={self.state_labels})"
        )
