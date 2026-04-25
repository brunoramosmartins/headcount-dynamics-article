"""Discrete-time Markov chain (DTMC) module.

Wraps a row-stochastic transition matrix P and provides:

- n-step transitions P^n
- a placeholder for the stationary distribution (full solver in Phase 2)
- structural queries: irreducibility, absorption, communication classes

Conventions:

- States are indexed 0, 1, ..., K-1.
- Distributions are row vectors. They evolve by right-multiplication: pi_{n+1} = pi_n P.
- Each row of P is a conditional distribution and sums to 1.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

ArrayF = NDArray[np.float64]


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
    # Stationary distribution (placeholder; full Phase 2 solver later)
    # ──────────────────────────────────────────────────────────────────

    def stationary(self) -> ArrayF:
        """Return a stationary distribution pi such that pi P = pi.

        Phase 1 placeholder: solves the linear system using the left
        eigenvector of P with eigenvalue 1. For multiple stationary
        distributions (reducible chains) the result is one valid solution
        with no guarantee of uniqueness; Phase 2 develops the full theory.

        Returns:
            Row vector of length K, non-negative, summing to 1.
        """
        # Solve pi (P - I) = 0 with sum(pi) = 1 via least squares.
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
