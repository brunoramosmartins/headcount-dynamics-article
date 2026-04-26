"""Absorption analysis for discrete-time Markov chains.

Tools for absorbing chains and first-passage analysis on general chains:

- canonical-form decomposition into transient/absorbing blocks
- the fundamental matrix ``N = (I - Q)^{-1}``
- expected absorption time ``t = N 1``
- absorption probabilities ``B = N R``
- variance of absorption time
- first-passage time and hitting probability between any two states

All functions accept either a raw transition matrix ``P`` (numpy array) or a
``MarkovChain`` instance from :mod:`src.dtmc`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .dtmc import MarkovChain

ArrayF = NDArray[np.float64]
ArrayI = NDArray[np.int64]


@dataclass(frozen=True)
class CanonicalForm:
    """Canonical decomposition of an absorbing chain.

    Attributes:
        Q: Transient-to-transient submatrix, shape (r, r).
        R: Transient-to-absorbing submatrix, shape (r, K - r).
        transient: Indices of transient states (in original ordering).
        absorbing: Indices of absorbing states (in original ordering).
    """

    Q: ArrayF
    R: ArrayF
    transient: ArrayI
    absorbing: ArrayI


# ──────────────────────────────────────────────────────────────────────────
# Canonical form
# ──────────────────────────────────────────────────────────────────────────


def _as_matrix(P: ArrayF | MarkovChain) -> tuple[ArrayF, float]:
    """Return (P_array, atol) accepting either form of input."""
    if isinstance(P, MarkovChain):
        return P.P, P.atol
    arr = np.asarray(P, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
        raise ValueError(f"P must be a square 2D array; got shape {arr.shape}.")
    if not np.allclose(arr.sum(axis=1), 1.0, atol=1e-9):
        raise ValueError("Rows of P must sum to 1.")
    return arr, 1e-9


def canonical_form(P: ArrayF | MarkovChain) -> CanonicalForm:
    """Split ``P`` into transient and absorbing blocks.

    A state ``i`` is absorbing if ``P[i, i] == 1``. Every other state is
    treated as transient (the function does not check that transient states
    actually reach an absorbing state — see :func:`is_absorbing_chain`).

    Args:
        P: Square row-stochastic matrix or :class:`MarkovChain`.

    Returns:
        :class:`CanonicalForm` with the four pieces.

    Raises:
        ValueError: if no absorbing state is found.
    """
    arr, atol = _as_matrix(P)
    diag = np.diag(arr)
    absorbing_mask = np.isclose(diag, 1.0, atol=atol)
    absorbing = np.where(absorbing_mask)[0].astype(np.int64)
    transient = np.where(~absorbing_mask)[0].astype(np.int64)
    if absorbing.size == 0:
        raise ValueError("P has no absorbing state; absorption analysis does not apply.")
    if transient.size == 0:
        # Degenerate but valid: every state is absorbing.
        return CanonicalForm(
            Q=np.empty((0, 0)),
            R=np.empty((0, absorbing.size)),
            transient=transient,
            absorbing=absorbing,
        )
    Q = arr[np.ix_(transient, transient)]
    R = arr[np.ix_(transient, absorbing)]
    return CanonicalForm(Q=Q, R=R, transient=transient, absorbing=absorbing)


def is_absorbing_chain(P: ArrayF | MarkovChain) -> bool:
    """True iff every transient state can reach an absorbing state.

    Returns ``False`` when ``P`` has no absorbing state at all.
    """
    try:
        cf = canonical_form(P)
    except ValueError:
        return False
    if cf.transient.size == 0:
        return True
    # Reachability within the transient block plus a step to absorbing.
    arr, atol = _as_matrix(P)
    R_mat = arr > atol
    # Warshall transitive closure on full state space, then check that from
    # every transient state there's a path ending in some absorbing state.
    K = arr.shape[0]
    np.fill_diagonal(R_mat, True)
    for k in range(K):
        R_mat = R_mat | (R_mat[:, [k]] & R_mat[[k], :])
    for i in cf.transient:
        if not any(R_mat[i, a] for a in cf.absorbing):
            return False
    return True


# ──────────────────────────────────────────────────────────────────────────
# Core quantities
# ──────────────────────────────────────────────────────────────────────────


def fundamental_matrix(P: ArrayF | MarkovChain) -> ArrayF:
    """Return ``N = (I - Q)^{-1}`` for an absorbing chain.

    ``N[i, j]`` is the expected number of visits to transient state ``j``
    starting from transient state ``i`` before absorption (transient indexing
    follows :func:`canonical_form`).
    """
    cf = canonical_form(P)
    if cf.transient.size == 0:
        return np.empty((0, 0))
    r = cf.Q.shape[0]
    return np.linalg.solve(np.eye(r) - cf.Q, np.eye(r))


def absorption_times(P: ArrayF | MarkovChain) -> ArrayF:
    """Expected time to absorption from each transient state: ``t = N 1``."""
    N = fundamental_matrix(P)
    if N.size == 0:
        return np.empty(0)
    return N.sum(axis=1)


def absorption_variances(P: ArrayF | MarkovChain) -> ArrayF:
    """Variance of absorption time per transient state.

    Computed as ``(2N - I) t - t^2`` element-wise. See Grinstead & Snell §11.2.
    """
    N = fundamental_matrix(P)
    if N.size == 0:
        return np.empty(0)
    t = N.sum(axis=1)
    r = N.shape[0]
    return (2 * N - np.eye(r)) @ t - t * t


def absorption_probabilities(P: ArrayF | MarkovChain) -> ArrayF:
    """Absorption-probability matrix ``B = N R``.

    ``B[i, a]`` is the probability of being absorbed at the ``a``-th absorbing
    state given start at the ``i``-th transient state.
    """
    cf = canonical_form(P)
    if cf.transient.size == 0:
        # Already absorbed — probability 1 at the starting state, 0 elsewhere.
        return np.eye(cf.absorbing.size)
    N = fundamental_matrix(P)
    return N @ cf.R


# ──────────────────────────────────────────────────────────────────────────
# First passage / hitting on general chains
# ──────────────────────────────────────────────────────────────────────────


def _make_target_absorbing(P: ArrayF, target: int) -> ArrayF:
    """Return a copy of P with row ``target`` replaced by the indicator of ``target``."""
    Q = P.copy()
    Q[target, :] = 0.0
    Q[target, target] = 1.0
    return Q


def first_passage_time(
    P: ArrayF | MarkovChain, source: int, target: int
) -> float:
    """Expected number of steps for the chain to first reach ``target``
    starting from ``source``.

    Implementation: turn ``target`` into the only absorbing state and read off
    the expected absorption time from the corresponding entry of ``t``.

    Args:
        P: Transition matrix (or :class:`MarkovChain`).
        source: Starting state index.
        target: State whose first passage time is wanted.

    Returns:
        Expected first passage time. Returns ``0.0`` when ``source == target``
        (interpreting the start as an immediate hit).
    """
    arr, _ = _as_matrix(P)
    K = arr.shape[0]
    if not (0 <= source < K) or not (0 <= target < K):
        raise ValueError("source and target must be valid state indices.")
    if source == target:
        return 0.0
    # Make target absorbing AND remove other absorbing states' self-loops by
    # routing them back into the chain via a uniform restart — but that would
    # change the dynamics. The right model is: keep the other absorbing
    # states as they are; if the chain hits one, the first passage time from
    # source to target is infinite for those paths. To get a well-defined
    # finite expectation, we require source can reach target with probability
    # 1. Otherwise we return +inf.
    modified = _make_target_absorbing(arr, target)
    # If any pre-existing absorbing state (other than target) is reachable
    # from source, the expectation is infinite.
    other_absorbing = [
        i for i in range(K) if i != target and np.isclose(arr[i, i], 1.0)
    ]
    if other_absorbing:
        # Reachability of those states from source.
        reach = arr > 0
        np.fill_diagonal(reach, True)
        for k in range(K):
            reach = reach | (reach[:, [k]] & reach[[k], :])
        if any(reach[source, a] for a in other_absorbing):
            return float("inf")
    cf = canonical_form(modified)
    # Locate source in the transient ordering.
    pos = int(np.where(cf.transient == source)[0][0])
    t = absorption_times(modified)
    return float(t[pos])


def hitting_probability(
    P: ArrayF | MarkovChain, source: int, target: int
) -> float:
    """Probability the chain ever visits ``target`` starting from ``source``.

    Computed by making ``target`` and every other absorbing state remain
    absorbing, then reading the column of ``B`` corresponding to ``target``.
    """
    arr, _ = _as_matrix(P)
    K = arr.shape[0]
    if not (0 <= source < K) or not (0 <= target < K):
        raise ValueError("source and target must be valid state indices.")
    if source == target:
        return 1.0
    modified = _make_target_absorbing(arr, target)
    cf = canonical_form(modified)
    if source in cf.absorbing:
        # Source itself is absorbing in the modified chain: it cannot move,
        # so it never reaches target.
        return 0.0
    B = absorption_probabilities(modified)
    pos_t = int(np.where(cf.transient == source)[0][0])
    pos_a = int(np.where(cf.absorbing == target)[0][0])
    return float(B[pos_t, pos_a])


# ──────────────────────────────────────────────────────────────────────────
# Mean return time from the stationary distribution
# ──────────────────────────────────────────────────────────────────────────


def mean_return_time(mc: MarkovChain, i: int) -> float:
    """Mean return time to state ``i`` for an irreducible chain: ``1 / pi_i``.

    Wrapper around :meth:`MarkovChain.stationary` that returns the recurrence
    interpretation rather than the distribution.
    """
    if not (0 <= i < mc.n_states):
        raise ValueError("i must be a valid state index.")
    pi = mc.stationary()
    if pi[i] <= 0:
        return float("inf")
    return float(1.0 / pi[i])
