"""Path simulation for discrete-time Markov chains.

Generates trajectories by sampling next states from the rows of the transition
matrix. All randomness flows through an explicit ``numpy.random.Generator`` so
results are reproducible from a single seed.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from .dtmc import MarkovChain

ArrayI = NDArray[np.int64]
ArrayF = NDArray[np.float64]


def simulate_path(
    mc: MarkovChain,
    initial: int,
    n_steps: int,
    seed: int | None = None,
) -> ArrayI:
    """Simulate a single path of length ``n_steps + 1`` starting at ``initial``.

    Args:
        mc: A MarkovChain.
        initial: Starting state index in [0, K).
        n_steps: Number of transitions to sample. The returned path has
            ``n_steps + 1`` entries (initial state plus ``n_steps`` next states).
        seed: Optional seed for the random generator.

    Returns:
        Integer array of length ``n_steps + 1`` containing the visited states.
    """
    if not (0 <= initial < mc.n_states):
        raise ValueError(f"initial must be in [0, {mc.n_states}); got {initial}.")
    if n_steps < 0:
        raise ValueError("n_steps must be non-negative.")

    rng = np.random.default_rng(seed)
    path = np.empty(n_steps + 1, dtype=np.int64)
    path[0] = initial
    states = np.arange(mc.n_states)
    for t in range(n_steps):
        path[t + 1] = rng.choice(states, p=mc.P[path[t]])
    return path


def simulate_paths(
    mc: MarkovChain,
    initial: int,
    n_steps: int,
    n_paths: int,
    seed: int | None = None,
) -> ArrayI:
    """Simulate ``n_paths`` independent trajectories.

    Args:
        mc: A MarkovChain.
        initial: Starting state index, identical across paths.
        n_steps: Number of transitions per path.
        n_paths: Number of paths.
        seed: Optional seed; uses a single ``Generator`` to sample all paths.

    Returns:
        Integer array of shape ``(n_paths, n_steps + 1)``.
    """
    if n_paths <= 0:
        raise ValueError("n_paths must be positive.")
    if not (0 <= initial < mc.n_states):
        raise ValueError(f"initial must be in [0, {mc.n_states}); got {initial}.")
    if n_steps < 0:
        raise ValueError("n_steps must be non-negative.")

    rng = np.random.default_rng(seed)
    paths = np.empty((n_paths, n_steps + 1), dtype=np.int64)
    paths[:, 0] = initial
    states = np.arange(mc.n_states)

    for t in range(n_steps):
        for k in range(n_paths):
            paths[k, t + 1] = rng.choice(states, p=mc.P[paths[k, t]])
    return paths


def state_distribution_at(paths: ArrayI, step: int, n_states: int) -> ArrayF:
    """Empirical distribution over states at a given time step.

    Args:
        paths: Array of shape ``(n_paths, n_steps + 1)``.
        step: Time index in ``[0, n_steps]``.
        n_states: Total number of states K (paths may not cover all of them).

    Returns:
        Float array of length K with empirical proportions, summing to 1.
    """
    if paths.ndim != 2:
        raise ValueError(f"paths must be 2D; got shape {paths.shape}.")
    if not (0 <= step < paths.shape[1]):
        raise ValueError(f"step must be in [0, {paths.shape[1]}); got {step}.")

    counts = np.bincount(paths[:, step], minlength=n_states).astype(np.float64)
    total = counts.sum()
    if total == 0:
        raise ValueError("No samples found at the requested step.")
    return counts / total
