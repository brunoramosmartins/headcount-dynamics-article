# Exercise Set 3 — Absorption, Hitting Times, and the Fundamental Matrix

## Proofs (paper)

**1.** Prove that if $Q$ is the transient-to-transient submatrix of an absorbing
chain, then all eigenvalues of $Q$ satisfy $|\lambda| < 1$.
*Hint:* $Q^n \to 0$ as $n \to \infty$. What does this imply about the spectral
radius?

**2.** Prove that $N = (I - Q)^{-1}$ exists and equals $\sum_{k=0}^{\infty} Q^k$.
*Technique:* Neumann series, condition $\rho(Q) < 1$.

**3.** Prove that $N_{ij} = \mathbb{E}[\text{number of visits to } j \text{ before absorption} \mid X_0 = i]$.
*Hint:* $N = \sum_{k=0}^{\infty} Q^k$, and
$(Q^k)_{ij} = P(X_k = j, \text{not absorbed before } k \mid X_0 = i)$.

**4.** Prove the mean return time theorem: for an irreducible positive
recurrent chain, $\mathbb{E}[T_i \mid X_0 = i] = 1/\pi_i$.
*Hint:* consider the renewal process formed by returns to state $i$.

---

## Computations (paper)

**5.** A simplified headcount model has 3 states: Working (W), Promoted (P,
absorbing), Exit (E, absorbing). The transition matrix is:

$$
P = \begin{pmatrix}
0.90 & 0.05 & 0.05 \\
0    & 1    & 0    \\
0    & 0    & 1
\end{pmatrix}
$$

Compute the fundamental matrix $N$, expected absorption time, and probability
of being promoted vs exiting.
*Note:* 1×1 transient submatrix, so $N$ is a scalar.

**6.** For the full 4-state headcount model (Exit absorbing): write the
canonical form, identify $Q$ and $R$, and compute $N = (I - Q)^{-1}$.
*Requires inverting a 3×3 matrix.* Compute the expected time to exit for each
starting level.

**7.** Compute the first passage time from Junior to Senior in the full model
(with Exit absorbing). What fraction of Juniors eventually reach Senior?
*Use $B = NR$.*

---

## Self-check

- [ ] $N$ has all entries $\geq 0$
- [ ] $t_i = \sum_j N_{ij}$ is strictly positive and finite
- [ ] Rows of $B$ sum to 1 (every transient state ends up absorbed with probability 1)
- [ ] The first passage probability Junior → Senior matches the Senior column of $B$
      for the variant that treats Senior and Exit as competing absorbing states
