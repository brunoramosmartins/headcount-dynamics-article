# Exercise Set 2 — Stationary Distributions and Convergence

## Proofs (paper)

**1.** Prove that for a finite irreducible Markov chain, the stationary
distribution $\pi$ exists.
*Hint:* consider the Cesàro average $\bar{\pi}_n = \frac{1}{n}\sum_{k=0}^{n-1} e_i P^k$
and show that any limit point satisfies $\pi P = \pi$. Use compactness of the
probability simplex.

**2.** Prove that if $P$ is irreducible and aperiodic, then $\pi$ is unique.
*Hint:* suppose $\pi'$ is another stationary distribution. Show
$\pi' = \pi' P^n \to \pi' \cdot \mathbf{1} \pi^T = \pi$.

**3.** Prove that the eigenvalues of a row-stochastic matrix satisfy
$|\lambda| \leq 1$.
*Hint:* $\|P\|_\infty = 1$ (maximum absolute row sum). If $Pv = \lambda v$, take
$\|\cdot\|_\infty$.

**4.** Prove that the rate of convergence
$\|P^n - \mathbf{1}\pi^T\|$ is governed by $|\lambda_2|^n$ where $\lambda_2$ is
the second-largest eigenvalue in modulus.
*Hint:* decompose $P$ spectrally; all terms except the rank-1 projection onto
$\pi$ decay exponentially.

---

## Computations (paper)

**5.** For the transition matrix of Exercise Set 1 (Exercise 5), solve
$\pi P = \pi$ with $\sum \pi_i = 1$ by hand. Set up 4 equations, note one is
redundant, substitute the normalization.

**6.** Compute the eigenvalues of the same $P$.
*Hint:* $\lambda = 1$ is always an eigenvalue. For a $4 \times 4$ matrix, use
the characteristic polynomial or numerical estimation. What is the spectral
gap $1 - |\lambda_2|$? Estimate the mixing time
$t_{mix} \approx \frac{1}{1 - |\lambda_2|}$.

**7.** If the Junior attrition rate $p_{14}$ increases from 0.04 to 0.08, what
happens to the stationary distribution? Recompute $\pi$ and compare. What is
the business interpretation?

---

## Self-check

- [ ] $\pi$ satisfies $\pi_i \geq 0$ and $\sum \pi_i = 1$
- [ ] Eigenvalue 1 has left eigenvector proportional to $\pi$
- [ ] Spectral gap > 0 means the chain mixes
- [ ] Business interpretation of (7): increasing junior churn shifts mass toward
      recycling-through-Exit — fewer Seniors in equilibrium
