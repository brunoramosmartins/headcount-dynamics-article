# Exercise Set 4 — Continuous-Time and Birth-Death

## Proofs (paper)

**1.** Prove that the generator matrix $Q$ of a CTMC has row sums equal to
zero.
*Hint:* from the definition of transition rates and the requirement $P(t)$ is
row-stochastic for all $t \geq 0$. Differentiate $P(t)\mathbf{1} = \mathbf{1}$
at $t = 0$.

**2.** Derive the Kolmogorov forward equation $P'(t) = P(t)Q$ from the
Chapman-Kolmogorov equation $P(t+h) = P(t)P(h)$.
*Technique:* take the limit $h \to 0$, use $P(h) = I + Qh + o(h)$.

**3.** Prove that $P(t) = e^{Qt}$ satisfies both forward and backward
Kolmogorov equations.
*Technique:* differentiate the matrix exponential series.

**4.** Derive the detailed balance equations for a birth-death process:
$\pi_n \lambda_n = \pi_{n+1} \mu_{n+1}$.
*Hint:* start from $\pi Q = 0$ and use the tridiagonal structure of $Q$.

**5.** Prove that for a birth-death process with $\lambda_n = \lambda$
(constant) and $\mu_n = n\mu$ (linear), the stationary distribution is
$\pi_n = e^{-\rho} \rho^n / n!$ where $\rho = \lambda/\mu$.
*Technique:* solve the detailed balance recurrence and identify the Poisson PMF.

---

## Computations (paper)

**6.** A team has constant hiring rate $\lambda = 3$ (people/month) and
per-capita attrition $\mu = 0.04$ (per person/month). Compute the steady-state
mean team size. Compute $P(n \geq 80)$ and $P(n \leq 60)$ at steady state.

**7.** Same parameters, starting at $n_0 = 45$: the expected trajectory
satisfies $\frac{dn}{dt} = \lambda - \mu n$. Solve this ODE. What is
$\mathbb{E}[n(6)]$? $\mathbb{E}[n(12)]$? When does $\mathbb{E}[n(t)]$ reach 70?

**8.** Birth-death with $\lambda_n = \lambda$ and $\mu_n = \mu$ (both
constant). Is the chain positive recurrent? Compare with a random walk. Under
what condition on $\lambda/\mu$ does a stationary distribution exist?

---

## Self-check

- [ ] Rows of $Q$ sum to 0 (not 1 as in DTMC)
- [ ] $e^{Qt}$ is always row-stochastic when $Q$ is a generator
- [ ] Poisson $\pi_n$ sums to 1 (use the series $\sum \rho^n / n! = e^\rho$)
- [ ] ODE solution: $n(t) = \lambda/\mu + (n_0 - \lambda/\mu) e^{-\mu t}$
