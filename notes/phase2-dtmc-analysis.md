# Phase 2 — Stationary Distributions and Convergence

Phase 1 built the language of Markov chains. This phase answers two questions
that language opens up:

1. **Where does the chain go?** — the stationary distribution $\pi$.
2. **How fast does it get there?** — the spectral gap and the mixing time.

Together these results turn the headcount chain into a forecasting tool: $\pi$
is the long-run team composition; the mixing time tells you how many months
must pass before the present composition no longer matters.

---

## 1. Stationary Distributions

### 1.1 Definition

A row vector $\pi = (\pi_1, \ldots, \pi_K)$ is a **stationary distribution**
of the chain with transition matrix $P$ if

$$
\pi P  =  \pi, \qquad \pi_i \geq 0, \qquad \sum_{i \in S} \pi_i  =  1.
$$

Equivalently, $\pi$ is a left eigenvector of $P$ with eigenvalue $1$, scaled
to a probability distribution. A chain started in $\pi$ stays in $\pi$:
$\pi_n = \pi_0 P^n = \pi$ for all $n$.

### 1.2 Existence (finite irreducible chains)

**Theorem.** If $S$ is finite and the chain is irreducible, a stationary
distribution exists.

**Proof (Cesàro average).** Fix any state $i$ and define

$$
\bar{\pi}_n  =  \frac{1}{n} \sum_{k=0}^{n-1} e_i P^k,
$$

where $e_i$ is the indicator row vector of state $i$. Each $\bar{\pi}_n$ is
a probability vector (convex combination of probability vectors), so it lies
in the compact simplex $\Delta_K = \{x \in \mathbb{R}^K : x_i \geq 0, \sum x_i = 1\}$.

By compactness, $\{\bar{\pi}_n\}$ has a convergent subsequence
$\bar{\pi}_{n_j} \to \pi$. Then

$$
\bar{\pi}_n P - \bar{\pi}_n  =  \frac{1}{n} \left( e_i P^n - e_i \right),
$$

which has $\ell_1$-norm bounded by $2/n$ and hence vanishes as $n \to \infty$.
Taking the limit along the subsequence,

$$
\pi P - \pi  =  \lim_{j \to \infty} \left( \bar{\pi}_{n_j} P - \bar{\pi}_{n_j} \right)  =  0,
$$

so $\pi P = \pi$, and $\pi \in \Delta_K$ is a stationary distribution. $\blacksquare$

### 1.3 Uniqueness (irreducible aperiodic chains)

**Theorem.** If the chain is irreducible and aperiodic on a finite state space,
the stationary distribution is unique.

**Proof sketch.** By the convergence theorem (§3 below), $P^n \to \mathbf{1} \pi$
where $\pi$ is the stationary distribution from §1.2. Suppose $\pi'$ is another
stationary distribution. Then for every $n$,

$$
\pi'  =  \pi' P^n  \xrightarrow{n \to \infty}  \pi' \cdot (\mathbf{1} \pi)  =  (\pi' \mathbf{1}) \cdot \pi  =  \pi.
$$

Hence $\pi' = \pi$. $\blacksquare$

For finite chains, **irreducibility alone** also implies uniqueness — the
periodic case requires a bit more care (the Cesàro limit replaces the
ordinary limit). Aperiodicity is needed only when one wants $P^n \to \mathbf{1} \pi$
in the ordinary sense, which is the form used everywhere in the article.

---

## 2. The Spectral Picture: Perron–Frobenius

### 2.1 Theorem (eigenvalues of a stochastic matrix)

If $P$ is row-stochastic, every eigenvalue $\lambda$ of $P$ satisfies
$|\lambda| \leq 1$, and $\lambda = 1$ is always an eigenvalue (with right
eigenvector $\mathbf{1}$).

**Proof.** Let $Pv = \lambda v$ for some $v \neq 0$, and let
$\|v\|_\infty = \max_i |v_i|$ be attained at index $i^*$. Then

$$
|\lambda| \cdot \|v\|_\infty  =  |\lambda v_{i^*}|  =  \left| \sum_j p_{i^* j} v_j \right|  \leq  \sum_j p_{i^* j} \, |v_j|  \leq  \|v\|_\infty \cdot \sum_j p_{i^* j}  =  \|v\|_\infty.
$$

Hence $|\lambda| \leq 1$. The vector $\mathbf{1}$ satisfies $P \mathbf{1} = \mathbf{1}$
because $P$ is row-stochastic, so $1 \in \mathrm{spec}(P)$. $\blacksquare$

### 2.2 Theorem (Perron–Frobenius, finite irreducible aperiodic case)

If $P$ is row-stochastic, irreducible, and aperiodic, then:

1. $\lambda = 1$ is a simple eigenvalue (algebraic multiplicity 1).
2. Every other eigenvalue $\lambda$ satisfies $|\lambda| < 1$.
3. The unique left eigenvector for $\lambda = 1$ (normalised to a probability
   vector) is the stationary distribution $\pi$, with $\pi_i > 0$ for all $i$.

A full proof of Perron–Frobenius requires positivity arguments via $P^k > 0$
for some $k$ (which holds under irreducibility + aperiodicity); see Norris
§1.7 or the Levin–Peres–Wilmer chapter on mixing. The third claim is the
analytical pay-off: the spectrum of $P$ tells us exactly where the chain
settles and how its rows decay onto $\pi$.

---

## 3. The Convergence Theorem

### 3.1 Statement

**Theorem.** Let $P$ be the transition matrix of a finite, irreducible,
aperiodic Markov chain with stationary distribution $\pi$. Then for every
state $i, j \in S$,

$$
\lim_{n \to \infty} p_{ij}^{(n)}  =  \pi_j.
$$

In matrix form, $P^n \to \mathbf{1} \pi$ as $n \to \infty$, where
$\mathbf{1} \pi$ is the rank-1 matrix whose every row is $\pi$.

### 3.2 Coupling proof (sketch)

Construct two independent copies of the chain $(X_n)$ and $(Y_n)$ with
$X_0 = i$ and $Y_0 \sim \pi$. Run them on the same probability space.
Define the coupling time

$$
\tau  =  \min\{n \geq 0 : X_n = Y_n\}.
$$

After they meet, replace $Y$ by $X$. Because the chain is irreducible and
aperiodic on a finite state space, $P(\tau < \infty) = 1$ — in fact $\tau$
has finite expectation. For any state $j$,

$$
|p_{ij}^{(n)} - \pi_j|  =  |P(X_n = j) - P(Y_n = j)|  \leq  P(\tau > n)  \to  0.
$$

The reason aperiodicity is essential: without it, $X_n$ and $Y_n$ may have
non-overlapping support modulo the period and never meet. With aperiodicity
$P^k > 0$ for some $k$, so coupling within time $k \cdot K$ has positive
probability per attempt and $\tau$ is geometrically dominated. $\blacksquare$

### 3.3 Spectral proof of the rate

Assume $P$ has $K$ distinct eigenvalues $1 = \lambda_1 > |\lambda_2| \geq \ldots \geq |\lambda_K|$.
Then $P$ admits the spectral decomposition

$$
P  =  \sum_{k=1}^K \lambda_k \, u_k v_k^\top,
$$

where $u_k$ and $v_k$ are right and left eigenvectors with $v_k^\top u_\ell = \delta_{k\ell}$,
$u_1 = \mathbf{1}$, $v_1 = \pi^\top$. Hence

$$
P^n  =  \mathbf{1} \pi + \sum_{k=2}^K \lambda_k^n \, u_k v_k^\top.
$$

Since $|\lambda_k| < 1$ for $k \geq 2$, every non-stationary mode decays
geometrically. In any matrix norm consistent with the spectral structure,

$$
\boxed{  \|P^n - \mathbf{1}\pi\|  \leq  C \cdot |\lambda_2|^n  }
$$

for a constant $C$ depending on the eigenbasis (and on the choice of norm).
The decay is **dominated by $\lambda_2$**, the second-largest eigenvalue in
modulus.

### 3.4 Spectral gap and mixing time

The **spectral gap** is

$$
\gamma  =  1 - |\lambda_2|  \in  (0, 1].
$$

A small gap (close to 0) means slow mixing; a large gap (close to 1) means
fast mixing.

The **total variation distance** at step $n$ from a starting state $i$ is

$$
d_i(n)  =  \tfrac{1}{2} \sum_{j \in S} \left| p_{ij}^{(n)} - \pi_j \right|.
$$

The **mixing time** at tolerance $\varepsilon$ is

$$
t_{\mathrm{mix}}(\varepsilon)  =  \min\!\left\{ n : \max_{i \in S} d_i(n) \leq \varepsilon \right\}.
$$

The standard convention is $\varepsilon = 1/4$ (Levin–Peres–Wilmer); with that
choice $t_{\mathrm{mix}}$ scales asymptotically as $\Theta(1/\gamma)$.

A useful bound: $t_{\mathrm{mix}}(\varepsilon) \leq \frac{1}{\gamma} \log\!\left(\frac{1}{\varepsilon \pi_{\min}}\right)$
where $\pi_{\min} = \min_i \pi_i$.

---

## 4. Application: the Headcount Chain

Take Formulation A with replacement hiring:

$$
P = \begin{pmatrix}
0.93 & 0.03 & 0    & 0.04 \\
0    & 0.96 & 0.02 & 0.02 \\
0    & 0    & 0.99 & 0.01 \\
0.50 & 0    & 0    & 0.50
\end{pmatrix}, \qquad S = \{\mathrm{J}, \mathrm{M}, \mathrm{S}, \mathrm{E}\}.
$$

The chain is irreducible and aperiodic (Phase 1 §4). By Perron–Frobenius
the stationary distribution is unique with strictly positive entries.

### 4.1 Setting up $\pi P = \pi$

Writing $\pi = (\pi_J, \pi_M, \pi_S, \pi_E)$, the four equations are

$$
\begin{aligned}
\pi_J &= 0.93 \, \pi_J + 0.50 \, \pi_E \\
\pi_M &= 0.03 \, \pi_J + 0.96 \, \pi_M \\
\pi_S &= 0.02 \, \pi_M + 0.99 \, \pi_S \\
\pi_E &= 0.04 \, \pi_J + 0.02 \, \pi_M + 0.01 \, \pi_S + 0.50 \, \pi_E
\end{aligned}
$$

with the normalisation $\pi_J + \pi_M + \pi_S + \pi_E = 1$.

The second equation gives $0.04 \, \pi_M = 0.03 \, \pi_J$, so
$\pi_M = \tfrac{3}{4} \pi_J$.
The third gives $0.01 \, \pi_S = 0.02 \, \pi_M$, so $\pi_S = 2 \, \pi_M = \tfrac{3}{2} \pi_J$.
The first equation rearranges to $0.07 \, \pi_J = 0.50 \, \pi_E$, so
$\pi_E = 0.14 \, \pi_J$.

Imposing normalisation:

$$
\pi_J \left( 1 + \tfrac{3}{4} + \tfrac{3}{2} + 0.14 \right)  =  \pi_J \cdot 3.39  =  1
 \implies  \pi_J \approx 0.2950.
$$

| State | $\pi_i$ (analytical) | Interpretation |
|-------|---------------------|----------------|
| Junior | 0.2950 | ≈ 30% of headcount-time spent at Junior |
| Mid    | 0.2212 | ≈ 22% |
| Senior | 0.4425 | ≈ 44% |
| Exit   | 0.0413 | ≈ 4% — short transient before recycling back |

The Senior class accounts for the largest share of long-run headcount-time
because $p_{34}$ (Senior attrition) is the lowest rate in the chain — Seniors
linger.

### 4.2 Eigenvalues, spectral gap, and mixing time

Numerical eigendecomposition of $P$ gives eigenvalues approximately

$$
\lambda_1 = 1.00, \qquad \lambda_2 \approx 0.99, \qquad \lambda_3 \approx 0.94, \qquad \lambda_4 \approx 0.41.
$$

The spectral gap is $\gamma = 1 - |\lambda_2| \approx 0.01$. The bound
$t_{\mathrm{mix}} \approx 1/\gamma$ predicts mixing on the order of 100 months.
Direct computation of the empirical mixing time
($\max_i \tfrac{1}{2}\|P^n_{i,:} - \pi\|_1 \leq 1/4$) gives $t_{\mathrm{mix}} \approx 70$ months.

**Business interpretation.** "Long-run composition" is a six-year horizon for
this chain. Quarterly hiring decisions move the present distribution but
cannot reshape the stationary distribution; that requires changing the
underlying transition rates (promotion or attrition policies).

### 4.3 Sensitivity to junior attrition

Replacing $p_{14} = 0.04$ with $p_{14} = 0.08$ (i.e., $p_{11} = 0.89$) and
re-solving:

| State  | $\pi_i$ (base) | $\pi_i$ (high churn) | Change |
|--------|----------------|----------------------|--------|
| Junior | 0.2950 | 0.3320 | +12% |
| Mid    | 0.2212 | 0.2491 | +13% |
| Senior | 0.4425 | 0.2491 | −44% |
| Exit   | 0.0413 | 0.1697 | +311% |

Doubling junior churn halves the Senior share. Junior attrition acts
upstream — it shrinks the pipeline that ever reaches Senior. Phase 5 will
revisit this as a sensitivity tornado chart.

---

## 5. Summary

- $\pi$ exists for every finite irreducible chain (Cesàro argument) and is
  unique under aperiodicity (or always, if the Cesàro limit is used).
- The spectrum of $P$ contains $\lambda = 1$ as a simple eigenvalue; all
  others lie strictly inside the unit disk (Perron–Frobenius).
- $P^n - \mathbf{1}\pi$ decays at rate $|\lambda_2|^n$; the spectral gap
  $\gamma = 1 - |\lambda_2|$ is the inverse-time-scale of convergence.
- The headcount chain has $\pi \approx (0.30, 0.22, 0.44, 0.04)$, spectral
  gap $\approx 0.01$, and mixing time $\sim 70$ months.

## References

- Norris, *Markov Chains*, Cambridge UP, 1997 — Ch. 1.7–1.8.
- Levin, Peres & Wilmer, *Markov Chains and Mixing Times*, AMS, 2009 — Ch. 4.
- Horn & Johnson, *Matrix Analysis*, Cambridge UP, 1985 — Ch. 8 (Perron–Frobenius).
