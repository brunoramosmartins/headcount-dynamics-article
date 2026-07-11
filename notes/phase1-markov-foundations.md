# Phase 1 — Markov Chain Foundations

## 1. The Markov Property

Let $(X_n)_{n \geq 0}$ be a sequence of random variables taking values in a finite (or countable) state space $S = \{1, 2, \ldots, K\}$. The sequence is a **Markov chain** if, for every $n \geq 0$ and every $i_0, i_1, \ldots, i_{n+1} \in S$ such that $P(X_n = i_n, \ldots, X_0 = i_0) > 0$,

$$
P(X_{n+1} = j \mid X_n = i, X_{n-1} = i_{n-1}, \ldots, X_0 = i_0)  =  P(X_{n+1} = j \mid X_n = i)  =  p_{ij}.
$$

The conditional distribution of the next state depends on the present state alone — the past adds no information beyond what $X_n$ already carries. The numbers $p_{ij}$ are the **one-step transition probabilities**, assembled into the **transition matrix**

$$
P = (p_{ij})_{i,j \in S}.
$$

Throughout this article the chain is assumed **time-homogeneous**: $p_{ij}$ does not depend on $n$.

### 1.1 The transition matrix is row-stochastic

**Claim.** $p_{ij} \geq 0$ for all $i, j$, and $\sum_{j \in S} p_{ij} = 1$ for every $i$.

**Proof.** Non-negativity is immediate from the definition of probability. For each fixed $i$,

$$
\sum_{j \in S} p_{ij}  =  \sum_{j \in S} P(X_{n+1} = j \mid X_n = i)  =  P\!\left(\bigcup_{j \in S} \{X_{n+1} = j\} \,\Big|\, X_n = i\right)  =  1,
$$

since the events $\{X_{n+1} = j\}$ for $j \in S$ form a partition of the sample space (the chain must be in some state at time $n+1$). $\blacksquare$

Equivalently, $P \mathbf{1} = \mathbf{1}$ where $\mathbf{1} = (1, \ldots, 1)^T$. Each row of $P$ is a probability distribution over the next state.

---

## 2. $n$-Step Transitions and Chapman–Kolmogorov

### 2.1 Definition

The **$n$-step transition probability** is

$$
p_{ij}^{(n)}  =  P(X_n = j \mid X_0 = i).
$$

By convention $p_{ij}^{(0)} = \delta_{ij}$ (Kronecker delta) and $p_{ij}^{(1)} = p_{ij}$.

### 2.2 Theorem ($n$-step matrix)

**Claim.** $p_{ij}^{(n)} = (P^n)_{ij}$ for every $n \geq 0$.

**Proof (induction on $n$).** The base case $n = 0$ holds by convention. Assume $p_{ij}^{(n)} = (P^n)_{ij}$ for some $n \geq 0$. Conditioning on $X_n$,

$$
p_{ij}^{(n+1)}  =  P(X_{n+1} = j \mid X_0 = i)  =  \sum_{k \in S} P(X_{n+1} = j, X_n = k \mid X_0 = i).
$$

Apply the multiplication rule and the Markov property:

$$
P(X_{n+1} = j, X_n = k \mid X_0 = i)  =  P(X_{n+1} = j \mid X_n = k) \cdot P(X_n = k \mid X_0 = i)  =  p_{kj} \cdot p_{ik}^{(n)}.
$$

Hence

$$
p_{ij}^{(n+1)}  =  \sum_{k \in S} p_{ik}^{(n)} \, p_{kj}  =  \sum_{k \in S} (P^n)_{ik} \, P_{kj}  =  (P^{n+1})_{ij}. \quad\blacksquare
$$

### 2.3 Theorem (Chapman–Kolmogorov)

**Claim.** For every $m, n \geq 0$ and every $i, j \in S$,

$$
p_{ij}^{(m+n)}  =  \sum_{k \in S} p_{ik}^{(m)} \, p_{kj}^{(n)}.
$$

**Proof.** Conditioning on $X_m$, then applying the Markov property and time-homogeneity,

$$
p_{ij}^{(m+n)}
= P(X_{m+n} = j \mid X_0 = i)
= \sum_{k} P(X_{m+n} = j, X_m = k \mid X_0 = i)
$$

$$
= \sum_{k} P(X_{m+n} = j \mid X_m = k) \cdot P(X_m = k \mid X_0 = i)
= \sum_{k} p_{kj}^{(n)} \, p_{ik}^{(m)}. \quad\blacksquare
$$

In matrix form: $P^{m+n} = P^m P^n$. The Chapman–Kolmogorov equation is the probabilistic content of matrix multiplication associativity.

### 2.4 Corollary ($P^n$ is row-stochastic)

By induction, $(P^{n+1}) \mathbf{1} = P^n (P \mathbf{1}) = P^n \mathbf{1} = \mathbf{1}$. Every $n$-step transition matrix is itself row-stochastic.

### 2.5 Distribution evolution

If $\pi_n = (P(X_n = 1), \ldots, P(X_n = K))$ is the row vector of marginal probabilities at time $n$, then

$$
\pi_{n+1}  =  \pi_n \, P, \qquad \pi_n  =  \pi_0 \, P^n.
$$

The distribution flows by **right multiplication by $P$**. This is the convention used throughout the article.

---

## 3. State Classification

Fix the chain $\{X_n\}$ with transition matrix $P$.

### 3.1 Accessibility and communication

State $j$ is **accessible from** state $i$, written $i \to j$, if $p_{ij}^{(n)} > 0$ for some $n \geq 0$.

States $i$ and $j$ **communicate**, written $i \leftrightarrow j$, if both $i \to j$ and $j \to i$.

### 3.2 Theorem (communication is an equivalence relation)

**Claim.** The relation $\leftrightarrow$ on $S$ is reflexive, symmetric, and transitive.

**Proof.**

- **Reflexivity.** $p_{ii}^{(0)} = 1 > 0$, so $i \to i$. Hence $i \leftrightarrow i$.
- **Symmetry.** Immediate from the definition.
- **Transitivity.** Suppose $i \leftrightarrow j$ and $j \leftrightarrow k$. There exist $m_1, n_1 \geq 0$ with $p_{ij}^{(m_1)} > 0$ and $p_{ji}^{(n_1)} > 0$, and $m_2, n_2 \geq 0$ with $p_{jk}^{(m_2)} > 0$ and $p_{kj}^{(n_2)} > 0$. By Chapman–Kolmogorov,

  $$
  p_{ik}^{(m_1 + m_2)}  =  \sum_{\ell} p_{i\ell}^{(m_1)} \, p_{\ell k}^{(m_2)}  \geq  p_{ij}^{(m_1)} \, p_{jk}^{(m_2)}  >  0.
  $$

  Symmetrically $p_{ki}^{(n_2 + n_1)} \geq p_{kj}^{(n_2)} p_{ji}^{(n_1)} > 0$. Hence $i \leftrightarrow k$. $\blacksquare$

The equivalence classes of $\leftrightarrow$ are the **communication classes**. The chain is **irreducible** if there is exactly one communication class — every state is reachable from every other.

### 3.3 Recurrence and transience

For a state $i$, define the **return time** $T_i = \min\{n \geq 1 : X_n = i\}$ (with $T_i = \infty$ if no return). Let $f_i = P(T_i < \infty \mid X_0 = i)$.

State $i$ is **recurrent** if $f_i = 1$ and **transient** if $f_i < 1$.

Equivalent characterisation:

- $i$ is recurrent $\iff \sum_{n=0}^{\infty} p_{ii}^{(n)} = \infty$.
- $i$ is transient $\iff \sum_{n=0}^{\infty} p_{ii}^{(n)} < \infty$.

A useful fact (proof sketched in references): all states in the same communication class are of the same type. In particular, in an irreducible finite chain every state is recurrent.

### 3.4 Absorbing states

A state $i$ is **absorbing** if $p_{ii} = 1$, equivalently if $\{i\}$ is a communication class with $i \not\to j$ for any $j \neq i$. An **absorbing chain** is one in which every state is either absorbing or eventually reaches an absorbing state with positive probability.

### 3.5 Periodicity

The **period** of state $i$ is

$$
d(i)  =  \gcd\{n \geq 1 : p_{ii}^{(n)} > 0\}.
$$

State $i$ is **aperiodic** if $d(i) = 1$, **periodic** otherwise. Periodicity is a class property: if $i \leftrightarrow j$ then $d(i) = d(j)$.

---

## 4. Application: the Headcount Chain

Take Formulation A from the model design with replacement hiring, $p_{41} = 0.5$:

$$
P = \begin{pmatrix}
0.93 & 0.03 & 0    & 0.04 \\
0    & 0.96 & 0.02 & 0.02 \\
0    & 0    & 0.99 & 0.01 \\
0.50 & 0    & 0    & 0.50
\end{pmatrix}, \qquad S = \{\mathrm{J}, \mathrm{M}, \mathrm{S}, \mathrm{E}\}.
$$

### 4.1 Accessibility graph

| from \ to | J | M | S | E |
|-----------|---|---|---|---|
| J         | ✓ | ✓ |   | ✓ |
| M         |   | ✓ | ✓ | ✓ |
| S         |   |   | ✓ | ✓ |
| E         | ✓ |   |   | ✓ |

Direct arrows: $\mathrm{J} \to \mathrm{M}$, $\mathrm{M} \to \mathrm{S}$, $\mathrm{J,M,S} \to \mathrm{E}$, $\mathrm{E} \to \mathrm{J}$.

### 4.2 Communication classes

- $\mathrm{J} \to \mathrm{M}$ via $p_{12}$. Reverse: $\mathrm{M} \to \mathrm{S} \to \mathrm{E} \to \mathrm{J}$ has positive probability, so $\mathrm{M} \to \mathrm{J}$. Hence $\mathrm{J} \leftrightarrow \mathrm{M}$.
- $\mathrm{M} \to \mathrm{S}$ via $p_{23}$. Reverse: $\mathrm{S} \to \mathrm{E} \to \mathrm{J} \to \mathrm{M}$, so $\mathrm{S} \to \mathrm{M}$.
- $\mathrm{S} \to \mathrm{E}$ via $p_{34}$. Reverse: $\mathrm{E} \to \mathrm{J} \to \mathrm{M} \to \mathrm{S}$, so $\mathrm{E} \to \mathrm{S}$.

All four states communicate. The chain is **irreducible**. Every state is recurrent. Aperiodicity follows from $p_{ii} > 0$ for every $i$ (each state has a self-loop), so $d(i) = 1$.

### 4.3 The absorbing variant

Set $p_{41} = 0$, $p_{44} = 1$. Now $\mathrm{E}$ is absorbing. Communication classes become $\{\mathrm{J}, \mathrm{M}, \mathrm{S}\}$ no longer communicate (because $\mathrm{E}$ no longer routes back), and each transient class is in fact a singleton: $\{\mathrm{J}\}, \{\mathrm{M}\}, \{\mathrm{S}\}, \{\mathrm{E}\}$. The transient states $\mathrm{J}, \mathrm{M}, \mathrm{S}$ all reach $\mathrm{E}$ in finite expected time (Phase 3 makes this quantitative).

---

## 5. Worked Example — 3-state chain

Take a simplified chain with states $\{\mathrm{J}, \mathrm{M}, \mathrm{S}\}$ (absorbing variant collapsed by ignoring exits):

$$
P = \begin{pmatrix}
0.7 & 0.3 & 0   \\
0   & 0.8 & 0.2 \\
0   & 0   & 1
\end{pmatrix}.
$$

### 5.1 Compute $P^2$

$$
P^2 = \begin{pmatrix}
0.7 & 0.3 & 0 \\
0 & 0.8 & 0.2 \\
0 & 0 & 1
\end{pmatrix}
\begin{pmatrix}
0.7 & 0.3 & 0 \\
0 & 0.8 & 0.2 \\
0 & 0 & 1
\end{pmatrix}
= \begin{pmatrix}
0.49 & 0.45 & 0.06 \\
0    & 0.64 & 0.36 \\
0    & 0    & 1
\end{pmatrix}.
$$

Sanity checks: rows sum to $0.49 + 0.45 + 0.06 = 1$, $0.64 + 0.36 = 1$, $1$. ✓

Interpretation:

- $p_{\mathrm{JJ}}^{(2)} = 0.49$: a Junior is still Junior after 2 months with probability 49%.
- $p_{\mathrm{JM}}^{(2)} = 0.45$: a Junior is Mid after 2 months with probability 45%.
- $p_{\mathrm{JS}}^{(2)} = 0.06$: a Junior is Senior after 2 months with probability 6%.

### 5.2 Compute $P^3 = P^2 \cdot P$

$$
P^3 = \begin{pmatrix}
0.49 & 0.45 & 0.06 \\
0    & 0.64 & 0.36 \\
0    & 0    & 1
\end{pmatrix}
\begin{pmatrix}
0.7 & 0.3 & 0 \\
0 & 0.8 & 0.2 \\
0 & 0 & 1
\end{pmatrix}
= \begin{pmatrix}
0.343 & 0.507 & 0.150 \\
0     & 0.512 & 0.488 \\
0     & 0     & 1
\end{pmatrix}.
$$

After 3 months the Junior probability has fallen from 0.7 to 0.343, while the Senior probability has grown from 0 to 0.150.

### 5.3 State classification of the 3-state chain

- $\mathrm{S}$ is absorbing: $p_{\mathrm{SS}} = 1$.
- $\mathrm{J}, \mathrm{M}$ are transient: from either, the chain reaches $\mathrm{S}$ with positive probability and never returns.

The chain is **not irreducible**; it has three communication classes: $\{\mathrm{J}\}$, $\{\mathrm{M}\}$, $\{\mathrm{S}\}$.

---

## 6. Summary

- A Markov chain is fully specified by its state space $S$ and transition matrix $P$.
- $P$ is row-stochastic; $P^n$ is row-stochastic for all $n \geq 0$.
- Chapman–Kolmogorov: $P^{m+n} = P^m P^n$. Distributions evolve via $\pi_n = \pi_0 P^n$.
- Communication is an equivalence relation. Recurrence and periodicity are class properties.
- The headcount chain with replacement hiring is irreducible, aperiodic, finite — and therefore positive recurrent (Phase 2 will compute its stationary distribution).

## References

- Norris, *Markov Chains*, Cambridge UP, 1997 — Ch. 1.
- Ross, *Introduction to Probability Models*, Academic Press, 11th ed., 2014 — Ch. 4.
- Grinstead & Snell, *Introduction to Probability*, AMS, 1997 — Ch. 11.
