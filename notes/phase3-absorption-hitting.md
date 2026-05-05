# Phase 3 — Absorption, Hitting Times, and the Fundamental Matrix

Phases 1 and 2 dealt with the *recurrent* picture: stationary distributions
and how fast a chain forgets its initial state. This phase deals with the
*transient* picture: how long does it take a chain to leave the transient
states for good, and what is the probability of leaving through a particular
absorbing state?

The single object that answers both questions is the **fundamental matrix**

$$
N  =  (I - Q)^{-1},
$$

where $Q$ is the transient-to-transient submatrix of $P$. The whole phase is
an unfolding of $N$.

---

## 1. Absorbing chains and the canonical form

### 1.1 Definitions

A state $i$ is **absorbing** if $p_{ii} = 1$. A chain is an **absorbing chain**
if (a) it has at least one absorbing state and (b) every transient state can
reach some absorbing state in finitely many steps with positive probability.

Order the states so that the $r$ transient states come first and the $K - r$
absorbing states come last. The transition matrix takes the **canonical form**

$$
P  =  \begin{pmatrix} Q & R \\ 0 & I \end{pmatrix},
$$

where

- $Q \in \mathbb{R}^{r \times r}$ is the transient-to-transient block,
- $R \in \mathbb{R}^{r \times (K-r)}$ is the transient-to-absorbing block,
- $0$ is the $(K-r) \times r$ zero block,
- $I$ is the $(K-r) \times (K-r)$ identity (absorbing states stay absorbing).

Each row of $P$ sums to 1, so each row of $(Q  R)$ also sums to 1.

### 1.2 Powers of $P$

Powers of the canonical form have a clean structure:

$$
P^n  =  \begin{pmatrix} Q^n & R_n \\ 0 & I \end{pmatrix},
\qquad R_n  =  \sum_{k=0}^{n-1} Q^k R.
$$

In particular, $(P^n)_{ij}$ for $i, j$ both transient is exactly $(Q^n)_{ij}$ —
the probability that the chain is in transient state $j$ at time $n$ given it
started in transient state $i$ **and has not yet been absorbed**.

---

## 2. The Fundamental Matrix

### 2.1 Theorem ($Q^n \to 0$)

**Claim.** Every eigenvalue of $Q$ satisfies $|\lambda| < 1$. Consequently
$Q^n \to 0$ as $n \to \infty$.

**Proof.** For a transient state $i$, the probability the chain has not been
absorbed by time $n$ is $\sum_{j \text{ transient}} (Q^n)_{ij}$. Because the
chain reaches an absorbing state with probability 1, this row sum tends to 0;
in particular every entry tends to 0. Hence $Q^n \to 0$ entrywise.

For the spectral statement: if $Qv = \lambda v$ with $|\lambda| \geq 1$, then
$Q^n v = \lambda^n v$ does not vanish, contradicting $Q^n \to 0$. So
$|\lambda| < 1$ for every eigenvalue, equivalently $\rho(Q) < 1$. $\blacksquare$

### 2.2 Theorem (Neumann series)

**Claim.** $I - Q$ is invertible and

$$
N  =  (I - Q)^{-1}  =  \sum_{k=0}^{\infty} Q^k.
$$

**Proof.** Since $\rho(Q) < 1$, the series $\sum_k Q^k$ converges absolutely
(in any matrix norm). Compute

$$
(I - Q) \sum_{k=0}^{n} Q^k  =  \sum_{k=0}^{n} Q^k - \sum_{k=1}^{n+1} Q^k  =  I - Q^{n+1}  \xrightarrow{n \to \infty}  I.
$$

Hence $\sum_{k=0}^\infty Q^k$ is a left inverse of $I - Q$. The same calculation
on the right gives a right inverse. Inverses agree, so
$(I-Q)^{-1} = \sum_{k=0}^\infty Q^k$. $\blacksquare$

### 2.3 Theorem ($N_{ij}$ counts expected visits)

**Claim.** For transient states $i, j$,

$$
N_{ij}  =  \mathbb{E}\!\left[ \,\#\{n \geq 0 : X_n = j\} \,\big|\, X_0 = i \right].
$$

**Proof.** Let $V_j = \sum_{n=0}^\infty \mathbf{1}\{X_n = j\}$ count the total
visits to $j$. Linearity of expectation gives

$$
\mathbb{E}[V_j \mid X_0 = i]  =  \sum_{n=0}^\infty P(X_n = j \mid X_0 = i)  =  \sum_{n=0}^\infty (Q^n)_{ij}  =  N_{ij}.
$$

The middle step uses the canonical-form computation $P(X_n = j \mid X_0 = i) = (Q^n)_{ij}$
for transient $i, j$ (the chain has not been absorbed). $\blacksquare$

The diagonal entry $N_{ii}$ is the expected number of times the chain visits
$i$ starting from $i$, including the initial visit. Hence $N_{ii} \geq 1$.

### 2.4 Expected time to absorption

**Theorem.** Let $T = \min\{n \geq 0 : X_n \text{ is absorbing}\}$. Then for
each transient state $i$,

$$
t_i  :=  \mathbb{E}[T \mid X_0 = i]  =  \sum_{j \text{ transient}} N_{ij},
$$

i.e. $t = N \mathbf{1}$.

**Proof.** Total time before absorption equals the total number of visits to
all transient states. Linearity:

$$
t_i  =  \mathbb{E}\!\left[ \sum_{j \text{ transient}} V_j \,\big|\, X_0 = i \right]  =  \sum_j N_{ij}. \quad \blacksquare
$$

### 2.5 Absorption probabilities

**Theorem.** For transient $i$ and absorbing $a$, let $B_{ia}$ be the
probability the chain is eventually absorbed at state $a$ given $X_0 = i$. Then

$$
B  =  N R.
$$

**Proof.** Decompose $B_{ia}$ by the time of first transition out of the
transient states. The chain visits transient $j$ a total of $N_{ij}$ times in
expectation (Theorem 2.3); each visit, it leaves to absorbing state $a$ with
probability $R_{ja}$. The events at distinct visits are disjoint, so

$$
B_{ia}  =  \sum_{j \text{ transient}} N_{ij} \, R_{ja}  =  (NR)_{ia}. \quad\blacksquare
$$

Each row of $B$ sums to 1: starting from a transient state, the chain is
absorbed somewhere with probability 1.

### 2.6 Variance of absorption time

A useful by-product (proof in Grinstead & Snell §11.2): the variance of $T$
given $X_0 = i$ equals the $i$-th entry of

$$
\mathrm{Var}(T)  =  (2N - I) \, t  -  t \odot t,
$$

where $\odot$ is element-wise product. We will use this to put confidence
bands on absorption-time forecasts.

---

## 3. Hitting Times for General Chains

### 3.1 Definitions

The **first passage time** from $i$ to $j$ is

$$
T_j  =  \min\{n \geq 1 : X_n = j\}.
$$

Define $h_j(i) = \mathbb{E}[T_j \mid X_0 = i]$ (with $h_j(j) = 1/\pi_j$ by the
mean-return-time theorem below).

### 3.2 Theorem (linear system for hitting times)

**Claim.** For $i \neq j$,

$$
h_j(i)  =  1 + \sum_{k \neq j} p_{ik} \, h_j(k),
$$

with $h_j(j) = 0$ when interpreted as the *time to first reach $j$ from $j$*
(treating the start at $j$ as an immediate hit). The boundary convention
matters; see the discussion below.

**Proof.** Condition on the first step $X_1$:

$$
h_j(i)  =  \mathbb{E}[T_j \mid X_0 = i]  =  \sum_{k} p_{ik} \, \mathbb{E}[T_j \mid X_0 = i, X_1 = k].
$$

By the Markov property, $\mathbb{E}[T_j \mid X_1 = k]$ equals $1$ if $k = j$
(the chain has already arrived) and $1 + h_j(k)$ otherwise. Substituting,

$$
h_j(i)  =  p_{ij} \cdot 1 + \sum_{k \neq j} p_{ik} (1 + h_j(k))  =  1 + \sum_{k \neq j} p_{ik} h_j(k). \quad\blacksquare
$$

Treating $j$ as absorbing (set its row to $\delta_j$) reduces hitting-time
computation to absorption analysis: the expected time to hit $j$ is the
expected absorption time of the modified chain.

### 3.3 Theorem (mean return time)

**Claim.** For an irreducible positive recurrent chain with stationary
distribution $\pi$,

$$
\mathbb{E}[T_i \mid X_0 = i]  =  \frac{1}{\pi_i}.
$$

**Proof sketch.** Let $T_i^{(1)} = T_i$ and $T_i^{(k+1)}$ be the time of the
$(k+1)$-th return to state $i$. The differences $T_i^{(k+1)} - T_i^{(k)}$ are
i.i.d. by the strong Markov property, with mean $\mu = \mathbb{E}[T_i \mid X_0 = i]$.
By the renewal theorem (or the strong law of large numbers applied to the
counting process of returns), the long-run fraction of time the chain spends
in $i$ is

$$
\lim_{n \to \infty} \frac{1}{n} \sum_{k=0}^{n-1} \mathbf{1}\{X_k = i\}  =  \frac{1}{\mu} \quad \text{a.s.}
$$

But by the ergodic theorem this limit also equals $\pi_i$. Hence
$\pi_i = 1/\mu$, i.e. $\mu = 1/\pi_i$. $\blacksquare$

This is the bridge between Phases 2 and 3. The stationary distribution and
the mean return times are *the same data*, expressed in two different
languages.

---

## 4. Worked Example — three-state absorbing chain

Take the simplified model with states $\{W, P, E\}$ where $W$ is Working
(transient) and $P, E$ are Promoted and Exit (absorbing):

$$
P = \begin{pmatrix}
0.90 & 0.05 & 0.05 \\
0    & 1    & 0    \\
0    & 0    & 1
\end{pmatrix}.
$$

### 4.1 Canonical form

$Q = (0.90)$ (a $1 \times 1$ scalar block), $R = (0.05    0.05)$.

### 4.2 Fundamental matrix

$$
N  =  (I - Q)^{-1}  =  (1 - 0.90)^{-1}  =  (0.10)^{-1}  =  10.
$$

A Working employee is in the Working state for 10 months in expectation
before being absorbed.

### 4.3 Absorption time

$t = N \mathbf{1} = 10$ months.

### 4.4 Absorption probabilities

$$
B  =  N R  =  10 \cdot (0.05    0.05)  =  (0.50    0.50).
$$

A Working employee is equally likely to end in Promoted or Exit. The 50–50
result is sensitive to the equal weights $0.05 = 0.05$; doubling either rate
shifts the probability proportionally.

---

## 5. Application — full headcount chain (Exit absorbing)

Set $p_{44} = 1$ (no replacement hiring). The states are $\{J, M, S, E\}$
with $E$ absorbing. The canonical form has

$$
Q  =  \begin{pmatrix} 0.93 & 0.03 & 0    \\ 0    & 0.96 & 0.02 \\ 0    & 0    & 0.99 \end{pmatrix}, \qquad
R  =  \begin{pmatrix} 0.04 \\ 0.02 \\ 0.01 \end{pmatrix}.
$$

### 5.1 Fundamental matrix

Inverting $I - Q$ (upper-triangular, easy by back-substitution):

$$
I - Q  =  \begin{pmatrix} 0.07 & -0.03 & 0    \\ 0    & 0.04 & -0.02 \\ 0    & 0    & 0.01 \end{pmatrix},
\qquad
N  =  \begin{pmatrix} 14.286 & 10.714 & 21.429 \\ 0     & 25.000 & 50.000 \\ 0     & 0     & 100.000 \end{pmatrix}.
$$

Each entry $N_{ij}$ is the expected number of months a starting employee
spends in level $j$ before exiting.

### 5.2 Expected time to exit

$$
t  =  N \mathbf{1}  =  \begin{pmatrix} 46.43 \\ 75.00 \\ 100.00 \end{pmatrix}.
$$

| Starting level | Expected months until exit |
|----------------|---------------------------:|
| Junior         | 46.4 |
| Mid            | 75.0 |
| Senior         | 100.0 |

A Senior expects roughly twice the tenure of a Junior. The Senior's $1/\mu$
intuition matches: Senior attrition is 1% per month, so expected time to exit
from Senior is $1/0.01 = 100$ months. Mid takes longer than Junior because
Mid has more time-to-decay before exit.

### 5.3 Absorption probabilities

Since $E$ is the only absorbing state, $B = NR = \mathbf{1}$ — every transient
state is absorbed at $E$ with probability 1. To say something interesting, we
need to make Senior also absorbing.

### 5.4 Junior → Senior: first passage and reach probability

Modify $P$ by treating both Senior and Exit as absorbing. New canonical form:

$$
Q'  =  \begin{pmatrix} 0.93 & 0.03 \\ 0 & 0.96 \end{pmatrix}, \qquad
R'  =  \begin{pmatrix} 0 & 0.04 \\ 0.02 & 0.02 \end{pmatrix}, \qquad \text{columns} = (S, E).
$$

Then

$$
N' = (I - Q')^{-1} = \begin{pmatrix} 14.286 & 10.714 \\ 0 & 25.000 \end{pmatrix}, \qquad
B' = N' R' = \begin{pmatrix} 0.214 & 0.786 \\ 0.500 & 0.500 \end{pmatrix}.
$$

| Starting | $P(\text{reach Senior first})$ | $P(\text{Exit first})$ |
|----------|-------------------------------:|-----------------------:|
| Junior   | 0.214 | 0.786 |
| Mid      | 0.500 | 0.500 |

Only ~21% of Juniors ever reach Senior under these rates; the rest exit first.
Of those that do reach Senior, the conditional expected first-passage time is
$t'_J / B'_{J,S} \approx 25/0.214 \approx 117$ months — a different question
that requires conditional expectation tools.

The unconditional expected first passage from Junior to "Senior or Exit" is
$t' = N'\mathbf{1} = (25.0, 25.0)^\top$ months.

### 5.5 Linking back to Phase 2

Phase 2 gave the stationary distribution under recycling: $\pi \approx
(0.30, 0.22, 0.44, 0.04)$. The mean return time to Senior under the same
recycling chain is $1/\pi_S \approx 2.26$ months in *Phase 2 sense* —
confusingly small because the recycling chain mixes. The interpretation is:
once the chain has reached its stationary distribution, gaps between Senior
visits are short on average; that does **not** contradict the 100-month
absorption time we just computed, which assumes no recycling.

---

## 6. Summary

- The fundamental matrix $N = (I - Q)^{-1}$ counts expected visits to
  transient states.
- Expected absorption time: $t = N \mathbf{1}$.
- Absorption probabilities: $B = N R$.
- First passage to a target state reduces to absorption analysis after
  marking the target as absorbing.
- Mean return time to state $i$ in a positive recurrent chain equals
  $1/\pi_i$ — the stationary distribution and return times are dual data.
- For the headcount chain with Exit absorbing: 46 / 75 / 100 expected months
  to exit from Junior / Mid / Senior; only 21% of Juniors ever reach Senior.

## References

- Grinstead & Snell, *Introduction to Probability*, AMS, 1997 — Ch. 11.
- Norris, *Markov Chains*, Cambridge UP, 1997 — Ch. 1.5.
- Kemeny & Snell, *Finite Markov Chains*, Springer, 1976 — Ch. 3.
