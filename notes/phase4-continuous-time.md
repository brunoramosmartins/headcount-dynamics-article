# Phase 4 — Continuous-Time Markov Chains and Birth-Death Processes

Phases 1–3 modelled the headcount as a sequence of monthly snapshots. Phase 4
removes the snapshots: hiring and attrition are **flows**, not events that
land on the first of each month. The mathematical object is a
*continuous-time Markov chain* (CTMC). The applied object is the
*birth-death process*. The pay-off is a closed-form stationary distribution
for team size: a Poisson with mean $\lambda / \mu$.

---

## 1. Continuous-Time Markov Chains

### 1.1 Definition

Let $\{X_t\}_{t \geq 0}$ be a stochastic process on a finite (or countable)
state space $S$. It is a **time-homogeneous CTMC** if for every
$0 \leq s_1 < \cdots < s_k < s$ and states $i_1, \ldots, i_k, i, j$,

$$
P(X_{s+t} = j \mid X_s = i, X_{s_k} = i_k, \ldots, X_{s_1} = i_1)
 =  P(X_{s+t} = j \mid X_s = i)
 =  p_{ij}(t),
$$

where $p_{ij}(t)$ depends only on the elapsed time $t$. The function
$P(t) = (p_{ij}(t))$ is the **transition probability matrix at time $t$**.

Two structural facts are immediate:

- $P(0) = I$.
- $P(t)$ is row-stochastic for every $t \geq 0$.
- Chapman–Kolmogorov: $P(s + t) = P(s) P(t)$ for $s, t \geq 0$.

### 1.2 The generator matrix

Assume $P(t)$ is differentiable at $t = 0$. The **generator** is

$$
Q  =  P'(0)  =  \lim_{h \downarrow 0} \frac{P(h) - I}{h}.
$$

For $i \neq j$, $Q_{ij}$ is the **rate** at which the chain transitions from
$i$ to $j$. Specifically, given $X_t = i$, the probability of moving to $j$
during a small interval $(t, t + h]$ is $Q_{ij} h + o(h)$.

The diagonal entry $-Q_{ii}$ is the **total rate of leaving** state $i$.

#### Theorem (row sums of $Q$)

**Claim.** $\sum_j Q_{ij} = 0$ for every $i$.

**Proof.** Differentiate $P(t) \mathbf{1} = \mathbf{1}$ at $t = 0$:

$$
0  =  \frac{d}{dt} P(t) \mathbf{1} \,\Big|_{t=0}  =  P'(0) \mathbf{1}  =  Q \mathbf{1}. \quad\blacksquare
$$

So a generator $Q$ satisfies:

- $Q_{ij} \geq 0$ for $i \neq j$ (non-negative rates),
- $Q_{ii} = -\sum_{j \neq i} Q_{ij}$ (so each row sums to 0).

### 1.3 Holding times and embedded jumps

A useful intuitive picture: in state $i$, the chain waits an
$\mathrm{Exponential}(-Q_{ii})$ amount of time, then jumps to state $j$ with
probability $Q_{ij} / (-Q_{ii})$, independently of the holding time and the
past. Holding times are memoryless because the chain itself is.

---

## 2. Kolmogorov Equations

The matrix-valued ODE that governs $P(t)$.

### 2.1 The forward equation

#### Theorem ($P'(t) = P(t) Q$)

**Claim.** For all $t \geq 0$, $P(t)$ satisfies the **forward Kolmogorov equation**

$$
\frac{d}{dt} P(t)  =  P(t) \, Q.
$$

**Proof.** Apply Chapman–Kolmogorov at $s = t$, $u = h$:

$$
P(t + h) - P(t)  =  P(t) P(h) - P(t)  =  P(t) [P(h) - I].
$$

Divide by $h$ and let $h \downarrow 0$, using $P(h) = I + Qh + o(h)$:

$$
\frac{P(t + h) - P(t)}{h}  \to  P(t) Q. \quad\blacksquare
$$

### 2.2 The backward equation

By analogous reasoning (decompose $P(t+h) = P(h) P(t)$),

$$
\frac{d}{dt} P(t)  =  Q \, P(t).
$$

**Same equation, different reading.** The forward form is "given the chain
is at $i$ now and will be at $j$ at $t + h$, what happens in the last $h$?";
the backward form is "what happens in the first $h$?". The two derivations
are not automatic in the infinite case but coincide for finite $S$.

### 2.3 The matrix exponential solution

#### Theorem ($P(t) = e^{Qt}$)

**Claim.** With initial condition $P(0) = I$, the unique solution to either
Kolmogorov equation is

$$
P(t)  =  e^{Qt}  =  \sum_{k=0}^{\infty} \frac{(Qt)^k}{k!}.
$$

**Proof.** The series converges absolutely on every compact interval (since
$Q$ is bounded for finite $S$). Differentiating term-by-term,

$$
\frac{d}{dt} e^{Qt}
 =  \sum_{k=1}^\infty \frac{k \, Q^k t^{k-1}}{k!}
 =  Q \cdot \sum_{k=0}^\infty \frac{(Qt)^k}{k!}
 =  Q \, e^{Qt}
 =  e^{Qt} \, Q,
$$

where the last equality uses $Q$ commutes with its powers. So $e^{Qt}$ solves
both equations. Uniqueness follows from standard Picard–Lindelöf for matrix
ODEs. $\blacksquare$

A standard numerical technique is `scipy.linalg.expm`, which uses Padé
approximation with scaling-and-squaring; numerically robust for moderate
state spaces.

### 2.4 Stationary distribution

A row vector $\pi$ is **stationary** for the CTMC if

$$
\pi \, Q  =  0, \qquad \pi_i \geq 0, \qquad \sum_i \pi_i  =  1.
$$

Equivalently, $\pi e^{Qt} = \pi$ for all $t$ — once started in $\pi$, the
chain stays in $\pi$. Existence and uniqueness mirror the discrete case: a
finite irreducible CTMC has a unique stationary distribution. The convergence
theorem becomes $P(t) \to \mathbf{1} \pi$ as $t \to \infty$ at rate $e^{-\gamma t}$
where $\gamma$ is the spectral gap of $-Q$.

---

## 3. Birth-Death Processes

A **birth-death process** is a CTMC on $S = \{0, 1, 2, \ldots\}$ in which
transitions go only to neighbouring states. Concretely, given $X_t = n$:

- a **birth** ($n \to n + 1$) occurs at rate $\lambda_n \geq 0$;
- a **death** ($n \to n - 1$) occurs at rate $\mu_n \geq 0$ (with $\mu_0 = 0$).

The generator is tridiagonal:

$$
Q  = 
\begin{pmatrix}
-\lambda_0 & \lambda_0 & 0 & 0 & \cdots \\
\mu_1 & -(\lambda_1 + \mu_1) & \lambda_1 & 0 & \cdots \\
0 & \mu_2 & -(\lambda_2 + \mu_2) & \lambda_2 & \cdots \\
0 & 0 & \mu_3 & -(\lambda_3 + \mu_3) & \cdots \\
\vdots & & & & \ddots
\end{pmatrix}.
$$

### 3.1 Detailed balance

#### Theorem (detailed balance)

**Claim.** A stationary distribution $\pi$ of a birth-death process satisfies

$$
\pi_n \, \lambda_n  =  \pi_{n+1} \, \mu_{n+1} \qquad \text{for all } n \geq 0.
$$

**Proof.** Write out $\pi Q = 0$ component-wise. The $n$-th equation
($n \geq 1$) is

$$
\pi_{n-1} \lambda_{n-1} - \pi_n (\lambda_n + \mu_n) + \pi_{n+1} \mu_{n+1}  =  0.
$$

Rearrange:

$$
\bigl( \pi_{n+1} \mu_{n+1} - \pi_n \lambda_n \bigr) - \bigl( \pi_n \mu_n - \pi_{n-1} \lambda_{n-1} \bigr)  =  0.
$$

So the quantity $J_n = \pi_n \lambda_n - \pi_{n+1} \mu_{n+1}$ is constant in
$n$. The boundary equation at $n = 0$ ($-\pi_0 \lambda_0 + \pi_1 \mu_1 = 0$)
gives $J_0 = 0$. Hence $J_n = 0$ for every $n$. $\blacksquare$

**Reading.** "The probability flow from $n$ to $n + 1$ equals the flow from
$n + 1$ to $n$." Detailed balance is the continuous-time analog of the time-
reversibility argument: the rate of births out of $n$ exactly balances the
rate of deaths into $n$.

### 3.2 Solving the recurrence

From $\pi_{n+1} = (\lambda_n / \mu_{n+1}) \, \pi_n$, by induction:

$$
\boxed{  
\pi_n  =  \pi_0 \cdot \prod_{k=0}^{n-1} \frac{\lambda_k}{\mu_{k+1}}
  }
$$

The constant $\pi_0$ is fixed by the normalisation $\sum_{n \geq 0} \pi_n = 1$:

$$
\pi_0  =  \left( 1 + \sum_{n = 1}^{\infty} \prod_{k=0}^{n-1} \frac{\lambda_k}{\mu_{k+1}} \right)^{-1}.
$$

The chain is **positive recurrent** (a stationary distribution exists) iff
this sum converges.

### 3.3 The Poisson steady-state — M/M/$\infty$ queue

A specially clean case: $\lambda_n = \lambda$ (constant hiring rate) and
$\mu_n = n \mu$ (per-capita attrition with rate $\mu$ each). Then

$$
\prod_{k=0}^{n-1} \frac{\lambda_k}{\mu_{k+1}}  =  \prod_{k=0}^{n-1} \frac{\lambda}{(k+1)\mu}  =  \frac{\lambda^n}{n! \, \mu^n}  =  \frac{\rho^n}{n!}, \qquad \rho := \frac{\lambda}{\mu}.
$$

So $\pi_n = \pi_0 \cdot \rho^n / n!$. Normalisation gives
$\pi_0 \sum_{n \geq 0} \rho^n / n! = \pi_0 e^\rho = 1$, hence $\pi_0 = e^{-\rho}$
and

$$
\boxed{  \pi_n  =  e^{-\rho} \, \frac{\rho^n}{n!}  \sim  \mathrm{Poisson}(\rho).  }
$$

**Headcount interpretation.** Constant hiring at rate $\lambda$ per month and
proportional attrition at rate $\mu$ per person per month produce a team size
that, in the long run, is Poisson-distributed with mean $\rho = \lambda / \mu$
and variance also $\rho$. The standard deviation grows like $\sqrt{\rho}$, so
the *coefficient of variation* falls like $1/\sqrt{\rho}$ — bigger teams
fluctuate proportionally less. This is the same reason aggregating insurance
risks reduces variance.

### 3.4 ODE for the expected trajectory

Take expectations on both sides of the rate equation:

$$
\frac{d}{dt} \mathbb{E}[X_t]  =  \lambda - \mu \, \mathbb{E}[X_t].
$$

Solving with initial condition $\mathbb{E}[X_0] = n_0$:

$$
\boxed{  \mathbb{E}[X_t]  =  \rho + (n_0 - \rho) \, e^{-\mu t}.  }
$$

The trajectory exponentially relaxes to the steady-state mean $\rho$ with
**half-life** $\ln 2 / \mu$. With $\mu = 0.025$, the half-life is $\ln 2 / 0.025 \approx 27.7$
months — over two years. The team forgets its initial size at this rate.

---

## 4. Application to the Headcount Chain

### 4.1 Default parameters

| Parameter | Value | Interpretation |
|-----------|-------|----------------|
| $\lambda$ | 2 / month | Hiring rate |
| $\mu$ | 0.025 per person per month | Per-capita attrition rate |
| $\rho = \lambda / \mu$ | 80 | Steady-state mean team size |
| $n_0$ | 45 | Current team size |

### 4.2 Steady-state distribution

$\pi_n = \mathrm{Poisson}(80)$. Standard deviation $\sqrt{80} \approx 8.94$.
A 95% interval is roughly $[63, 98]$.

| Quantity | Value |
|----------|-------|
| $\mathbb{E}[X_\infty]$ | 80 |
| $\mathrm{sd}(X_\infty)$ | 8.94 |
| $P(X_\infty \geq 80)$ | $\approx 0.52$ (Poisson(80) is nearly symmetric) |
| $P(X_\infty < 60)$ | $\approx 0.013$ |

### 4.3 Expected trajectory

With $n_0 = 45$:

$$
\mathbb{E}[X_t]  =  80 + (45 - 80) \, e^{-0.025 t}  =  80 - 35 \, e^{-0.025 t}.
$$

| $t$ (months) | $\mathbb{E}[X_t]$ |
|-------------:|------------------:|
| 0 | 45 |
| 6 | 49.9 |
| 12 | 54.1 |
| 24 | 60.8 |
| 36 | 65.8 |
| $\infty$ | 80 |

Time to reach mean 70: solve $80 - 35 e^{-0.025 t} = 70 \Rightarrow t = \ln(35/10) / 0.025 \approx 50.1$ months.

### 4.4 $P(X_{12} \geq 50)$

This is not a steady-state question — at $t = 12$ the distribution has not
converged. To answer it, build the truncated generator $Q$ on
$\{0, 1, \ldots, N\}$ for some large $N$ (e.g. $N = 200$), compute
$P(t) = e^{Qt}$, and read off

$$
P(X_{12} \geq 50 \mid X_0 = 45)  =  \sum_{j \geq 50} (e^{Qt})_{45, j}.
$$

The implementation in `src/ctmc.py` does exactly this. Numerically (with
$N = 200$) we find $P(X_{12} \geq 50) \approx 0.80$ — that is the actual
probability of "having at least 50 people one year out" given today's 45.
A deterministic forecast that simply projects $\mathbb{E}[X_{12}] \approx 54$
hides the 20% chance of finishing below 50.

### 4.5 Connection to Phase 1 (uniformization)

A CTMC with bounded rates can be related to a DTMC via **uniformization**:
choose $\nu \geq \max_i (-Q_{ii})$ and define

$$
\tilde{P}  =  I + \frac{1}{\nu} Q.
$$

Then $\tilde{P}$ is row-stochastic, and

$$
e^{Qt}  =  e^{-\nu t} \sum_{k \geq 0} \frac{(\nu t)^k}{k!} \, \tilde{P}^k.
$$

The CTMC at time $t$ equals the DTMC $\tilde{P}$ run for a $\mathrm{Poisson}(\nu t)$
number of steps. This is a useful simulation device and a conceptual bridge
back to Phase 1.

---

## 5. Summary

- A CTMC is described by a generator $Q$ with non-negative off-diagonal
  entries and zero row sums.
- Transition probabilities solve $P'(t) = P(t) Q$ with closed-form
  $P(t) = e^{Qt}$.
- Stationary distributions satisfy $\pi Q = 0$.
- A birth-death process has tridiagonal $Q$; detailed balance gives
  $\pi_n = \pi_0 \prod_{k < n} \lambda_k / \mu_{k+1}$.
- The M/M/$\infty$ specialisation ($\lambda$ constant, $\mu_n = n\mu$) yields
  $\pi_n = \mathrm{Poisson}(\rho)$ with $\rho = \lambda / \mu$.
- Headcount with $\lambda = 2$, $\mu = 0.025$ has steady-state mean 80 and
  reaches mean 70 in about 50 months from a starting size of 45.

## References

- Norris, *Markov Chains*, Cambridge UP, 1997 — Ch. 2.
- Ross, *Introduction to Probability Models*, Academic Press, 11th ed., 2014 — Ch. 6.
- Taylor & Karlin, *An Introduction to Stochastic Modeling*, Academic Press, 3rd ed., 1998 — Ch. 6.
