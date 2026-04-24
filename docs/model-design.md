# Headcount Model Design

Two complementary formulations model the same workforce from different angles.
Together they answer questions about composition, growth, absorption, and
time-to-target.

---

## Formulation A — Career Level States (DTMC)

Each employee is, at the start of each month, in exactly one state.

| State | Label  | Description              |
|-------|--------|--------------------------|
| 1     | Junior | Entry-level              |
| 2     | Mid    | Intermediate             |
| 3     | Senior | Experienced              |
| 4     | Exit   | Left the company         |

Transitions happen monthly: stay, promote, or exit. Optionally, `Exit → Junior`
models replacement hiring (recycling). Without recycling, `Exit` is absorbing.

### Transition Matrix

$$
P = \begin{pmatrix}
p_{11} & p_{12} & 0      & p_{14} \\
0      & p_{22} & p_{23} & p_{24} \\
0      & 0      & p_{33} & p_{34} \\
p_{41} & 0      & 0      & p_{44}
\end{pmatrix}
$$

- $p_{12}$: promotion Junior → Mid
- $p_{23}$: promotion Mid → Senior
- $p_{i4}$: attrition from level $i$
- $p_{41}$: replacement hiring (Exit → Junior); absorbing variant sets $p_{44}=1$
- No demotions (upper-triangular on career progression)

### Default Parameters

| Parameter               | Value        | Rationale                           |
|-------------------------|--------------|-------------------------------------|
| $n_{\text{initial}}$    | 45           | Current team size                   |
| $p_{12}$ (promo J→M)    | 0.03/month   | ~1 promotion per year from pool     |
| $p_{23}$ (promo M→S)    | 0.02/month   | Slower senior promotion             |
| $p_{14}$ (attrition J)  | 0.04/month   | Higher junior turnover              |
| $p_{24}$ (attrition M)  | 0.02/month   | Moderate mid attrition              |
| $p_{34}$ (attrition S)  | 0.01/month   | Low senior attrition                |

---

## Formulation B — Team Size States (Birth-Death)

The state is the total team size $n \in \{0, 1, 2, \ldots, N_{max}\}$.

| Event              | Rate          | Description                         |
|--------------------|---------------|-------------------------------------|
| Birth (hiring)     | $\lambda_n$   | New employee joins                  |
| Death (attrition)  | $\mu_n$       | Employee leaves                     |

### Generator Matrix

$$
Q = \begin{pmatrix}
-\lambda_0 & \lambda_0 & 0 & \cdots \\
\mu_1 & -(\lambda_1 + \mu_1) & \lambda_1 & \cdots \\
0 & \mu_2 & -(\lambda_2 + \mu_2) & \cdots \\
\vdots & & & \ddots
\end{pmatrix}
$$

### Default Parameters

| Parameter      | Value              | Rationale                              |
|----------------|--------------------|----------------------------------------|
| $\lambda$      | 2/month            | Steady hiring rate                     |
| $\mu$          | 0.025 × $n$        | Size-proportional attrition            |

Steady-state mean team size: $\lambda / \mu = 80$.
The stationary distribution is Poisson with mean $\lambda / \mu$.

---

## Five Scenarios

| ID | Name              | Description                                          |
|----|-------------------|------------------------------------------------------|
| S1 | Steady Growth     | Team grows from 45 to 60 over 12 months              |
| S2 | Hiring Freeze     | All hiring stops; model attrition trajectory         |
| S3 | Layoff Recovery   | Sudden reduction to 35; model recovery time          |
| S4 | Composition Shift | Increase senior ratio from 20% to 30%                |
| S5 | Seasonal Hiring   | $\lambda$ varies by quarter (high Q1/Q3, low Q2/Q4)  |

---

## Expansion Points

- **Multi-role model**: extend states to (level × role), e.g., `(Junior, Backend)`
- **Time-varying rates**: make $p_{ij}(t)$ or $\lambda(t)$, $\mu(t)$ functions of $t$
- **Coupled chains**: one employee's state affects hiring rate (budget-constrained)
- **Cost overlay**: assign cost to each state for direct budget linkage (Article 1)

---

## Questions the Model Answers

1. **Steady state** — What team composition does the system converge to?
2. **Hitting time** — Expected time to grow from 45 to 55 employees?
3. **Absorption** — If hiring freezes, expected time until team below 40?
4. **Sensitivity** — Which rate most impacts steady-state size?
5. **Mixing time** — How many months until initial composition is "forgotten"?
