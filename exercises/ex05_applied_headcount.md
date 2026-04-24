# Exercise Set 5 — Applied Headcount Model

## Computations (paper)

**1.** A team has 45 employees: 15 Juniors, 20 Mids, 10 Seniors. Monthly
transition rates: Junior promotion 3%, Junior attrition 4%, Mid promotion 2%,
Mid attrition 2%, Senior attrition 1%. Replacement hiring fills Junior
positions. Compute the expected team composition after 1, 3, and 12 months.
*Use $\pi_n = \pi_0 P^n$.*

**2. Hiring freeze scenario.** Using the same rates but setting all hiring to
zero (Exit becomes absorbing), compute:

- Expected time until team drops below 40
- Expected time until all employees have left
- Probability the last person standing is a Senior

**3. Growth target.** Birth-death model with $\lambda = 3$, $\mu = 0.04$ per
person. Team starts at $n = 45$, target is $n = 55$. Solve the ODE
$dn/dt = \lambda - \mu n$ to find $\mathbb{E}[T_{\text{target}}]$. Compare with
steady-state mean $\lambda/\mu$. Is the target above or below steady state?

**4. Budget bridge.** Headcount at month 12 is approximately
$n_{12} \sim \text{Poisson}(75)$. Average salary
$S \sim \text{LogNormal}(9.2, 0.09)$. Using the law of total expectation and
total variance:

- Compute $\mathbb{E}[\text{Total Salary}] = \mathbb{E}[n_{12}] \cdot \mathbb{E}[S] \cdot 12$
- Compute $\mathrm{Var}[\text{Total Salary}]$
- *This connects Article 3 (headcount process) to Article 1 (budget simulation).*

---

## Interpretation prompts

**A.** In exercise 2, if replacement hiring at rate $p_{41} = 0.5$ is
reintroduced, what changes qualitatively? *Hint: the chain becomes irreducible
and the absorption analysis no longer applies directly — the stationary
analysis of Set 2 does.*

**B.** In exercise 3, what is the hitting-time interpretation of the
"half-life" $\ln 2 / \mu$? How many months does the gap $|n(t) - \lambda/\mu|$
take to halve?

**C.** Sensitivity: in exercise 1, which of $\{p_{12}, p_{14}, p_{24}, p_{34}\}$
has the biggest impact on the month-12 Senior count? Sketch a tornado chart
qualitatively before computing.

---

## Self-check

- [ ] Total composition after month 1 sums to 45 only if replacement hiring
      exactly compensates attrition — otherwise it drifts
- [ ] In exercise 2, probabilities for "last person is Senior" come from
      absorption probabilities conditioned on being the last non-Exit state
- [ ] ODE solution from Set 4 re-used in exercise 3 — check you get the same answer
- [ ] $\mathrm{Var}[\text{Total}] = 144 \cdot (\mathbb{E}[n]^2 \mathrm{Var}[S] + \mathbb{E}[S]^2 \mathrm{Var}[n] + \mathrm{Var}[n]\mathrm{Var}[S])$
