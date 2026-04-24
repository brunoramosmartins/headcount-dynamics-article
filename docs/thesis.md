# Thesis — The Team That Replaces Itself

## Central Claim (v0.1)

> A headcount plan that says "we will have 50 people by December" is a point
> forecast of a stochastic process — it ignores the probability of actually
> reaching that state, the expected time to get there, and the steady-state
> the system naturally tends toward. Markov chain theory provides exact answers
> to all three questions: the stationary distribution tells you where the team
> will settle, the fundamental matrix tells you how long transitions take, and
> the birth-death extension captures the continuous flow of hiring and attrition.

## Central Axis

Headcount is not a number — it is a stochastic process. Planning should work
*with* the process, not against it.

```
Deterministic plan:  "We will have 50 people in Q4."
Stochastic model:    "P(n ≥ 50 in Q4) = 0.73, E[T to reach 50] = 4.2 months."
```

## Scope

- Discrete-time Markov chains (DTMC): month-to-month career-level transitions
- Stationary distributions, convergence, spectral analysis
- Absorption analysis: fundamental matrix, hitting times, absorption probabilities
- Continuous-time Markov chains (CTMC) and Kolmogorov equations
- Birth-death processes as a model of team size evolution
- Applied headcount scenarios: growth, freeze, layoff, composition, seasonal hiring
- Connection to budget simulation (link to Article 1 and Article 2)

## Anti-Scope

- Queueing networks (Jackson, BCMP) — beyond a single birth-death queue
- Semi-Markov processes
- Hidden Markov models
- Markov Chain Monte Carlo (MCMC) — the chain is the model, not a sampler
- Real company data — all parameters are synthetic and plausible

## Audience

- Analytics engineers and data scientists who want a rigorous but applied
  introduction to stochastic processes
- IT / FP&A planners who already know forecasting and want probabilistic tools
- Readers of Articles 1 (Monte Carlo) and 2 (Distributions) completing the
  trilogy

## Abstract

Workforce headcount is usually planned as a single number per period. This
article argues the number is a thin projection of a richer object — a
stochastic process — and that working directly with the process yields answers
that deterministic planning cannot provide: the probability of hitting a
target, the expected time to reach it, the long-run team composition, and the
sensitivity of all three to transition parameters. Starting from the
definition of a Markov chain and proving the Chapman-Kolmogorov equation, the
article builds up to the fundamental matrix for absorption analysis, then
extends to continuous-time birth-death processes as a model of hiring and
attrition flows. Five applied scenarios (steady growth, hiring freeze, layoff
recovery, composition shift, seasonal hiring) illustrate the practical payoff.
The final section shows how the headcount distribution feeds the budget
simulation of Article 1, completing the probabilistic toolkit for IT budget
planning.
