# TIL — A Headcount Plan Is a Point Estimate of a Distribution

**Phase:** 5 · **Topic:** applied modelling · **Domain:** workforce planning

## Hook

"We will have 50 people in Q4" is a useful answer. It is also a lossy
compression of the real object the planner cares about.

## Insight

The real object is a **distribution** over team sizes at Q4. It has a mean
(which the plan sort of tracks), a variance (which the plan hides entirely),
and a shape (which tells you the probability of over- or under-shooting).
Three questions the single number cannot answer:

1. What is $P(n_{Q4} \geq 50)$?
2. What is $\mathbb{E}[T \mid \text{first hit 50}]$?
3. Where does the team settle if current rates persist?

All three have closed-form answers in the stochastic model. None have answers
in the deterministic one.

## Example

Under the default birth-death parameters, the model says the team follows
$n(t) = 80 + (45 - 80)e^{-0.025 t}$ in expectation, with Poisson spread around
it. At $t = 12$ months, $\mathbb{E}[n] \approx 54$, $\sigma \approx 7$.
$P(n_{12} \geq 50) \approx 0.72$. A plan that says "50 people by December" is
betting on a 72% probability — which you can now size, hedge, or accept.

## Takeaway

When you model the process, the "plan" becomes one statistic among many you
can compute. The probability of hitting the plan becomes as important as the
plan itself. That is what probabilistic planning buys you.
