# TIL — A Team with Constant Hiring and Proportional Attrition Is Poisson-Sized

**Phase:** 4 · **Topic:** birth-death processes · **Domain:** workforce modelling

## Hook

Here is a surprising fact. If you hire at constant rate $\lambda$ per month
and each employee has an independent probability $\mu$ per month of leaving,
the steady-state team size is Poisson-distributed with mean $\lambda / \mu$.
The same distribution you saw for arrivals in a call center.

## Insight

This is the M/M/∞ queue in disguise. Detailed balance gives

$$
\pi_n \lambda_n = \pi_{n+1} \mu_{n+1} \implies \pi_n \lambda = \pi_{n+1} (n+1) \mu
$$

Solving the recurrence with $\rho = \lambda/\mu$:

$$
\pi_n = e^{-\rho} \frac{\rho^n}{n!}
$$

— the Poisson PMF.

## Example

Hiring rate $\lambda = 2$ per month, per-capita attrition $\mu = 0.025$.
Steady-state mean team size is $\lambda/\mu = 80$. The distribution is
$\text{Poisson}(80)$, so the standard deviation is $\sqrt{80} \approx 9$.
"We aim for 80" really means "there is a ~95% chance the real number is
between 62 and 98."

## Takeaway

The Poisson distribution is not a model for counts that happen to look
Poisson. It is the *consequence* of a simple generative mechanism: independent
arrivals at constant rate, independent decays at individual rate. Once you
recognize the mechanism, you get the distribution for free.
