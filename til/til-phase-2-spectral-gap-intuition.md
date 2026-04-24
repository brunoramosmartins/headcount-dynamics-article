# TIL — The Spectral Gap Is the Speed of Forgetting

**Phase:** 2 · **Topic:** convergence & mixing · **Domain:** workforce modelling

## Hook

Every Markov chain converges to its stationary distribution. Less obvious: the
**rate** of convergence is one number — the spectral gap.

## Insight

If $\lambda_1 = 1 > |\lambda_2| \geq |\lambda_3| \geq \ldots$ are the
eigenvalues of $P$ in modulus, then

$$
\|P^n - \mathbf{1}\pi^T\| \leq C \cdot |\lambda_2|^n
$$

The **spectral gap** $1 - |\lambda_2|$ controls how fast the chain mixes. A
chain with $\lambda_2 = 0.99$ takes $\sim 100$ steps to forget its start; one
with $\lambda_2 = 0.5$ takes $\sim 2$.

## Example

For the headcount transition matrix with default rates, numerical
eigendecomposition gives roughly $\lambda_2 \approx 0.98$. Spectral gap
$\approx 0.02$. Mixing time $\sim 50$ months — more than four years before the
team "forgets" its initial composition. That is long enough that short-term
hiring changes can reshape the near-term distribution without the stationary
distribution ever being reached.

## Takeaway

"Long-run" is a verb, not a noun. The spectral gap tells you **how long** the
long-run actually takes. If it is larger than your planning horizon, the
stationary distribution is a useful compass but not a destination.
