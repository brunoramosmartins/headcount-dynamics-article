# TIL — $N = (I - Q)^{-1}$ Counts Visits, Not Steps

**Phase:** 3 · **Topic:** absorbing chains · **Domain:** workforce modelling

## Hook

The fundamental matrix looks like a formal trick from linear algebra. It is
actually a **counter**: entry $N_{ij}$ is the expected number of times the
chain visits state $j$ before being absorbed, starting from $i$.

## Insight

$$
N = (I - Q)^{-1} = \sum_{k=0}^{\infty} Q^k
$$

Each $(Q^k)_{ij}$ is the probability that the chain is in $j$ at step $k$
*without yet being absorbed*. Sum over all $k$: expected number of visits.
That single matrix then gives you almost everything:

- Expected absorption time: $t = N \mathbf{1}$ (sum visits across all $j$)
- Absorption probabilities: $B = N R$
- Variance of absorption time: $(2N - I)t - t \odot t$

## Example

Headcount model with Exit absorbing. The transient submatrix $Q$ is
$3 \times 3$ (Junior, Mid, Senior). Inverting $I - Q$ gives expected time-to-
exit per starting level — for the default parameters, Junior expects $\sim 25$
months, Mid $\sim 50$ months, Senior $\sim 100$ months. Senior retention is
worth more than double mid retention.

## Takeaway

Every time you invert $I - Q$ for an absorbing chain, you are secretly
summing an infinite series of "probability of still being alive at step $k$".
The matrix inverse is just the closed form.
