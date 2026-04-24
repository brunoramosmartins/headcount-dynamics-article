# TIL — A Transition Matrix Sums to 1 by Row, Not by Column (Here Is Why)

**Phase:** 1 · **Topic:** Markov chains · **Domain:** workforce modelling

## Hook

Most people meet matrices through systems of equations, where we think
column-by-column (coefficients of each variable). Transition matrices flip
this: the meaningful direction is the **row**.

## Insight

Row $i$ of $P$ is the conditional distribution
$P(X_{n+1} = \cdot \mid X_n = i)$. A probability distribution must sum to 1,
so the row sums to 1. Columns have no such constraint — a column sums to the
probability that *anything* transitions into state $j$, which depends on where
mass currently lives, not on the laws of probability.

This is why the update is $\pi_{n+1} = \pi_n P$ (a row vector multiplied from
the **left**), not $P \pi$. The distribution flows through the rows.

## Example

In the headcount model, the Junior row is `(0.93, 0.03, 0, 0.04)` — a Junior
stays Junior with 93% probability, gets promoted with 3%, cannot jump straight
to Senior, exits with 4%. The row sums to 1. The Exit column sums to
$0.04 + 0.02 + 0.01 + 0.50 = 0.57$ — meaningless in isolation.

## Takeaway

If your transition matrix has columns that sum to 1, you either (a) wrote the
left-stochastic convention (common in some physics textbooks) or (b) made a
mistake. Pick a convention, write it down, and multiply accordingly.
