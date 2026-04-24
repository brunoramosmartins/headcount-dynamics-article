# Exercise Set 1 — Markov Chain Foundations

> Work everything by hand first. Only check with code after you have a full
> solution on paper. Each exercise notes the theorem or technique it exercises.

## Proofs (paper)

**1.** Prove that if $P$ is a row-stochastic matrix, then $P^n$ is also
row-stochastic for all $n \geq 1$.
*Technique:* induction on $n$. Show that $(P^{n+1}\mathbf{1}) = P^n (P\mathbf{1}) = P^n \mathbf{1} = \mathbf{1}$.

**2.** Prove the Chapman-Kolmogorov equation
$$p_{ij}^{(m+n)} = \sum_{k} p_{ik}^{(m)} p_{kj}^{(n)}$$
directly from the law of total probability and the Markov property.
*Technique:* condition on $X_m$, apply the Markov property to split the two halves.

**3.** Prove that "communicates with" is an equivalence relation on the state
space: reflexive ($i \leftrightarrow i$), symmetric
($i \leftrightarrow j \Rightarrow j \leftrightarrow i$), and transitive.

**4.** Prove that all states in the same communication class are of the same
type (all transient or all recurrent).
*Hint:* show that if $i$ is recurrent and $i \leftrightarrow j$, then $j$ is
recurrent. Use the fact that recurrence ↔ $\sum_n p_{ii}^{(n)} = \infty$.

---

## Computations (paper)

**5.** Consider the headcount model with 3 career levels + Exit:

$$
P = \begin{pmatrix}
0.93 & 0.03 & 0    & 0.04 \\
0    & 0.96 & 0.02 & 0.02 \\
0    & 0    & 0.99 & 0.01 \\
0.50 & 0    & 0    & 0.50
\end{pmatrix}
$$

Compute $P^2$ by hand. What is the probability that a Junior employee is still
Junior after 2 months? What is the probability they have been promoted to Mid?

**6.** For the same $P$, classify each state as transient or recurrent. Is the
chain irreducible? Identify the communication classes.

**7.** Starting with the initial distribution $\pi_0 = (1, 0, 0, 0)$ (all
Juniors), compute $\pi_1 = \pi_0 P$ and $\pi_2 = \pi_1 P$ (or $\pi_0 P^2$).
Verify that $\pi_2$ matches your $P^2$ computation from exercise 5.

---

## Self-check

- [ ] Exercise 1: proof uses induction and the identity $P\mathbf{1} = \mathbf{1}$
- [ ] Exercise 2: the summation index is the intermediate state at step $m$
- [ ] Exercises 5–7: all row sums of $P^2$ equal 1 (sanity check)
- [ ] Exercise 7: $\pi_2$ from two one-step iterations equals the first row of $P^2$
