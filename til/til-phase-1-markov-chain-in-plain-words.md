# TIL — A Markov Chain Is a Table of "What Happens Next"

**Phase:** 1 · **Topic:** Markov chains, intuition · **Domain:** workforce modelling

## Hook

You don't need a doctorate to grasp the core idea. A Markov chain is a way
of describing a system that keeps changing state, where the only thing
needed to predict the next move is *the current state* — not how you got
there.

## Insight

Think of an employee in a software team. At any month they are in one of a
small number of buckets — Junior, Mid, Senior, or *gone*. Next month, they
either stay, move up, or leave. If you know the *probabilities* of each
move from each bucket, you have all you need to predict what happens.
Those probabilities live in a table — one row per starting bucket, one
column per destination bucket — and the rows sum to 1 because *something*
has to happen.

That's it. That table is the entire model. No history of past job moves,
no individual variation, no calendar effects — just "if you are here now,
where do you end up next month?". Everything else in the article — the
long-run team composition, the probability of hitting a target headcount,
the expected time until everyone has left under a hiring freeze — is a
calculation on that one table.

## Example

For a junior developer, the next-month transition probabilities might be:
93% stay Junior, 3% promoted to Mid, 4% leave the company. Three numbers.
Same for Mid (96% stay, 2% promote, 2% leave) and Senior (99% stay, 1%
leave). That entire team can now be analysed without any further input.

## Takeaway

If a system can be described by "what state am I in, and what are the
probabilities of moving to each other state next?", it is a Markov chain
— and you can compute its long-run behaviour, its hitting times, and its
absorption probabilities with linear algebra. The hard part is naming the
states and estimating the probabilities. After that, it is bookkeeping.
