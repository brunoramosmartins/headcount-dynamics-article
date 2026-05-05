# TIL — Today I Learned

Short, portfolio-ready notes capturing one non-obvious insight from each phase.
Each TIL is a **skeleton**: hook → insight → example → takeaway. The author
writes the final text in their own voice.

## Why TILs?

The article publishes once. The TILs publish continuously — one per phase on
LinkedIn, the personal blog, or Medium. They keep the work visible while the
article matures.

## Format

Each file:

- Title that states the insight, not the topic
- 100–300 words
- One concrete example (numbers or a small matrix)
- One-line takeaway
- Tags: the phase, the math branch, the applied domain

## Index

Two TILs are tied to Phase 1. The "plain-words" entry is a beginner-friendly
overview of the Markov chain idea — useful as a first post for an audience
unfamiliar with the topic. The "row-stochastic" entry digs into a single
non-obvious technical detail. Pick whichever fits the audience of the
moment, or run them as a two-post series.

| Phase | TIL file                                               | Insight                                       |
|-------|--------------------------------------------------------|-----------------------------------------------|
| 1     | [til-phase-1-markov-chain-in-plain-words.md](til-phase-1-markov-chain-in-plain-words.md) | What a Markov chain *is*, in plain language  |
| 1     | [til-phase-1-why-row-stochastic.md](til-phase-1-why-row-stochastic.md) | Why rows of $P$ sum to 1 (not columns)        |
| 2     | [til-phase-2-spectral-gap-intuition.md](til-phase-2-spectral-gap-intuition.md) | Spectral gap = mixing speed                   |
| 3     | [til-phase-3-fundamental-matrix.md](til-phase-3-fundamental-matrix.md) | $N = (I-Q)^{-1}$ counts visits, not steps     |
| 4     | [til-phase-4-poisson-from-birth-death.md](til-phase-4-poisson-from-birth-death.md) | A team with constant hiring is Poisson-sized  |
| 5     | [til-phase-5-headcount-as-process.md](til-phase-5-headcount-as-process.md) | A plan is a point estimate of a distribution  |
| 6     | [til-phase-6-fixed-seeds-repro.md](til-phase-6-fixed-seeds-repro.md) | Why every figure gets a seed                  |

## Publishing order

Draft in this folder during the phase, publish after the phase PR is merged.
Order on LinkedIn matches the index. Each post links to the GitHub folder so
the audience can follow progress.
