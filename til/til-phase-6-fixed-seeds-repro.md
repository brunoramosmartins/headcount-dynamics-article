# TIL — Every Figure Deserves a Seed

**Phase:** 6 · **Topic:** experiment reproducibility · **Domain:** scientific communication

## Hook

Run a Monte Carlo experiment twice and get two slightly different figures. Run
it three years later and you cannot reproduce either. The fix is one line.

## Insight

At the top of every experiment script:

```python
SEED = 42
rng = numpy.random.default_rng(SEED)
```

Then pass `rng` into every sampler. Three habits follow:

1. The figure is **deterministic** given the code — anyone who re-runs the
   script gets the same picture.
2. The random-number generator is **explicit** — no hidden global state.
3. The seed sits in the **script**, not in a config file or environment
   variable — one fewer thing to lose.

## Example

For the animated trajectory GIF in the headcount article, 100 simulated paths
of 24 months each are sampled. With `default_rng(42)`, the first path is
always the same `[45, 46, 47, 46, 48, ...]`. Reviewers can flag a specific
trajectory by index; I can reproduce it without guessing.

## Takeaway

The computation cost of setting a seed is zero. The communication cost of
**not** setting one — a reviewer saying "I cannot reproduce your figure 4" —
can be hours. Always seed.
