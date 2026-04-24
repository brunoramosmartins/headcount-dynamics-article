# Headcount Dynamics — The Team That Replaces Itself

**Headcount and turnover as stochastic processes — from Markov chains to birth-death models.**

A portfolio-grade technical article that models workforce dynamics — hiring,
attrition, promotions, and steady-state composition — as stochastic processes.
The article starts with discrete-time Markov chains (month-to-month
transitions) and extends to continuous-time birth-death processes (differential
equations for team size evolution).

**Key insight:** headcount is not a number you plan — it is a random process
with a stationary distribution, absorption probabilities, and expected hitting
times that can be computed analytically.

---

## The Trilogy

| Aspect            | Article 1 (Monte Carlo)     | Article 2 (Distributions) | Article 3 (This one) |
|-------------------|------------------------------|---------------------------|----------------------|
| Core question     | How to simulate total budget | Which distributions model each cost | How headcount evolves over time |
| Math branch       | Simulation theory            | Statistical inference     | Stochastic processes |
| Time dimension    | Static                       | Static                    | **Dynamic**          |

---

## Repository Structure

```
headcount-dynamics-article/
├── .claude/CLAUDE.md                 # Claude Code project rules
├── .github/                          # Templates and setup scripts
│   ├── ISSUE_TEMPLATE/
│   ├── PULL_REQUEST_TEMPLATE.md
│   └── setup/                        # bash scripts for labels, milestones, issues
├── article/                          # Final article (English + Portuguese)
├── docs/                             # Thesis, model design, outline
├── src/                              # Core modules (dtmc, ctmc, absorption, ...)
├── scripts/                          # Experiment scripts
├── notebooks/                        # Jupyter notebooks
├── exercises/                        # Paper exercises (one file per phase)
├── til/                              # Today-I-Learned posts (portfolio)
├── notes/                            # Theory notes per phase
├── tests/                            # Unit tests
├── figures/                          # Generated figures
├── pyproject.toml
└── LICENSE
```

---

## Setup

```bash
# Clone and install in editable mode
pip install -e ".[dev]"
```

---

## GitHub Bootstrap (one-time, by author)

After pushing the first commit and creating the repo on GitHub:

```bash
# 1. Labels
bash .github/setup/labels.sh brunoramosmartins/headcount-dynamics-article

# 2. Milestones (one per phase)
bash .github/setup/milestones.sh brunoramosmartins/headcount-dynamics-article

# 3. Issues (all phases)
bash .github/setup/issues.sh brunoramosmartins/headcount-dynamics-article
```

---

## Workflow

Each phase:

1. **Paper exercises** (`exercises/exNN_*.md`) — proofs and numerical computations by hand.
2. **Theory notes** (`notes/phaseN-*.md`) — rigorous derivations.
3. **Implementation** (`src/*.py`) — with type hints, docstrings, and tests.
4. **Experiments** (`scripts/exp_*.py`) — produce figures in `figures/`.
5. **TIL post** (`til/til-phaseN-*.md`) — portfolio-ready takeaway.
6. **Commit, PR, tag** — following conventional commits and the roadmap.

Author runs all validation locally:

```bash
pytest tests/
ruff check .
```

---

## Development Commands

| Task                  | Command                                     |
|-----------------------|---------------------------------------------|
| Install in dev mode   | `pip install -e ".[dev]"`                   |
| Run tests             | `pytest tests/`                             |
| Lint                  | `ruff check .`                              |
| Auto-fix              | `ruff check . --fix`                        |
| Run an experiment     | `python scripts/exp_transition_convergence.py` |

---

## Roadmap

See [`roadmap-headcount-dynamics-v1.md`](./roadmap-headcount-dynamics-v1.md) for
the full phase-by-phase plan.

| Phase | Focus                                    | Tag                       |
|-------|------------------------------------------|---------------------------|
| 0     | Foundation                               | `v0.1-foundation`         |
| 1     | Markov chain foundations                 | `v0.2-markov-foundations` |
| 2     | Stationary distributions & convergence   | `v0.3-dtmc-analysis`      |
| 3     | Absorption & hitting times               | `v0.4-absorption-hitting` |
| 4     | Continuous-time, birth-death             | `v0.5-continuous-time`    |
| 5     | Applied headcount model                  | `v0.6-headcount-model`    |
| 6     | Experiments & visualizations             | `v0.7-experiments`        |
| 7     | Article writing                          | `v0.8-article-draft`      |
| 8     | Review & publish                         | `v1.0.0`                  |

---

## License

MIT — see [LICENSE](./LICENSE).
