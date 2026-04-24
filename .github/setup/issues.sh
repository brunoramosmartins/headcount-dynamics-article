#!/bin/bash
# Usage: bash .github/setup/issues.sh owner/repo
# Run AFTER labels.sh and milestones.sh

set -euo pipefail
REPO="${1:?Usage: bash issues.sh owner/repo}"

echo "Creating issues for $REPO..."

# ── PHASE 0 ──

gh issue create --repo "$REPO" \
  --title "[Phase 0] Write thesis and define scope" \
  --label "phase:0,type:documentation,priority:high" \
  --milestone "Phase 0 — Foundation" \
  --body "## Context
The thesis must argue that headcount is a stochastic process, not a deterministic plan.

## Tasks
- [ ] Draft central claim (v0.1)
- [ ] Define scope: DTMC, stationary distributions, absorption, CTMC, birth-death, applied scenarios
- [ ] Define anti-scope: queueing networks, semi-Markov, hidden Markov, MCMC, real company data
- [ ] Write 1-paragraph abstract

## Definition of Done
- [ ] \`docs/thesis.md\` with thesis, scope, anti-scope, audience, abstract

## References
- Project motivation (IT headcount budgeting)"

gh issue create --repo "$REPO" \
  --title "[Phase 0] Design headcount model and state space" \
  --label "phase:0,type:documentation,priority:high" \
  --milestone "Phase 0 — Foundation" \
  --body "## Context
The headcount model defines the state space, transition rates, and both
formulations (career levels DTMC + team size birth-death).

## Tasks
- [ ] Define Formulation A: career level states (Junior, Mid, Senior, Exit)
- [ ] Define Formulation B: team size birth-death (n = 0, 1, ..., N_max)
- [ ] Specify default parameters for both formulations
- [ ] Define 5 scenarios (growth, freeze, layoff, composition, seasonal)
- [ ] Document expansion points

## Definition of Done
- [ ] \`docs/model-design.md\` with both formulations and all parameters

## References
- Headcount Model Design section of roadmap"

gh issue create --repo "$REPO" \
  --title "[Phase 0] Configure repository, GitHub setup scripts, and Claude Code rules" \
  --label "phase:0,type:infrastructure,priority:high" \
  --milestone "Phase 0 — Foundation" \
  --body "## Context
Professional repo setup. GitHub configuration via bash scripts. Claude Code
rules enforce: no branches, no commits, no PRs, no co-author, no requirements.txt.

## Tasks
- [ ] Initialize all directories
- [ ] Create \`.claude/CLAUDE.md\`
- [ ] Create \`.github/ISSUE_TEMPLATE/task.md\` and \`bug.md\`
- [ ] Create \`.github/PULL_REQUEST_TEMPLATE.md\`
- [ ] Create \`.github/setup/labels.sh\`, \`milestones.sh\`, \`issues.sh\`
- [ ] Write \`pyproject.toml\` with deps + ruff config (no requirements.txt)
- [ ] Write initial \`README.md\`

## Definition of Done
- [ ] All scripts run successfully against the repo
- [ ] \`pip install -e .\` works
- [ ] Claude Code loads rules correctly

## References
- Repository structure and Claude Code config in roadmap"

# ── PHASE 1 ──

gh issue create --repo "$REPO" \
  --title "[Phase 1] Derive Markov chain foundations and Chapman-Kolmogorov" \
  --label "phase:1,type:theory,priority:critical" \
  --milestone "Phase 1 — Markov Foundations" \
  --body "## Context
The Markov property and Chapman-Kolmogorov equation are the two pillars.
State classification determines which tools apply later.

## Tasks
- [ ] Define Markov property formally
- [ ] Prove transition matrix is row-stochastic
- [ ] Derive n-step transitions: P^(n) = P^n
- [ ] Prove Chapman-Kolmogorov
- [ ] Classify states: transient, recurrent, absorbing
- [ ] Prove communication is equivalence relation
- [ ] 3-state worked example

## Definition of Done
- [ ] All proofs step-by-step in \`notes/phase1-markov-foundations.md\`
- [ ] State classification applied to headcount model

## References
- Norris, Markov Chains, Ch. 1
- Ross, Probability Models, Ch. 4"

gh issue create --repo "$REPO" \
  --title "[Phase 1] Implement DTMC module and path simulator" \
  --label "phase:1,type:code,priority:high" \
  --milestone "Phase 1 — Markov Foundations" \
  --body "## Context
The DTMC module wraps transition matrix operations and provides path simulation.

## Tasks
- [ ] Implement \`src/dtmc.py\`: MarkovChain class
- [ ] Implement \`src/simulation.py\`: path simulation
- [ ] Create tests (DO NOT RUN)

## Definition of Done
- [ ] Both modules implemented with type hints and docstrings
- [ ] Tests created
- [ ] \"Please run \`pytest tests/\` and \`ruff check .\`\"

## References
- Theory from Phase 1 theory issue"

# ── PHASE 2 ──

gh issue create --repo "$REPO" \
  --title "[Phase 2] Prove stationary distribution existence, uniqueness, and convergence" \
  --label "phase:2,type:theory,priority:critical" \
  --milestone "Phase 2 — Stationary Distributions" \
  --body "## Context
The stationary distribution tells the long-run team composition. Convergence
rate (spectral gap) tells how fast the team forgets its initial state.

## Tasks
- [ ] Define πP = π
- [ ] Prove existence (finite irreducible)
- [ ] Prove uniqueness (irreducible positive recurrent)
- [ ] State Perron-Frobenius
- [ ] Prove convergence: P^n → 1π^T
- [ ] Derive rate: ||P^n − 1π^T|| ≤ C|λ₂|^n
- [ ] Define spectral gap and mixing time

## Definition of Done
- [ ] Complete proof chain in \`notes/phase2-dtmc-analysis.md\`

## References
- Norris Ch. 1.7–1.8
- Levin, Peres & Wilmer, Mixing Times"

gh issue create --repo "$REPO" \
  --title "[Phase 2] Implement stationary computation and spectral analysis" \
  --label "phase:2,type:code,priority:high" \
  --milestone "Phase 2 — Stationary Distributions" \
  --body "## Context
Apply theory to headcount model: solve πP = π, compute spectral gap.

## Tasks
- [ ] Add \`.stationary()\`, \`.eigenvalues()\`, \`.spectral_gap()\`, \`.mixing_time()\` to MarkovChain
- [ ] Solve πP = π for headcount model
- [ ] Create convergence heatmap and eigenvalue spectrum
- [ ] Tests (DO NOT RUN)

## Definition of Done
- [ ] Analytical π matches numerical
- [ ] Figures generated

## References
- Theory from Phase 2 theory issue"

# ── PHASE 3 ──

gh issue create --repo "$REPO" \
  --title "[Phase 3] Derive fundamental matrix, absorption, and hitting times" \
  --label "phase:3,type:theory,priority:critical" \
  --milestone "Phase 3 — Absorption & Hitting Times" \
  --body "## Context
The fundamental matrix N = (I − Q)⁻¹ encodes absorption times and probabilities.

## Tasks
- [ ] Define absorbing chains and canonical form
- [ ] Prove I − Q is invertible
- [ ] Derive N = (I − Q)⁻¹
- [ ] Prove N_ij = expected visits to j from i
- [ ] Derive t = N·1 and B = NR
- [ ] Derive first passage times
- [ ] Prove mean return time = 1/π_i

## Definition of Done
- [ ] Complete derivation in \`notes/phase3-absorption-hitting.md\`

## References
- Grinstead & Snell Ch. 11, Norris Ch. 1.5"

gh issue create --repo "$REPO" \
  --title "[Phase 3] Implement absorption module and headcount scenarios" \
  --label "phase:3,type:code,priority:high" \
  --milestone "Phase 3 — Absorption & Hitting Times" \
  --body "## Context
Answers: how long until empty team (freeze), time Junior→Senior, P(Senior before Exit).

## Tasks
- [ ] Implement \`src/absorption.py\`
- [ ] Apply to 3 headcount scenarios
- [ ] Create tests (DO NOT RUN)

## Definition of Done
- [ ] All functions match hand-computed examples
- [ ] Three scenarios computed

## References
- Theory from Phase 3 theory issue"

# ── PHASE 4 ──

gh issue create --repo "$REPO" \
  --title "[Phase 4] Derive CTMC and birth-death process theory" \
  --label "phase:4,type:theory,priority:critical" \
  --milestone "Phase 4 — Continuous-Time" \
  --body "## Context
Birth-death processes model hiring/attrition as continuous flows.

## Tasks
- [ ] Define generator Q
- [ ] Derive Kolmogorov forward equation
- [ ] Prove P(t) = e^{Qt}
- [ ] Derive detailed balance for birth-death
- [ ] Prove Poisson steady-state for λ/nμ case

## Definition of Done
- [ ] All in \`notes/phase4-continuous-time.md\`

## References
- Norris Ch. 2, Ross Ch. 6"

gh issue create --repo "$REPO" \
  --title "[Phase 4] Implement CTMC module and birth-death headcount model" \
  --label "phase:4,type:code,priority:high" \
  --milestone "Phase 4 — Continuous-Time" \
  --body "## Context
Matrix exponential and birth-death for team size evolution.

## Tasks
- [ ] Implement \`src/ctmc.py\`: CTMC and BirthDeathProcess classes
- [ ] Apply to headcount: λ=2, μ=0.025/person
- [ ] Create tests (DO NOT RUN)

## Definition of Done
- [ ] Matrix exponential matches scipy
- [ ] Birth-death steady state matches Poisson

## References
- Theory from Phase 4 theory issue"

# ── PHASE 5 ──

gh issue create --repo "$REPO" \
  --title "[Phase 5] Build applied headcount model and run 5 scenarios" \
  --label "phase:5,type:code,type:experiment,priority:critical" \
  --milestone "Phase 5 — Applied Headcount Model" \
  --body "## Context
Full headcount model combining DTMC and birth-death. Five scenarios.

## Tasks
- [ ] Implement \`src/headcount.py\`
- [ ] Run scenarios S1–S5
- [ ] Create figures
- [ ] Connect to budget (link to Articles 1 & 2)
- [ ] Tests (DO NOT RUN)

## Definition of Done
- [ ] All 5 scenarios analysed with figures
- [ ] Budget connection explicit

## References
- All theory from Phases 1–4"

gh issue create --repo "$REPO" \
  --title "[Phase 5] Sensitivity analysis on transition rates" \
  --label "phase:5,type:experiment,priority:high" \
  --milestone "Phase 5 — Applied Headcount Model" \
  --body "## Context
Which transition rates matter most for steady-state team size?

## Tasks
- [ ] Vary each rate ±50%
- [ ] Compute impact on steady-state
- [ ] Create tornado chart

## Definition of Done
- [ ] Tornado chart with business interpretation

## References
- Default parameters from model design"

# ── PHASE 6 ──

gh issue create --repo "$REPO" \
  --title "[Phase 6] Run all experiments and create publication figures" \
  --label "phase:6,type:experiment,priority:critical" \
  --milestone "Phase 6 — Experiments & Visualizations" \
  --body "## Context
All 7 experiments with publication-quality figures.

## Tasks
- [ ] Experiments A–G (see roadmap for details)
- [ ] 300 DPI PNGs, fixed seeds
- [ ] Animated trajectory GIF

## Definition of Done
- [ ] All figures in \`figures/\`
- [ ] GIF < 5 MB

## References
- All theory from Phases 1–5"

# ── PHASE 7 ──

gh issue create --repo "$REPO" \
  --title "[Phase 7] Write article sections 1–5 (English)" \
  --label "phase:7,type:writing,priority:high" \
  --milestone "Phase 7 — Article Writing" \
  --body "## Context
First half: Markov chains through birth-death processes.

## Tasks
- [ ] Write sections 1–5 in \`article/headcount-dynamics.md\`
- [ ] Self-contained derivations, consistent notation

## Definition of Done
- [ ] Sections 1–5 complete in English

## References
- Theory notes Phases 1–4"

gh issue create --repo "$REPO" \
  --title "[Phase 7] Write article sections 6–10 (English)" \
  --label "phase:7,type:writing,priority:high" \
  --milestone "Phase 7 — Article Writing" \
  --body "## Context
Second half: applied model, experiments, budget bridge, conclusion.

## Tasks
- [ ] Write sections 6–10
- [ ] Budget connection referencing Articles 1 & 2
- [ ] Practical framework

## Definition of Done
- [ ] Full article complete in English

## References
- Phase 5 model, Phase 6 figures"

gh issue create --repo "$REPO" \
  --title "[Phase 7] Translate article to Portuguese" \
  --label "phase:7,type:writing,priority:medium" \
  --milestone "Phase 7 — Article Writing" \
  --body "## Context
Portuguese version for local review. Translation, not rewrite.

## Tasks
- [ ] Translate \`headcount-dynamics.md\` → \`headcount-dynamics-ptbr.md\`
- [ ] Keep LaTeX and figure references unchanged

## Definition of Done
- [ ] Full Portuguese translation

## References
- English article"

# ── PHASE 8 ──

gh issue create --repo "$REPO" \
  --title "[Phase 8] Mathematical validation and code reproducibility" \
  --label "phase:8,type:review,priority:critical" \
  --milestone "Phase 8 — Review & Publish" \
  --body "## Context
Final quality gate.

## Tasks
- [ ] Review all derivations
- [ ] Author runs: \`pip install -e .\` → scripts → tests → ruff
- [ ] Fix discrepancies

## Definition of Done
- [ ] Zero errors, all tests pass, ruff clean

## References
- Full article"

gh issue create --repo "$REPO" \
  --title "[Phase 8] Publish to GitHub Pages and Medium" \
  --label "phase:8,type:writing,type:infrastructure,priority:high" \
  --milestone "Phase 8 — Review & Publish" \
  --body "## Context
Publication and distribution.

## Tasks
- [ ] Copy English MD to github.io repo, run pipeline, verify
- [ ] Medium cross-post with canonical link
- [ ] LinkedIn post
- [ ] Update README

## Definition of Done
- [ ] Article live, Medium published, LinkedIn drafted

## References
- Author's github.io pipeline"

echo "All issues created."
echo ""
echo "Verify milestone IDs:"
echo "  gh api repos/$REPO/milestones --jq '.[] | {number, title}'"
