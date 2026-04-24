#!/bin/bash
# Usage: bash .github/setup/milestones.sh owner/repo

set -euo pipefail
REPO="${1:?Usage: bash milestones.sh owner/repo}"

echo "Creating milestones for $REPO..."

gh api "repos/$REPO/milestones" -f title="Phase 0 — Foundation" \
  -f description="Thesis, model design, project scaffold, GitHub configuration." --silent
gh api "repos/$REPO/milestones" -f title="Phase 1 — Markov Foundations" \
  -f description="Markov chain definition, Chapman-Kolmogorov, state classification, DTMC module." --silent
gh api "repos/$REPO/milestones" -f title="Phase 2 — Stationary Distributions" \
  -f description="Stationary distribution existence/uniqueness, convergence, spectral analysis." --silent
gh api "repos/$REPO/milestones" -f title="Phase 3 — Absorption & Hitting Times" \
  -f description="Fundamental matrix, absorption probabilities, first passage times." --silent
gh api "repos/$REPO/milestones" -f title="Phase 4 — Continuous-Time" \
  -f description="CTMC, Kolmogorov equations, birth-death processes." --silent
gh api "repos/$REPO/milestones" -f title="Phase 5 — Applied Headcount Model" \
  -f description="Full headcount model, 5 scenarios, sensitivity analysis." --silent
gh api "repos/$REPO/milestones" -f title="Phase 6 — Experiments & Visualizations" \
  -f description="All experiments, publication figures, animated trajectories." --silent
gh api "repos/$REPO/milestones" -f title="Phase 7 — Article Writing" \
  -f description="Full article in English and Portuguese." --silent
gh api "repos/$REPO/milestones" -f title="Phase 8 — Review & Publish" \
  -f description="Validation, reproducibility, publication." --silent

echo "All milestones created."
echo ""
echo "Verify milestone IDs:"
echo "  gh api repos/$REPO/milestones --jq '.[] | {number, title}'"
