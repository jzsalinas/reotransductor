# Reotransductor — Agent Operating Rules

## Purpose

Reotransductor is a scientific hypothesis-testing project.

The objective of the project is NOT to confirm the Reotransductor hypothesis.
The objective is to determine rigorously whether the hypothesis is mathematically consistent, physically meaningful, numerically reproducible, falsifiable, and compatible or incompatible with evidence.

A negative result is a valid and scientifically valuable outcome.

## Communication

- Communicate with the project owner in Spanish unless explicitly requested otherwise.
- Write source-code comments, docstrings, identifiers, and technical code documentation in English.
- Use precise, neutral, scientific language.
- Avoid rhetorical claims, promotional language, and decorative icons in source code or scientific technical documentation.

## Scientific integrity

Never optimize implementation, equations, constants, initial conditions, numerical cutoffs, or interpretation merely to obtain a desired observational result.

Never transform:
- a hypothesis into an established fact;
- a numerical result into an observational validation;
- a phenomenological calibration into a first-principles derivation;
- implementation success into scientific confirmation;
- correlation into causation.

Synthetic data, bootstrap fields, phenomenological constants, numerical regulators, sub-grid models, fitted parameters, and approximations must always be explicitly identified as such.

Never invent mathematical or physical justification for an existing implementation.

When physical justification is uncertain, report the uncertainty instead of constructing a justification.

## Model governance

Before modifying project-specific physics, equations, constants, cosmological assumptions, or phenomenology, consult:

1. GENESIS_RECONCILIATION.md
2. FOUNDATIONAL_HYPOTHESIS.md
3. SCIENTIFIC_CONTRACT.md
4. MODEL_SPECIFICATION.md

If any of these files is missing or does not define the relevant scientific decision, do NOT make an autonomous physics change.

Instead, report:

SCIENTIFIC DECISION REQUIRED

and explain the unresolved choice.

Established physics, physical constants, observational values, datasets, and literature claims must be verified against appropriate primary sources or provenance records when relevant.

## Paper / code consistency

Paper I describes the theoretical hypothesis and its explicitly declared extensions.

Paper II must describe the numerical model actually implemented.

If a manuscript and the executable implementation disagree, do NOT silently change either one.

Report the discrepancy and classify it as one of:

- IMPLEMENTATION BUG
- DOCUMENTATION MISMATCH
- NUMERICAL APPROXIMATION
- PHENOMENOLOGICAL ASSUMPTION
- UNRESOLVED SCIENTIFIC CHOICE

Paper III and observational results must never be used retroactively to tune the foundational theory or numerical equations merely to improve agreement with observations.

## Code changes

For every scientifically relevant code change:

- identify the governing equation or model assumption affected;
- identify whether the change is physical, phenomenological, or purely numerical;
- perform dimensional/unit consistency checks where applicable;
- add or update tests;
- report expected scientific consequences;
- report whether previous simulation results may have become invalid or non-comparable.

Do not introduce hidden clamps, floors, ceilings, boosts, normalization constants, or empirical factors.

## Tests

Tests must exercise production configuration whenever they are intended to validate production physics.

A test that only reproduces current implementation behavior is not sufficient evidence of physical correctness.

Numerical convergence tests must compare physically equivalent states, not merely equal iteration counts, unless equal iteration count is itself the intended benchmark.

## Git safety

Never run git add, git commit, git push, git reset, git checkout, git clean, or destructive Git operations unless explicitly requested by the project owner.

Never stage unrelated user changes.

Never use `git add .` by default.

When changes are requested, show the relevant diff/status before any commit.

## Repository safety

Do not delete checkpoints, datasets, results, logs, local configuration, credentials, or generated scientific artifacts unless explicitly requested.

Do not modify running-production checkpoints in place.

Preserve reproducibility and provenance.

## Default behavior

When asked to audit:
- inspect first;
- modify nothing unless explicitly requested;
- distinguish code behavior from manuscript claims;
- distinguish established physics from Reotransductor hypotheses;
- prefer reporting a contradiction over resolving it through assumption.

When asked to implement:
- implement only the explicitly approved scientific decision;
- do not broaden scope;
- run relevant tests;
- report exactly what changed.
