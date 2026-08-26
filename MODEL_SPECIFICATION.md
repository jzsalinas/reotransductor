# Reotransductor — Model Specification

## 1. Status and purpose

This document defines the scientific model hierarchy of Reotransductor.

It does NOT declare the current cosmological engine to be physically correct.

It does NOT automatically approve a feature merely because that feature
already exists in source code, manuscripts, historical documentation, or
previous simulation results.

The purpose of this specification is to distinguish:

1. established physical baseline;
2. foundational Reotransductor candidates under review;
3. approved model extensions;
4. provisional numerical choices;
5. experimental extensions;
6. legacy or unresolved assumptions.

Until an item is explicitly promoted to APPROVED, an AI agent must not
treat its current implementation as scientific authority.

---

## 2. Specification status

Current specification version:

`MODEL-SPEC-0.2`

Current overall model status:

**FOUNDATIONAL RECONCILIATION IN PROGRESS —
COSMOLOGICAL MODEL NOT FROZEN**

This means that the historical and later formulations of the Reotransductor
clock are being reconciled. The full cosmological implementation also contains
scientific and numerical choices requiring independent review.

Paper III and observational claims are outside the scope of the present
model-freezing phase.

---

## 3. Status vocabulary

Every model component must have one of the following statuses.

### ESTABLISHED_BASELINE

Accepted external physics used as a reference framework.

This does not mean that every numerical approximation used to represent it
is automatically valid.

### FROZEN_FOUNDATIONAL

Part of the minimal Reotransductor hypothesis defined in
`FOUNDATIONAL_HYPOTHESIS.md`.

Changing it constitutes a new foundational model version.

### APPROVED_EXTENSION

A project-specific extension that has received explicit scientific approval
after its assumptions, equations, units, consequences, and falsifiability
were reviewed.

### PROVISIONAL_NUMERICAL

A numerical mechanism temporarily permitted for stability, regularization,
finite precision, or computational practicality.

It must not be interpreted as physical law.

### EXPERIMENTAL_EXTENSION

A hypothesis being explored but not required by the foundational model.

### UNDER_REVIEW

An existing or proposed component whose scientific formulation has not yet
been accepted.

### LEGACY

A historical implementation or interpretation that may remain in the
repository but is not part of the approved scientific model.

### NOT_IMPLEMENTED

A defined alternative or future branch that is not currently executable
physics.

---

## 4. Layer 0 — Established physical baseline

The project accepts the following only within the domains where established
physics supports them.

### B0.1 Geometric proper time

General Relativity provides geometric proper time:

\[
d\tau_{\mathrm{geom}}
=
\frac{1}{c}
\sqrt{-g_{\mu\nu}dx^\mu dx^\nu}.
\]

Status:

`ESTABLISHED_BASELINE`

Reotransductor does not replace this relation at the foundational level.

### B0.2 Local entropy production

Relativistic non-equilibrium thermodynamics permits description through an
entropy current \(s^\mu\) satisfying, for physically admissible irreversible
processes,

\[
\Sigma = \nabla_\mu s^\mu \ge 0.
\]

Status:

`ESTABLISHED_BASELINE`

The inequality is established baseline.

A particular Reotransductor formula chosen to approximate \(\Sigma\) is NOT
automatically established physics.

### B0.3 Standard constants

Quantities such as

\[
c,\quad G,\quad \hbar,\quad k_B
\]

must use traceable accepted values when represented as fundamental
constants.

Status:

`ESTABLISHED_BASELINE`

Their presence in an equation does not make that equation first-principles.

---

## 5. Layer 1 — Foundational Reotransductor candidates under review

No project-specific clock interpretation is currently frozen. The additive
formulation recorded below is candidate FH-001B and remains under review
pending resolution of FH-001.

### C1.1 Candidate thermodynamic clock

Introduce a candidate scalar field

\[
\tau_R(x^\mu).
\]

Status:

`UNDER_REVIEW`

### C1.2 Candidate additive dissipative-clock hypothesis (FH-001B)

The foundational postulate is

\[
u^\mu\nabla_\mu\tau_R
=
1+\kappa\Sigma.
\]

Equivalently,

\[
\frac{d\Delta\tau_{\mathrm{diss}}}{dt}
=
\kappa\Sigma.
\]

Status:

`UNDER_REVIEW`

This equation is a hypothesis.

It is not derived from the cited intellectual influences.

### C1.3 Null model

The model must always permit

\[
\kappa=0.
\]

Status:

`UNDER_REVIEW`

### C1.4 No foundational gravitational backreaction

The foundational hypothesis does NOT require \(\tau_R\) or
\(\Delta\tau_{\mathrm{diss}}\) to gravitate.

Any gravitational backreaction of the thermodynamic clock is therefore an
extension and requires independent justification.

Status:

`UNDER_REVIEW`

### C1.5 No foundational cosmological-cycle requirement

The foundational hypothesis does NOT require:

- Conformal Cyclic Cosmology;
- cosmological eons;
- holographic memory;
- black-hole endpoint physics;
- dark-matter replacement;
- cosmological structure formation.

Status:

`UNDER_REVIEW`

---

## 6. Layer 2 — Numerical hypothesis-testing framework

A numerical PDE framework may be used to investigate consequences of the
foundational hypothesis.

However, the current production engine is not yet declared an approved
scientific realization of the full cosmological model.

Status:

`UNDER_REVIEW`

The engine must eventually satisfy:

1. equations in code match the approved specification;
2. code units map consistently to physical units;
3. regulators are explicit;
4. production parameters are tested;
5. convergence is demonstrated at physically equivalent states;
6. checkpoints preserve complete scientific configuration;
7. Paper II documents the executable model exactly.

---

## 7. Extension registry

### E1 — Passive intra-eon dissipative tracer

Description:

During ordinary evolution, newly accumulated
\(\Delta\tau_{\mathrm{diss}}\) may be treated as a passive diagnostic field
without immediate gravitational backreaction.

Status:

`UNDER_REVIEW`

This architecture is compatible with the foundational model but has not yet
been frozen as the unique physical interpretation.

---

### E2 — Gravitational coupling of inherited temporal memory

Description:

A fossil temporal-memory field from a previous cosmological cycle is used
as an effective source in the gravitational potential.

Representative form:

\[
\rho_{\mathrm{eff}}
=
\rho+\mathcal{H}\tau_{\mathrm{prior}}.
\]

Status:

`EXPERIMENTAL_EXTENSION`

The value, dimensions, interpretation, and conservation implications of
\(\mathcal{H}\) require independent derivation or explicit phenomenological
classification.

This mechanism is NOT part of the foundational hypothesis.

---

### E3 — Holographic 1:1 gravitational parity

Description:

The specific choice

\[
\mathcal{H}=1
\]

interpreting inherited temporal memory as gravitationally equivalent to
matter.

Status:

`UNDER_REVIEW`

No first-principles status is currently assigned.

---

### E4 — Synthetic Eon-0 bootstrap

Description:

A synthetic initial memory field correlated with primordial density is used
to initialize the first numerically represented eon.

Status:

`UNDER_REVIEW`

Classification:

`INITIALIZATION / PHENOMENOLOGICAL`

The existence of an infinite cyclic model does not by itself determine the
numerical bootstrap amplitude.

The bootstrap must not be used as evidence for the foundational hypothesis.

---

### E5 — Inter-eon Fourier phase memory

Description:

Primordial fluctuations of a later eon may inherit Fourier phase
information from the predecessor dissipative-memory field.

Status:

`EXPERIMENTAL_EXTENSION`

Any memory coefficient such as \(\alpha_{\mathrm{mem}}\) is phenomenological
unless independently derived.

The exact circular-phase interpolation used by the implementation must be
documented if this extension is retained.

---

### E6 — Conformal Cyclic Cosmology embedding

Description:

Reotransductor may be embedded experimentally within a Penrose-inspired
Conformal Cyclic Cosmology architecture.

Status:

`EXPERIMENTAL_EXTENSION`

CCC is not implied by the foundational Reotransductor hypothesis.

Likewise, Reotransductor is not implied by CCC.

---

### E7 — Finite numerical conformal cutoff

Description:

A finite expansion-coordinate threshold may be used as a computational
proxy for an asymptotic cosmological boundary.

Status:

`PROVISIONAL_NUMERICAL`

A finite cutoff such as a value historically represented by `7.0` must not
be described as the physical CCC boundary itself.

Its effect must be studied through cutoff sensitivity.

---

### E8 — Stable black-hole-era memory saturation

Description:

During a resolved astrophysical era, the dissipative excess may be bounded
numerically to prevent finite-precision runaway while compact objects are
treated as effectively stable over the simulated interval.

Status:

`PROVISIONAL_NUMERICAL`

The numerical saturation of \(\Delta\tau_{\mathrm{diss}}\) is NOT equivalent
to physical cessation of proper time.

A fixed ceiling such as \(10^7\) has no fundamental status unless separately
derived.

---

### E9 — Hawking evaporation endpoint

Description:

A remote-future branch in which astrophysical black holes lose mass through
semiclassical Hawking evaporation before a possible conformal crossover.

Status:

`NOT_IMPLEMENTED`

This branch must not be claimed as part of executable results until mass,
energy, radiation, and timescale evolution are implemented consistently.

---

### E10 — Quantum white-hole endpoint

Description:

A local black-to-white-hole transition or Planck-star-inspired quantum
gravity branch.

Status:

`NOT_IMPLEMENTED`

If ever implemented, it must remain a local physical event and must not
automatically trigger a global cosmological reset unless separately
justified.

---

## 8. Foundational decision register

The recovered historical genesis requires four foundational decisions before
cosmological model freezing can resume. No AI agent may resolve these choices
autonomously.

### FH-001 — Meaning of the Reotransductor clock

Alternatives:

- **A — Dissipative progress variable:** a scalar measuring accumulated
  irreversible progress without being identified as physical time.
- **B — Additive dissipative clock:** geometric proper time remains the
  baseline and dissipation adds a clock-like contribution.
- **C — Strong emergent PLAY-time:** irreversible processes generate the
  physical progression represented by the candidate clock, with zero rate
  when the relevant irreversible production vanishes.

Status:

`SCIENTIFIC DECISION REQUIRED`

### FH-002 — Global playhead versus local scalar field

The project must choose whether the candidate clock is a spatially integrated
global playhead, a local scalar field along material worldlines, or a precisely
defined relationship between both.

Status:

`SCIENTIFIC DECISION REQUIRED`

### FH-003 — Physical status of Index I

The project must determine whether the Index is historical motivation, a
rigorously defined observable, an independent dynamical field, or a redundant
quantity expressible through established entropy or statistical observables.

Status:

`SCIENTIFIC DECISION REQUIRED`

### FH-004 — Physical status of the mesh

The project must determine whether the historical metaphor has any precise
physical referent. Absolute zero, the Cosmic Microwave Background, quantum
vacuum, spacetime, and the electromagnetic field must not be identified with
one another.

Status:

`SCIENTIFIC DECISION REQUIRED`

SD-001 through SD-012 remain valid audit findings. Cosmological resolution is
**PAUSED** until FH-001 and FH-002 are resolved.

---

## 9. Known cosmological model decisions requiring review

The following scientific decisions are intentionally unresolved.

An AI agent must NOT choose among alternatives autonomously.

---

### SD-001 — Entropy-production closure

The foundational model requires a scalar

\[
\Sigma\ge0,
\]

but does not yet approve a unique computational decomposition.

Possible thermal, viscous, gravitational, shock, information, or sub-grid
terms must each receive independent provenance and dimensional analysis.

In particular, any gravitational contribution proportional to quantities
such as

\[
\rho |\nabla\Phi|^2 / T
\]

must be classified as a Reotransductor-specific phenomenological ansatz
unless independently derived.

Status:

`SCIENTIFIC DECISION REQUIRED`

---

### SD-002 — Effective coupling \(\kappa_{\mathrm{eff}}\)

The natural dimensional scale

\[
\kappa_0
=
\frac{\hbar^2G^2}{c^7k_B}
\]

does not by itself determine the macroscopic simulation coupling.

Any:

- coarse-graining factor;
- gauge normalization;
- degree-of-freedom scaling;
- empirical multiplier;
- lower or upper clipping;

must be explicitly classified.

Status:

`SCIENTIFIC DECISION REQUIRED`

---

### SD-003 — Expansion variable

The cosmological engine historically uses a variable called `scale_factor`
with landmark values that do not trivially correspond to the conventional
FLRW normalization \(a(z=0)=1\).

The project must decide whether the numerical variable represents:

1. the physical FLRW scale factor; or
2. a distinct numerical expansion coordinate \(A\).

If it is \(A\), the mapping

\[
a_{\mathrm{FLRW}} = f(A)
\]

or equivalently

\[
z=f(A)
\]

must be explicitly specified.

Status:

`SCIENTIFIC DECISION REQUIRED`

---

### SD-004 — Cosmological expansion law

The exact production-engine equations governing:

- \(H\);
- scale-factor evolution;
- Hubble damping;
- cosmological cooling;
- Poisson scaling;

must be selected and then documented identically in Paper II and code.

Status:

`SCIENTIFIC DECISION REQUIRED`

---

### SD-005 — Sound-speed prescription

The production sound-speed law must be derived or explicitly classified.

Its:

- temperature dependence;
- floor;
- ceiling;
- code-unit mapping;
- physical velocity;
- resulting Jeans scale;

must be mutually consistent.

Historical values or numerical clamps are not automatically approved.

Status:

`SCIENTIFIC DECISION REQUIRED`

---

### SD-006 — Gravitational code-unit normalization

Any screened or effective value used for the Poisson coupling in code units
must have an explicit nondimensional derivation or phenomenological
classification.

It must not be described simply as Newton's \(G\) if it is not numerically
the nondimensional image of the physical constant.

Status:

`SCIENTIFIC DECISION REQUIRED`

---

### SD-007 — Temperature representation

The project must define one authoritative relation between the numerical
temperature field and physical Kelvin temperature.

Multiple incompatible mappings are not permitted in the frozen model.

Status:

`SCIENTIFIC DECISION REQUIRED`

---

### SD-008 — Sub-grid virial shock model

Any enhancement intended to recover unresolved entropy production must
define:

- activation criterion;
- resolution dependence;
- mathematical form;
- physical interpretation;
- calibration source;
- sensitivity;
- convergence behavior.

It must be labeled as a sub-grid closure, not a first-principles result.

Status:

`SCIENTIFIC DECISION REQUIRED`

---

### SD-009 — Initial-condition spectrum

The primordial density and velocity initial conditions, power spectrum,
normalization, random-phase procedure, and any synthetic fossil memory must
be specified before cosmological prediction claims are made.

Status:

`SCIENTIFIC DECISION REQUIRED`

---

### SD-010 — Information field \(I\)

Any independently evolved information variable and its Landauer-inspired
decay/sustenance dynamics must demonstrate what physical observable it
represents and whether it influences the foundational clock model.

If it has no required role, it may remain an experimental diagnostic rather
than part of the frozen cosmological model.

Status:

`SCIENTIFIC DECISION REQUIRED`

---

### SD-011 — Black-hole identification versus memory saturation

A numerical threshold on \(\Delta\tau_{\mathrm{diss}}\) must not by itself be
used as proof that a numerical cell physically represents a black hole.

If the simulation intends to identify black holes, an independent compactness
or horizon criterion must be defined.

Status:

`SCIENTIFIC DECISION REQUIRED`

---

### SD-012 — Physical time and grid resolution

The mapping between:

- code timestep;
- cell size;
- physical elapsed time;
- expansion evolution;

must be fixed before physical convergence studies are interpreted.

Status:

`SCIENTIFIC DECISION REQUIRED`

---

## 10. Legacy interpretations currently not approved

The following types of statement are explicitly not scientific authority
merely because they may appear in historical files.

### L1

"Numerical saturation means physical time freezes."

Status:

`LEGACY`

### L2

"A local black-hole/Planck-star bounce triggers a global eon transition."

Status:

`LEGACY`

### L3

"An observational fit proves or validates Reotransductor."

Status:

`LEGACY`

### L4

"All implemented cosmological terms are first-principles."

Status:

`LEGACY`

### L5

"A parameter is physically derived because its value yields the expected
astrophysical structure."

Status:

`LEGACY`

### L6

"Passing software tests establishes physical correctness."

Status:

`LEGACY`

---

## 11. Current approved scientific model

At `MODEL-SPEC-0.2`, no project-specific clock interpretation is frozen.
The additive formulation

\[
u^\mu\nabla_\mu\tau_R
=
1+\kappa\Sigma,
\qquad
\Sigma\ge0,
\]

is retained as candidate FH-001B together with the explicit null hypothesis

\[
\kappa=0.
\]

Candidates FH-001A and FH-001C remain open, and the global-versus-local choice
in FH-002 is unresolved. No unique foundational or full cosmological
realization has been frozen.

This limitation is intentional.

The existing cosmological engine is preserved as a valuable historical and
experimental implementation to be audited component by component.

---

## 12. Promotion rule

A component may move from `UNDER_REVIEW` or `EXPERIMENTAL_EXTENSION` to
`APPROVED_EXTENSION` only after documenting:

1. exact equations;
2. variable definitions;
3. dimensional consistency;
4. physical or phenomenological provenance;
5. numerical implementation;
6. tests;
7. sensitivity to free parameters and regulators;
8. expected falsifiable consequence;
9. compatibility with the rest of the approved model;
10. explicit project-owner approval.

Implementation existence alone is never sufficient.

---

## 13. Model freeze criteria

A cosmological model version may be declared frozen only when:

1. FH-001 and FH-002 have been explicitly resolved;
2. all governing equations are explicitly specified;
3. all scientifically relevant parameters are classified;
4. unit mappings are unique and documented;
5. unresolved `SCIENTIFIC DECISION REQUIRED` items affecting the production
   model have been resolved;
6. Paper II matches the executable equations;
7. production tests use production parameters;
8. numerical convergence has been tested at physically equivalent states;
9. checkpoint manifests reproduce the full scientific configuration;
10. no observational result is being used to retroactively tune the frozen
   equations;
11. the exact model receives a version identifier.

Only after this freeze should downstream observational validation be treated
as a test of that model version.

---

## 14. Required audit output before code remediation

Before scientifically modifying the current cosmological engine, perform a
read-only audit producing a matrix with columns:

| ID | Component | Foundational requirement | Current code behavior | Paper I | Paper II | Status | Required decision |
|---|---|---|---|---|---|---|---|

Every discrepancy must be classified according to `AGENTS.md`.

The audit must not fix the discrepancies.

Its purpose is to establish the starting point for finite, traceable model
decisions.

---

## 15. Relationship to Paper III

Paper III is intentionally excluded from the current model specification
phase.

No observational agreement or disagreement from Paper III may determine the
resolution of SD-001 through SD-012 unless the modification is explicitly
classified as post-hoc exploratory model development.

After a model version is frozen, observational testing may resume using
pre-specified metrics and independent evidence.

---

## 16. Current epistemic statement

Reotransductor currently contains:

- a reopened foundational decision register;
- a substantial experimental cosmological codebase;
- multiple physically interesting extensions;
- unresolved model choices;
- numerical and theoretical work requiring verification.

It does NOT yet contain a frozen, fully validated cosmological theory.

Current status:

**MODEL-SPEC-0.2 — FOUNDATIONAL RECONCILIATION IN PROGRESS /
COSMOLOGICAL MODEL NOT FROZEN**
