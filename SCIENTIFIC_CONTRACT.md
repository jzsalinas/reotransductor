# Reotransductor — Scientific Contract

## 1. Purpose

This document defines the methodological rules governing scientific work in the Reotransductor project.

These rules apply to humans, AI agents, numerical experiments, manuscripts, and interpretation of results.

The purpose of this contract is to prevent the project from drifting from hypothesis testing into hypothesis confirmation.

The foundational question is defined separately in:

`FOUNDATIONAL_HYPOTHESIS.md`

This document does not establish whether that hypothesis is correct.

It defines how the hypothesis must be investigated.

---

## 2. Primary scientific principle

The project must attempt to determine whether the Reotransductor hypothesis is false as seriously as it attempts to determine whether it may be viable.

The desired outcome of an experiment must never influence the governing equations used to produce that experiment.

A failed prediction, null result, numerical instability, incompatibility with observation, or exclusion of the model is a scientifically valid outcome.

No code change may be introduced solely to rescue the hypothesis from an unfavorable result.

---

## 3. Pre-specification before observation

Whenever reasonably possible, the following must be specified before examining the observational outcome of a scientific test:

- governing equations;
- free parameters;
- numerical regulators;
- initial conditions;
- priors;
- selection criteria;
- fitting procedure;
- comparison metric;
- uncertainty model;
- success/failure criterion.

If any of these are changed after examining results, the change must be explicitly classified as post-hoc.

The same result used to motivate a post-hoc modification must not subsequently be presented as an independent prediction of the modified model.

A new independent test or held-out dataset is required for such a claim.

---

## 4. Separation of development and validation

Reotransductor work must distinguish at least three phases.

### 4.1 Development

Used to construct, debug, calibrate, or explore the model.

Results obtained during development may guide engineering decisions but are not independent validation evidence.

### 4.2 Verification

Used to determine whether the numerical implementation correctly solves the model that was specified.

Verification asks:

> Are we solving the equations right?

Examples include:

- unit tests;
- dimensional checks;
- conservation tests;
- manufactured solutions;
- convergence studies;
- CPU/GPU equivalence;
- checkpoint reproducibility.

Verification does NOT establish that the underlying physical model is true.

### 4.3 Validation

Used to determine whether the specified model adequately describes physical reality.

Validation asks:

> Are we solving the right equations?

Validation requires comparison with independent experimental or observational evidence, including uncertainty and appropriate baseline models.

Software verification must never be reported as physical validation.

---

## 5. Hierarchy of scientific claims

Every scientific statement must be classified implicitly or explicitly according to its evidentiary level.

From weakest to strongest:

1. IDEA
2. HYPOTHESIS
3. POSTULATE
4. PHENOMENOLOGICAL MODEL
5. NUMERICAL CONSEQUENCE
6. PREDICTION
7. OBSERVATIONAL COMPATIBILITY
8. EMPIRICAL SUPPORT
9. ESTABLISHED RESULT

A claim may only move upward when the corresponding evidence exists.

In particular:

- a postulate is not a derivation;
- a fit is not a prediction;
- compatibility is not confirmation;
- lower chi-square is not automatically evidence for new physics;
- successful code execution is not empirical support;
- reproduction of a known pattern is not proof of the proposed mechanism.

---

## 6. Classification of quantities and parameters

Every scientifically relevant parameter or term must be classified as one of:

### FUNDAMENTAL

Directly corresponding to an established physical constant or experimentally measured quantity.

Examples may include \(c\), \(G\), \(k_B\), or \(\hbar\).

### DERIVED

Obtained mathematically from clearly stated assumptions and quantities without empirical tuning.

### PHENOMENOLOGICAL

Introduced to represent unresolved or hypothesized physics and constrained or calibrated empirically.

### NUMERICAL

Introduced for discretization, stability, regularization, finite precision, convergence, runtime, or computational feasibility.

### INITIALIZATION

Used to construct initial or boundary conditions but not derived from the dynamical equations.

### FITTED

Estimated directly from the data being analyzed.

No phenomenological, numerical, initialization, or fitted quantity may be presented as fundamental or first-principles.

---

## 7. No hidden numerical physics

The following numerical mechanisms must never be introduced silently:

- floors;
- ceilings;
- clipping;
- saturation;
- boosts;
- damping factors;
- normalization constants;
- smoothing kernels;
- sub-grid enhancements;
- artificial diffusion;
- artificial viscosity;
- threshold triggers;
- empirical rescaling.

Every such mechanism must document:

1. why it exists;
2. whether it is numerical or physical;
3. its units;
4. its effect on the equations;
5. its sensitivity range;
6. whether results survive reasonable variation or removal.

If a regulator materially changes a scientific observable, that dependence must be reported.

---

## 8. Dimensional and unit discipline

Every new physical equation must pass dimensional consistency checks before implementation.

Code units must have an explicit mapping to physical units.

Dimensionless variables must be identified as dimensionless.

A quantity must not be assigned a physical interpretation merely because its numerical magnitude resembles an astrophysical scale.

Conversion factors must be traceable.

No undocumented unit conversion constant may be introduced.

If several mappings from code units to physical units are possible, the ambiguity must be reported as:

`SCIENTIFIC DECISION REQUIRED`

---

## 9. Established physics versus Reotransductor physics

Established physical laws must be separated from project-specific hypotheses.

When an equation is standard, its provenance should be identified where relevant.

When a term is specific to Reotransductor, it must be described as such.

A citation to established physics must not be used to imply that the cited source derives a Reotransductor-specific equation unless it actually does.

In particular, intellectual inspiration does not constitute derivation.

---

## 10. No retrospective theoretical justification

If the code contains an unexplained constant, equation, threshold, or behavior, the agent must not invent a physical derivation for it.

The correct response is to classify the item as one of:

- undocumented numerical choice;
- phenomenological assumption;
- legacy implementation;
- unresolved scientific choice;
- implementation bug.

A physical justification may only be added after it has been independently established and explicitly approved.

---

## 11. Baseline models

Every claimed new physical effect must be compared against an appropriate null or standard baseline.

For the foundational Reotransductor hypothesis, the primary null model is:

\[
\kappa = 0.
\]

When applicable, comparisons should also include relevant conventional physical models.

A new model must not be declared successful merely because it qualitatively resembles observations.

The comparison must establish what additional explanatory or predictive content the new hypothesis provides.

---

## 12. Parameter fitting and predictive independence

Parameters fitted using a dataset cannot be treated as independent predictions for that same dataset.

The project must distinguish between:

- calibration data;
- training/development data;
- validation data;
- held-out or genuinely predictive data.

If the same dataset is reused, the manuscript must state this explicitly.

Degrees of freedom introduced through fitting must be included when comparing models.

---

## 13. Statistical discipline

Scientific comparisons must report appropriate uncertainty.

When applicable, analyses should consider:

- observational errors;
- covariance;
- numerical uncertainty;
- finite-volume uncertainty;
- resolution dependence;
- parameter uncertainty;
- model complexity.

Reduced chi-square, likelihood, residuals, information criteria, posterior distributions, or other statistics must be interpreted according to their assumptions.

Statistical significance must not be inferred from a metric that does not support such an interpretation.

---

## 14. Numerical convergence

A scientifically relevant numerical result must not be considered robust solely because the simulation remains numerically stable.

Convergence must compare physically equivalent states.

Grid-resolution studies must account for resolution-dependent physical time steps, spatial scale, regulator behavior, and sub-grid prescriptions.

Where possible, convergence order or quantitative error reduction should be estimated.

If the observable fails to converge, this must be reported.

---

## 15. Reproducibility

A scientific simulation result intended for publication must be reproducible from a documented configuration.

The reproducibility record should include, where applicable:

- source-code revision;
- model specification revision;
- random seed;
- grid resolution;
- physical box size;
- timestep;
- physical constants;
- free parameters;
- numerical regulators;
- hardware/backend;
- initial-condition generator;
- checkpoint provenance.

Incompatible checkpoints must fail explicitly rather than resume silently.

---

## 16. Preservation of negative results

Negative and null results must not be deleted merely because they contradict expectations.

Failed experiments that materially influence model development should be preserved in appropriate logs, reports, or version history.

The scientific record should make it possible to distinguish:

- hypotheses considered;
- tests attempted;
- failures observed;
- changes subsequently introduced.

---

## 17. Paper governance

### Paper I

Must define the foundational theoretical hypothesis and clearly labeled extensions.

It must distinguish established physics from Reotransductor postulates.

### Paper II

Must document the numerical equations actually executed by the production engine.

Idealized equations may be included for explanation, but any difference from the implemented system must be stated.

### Paper III

Must be treated as downstream empirical testing.

Observational outcomes from Paper III must not retroactively alter Papers I or II merely to improve agreement.

If Paper III motivates a new model version, the revised model must be frozen and tested again using independent evidence.

---

## 18. Model versioning

A scientifically meaningful change to governing equations, physical assumptions, or phenomenological parameters defines a new model version.

Results from different model versions must not be mixed without explicit labeling.

A numerical bug fix that changes previous scientific outputs must document which previous results are invalidated.

A purely engineering change that leaves numerical results unchanged should be distinguished from a scientific model revision.

---

## 19. Stop conditions for AI agents

An AI agent must stop and report:

`SCIENTIFIC DECISION REQUIRED`

rather than modifying physics autonomously when:

- manuscript and code disagree about the governing model;
- multiple physically plausible interpretations exist;
- a constant lacks provenance;
- a unit mapping is ambiguous;
- a proposed fix would change scientific predictions;
- a requested implementation requires choosing between competing physical models;
- observational agreement is being used as the primary justification for changing an equation.

The report must explain the alternatives without selecting the scientifically preferred result on behalf of the project owner.

---

## 20. Scientific audit rule

Audits should preferentially be differential after a model version has been frozen.

A differential audit asks:

> What scientifically relevant behavior changed since the last approved model state?

Full-repository audits remain appropriate after major architectural changes or when provenance is uncertain.

The objective is to achieve finite, traceable scientific iterations rather than repeated uncontrolled reinterpretation of the entire project.

---

## 21. Publication standard

Any future manuscript must permit a skeptical independent researcher to distinguish:

- established theory;
- new hypothesis;
- phenomenological closure;
- numerical approximation;
- parameter calibration;
- model prediction;
- observational comparison.

Important limitations must appear in the main scientific narrative rather than being hidden only in supplementary documentation or source code.

---

## 22. Final rule

The project must never ask:

> What modification makes Reotransductor agree with the data?

without first asking:

> Was that modification independently justified before seeing the result?

When the answer is no, the change must be classified as exploratory or post-hoc.

Scientific credibility takes priority over preserving the hypothesis.
