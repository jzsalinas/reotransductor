# Cusp-Core Laboratory — Experiment 002 Preregistration

## 1. Experiment identity

**Experiment:** Cusp-Core Laboratory — Experiment 002

**Relationship:** Fresh repetition of the same frozen physical experiment
after correction of independently diagnosed numerical defects in CCL-NUM-1.

**Physical model:** Unchanged from Experiment 001.

**Numerical version:** `CCL-NUM-2`.

Experiment 001 remains permanently classified as:

```text
CCL-NUM-1
Experiment 001
Outcome = SCIENTIFIC EXPERIMENT INVALID — NUMERICAL FAILURE
```

Experiment 001 is not reinterpreted by this preregistration.

## 2. Frozen experiment configuration

Frozen conductivity:

```text
K_hat = 1.5
```

Frozen resolutions:

```text
N = 512
N = 1024
N = 2048
```

Frozen output times:

```text
t/t0 = 0
       0.25
       0.5
       1.0
       2.0
```

Frozen physical scales:

```text
r_s = 1 kpc
rho_s = 1 Msun pc^-3
```

These scale definitions fix the physical unit mapping. They are not an
observational normalization and do not constitute a first-principles
prediction of absolute galaxy rotation speeds.

## 3. Frozen scientific criteria

```text
mass relative error <= 1e-12
total-energy relative error <= 1e-6

r_min = 8 dr

resolved core requires:
    r_core / dr >= 16

core boundary:
    gamma(r_core) = -0.5

convergence criteria exactly as already specified in SPECIFICATION.md

flat rotation-curve region requires:
    |beta| < 0.1
continuously for at least a factor two in radius.
```

Here `gamma=d ln(rho)/d ln(r)` and `beta=d ln(v_c)/d ln(r)`. The rotation
curve is calculated only from the simulated mass distribution through

\[
v_c(r)=\sqrt{G M(<r)/r}.
\]

No observational normalization is permitted. No SPARC data may be used
during the experiment.

## 4. Pre-registered outcome classification

No additional outcome category may be invented after seeing the result.

### A. CONDUCTIVE CORE FORMATION SUPPORTED

Requires all of the following:

- the `K_hat=0` control remains cuspy;
- the conductive branch develops a resolved core;
- `r_core/dr >= 16`;
- `r_core` approaches a nonzero physical value with resolution;
- gamma and density profiles converge;
- mass conservation passes;
- total-energy conservation passes;
- conductive energy redistribution causally accounts for the structural
  evolution.

### B. NO CORE IN THIS PREREGISTERED CONDUCTIVE REGIME

Use if all numerical requirements pass but `K_hat=1.5` does not produce a
resolved convergent core through `2 t0`.

This rejects only this exact preregistered regime.

### C. APPARENT CORE IS NUMERICAL / NON-CONVERGENT

Use if apparent flattening disappears under refinement, scales with `dr`, or
fails the frozen convergence criteria.

### D. SCIENTIFIC EXPERIMENT INVALID — CONSERVATION FAILURE

Use if either frozen conservation requirement fails.

### E. SCIENTIFIC EXPERIMENT INVALID — NUMERICAL FAILURE

Use if a non-conservation numerical failure prevents a valid scientific
classification under A, B, or C.

## 5. Experiment 001 firewall

The transient result from Experiment 001,

```text
gamma(r_min) ~ -0.404
central-density reduction ~65.8%
```

must not be used as a target, benchmark, expected result, acceptance
criterion, or tuning reference. Experiment 002 must be evaluated
independently. No attempt may be made to reproduce that transient.

## 6. Verified instrument status before execution

The numerical instrument was verified before this preregistration:

```text
CCL-NUM-2
29/29 tests passed
K_hat=0 control validated at N=512,1024,2048 through 2 t0
maximum reported total-energy error = 1.54e-15
maximum reported mass error = 0
no resolved core
positivity limiter activations = 0
```

This is numerical verification and control validation, not evidence for
conductive core formation. Experiment 002 has not been run at the time of
this preregistration.
