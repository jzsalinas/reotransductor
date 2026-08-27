# Frozen Physical and CCL-NUM-2 Numerical Specification — Independent Cusp-Core Laboratory

Numerical version: `CCL-NUM-2`

`Experiment 001` remains permanently identified with `CCL-NUM-1` and
`OUTCOME = SCIENTIFIC EXPERIMENT INVALID — NUMERICAL FAILURE`. Its artifacts
under `results/8786e1f_khat1p5/` are immutable historical evidence and are not
reinterpreted by this revision.

## 1. Scientific scope

The laboratory contains one Newtonian, collisional, calorically perfect,
monatomic ideal gas in spherical symmetry. Its established interactions are
self-gravity and Fourier heat conduction. It contains no physical viscosity,
radiative cooling, external heating, cosmological expansion, observational
normalization, synthetic gravitational source, or Reotransductor-specific
physics.

The conductive branch asks only whether this specified mathematical fluid and
transport regime can transform the specified cusp into a numerically resolved
core. It does not identify the fluid with SIDM, baryonic plasma, or dark matter.

## 2. Exact PDE system

For density `rho`, radial velocity `u`, pressure `P`, specific internal energy
`e`, gas energy density `E`, potential `Phi`, temperature-like variable
`Theta`, and radial heat flux `q`,

\[
 P=(\gamma-1)\rho e,\qquad
 E=\rho e+\tfrac12\rho u^2,\qquad
 \Theta=P/\rho,\qquad \gamma=5/3,
\]

\[
 \partial_t\rho+\frac1{r^2}\partial_r(r^2\rho u)=0,
\]

\[
 \partial_t(\rho u)
 +\frac1{r^2}\partial_r\left[r^2(\rho u^2+P)\right]
 =\frac{2P}{r}-\rho\,\partial_r\Phi,
\]

\[
 \partial_t E
 +\frac1{r^2}\partial_r\left[r^2u(E+P)\right]
 =-\rho u\,\partial_r\Phi
 -\frac1{r^2}\partial_r(r^2q),
\]

\[
 \frac1{r^2}\partial_r(r^2\partial_r\Phi)=4\pi G\rho,
 \qquad
 \partial_r\Phi=\frac{G M(<r)}{r^2},
 \qquad
 M(<r)=4\pi\int_0^r\rho(s)s^2ds,
\]

\[
 q=-K\,\partial_r\Theta,\qquad
 L_{\rm cond}=4\pi r^2q.
\]

All conductive updates use conservative face luminosities. No density floor,
pressure floor, energy floor, clamp, ceiling, velocity clipping, artificial
density diffusion, temperature relaxation, or hidden source is permitted. A
non-positive density or internal energy is an explicit numerical failure.

## 3. Units and parameter classification

The code is nondimensional with

\[
 L_0=r_s,\quad \rho_0=\rho_s,\quad M_0=\rho_s r_s^3,
 \quad v_0=\sqrt{G\rho_s r_s^2},\quad t_0=(G\rho_s)^{-1/2}.
\]

Thus `G_hat=1`. The reference physical realization is
`r_s=1 kpc`, `rho_s=1 Msun pc^-3`. These are `INITIALIZATION` scale
definitions, not observational fits. They determine `v0` and therefore do not
constitute a first-principles prediction of absolute galaxy rotation speeds.

The controlled conductivity is

\[
 \widehat K=\frac{K}{\rho_s r_s v_0},\qquad
 \widehat K_\star=1.5.
\]

It is a pre-registered controlled transport parameter, not a microphysical
conductivity. The reference-scale time is

\[
 t_{\rm cond,0}=\frac{t_0}{(\gamma-1)\widehat K}.
\]

At `K_hat_star=1.5`, `t_cond,0=t0`. This is not the local conduction time at
`r_s`; using `rho(r_s)=rho_s/4` and length `r_s` gives the approximate local
scale `t_cond(r_s;r_s)=t_cond,0/4`.

## 4. Initial equilibrium and weak central regularity

On `0 <= r <= R=10 r_s`, with `x=r/r_s`,

\[
 \rho(r)=\frac{\rho_s}{x(1+x)^2},\qquad u(r)=0,
\]

\[
 \frac{dP}{dr}=-\rho\frac{GM(<r)}{r^2}.
\]

The outer pressure is fixed by `dTheta/dr=0` at `R`. Explicitly,

\[
 \Theta(R)=-\frac{g(R)}{d\ln\rho/dr|_R},\qquad P(R)=\rho(R)\Theta(R),
\]

and `P(r)=P(R)+integral_r^R rho(s)g(s) ds`.

As `r -> 0`, `rho ~ rho_s r_s/r`, `M ~ 2 pi rho_s r_s r^2`,
`g -> 2 pi G rho_s r_s`, `P ~ C ln(r_s/r)`, and
`Theta ~ g(0) r ln(r_s/r)`. Consequently `dTheta/dr` and `q` diverge
logarithmically but `L_cond ~ r^2 ln(r_s/r) -> 0`. The initial conductive
problem is weakly regular/integrable, not classically pointwise regular.

Every initial conserved quantity is a finite-volume cell average of the
continuous singular state. Density averages use exact NFW mass differences.
Pressure averages use the identity

\[
 \int_a^b r^2P\,dr
 =\frac{[r^3P]_a^b}{3}+\frac13\int_a^b r^3\rho g\,dr,
\]

so `P(0)` is never evaluated. Neither `rho(0)` nor `q(0)` is evaluated.
The central boundary condition is `L_cond(0)=0`; it is not the false
pointwise assertion `lim q(r)=0`.

## 5. Boundary conditions

Both radial boundaries are impermeable to mass and advective energy. The
origin uses spherical reflection through its zero-area face. The outer face is
a rigid reflecting wall at `R=10`. Conductive luminosity is zero at both
boundaries, so conduction only redistributes internal energy inside the
domain.

## 6. Numerical method

- Uniform radial spherical finite volumes in `float64`.
- Exact cell volumes `4 pi (r_+^3-r_-^3)/3` and face areas `4 pi r_f^2`.
- Piecewise-linear MUSCL reconstruction of conserved-variable perturbations
  about the exact NFW equilibrium. Linear perturbations use the exact spherical
  volume centroid
  `r_bar,V=(3/4)(r_+^4-r_-^4)/(r_+^3-r_-^3)`.
- Monotonized-central (MC) slopes for all conserved perturbations.
- Conservative positivity scaling at every used reconstruction point,
  `U_limited=U_bar+alpha(U_raw-U_bar)`, using the largest common
  `0<alpha<=1` that gives `rho>0` and `rho E-(rho u)^2/2>0`. The same `alpha`
  scales every conserved component. Cell averages are unchanged. A
  nonphysical cell average is an immediate numerical failure.
- HLLC advective Riemann flux.
- Well-balanced equilibrium subtraction: numerical fluxes and source terms
  are evolved relative to the same hydrostatic NFW equilibrium. The exact
  initialized state therefore has a zero semidiscrete residual while the
  correction vanishes under spatial refinement for general states.
- Self-gravity from the conjugate potential of the same exact discrete
  piecewise-constant-shell energy functional used in conservation accounting;
  no external or synthetic source.
- Second-order SSP-RK2 for the Euler-Poisson operator.
- Conservative matrix-exponential conduction for `H dTheta/dt=A Theta`.
  Production actions use the symmetric similar tridiagonal operator
  `B=H^-1/2 A H^-1/2` and adaptive, fully reorthogonalized Lanczos projection.
  Only the small Krylov projection is dense; no production-sized dense
  exponential is constructed.
- Symmetric conduction-half / hydro-full / conduction-half Strang splitting.

The weighted constant thermal mode is split explicitly before the Lanczos
action. The fixed float64 exponential-action contract is:

- relative a-posteriori action tolerance `5e-13`;
- absolute a-posteriori action tolerance `5e-15`;
- weighted thermal-invariant relative tolerance `5e-13`;
- constant-mode preservation relative tolerance `5e-14`;
- positivity numerical allowance `5e-13 * max(1,max|Theta_initial|)`;
- discrete maximum-principle allowance
  `5e-13 * max(1,max|Theta_initial|)`;
- initial/increment/maximum Krylov dimensions `16/16/256`.

The action is accepted only when the successive-projection a-posteriori
difference and the Lanczos residual estimate meet the contract (or an
invariant Krylov subspace terminates exactly), and all invariant, constant,
positivity, and maximum-principle checks pass. Temperature is never clipped.
Negativity beyond the certified allowance is `NUMERICAL FAILURE`; any
non-positive state, even within an allowance, is unusable and also stops.

## 7. Conservation accounting

The implementation records independently

\[
 M=\sum_i\rho_iV_i,
\quad U=\sum_i\rho_i e_iV_i,
\quad T=\sum_i\tfrac12\rho_i u_i^2V_i,
\]

\[
 W=-\int_0^R\frac{GM(r)}r\,dM(r),\qquad E_{\rm total}=U+T+W,
\]

and every conductive face luminosity. `W` is integrated exactly for the
piecewise-constant finite-volume density. The required tolerances are relative
mass error `<=1e-12` and relative total-energy error `<=1e-6`. They may not be
relaxed based on development outcomes without a new scientific/numerical
decision.

Writing cell masses as `m_i=rho_i V_i`, the implemented quadratic functional
is

\[
 W_h=-\sum_i\left(a_i m_i M_i+b_i m_i^2\right),\qquad
 M_i=\sum_{j<i}m_j,
\]

where

\[
 a_i={4\pi\over V_i}{r_{i+1}^2-r_i^2\over2},\qquad
 b_i={1\over3}\left({4\pi\over V_i}\right)^2
 \left[{r_{i+1}^5-r_i^5\over5}
 -{r_i^3(r_{i+1}^2-r_i^2)\over2}\right].
\]

Its conjugate potential is

\[
 \psi_k={\partial W_h\over\partial m_k}
 =-\left(a_kM_k+2b_km_k+\sum_{i>k}a_i m_i\right).
\]

For the RK2-integrated face mass transfer `F_m,f`, `psi` is evaluated at the
midpoint cell mass. The quadratic identity is then exact:

\[
 \Delta W_h=\sum_f F_{m,f}(\psi_R-\psi_L).
\]

The gas receives the local face work
`-F_m,f(psi_R-psi_L)`, shared equally by its adjacent finite volumes. This is
part of the timestep construction, not a post-step correction or energy
renormalization. The same first-stage construction supplies the RK2 predictor.
The gravitational momentum source uses the gradient of this conjugate
potential with the frozen equilibrium subtraction.

## 8. Controlled experiment and protected execution

The branches have identical equations, initial state, grid, boundary
conditions, solver, and diagnostics:

- `CONTROL`: `K_hat=0`.
- `CONDUCTIVE`: `K_hat=K_hat_star=1.5`.

The control is evaluated at `t/t0={0.25,0.5,1.0,2.0}` on
`N_r={512,1024,2048}`. It must satisfy the conservation tolerances, have
`max|u|/v0<=1e-10`, volume-weighted relative `L1` density change `<=1e-10`,
maximum resolved slope change `<=1e-8`, and form no resolved core.

The final conductive experiment at `N_r={512,1024,2048}` is protected and is
not executed as part of implementation or control validation.

## 8.1 CCL-NUM-2 verification gate

Before freezing this numerical version, the implementation must pass the
pre-registered tests for constant temperature, positive Gaussian diffusion,
spherical analytic diffusion, thermal conservation, zero origin luminosity,
stiffness through `mu=1e5`, the central NFW weak state, temporal and spatial
convergence, coupled hydro-gravity energy, extreme positive reconstruction,
exact NFW equilibrium, locally convergent collapse/expansion work, absence of
silent repair, dense-reference comparison, and thermal invariant/constant
mode preservation. Verification is software evidence only, not physical
validation.

## 9. Resolution, slope, and core protocol

All comparisons use identical dimensionless physical times, never equal step
counts. The minimum analyzed radius is

\[
 r_{\min}=8\,\Delta r.
\]

The logarithmic slope estimator is a second-order finite difference in
`ln(rho)` versus `ln(r)`:

\[
 \gamma(r)=\frac{d\ln\rho}{d\ln r}.
\]

The core radius is the first outward interpolated crossing
`gamma(r_core)=-0.5`, starting outside `r_min`. A claimed core additionally
requires `r_core/Delta r >= 16` at every resolution used for that claim.

For the future conductive experiment, convergence requires between the two
finest grids at identical times: relative `r_core` difference `<=10%` and
absolute slope difference `<=0.05` at matched radii, in addition to both
conservation tolerances. Failure of an estimator to exist at one resolution is
non-convergence, not a zero core radius.

## 10. Rotation curve

No normalization to data is allowed. The prediction conditional on the
simulated mass distribution is

\[
 v_c(r)=\sqrt{\frac{G M(<r)}{r}}.
\]

Physical radii and speeds are obtained only through the fixed scale mapping in
Section 3. No SPARC curve or other observational input enters this calculation.

## 11. Falsification scope

The pre-registered conductive regime fails if it does not produce a resolved,
converged core; if the control also forms a core; if control stationarity
fails; if mass or energy tolerances fail; if inferred core size follows grid
scale; or if density/internal energy becomes non-positive.

Failure at `K_hat=1.5` rejects this pre-specified transport regime. It does not
prove that all conducting self-gravitating fluids cannot form cores. Conversely,
a converged core would be a numerical consequence of this generic model, not
evidence that a specific material or dark-matter microphysics has been found.
