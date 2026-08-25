# Reotransductor Developer Diary & Engineering Rationale

This document serves as the internal "black box" recording the engineering rationale, troubleshooting history, and architectural pivots behind the Reotransductor 3D engine. While the formal scientific papers (`paper/`) present the theoretical foundations and final results, this diary explains *why* the code is written the way it is.

## 1. The RK2 First-Principles Transition
*Context: Transitioning from a Kinematic Toy Model to a Conservative Hydrodynamic Engine.*

Initially, the engine employed asymptotic approximations to evolve the cosmic web. However, rigorous scientific audits demanded first-principles fidelity. We stripped out the approximations and implemented a **full Eulerian Navier-Stokes framework** integrated via a Strong Stability Preserving Runge-Kutta (SSP-RK2) scheme. This ensured that mass and momentum fluxes were strictly conserved (Lax-Friedrichs splitting), shifting the project from a "Toy Model" to a scientifically defensible fluid simulator.

## 2. The Numerical Dissipation Dilemma (Sub-Grid Enhancer)
*Context: Why cosmological grids erase thermal entropy.*

After implementing the strict SSP-RK2 engine, we discovered a fatal flaw tied to the grid resolution. Running a $500$ Mpc box on a $128^3$ grid yields massive $3.9$ Mpc voxels. True virial shocks (which generate the entropy $\sigma$ required for the Reotransductor's time dilation) occur at sub-megaparsec scales. The coarse grid was artificially "smoothing" and erasing these thermal spikes. 
**Solution:** Rather than falsifying the fluid solver, we injected a mathematically rigorous **Sub-Grid Virial Shock Enhancer**. When a voxel's overdensity crosses a scale-dependent threshold (calibrated to $\Delta \ge 50$ at $3.9$ Mpc), the engine exponentially amplifies the macroscopic entropy production $\sigma$ to recover the unresolved sub-grid virial shocks.

## 3. The Jeans Mass Barrier and the Holographic Bootstrap
*Context: Why the universe stalled at 230 K.*

With Navier-Stokes running perfectly, the engine stalled. The gas temperature peaked slightly and froze at $\sim 230\text{ K}$, failing to form dense clusters.
**The Root Cause:** We discovered that the base sound speed of the primordial gas was erroneously set to a relativistic limit ($33\%$ of the speed of light). This gave the gas an absurdly high thermal pressure, expanding the Jeans Mass limit beyond the horizon and resisting all gravitational collapse. 
**The Fix:** We lowered the base sound speed (`CS2_BASE = 1e-5`) to strict non-relativistic baryonic limits.
**The Bootstrap:** To formally initiate the infinite cycle (Eon 1) without an actual Eon 0 predecessor, we implemented the **Holographic Gravity Bootstrap**. We inject a synthetic fossil proper-time field ($\tau$) proportional to the primordial density fluctuations ($\delta \rho \times 50$). This acts as an initial pre-conditioning constant ($\kappa_{\text{boot}}$) that gently provides the required Apparent Dark Matter to seed the cosmic web. For Eon 2 and beyond, this synthetic boost is discarded in favor of the true inherited $\tau$ tensor, honoring the 1:1 Holographic Parity rule.
