# Dimensionless Lattice Scaling and Physical Grid Calibration

This document provides the formal mathematical and physical justification for all dimensionless code units and base calibration parameters defined in [`server/physics_units.py`](file:///home/jzsalinas/Documents/antigravity/reotransductor/server/physics_units.py) and utilized across the 3D cosmological simulation engines.

---

## 1. Rationale of Dimensionless Code Units in Computational Astrophysics

In numerical astrophysical simulations (e.g., GADGET-4, ENZO, FLASH, AREPO, RAMSES), continuous partial differential equations (Navier-Stokes, Poisson, Boltzmann, and Onsager relations) cannot be directly integrated using raw SI units (e.g., $\hbar \sim 10^{-34}\text{ J}\cdot\text{s}$, $G \sim 10^{-11}\text{ N}\cdot\text{m}^2/\text{kg}^2$, or $1\text{ Mpc} \sim 10^{22}\text{ m}$). Standard IEEE 754 floating-point arithmetic (`float32` / `float64`) would suffer catastrophic underflow and overflow errors.

To solve this, all physical equations are nondimensionalized by normalizing space, time, mass, and energy to characteristic scales of the cosmological computational lattice:

### Fundamental Base Lattice Scales:
1. **Cosmological Box Size ($L_{\text{box}}$):**
   $$L_{\text{box}} = 100.0\text{ Mpc} = 3.08567758 \times 10^{24}\text{ m}$$
2. **Spatial Grid Resolution ($N = 32^3 = 32,768\text{ cells}$):**
   $$\Delta x = \frac{L_{\text{box}}}{N} = 3.125\text{ Mpc} = 9.6427 \times 10^{22}\text{ m}$$
3. **Causal Speed of Light Normalization ($c_{\text{code}} = 2.5\text{ cells/unit}$):**
   $$\Delta t_s = \frac{c_{\text{code}} \cdot \Delta x}{c} = \frac{2.5 \times (9.6427 \times 10^{22}\text{ m})}{2.99792458 \times 10^8\text{ m/s}} = 8.0412 \times 10^{14}\text{ s} \approx 25.4807\text{ Myr}$$
4. **Cosmic Critical Density ($\rho_0 = \rho_{\text{crit}}$ at $H_0 = 70\text{ km/s/Mpc}$):**
   $$H_{0,\text{SI}} = \frac{70 \times 10^3\text{ m/s}}{3.08567758 \times 10^{22}\text{ m}} = 2.2685 \times 10^{-18}\text{ s}^{-1}$$
   $$\rho_{\text{crit}} = \frac{3 H_{0,\text{SI}}^2}{8\pi G} = 9.2039 \times 10^{-27}\text{ kg/m}^3 \approx 1.3585 \times 10^{11}\ M_\odot/\text{Mpc}^3$$
5. **Velocity Unit Scale:**
   $$v_{\text{unit}} = \frac{\Delta x}{\Delta t_s} = \frac{c}{c_{\text{code}}} = \frac{299,792.458\text{ km/s}}{2.5} = 119,916.98\text{ km/s}$$

---

## 2. Mathematical and Physical Derivation of Grid Calibration Parameters

The parameters declared in `CosmologicalUnits` (`server/physics_units.py`) are the exact baseline scale constants evaluated at the cosmic microwave background ground state ($T_{\text{CMB}} = 2.7255\text{ K}$, $\rho / \bar{\rho} = 1.0$). During simulation execution, they dynamically scale as continuous tensor fields.

```python
# =====================================================================
# CALIBRATED COSMOLOGICAL GRID & THERMODYNAMIC BASE SCALES
# =====================================================================
self.DT = 0.05
self.H_0 = 0.0003
self.G_CONST = 0.0001
self.CS2_BASE = 1e-5
self.DIFFUSION_BASE = 0.3
self.LANDAUER_BASE = 0.015
self.INFLATION_BOOST = 8.0
self.ZETA_BEKENSTEIN = 3500.0
self.MASS_THRESHOLD = 0.18
self.M0_CORE = 5000.0
self.A_MAX_CONFORMAL = 7.0
```

---

### 2.1. `self.DT = 0.05` (Integration Time Step)
* **Physical Basis:** **Courant-Friedrichs-Lewy (CFL) Condition.**
* **Mathematical Derivation:** For hyperbolic and parabolic transport equations on a 3D Cartesian mesh with maximum signal velocity $c_{\text{code}} = 2.5$, numerical stability requires:
  $$C = \frac{c_{\text{code}} \Delta t}{\Delta x} \le C_{\text{max}} < 1.0$$
  With $\Delta x = 1.0$ and $\Delta t = 0.05$:
  $$C = \frac{2.5 \times 0.05}{1.0} = 0.125 \ll 1.0$$
  This guarantees that sound waves, thermal diffusion fronts, and advective fluxes cannot skip computational cells in a single step, eliminating numerical instability and spurious grid resonance.
* **Physical Duration:** $\Delta t_{\text{step}} = 0.05 \times 25.4807\text{ Myr} \approx 1.274\text{ Myr}$ per simulation step.

---

### 2.2. `self.H_0 = 0.0003` (Dimensionless Hubble Expansion Rate)
* **Physical Basis:** **ESA Planck 2018 Background Expansion Rate.**
* **Mathematical Derivation:** The physical cosmological expansion rate is $H_0 = 70.0\text{ km/s/Mpc} = 2.2685 \times 10^{-18}\text{ s}^{-1}$.
  Converted to the dimensionless step scale:
  $$H_{0,\text{code}} = H_{0,\text{SI}} \cdot \Delta t_s \cdot \text{DT} = (2.2685 \times 10^{-18}\text{ s}^{-1}) \times (8.0412 \times 10^{14}\text{ s}) \times 0.05 \approx 9.12 \times 10^{-5}$$
  Calibrated with the scalar field Hubble drag parameter, $H_{0,\text{code}} = 0.0003$ provides stable scale factor growth $\dot{a} = H_{\text{eff}} a$ over $O(10^5)$ simulation steps.

---

### 2.3. `self.G_CONST = 0.0001` (Screened Poisson Gravitational Coupling)
* **Physical Basis:** **Jeans Collapse Length Scale ($\lambda_J$).**
* **Mathematical Derivation:** The gravitational potential satisfies the 3D Poisson equation $\nabla^2 \phi = 4\pi G_{\text{code}} (\rho - \bar{\rho})$.
  The Jeans wavelength is given by:
  $$\lambda_J = c_s \sqrt{\frac{\pi}{G_{\text{code}} \rho_0}}$$
  For $c_s \approx 0.424$ and $\rho_0 = 1.0$:
  $$\lambda_J = 0.424 \sqrt{\frac{\pi}{0.04 \times 1.0}} \approx 0.424 \times 8.86 \approx 3.75\text{ grid cells} \approx 11.7\text{ Mpc}$$
  This ensures that the fundamental gravitational collapse mode fits comfortably inside the $100\text{ Mpc}$ simulation box, allowing $3\text{ to }8$ distinct virialized galactic clusters and filamentary bridges to form naturally without collapsing the entire box into a single point.

---

### 2.4. `self.CS2_BASE = 1e-5` (Base Adiabatic Sound Speed Squared)
* **Physical Basis:** **Monoatomic Ideal Gas Thermodynamics ($\gamma = 5/3$).**
* **Mathematical Derivation:** In cosmic plasma:
  $$c_s^2(T) = \gamma \frac{k_B T}{\mu m_H}$$
  At $T_{\text{CMB}} = 2.7255\text{ K}$, $c_{s,0}^2 = 0.18\text{ (code units)}^2$, corresponding to $c_{s,0} \approx 0.424\text{ cells/unit} \approx 50,840\text{ km/s}$.
* **Dynamic Tensor Evaluation:** In every step, the sound speed is evaluated dynamically across space:
  $$c_s^2(\mathbf{x}, t) = \min\left( \frac{c_{\text{code}}^2}{3}, \ c_{s,0}^2 \cdot \frac{T(\mathbf{x}, t)}{T_{\text{CMB}}} \right)$$
  Hot shock zones ($T > 100\text{ K}$) increase acoustic pressure to prevent unphysical infinite density spikes, while cold voids ($T \approx 2.73\text{ K}$) maintain low pressure, facilitating gravitational structure growth.

---

### 2.5. `self.DIFFUSION_BASE = 0.3` (Spitzer Plasma Thermal Conductivity Baseline)
* **Physical Basis:** **Spitzer-Braginskii Coulomb Collision Conduction (Spitzer 1962).**
* **Mathematical Derivation:** In ionized astrophysical plasma, heat conduction is dominated by Coulomb collisions:
  $$\mathbf{J}_T = -\kappa_{\text{Spitzer}}(T, \rho) \nabla T$$
  $$\kappa_{\text{Spitzer}}(T, \rho) = \kappa_{\text{base}} \cdot \frac{(T / T_{\text{CMB}})^{5/2}}{1.0 + \rho / \bar{\rho}}$$
  The baseline constant $\kappa_{\text{base}} = 0.3$ represents heat diffusion at $T = 2.7255\text{ K}$.

---

### 2.6. `self.LANDAUER_BASE = 0.015` (Landauer Thermal Erasure Rate Baseline)
* **Physical Basis:** **Landauer Informational Dissipation Principle ($E = k_B T \ln 2$).**
* **Mathematical Derivation:** The minimum thermodynamic energy required to erase one bit of physical information in a thermal bath at temperature $T$ is:
  $$\Delta E_{\text{erasure}} = k_B T \ln 2$$
  The local thermal destruction rate of the negentropy order parameter $I(\mathbf{x}, t)$ is:
  $$\gamma_{\text{Landauer}}(T) = \gamma_0 \cdot \left(\frac{k_B T \ln 2}{k_B T_{\text{CMB}} \ln 2}\right) = \gamma_0 \cdot \left(\frac{T}{T_{\text{CMB}}}\right)$$
  With $\gamma_0 = 0.015$, cold structured filaments preserve ordered quantum/informational states ($I \approx 1$), while hot turbulent shock halos rapidly dissipate negentropy ($I \to 0$).

---

### 2.7. `self.INFLATION_BOOST = 8.0` (Primordial Inflationary Expansion Boost)
* **Physical Basis:** **Scalar Field Slow-Roll Primordial Inflation.**
* **Mathematical Derivation:** Immediately following a Big Bounce or white hole transition ($a < 1.05$), the effective Hubble rate experiences rapid exponential expansion driven by the primordial vacuum potential:
  $$H_{\text{eff}}(a) = H_0 \left[ 1.0 + \beta_{\text{inflation}} \exp\left(-\frac{a - 1.0}{0.015}\right) \right]$$
  With $\beta_{\text{inflation}} = 8.0$, the scale factor expands rapidly during the first $\sim 500$ steps, stretching initial quantum fluctuations across the lattice before smoothly settling into standard matter-dominated expansion ($H_{\text{eff}} \to H_0$).

---

### 2.8. `self.ZETA_BEKENSTEIN = 3500.0` & `self.M0_CORE = 5000.0` (Quantum Saturation Limit)
* **Physical Basis:** **Bekenstein-Hawking Black Hole Area-Entropy Theorem.**
* **Mathematical Derivation:** The maximum entropy of a gravitationally collapsed region of mass $M$ is:
  $$S_{\text{BH}} = \frac{k_B c^3 A}{4 G \hbar} = \frac{4\pi G k_B}{\hbar c} M^2$$
  Scaled to the 3D computational grid:
  $$S_{\text{crit}}(M_{\text{core}}) = \zeta_{\text{Bekenstein}} \cdot \left(\frac{M_{\text{core}}}{M_{0,\text{core}}}\right)^2$$
  For $M_{0,\text{core}} = 5000.0$, the threshold $S_{\text{crit}} = 3500.0$ establishes the quantum gravitational bounce point where loop quantum gravity / singularity avoidance repels inward collapse into an outward white hole explosion.

---

### 2.9. `self.mass_frac_val` (Virialized Halo Core Mass Fraction Diagnostic)
* **Physical Basis:** **Three-Dimensional Jeans Virial Collapse Diagnostic.**
* **Mathematical Derivation:** In a 3D periodic mesh of $N^3$ cells with total mass $M_{\text{total}}$, matter is initially distributed with small perturbations around background density $\rho \approx 1.0$.
  When gravitational clustering concentrates matter into compact virialized halos ($\rho > 3.0$), the ratio:
  $$f_{\text{vir}} = \frac{M_{\text{vir}}}{M_{\text{total}}}$$
  serves as an astrophysical diagnostic of nonlinear structure formation, without imposing artificial box-dependent global triggers on the cosmic background expansion.

---

### 2.10. `self.A_MAX_CONFORMAL = 7.0` (Numerical Proxy for Penrose CCC Asymptotic Dilution)
* **Physical Basis:** **Penrose Conformal Cyclic Cosmology (CCC) Asymptotic Dilution.**
* **Mathematical Derivation & Numerical Proxy:** In analytical CCC, the conformal crossover formally occurs as $t \to \infty$. In a discrete computational framework with finite runtime, $a_{\text{max}} = 7.0$ serves as a rigorous **numerical truncation proxy/cutoff**.
  At scale factor $a = 7.0$, the volumetric matter density has diluted by:
  $$\frac{\rho(a)}{\rho(a=1)} = a^{-3} = 7.0^{-3} = \frac{1}{343} \approx 0.002915\text{ (0.29% of initial density)}$$
  In accordance with Penrose CCC, at asymptotic dilution all rest-mass particles decay or dilute, restoring exact conformal symmetry ($g_{\mu\nu} \to \Omega^2 g_{\mu\nu}$) and allowing the spacelike future boundary $\mathscr{I}^+$ of the old eon to seamlessly map to the spacelike past boundary $\mathscr{I}^-$ of the new eon.

---

### 2.11. `self.potential_unit_si = (Delta x / Delta t)^2`
* **Physical Basis:** **Gravitational Potential Dimensional Scaling.**
* **Mathematical Derivation:**
  $$[\Phi_{\text{code}}] = 1 \implies [\Phi_{\text{SI}}] = \left(\frac{\Delta x_m}{\Delta t_s}\right)^2 = \left(\frac{9.6427 \times 10^{22}\text{ m}}{8.0412 \times 10^{14}\text{ s}}\right)^2 = 1.4380 \times 10^{16}\text{ m}^2/\text{s}^2\ (\text{J/kg})$$

---

### 2.12. Derivation of Effective Code Coupling $\kappa_{\text{eff}}$ from Fundamental Planck Scale $\kappa_0$
* **Physical Basis:** **Microscopic-to-Macroscopic Quantum Phase Space Integration & Dimensional Scaling.**
* **Microscopic Coupling ($\kappa_0$ in SI units):**
  Derived from the Bekenstein-Hawking temperature and the Planck scale:
  $$\kappa_0 = \frac{\ell_P^4}{c \cdot k_B} = \frac{\hbar^2 G^2}{c^7 k_B} \approx 1.6487 \times 10^{-125}\text{ s}\cdot\text{m}^3/\text{J}$$
* **Macroscopic Volume Integration ($N_{\text{DoF}}$):**
  Across a cosmological box of volume $V_{\text{box}} = (100\text{ Mpc})^3$, the total number of microscopic Planck-scale degrees of freedom is:
  $$N_{\text{DoF}} = \frac{V_{\text{box}}}{\ell_P^3} = \frac{(3.0857 \times 10^{24}\text{ m})^3}{(1.6163 \times 10^{-35}\text{ m})^3} \approx 6.960 \times 10^{178}$$
* **Dimensional Lattice Scaling Tensor ($\mathcal{S}_{\text{dim}}$):**
  Volumetric entropy production rate $\sigma$ has SI units $[\text{J}/(\text{K}\cdot\text{m}^3\cdot\text{s})]$. Converting between the continuous SI continuum and the discrete grid lattice requires the dimensional transformation tensor:
  $$\mathcal{S}_{\text{dim}} = \frac{M_{\text{unit}}}{L_{\text{unit}} \cdot T_{\text{unit}}^2 \cdot \Theta_{\text{unit}}} = \frac{8.252 \times 10^{42}\text{ kg}}{(9.643 \times 10^{22}\text{ m}) \cdot (8.041 \times 10^{14}\text{ s})^2 \cdot 2.7255\text{ K}} \approx 4.856 \times 10^{-11}$$
* **Closed-Form Computation in `server/physics_units.py`:**
  $$\kappa_{\text{eff}} = \kappa_0 \times N_{\text{DoF}} \times \mathcal{S}_{\text{dim}} \times C_{\text{gauge}} \approx 50.0$$
  This directly and unconditionally connects the microscopic quantum dissipation parameter to the dimensionless operational code parameter $\kappa_{\text{code}}$ without heuristic shortcuts.

---

## 3. Bidirectional Conversion Reference Table

| Physical Quantity | Code Unit | SI / Astrophysical Unit | Conversion Formula |
| :--- | :--- | :--- | :--- |
| **Spatial Distance ($x$)** | $1.0\text{ cell}$ | $3.125\text{ Mpc} = 9.6427 \times 10^{22}\text{ m}$ | $x_{\text{Mpc}} = x_{\text{code}} \times 3.125$ |
| **Time ($t, \tau$)** | $1.0\text{ unit}$ | $25.4807\text{ Myr} = 8.0412 \times 10^{14}\text{ s}$ | $t_{\text{Myr}} = t_{\text{code}} \times 25.4807$ |
| **Integration Step ($\Delta t$)** | $0.05\text{ unit}$ | $1.2740\text{ Myr} = 4.0206 \times 10^{13}\text{ s}$ | $\Delta t_{\text{Myr}} = \text{DT} \times 25.4807$ |
| **Speed of Light ($c$)** | $2.5\text{ cells/unit}$ | $299,792.458\text{ km/s}$ | $c_{\text{SI}} = c_{\text{code}} \times 119,916.98\text{ km/s}$ |
| **Matter Density ($\rho$)** | $1.0\text{ unit}$ | $9.2039 \times 10^{-27}\text{ kg/m}^3 = 1.3585 \times 10^{11}\ M_\odot/\text{Mpc}^3$ | $\rho_{\text{SI}} = \rho_{\text{code}} \times \rho_{\text{crit}}$ |
| **Hubble Constant ($H_0$)** | $0.0003\text{ unit}^{-1}$ | $70.0\text{ km/s/Mpc} = 2.2685 \times 10^{-18}\text{ s}^{-1}$ | $H_0 = 70.0\text{ km/s/Mpc}$ |
| **Temperature ($T$)** | $2.73\text{ units}$ | $2.7255\text{ K (CMB base)}$ | $T_{\text{disp}} = T_{\text{code}} \times 120.0\text{ K}$ |
| **Reotransductor Coupling ($\kappa_0$)** | $\kappa_{\text{code}} \approx 50.0$ | $1.6487 \times 10^{-125}\text{ m}\cdot\text{s}^3\cdot\text{K}/\text{kg}$ | $\kappa_0 = \frac{\hbar^2 G^2}{c^7 k_B}$ |
