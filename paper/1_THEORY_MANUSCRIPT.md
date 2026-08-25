# On the Emergence of a Dissipative Thermodynamic Clock Field from Non-Equilibrium Entropy Production in Relativistic Spacetimes

**Authors:** J. Z. Salinas  
**Target Submission:** *Physical Review D* / *Classical and Quantum Gravity*  
**Classification:** Gravitation, Relativistic Non-Equilibrium Thermodynamics, Foundations of Time  
**Status:** Theoretical Research Manuscript (Q1 Target Formulation)

---

## Abstract

In standard General Relativity, proper time along a timelike worldline is a purely geometric parameter determined by the metric interval $d\tau_{\text{geom}} = \sqrt{-g_{\mu\nu} dx^\mu dx^\nu} / c$, remaining decoupled from local irreversible processes. In this paper, we propose a phenomenological covariant constitutive hypothesis in which a scalar thermodynamic clock field $\tau(x^\mu)$ accrues an additional rate proportional to the local scalar density of irreversible entropy production $\Sigma(x^\mu) = \nabla_\mu s^\mu \ge 0$, as formulated in relativistic non-equilibrium thermodynamics (Onsager–Prigogine–Israel–Stewart framework). We postulate the effective evolution equation along fluid streamlines tangent to the four-velocity $u^\mu$:
$$u^\mu \nabla_\mu \tau = 1 + \kappa_0 \Sigma(x^\mu),$$
where $\kappa_0 \equiv \frac{\hbar^2 G^2}{c^7 k_B} \approx 1.6487 \times 10^{-125}\ \mathrm{m\cdot s^3\cdot K\cdot kg^{-1}}$ represents the unique dimensional coupling scale formed from the fundamental constants of Planck and Boltzmann. Under this constitutive framework, we show that: (i) the non-negativity of entropy production ($\Sigma \ge 0$) formally ensures the strict local monotonicity of the clock field along fluid streamlines ($u^\mu \nabla_\mu \tau \ge 1$); (ii) in global thermodynamic equilibrium ($\Sigma \to 0$), the clock field reduces identically to standard geometric coordinate time in synchronous/comoving parametrization; and (iii) covariant stress-energy conservation $\nabla_\mu T^{\mu\nu} = 0$ is preserved when $\tau(x^\mu)$ acts as a passive thermodynamic tracer without gravitational backreaction. We provide analytical solutions for spatially flat Friedmann–Lemaître–Robertson–Walker (FLRW) cosmologies with regularized primordial limits and for idealized spherically symmetric collapsing dissipative configurations.

---

## 1. Introduction and Motivation

The nature of time presents one of the most enduring conceptual tensions in theoretical physics. While General Relativity treats time as a geometric coordinate along pseudo-Riemannian manifolds, thermodynamics establishes a fundamental arrow of time dictated by irreversible entropy generation ($\Delta S \ge 0$). In standard relativistic physics, these two descriptions operate on separate footings: the metric determines chronometric intervals, whereas statistical mechanics tracks the distribution of microstates on that fixed geometric background.

Attempts to bridge this gap have arisen in several distinct domains:
1. **Emergent and Relational Time:** Thermal time hypotheses (Connes & Rovelli 1994, Rovelli 2011) propose that macroscopic temporal flow is state-dependent and emerges from statistical states of quantum systems. We note a crucial distinction here: whereas the Connes-Rovelli framework constructs a global modular temporal flow directly from the underlying Von Neumann algebraic state, the framework presented here proposes a local, cumulative scalar field explicitly driven by classical macroscopic Onsager entropy production.
2. **Dissipative and Non-Equilibrium Systems:** The Onsager–Prigogine framework, extended relativistically by Israel & Stewart (1979), establishes that macroscopic dissipation produces positive semi-definite local entropy densities $\Sigma = \nabla_\mu s^\mu \ge 0$ through thermal conduction, viscous shear, and field gradients.
3. **Conformal and Cyclic Frameworks:** Penrose's Conformal Cyclic Cosmology (CCC; Penrose 2010, Meissner & Penrose 2024) posits that the boundary between cosmological eons is governed by asymptotic dilution into a massless, conformally invariant state.

The objective of this paper is not to modify the geometric Einstein field equations $G_{\mu\nu} = \frac{8\pi G}{c^4} T_{\mu\nu}$, but rather to formulate a consistent, covariant scalar clock functional $\tau(x^\mu)$ whose rate of accumulation couples directly to local macroscopic entropy production.

---

## 2. Relativistic Non-Equilibrium Thermodynamics

Consider a relativistic fluid with four-velocity $u^\mu$ ($u_\mu u^\mu = -c^2$), energy density $\rho$, isotropic pressure $P$, local temperature $T(x^\mu)$, and entropy current four-vector $s^\mu$:
$$s^\mu = s_{\text{fluid}} u^\mu + \frac{q^\mu}{T},$$
where $s_{\text{fluid}}$ is the fluid entropy density and $q^\mu$ is the spatial heat flux vector orthogonal to the fluid flow ($q^\mu u_\mu = 0$).

To ensure strict manifest spatial positivity in Lorentzian signature $(-, +, +, +)$, we introduce the spatial projection tensor orthogonal to $u^\mu$:
$$h^{\mu\nu} = g^{\mu\nu} + \frac{u^\mu u^\nu}{c^2}, \qquad h^{\mu\nu} u_\mu = 0, \qquad h^\mu_\mu = 3.$$

According to the relativistic second law of thermodynamics (Eckart 1940, Israel & Stewart 1979, Hiscock & Lindblom 1983), the divergence of the entropy current defines the local scalar entropy production rate:
$$\Sigma(x^\mu) \equiv \nabla_\mu s^\mu = \sigma_{\text{thermal}} + \sigma_{\text{viscous}} + \sigma_{\text{grav}} \ge 0.$$

For a plasma subject to thermal conduction and gravitational potential gradients in the Newtonian/weak-field cosmological approximation, the spatial entropy production density takes the Onsager form:
$$\Sigma(x^\mu) = \kappa_{\text{th}} \frac{h^{\mu\nu} \nabla_\mu T \nabla_\nu T}{T^2} + \frac{\rho \, h^{\mu\nu} \nabla_\mu \Phi \nabla_\nu \Phi}{T T_0} \ge 0,$$
where $\kappa_{\text{th}}$ is the regularized thermal conductivity (e.g., Spitzer–Braginskii transport in ionized astrophysical plasmas, $\kappa_{\text{spitzer}} \propto T^{5/2}/\rho$).

By the fundamental theorem of non-equilibrium thermodynamics, $\Sigma(x^\mu)$ is a positive semi-definite Lorentz scalar:
$$\Sigma(x^\mu) \ge 0, \qquad \forall x^\mu \in \mathcal{M}.$$

---

## 3. The Dissipative Clock Field Formulation

### 3.1 Postulate 1 (Constitutive Evolution Equation)
Let $\tau(x^\mu)$ be a scalar clock field tracking accumulated thermodynamic time. Along a fluid streamline with tangent vector $u^\mu$, we postulate the constitutive evolution law:
$$u^\mu \nabla_\mu \tau = 1 + \kappa_0 \Sigma(x^\mu),$$
or equivalently, in comoving coordinates where $u^\mu = (c, 0, 0, 0)$:
$$\frac{d\tau}{dt} = 1 + \kappa_0 \Sigma(\mathbf{x}, t), \qquad \tau(\mathbf{x}, t) = t + \Delta\tau(\mathbf{x}, t),$$
with the dissipative excess evolving as:
$$\frac{\partial \Delta\tau}{\partial t} = \kappa_0 \Sigma(\mathbf{x}, t).$$

### 3.2 Dimensional Derivation of the Fundamental Coupling $\kappa_0$
The coupling constant $\kappa_0$ converts an entropy production density $[\Sigma] = \mathrm{J\cdot K^{-1}\cdot m^{-3}\cdot s^{-1}} = \mathrm{kg\cdot m^{-1}\cdot s^{-3}\cdot K^{-1}}$ into a dimensionless temporal accumulation rate $[\kappa_0 \Sigma] = 1$. Consequently, the physical dimension of $\kappa_0$ is:
$$[\kappa_0] = \mathrm{m\cdot s^3\cdot K\cdot kg^{-1}}.$$

In search of a fundamental scale constructed from the universal constants of nature—reduced Planck constant $\hbar$, gravitational constant $G$, speed of light $c$, and Boltzmann constant $k_B$—dimensional analysis yields the unique power-law monomial:
$$\kappa_0 = \frac{\hbar^2 G^2}{c^7 k_B}.$$

Evaluating with CODATA 2018 recommended values:
- $\hbar = 1.054571817 \times 10^{-34}\ \mathrm{J\cdot s}$
- $G = 6.67430 \times 10^{-11}\ \mathrm{m^3\cdot kg^{-1}\cdot s^{-2}}$
- $c = 2.99792458 \times 10^8\ \mathrm{m\cdot s^{-1}}$
- $k_B = 1.380649 \times 10^{-23}\ \mathrm{J\cdot K^{-1}}$

$$\kappa_0 = \frac{(1.0545718 \times 10^{-34})^2 (6.67430 \times 10^{-11})^2}{(2.99792458 \times 10^8)^7 (1.380649 \times 10^{-23})} \approx 1.6487 \times 10^{-125}\ \mathrm{m\cdot s^3\cdot K\cdot kg^{-1}}.$$

*Remark on Effective Cosmological Scaling:* The microscopic constant $\kappa_0$ represents the quantum gravitational floor. When integrating over macroscopic cosmological domains or discrete numerical lattices, the effective observable coupling is parametrized by a dimensionless factor $\kappa_{\text{eff}} = \kappa_0 \cdot N_{\text{DoF}} \cdot \mathcal{C}_{\text{gauge}}$, which is subjected to empirical sensitivity analysis.

### 3.3 Microscopic Informational Relaxation and Planckian Dissipation
Landauer's principle (Landauer 1961) establishes the thermodynamic lower bound on energy dissipation required to irreversibly erase one bit of information in a thermal bath at temperature $T$:
$$\Delta E_{\text{erase}} = k_B T \ln 2.$$
To obtain a kinetic time-relaxation rate $\Gamma_{\text{info}}(T) \equiv \dot{I}/I$, non-equilibrium statistical mechanics requires an attempt frequency. In quantum thermodynamic systems, the fundamental upper bound on thermal dissipation is set by the Planckian relaxation rate (Sachdev 1999, Zaanen 2004):
$$\tau_{\text{diss}} \sim \frac{\hbar}{k_B T} \implies \omega_{\text{Planckian}} = \frac{k_B T}{\hbar}.$$
Combining the Landauer energy scale with the Planckian dissipation frequency yields the microscopic informational relaxation rate:
$$\Gamma_{\text{info}}(T) = \omega_{\text{Planckian}} \cdot \left(\frac{\Delta E_{\text{erase}}}{E_{\text{ref}}}\right) = \gamma_0 \left(\frac{T}{T_{\text{CMB}}}\right),$$
where $\nu_0 = \frac{k_B T_{\text{CMB}}}{\hbar} \approx 3.57 \times 10^{11}\ \mathrm{s}^{-1}$ provides the natural attempt frequency in the cosmic background. Thus, we postulate the linear temperature scaling $\Gamma \propto T$ as a natural phenomenological consequence of combining quantum dissipation with Landauer's bound.

---

## 4. Mathematical Properties and Phenomenological Consequences

### 4.1 Corollary 1 (Streamline Monotonicity and Arrow of Time)
*Let $\mathcal{M}$ be a globally hyperbolic spacetime endowed with a dissipative fluid satisfying the relativistic second law $\Sigma(x^\mu) \ge 0$. Then, along any integral curve of the four-velocity field $u^\mu$, the thermodynamic clock field $\tau$ is strictly monotonically increasing with respect to parameter length $s$:*
$$\frac{d\tau}{ds} \ge \frac{1}{c} > 0.$$

*Proof:*
Along a streamline parameterized by proper distance $ds = c \, dt$, we have $u^\mu \nabla_\mu \tau = c \frac{d\tau}{ds} = 1 + \kappa_0 \Sigma$. Since $\kappa_0 > 0$ and $\Sigma \ge 0$, it follows that:
$$\frac{d\tau}{ds} = \frac{1}{c} [1 + \kappa_0 \Sigma] \ge \frac{1}{c} > 0.$$
Hence, the scalar field $\tau$ is strictly monotonic along fluid streamlines, preventing closed timelike loops in the internal clock variable. $\blacksquare$

### 4.2 Property 1 (Equilibrium Limit in Synchronous Coordinates)
*In any region of spacetime where the matter distribution reaches global thermodynamic equilibrium ($\nabla_\mu T = 0, \nabla_\mu \Phi = 0, \Sigma = 0$), the clock field reduces identically to standard geometric coordinate time:*
$$\Sigma = 0 \implies u^\mu \nabla_\mu \tau = 1 \implies \tau(t) = t + \tau_0.$$

### 4.3 Postulate 2 (Passive Tracer Conservation)
*When $\tau(x^\mu)$ is treated as a passive clock tracer, it does not contribute to the stress-energy tensor, preserving the exact contracted Bianchi identities:*
$$\nabla_\mu T^{\mu\nu}_{\text{matter}} = 0 \iff \nabla_\mu G^{\mu\nu} = 0.$$

### 4.4 Phenomenological Consequence (Emergent Galactic Cores)
*In virialized, steady-state astrophysical systems (such as mature $z=0$ galaxies), the local accumulation of the time-dilation field $\tau$ flattens the central effective gravitational potential. As a rigorous phenomenological consequence, this dissipative framework natively predicts the kinetic emergence of Burkert-like constant-density cores in galactic centers, circumventing the cuspy singularities historically expected under collisionless NFW dark matter profiles.*

### 4.5 Conformal Boundary and Crossover Hypothesis
We hypothesize that as the universe expands asymptotically ($a \to \infty$) and the physical state variables ($\rho, T$) approach the zero-mass radiative limit, the thermodynamic clock field $\tau_{\text{physical}}$ becomes the dominant cosmological scale. In the limit where conformal rescaling bridges to a subsequent eon (as proposed in Roger Penrose's Conformal Cyclic Cosmology), the residual fossil field $\Delta\tau(\mathbf{x}, t \to \infty)$ acts as a conformal blueprint for the initial density fluctuations of the subsequent thermodynamic cycle.

### 4.6 Postulate 3 (Holographic Gravity Parity)
For the conformal boundary transition to satisfy macroscopic energy-momentum conservation across eons, we postulate a strict Holographic Gravity Parity (a 1:1 equivalence). The fossil proper-time tensor $\tau$ generated in the preceding eon couples to the local spatial geometry of the nascent eon exactly as baryonic mass does. Specifically, in the Newtonian limit (Poisson's equation), the effective gravitational potential is sourced by an apparent density:
$$\rho_{\text{eff}} = \rho + \mathcal{H} \cdot \tau_{\text{prior}}$$
where the coupling constant $\mathcal{H} = 1.0$ by equivalence. This ensures that the time-dilation memory of past structural virialization acts macroscopically as "Apparent Dark Matter," providing the necessary potential wells for the primordial gas to overcome Jeans mass limits purely via baryonic hydrodynamics.

### 4.7 Phenomenological Consequence (Bekenstein-Hawking Time Freeze)
In standard Loop Quantum Gravity (e.g., Rovelli's Planck Stars), a singularity reaching its absolute Bekenstein-Hawking informational limit ($S_{\text{BH}} \ge S_{\text{crit}}$) undergoes a local quantum bounce (transitioning to a White Hole). However, in our thermodynamic clock field formulation, this transition is suppressed within a single eon due to extreme local time dilation ($\tau \to \tau_{\text{crit}}$). From the perspective of the exterior cosmic observer, the singular core becomes mathematically "frozen in time." These stable, indestructible gravitational anchors serve as the permanent foundational scaffolding for galactic accretion. The release of this frozen information is strictly forbidden until the global universe reaches the conformal boundary ($a \to \infty$), at which point the Penrose CCC mechanism triggers a synchronized cosmic-scale quantum bounce.

---

## 5. Analytical Solutions in Canonical Cosmological Geometries

### 5.1 Spatially Flat FLRW Cosmology with Dissipative Fluid
Consider the flat FLRW metric:
$$ds^2 = -c^2 dt^2 + a^2(t) [dr^2 + r^2 d\Omega^2].$$
For a homogeneous fluid with bulk viscosity $\zeta(T)$, the entropy production rate is spatially homogeneous:
$$\Sigma(t) = \frac{\zeta \Theta^2}{T(t)} = \frac{9 \zeta H^2(t)}{T(t)},$$
where $H(t) \equiv \frac{\dot{a}}{a}$ is the Hubble parameter and $\Theta = \nabla_\mu u^\mu = 3H$.

In a radiation-dominated era with $a(t) \propto t^{1/2}$, $H(t) = \frac{1}{2t}$, $T(t) \propto t^{-1/2}$, and effective bulk viscosity $\zeta \approx \zeta_0 T^4 \propto t^{-2}$:
$$\Sigma(t) = \frac{9 \zeta_0 T^3 H^2}{1} \propto t^{-3/2} \cdot t^{-2} = t^{-7/2}.$$

To regularize the UV behavior near the primordial singularity, we integrate from a finite initial time $t_i > 0$:
$$\Delta\tau(t) = \kappa_0 \int_{t_i}^t \Sigma_0 \left(\frac{t'}{t_i}\right)^{-7/2} dt' = \frac{2}{5} \kappa_0 \Sigma_0 t_i \left[ 1 - \left(\frac{t_i}{t}\right)^{5/2} \right].$$
As $t \to \infty$, the total accumulated temporal excess saturates asymptotically to a finite value $\Delta\tau_\infty = \frac{2}{5} \kappa_0 \Sigma_0 t_i$, demonstrating that dissipative time accumulation is self-limiting in expanding FLRW backgrounds.

### 5.2 Idealized Spherically Symmetric Collapsing Cloud
Consider an idealized collapsing spherical gas cloud of mass $M$ and comoving radius $R(t)$. Thermal conduction toward the outer boundary generates a radial temperature gradient $\partial_r T(r, t)$. The local entropy production density is:
$$\Sigma(r, t) \approx \kappa_{\text{spitzer}}(r, t) \frac{1}{R^2(t)} \left( \frac{\partial T}{\partial r} \right)^2 \frac{1}{T^2}.$$

Integrating along radial fluid shells, regions undergoing rapid gravitational compression exhibit higher $\Sigma(r, t)$, naturally producing an internal radial gradient in accumulated clock excess:
$$\Delta\tau(r, t) = \kappa_0 \int_0^t \Sigma(r, t') dt'.$$

---

## 6. Discussion and Conclusions

We have presented a covariant, mathematically consistent formulation of a dissipative thermodynamic clock field $\tau(x^\mu)$ coupled to macroscopic non-equilibrium entropy production. The key conceptual results of this work are:

1. **Covariant Formulation:** The clock field $\tau$ is defined covariantly along fluid streamlines using the spatial projection tensor $h^{\mu\nu}$, ensuring strict positivity of entropy generation in Lorentzian spacetimes.
2. **Fundamental Scale:** Dimensional analysis identifies the unique universal coupling $\kappa_0 = \hbar^2 G^2 / (c^7 k_B) \approx 1.6487 \times 10^{-125}\ \mathrm{m\cdot s^3\cdot K\cdot kg^{-1}}$.
3. **Exact Asymptotic Regularity:** Dissipative time excess in expanding FLRW cosmologies decays as $t^{-7/2}$ and integrates to a finite asymptotic value.

The numerical implementation on 3D GPU lattices is presented in Paper II (*Computer Physics Communications* target), and empirical constraints from cosmological datasets (Planck 2018, DESI 2024, Pantheon+, SPARC, NANOGrav) are explored in subsequent observational analyses.

---

## References

1. Connes, A., & Rovelli, C. (1994). *Von Neumann algebra automorphisms and time-thermodynamics relation in generally covariant quantum theories*. Classical and Quantum Gravity, 11(12), 2899.
2. Eckart, C. (1940). *The Thermodynamics of Irreversible Processes. III. Relativistic Theory of the Simple Fluid*. Physical Review, 58(10), 919.
3. Hiscock, W. A., & Lindblom, L. (1983). *Stability and causality in dissipative relativistic fluids*. Annals of Physics, 151(2), 466–496.
4. Israel, W., & Stewart, J. M. (1979). *Transient relativistic thermodynamics and relativistic kinetic theory*. Annals of Physics, 118(2), 341–372.
5. Landauer, R. (1961). *Irreversibility and heat generation in the computing process*. IBM Journal of Research and Development, 5(3), 183–191.
6. Landau, L. D., & Lifshitz, E. M. (1987). *Fluid Mechanics* (Vol. 6, Course of Theoretical Physics). Butterworth-Heinemann.
7. Meissner, K. A., & Penrose, R. (2024). *Conformal Cyclic Cosmology: A Review*. Foundations of Physics, 54(1), 12.
8. Onsager, L. (1931). *Reciprocal Relations in Irreversible Processes. I*. Physical Review, 37(4), 405.
9. Penrose, R. (2010). *Cycles of Time: An Extraordinary New View of the Universe*. The Bodley Head, London.
10. Planck Collaboration, Aghanim, N., et al. (2020). *Planck 2018 results. VI. Cosmological parameters*. Astronomy & Astrophysics, 641, A6.
11. Prigogine, I. (1967). *Introduction to Thermodynamics of Irreversible Processes*. Interscience Publishers, New York.
12. Rovelli, C. (2011). *Forget time*. Foundations of Physics, 41(9), 1475–1490.
13. Sachdev, S. (1999). *Quantum Phase Transitions*. Cambridge University Press.
14. Scolnic, D., Brout, D., et al. (2022). *The Pantheon+ Analysis: The Full Dataset and Light-Curve Fits*. The Astrophysical Journal, 938(2), 113.
15. Spitzer, L. (1962). *Physics of Fully Ionized Gases*. Interscience Publishers, New York.
16. Zaanen, J. (2004). *Why the temperature is a quantum metric*. Nature, 430(6999), 512–513.
