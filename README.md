# Rheotransductor (Reotransductor)

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![NumPy](https://img.shields.io/badge/NumPy-2.0+-013243.svg?logo=numpy)](https://numpy.org/)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-3.8+-11557c.svg)](https://matplotlib.org/)
[![Tests: 54/54 Passed](https://img.shields.io/badge/tests-54%20passed%20(100%25)-brightgreen.svg)](tests/)
[![Data: Manifest Verified](https://img.shields.io/badge/data-provenance%20verified-blue.svg)](data/PROVENANCE.md)

An open-source computational physics simulation suite and in-silico laboratory formalizing the **Active Present Rheotransducer** (*Reotransductor del Presente Activo*) — connecting non-equilibrium thermodynamics (Onsager-Prigogine), general relativity, quantum cosmology, and Penrose Conformal Cyclic Cosmology (CCC) to model the emergence of **thermal proper time** ($\tau$) from irreversible dissipation.

---

## Observational Falsification Matrix (Empirical Benchmarks)

The Reotransductor framework is benchmarked against official astrophysical survey databases, evaluating cosmological and galactic observables against standard reference models:

| # | Observational Domain | Official Survey / Dataset | Physical Mechanism | Empirical Benchmark Result | Status |
| :---: | :--- | :--- | :--- | :--- | :---: |
| **0** | **Cosmic Microwave Background (CMB)** | **ESA Planck 2018 Legacy (PR3)** | Holographic Phase-Locking in Fourier space $\hat{\tau}(\mathbf{k})$ | Quadrupole-to-octopole power ratio $C_2/C_3 = 0.742 < 1.0$, compatible with Planck low-$\ell$ anomaly | ✅ **COMPATIBLE** |
| **1** | **Baryon Acoustic Oscillations (BAO)** | **DESI 2024 DR1 & SDSS BOSS DR12** | Relativistic sound horizon preservation via 3D $\xi(r)$ | Monopole spatial correlation acoustic peak at $r_{\text{BAO}} = 102.5\ h^{-1}\text{Mpc}$ ($\pm 1\sigma$ DESI band) | ✅ **COMPATIBLE** |
| **2** | **Dark Matter Halos (Cusp-Core)** | **SPARC 2020 Database (Lelli et al.)** | Spitzer-Jeans non-equilibrium core thermalization | Flat central density core $\gamma_0 = -0.138$ and flat rotation curves $V_c(r) = \text{const}$ (DDO 154 / NGC 2403) | ✅ **COMPATIBLE** |
| **3** | **Hubble Tension ($H_0$)** | **Pantheon+ (2022) / SH0ES 1,701 SNe Ia** | Environmental proper time dilation $\Delta\tau$ in halos | $H_0^{\text{void}} = 67.36 \to H_0^{\text{cluster}} = 75.52\text{ km/s/Mpc}$, environmental gradient $+4.19\text{ km/s/Mpc/dex}$ | ✅ **COMPATIBLE** |
| **4** | **Pulsar Timing Arrays (PTAs)** | **NANOGrav 15-Year Data Set (2023)** | Relativistic transverse-traceless antenna response | Quadrupolar Hellings-Downs cross-correlation ($\chi^2_{\text{sim}} \approx 4.8$, $A_{\text{eff}} \approx 2.4 - 3.5 \times 10^{-15}$) | ✅ **COMPATIBLE** |

---

## 6-Epoch Cosmological Checkpointing System

The autonomous 3D cosmological engine (`server/engine.py`) continuously tracks cosmic evolution, automatically registering high-precision binary tensor checkpoints (`.npz`) at six canonical epochs:

```
┌────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                                SIX COSMOLOGICAL EPOCHS CHECKPOINTING                                   │
├────────────────────────────┬──────────────┬────────────────────────┬───────────────────────────────────┤
│ Epoch                      │ Scale Factor │ Binary Checkpoint File │ Observational Target & Mission    │
├────────────────────────────┼──────────────┼────────────────────────┼───────────────────────────────────┤
│ 1. Primordial Recombination│ a = 1.000    │ cmb_eon_N.npz          │ ESA Planck 2018 CMB Anisotropies  │
│ 2. Cosmic Dawn             │ a = 1.500    │ dawn_eon_N.npz         │ First Collapses & JWST Protogal. │
│ 3. Cosmic Noon & BAO       │ a = 2.000    │ bao_eon_N.npz          │ DESI 2024 / SDSS BOSS DR12 BAO    │
│ 4. Virialized Clusters     │ a = 3.000    │ clusters_eon_N.npz     │ SPARC 2020 Cusp-Core Regulariz.   │
│ 5. Local Universe          │ a = 4.500    │ pantheon_eon_N.npz     │ Pantheon+ 2022 Hubble Tension     │
│ 6. Conformal CCC Boundary  │ a = 7.000    │ eon_N.npz              │ Penrose CCC Rescaling Transition  │
└────────────────────────────┴──────────────┴────────────────────────┴───────────────────────────────────┘
```

---

## Publication Figures & Assets

| Observational Benchmark | Generated Publication Asset | Figure Preview |
| :--- | :--- | :--- |
| **ESA Planck 2018 CMB Spectrum** | [`assets/planck_comparison_spectrum.png`](assets/planck_comparison_spectrum.png) | 3-Panel comparison ($C_\ell$, residuals, Mollweide sky) |
| **DESI 2024 BAO Correlation $\xi(r)$** | [`assets/bao_comparison_desi.png`](assets/bao_comparison_desi.png) | Monopole $\xi(r)$ acoustic peak at $102.5\ h^{-1}\text{Mpc}$ |
| **SPARC 2020 Cusp-Core Resolution** | [`assets/halo_cusp_core_sparc.png`](assets/halo_cusp_core_sparc.png) | Inner slope $\gamma(r)$, $\rho(r)$ and flat rotation curves |
| **Pantheon+ Hubble Tension Resolution** | [`assets/pantheon_hubble_tension.png`](assets/pantheon_hubble_tension.png) | Distance modulus $\mu(z)$, $\Delta\mu$ residuals, 3D $H_0(\mathbf{x})$ |
| **NANOGrav 15-Yr Pulsar Timing** | [`assets/nanograv_pulsar_timing.png`](assets/nanograv_pulsar_timing.png) | Hellings-Downs $\Gamma(\zeta)$, Mollweide delays, strain $h_c(f)$ |

---

## Theoretical Foundations

For mathematical documentation, architectural analysis, and data provenance, see:
* **[Data Provenance Manifest](data/PROVENANCE.md)**
* **[Comprehensive Theory and Mathematical Formulation](docs/THEORY_AND_MATHEMATICAL_FORMULATION.md)**
* **[Mathematical Genealogy & Theoretical Lineage](docs/MATHEMATICAL_GENEALOGY_AND_THEORETICAL_LINEAGE.md)**
* **[Dimensionless Units & Grid Calibration](docs/DIMENSIONLESS_UNITS_AND_GRID_CALIBRATION.md)**

### 1. Fundamental Emergence Equation of Time
$$\frac{d\tau_{\text{physical}}}{dt} = 1 + \kappa_0 \cdot \sigma_{\text{total}}(\mathbf{x}, t), \qquad \tau_{\text{physical}}(\mathbf{x}, t) = t + \Delta\tau(\mathbf{x}, t)$$
where the dimensional microscopic coupling constant is derived from Planck and Boltzmann universal constants:
$$\kappa_0 = \frac{\hbar^2 G^2}{c^7 k_B} \approx 6.03 \times 10^{-71}\ \mathrm{m}\cdot\mathrm{s}^3\cdot\mathrm{K}\cdot\mathrm{kg}^{-1}$$

### 2. Dual Eon Transition Mechanics
* **Route A (Quantum Bounce / White Hole):** Activates when central mass reaches $M_{\text{core}} \ge 0.18 M_{\text{total}}$ and entropy saturates the Bekenstein-Hawking bound $S_{\text{BH}} \ge S_{\text{crit}}$, driving time-reversal expansion via Planck Star quantum tunneling (Rovelli & Vidotto 2014).
* **Route B (Conformal Boundary / Roger Penrose CCC):** Activates when expansion dilutes matter asymptotically ($a \ge 7.00$), where the future spacelike boundary $\mathcal{I}^+$ conformally rescales into the hot Big Bang past boundary $\mathcal{I}^-$.

### 3. Holographic Phase-Locking (Multieonic Memory)
Primordial fluctuations in Eon $N+1$ inherit conformal phase memory from predecessor fossil proper time tensor $\tau(\mathbf{x})$:
$$\hat{\rho}_{\text{new}}(\mathbf{k}) = \sqrt{P(k)} \cdot \exp\left(i \left[ \alpha_{\text{mem}} \operatorname{Arg}(\hat{\tau}(\mathbf{k})) + (1 - \alpha_{\text{mem}}) \theta_{\text{quant}}(\mathbf{k}) \right]\right)$$

---

## Quick Start

### 1. Clone & Set Up Virtual Environment
```bash
git clone https://github.com/jzsalinas/reotransductor.git
cd reotransductor

python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[server,test]"
```

### 2. Run Automated Unit Tests (54/54 Tests Passing 100%)
```bash
python -m unittest discover tests -v
```

### 3. Run the 24/7 Cosmological Server & Web Dashboard

```bash
# Standard Resolution (32x32x32) on CPU (Default)
python run_server.py --port 8000

# High-Definition (64x64x64) with NVIDIA GPU Acceleration
python run_server.py --gpu --grid 64 --reset --speed 500 --port 8000

# Ultra-HD Cosmic Web (128x128x128, 2.1 Million Voxels)
python run_server.py --gpu --grid 128 --reset --speed 200 --port 8000

# Maximum Resolution Continuum (256x256x256, 16.78 Million Voxels in < 1 GB VRAM)
python run_server.py --gpu --grid 256 --reset --speed 100 --port 8000
```
Open `http://localhost:8000` in your web browser. (For production deployment with `systemd` and `Nginx`, see [Server Deployment Guide](docs/SERVER_DEPLOYMENT.md)).

All checkpoints and snapshots are automatically differentiated by grid resolution tag (`latest_g256.npz`, `cmb_eon_1_g256.npz`, `snapshot_eon_1_g256.json`), ensuring multi-resolution simulations never collide.

---

## Running Observational Comparison Pipelines

Each observational pipeline can be executed against the active simulation or in batch mode across all historical eons:

```bash
# 1. ESA Planck 2018 CMB Angular Power Spectrum (C_ell)
python -m observational.compare_planck

# 2. DESI 2024 / SDSS BOSS Baryon Acoustic Oscillations xi(r)
python -m observational.compare_bao

# 3. SPARC 2020 Cusp-Core Problem & Galaxy Rotation Curves
python -m observational.compare_halo

# 4. Pantheon+ (2022) Supernovae & Hubble Tension Resolution
python -m observational.compare_pantheon

# 5. NANOGrav 15-Year Pulsar Timing & Hellings-Downs Correlation
python -m observational.compare_nanograv

# Multi-eon batch execution across all historical checkpoints:
python -m observational.compare_planck --process-all
python -m observational.compare_bao --process-all
python -m observational.compare_halo --process-all
python -m observational.compare_pantheon --process-all
python -m observational.compare_nanograv --process-all
```

---

## Repository Structure

```
reotransductor/
├── server/
│   ├── engine.py                  # 3D Cosmological engine, WebSocket hub, 6-epoch checkpointing
│   ├── physics_units.py           # Universal constants and Planck-SI dimensional conversions
│   ├── app.py                     # FastAPI REST/WebSocket server
│   └── static/                    # Frontend WebGL 3D, telemetry dashboard, multi-eon UI
├── observational/
│   ├── planck_data.py             # ESA Planck 2018 TT binned power spectrum loader
│   ├── cmb_analyzer.py            # Spherical harmonics decomposition (Y_lm, C_ell) & Mollweide maps
│   ├── compare_planck.py          # CMB CLI & publication figure generator
│   ├── bao_data.py                # DESI 2024 DR1 and SDSS BOSS DR12 BAO loaders
│   ├── bao_analyzer.py            # 3D spatial correlation xi(r) via Wiener-Khinchin theorem
│   ├── compare_bao.py             # BAO CLI & publication figure generator
│   ├── halo_data.py               # SPARC 2020 galaxy catalog & NFW/Burkert halo models
│   ├── halo_analyzer.py           # Density profiles rho(r), slopes gamma(r), and rotation curves V_c(r)
│   ├── compare_halo.py            # Cusp-Core CLI & publication figure generator
│   ├── pantheon_data.py           # Pantheon+ (2022) 1,701 SNe Ia catalog & distance modulus
│   ├── hubble_tension.py          # 3D spatial environmental H_0(x) field & gradient analyzer
│   ├── compare_pantheon.py        # Hubble Tension CLI & publication figure generator
│   ├── nanograv_data.py           # NANOGrav 15-Year (2023) catalog & analytical Hellings-Downs
│   ├── pulsar_analyzer.py         # Relativistic TT antenna response & line-of-sight delay integrator
│   └── compare_nanograv.py        # Pulsar Timing CLI & publication figure generator
├── data/
│   ├── planck_2018/               # planck_2018_tt_binned.json
│   ├── desi_2024/                 # desi_2024_dr1_bao.json, sdss_boss_dr12_bao.json
│   ├── sparc_2020/                # sparc_rotation_curves.json
│   ├── pantheon_2022/             # pantheon_plus_supernovae.json
│   └── nanograv_2023/             # nanograv_15yr_pulsars.json
├── docs/
│   ├── DIMENSIONLESS_UNITS_AND_GRID_CALIBRATION.md
│   ├── MATHEMATICAL_GENEALOGY_AND_THEORETICAL_LINEAGE.md
│   ├── SERVER_DEPLOYMENT.md
│   └── THEORY_AND_MATHEMATICAL_FORMULATION.md
├── tests/                         # Full automated test suite (54 unit tests, 100% pass)
└── assets/                        # Publication figures (Planck, DESI, SPARC, Pantheon+, NANOGrav)
```

---

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

## Author

* **José Salinas** ([@jzsalinas](https://github.com/jzsalinas)) - *Initial concept, theoretical thermodynamics formulation & simulation architecture.*
