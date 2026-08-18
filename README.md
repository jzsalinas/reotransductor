# Rheotransductor (Reotransductor)

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![NumPy](https://img.shields.io/badge/NumPy-2.0+-013243.svg?logo=numpy)](https://numpy.org/)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-3.8+-11557c.svg)](https://matplotlib.org/)

A suite of physical and numerical simulations exploring the **Active Present Rheotransducer** (*Reotransductor del Presente Activo*) — modeling the emergence of local **thermal time** ($\tau$) driven by Onsager entropy production, non-equilibrium dissipative heat flows, finite energy reservoirs, and Landauer's informational limit.

---

## Available Simulations

| Simulation | File | Description |
|---|---|---|
| **24/7 Cosmological Server & Web Dashboard** | [`run_server.py`](run_server.py) | Autonomous headless multi-core engine with automated checkpointing, multi-eon history recorder, and real-time WebSocket Web Dashboard (Dell PowerEdge R820 / Nginx ready). |
| **3D Cosmological Model (GPU CUDA)** | [`simulador_cosmologico_3d_gpu.py`](simulador_cosmologico_3d_gpu.py) | Desktop GPU-accelerated 3D simulation with CuPy in-VRAM integration, Planck CMB sky projection, and Matplotlib GUI. |
| **3D Cosmological Model (CPU NumPy)** | [`simulador_cosmologico_3d.py`](simulador_cosmologico_3d.py) | Desktop 3D simulation with pure NumPy, 3D Poisson gravity, Planck CMB sky projection, and Matplotlib GUI. |
| **2D Cosmological Model** | [`simulador_cosmologico.py`](simulador_cosmologico.py) | 2D cosmological model with cellular time emergence, Bekenstein entropy saturation, White Hole bounce, and fossil memory. |
| **Core Thermodynamic Rheotransducer** | [`simulador_reotransductor.py`](simulador_reotransductor.py) | Fundamental non-equilibrium thermodynamic simulation with finite boilers, Landauer negentropy decay, and Onsager dissipation. |

---

## Simulation Dashboards

### 1. 3D Cosmological Model & CMB Sky Projection (Mollweide)
![3D Cosmological Rheotransducer Dashboard](assets/preview_cosmology_3d.png)

### 2. 2D Cosmological Rheotransducer Dashboard
![Cosmological Rheotransducer Dashboard](assets/preview_cosmology.png)

### 3. Core Thermodynamic Dashboard
![Rheotransductor Simulation Dashboard](assets/preview.png)

---

## Theoretical Framework & Mathematical Formulation

For the complete, rigorous mathematical derivations, differential equations, and code implementation mapping, see the dedicated document:
👉 **[Comprehensive Theory and Mathematical Formulation](docs/THEORY_AND_MATHEMATICAL_FORMULATION.md)**

The project integrates non-equilibrium thermodynamics, gravitational collapse, and information theory across spatial continua:

### 1. Thermal Diffusion & Fluid Dynamics
* **Thermal Conduction**: Heat flows according to Fourier's equation:
  $$\frac{\partial T}{\partial t} = k \, \nabla^2 T$$
* **Cosmological Matter Advection**: In the cosmological model, matter infall is driven by the gravitational potential $\Phi$ solved exactly via 2D Fast Fourier Transform ($\nabla^2 \Phi = 4\pi G (\rho - \bar{\rho})$):
  $$\mathbf{v} = -\mu \nabla \Phi, \quad \frac{\partial \rho}{\partial t} = -\nabla \cdot (\rho \mathbf{v}) + D_\rho \nabla^2 \rho$$

### 2. Onsager Dissipation & Entropy Production
The thermodynamic driving forces are the inverse temperature gradient and gravitational flux:
$$\mathbf{X}_T = \nabla\left(\frac{1}{T}\right), \quad \sigma = \mathbf{J}_T \cdot \nabla\left(\frac{1}{T}\right) + \frac{\rho \|\nabla \Phi\|^2}{T} \ge 0$$
where $\mathbf{J}_T = -k \nabla T$ is the heat flux.

### 3. The Rheotransducer: Emergent Cellular Proper Time ($\tau$)
Classical cosmology (e.g., Friedman-Lemaître-Robertson-Walker metrics) relies on an idealized, global cosmic clock. In reality, **proper time is local and emergent**:
$$\frac{d\tau_i}{dt} = \kappa \cdot \sigma_i(x, y)$$
* **Collapsing Clusters & Superclusters**: High matter accretion and thermal dissipation produce high entropy production $\sigma \implies$ time ticks vigorously ($d\tau/dt \gg 0$).
* **Cosmic Voids**: In near-uniform underdense voids ($\nabla T \to 0, \nabla \Phi \to 0$), $\sigma \to 0 \implies$ **local time practically freezes**.

### 4. Dynamic Informational Field & Landauer's Limit ($I$)
Low-entropy structures (coherent states, complex matter, biological systems) require continuous negentropic consumption to counteract thermal noise and Landauer decay:
* **Core Thermodynamic Model**: Models discrete bounded index islands ($I_A, I_B, I_C$) fed by local thermal flux $\mathbf{J}_T$.
* **Cosmological Model (Continuous Field $I(\mathbf{r}, t)$)**: In accordance with Schrödinger, Prigogine, and Landauer, order is represented as a **continuous dynamic scalar field** that advects with matter velocity $\mathbf{v}$ and self-organizes in dissipative filamentary nodes:
  $$\frac{\partial I}{\partial t} + \nabla \cdot (I \mathbf{v}) = D_I \nabla^2 I + \alpha \sigma \left(\frac{\rho}{\bar{\rho}}\right) - \beta T - \gamma_{\text{Landauer}} I$$
  When galaxies and filaments merge, the informational field naturally coalesces into unified negentropic superclusters, while expanding cosmic voids smoothly decay toward informational erasure ($I \to 0$).

### 5. Dual Cosmological Eon Transitions: Bekenstein Singularity & Penrose CCC Conformal Boundary
In standard general relativity, gravitational collapse leads to non-physical mathematical singularities ($\rho \to \infty$). The Reotransductor cosmology implements two complementary, physically grounded transition mechanisms:
* **Route A — Loop Quantum Gravity Singularity Bounce (Bekenstein Saturation)**: Grounded in **Loop Quantum Gravity (LQG/LQC)** (Ashtekar, Rovelli, Vidotto, Christodoulou), a collapsed black hole core ($\rho > 1.0$) saturates when its accumulated internal informational entropy exceeds the **Bekenstein-Hawking bound** ($S_{\text{max}} \propto M_{\text{BH}}^2$). Upon saturation, it undergoes quantum tunneling into an explosive **White Hole Blast**.
* **Route B — Penrose Conformal Cyclic Cosmology (CCC Heat-Death Crossover)**: In accordance with **Sir Roger Penrose's Conformal Cyclic Cosmology (Nobel Prize 2020)**, if cosmological expansion disperses matter before a singular black hole forms ($a \to a_{\text{max}} \ge 7.0$), the universe dilutes into asymptotic thermal heat death ($T \to 2.73\text{ K}, \nabla\Phi \to 0$). In this scale-invariant conformal regime, the cold future conformal boundary ($\mathcal{I}^+$) seamlessly rescales into the hot Big Bang past boundary ($\mathcal{I}^-$) of Eon $N+1$.
* **Conformal Memory Carrier ($\tau$)**: While scale factor $a$ and temperature $T$ undergo quantum reheating, the **Rheotransducer proper time field $\tau(\mathbf{x})$ is monotonically preserved across eons**, acting as the conformal memory tensor that seeds future structure formation.
* **Cosmic Inflation & Causal Speed Limit ($c$)**: Each new eon initiates with a brief **quantum inflationary super-expansion** ($a < 1.05$) that homogeneously stretches primordial perturbation modes across cosmological scales, strictly obeying relativistic causal advection ($\|\mathbf{v}\| \le c$).

---

## Quick Start

### Prerequisites
* **Python 3.10+**
* `pip` or [`uv`](https://github.com/astral-sh/uv)

### 1. Clone the repository
```bash
git clone https://github.com/jzsalinas/reotransductor.git
cd reotransductor
```

### 2. Set up virtual environment & install dependencies

Using standard `venv` and `pip`:
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Or using `uv` (fast):
```bash
uv venv .venv
source .venv/bin/activate
uv pip install -r requirements.txt
```

### 3. Run the simulations

**To run the 24/7 Cosmological Server & Live Web Dashboard (Dell PowerEdge R820 / Headless / Nginx):**
```bash
python run_server.py --port 8000
```
Open `http://localhost:8000` in your web browser. For full production deployment with `systemd` and `Nginx` reverse proxy, see 👉 **[Server Deployment Guide](docs/SERVER_DEPLOYMENT.md)**.

**To run the 3D Cosmological Model on GPU (NVIDIA CUDA / CuPy - Ultra-Fast):**
```bash
python simulador_cosmologico_3d_gpu.py
```

**To run the 3D Cosmological Model on CPU (pure NumPy):**
```bash
python simulador_cosmologico_3d.py
```

**To run the 2D Cosmological Model:**
```bash
python simulador_cosmologico.py
```

**To run the Core Thermodynamic Simulation:**
```bash
python simulador_reotransductor.py
```

---

## Interactive Controls

* **Web Dashboard**: Real-time WebSocket streaming, speed slider (1x - 500x), pause/resume, manual checkpoint triggers, and multi-eon history log.
* **Desktop Desktop Toolbar**: Pan, zoom, and export vector/raster snapshots directly from Matplotlib.

---

## Repository Structure

```
reotransductor/
├── assets/
│   ├── preview.png                 # Core simulation snapshot
│   ├── preview_cosmology.png       # 2D Cosmological simulation snapshot
│   └── preview_cosmology_3d.png    # 3D Cosmological simulation snapshot
├── docs/
│   ├── SERVER_DEPLOYMENT.md        # 24/7 Server deployment guide (Dell R820, systemd, Nginx)
│   └── THEORY_AND_MATHEMATICAL_FORMULATION.md # Complete physical foundations and equations
├── server/
│   ├── __init__.py
│   ├── engine.py                   # Autonomous headless physics engine with auto-checkpointing
│   ├── app.py                      # FastAPI server with WebSocket real-time streaming hub
│   └── static/                     # Web Dashboard frontend (HTML5, CSS3, ES6 Canvas/WebGL)
├── run_server.py                   # 24/7 Server entrypoint
├── simulador_cosmologico_3d_gpu.py # Desktop 3D GPU CUDA simulator
├── simulador_cosmologico_3d.py     # Desktop 3D CPU NumPy simulator
├── simulador_cosmologico.py        # Desktop 2D Cosmological simulator
├── simulador_reotransductor.py     # Core thermodynamic simulator
├── pyproject.toml                  # Project metadata
├── requirements.txt                # Production dependencies
└── README.md                       # Documentation & guide
```

---

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

## Author

* **José Salinas** ([@jzsalinas](https://github.com/jzsalinas)) - *Initial concept, thermodynamics formulation & simulation engine.*
