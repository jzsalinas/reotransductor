# Rheotransductor (Reotransductor)

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![NumPy](https://img.shields.io/badge/NumPy-2.0+-013243.svg?logo=numpy)](https://numpy.org/)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-3.8+-11557c.svg)](https://matplotlib.org/)

A suite of physical and numerical simulations exploring the **Active Present Rheotransducer** (*Reotransductor del Presente Activo*) — modeling the emergence of local **thermal time** ($\tau$) driven by Onsager entropy production, non-equilibrium dissipative heat flows, finite energy reservoirs, and Landauer's informational limit.

---

## 🌟 Available Simulations

| Simulation | File | Description |
|---|---|---|
| **Core Rheotransducer** | [`simulador_reotransductor.py`](simulador_reotransductor.py) | Fundamental non-equilibrium thermodynamic simulation with finite boilers, Landauer negentropy decay, and Onsager dissipation. |
| **Cosmological Rheotransducer** | [`simulador_cosmologico.py`](simulador_cosmologico.py) | Cosmological model replacing the unphysical global FLRW time with **cellular emergent proper time** ($\tau_i$) across a collapsing cosmic web and freezing voids. |

---

## 📸 Simulation Dashboards

### 1. Cosmological Rheotransducer Dashboard
![Cosmological Rheotransducer Dashboard](assets/preview_cosmology.png)

### 2. Core Thermodynamic Dashboard
![Rheotransductor Simulation Dashboard](assets/preview.png)

---

## 🌌 Theoretical Framework

The project integrates non-equilibrium thermodynamics, gravitational collapse, and information theory across a 2D spatial continuum:

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

### 4. Negentropy Indices & Landauer's Limit ($I$)
Low-entropy structures (galaxies, biological organisms, coherent states) require continuous negentropic consumption to counteract thermal noise and Landauer decay:
$$\frac{dI}{dt} = \alpha \|\mathbf{J}_{\text{total}}\| - \beta T - \gamma_{\text{Landauer}}$$
* **`A` (Cyan - Central Supercluster)**: Sustained by intense gravitational accretion and dissipative negentropy.
* **`B` (Lime - Filament Node)**: Moderately sustained by filamentary inflow.
* **`C` (Red - Deep Void Dwarf)**: Decays monotonically due to absence of sustaining dissipative energy flux.

---

## ⚡ Quick Start

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

**To run the Cosmological Model:**
```bash
python simulador_cosmologico.py
```

**To run the Core Thermodynamic Simulation:**
```bash
python simulador_reotransductor.py
```

---

## 🎮 Interactive Controls

* **Acceleration Slider**: Dynamically adjust execution from **1x** up to **50x** physics steps per frame.
* **Matplotlib Navigation Toolbar**: Pan, zoom, and export vector or raster snapshots.

---

## 📁 Repository Structure

```
reotransductor/
├── assets/
│   ├── preview.png                 # Core simulation snapshot
│   └── preview_cosmology.png       # Cosmological simulation snapshot
├── .gitignore                      # Python & virtual environment ignore rules
├── LICENSE                         # MIT License
├── pyproject.toml                  # Project metadata & packaging configuration
├── requirements.txt                # Production dependencies
├── README.md                       # Documentation & mathematical foundations
├── simulador_cosmologico.py        # Cosmological model with cellular time emergence
└── simulador_reotransductor.py      # Core thermodynamic Rheotransducer simulation
```

---

## 📜 License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

## 👤 Author

* **José Salinas** ([@jzsalinas](https://github.com/jzsalinas)) - *Initial concept, thermodynamics formulation & simulation engine.*
