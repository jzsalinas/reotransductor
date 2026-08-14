# Rheotransductor (Reotransductor)

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![NumPy](https://img.shields.io/badge/NumPy-2.0+-013243.svg?logo=numpy)](https://numpy.org/)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-3.8+-11557c.svg)](https://matplotlib.org/)

A physical and numerical simulation exploring the **Active Present Rheotransducer** (*Reotransductor del Presente Activo*) — modeling the emergence of local **thermal time** ($\tau$) driven by Onsager entropy production, non-equilibrium dissipative heat flows, finite energy reservoirs, and Landauer's informational limit.

---

## 📸 Simulation Dashboard

![Rheotransductor Simulation Dashboard](assets/preview.png)

---

## 🌌 Theoretical Framework

The simulation integrates non-equilibrium thermodynamics, information theory, and the thermal time hypothesis across a 2D spatial continuum:

### 1. Thermal Diffusion (Fourier Law)
Heat propagates through the medium according to Fourier's heat equation without boundary wrap-around:
$$\frac{\partial T}{\partial t} = k \, \nabla^2 T$$
where $k$ is the thermal conductivity (`DIFFUSION_COEFF = 0.5`) and outer boundaries are fixed to cosmic background temperature ($T_{\text{vacuum}} = 2.73\text{ K}$).

### 2. Onsager Dissipation & Entropy Production
The thermodynamic force driving dissipation is the gradient of inverse temperature:
$$\mathbf{X} = \nabla\left(\frac{1}{T}\right)$$
Together with the heat flux $\mathbf{J} = -k \nabla T$, the local volumetric entropy production rate $\sigma$ is strictly non-negative:
$$\sigma = \mathbf{J} \cdot \nabla\left(\frac{1}{T}\right) \ge 0$$

### 3. The Rheotransducer: Emergent Thermal Time ($\tau$)
Local proper time does not tick uniformly; it emerges from dissipative irreversibility:
$$\frac{d\tau}{dt} = \kappa \cdot \max\left(0,\, \mathbf{J} \cdot \nabla\left(\frac{1}{T}\right)\right)$$
where $\kappa$ is the Rheotransducer coupling constant (`KAPPA = 50.0`). In thermal equilibrium ($\nabla T = 0$), $d\tau/dt = 0$ — **time ceases to flow**.

### 4. Negentropy Indices & Landauer's Limit ($I$)
Low-entropy structures (biological organisms, memory registers, coherent states) require continuous negentropic consumption to counteract thermal noise and Landauer decay:
$$\frac{dI}{dt} = \alpha \|\mathbf{J}\| - \beta T - \gamma_{\text{Landauer}}$$
* **Sustenance**: Fed by local dissipative energy flux $\alpha \|\mathbf{J}\|$.
* **Thermal Noise**: Thermal agitation $\beta T$ disrupting order.
* **Landauer Decay**: Intrinsic erasure cost $\gamma_{\text{Landauer}} = 0.02$.

### 5. Finite Energy Reservoirs & Asymptotic Cooling
Boilers are modeled as finite sensible energy reservoirs ($E = C \cdot \Delta T$). Temperature drops smoothly proportional to remaining fuel:
$$T_{\text{boiler}}(t) = T_{\text{vacuum}} + (T_{\text{max}} - T_{\text{vacuum}}) \times \left(\frac{\text{Fuel}(t)}{\text{Fuel}_{\text{init}}}\right)$$
In accordance with Newton's law of cooling, heat dissipation naturally exhibits asymptotic exponential decay.

---

## 🎛️ Dashboard Overview

The interactive interface is divided into 4 real-time synchronized monitors:

| Panel | Metric | Description |
|---|---|---|
| **Top-Left** | **Temperature ($T$)** | Real-time thermal grid with dynamic fuel status and live boiler temperatures in Kelvins. |
| **Top-Right** | **Rheotransducer Velocity ($d\tau/dt$)** | Local entropy production density representing the current velocity of time's flow. |
| **Bottom-Left** | **Accumulated Time ($\tau$)** | Total emergent thermal seconds accumulated per spatial coordinate. |
| **Bottom-Right** | **Negentropy Index ($I$)** | Structural integrity of bounded low-entropy islands (**A**, **B**, **C**). |

### Monitored Negentropy Islands:
* **`A` (Cyan)**: Positioned close to Boiler 1 ($1000\text{ K}$) — Sustained by high dissipative energy flux.
* **`B` (Lime)**: Positioned near Boiler 2 ($800\text{ K}$) — Moderately sustained with steady integrity.
* **`C` (Red)**: Isolated in cold vacuum ($2.73\text{ K}$) — Decays monotonically due to absence of sustaining heat flux.

---

## ⚡ Quick Start

### Prerequisites
* **Python 3.10+**
* `pip` or [`uv`](https://github.com/astral-sh/uv)

### 1. Clone the repository
```bash
git clone https://github.com/your-username/reotransductor.git
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

### 3. Run the simulation
```bash
python simulador_reotransductor.py
```

---

## 🎮 Interactive Controls

* **Acceleration Slider**: Adjust the simulation speed dynamically from **1x** up to **50x** physics steps per frame using the slider at the bottom of the dashboard.
* **Matplotlib Navigation Toolbar**: Pan, zoom into specific grid regions, or export snapshot figures directly.

---

## 📁 Repository Structure

```
reotransductor/
├── assets/
│   └── preview.png               # High-resolution simulation snapshot
├── .gitignore                    # Python & virtual environment ignore rules
├── LICENSE                       # MIT License
├── pyproject.toml                # Project metadata & packaging configuration
├── requirements.txt              # Production dependencies
├── README.md                     # Project documentation
└── simulador_reotransductor.py    # Main simulation engine and UI
```

---

## 📜 License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

## 👤 Author

* **José Salinas** ([@jzsalinas](https://github.com/jzsalinas)) - *Initial concept, thermodynamics formulation & simulation engine.*
