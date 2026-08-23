"""
Multi-Resolution Numerical Convergence and Mass Conservation Study.
Evaluates simulation fidelity, mass drift, execution throughput, and observable invariance
across spatial grid resolutions (N = 16^3, 32^3, 64^3, 128^3, 256^3).
"""

import os
import sys
import time
import json
import numpy as np

# Ensure project root is in Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from server.engine import CosmologicalEngine
from observational.halo_analyzer import HaloRadialProfileAnalyzer
from observational.bao_analyzer import BAOSpatialCorrelationAnalyzer
from observational.cmb_analyzer import CMBSphericalHarmonicsAnalyzer
from observational.hubble_tension import HubbleTensionAnalyzer


def run_convergence_benchmark(resolutions=[16, 32, 64, 128, 256], n_steps=25, box_size_mpc=500.0):
    """
    Executes standardized multi-resolution hydrodynamic integration benchmark.
    Returns structured results dictionary.
    """
    results = []
    print("=" * 80)
    print("  🔬 REOTRANSDUCTOR: MULTI-RESOLUTION NUMERICAL CONVERGENCE STUDY")
    print("=" * 80)
    print(f"• Physical Box Size: {box_size_mpc} Mpc | Benchmark Steps: {n_steps}")
    print("-" * 80)
    
    for N in resolutions:
        print(f"\n[Benchmarking Grid N = {N}^3 ({N**3:,} voxels)]...")
        dx_mpc = box_size_mpc / N
        
        # Initialize engine in non-resuming clean mode
        engine = CosmologicalEngine(grid_size=N, box_size_mpc=box_size_mpc, auto_resume=False)
        
        # Initial mass
        m_initial = float(np.sum(engine.to_cpu(engine.rho)) * (dx_mpc**3))
        
        # Timed execution loop
        t0 = time.perf_counter()
        for _ in range(n_steps):
            engine.step()
        t1 = time.perf_counter()
        wall_time_per_step_ms = ((t1 - t0) / n_steps) * 1000.0
        
        # Final mass and mass drift
        m_final = float(np.sum(engine.to_cpu(engine.rho)) * (dx_mpc**3))
        mass_drift_rel = abs(m_final - m_initial) / max(1e-10, m_initial)
        
        # Compute primary physical observables on current state
        rho_cpu = engine.to_cpu(engine.rho).astype(np.float64)
        tau_cpu = engine.to_cpu(engine.tau).astype(np.float64)
        t_cpu = engine.to_cpu(engine.T).astype(np.float64)
        
        # 1. CMB Low-ell power
        cmb_analyzer = CMBSphericalHarmonicsAnalyzer(n_theta=32, n_phi=64, ell_max=8)
        c_map = cmb_analyzer.extract_celestial_sphere_from_grid(tau_cpu, T_3d=t_cpu, rho_3d=rho_cpu)
        cmb_metrics = cmb_analyzer.compute_angular_power_spectrum(c_map)
        c2_c3_ratio = float(cmb_metrics["ratio_C2_C3"])
        
        # 2. Halo Core Inner Slope
        halo_analyzer = HaloRadialProfileAnalyzer(grid_size=N, box_size_mpc=box_size_mpc, n_shells=min(16, N // 2))
        halo_metrics = halo_analyzer.evaluate_halo_diagnostics(rho_cpu)
        gamma_0 = halo_metrics["inner_slope_gamma0"]
        
        # 3. BAO Correlation Peak
        bao_analyzer = BAOSpatialCorrelationAnalyzer(grid_size=N, box_size_mpc=box_size_mpc, n_bins=16)
        bao_metrics = bao_analyzer.evaluate_cosmological_fields(rho_cpu, scale_factor=engine.scale_factor)
        r_bao = bao_metrics["bao_peak_radius_mpc"]
        
        # 4. Environmental Hubble Shift
        hubble_analyzer = HubbleTensionAnalyzer()
        hubble_metrics = hubble_analyzer.compute_3d_environmental_h0_field(rho_cpu, tau_cpu)
        delta_h0 = hubble_metrics["h0_cluster"] - hubble_metrics["h0_void"]
        
        # Estimate GPU VRAM footprint in MB
        # Total state tensors: rho, vx, vy, vz, T, Phi, I, tau (8 tensors * 4 bytes) + scratch
        vram_mb = (N**3 * 4 * 12) / (1024 * 1024)
        
        entry = {
            "grid_size": N,
            "total_voxels": N**3,
            "cell_size_mpc": round(dx_mpc, 3),
            "wall_time_ms": round(wall_time_per_step_ms, 2),
            "estimated_vram_mb": round(vram_mb, 1),
            "mass_drift": float(mass_drift_rel),
            "gamma_0": round(gamma_0, 3),
            "r_bao_mpc": round(r_bao, 2),
            "c2_c3_ratio": round(c2_c3_ratio, 3),
            "delta_h0": round(delta_h0, 3)
        }
        results.append(entry)
        
        print(f"  • dx = {dx_mpc:.3f} Mpc | Step: {wall_time_per_step_ms:.2f} ms | VRAM: {vram_mb:.1f} MB")
        print(f"  • Mass Drift: {mass_drift_rel:.2e} | Halo Slope γ₀: {gamma_0:.3f} | BAO r_peak: {r_bao:.2f} Mpc")
        print(f"  • CMB C₂/C₃: {c2_c3_ratio:.3f} | Environmental ΔH₀: +{delta_h0:.2f} km/s/Mpc")
        
    return results


def plot_convergence_results(results, output_path: str = "assets/numerical_convergence_study.png"):
    """Generates 4-panel publication-grade white paper figure showing numerical convergence."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    grids = [r["grid_size"] for r in results]
    dx = [r["cell_size_mpc"] for r in results]
    mass_drift = [r["mass_drift"] for r in results]
    gamma_0 = [r["gamma_0"] for r in results]
    r_bao = [r["r_bao_mpc"] for r in results]
    c2_c3 = [r["c2_c3_ratio"] for r in results]
    time_ms = [r["wall_time_ms"] for r in results]
    
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.5), dpi=300)
    fig.patch.set_facecolor('white')
    
    # -------------------------------------------------------------------------
    # Panel 1: Mass Conservation & Numerical Drift
    # -------------------------------------------------------------------------
    ax1 = axes[0, 0]
    ax1.set_facecolor('white')
    ax1.plot(grids, mass_drift, marker='o', color='#1d4ed8', linewidth=2.2, markersize=7, label=r'Mass Drift $|\Delta M / M_0|$')
    ax1.axhline(y=1e-4, color='#dc2626', linestyle='--', linewidth=1.5, label='Conservation Bound ($10^{-4}$)')
    ax1.set_yscale('log')
    ax1.set_xlabel('Grid Resolution $N$', color='#0f172a', fontsize=11, fontweight='bold')
    ax1.set_ylabel(r'Fractional Mass Drift $|\Delta M / M_0|$', color='#0f172a', fontsize=11, fontweight='bold')
    ax1.set_title('(A) Global Mass Conservation', color='#0f172a', fontsize=12, fontweight='bold')
    ax1.tick_params(colors='#334155', labelsize=9.5)
    ax1.grid(True, linestyle='--', alpha=0.5, color='#e2e8f0')
    ax1.legend(loc='upper right', facecolor='#ffffff', edgecolor='#cbd5e1', fontsize=8.5)
    
    # -------------------------------------------------------------------------
    # Panel 2: Physical Observable Invariance (Halo Core Slope gamma_0)
    # -------------------------------------------------------------------------
    ax2 = axes[0, 1]
    ax2.set_facecolor('white')
    ax2.plot(grids, gamma_0, marker='s', color='#15803d', linewidth=2.2, markersize=7, label=r'Inner Slope $\gamma_0(N)$')
    ax2.axhline(y=-0.12, color='#0284c7', linestyle='--', linewidth=1.5, label=r'SPARC DDO 154 Observed ($\gamma = -0.12$)')
    ax2.axhline(y=-1.00, color='#b91c1c', linestyle=':', linewidth=1.5, label=r'NFW Cuspy Model ($\gamma = -1.00$)')
    ax2.set_xlabel('Grid Resolution $N$', color='#0f172a', fontsize=11, fontweight='bold')
    ax2.set_ylabel(r'Inner Logarithmic Slope $\gamma_0 = d\ln\rho / d\ln r$', color='#0f172a', fontsize=11, fontweight='bold')
    ax2.set_title('(B) Dark Matter Halo Core Invariance', color='#0f172a', fontsize=12, fontweight='bold')
    ax2.set_ylim(-1.1, 0.1)
    ax2.tick_params(colors='#334155', labelsize=9.5)
    ax2.grid(True, linestyle='--', alpha=0.5, color='#e2e8f0')
    ax2.legend(loc='lower right', facecolor='#ffffff', edgecolor='#cbd5e1', fontsize=8.5)
    
    # -------------------------------------------------------------------------
    # Panel 3: BAO Acoustic Scale Invariance
    # -------------------------------------------------------------------------
    ax3 = axes[1, 0]
    ax3.set_facecolor('white')
    ax3.plot(grids, r_bao, marker='^', color='#ea580c', linewidth=2.2, markersize=7, label=r'Simulation Peak $r_{\mathrm{BAO}}(N)$')
    ax3.axhline(y=101.4, color='#0284c7', linestyle='--', linewidth=1.5, label=r'DESI DR2 / BOSS Fiducial ($101.4\ \mathrm{Mpc}/h$)')
    ax3.set_xlabel('Grid Resolution $N$', color='#0f172a', fontsize=11, fontweight='bold')
    ax3.set_ylabel(r'Acoustic Peak Scale $r_{\mathrm{BAO}}$ (Mpc)', color='#0f172a', fontsize=11, fontweight='bold')
    ax3.set_title('(C) BAO Acoustic Scale Convergence', color='#0f172a', fontsize=12, fontweight='bold')
    ax3.tick_params(colors='#334155', labelsize=9.5)
    ax3.grid(True, linestyle='--', alpha=0.5, color='#e2e8f0')
    ax3.legend(loc='lower right', facecolor='#ffffff', edgecolor='#cbd5e1', fontsize=8.5)
    
    # -------------------------------------------------------------------------
    # Panel 4: Computational Scaling (Wall Time vs Grid Size)
    # -------------------------------------------------------------------------
    ax4 = axes[1, 1]
    ax4.set_facecolor('white')
    ax4.plot(grids, time_ms, marker='d', color='#7c3aed', linewidth=2.2, markersize=7, label='GPU Wall Time per Step')
    ax4.set_yscale('log')
    ax4.set_xlabel('Grid Resolution $N$', color='#0f172a', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Execution Time per Step (ms)', color='#0f172a', fontsize=11, fontweight='bold')
    ax4.set_title('(D) Computational Throughput Scaling', color='#0f172a', fontsize=12, fontweight='bold')
    ax4.tick_params(colors='#334155', labelsize=9.5)
    ax4.grid(True, linestyle='--', alpha=0.5, color='#e2e8f0')
    ax4.legend(loc='upper left', facecolor='#ffffff', edgecolor='#cbd5e1', fontsize=8.5)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()
    print(f"\n✅ Convergence Study Figure saved to: {output_path}")


if __name__ == "__main__":
    results = run_convergence_benchmark(resolutions=[16, 32, 64, 128, 256], n_steps=20, box_size_mpc=500.0)
    plot_convergence_results(results)
    
    # Output formatted Markdown table
    print("\n" + "=" * 80)
    print("  📋 NUMERICAL CONVERGENCE & SCALING SUMMARY TABLE")
    print("=" * 80)
    print("| Resolution $N$ | Total Cells | $\\Delta x$ (Mpc) | Mass Drift $|\\Delta M/M_0|$ | Halo Slope $\\gamma_0$ | BAO $r_{\\mathrm{BAO}}$ | Step Time (ms) | Peak VRAM |")
    print("|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|")
    for r in results:
        print(f"| ${r['grid_size']}^3$ | {r['total_voxels']:,} | {r['cell_size_mpc']:.2f} | {r['mass_drift']:.2e} | {r['gamma_0']:.3f} | {r['r_bao_mpc']:.1f} Mpc | {r['wall_time_ms']:.1f} ms | {r['estimated_vram_mb']:.1f} MB |")
    print("=" * 80)
