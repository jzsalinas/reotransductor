"""Pre-registered control-only validation protocol."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from .config import LaboratoryConfig
from .diagnostics import core_radius, logarithmic_slope, relative_error
from .hydrodynamics import primitive
from .solver import CuspCoreSolver


CONTROL_RESOLUTIONS = (512, 1024, 2048)
CONTROL_TIMES = (0.25, 0.5, 1.0, 2.0)
SLOPE_SAMPLE_RADII = (0.25, 0.5, 1.0, 2.0, 5.0, 9.5)


@dataclass(frozen=True)
class ControlMeasurement:
    cells: int
    time: float
    max_velocity: float
    density_l1_change: float
    max_slope_change: float
    mass_error: float
    total_energy_error: float
    core_radius: float | None
    slope_samples: dict[str, float]
    positivity_limiter_activation_frequency: float
    passed: bool


def validate_control_resolution(
    cells: int,
    times: tuple[float, ...] = CONTROL_TIMES,
    config: LaboratoryConfig | None = None,
) -> list[ControlMeasurement]:
    """Run K_hat=0 at fixed dimensionless physical output times."""
    cfg = config or LaboratoryConfig()
    solver = CuspCoreSolver(cells, conductivity_hat=0.0, config=cfg)
    initial = solver.snapshot()
    rho_initial, _, _ = primitive(initial.state, cfg.gamma)
    slope_initial = logarithmic_slope(solver.grid.centers, rho_initial)
    resolved = solver.grid.centers >= cfg.minimum_resolved_cells * solver.grid.dr
    measurements: list[ControlMeasurement] = []
    for output_time in times:
        snapshot = solver.evolve_to(output_time)
        rho, velocity, _ = primitive(snapshot.state, cfg.gamma)
        slope = logarithmic_slope(solver.grid.centers, rho)
        density_l1 = float(
            np.sum(np.abs(rho - rho_initial) * solver.grid.volumes)
            / initial.energy.mass
        )
        slope_change = float(np.max(np.abs(slope[resolved] - slope_initial[resolved])))
        mass_error = relative_error(snapshot.energy.mass, initial.energy.mass)
        energy_error = relative_error(snapshot.energy.total, initial.energy.total)
        estimated_core = core_radius(
            solver.grid.centers,
            rho,
            solver.grid.dr,
            cfg.minimum_resolved_cells,
            cfg.required_core_cells,
        )
        max_velocity = float(np.max(np.abs(velocity)))
        slope_samples = {
            f"r={radius:g}": float(
                np.interp(np.log(radius), np.log(solver.grid.centers), slope)
            )
            for radius in SLOPE_SAMPLE_RADII
        }
        limiter_frequency = solver.hydrodynamics.limiter_activation_frequency
        passed = (
            max_velocity <= cfg.control_velocity_tolerance
            and density_l1 <= cfg.control_density_l1_tolerance
            and slope_change <= cfg.control_slope_tolerance
            and mass_error <= cfg.mass_tolerance
            and energy_error <= cfg.energy_tolerance
            and estimated_core is None
        )
        measurements.append(
            ControlMeasurement(
                cells=cells,
                time=output_time,
                max_velocity=max_velocity,
                density_l1_change=density_l1,
                max_slope_change=slope_change,
                mass_error=mass_error,
                total_energy_error=energy_error,
                core_radius=estimated_core,
                slope_samples=slope_samples,
                positivity_limiter_activation_frequency=limiter_frequency,
                passed=passed,
            )
        )
    return measurements


def run_control_validation(
    resolutions: tuple[int, ...] = CONTROL_RESOLUTIONS,
    times: tuple[float, ...] = CONTROL_TIMES,
) -> tuple[bool, list[dict[str, object]]]:
    records: list[ControlMeasurement] = []
    for cells in resolutions:
        records.extend(validate_control_resolution(cells, times))
    return all(record.passed for record in records), [asdict(record) for record in records]
