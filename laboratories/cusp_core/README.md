# Independent Cusp-Core Laboratory

This laboratory tests a generic conducting self-gravitating ideal gas. It does
not currently represent SIDM, baryonic plasma, or a physical dark-matter
model.

It is deliberately isolated from the Reotransductor cosmological engine,
observational analysis, and experiments. It contains no `tau` field, `I`
field, eons, CCC, holographic memory, cosmological expansion, synthetic mass
source, observational calibration, or SPARC input.

The sole question reserved for the later protected scientific experiment is:

> Can conductive energy transport causally transform a resolved cuspy
> self-gravitating halo into a resolved cored halo?

Passing the software tests or the nonconducting control establishes neither
that a physical dark-matter candidate has been specified nor that conduction
will form a core. The fixed conductivity is a controlled transport parameter,
not a microphysical derivation.

## Package map

- `config.py`: frozen model values, numerical tolerances, and physical scales.
- `grid.py`: exact spherical face areas and cell volumes.
- `equilibrium.py`: NFW mass, exact density cell averages, and hydrostatic
  pressure cell averages without evaluating the central singularity.
- `hydrodynamics.py`: MUSCL-MC/HLLC Euler-Poisson operator with equilibrium
  subtraction and SSP-RK2 time integration.
- `conduction.py`: conservative Crank-Nicolson Fourier transport through face
  luminosities.
- `diagnostics.py`: mass, energy, potential, slope, core radius, and rotation
  curve.
- `solver.py`: symmetric Strang-split evolution.
- `validation.py` and `run_control.py`: protected `K_hat=0` validation only.
- `tests/`: isolated verification suite.

## Verification and control

Run the complete local suite with the repository virtual environment:

```bash
.venv/bin/python -m unittest discover \
  -s laboratories/cusp_core/tests -t . -v
```

Run only the pre-registered nonconducting control:

```bash
.venv/bin/python -m laboratories.cusp_core.run_control
```

The control command does not expose or execute the `K_hat=1.5` scientific
branch. The protected conductive experiment must only be added or invoked
after explicit project-owner approval.

The reference realization `r_s=1 kpc`, `rho_s=1 Msun pc^-3` fixes the unit
mapping. Rotation speeds emitted in km/s are conditional on that mapping and
are not first-principles predictions of absolute galaxy velocities.
