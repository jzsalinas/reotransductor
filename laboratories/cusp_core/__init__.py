"""Independent spherical conducting self-gravitating gas laboratory."""

from .config import LaboratoryConfig, PhysicalScales
from .solver import CuspCoreSolver, NumericalFailure

__all__ = [
    "CuspCoreSolver",
    "LaboratoryConfig",
    "NumericalFailure",
    "PhysicalScales",
]
