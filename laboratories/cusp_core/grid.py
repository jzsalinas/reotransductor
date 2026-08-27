"""Uniform radial finite-volume mesh with exact spherical geometry."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SphericalGrid:
    faces: np.ndarray
    centers: np.ndarray
    volumes: np.ndarray
    areas: np.ndarray
    dr: float

    @classmethod
    def uniform(cls, cells: int, radius: float) -> "SphericalGrid":
        if cells < 4:
            raise ValueError("At least four radial cells are required")
        if radius <= 0.0:
            raise ValueError("The outer radius must be positive")
        faces = np.linspace(0.0, radius, cells + 1, dtype=np.float64)
        centers = 0.5 * (faces[:-1] + faces[1:])
        volumes = (4.0 * np.pi / 3.0) * (faces[1:] ** 3 - faces[:-1] ** 3)
        areas = 4.0 * np.pi * faces**2
        return cls(faces, centers, volumes, areas, float(radius / cells))
