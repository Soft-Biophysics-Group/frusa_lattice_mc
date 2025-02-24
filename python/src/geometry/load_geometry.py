"""
Vincent Ouazan-Reboul, 2025

Functions to easily load geometry classes from model files.
Essentially top-level wrappers for geometry to be mostly used in other codes.
"""

import config as cfg
from geometry.cubic import CubicLattice, CubicParticle
from geometry.lattice_geometry import LatticeGeometry
from geometry.particle_geometry import ParticleGeometry
from geometry.triangular import TriangularLattice, TriangularParticle
from pathlib import Path


def geometry_from_model_file(
    model_file: str | Path = cfg.default_model_params_file,
) -> tuple[LatticeGeometry, ParticleGeometry]:
    """
    Returns a (LatticeGeometry, ParticleGeometry) tuple from a model file at `model_file`.
    If `model_file` is not specified, looks for it in `input/model_params.json`.

    Only lattices supported so far are triangular and cubic.
    """
    model_params = cfg.load_model_file(model_file)

    if model_params["lattice_name"] == "triangular":
        return TriangularLattice.from_model_file(model_file), TriangularParticle()
    elif model_params["lattice_name"] == "cubic":
        return CubicLattice.from_model_file(model_file), CubicParticle()
    else:
        print("Error: lattice not supported!")

    return LatticeGeometry(), ParticleGeometry()
