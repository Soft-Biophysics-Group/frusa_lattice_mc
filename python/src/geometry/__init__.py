from .lattice_geometry import LatticeGeometry
from .particle_geometry import ParticleGeometry
from pathlib import Path
from config import load_model_file, default_model_params_file


def get_particle_lattice(model_file: str | Path = default_model_params_file):
    return (
        ParticleGeometry.from_model_file(model_file),
        LatticeGeometry.from_model_file(model_file),
    )


__all_ = ["LatticeGeometry", "ParticleGeometry", "get_particle_lattice"]
