# pyright basic
"""
Vincent Ouazan-Reboul, 2025
Tools to implement geometry of cubic particles, in order to easily generate contact maps and plot
simulation results

TODOS:
- Find a better way to explain the docstrings
"""

# mathutils is provided by bpy
import numpy as np
from scipy.spatial.transform import Rotation as R
from geometry.particle_geometry import ParticleGeometry
from geometry.lattice_geometry import LatticeGeometry
import config as cfg
from pathlib import Path

# Globals
ID_ROT = R.identity()

# PARTICLE PARAMETERS
# scipy convention is positive = counter-clockwise when looking at axis through + sign
# At this point I encourage you to play with the cube to understand why we define rotations like
# this
C4XM = R.from_euler("x", -90, degrees=True)
C4Y = R.from_euler("y", 90, degrees=True)
C4ZM = R.from_euler("z", -90, degrees=True)
C2Z = C4ZM * C4ZM
C4XM_POWERS = [C4XM**i for i in range(4)]
BOND_ORIENTATIONS_POSITIVE = [ID_ROT, C4ZM, C4Y]  # , C2Z, C2Z * C4Z, C2Z * C4Y]
ALL_BOND_ORIENTATIONS = [*BOND_ORIENTATIONS_POSITIVE, C2Z, C2Z * C4ZM, C2Z * C4Y]
ORIENTATION_0_VECTORS = np.array([[1, 0, 0], [0, 1, 0]])
BONDS = [
    ( 0,  0,  1),
    ( 0,  1,  0),
    ( 1,  0,  0),
    (-1,  0,  0),
    ( 0, -1,  0),
    ( 0,  0, -1),
]
# LATTICE PARAMETERS
BASIS_VECTORS = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])


class CubicLattice(LatticeGeometry):
    def __init__(
        self, lx: int = 1, ly: int = 1, lz: int = 1, lattice_spacing: float = 1.0
    ):
        super().__init__(BASIS_VECTORS, BONDS, lx, ly, lz, lattice_spacing)

    @classmethod
    def from_model_file(cls, model_file: str | Path = cfg.default_model_params_file,
                        lattice_spacing: float = 1.0):
        model_params = cfg.load_model_file(model_file)
        lx = model_params["lx"]
        ly = model_params["ly"]
        lz = model_params["lz"]

        return cls(lx, ly, lz, lattice_spacing)


class CubicParticle(ParticleGeometry):
    def __init__(self):
        super().__init__(
            orientation_0_vectors=ORIENTATION_0_VECTORS,
            face_0_permutation_rotations=BOND_ORIENTATIONS_POSITIVE,
            rotations_around_face_0=C4XM_POWERS,
            opposite_face_rotation=C2Z,
            bond_rotations=ALL_BOND_ORIENTATIONS,
        )


class CubicGeometry:
    particle = CubicParticle()

    def __init__(
        self, lx: int = 1, ly: int = 1, lz: int = 1, lattice_spacing: float = 1.0
    ):
        # We keep track of only x and y internal vectors, z is redundant
        self.lattice = CubicLattice(lx, ly, lz, lattice_spacing)
