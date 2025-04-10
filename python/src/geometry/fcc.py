# pyright basic
"""
Vincent Ouazan-Reboul, 2025
Tools to implement geometry of rhombic-dodecahedric particles, in order to easily generate
contact maps and plot simulation results

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
# Base rotations
ID_ROT = R.identity()
# scipy convention is positive = counter-clockwise when looking at axis through + sign
# At this point I encourage you to play with the dodecahedron to understand why we define
# rotations like this
C2Z = R.from_euler("z", 180, degrees=True)
# 2pi/3 clockwise rotation around axis (+1, 0, +1)
C3XZM = R.from_euler("xy", [-90, -90], degrees=True)
# 2pi/3 counterclockwise rotation around axis (+1, 0, -1)
C3XMZ = R.from_euler("xy", [90, 90], degrees=True)
FACE_ROTATION = R.from_euler("zx", [-90, -90], degrees=True)

BOND_ORIENTATIONS_POSITIVE = [ID_ROT, C3XZM, C3XZM**2, C3XMZ, C3XMZ**2, C2Z]
BOND_ORIENTATIONS_NEGATIVE = [
    C2Z * orientation for orientation in BOND_ORIENTATIONS_POSITIVE
]
ALL_BOND_ORIENTATIONS = BOND_ORIENTATIONS_POSITIVE + BOND_ORIENTATIONS_NEGATIVE
ORIENTATION_0_VECTORS = np.array([[1, 0, 0], [0, 1, 0]])

# PARTICLE PARAMETERS
BONDS = [
    (  1,   0,   0),
    (  0,   1,   0),
    (  0,   0,   1),
    (- 1,   0,   0),
    (  1,   0, - 1),
    (  1, - 1,   0),
    ( -1,   0,   0),
    (  0,  -1,   0),
    (  0,   0,  -1),
    (  1,   0,   0),
    (- 1,   0,   1),
    (- 1,   1,   0)
]

# LATTICE PARAMETERS
BASIS_VECTORS = np.array([[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]])


class FccLattice(LatticeGeometry):
    def __init__(
        self, lx: int = 1, ly: int = 1, lz: int = 1, lattice_spacing: float = 1.0
    ):
        super().__init__(BASIS_VECTORS, BONDS, lx, ly, lz, lattice_spacing)

    @classmethod
    def from_model_file(
        cls,
        model_file: str | Path = cfg.default_model_params_file,
        lattice_spacing: float = 1.0,
    ):
        model_params = cfg.load_model_file(model_file)
        lx = model_params["lx"]
        ly = model_params["ly"]
        lz = model_params["lz"]

        return cls(lx, ly, lz, lattice_spacing)

class FccParticle(ParticleGeometry):
    def __init__(self):
        super().__init__(
            orientation_0_vectors=ORIENTATION_0_VECTORS,
            face_0_permutation_rotations=BOND_ORIENTATIONS_POSITIVE,
            rotations_around_face_0=[FACE_ROTATION,],
            opposite_face_rotation=C2Z,
            bond_rotations=ALL_BOND_ORIENTATIONS,
        )

