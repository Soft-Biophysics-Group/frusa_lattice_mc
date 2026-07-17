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

# Globals
# Base rotations
ID_ROT = R.identity()
# scipy convention is positive = counter-clockwise when looking at axis through + sign
# At this point I encourage you to play with the dodecahedron to understand why we define
# rotations like this
C4ZM = R.from_euler("z", -90, degrees=True)
C2Z = R.from_euler("z", 180, degrees=True)
# 2pi/3 clockwise rotation around axis (+1, 0, +1)
C3XZM = R.from_euler("xy", [-90, -90], degrees=True)
# 2pi/3 counterclockwise rotation around axis (+1, 0, -1)
C3XMZ = R.from_euler("xy", [90, 90], degrees=True)
FACE_ROTATION = R.from_euler("xz", [180, 90], degrees=True)

BOND_ORIENTATIONS_POSITIVE = [ID_ROT, C3XZM, C3XZM**2, C3XMZ, C3XMZ**2, C4ZM]
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
    (  1,   0, - 1),
    (  1, - 1,   0),
    (  0,   1, - 1),
    ( -1,   0,   0),
    (  0, - 1,   0),
    (  0,   0,  -1),
    (- 1,   0,   1),
    (- 1,   1,   0),
    (  0, - 1,   1),
]


# LATTICE PARAMETERS
# Careful: these are column vectors!
BASIS_VECTORS = np.array([[1.0, 0.0, 1.0], [1.0, 1.0, 0.0], [0.0, 1.0, 1.0]])
