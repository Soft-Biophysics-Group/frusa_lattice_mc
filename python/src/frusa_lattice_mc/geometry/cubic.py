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
