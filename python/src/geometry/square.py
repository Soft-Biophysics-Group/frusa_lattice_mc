# pyright basic
"""
Vincent Ouazan-Reboul, 2025
Tools to implement geometry of hexagonal particles, in order to easily generate contact maps and plot
simulation results

TODOS:
- Find a better way to explain the docstrings
"""

# mathutils is provided by bpy
import numpy as np
from scipy.spatial.transform import Rotation as R

# Globals
ID_ROT = R.identity()
C4Z = R.from_euler("z", 90, degrees=True)
C2Z = R.from_euler("z", 180, degrees=True)
ORIENTATION_0_VEC = np.array([[1, 0, 0]])
BONDS = [
    ( 1,  0, 0),
    ( 0,  1, 0),
    (-1,  0, 0),
    ( 0, -1, 0),
]
BOND_ORIENTATIONS_POSITIVE = [ID_ROT, C4Z]
BOND_ROTATIONS = [
    *BOND_ORIENTATIONS_POSITIVE,
    *[C2Z * rot for rot in BOND_ORIENTATIONS_POSITIVE],
]
BASIS_VECTORS = np.array([[1., 0., 0.], [0., 1., 0.]])
