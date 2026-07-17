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
C6Z = R.from_euler("z", 60, degrees=True)
C2Z = R.from_euler("z", 180, degrees=True)
ORIENTATION_0_VEC = np.array([[1, 0, 0]])
BONDS = [
    ( 1,  0, 0),
    ( 0,  1, 0),
    (-1,  1, 0),
    (-1,  0, 0),
    ( 0, -1, 0),
    ( 1, -1, 0),
]
BOND_ORIENTATIONS_POSITIVE = [ID_ROT, C6Z, C6Z**2]
BOND_ROTATIONS = [
    *BOND_ORIENTATIONS_POSITIVE,
    *[C2Z * rot for rot in BOND_ORIENTATIONS_POSITIVE],
]
SR32: float = np.sqrt(3) / 2
BASIS_VECTORS = np.array([[1, 0.5, 0.], [0, SR32, 0.]])
