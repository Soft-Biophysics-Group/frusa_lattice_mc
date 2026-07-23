# pyright basic
"""
Vincent Ouazan-Reboul, 2026
Tools to implement geometry of square particles, in order to easily generate contact maps and plot
simulation results
"""

# mathutils is provided by bpy
import numpy as np
from scipy.spatial.transform import Rotation as R
from ..analysis.clusters import Cluster

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

# Function to measure dimensions of a rectangular cluster
def measure_effective_length_width(cluster:Cluster) -> tuple[int, int]:
    """See calculation AAH, 26-07-22 for derivation. Returns the largest dimension first,
    smallest second."""
    area = cluster.size
    perimeter = cluster.n_interfaces_with_exterior

    # Convention: length = largest dimension, width = smallest
    discriminant = np.sqrt(1 - 16 * area / perimeter**2)
    length = perimeter / 6 * (1 + discriminant)
    width  = perimeter / 6 * (1 - discriminant)

    return (length, width)
