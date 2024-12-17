from typing import Iterator, List
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.axes import Axes
import numpy as np


"""
Code is strongly inspired from Lara's; see 
/Users/vincent/research/projects/23_frustratedSelfAssembly/2311_laraSimulationCode/mySelfAssembly2/src/latticeparticles/LatticeTools.py
"""

"""
For now: we get by with only one particle type. If we want to make
fancier figures, or to have several particle types, we will need to add several
particle representations. This is left as an exercise to the user.
"""

"""
TODOS:
    - Wrap the objects reused for several lattices in a plotting_utils.py file.
    - Add support for several particle types: make the ParticleRepresentation
      class accept several colors for instance
"""

# To speed up calculations
sr32 = np.sqrt(3) / 2

def plot_chain(orientations: List[int],
               ax: Axes):
    n_sites = len(orientations)
    for site in range(n_sites):
        # Only plot non-empty sites
        orientation = orientations[site]
        if orientation != -1:
            ParticleRepresentation().plot(site, orientation, ax)
    ax.set_xlim(-1, n_sites+1)
    ax.set_ylim(-0.3, 0.3)


# Strongly inspired by Lara's code
class ParticleRepresentation:
    # 2 vertices, one in front and one in back
    # Here particle length is normalized to 1
    def __init__(self, side_length = 1.0):
        """
        How the faces permute when we change the particle orientation.
        Directly lifted from triangular.h
        """
        self.permutations = np.array([
            [0, 1, 2, 3, 4, 5],
            [1, 2, 3, 4, 5, 0],
            [2, 3, 4, 5, 0, 1],
            [3, 4, 5, 0, 1, 2],
            [4, 5, 0, 1, 2, 3],
            [5, 0, 1, 2, 3, 4],
        ])
        # Face colors in the reference orientation (0)
        self.colors = ["bf9c76ff", "cf938dff", "b996c1ff", "7fa5d3ff",
                       "67b0b0ff", "92ab7dff"]
        self.border_color = "black"
        self.side_length = side_length
        self.n_faces = 6
        self.rotation_matrix = np.array([
            [1/2, -sr32],
            [sr32, 1/2]
        ])
        self.face_0_corners = [np.array([0, 0]),
                               np.array([sr32 * self.side_length,
                                         -self.side_length/2]),
                               np.array([sr32 * self.side_length,
                                         self.side_length/2])]
        self.all_face_corners = self.init_face_coords()

    def rotate_point(self, coords, n_steps):
        for i in range(n_steps):
            coords = np.matmul(self.rotation_matrix, coords)

    def init_face_coords(self):
        all_faces_corners = []
        for face_nr in range(self.n_faces):
            face_corners = np.copy(self.face_0_corners)
            for corner in face_corners:
                self.rotate_point(corner, face_nr)
            all_faces_corners.append(face_corners)
        return all_faces_corners


    # V: I stopped before working on this function
    def plot(self, x_center: int, y_center:int, orientation: int,
             ax: Axes) -> None:
        permuted_faces = self.permutations[orientation]
        for i in range(self.n_faces):
            coords = self.all_face_corners[i]
            # Outline: use a patch_collection once I figure out how to create
            # them
