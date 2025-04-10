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


def plot_chain(orientations: List[int], ax: Axes):
    n_sites = len(orientations)
    for site in range(n_sites):
        # Only plot non-empty sites
        orientation = orientations[site]
        if orientation != -1:
            ParticleRepresentation().plot(site, orientation, ax)
    ax.set_xlim(-1, n_sites + 1)
    ax.set_ylim(-0.3, 0.3)


# Strongly inspired by Lara's code
class ParticleRepresentation:
    # 2 vertices, one in front and one in back
    # Here particle length is normalized to 1
    def __init__(self, width:float=1.0, height:float=0.5):
        """
        How the faces permute when we change the particle orientation.
        Here it's trivial: we have 2 faces, 0 and 1, which are inverted when
        the orientation goes from 0 to 1!
        """
        self.permutations = np.array([[0, 1],[1,0]])
        # Face colors in the reference orientation (0)
        self.colors = ["black", "red"]
        self.width = width
        self.height = height
        self.face_coords = np.array([[-width/2, -height/2], [0, -height/2]])
        self.n_faces = 2

    def plot(self, x_center: int, orientation: int,
             ax: Axes) -> None:
        permuted_faces = self.permutations[orientation]
        for i in range(self.n_faces):
            bottom_left_coords = (x_center+self.face_coords[i, 0],
                                  self.face_coords[i, 1])
            ax.add_patch(Rectangle(bottom_left_coords, self.width/2,
                                   self.height,
                                   facecolor=self.colors[permuted_faces[i]]))
