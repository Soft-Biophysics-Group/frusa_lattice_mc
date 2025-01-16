# pyright: basic
from typing import Iterator, List
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import matplotlib.patches as mpatches
import numpy as np
from lattice_utils import (
    lattice_coords_to_lattice_site_2d,
    lattice_site_to_lattice_coords_2d,
    get_full_sites_characteristics,
)
from pathlib import Path
import config as cfg
from matplotlib.axes import Axes
from contact_utils import ContactMapWrapper

# Global color sets
ARROW_COLORS = ["blue", "red", "green"]
# colormap of camembert contacts: forbidden contacts are in red, crystal in cyan, line in blue,
# nothing in orange
CAMEMBERT_CONTACTS_CMAP = ["red", "cyan", "blue", "orange"]
CAMEMBERT_CONTACTS = np.zeros((6,6), dtype = int)
for i in range(6):
    CAMEMBERT_CONTACTS[i, (i+3)%6] = 1
CAMEMBERT_CONTACTS[0, 2] = 2
CAMEMBERT_CONTACTS[2, 0] = 2
CAMEMBERT_CONTACTS[1, 5] = 2
CAMEMBERT_CONTACTS[5, 1] = 2

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


# Strongly inspired by Lara's code
class ParticleRepresentation:
    # 2 vertices, one in front and one in back
    # Here particle length is normalized to 1
    def __init__(self, lattice_spacing=1.0):
        """
        How the faces permute when we change the particle orientation.
        Directly lifted from triangular.h
        """
        self.lattice_spacing = lattice_spacing
        self.radius = 0.5 / sr32 * self.lattice_spacing
        self.permutations = np.array(
            [
                [0, 1, 2, 3, 4, 5],
                [1, 2, 3, 4, 5, 0],
                [2, 3, 4, 5, 0, 1],
                [3, 4, 5, 0, 1, 2],
                [4, 5, 0, 1, 2, 3],
                [5, 0, 1, 2, 3, 4],
            ]
        )
        self.opposite_bonds = np.array([3, 4, 5, 0, 1, 2])
        self.bonds = np.array([
            [1, 0],
            [1, 0],
            [-1, 1],
            [-1, 0],
            [0, -1],
            [1, -1]
        ])

        # Face colors in the reference orientation (0)
        self.colors = [
            "bf9c76ff",
            "cf938dff",
            "b996c1ff",
            "7fa5d3ff",
            "67b0b0ff",
            "92ab7dff",
        ]
        self.border_color = "black"
        self.side_length = self.radius / np.sqrt(3)
        self.n_faces = 6
        self.rotation_matrix = np.array([[1 / 2, -sr32], [sr32, 1 / 2]])
        self.face_0_corners = np.array(
            [
                # [[0.5 * self.lattice_spacing, 0.5 * self.lattice_spacing], [-sr32, sr32]]
                [0.5 * self.lattice_spacing, 0.5 * self.lattice_spacing],
                [-0.25 / sr32 * self.lattice_spacing, 0.25 / sr32 * self.lattice_spacing],
            ]
        )
        self.all_face_corners = self.init_face_coords()
        self.lattice_vectors_in_cartesian = np.array([[1, 0.5], [0, sr32]])

    def rotate_point(self, coords, n_steps=1):
        for i in range(n_steps):
            coords = np.matmul(self.rotation_matrix, coords)
        return coords

    def init_face_coords(self):
        all_faces_corners = []
        for face_nr in range(self.n_faces):
            face_corners = np.copy(self.face_0_corners)
            for i in range(2):
                this_corner = face_corners[:, i]
                face_corners[:, i] = self.rotate_point(this_corner, face_nr)
            all_faces_corners.append(face_corners)
        return all_faces_corners

    def get_bond(self, bond_vector):
        matching_coeffs = self.bonds == bond_vector
        matching_vectors = np.prod(matching_coeffs, axis=1)
        if matching_vectors.sum() != 1:
            print("The vector does not correspond to a bond!")
            return -1
        else:
            return np.argmax(matching_vectors)

    def get_neighbour_sites(self, site_index, lx, ly):
        x, y = lattice_site_to_lattice_coords_2d(site_index, lx)
        neighbours = []
        for bond in self.bonds:
            x_neighbour, y_neighbour = np.array([x, y]) + bond
            # Quick and dirty implementation of periodic boundary conditions
            if x_neighbour >= lx:
                x_neighbour -= lx
            if x_neighbour < 0:
                x_neighbour += lx
            if y_neighbour >= ly:
                y_neighbour -= ly
            if y_neighbour < 0:
                y_neighbour += ly

            neighbour_index = lattice_coords_to_lattice_site_2d(x_neighbour, y_neighbour, lx)
            neighbours.append(neighbour_index)
        return neighbours

    def get_faces_in_contact(
        self, site_1, site_2, orientation_1, orientation_2, lx, ly
    ):
        neighbours = self.get_neighbour_sites(site_1, lx, ly)
        if site_2 in neighbours:
            bond_index = np.where(np.array(neighbours) == site_2)[0][0]
        else:
            print("The sites are not neighbours!")
            return -1, -1

        if orientation_1 == -1:
            face_1 = -1
        else:
            face_1 = self.permutations[bond_index, orientation_1]
        if orientation_2 == -1:
            face_2 = -1
        else:
            face_2 = self.permutations[(bond_index+3)%6, orientation_2]

        return face_1, face_2

    def plot_particle_outline(
        self,
        x_center_lattice: int,
        y_center_lattice: int,
        lx: int,
        ax: Axes,
        squared: bool = False,
    ) -> None:
        x_center, y_center = self.lattice_to_cartesian(
            x_center_lattice, y_center_lattice
        )
        if squared:
            x_center = square_coordinates(x_center, lx)
        h = mpatches.RegularPolygon(
            (x_center, y_center),
            self.n_faces,
            radius=self.radius,
            facecolor=None,
            fill=False,
            edgecolor="black",
        )
        ax.add_artist(h)

    def plot_particle_orientation(
        self,
        x_center_lattice: int,
        y_center_lattice: int,
        lx:int,
        orientation: int,
        ax: Axes,
        squared:bool = False,
        color="blue",
    ) -> None:
        """
        Adds an arrow at the center of the particle indicating the particle orientation.
        The arrow is of length side_length/2, and points towards face 0 in the particle's
        current orientation.
        """
        a_length = self.radius / 2
        # Get the edge corresponding to face 0 in current orientation
        # And deduce the arrow orientation in lattice coordinates
        edge = np.where(self.permutations[orientation] == 0)[0][0]
        arrow_vector = np.array([a_length, 0])
        arrow_vector = self.rotate_point(arrow_vector, edge)
        # FancyArrow uses arrow base location as input, so we'll need that
        x_center_cartesian, y_center_cartesian = self.lattice_to_cartesian(
            x_center_lattice, y_center_lattice
        )
        if squared:
            x_center_cartesian = square_coordinates(x_center_cartesian, lx)
        x_arrow_base = x_center_cartesian - arrow_vector[0] / 2
        y_arrow_base = y_center_cartesian - arrow_vector[1] / 2

        arrow = mpatches.FancyArrow(
            x_arrow_base,
            y_arrow_base,
            arrow_vector[0],
            arrow_vector[1],
            width=a_length / 5,
            color=color,
        )
        ax.add_artist(arrow)

    def lattice_to_cartesian(self, x_lattice: int, y_lattice: int) -> int:
        return (
            x_lattice * self.lattice_vectors_in_cartesian[:, 0]
            + y_lattice * self.lattice_vectors_in_cartesian[:, 1]
        )

    def plot_contact_camembert(
        self,
        x_center,
        y_center,
        face_1,
        face_2,
        bond,
        lx,
        ax,
        squared:bool = False,
        cmap=CAMEMBERT_CONTACTS_CMAP,
    ):
        centered_face_corners = self.all_face_corners[bond]
        x_center_cartesian, y_center_cartesian = self.lattice_to_cartesian(x_center, y_center)
        if squared:
            x_center_cartesian = square_coordinates(x_center_cartesian, lx)
        particles_face_corners = centered_face_corners + np.array(
            [
                [x_center_cartesian, x_center_cartesian],
                [y_center_cartesian, y_center_cartesian],
            ]
        )
        if face_1 == -1 or face_2 == -1:
            color = cmap[-1]
        else:
            color = cmap[CAMEMBERT_CONTACTS[face_1, face_2]]

        ax.plot(particles_face_corners[0, :], particles_face_corners[1, :], color=color)
        return


    def plot_all_contacts_camembert(
        self,
        ax: Axes,
        results_index: int = -1,
        results_folder: str | Path = "",
        results_file: str | Path = "",
        squared = False
    ):
        results = cfg.load_structure(results_index, results_folder, results_file)
        lx = cfg.load_model_file()["lx"]
        ly = cfg.load_model_file()["ly"]

        for (
            site,
            _,
            orientation
        ) in get_full_sites_characteristics(results):
            x_1, y_1 = lattice_site_to_lattice_coords_2d(site, lx)
            neighbours = self.get_neighbour_sites(site, lx, ly)

            for bond, neighbour in enumerate(neighbours):
                # Let's avoid plotting the same contact twice
                face_1, face_2 = self.get_faces_in_contact(
                    site, neighbour, orientation, results[1, neighbour], lx, ly
                )
                self.plot_contact_camembert(
                    x_1, y_1, face_1, face_2, bond, lx, ax, squared
                )
        return

    def plot_result_outlines(
        self,
        results_index: int = -1,
        results_file: str | Path = "",
        results_folder: str | Path = "",
        ax: Axes | None = None,
        squared: bool = False,
        **kwargs,
    ) -> tuple[Figure, Axes]:
        if ax is None:
            fig, ax = plt.subplots(**kwargs)
        else:
            fig = ax.get_figure()
        results = cfg.load_structure(results_index, results_folder, results_file)
        lx = cfg.load_model_file()["lx"]
        for (
            site,
            _,
            _,
        ) in get_full_sites_characteristics(results):
            x_lattice, y_lattice = lattice_site_to_lattice_coords_2d(site, lx)
            self.plot_particle_outline(x_lattice, y_lattice, lx, ax, squared = squared)

        return fig, ax

    def plot_results_arrows(
        self,
        results_index: int = -1,
        results_folder: str | Path = "",
        results_file: str | Path = "",
        squared: bool = False,
        **kwargs,
    ) -> tuple[Figure, Axes]:
        fig, ax = plt.subplots()
        results = cfg.load_structure(results_index, results_folder,  results_file)
        lx = cfg.load_model_file()["lx"]
        ly = cfg.load_model_file()["ly"]

        for (
            site,
            ptype,
            orientation,
        ) in get_full_sites_characteristics(results):
            x_lattice, y_lattice = lattice_site_to_lattice_coords_2d(site, lx)
            color = ARROW_COLORS[ptype]
            self.plot_particle_outline(x_lattice, y_lattice, lx, ax, squared)
            self.plot_particle_orientation(
                x_lattice, y_lattice, lx, orientation, ax, color=color, squared=squared
            )

        # Adjusting the viewing window
        x_min = -self.radius * 1.1
        y_min = -self.radius * 1.1
        if squared:
            x_max = lx + self.radius
            y_max = ly * sr32 + self.radius
        else:
            x_max, y_max = self.lattice_to_cartesian(lx, ly) + self.side_length
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

        return fig, ax

def square_coordinates(x_cartesian, lx):
    if x_cartesian >= lx:
        return x_cartesian - lx
    else:
        return x_cartesian
