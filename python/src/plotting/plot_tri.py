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

ARROW_COLORS = ["blue", "red", "green"]

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
        self.side_length = lattice_spacing * np.sqrt(3) / 3
        self.n_faces = 6
        self.rotation_matrix = np.array([[1 / 2, sr32], [-sr32, 1 / 2]])
        self.face_0_corners = [
            np.array([0, 0]),
            np.array([sr32 * self.side_length, -self.side_length / 2]),
            np.array([sr32 * self.side_length, self.side_length / 2]),
        ]
        self.all_face_corners = self.init_face_coords()
        self.lattice_vectors_in_cartesian = np.array([[1, 0], [1/2, sr32]])

    def rotate_point(self, coords, n_steps):
        for i in range(n_steps):
            coords[...] = np.matmul(self.rotation_matrix, coords)

    def init_face_coords(self):
        all_faces_corners = []
        for face_nr in range(self.n_faces):
            face_corners = np.copy(self.face_0_corners)
            for corner in face_corners:
                self.rotate_point(corner, face_nr)
            all_faces_corners.append(face_corners)
        return all_faces_corners

    def plot_particle_outline(
        self, x_center_lattice: int, y_center_lattice: int, ax: Axes
    ) -> None:
        x_center, y_center = self.lattice_to_cartesian(
            x_center_lattice, y_center_lattice
        )
        h = mpatches.RegularPolygon(
            (x_center, y_center),
            self.n_faces,
            radius=self.side_length,
            facecolor=None,
            fill=False,
            edgecolor="black",
        )
        ax.add_artist(h)

    def plot_particle_orientation(
        self,
        x_center_lattice: int,
        y_center_lattice: int,
        orientation: int,
        ax: Axes,
        color="blue",
    ) -> None:
        """
        Adds an arrow at the center of the particle indicating the particle orientation.
        The arrow is of length side_length/2, and points towards face 0 in the particle's
        current orientation.
        """
        a_length = self.side_length / 2
        # Get the edge corresponding to face 0 in current orientation
        # And deduce the arrow orientation in lattice coordinates
        edge = np.where(self.permutations[orientation] == 0)[0][0]
        arrow_vector = np.array([a_length, 0])
        self.rotate_point(arrow_vector, edge)
        # FancyArrow uses arrow base location as input, so we'll need that
        x_center_cartesian, y_center_cartesian = self.lattice_to_cartesian(
            x_center_lattice, y_center_lattice
        )
        x_arrow_base = x_center_cartesian - arrow_vector[0] / 2
        y_arrow_base = y_center_cartesian - arrow_vector[1] / 2

        arrow = mpatches.FancyArrow(
            x_arrow_base,
            y_arrow_base,
            arrow_vector[0],
            arrow_vector[1],
            width=a_length / 10,
            color=color,
        )
        ax.add_artist(arrow)

    def lattice_to_cartesian(self, x_lattice: int, y_lattice: int) -> int:
        return (
            x_lattice * self.lattice_vectors_in_cartesian[0, :]
            + y_lattice * self.lattice_vectors_in_cartesian[1, :]
        )

    def plot_result_outlines(
        self,
        results_index: int = -1,
        results_file: str | Path = "",
        ax: Axes | None = None,
        **kwargs
    ) -> tuple[Figure, Axes]:
        if ax is None:
            fig, ax = plt.subplots(**kwargs)
        else:
            fig = ax.get_figure()
        results = cfg.load_structure(results_index, results_file)
        lx = cfg.load_model_file()["lx"]
        for (
            site,
            _,
            _,
        ) in get_full_sites_characteristics(results):
            x_lattice, y_lattice = lattice_site_to_lattice_coords_2d(site, lx)
            self.plot_particle_outline(x_lattice, y_lattice, ax)

        return

    def plot_results_arrows(
        self, results_index: int = -1, results_file: str | Path = "", **kwargs
    ) -> tuple[Figure, Axes]:
        fig, ax = plt.subplots()
        results = cfg.load_structure(results_index, results_file)
        lx = cfg.load_model_file()["lx"]
        ly = cfg.load_model_file()["ly"]

        for (
            site,
            ptype,
            orientation,
        ) in get_full_sites_characteristics(results):
            x_lattice, y_lattice = lattice_site_to_lattice_coords_2d(site, lx)
            color = ARROW_COLORS[ptype]
            self.plot_particle_outline(x_lattice, y_lattice, ax)
            self.plot_particle_orientation(x_lattice, y_lattice, orientation, ax, color = color)

        # Adjusting the viewing window
        x_min = -self.side_length
        y_min = -self.side_length
        x_max, y_max = self.lattice_to_cartesian(lx, ly) + self.side_length
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

        return fig, ax
