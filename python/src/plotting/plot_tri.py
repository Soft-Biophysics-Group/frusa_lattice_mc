# pyright: basic
from typing import Iterator, List
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import matplotlib.patches as mpatches
import numpy as np
from geometry.triangular import TriangularLattice, TriangularParticle
from pathlib import Path
import config as cfg
from matplotlib.axes import Axes
from contact_utils import ContactMapWrapper

# Global color sets
ARROW_COLORS = ["black", "blue", "red", "green"]
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

Some important general notes:
- A (x, y) pair in lattice coordinates is hashed into a unique site index i using the formula
    i = x + lx * y, where lx is the number of sites in the x direction.
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
    """
    Class serving as an intermediary to plot hexagonal particles.
    """
    # 2 vertices, one in front and one in back
    # Here particle length is normalized to 1
    def __init__(self, lx:int = 1, ly:int = 1, lattice_spacing=1.0):
        """
        How the faces permute when we change the particle orientation.
        Directly lifted from triangular.h
        """
        self.lattice_spacing = lattice_spacing
        self.radius = 0.5 / sr32 * self.lattice_spacing
        self.lattice = TriangularLattice(lx, ly, lattice_spacing)
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
                [0.5 * self.lattice_spacing         , 0.5 * self.lattice_spacing        ],
                [-0.25 / sr32 * self.lattice_spacing, 0.25 / sr32 * self.lattice_spacing],
                [0                                  , 0                                 ]
            ]
        )
        self.all_face_corners = self.init_face_coords()

    @classmethod
    def from_model_file(
        cls, model_file=cfg.default_model_params_file, lattice_spacing=1.0
    ):
        model_dict = cfg.load_model_file(model_file)
        lx = model_dict["lx"]
        ly = model_dict["ly"]

        return cls(lx, ly, lattice_spacing)

    def init_face_coords(self):
        """
        Creates the coordinates of all faces' corners
        """
        all_faces_corners = []
        for face_nr in range(self.n_faces):
            face_corners = np.copy(self.face_0_corners)
            rotation = TriangularParticle().orientation_rotations[face_nr]
            for i in range(2):
                this_corner = face_corners[:, i]
                face_corners[:, i] = rotation.apply(this_corner)
            all_faces_corners.append(face_corners)
        return all_faces_corners

    # ----- SINGLE-PARTICLE PLOTTING ----- 
    def plot_particle_outline(
        self,
        x_center_lattice: int,
        y_center_lattice: int,
        ax: Axes,
        squared: bool = False,
        fill_color: str = ""
    ) -> None:
        """
        Plots the outline of a particle, optionally filled with color `fill_color`
        """
        x_center, y_center = self.lattice.lattice_to_cartesian(
            x_center_lattice, y_center_lattice
        )
        if fill_color == "":
            fill = False
        if squared:
            x_center = self.square_coordinates(x_center)
        h = mpatches.RegularPolygon(
            (x_center, y_center),
            self.n_faces,
            radius=self.radius,
            facecolor=None,
            fill=fill,
            edgecolor="black",
        )
        ax.add_artist(h)

    def plot_particle_orientation(
        self,
        x_center_lattice: int,
        y_center_lattice: int,
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
        arrow_vector = np.array([a_length, 0, 0])
        orientation_rotation = TriangularParticle().orientation_rotations[orientation]
        arrow_vector = orientation_rotation.inv().apply(arrow_vector)
        # FancyArrow uses arrow base location as input, so we'll need that
        x_center_cartesian, y_center_cartesian = self.lattice.lattice_to_cartesian(
            x_center_lattice, y_center_lattice
        )
        if squared:
            x_center_cartesian = self.square_coordinates(x_center_cartesian)
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

    def plot_contact(
        self,
        x_center,
        y_center,
        bond,
        ax,
        color,
        squared:bool = False,
    ):
        """
        Plots the contact between particles into different colors depending on their relative
        orientations.
        This is mostly used when trying to make "camembert" aggregates: we want to easily
        visualize crystalline domains and the defect lines between them.

        contact_to_colors maps every possible pair of faces to a contact color

        Note that this sometimes displays incorrect colors, but I have not had the time to debug
        it.
        """
        centered_face_corners = self.all_face_corners[bond]
        x_center_cartesian, y_center_cartesian = self.lattice.lattice_to_cartesian(
            x_center, y_center
        )
        if squared:
            x_center_cartesian = self.square_coordinates(x_center_cartesian)
        particles_face_corners = centered_face_corners + np.array(
            [
                [x_center_cartesian, x_center_cartesian],
                [y_center_cartesian, y_center_cartesian],
                [0                 , 0                 ]
            ]
        )

        # print(particles_face_corners)
        ax.plot(
            particles_face_corners[0, :], particles_face_corners[1, :], color=color
        )
        return

    # ----- PLOTTING SIMULATION RESULTS -----
    def plot_result_outlines(
        self,
        results_index: int = -1,
        results_file: str | Path = "",
        results_folder: str | Path = "",
        ax: Axes | None = None,
        squared: bool = False,
        **kwargs,
    ) -> tuple[Figure, Axes]:
        """
        Plots the outline of the particles contained in a structure file.

        ## Where this looks for structures

        If results_file points to a structure.dat file, plots the corresponding results.
        If `results_folder` and `results_index` are specified, plots the
        structure_results_index.dat file contained in results_folder.
        If only `results_folder` is specified or `results_index==-1`, plots the final structure
        contained in `results_folder`.
        If `results_folder` is not specified, looks for the results in `<root
        folder>/data/structures` folder.

        ## Extra parameters:
        - `squared` is a boolean which, if set to True, will use the periodic boundary
          conditions to wrap the lattice into a square window rather.
        - Any additional keyword arguments will be passed to matplotlib to create the figure.
        """
        if ax is None:
            fig, ax = plt.subplots(**kwargs)
        else:
            fig = ax.get_figure()
        results = cfg.load_structure(results_index, results_folder, results_file)
        for site in cfg.get_full_sites(results):
            x_lattice, y_lattice, _ = self.lattice.lattice_site_to_lattice_coords(site)
            self.plot_particle_outline(x_lattice, y_lattice, ax, squared = squared)

        return fig, ax

    def plot_results_arrows(
        self,
        results_index: int = -1,
        results_folder: str | Path = "",
        results_file: str | Path = "",
        squared: bool = False,
        **kwargs,
    ) -> tuple[Figure, Axes]:
        """
        Plots the particles contained in a structure file, with an arrow indicating the
        orientation of the particle (pointing towards face 0).
        The color of the arrow depends on the particle species, in an order ginve by the
        ARROW_COLORS global.

        ## Where this looks for structures

        If results_file points to a structure.dat file, plots the corresponding results.
        If `results_folder` and `results_index` are specified, plots the
        structure_results_index.dat file contained in results_folder.
        If only `results_folder` is specified or `results_index==-1`, plots the final structure
        contained in `results_folder`.
        If `results_folder` is not specified, looks for the results in `<root
        folder>/data/structures` folder.

        ## Extra parameters:
        - `squared` is a boolean which, if set to True, will use the periodic boundary
          conditions to wrap the lattice into a square window rather.
        - Any additional keyword arguments will be passed to matplotlib to create the figure.
        """
        fig, ax = plt.subplots()
        results = cfg.load_structure(results_index, results_folder,  results_file)

        for (
            site,
            ptype,
            orientation,
        ) in cfg.get_full_sites_characteristics(results):
            x_lattice, y_lattice, _ = self.lattice.lattice_site_to_lattice_coords(site)
            color = ARROW_COLORS[ptype]
            self.plot_particle_outline(x_lattice, y_lattice, ax, squared)
            self.plot_particle_orientation(
                x_lattice, y_lattice, orientation, ax, color=color, squared=squared
            )

        # Adjusting the viewing window
        x_min = -self.radius * 1.1
        y_min = -self.radius * 1.1
        if squared:
            x_max = self.lattice.lx + self.radius
            y_max = self.lattice.ly * sr32 + self.radius
        else:
            x_max, y_max = (
                self.lattice.lattice_to_cartesian(self.lattice.lx, self.lattice.ly)
                + self.side_length
            )
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

        return fig, ax

    def plot_contacts(
        self,
        ax: Axes,
        contact_to_color,
        results_index: int = -1,
        results_folder: str | Path = "",
        results_file: str | Path = "",
        squared=False,
    ):
        """
        Plots the contacts between particles to visualise crystalline domains and defects lines.
        Needs to be put on top of an existing `ax` Axes objects, because this will never be
        useful on its own.

        `contact_to_color` is a mappable mapping a pair of faces to a color numpy accepts
        (typically a color name or a hex code)

        ## Where this looks for structures

        If results_file points to a structure.dat file, plots the corresponding results.
        If `results_folder` and `results_index` are specified, plots the
        structure_results_index.dat file contained in results_folder.
        If only `results_folder` is specified or `results_index==-1`, plots the final structure
        contained in `results_folder`.
        If `results_folder` is not specified, looks for the results in `<root
        folder>/data/structures` folder.

        ## Extra parameters:
        - `squared` is a boolean which, if set to True, will use the periodic boundary
          conditions to wrap the lattice into a square window rather.
        - Any additional keyword arguments will be passed to matplotlib to create the figure.
        """
        results = cfg.load_structure(results_index, results_folder, results_file)

        for site, _, orientation in cfg.get_full_sites_characteristics(results):
            x_1, y_1, _= self.lattice.lattice_site_to_lattice_coords(site)
            # print(f"site {site}, coordinates {x_1, y_1}")
            neighbours = self.lattice.get_neighbour_sites(site)

            for bond, neighbour in enumerate(neighbours):
                x_2, y_2, _= self.lattice.lattice_site_to_lattice_coords(neighbour)
                # print( f"\tbond {bond}, neighbour {neighbour} ({x_2, y_2})")
                neighbour_orientation = results[1, neighbour]
                    # f"\tbond {bond}, neighbour {neighbour} ({x_2, y_2}), with orientation {neighbour_orientation}"
                # Let's avoid plotting the same contact twice
                face_1, face_2 = TriangularParticle().get_faces_in_contact(
                    orientation, neighbour_orientation, bond
                )
                # print(f"\t\tfaces:{face_1}, {face_2}")
                color = contact_to_color[face_1, face_2]
                # print(f"\t\t\t{color}")
                self.plot_contact(x_1, y_1, bond, ax, color, squared)
        return

    def square_coordinates(self, x_cartesian):
        if x_cartesian >= self.lattice.lx:
            return x_cartesian - self.lattice.lx
        else:
            return x_cartesian
