# pyright: basic
"""
Lattice-independent machinery for 2D plotting of lattice particles.

`ParticleRepresentation2D` holds all the plotting logic that does not depend on the specific
lattice: particles are drawn as space-filling regular polygons (the Voronoi cell of the site),
with faces/orientations resolved through the generic `LatticeGeometry` / `ParticleGeometry`
abstractions. Concrete lattices subclass it and only specify `lattice_name`, `n_faces` and
`colors`.

Code is strongly inspired from Lara's; see
/Users/vincent/research/projects/23_frustratedSelfAssembly/2311_laraSimulationCode/mySelfAssembly2/src/latticeparticles/LatticeTools.py
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import matplotlib.patches as mpatches
import numpy as np
from geometry import LatticeGeometry
from pathlib import Path
import config as cfg

# Arrow color per particle species
ARROW_COLORS = ["black", "blue", "red", "green"]
# Contact colors
DEFAULT_CONTACT_COLORS = {
    "mismatch": "#8b0000ff",
    "crystal": "#66cdaaff",
    "defect": "#0000cdff",
    "empty": "#ff856fff",
}


def load_structure_from_args(
    results_index: int | None = None,
    results_folder: str | Path | None = None,
    results_file: str | Path | None = None,
) -> np.ndarray:
    """
    Resolves structure location from the various argument combinations
    and returns the loaded structure array.
    """
    return cfg.load_structure(results_index, results_folder, results_file)


class ParticleRepresentation2D:
    """
    Base class plotting 2D lattice particles as space-filling regular polygons.

    Subclasses set the class attributes:
        lattice_name: name registered in `LatticeGeometry` / `ParticleGeometry` (e.g. "square").
        n_faces: number of faces == number of nearest-neighbour bonds == polygon sides.
        colors: face colors in the reference orientation (0), one per face.
    """

    lattice_name: str
    n_faces: int
    colors: list[str]

    # Concrete subclasses register here by lattice_name (see __init_subclass__)
    _representations: dict[str, type["ParticleRepresentation2D"]] = {}

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if getattr(cls, "lattice_name", None):
            ParticleRepresentation2D._representations[cls.lattice_name] = cls

    def __init__(self, lx: int = 1, ly: int = 1, lattice_spacing: float = 1.0):
        self.lattice_spacing = lattice_spacing
        # The polygon is the Voronoi cell of the site: apothem (centre->edge) is half the
        # nearest-neighbour distance, and circumradius follows from the number of sides.
        apothem = 0.5 * lattice_spacing
        self.radius = apothem / np.cos(np.pi / self.n_faces)
        self.side_length = self.radius / np.sqrt(3)
        self.lattice = LatticeGeometry.from_lattice_name(
            self.lattice_name, lx, ly, lattice_spacing=lattice_spacing
        )
        # Reuse the particle geometry attached to the lattice rather than rebuilding it
        self.particle = self.lattice.particle_geometry
        self.border_color = "black"
        # Orient the RegularPolygon so one edge (face 0) faces +x, i.e. bond 0
        self.polygon_orientation = -np.pi / 2 - np.pi / self.n_faces
        # Face 0 is the edge facing +x: it sits at x = apothem, spanning +-half_edge in y
        half_edge = self.radius * np.sin(np.pi / self.n_faces)
        self.face_0_corners = np.array(
            [
                [apothem, apothem],
                [-half_edge, half_edge],
                [0, 0],
            ]
        )
        self.all_face_corners = self.init_face_coords()

    @classmethod
    def from_lattice_name(
        cls, lattice_name: str, lx: int = 1, ly: int = 1, lattice_spacing: float = 1.0
    ) -> "ParticleRepresentation2D":
        """Returns the representation for `lattice_name` (e.g. "triangular", "square")."""
        # Import the concrete modules so their subclasses register themselves
        from . import triangular, square  # noqa: F401

        if lattice_name not in cls._representations:
            raise ValueError(
                f"No 2D plotting for lattice '{lattice_name}'. "
                f"Available: {sorted(cls._representations)}"
            )
        return cls._representations[lattice_name](lx, ly, lattice_spacing)

    @classmethod
    def from_model_file(
        cls, model_file=cfg.default_model_params_file, lattice_spacing=1.0
    ) -> "ParticleRepresentation2D":
        """Builds the representation matching the lattice recorded in the model file."""
        model_dict = cfg.load_model_file(model_file)
        # These keys are optional in the ModelParams TypedDict but required here
        lattice_name = model_dict.get("lattice_name")
        lx, ly = model_dict.get("lx"), model_dict.get("ly")
        if lattice_name is None or lx is None or ly is None:
            raise KeyError("model file must define 'lattice_name', 'lx' and 'ly'")
        return cls.from_lattice_name(lattice_name, lx, ly, lattice_spacing)

    def _resolve_fig_ax(self, ax: Axes | None, **kwargs) -> tuple[Figure, Axes]:
        """Returns the (figure, axes) to draw on, creating a new figure when ax is None."""
        if ax is not None:
            fig = ax.get_figure()
            if not isinstance(fig, Figure):
                raise ValueError("The provided Axes is not attached to a Figure.")
            return fig, ax
        # Fresh names: reusing `ax` would keep its declared Axes | None type
        new_fig, new_ax = plt.subplots(**kwargs)
        return new_fig, new_ax

    def init_face_coords(self):
        """Creates the coordinates of all faces' corners by rotating face 0."""
        all_faces_corners = []
        for face_nr in range(self.n_faces):
            face_corners = np.copy(self.face_0_corners)
            rotation = self.particle.orientation_rotations[face_nr]
            for i in range(2):
                face_corners[:, i] = rotation.apply(face_corners[:, i])
            all_faces_corners.append(face_corners)
        return all_faces_corners

    # ----- SINGLE-PARTICLE PLOTTING -----
    def plot_particle_outline(
        self,
        x_center_lattice: int,
        y_center_lattice: int,
        ax: Axes,
        squared: bool = False,
        fill_color: str = "",
    ) -> None:
        """Plots the outline of a particle, optionally filled with color `fill_color`."""
        x_center, y_center = self.lattice.lattice_to_cartesian(
            x_center_lattice, y_center_lattice
        )
        if fill_color == "":
            fill = False
            facecolor = None
        else:
            fill = True
            facecolor = fill_color
        if squared:
            x_center = self.square_coordinates(x_center)
        h = mpatches.RegularPolygon(
            (x_center, y_center),
            self.n_faces,
            radius=self.radius,
            orientation=self.polygon_orientation,
            facecolor=facecolor,
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
        squared: bool = False,
        color="blue",
    ) -> None:
        """
        Adds an arrow at the center of the particle pointing towards face 0 in the particle's
        current orientation.
        """
        a_length = self.radius / 2
        # In orientation o, face 0 points along bond -o, hence the inverse rotation
        arrow_vector = np.array([a_length, 0, 0])
        orientation_rotation = self.particle.orientation_rotations[orientation]
        arrow_vector = orientation_rotation.inv().apply(arrow_vector)
        # FancyArrow takes the arrow base, so centre the arrow on the particle
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
        squared: bool = False,
    ):
        """
        Plots the edge shared with the neighbour along `bond`, colored by the contact type.
        Used to visualise crystalline domains and the defect lines between them.
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
                [0, 0],
            ]
        )
        ax.plot(particles_face_corners[0, :], particles_face_corners[1, :], color=color)
        return

    # ----- PLOTTING SIMULATION RESULTS -----
    def plot_result_outlines(
        self,
        results: np.ndarray | None = None,
        *,
        results_index: int | None = None,
        results_file: str | Path | None = None,
        results_folder: str | Path | None = None,
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
        - `ax` is an optional Axes to draw onto; if None, a new figure is created.
        - `squared` is a boolean which, if set to True, will use the periodic boundary
          conditions to wrap the lattice into a square window rather.
        - Any additional keyword arguments will be passed to matplotlib to create the figure.
        """
        fig, ax = self._resolve_fig_ax(ax, **kwargs)
        results = cfg.load_structure(results_index, results_folder, results_file)
        for site in cfg.get_full_sites(results):
            x_lattice, y_lattice, _ = self.lattice.lattice_site_to_lattice_coords(site)
            self.plot_particle_outline(x_lattice, y_lattice, ax, squared=squared)

        return fig, ax

    def plot_results_arrows(
        self,
        results_index: int | None = None,
        results_folder: str | Path | None = None,
        results_file: str | Path | None = None,
        ax: Axes | None = None,
        squared: bool = False,
        **kwargs,
    ) -> tuple[Figure, Axes]:
        """
        Plots the particles contained in a structure file, with an arrow indicating the
        orientation of the particle (pointing towards face 0).
        The color of the arrow depends on the particle species, in an order given by the
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
        - `ax` is an optional Axes to draw onto; if None, a new figure is created.
        - `squared` is a boolean which, if set to True, will use the periodic boundary
          conditions to wrap the lattice into a square window rather.
        - Any additional keyword arguments will be passed to matplotlib to create the figure.
        """
        fig, ax = self._resolve_fig_ax(ax, **kwargs)
        results = cfg.load_structure(results_index, results_folder, results_file)

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
            y_max = self.lattice.ly * (np.sqrt(3) / 2) + self.radius
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
        contacts: list[tuple[int, int]],
        color: str,
        results_index: int | None = None,
        results_folder: str | Path = "",
        results_file: str | Path | None = None,
        squared=False,
    ):
        """
        Plots the contacts between particles to visualise crystalline domains and defect lines.
        Needs to be put on top of an existing `ax` Axes object.

        `contact_to_color` is a mappable mapping a pair of faces to a color numpy accepts
        (typically a color name or a hex code).

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
        """
        if results is None:
            results = cfg.load_structure(results_index, results_folder, results_file)

        for site, _, orientation in cfg.get_full_sites_characteristics(results):
            x_1, y_1, _ = self.lattice.lattice_site_to_lattice_coords(site)
            neighbours = self.lattice.get_neighbour_sites(site)

            for bond, neighbour in enumerate(neighbours):
                neighbour_orientation = results[1, neighbour]
                face_1, face_2 = self.particle.get_faces_in_contact(
                    orientation, neighbour_orientation, bond
                )
                if (face_1, face_2) in contacts or (face_2, face_1) in contacts:
                    self.plot_contact(x_1, y_1, bond, ax, color, squared)
        return

    def plot_other_contacts(
        self,
        ax: Axes,
        contacts: list[tuple[int, int]],
        color: str,
        results_index: int | None = None,
        results_folder: str | Path = "",
        results_file: str | Path | None = None,
        squared=False,
    ):
        results = cfg.load_structure(results_index, results_folder, results_file)

        for site, _, orientation in cfg.get_full_sites_characteristics(results):
            x_1, y_1, _ = self.lattice.lattice_site_to_lattice_coords(site)
            neighbours = self.lattice.get_neighbour_sites(site)

            for bond, neighbour in enumerate(neighbours):
                neighbour_orientation = results[1, neighbour]
                face_1, face_2 = self.particle.get_faces_in_contact(
                    orientation, neighbour_orientation, bond
                )
                not_in_contacts = (face_1, face_2) not in contacts and (
                    face_2,
                    face_1,
                ) not in contacts
                not_empty = neighbour_orientation != -1
                if not_in_contacts and not_empty:
                    self.plot_contact(x_1, y_1, bond, ax, color, squared)
        return

    def plot_contacts_w_empty(
        self,
        ax: Axes,
        color: str,
        results_index: int | None = None,
        results_folder: str | Path = "",
        results_file: str | Path | None = None,
        squared=False,
    ):
        results = cfg.load_structure(results_index, results_folder, results_file)

        for site, _, orientation in cfg.get_full_sites_characteristics(results):
            x_1, y_1, _ = self.lattice.lattice_site_to_lattice_coords(site)
            neighbours = self.lattice.get_neighbour_sites(site)

            for bond, neighbour in enumerate(neighbours):
                neighbour_orientation = results[1, neighbour]
                is_empty = neighbour_orientation == -1
                if is_empty:
                    self.plot_contact(x_1, y_1, bond, ax, color, squared)
        return

    def square_coordinates(self, x_cartesian):
        # Wrap x into [0, lx): shear offsets top rows by up to 0.5*ly, so use full period
        return np.mod(x_cartesian, self.lattice.lx)
