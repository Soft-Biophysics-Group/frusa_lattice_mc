"""Vincent Ouazan-Reboul, 2025
Generic functions to plot the simulation results of 3D particles using Blender.
"""

import config as cfg
from pathlib import Path
import bpy
import bmesh
import shutil

from geometry import LatticeGeometry, ParticleGeometry, get_particle_lattice
import mathutils
import numpy as np

from collections.abc import Callable
from typing import Any
from numpy.typing import NDArray

from . import cubic, fcc

from scipy.spatial.transform import Rotation as R

from designs.fcc import CRYSTAL_INTS_CORNER

from .build_grain_boundary_geo_nodes import geometry_nodes_node_group
from .build_composition_nodes import build_composition_nodes

from analysis.analyze_3d_size import get_aggregates

DEFAULT_MATERIAL = (14, 0, 255, 1.0)

class BlenderPlot:
    """Class to generate .blend files from simulation results.

    Handles loading the geometry associated with a particular lattice, loading the model of the
    corresponding particle, plotting the simulation results from an `..._model_file.json` and a
    `structure.dat` files.

    For a given lattice, if you want to list the different particle models, create a class
    instance using `bp = BlenderPlot.from_lattice_name(lattice_name)`, and then
    `print(bp.particle_paths)`.

    Don't forget to add new particle paths as you create new models!

    Typical use:
    ```python
    bp = BlenderPlot.from_model_file(model_file)
    bp.load_particle(path_to_particle)
    ```

    To make a simple plot from simulation results, use the
    `plot_particles_from_simulation_results` method.

    For a fancier, "camembert_style" plot where you plot grain boundaries for different
    contacts, use the `plot_canonical_style` method
    """
    _lattices: dict[str, type["BlenderPlot"]] = {}

    def __init__(
        self,
        particle_geometry: ParticleGeometry,
        lattice_geometry: LatticeGeometry,
        gen_face_outline_method,
        particle_paths: dict[str, Path],
        default_particle_path: Path | None = None,
    ) -> None:
        # Look for a Blender executable
        blender_bin = shutil.which("Blender")
        if blender_bin:
            print("Found:", blender_bin)
            bpy.app.binary_path = blender_bin
        else:
            print("Unable to find blender!")
        self.clear()

        self.particle_geometry = particle_geometry
        self.lattice_geometry = lattice_geometry
        self.plot_grain_boundary_method = gen_face_outline_method

        # self.load_particle(default_particle_path)


        self.particle_paths = particle_paths
        if default_particle_path is None:
            default_particle_path = list(self.particle_paths.values())[0]
        self.default_particle_path = default_particle_path

        self.particle_coords_cart = {}
        self.particle_objs = {}

        self.orientations: NDArray[np.int_]
        self.full_sites: NDArray[np.int_]

        self.fig_file: str | Path

        self.model_file : str | Path
        self.n_aggregates: int

        return

    @classmethod
    def from_lattice_name(
        cls,
        lattice_name: str,
        lx: int,
        ly: int,
        lz: int,
    ):
        particle_geometry = ParticleGeometry.from_lattice_name(lattice_name)
        lattice_geometry = LatticeGeometry.from_lattice_name(lattice_name, lx, ly, lz)
        boundary_plot_method = load_boundary_method(lattice_name)
        particle_paths, default_particle_path = import_particle_paths(lattice_name)


        return cls(
            particle_geometry,
            lattice_geometry,
            boundary_plot_method,
            particle_paths,
            default_particle_path,
        )

    @classmethod
    def from_model_file(
        cls,
        model_file: str | Path = cfg.default_model_params_file,
    ):
        model_params = cfg.load_model_file(Path(model_file))
        lattice_name = model_params["lattice_name"]
        lx = model_params["lx"]
        ly = model_params["ly"]
        lz = model_params["lz"]


        this_instance =  cls.from_lattice_name(lattice_name, lx, ly, lz)
        this_instance.model_file = model_file

        return this_instance

    def load_particle(self, path_to_particle: Path | None = None):
        # Load the particle and its material
        if path_to_particle is None:
            path_to_particle = self.default_particle_path

        bpy.ops.file.pack_all()
        bpy.ops.wm.obj_import(filepath=str(path_to_particle))
        self.obj = bpy.context.selected_objects[0]
        # Hide the original object from view and selection
        self.obj.hide_viewport = True
        self.obj.hide_select = True
        self.original_materials = list(self.obj.data.materials)

        return

    def cubify(self, site_coords_cartesian):
        x, y, z = (
            site_coords_cartesian[0],
            site_coords_cartesian[1],
            site_coords_cartesian[2],
        )
        cube_lx = 2 ** (1 / 3) * self.lattice_geometry.lx
        cube_ly = 2 ** (1 / 3) * self.lattice_geometry.ly
        cube_lz = 2 ** (1 / 3) * self.lattice_geometry.lz

        if x >= cube_lx:
            x -= cube_lx
        if y >= cube_ly:
            y -= cube_ly
        if z >= cube_lz:
            z -= cube_lz

        return np.array([x, y, z])

    def calc_barycenter_coords_cartesian(self, coords):
        lx_cart, ly_cart, lz_cart = (
            self.lattice_geometry.lx * self.lattice_geometry.lattice_spacing,
            self.lattice_geometry.ly * self.lattice_geometry.lattice_spacing,
            self.lattice_geometry.lz * self.lattice_geometry.lattice_spacing
        )

        theta_coords = np.zeros_like(coords)
        theta_coords[:, 0] = 2 * np.pi * coords[:, 0] / lx_cart
        theta_coords[:, 1] = 2 * np.pi * coords[:, 1] / ly_cart
        theta_coords[:, 2] = 2 * np.pi * coords[:, 2] / lz_cart

        complex_coords = np.exp(1j * theta_coords)
        barycenter_complex_coords = np.mean(complex_coords, axis=0)

        # Go back to (non-integer) barycenter coordinates
        x_barycenter = (
            lx_cart / (2 * np.pi) * np.angle(barycenter_complex_coords[0])
        ) % lx_cart
        y_barycenter = (
            ly_cart / (2 * np.pi) * np.angle(barycenter_complex_coords[1])
        ) % ly_cart
        z_barycenter = (
            lz_cart / (2 * np.pi) * np.angle(barycenter_complex_coords[2])
        ) % lz_cart


        return x_barycenter, y_barycenter, z_barycenter

    def draw_line_between_two_points(
        self,
        point1: tuple[float, float, float],
        point2: tuple[float, float, float],
        color: tuple[int, int, int, int] = (0, 0, 0, 1),
        strength: int = 1,
        thickness: float = 0.1
    ):
        # Code generated using chatgpt. Use with caution.

        # Create a new curve data block
        curve_collection = get_or_create_collection("Box")
        curve_data = bpy.data.curves.new(name="BoxEdge", type='CURVE')
        curve_data.dimensions = '3D'

        # Create a polyline spline
        spline = curve_data.splines.new(type='POLY')
        spline.points.add(1)  # Already has 1 point, needs to join 2

        for i, pt in enumerate((point1, point2)):
            spline.points[i].co = (*pt, 1)  # (x, y, z, w)

        # Add bevel for visible thickness
        curve_data.bevel_depth = thickness      # Adjust thickness
        curve_data.bevel_resolution = 4

        # Create associated object and associated material we want
        curve_obj = bpy.data.objects.new("BoxLine", curve_data)

        collection = get_or_create_collection("Box")
        collection.objects.link(curve_obj)

        mat = create_line_material(color, strength)
        curve_obj.data.materials.append(mat)

        return curve_obj

    def draw_cube(
        self,
        cube_length: float | None = None,
        color: tuple[int, int, int, int] | None = None,
        strength=1,
        thickness=0.2,
    ):
        if color is None:
            color = (0, 0, 0, 1)
        if cube_length is None:
            lx_cart, ly_cart, lz_cart = (
                self.lattice_geometry.lx * self.lattice_geometry.lattice_spacing,
                self.lattice_geometry.ly * self.lattice_geometry.lattice_spacing,
                self.lattice_geometry.lz * self.lattice_geometry.lattice_spacing
            )
            cube_edges = [
                ((0, 0, 0), (lx_cart, 0, 0)),
                ((0, 0, 0), (0, ly_cart, 0)),
                ((0, 0, 0), (0, 0, lz_cart)),
                ((lx_cart, ly_cart, 0), (0, ly_cart, 0)),
                ((lx_cart, ly_cart, 0), (lx_cart, 0, 0)),
                ((lx_cart, ly_cart, 0), (lx_cart, ly_cart, lz_cart)),
                ((lx_cart, 0, lz_cart), (lx_cart, ly_cart, lz_cart)),
                ((lx_cart, 0, lz_cart), (0, 0, lz_cart)),
                ((lx_cart, 0, lz_cart), (lx_cart, 0, 0)),
                ((0, ly_cart, lz_cart), (lx_cart, ly_cart, lz_cart)),
                ((0, ly_cart, lz_cart), (0, 0, lz_cart)),
                ((0, ly_cart, lz_cart), (0, ly_cart, 0)),
            ]

        else:
            clp1 = cube_length + 1
            cube_edges = [
                ((-1  , -1  , -1)  , (clp1, -1  , -1))  ,
                ((-1  , -1  , -1)  , (-1  , clp1, -1))  ,
                ((-1  , -1  , -1)  , (-1  , -1  , clp1)),
                ((clp1, clp1, -1)  , (-1  , clp1, -1))  ,
                ((clp1, clp1, -1)  , (clp1, -1  , -1))  ,
                ((clp1, clp1, -1)  , (clp1, clp1, clp1)),
                ((clp1, -1  , clp1), (clp1, clp1, clp1)),
                ((clp1, -1  , clp1), (-1  , -1  , clp1)),
                ((clp1, -1  , clp1), (clp1, -1  , -1))  ,
                ((-1  , clp1, clp1), (clp1, clp1, clp1)),
                ((-1  , clp1, clp1), (-1  , -1  , clp1)),
                ((-1  , clp1, clp1), (-1  , clp1, -1))  ,
            ]

        for edge_point_1, edge_point_2 in cube_edges:
            self.draw_line_between_two_points(
                edge_point_1,
                edge_point_2,
                color=color,
                strength=strength,
                thickness=thickness,
            )

        return

    def center_cubic(self, site_coords_cartesian, objs):
        barycenter_cartesian = self.calc_barycenter_coords_cartesian(
            site_coords_cartesian
        )
        lx_cart, ly_cart, lz_cart = (
            self.lattice_geometry.lx * self.lattice_geometry.lattice_spacing,
            self.lattice_geometry.ly * self.lattice_geometry.lattice_spacing,
            self.lattice_geometry.lz * self.lattice_geometry.lattice_spacing
        )
        # Astype(int) is so that we get integer shift, which makes the figure cleaner
        shift_vector = -np.array(barycenter_cartesian).astype(int) + np.array(
            [lx_cart // 2, ly_cart // 2, lz_cart // 2]
        )

        for obj, site_coord in zip(objs, site_coords_cartesian):
            obj.location = (np.array(obj.location) + shift_vector) % np.array(
                [lx_cart, ly_cart, lz_cart]
            )

        return


    def place_obj_copy_from_site_orientation(
        self,
        site_index: int,
        orientation: int,
        barycenter: tuple[int, int, int] | None = None,
        offset_cartesian: list[int] | NDArray[np.int_] | None = None,
        offset_lattice: list[float] | NDArray[np.float_] | None = None
    ):
        lattice_coords = self.lattice_geometry.lattice_site_to_lattice_coords(
            site_index
        )
        if barycenter is not None:
            lattice_coords -= barycenter - np.array(
                [
                    self.lattice_geometry.lx//2,
                    self.lattice_geometry.ly//2,
                    self.lattice_geometry.lz//2,
                ]
            )
            lattice_coords = self.lattice_geometry.apply_pbc(*lattice_coords)

        if offset_lattice is not None:
            lattice_coords = np.array(lattice_coords) + offset_lattice

        site_coords_cartesian = self.lattice_geometry.lattice_to_cartesian(
            *lattice_coords
        )

        if offset_cartesian is not None:
            site_coords_cartesian += offset_cartesian

        rotation = self.particle_geometry.orientation_rotations[orientation]
        quaternion_coeffs = rotation.as_quat(canonical=False, scalar_first=True)
        euler_xyz_angles = rotation.as_euler("xyz", degrees = True)

        obj_copy = duplicate_shift_rotate_obj(
            self.obj, self.original_materials, site_coords_cartesian, quaternion_coeffs
        )

        return obj_copy, site_coords_cartesian


    def calc_barycenter_coords(self, full_sites):
        # Done using the circular mean, https://en.wikipedia.org/wiki/Circular_mean

        # Normalize the coordinates to [0, 2pi] and take their complex exponential
        barycenter_coords_theta = np.zeros(3, dtype=complex)
        for site in full_sites:
            x_lattice, y_lattice, z_lattice = (
                    self.lattice_geometry.lattice_site_to_lattice_coords(site)
            )
            theta_x = 2 * np.pi * x_lattice / self.lattice_geometry.lx
            theta_y = 2 * np.pi * y_lattice / self.lattice_geometry.ly
            theta_z = 2 * np.pi * z_lattice / self.lattice_geometry.lz
            barycenter_coords_theta[0] += np.exp(1j * theta_x)
            barycenter_coords_theta[1] += np.exp(1j * theta_y)
            barycenter_coords_theta[2] += np.exp(1j * theta_z)

        barycenter_coords_theta /= len(full_sites)

        # Go back to (non-integer) barycenter coordinates
        x_bar_nonint = (
            self.lattice_geometry.lx
            / (2 * np.pi)
            * np.angle(barycenter_coords_theta[0])
        ) % self.lattice_geometry.lx
        y_bar_nonint = (
            self.lattice_geometry.ly
            / (2 * np.pi)
            * np.angle(barycenter_coords_theta[1])
        ) % self.lattice_geometry.ly
        z_bar_nonint = (
            self.lattice_geometry.lz
            / (2 * np.pi)
            * np.angle(barycenter_coords_theta[2])
        ) % self.lattice_geometry.lz

        x_barycenter = int(np.floor(x_bar_nonint))
        y_barycenter = int(np.floor(y_bar_nonint))
        z_barycenter = int(np.floor(z_bar_nonint))

        return x_barycenter, y_barycenter, z_barycenter

    def draw_fcc_cell(
        self, color: tuple[int, int, int, int] | None = None, strength=1, thickness=0.1
    ):
        if color is None:
            color = (0, 0, 0, 1)
        # For brevity
        lx, ly, lz = (
            self.lattice_geometry.lx,
            self.lattice_geometry.ly,
            self.lattice_geometry.lz,
        )
        edges_lattice = [
            ((0, 0, 0), (lx - 1, 0, 0)),
            ((0, 0, 0), (0, ly - 1, 0)),
            ((0, 0, 0), (0, 0, lz - 1)),
            ((lx - 1, ly - 1, 0), (0, ly - 1, 0)),
            ((lx - 1, ly - 1, 0), (lx - 1, 0, 0)),
            ((lx - 1, ly - 1, 0), (lx - 1, ly - 1, lz - 1)),
            ((lx - 1, 0, lz - 1), (lx - 1, ly - 1, lz - 1)),
            ((lx - 1, 0, lz - 1), (0, 0, lz - 1)),
            ((lx - 1, 0, lz - 1), (lx - 1, 0, 0)),
            ((0, ly - 1, lz - 1), (lx - 1, ly - 1, lz - 1)),
            ((0, ly - 1, lz - 1), (0, 0, lz - 1)),
            ((0, ly - 1, lz - 1), (0, ly - 1, 0)),
        ]

        edge_cartesian = [
            (
                tuple(self.lattice_geometry.lattice_to_cartesian(*point1)),
                tuple(self.lattice_geometry.lattice_to_cartesian(*point2)),
            )
            for point1, point2 in edges_lattice
        ]

        for edge_point_1, edge_point_2 in edge_cartesian:
            self.draw_line_between_two_points(
                edge_point_1,
                edge_point_2,
                color=color,
                strength=strength,
                thickness=thickness,
            )

        return

    def plot_particles_from_simulation_results(
        self,
        struct_index: int = -1,
        struct_folder: Path | str = "",
        struct_file: Path | str = "",
        fig_file: Path | str = "",
        centered: bool = False,
        cube_length: float | None = None,
        draw_box: bool = False,
        offset_cartesian: list[float] | NDArray[np.float_] | None = None,
        box_color: tuple[int, int, int, int] | None = None,
        apply_depth_effect: bool = True
    ):
        results = cfg.load_structure(
            struct_index=struct_index,
            struct_folder=struct_folder,
            struct_file=struct_file,
        )
        self.orientations = results[1, :]
        self.full_sites = cfg.get_full_sites(results)

        self.fig_file = fig_file

        if offset_cartesian is None:
            offset_cartesian = np.zeros(3)

        center_offset = np.zeros(3)
        if centered:
            barycenter_coords = self.calc_barycenter_coords(self.full_sites)
            if cube_length is not None:
                barycenter_after_centering_cart = np.array(
                    self.lattice_geometry.lattice_to_cartesian(
                        self.lattice_geometry.lx // 2,
                        self.lattice_geometry.ly // 2,
                        self.lattice_geometry.lz // 2,
                    )
                )
                new_center_cart = np.array(
                    [cube_length / 2, cube_length / 2, cube_length / 2]
                )
                center_offset = new_center_cart - barycenter_after_centering_cart

        if cube_length is not None:
            all_periodic_copies_vectors = np.vstack(
                [
                    [0.0, 0.0, 0.0],
                    self.lattice_geometry.lattice_vectors_in_cartesian,
                    -self.lattice_geometry.lattice_vectors_in_cartesian,
                    self.lattice_geometry.lattice_vectors_in_cartesian[:, 0]
                    - self.lattice_geometry.lattice_vectors_in_cartesian[:, 1],
                    self.lattice_geometry.lattice_vectors_in_cartesian[:, 0]
                    - self.lattice_geometry.lattice_vectors_in_cartesian[:, 2],
                    self.lattice_geometry.lattice_vectors_in_cartesian[:, 1]
                    - self.lattice_geometry.lattice_vectors_in_cartesian[:, 2],
                    self.lattice_geometry.lattice_vectors_in_cartesian[:, 1]
                    - self.lattice_geometry.lattice_vectors_in_cartesian[:, 0],
                    self.lattice_geometry.lattice_vectors_in_cartesian[:, 2]
                    - self.lattice_geometry.lattice_vectors_in_cartesian[:, 0],
                    self.lattice_geometry.lattice_vectors_in_cartesian[:, 2]
                    - self.lattice_geometry.lattice_vectors_in_cartesian[:, 1],
                ]
            )
            all_periodic_copies_vectors *= np.array(
                [
                    self.lattice_geometry.lx,
                    self.lattice_geometry.ly,
                    self.lattice_geometry.lz,
                ]
            ).T
        else:
            all_periodic_copies_vectors = np.array([[0.0, 0.0, 0.0]])

        barycenter_cartesian = np.zeros(3)
        n_particles = 0
        # Original cell
        for site in self.full_sites:
            orientation = self.orientations[site]
            self.particle_coords_cart[site] = []
            self.particle_objs[site] = []
            for periodic_vec in all_periodic_copies_vectors:
                obj_copy, site_coords_cartesian = (
                    self.place_obj_copy_from_site_orientation(
                        site,
                        orientation,
                        barycenter=barycenter_coords,
                        offset_cartesian=periodic_vec
                        + center_offset
                        + offset_cartesian,
                    )
                )
                if cube_length is not None and (
                    (site_coords_cartesian < 0).any()
                    or (site_coords_cartesian >= cube_length).any()
                ):
                    bpy.data.objects.remove(obj_copy, do_unlink=True)
                else:
                    n_particles += 1
                    self.particle_coords_cart[site].append(site_coords_cartesian)
                    self.particle_objs[site].append(obj_copy)
                    barycenter_cartesian += site_coords_cartesian
        barycenter_cartesian /= n_particles

        # Center the cubified simulation result
        if (cube_length is not None) and centered:
            new_barycenter_cartesian = np.zeros(3) + cube_length / 2
            offset_vector = new_barycenter_cartesian - barycenter_cartesian
            for obj_list in self.particle_objs.values():
                for obj in obj_list:
                    for i in range(3):
                        obj.location[i] += offset_vector[i]
                        obj.location[i] -= cube_length * (
                            obj.location[i] >= cube_length
                        )
                    # barycenter_cartesian += obj.location
        # barycenter_cartesian /= n_particles

        if draw_box:
            if cube_length is not None:
                self.draw_cube(cube_length, box_color)
            else:
                self.draw_fcc_cell(box_color)

        if fig_file:
            bpy.ops.wm.save_as_mainfile(filepath=self.fig_file)

        # Delete original cube: we most likely won't do anything with it anymore
        # bpy.data.objects.remove(self.obj, do_unlink=True)

        self.obj.hide_viewport = True
        self.obj.hide_render = True

        return


    def plot_particles_from_simulation_results_w_periodic_copies(
        self,
        struct_index: int = -1,
        struct_folder: Path | str = "",
        struct_file: Path | str = "",
        fig_file: Path | str = "",
        centered: bool = False,
        cube_length: float | None = None,
        draw_box: bool = False,
        offset_cartesian: list[float] | NDArray[np.float_] | None = None,
        box_color: tuple[int, int, int, int] | None = None
    ):
        results = cfg.load_structure(
            struct_index=struct_index,
            struct_folder=struct_folder,
            struct_file=struct_file,
        )
        self.orientations = results[1, :]
        self.full_sites = cfg.get_full_sites(results)

        self.fig_file = fig_file

        center_offset = np.zeros(3)
        if centered:
            barycenter_coords = self.calc_barycenter_coords(self.full_sites)
        all_periodic_copies_vectors = np.vstack(
            [
                [0.0, 0.0, 0.0],
                self.lattice_geometry.lattice_vectors_in_cartesian,
                -self.lattice_geometry.lattice_vectors_in_cartesian,
                self.lattice_geometry.lattice_vectors_in_cartesian[:, 0]
                + self.lattice_geometry.lattice_vectors_in_cartesian[:, 1],
                self.lattice_geometry.lattice_vectors_in_cartesian[:, 1]
                + self.lattice_geometry.lattice_vectors_in_cartesian[:, 2],
                self.lattice_geometry.lattice_vectors_in_cartesian[:, 2]
                + self.lattice_geometry.lattice_vectors_in_cartesian[:, 0],
            ]
        )
        all_periodic_copies_vectors *= np.array(
            [
                self.lattice_geometry.lx,
                self.lattice_geometry.ly,
                self.lattice_geometry.lz,
            ]
        ).T

        barycenter_cartesian = np.zeros(3)
        n_particles = 0
        # Original cell
        for site in self.full_sites:
            orientation = self.orientations[site]
            self.particle_coords_cart[site] = []
            self.particle_objs[site] = []
            for periodic_vec in all_periodic_copies_vectors:
                obj_copy, site_coords_cartesian = (
                    self.place_obj_copy_from_site_orientation(
                        site,
                        orientation,
                        offset_cartesian=periodic_vec + center_offset
                    )
                )
                n_particles += 1
                barycenter_cartesian += site_coords_cartesian
                self.particle_coords_cart[site].append(site_coords_cartesian)
                self.particle_objs[site].append(obj_copy)
        barycenter_cartesian /= n_particles

        if draw_box:
            if cube_length is not None:
                self.draw_cube(cube_length, box_color)
                self.draw_fcc_cell(box_color)
            else:
                self.draw_fcc_cell(box_color)

        if fig_file:
            bpy.ops.wm.save_as_mainfile(filepath=self.fig_file)

        # Delete original cube: we most likely won't do anything with it anymore
        # bpy.data.objects.remove(self.obj, do_unlink=True)

        self.obj.hide_viewport = True
        self.obj.hide_render = True

        return

    def gen_lozenge_boundary_mesh(
        self,
        site,
        bond_orientation,
        particle_center_cartesian_coords,
        boundary_offset=0.05,
    ):
        # Do nothing if the particles don't touch
        if bond_orientation == -1:
            return

        half_size = self.lattice_geometry.lattice_spacing / 2
        extra_offset = boundary_offset
        vertices = np.array(
            [
                [1 + extra_offset, -extra_offset, 0],
                [0.5, 0.5, 0.5 + extra_offset],
                [-extra_offset, 1 + extra_offset, 0],
                [0.5, 0.5, -0.5 - extra_offset],
            ]
        )

        rotation_to_apply = self.particle_geometry.bond_rotations[bond_orientation]
        vertices = rotation_to_apply.inv().apply(vertices)

        # Translate the lozenge to the right spot
        vertices += particle_center_cartesian_coords

        faces = [(0, 1, 2, 3)]
        edges = [(0, 1), (1, 2), (2, 3), (3, 0)]

        # Create a new mesh and object
        mesh = bpy.data.meshes.new(name="SquareMesh")
        mesh.from_pydata(vertices, edges, faces)
        mesh.update()

        return mesh

    def plot_all_boundaries_from_simulation_results(
        self,
        contacts_of_interest: list[tuple[int, int]],
        fig_file: Path | str = "",
        material = None,
        thickness = 0.2,
        boundary_offset = 0.1,
        offset = 0.0,
        use_rim_only = False
    ) -> None:
        all_meshes = []

        for i, site_1 in enumerate(self.full_sites):
            orientation_1 = self.orientations[site_1]
            for obj in self.particle_objs[site_1]:
                coords_cart_1 = obj.location
                neighbours_of_1 = self.lattice_geometry.get_neighbour_sites(site_1)
                full_neighbours_of_1 = set(neighbours_of_1) & set(self.full_sites)
                for neighbour in full_neighbours_of_1:
                    site_2 = neighbour
                    orientation_2 = self.orientations[site_2]
                    face_1, face_2, bond = (
                        self.lattice_geometry.get_faces_in_contact_and_bond(
                            site_1, orientation_1, site_2, orientation_2
                        )
                    )
                    for (
                        alt_face_1,
                        alt_face_2,
                    ) in self.particle_geometry.get_equivalent_face_pairs(
                        face_1, face_2
                    ):
                        if (alt_face_1, alt_face_2) in contacts_of_interest or (
                            alt_face_2,
                            alt_face_1,
                        ) in contacts_of_interest:
                            all_meshes.append(
                                self.gen_lozenge_boundary_mesh(
                                    site_1,
                                    bond,
                                    coords_cart_1,
                                    boundary_offset=boundary_offset,
                                )
                            )

        merge_boundaries_into_bmesh(
            all_meshes,
            thickness=thickness,
            material=material,
            use_rim_only=use_rim_only,
            offset=offset,
        )

        if self.fig_file:
            bpy.ops.wm.save_as_mainfile(filepath=self.fig_file)

        return

    def gen_lozenge_boundary_mesh(
        self,
        site,
        bond_orientation,
        particle_center_cartesian_coords,
        collection_name="Boundaries",
        boundary_offset = 0.05
    ):
        # Do nothing if the particles don't touch
        if bond_orientation == -1:
            return

        collection = get_or_create_collection(collection_name)

        half_size = self.lattice_geometry.lattice_spacing / 2
        extra_offset = boundary_offset
        vertices = np.array(
            [
                [1 + extra_offset, -extra_offset, 0],
                [0.5, 0.5, 0.5 + extra_offset],
                [-extra_offset, 1 + extra_offset, 0],
                [0.5, 0.5, -0.5 - extra_offset],
            ]
        )

        rotation_to_apply = self.particle_geometry.bond_rotations[bond_orientation]
        vertices = rotation_to_apply.inv().apply(vertices)

        # Translate the lozenge to the right spot
        vertices += particle_center_cartesian_coords

        faces = [(0, 1, 2, 3)]
        edges = [(0, 1), (1, 2), (2, 3), (3, 0)]

        # Create a new mesh and object
        mesh = bpy.data.meshes.new(name="SquareMesh")
        # obj = bpy.data.objects.new(name="Square", object_data=mesh)
        # collection.objects.link(obj)
        # bpy.context.collection.objects.link(obj)
        # solidify = obj.modifiers.new(name="Solidify", type='SOLIDIFY')
        # solidify.thickness = 0.2  # Adjust thickness
        # solidify.offset = 0  # 0 = centered, 1 = extrude outward only
        mesh.from_pydata(vertices, edges, faces)
        mesh.update()

        return mesh

    def clear(self):
        bpy.ops.object.select_all(action="SELECT")
        bpy.ops.object.delete()

        return

    def save(self, save_file):
        bpy.ops.wm.save_as_mainfile(filepath=str(save_file))

    def plot_canonical_style(
        self,
        mc_file: str | Path,
        fig_save_file: str | Path,
        defect_contacts: list[tuple[int, int]],
        centered: bool = True,
        cubified: bool = False,
        draw_box: bool = True,
        crystal_contacts: list[tuple[int, int]] | None = None,
        offset_cartesian: list[float] | NDArray[np.float_] | None = None,
        cube_length: float | None = None,
        box_color: tuple[int, int, int, int] | None = None,
        apply_depth_effect: bool = True,
        camera_location: tuple[float, float, float] | None = None,
        camera_rotation_deg: tuple[float, float, float] | None = None,
        domain_line_boundary_offset: float = 0.005
    ):
        mc_params = cfg.load_mc_file(mc_file)
        structure_folder = Path(mc_params["final_structure_address"])
        struct_file = structure_folder / "final_structure.dat"

        if struct_file.is_file():
            print(f"\nPlotting {struct_file}")
            if cubified and cube_length is None:
                cube_length = 2 ** (1 / 3) * self.lattice_geometry.lx
                print(self.lattice_geometry.lx)

            self.plot_particles_from_simulation_results(
                struct_file=struct_file,
                centered=centered,
                draw_box=draw_box,
                offset_cartesian=offset_cartesian,
                cube_length=cube_length,
                box_color=box_color,
            )
            self.plot_all_boundaries_from_simulation_results(
                defect_contacts,
                boundary_offset=domain_line_boundary_offset,
                thickness=0.06,
                offset=0.56,
                use_rim_only=True,
            )
            inner_mat = create_line_material((140 / 255, 175 / 255, 246 / 255, 1.0))
            self.plot_all_boundaries_from_simulation_results(
                defect_contacts,
                boundary_offset=-0.05,
                thickness=0.01,
                material=inner_mat,
            )

            if crystal_contacts is not None:
                outer_mat_crystal = create_line_material(
                    (0.262, 0.604, 0.402, 1.0)
                )
                self.plot_all_boundaries_from_simulation_results(
                    crystal_contacts,
                    boundary_offset=domain_line_boundary_offset,
                    thickness=0.06,
                    offset=0.56,
                    use_rim_only=True,
                    material = outer_mat_crystal
                )
                inner_mat_crystal = create_line_material(
                    (199 / 255, 229 / 255, 215 / 255, 1.0)
                )
                self.plot_all_boundaries_from_simulation_results(
                    crystal_contacts,
                    boundary_offset=-0.05,
                    thickness=0.01,
                    material=inner_mat_crystal,
                )

        # Fluff: use our preferred settings for convenience
        # Enable the cycles engine
        bpy.context.scene.render.engine = "CYCLES"
        bpy.context.preferences.addons[
            "cycles"
        ].preferences.compute_device_type = "METAL"
        bpy.context.scene.cycles.device = "GPU"
        # Make the film transparent
        bpy.context.scene.render.film_transparent = True

        # Make the world background white
        bpy.context.scene.world.use_nodes = True
        bg = bpy.context.scene.world.node_tree.nodes["Background"]
        bg.inputs[0].default_value = (1, 1, 1, 1)  # RGBA white

        # Set the render resolution to a square
        bpy.context.scene.render.resolution_x = 1080

        # Put the camera at our preferred location and angle
        # Create or get camera
        if camera_location is None:
            camera_location = (16.19, 76.6, 11.3)
        if not bpy.data.objects.get("Camera"):
            bpy.ops.object.camera_add(
                location=camera_location
            )  # Create new if not present
        camera = bpy.data.objects["Camera"]

        # Set camera location and rotation
        if camera_rotation_deg is None:
            camera_rotation_deg = (88.1825, 0, 173.5)
        deg_to_rad = 2 * np.pi / 360
        camera_rotation_rad = tuple(
            camera_rotation_deg_axis * deg_to_rad
            for camera_rotation_deg_axis in camera_rotation_deg
        )
        camera.rotation_euler = camera_rotation_rad
        # Put it in orthographic mode
        camera.data.type = "ORTHO"
        camera.data.ortho_scale = 28.8

        # Outline the clusters based on depth
        if apply_depth_effect:
           # Enable Z pass for rendering
            # Change the view layers accordingly
            bpy.context.view_layer.name = "WithBox"
            # Create new layer without the box drawn in
            new_layer = bpy.context.scene.view_layers.new(name="WithoutBox")
            # We finally need to select the box collection, which is a bit of a pain
            for coll in new_layer.layer_collection.children:
                if coll.name == "Box":
                    coll.exclude = True
            # And build the right node setup
            bpy.context.scene.view_layers["WithoutBox"].use_pass_z = True
            build_composition_nodes()

        self.save(fig_save_file)

        return

    def plot_canonical_style_from_mc_file(self):
        return

def import_particle_paths(lattice_name: str):
    if lattice_name == "cubic":
        return cubic.paths, cubic.paths["path_to_numbered_cube"]
    elif lattice_name == "fcc":
        return fcc.paths, fcc.paths["path_to_numbered_rhombic"]
    else:
        raise ValueError("No paths to particle models implemented for this lattice")


def load_boundary_method(lattice_name: str):
    if lattice_name == "cubic":
        return cubic.plot_grain_boundary
    elif lattice_name == "fcc":
        return fcc.plot_lozenge
    else:
        raise RuntimeError("No boundary plotting method implemented for this lattice")

    return


def duplicate_shift_rotate_obj(
    obj,
    original_materials,
    coords=[0.0, 0.0, 0.0],
    quats_coeffs_wxyz=[0.0, 0.0, 0.0, 0.0],
    collection_name="Particles",
):
    """Duplicates a loaded Blender object, placing it at location coords and rotating it using
    the angles euler_angles_xyz.
    Euler angles should be given in degrees.
    """
    collection = get_or_create_collection(collection_name)

    obj_copy = obj.copy()
    obj_copy.data = obj.data

    obj_copy.location = coords
    obj_copy.rotation_mode = "QUATERNION"
    rotation_quat = mathutils.Quaternion(quats_coeffs_wxyz)
    obj_copy.rotation_quaternion.rotate(rotation_quat)
    # Go back to xyz mode because it is easier to interpret
    obj_copy.rotation_mode = "XYZ"


    obj_copy.data.materials.clear()
    for mat in original_materials:
        obj_copy.data.materials.append(mat)  # Assign each material

    # Unhide obj_copy
    obj_copy.hide_viewport = False
    obj_copy.hide_select = False

    collection.objects.link(obj_copy)

    return obj_copy


# Code created with the help of ChatGPT. Use with caution
def get_or_create_collection(collection_name="Squares"):
    """Check if a collection exists, if not, create it."""
    if collection_name in bpy.data.collections:
        return bpy.data.collections[collection_name]

    # Create a new collection and link it to the scene
    new_collection = bpy.data.collections.new(collection_name)
    bpy.context.scene.collection.children.link(new_collection)
    return new_collection


def create_material(name="CustomMaterial", color=DEFAULT_MATERIAL):
    """Create a new material with a given RGBA color (supports transparency).
    Generated by ChatGPT (I know, I know. I was in a rush)"""
    mat = bpy.data.materials.new(name=name)
    mat.use_nodes = True  # Enable node-based material

    # Get the Principled BSDF node
    bsdf = mat.node_tree.nodes.get("Principled BSDF")

    if bsdf:
        bsdf.inputs["Base Color"].default_value = color  # Set color & alpha
        bsdf.inputs["Alpha"].default_value = color[3]  # Set transparency

    # Enable transparency settings
    mat.blend_method = "BLEND"  # Enables alpha blending (smooth transparency)
    # mat.shadow_method = "HASHED"  # Ensures shadows work with transparency
    mat.use_backface_culling = False  # Show both sides

    return mat


def merge_boundaries_into_bmesh(
    boundary_meshes, material=None, thickness=0.2, offset=0.0, use_rim_only=False
):
    # Generated using ChatGPT - beware!

    final_mesh = bpy.data.meshes.new("MergedMesh")
    final_obj = bpy.data.objects.new("Boundaries", final_mesh)

    collection = get_or_create_collection("Boundaries")
    collection.objects.link(final_obj)

    # Build one combined BMesh
    bm_main = bmesh.new()

    for mesh in boundary_meshes:
        # Copy into the main bmesh
        bm_main.from_mesh(mesh)

    bm_main.to_mesh(final_mesh)
    bm_main.free()

    if material is None:
        material = create_line_material(color=(12 / 256, 5 / 256, 216 / 256, 1))
    final_obj.data.materials.append(material)

    solidify = final_obj.modifiers.new(name="Solidify", type="SOLIDIFY")
    solidify.thickness = thickness  # Adjust thickness
    solidify.offset = offset  # 0 = centered, 1 = extrude outward only
    solidify.use_rim_only = use_rim_only

    if use_rim_only:
        geo_nodes = final_obj.modifiers.new(name="GeoModifier", type="NODES")
        geo_nodes.node_group = geometry_nodes_node_group(material.name)

    return


def create_line_material(
    color: tuple[float, float, float, float] = (0, 0, 0, 1),
    strength: int = 1,
    name="UnlitBlack",
):
    # Code generated using ChatGPT. Use with caution.
    # Create unlit black material
    mat = bpy.data.materials.new(name=name)
    mat.use_nodes = True
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links

    # Remove default nodes
    for node in nodes:
        nodes.remove(node)

    # Add emission shader (pure black)
    emission = nodes.new(type="ShaderNodeEmission")
    emission.inputs["Color"].default_value = color
    emission.inputs["Strength"].default_value = strength

    # Output node
    output = nodes.new(type="ShaderNodeOutputMaterial")
    links.new(emission.outputs["Emission"], output.inputs["Surface"])

    return mat


