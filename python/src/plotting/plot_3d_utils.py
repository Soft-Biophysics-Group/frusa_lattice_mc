"""
Vincent Ouazan-Reboul, 2025
Generic functions to plot the simulation results of 3D particles using Blender.
"""

import config as cfg
from pathlib import Path
import bpy
import shutil

from geometry import LatticeGeometry, ParticleGeometry, get_particle_lattice
import mathutils
import numpy as np

from collections.abc import Callable
from typing import Any

from . import plot_cubic, plot_fcc

DEFAULT_MATERIAL = (14, 0, 255, 1.0)

class BlenderPlot:
    _lattices: dict[str, type["BlenderPlot"]] = {}

    def __init__(
        self,
        particle_geometry: ParticleGeometry,
        lattice_geometry: LatticeGeometry,
        particle_model_file: str | Path,
        plot_boundary_method,
    ) -> None:
        # Look for a Blender executable
        blender_bin = shutil.which("Blender")
        if blender_bin:
            print("Found:", blender_bin)
            bpy.app.binary_path = blender_bin
        else:
            print("Unable to find blender!")

        self.particle_geometry = particle_geometry
        self.lattice_geometry = lattice_geometry
        self.particle_model_file = particle_model_file
        self.plot_grain_boundary_method = plot_boundary_method

        # Load the particle and its material
        bpy.ops.file.pack_all()
        bpy.ops.wm.obj_import(filepath=str(self.particle_model_file))
        self.obj = bpy.context.selected_objects[0]
        self.original_materials = list(obj.data.materials)

        return

    @classmethod
    def from_model_file(
        cls,
        model_file: str | Path = cfg.default_model_params_file,
        particle_model_file: str | Path = "",
    ):
        model_params = cfg.load_model_file(model_file)
        lattice_name = model_params["lattice_name"]
        particle_geometry, lattice_geometry = get_particle_lattice(model_params)

        # Pick the right boundary method
        if lattice_name == "cubic":
            plot_boundary_method = plot_cubic.plot_grain_boundary
        else:
            raise RuntimeError("No boundary plotting method implemented for this lattice")

        return cls(
            particle_geometry,
            lattice_geometry,
            particle_model_file,
            plot_boundary_method,
        )

    def place_obj_copy_from_site_orientation(self, site_index: int, orientation: int):
        # The factor of 2 is due to the size of the cubes
        site_coords = self.lattice_geometry.lattice_site_to_lattice_coords(site_index)
        rotation = self.particle_geometry.orientation_rotations[orientation]
        euler_angles = rotation.as_euler("xyz")

        duplicate_shift_rotate_obj(
            self.obj, self.original_materials, site_coords, euler_angles
        )

        return

    def plot_particles_from_simulation_results(
        self,
        struct_index: int = -1,
        struct_folder: Path | str = "",
        struct_file: Path | str = "",
        fig_file: Path | str = "",
        path_to_obj: Path | str = "",
    ):
        results = cfg.load_structure(
            struct_index=struct_index,
            struct_folder=struct_folder,
            struct_file=struct_file,
        )
        orientations = results[1, :]
        full_sites = cfg.get_full_sites(results)

        for site in full_sites:
            orientation = orientations[site]
            self.place_obj_copy_from_site_orientation(site, orientation)

        if fig_file:
            bpy.ops.wm.save_as_mainfile(filepath=fig_file)

        # Delete original cube
        bpy.data.objects.remove(self.obj, do_unlink=True)

        return

    # def plot_all_boundaries(
    #     self,
    #     boundary_contacts: list[tuple[int, int]],
    #     material=create_material(),
    #     struct_index=-1,
    #     struct_folder="",
    #     struct_file="",
    #     collection_name="Grain boundaries",
    # ):
    #     results = cfg.load_structure(
    #         struct_index=struct_index,
    #         struct_folder=struct_folder,
    #         struct_file=struct_file,
    #     )
    #     orientations = results[1, :]
    #     full_sites = cfg.get_full_sites(results)
    #     for site_1 in full_sites:
    #         orientation_1 = orientations[site_1]
    #         for site_2 in self.lattice_geometry.get_neighbour_sites(site_1):
    #             orientation_2 = orientations[site_2]
    #             self.plot_grain_boundary_method(
    #                 site_1,
    #                 site_2,
    #                 orientation_1,
    #                 orientation_2,
    #                 lattice,
    #                 boundary_contacts,
    #                 material,
    #                 collection_name=collection_name,
    #             )
    #
    #     return

    def clear(self):
        bpy.ops.object.select_all(action="SELECT")
        bpy.ops.object.delete()

        return

    def save(self, save_file):
        bpy.ops.wm.save_as_mainfile(filepath=str(save_file))


def duplicate_shift_rotate_obj(
    obj,
    original_materials,
    coords=[0.0, 0.0, 0.0],
    euler_angles_xyz=[0.0, 0.0, 0.0],
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
    new_rotation_euler = mathutils.Euler(euler_angles_xyz)
    obj_copy.rotation_euler.rotate(new_rotation_euler)

    obj_copy.data.materials.clear()
    for mat in original_materials:
        obj_copy.data.materials.append(mat)  # Assign each material

    collection.objects.link(obj_copy)

    return obj


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
