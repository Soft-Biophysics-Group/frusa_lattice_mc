"""
Vincent Ouazan-Reboul, 2025
Functions to plot the results of simulations of cubic particles.
So far only supports one type of cubes.
"""

# pyright: basic

import config as cfg
from pathlib import Path
from geometry.cubic import CubicGeometry
import bpy
import shutil
import sys
import mathutils
import numpy as np
from geometry.cubic import CubicGeometry

path_to_config = Path(__file__).parent.parent

# Getting blender executable
blender_bin = shutil.which("Blender")
if blender_bin:
    print("Found:", blender_bin)
    bpy.app.binary_path = blender_bin
else:
    print("Unable to find blender!")

# Start from an empty scene
bpy.ops.wm.read_factory_settings(use_empty=True)

# Make sure we get the path to the cube
path_to_numbered_cube = (
    Path(__file__).parent / "assets/oneCube/one_cube_numbered.obj"
).resolve()


def plot_cube(
    coords=[0.0, 0.0, 0.0],
    euler_angles_xyz=[0.0, 0.0, 0.0],
    path_to_cube: Path | str = path_to_numbered_cube,
):
    """
    Puts an object (typically a decorated cube) located at path_to_cube at coordinates coords.
    Euler angles should be given in degrees.
    """
    bpy.ops.wm.obj_import(filepath=str(path_to_cube))
    cube = bpy.context.selected_objects[0]
    cube.location = coords
    new_rotation_euler = mathutils.Euler(euler_angles_xyz)
    cube.rotation_euler.rotate(new_rotation_euler)

    return cube


def plot_cube_from_site_orientation(
    site_index: int, orientation: int, lx: int, ly: int
):
    # The factor of 2 is due to the size of the cubes
    cubic_lattice = CubicGeometry(lx=lx, ly=ly).lattice
    site_coords = 2 * cubic_lattice.lattice_site_to_lattice_coords(site_index)
    rotation = CubicGeometry.particle.orientation_rotations[orientation]
    euler_angles = rotation.as_euler("xyz")

    plot_cube(site_coords, euler_angles)
    return

def save_blender_fig(blend_file_path: str | Path):
    bpy.ops.wm.save_as_mainfile(filepath=str(blend_file_path))

def deleteAllObjects():
    """
    Deletes all objects in the current scene
    """
    deleteListObjects = [
        "MESH",
        "CURVE",
        "SURFACE",
        "META",
        "FONT",
        "HAIR",
        "POINTCLOUD",
        "VOLUME",
        "GPENCIL",
        "ARMATURE",
        "LATTICE",
        "EMPTY",
        "LIGHT",
        "LIGHT_PROBE",
        "CAMERA",
        "SPEAKER",
    ]

    # Select all objects in the scene to be deleted:

    for o in bpy.context.scene.objects:
        for i in deleteListObjects:
            if o.type == i:
                o.select_set(False)
            else:
                o.select_set(True)
    # Deletes all selected objects in the scene:

    bpy.ops.object.delete()


def plot_cubes_from_simulation_results(
    struct_index=-1,
    struct_folder="",
    struct_file="",
    model_file=cfg.default_model_params_file,
    fig_file = ""
):
    # Start by purging all objects we have in memory, otherwise we will plot bullshit
    bpy.ops.object.select_all(action = 'SELECT')
    bpy.ops.object.delete()

    params = cfg.load_model_file(model_file)
    lx = params["lx"]
    ly = params["ly"]
    results = cfg.load_structure(
        struct_index=struct_index, struct_folder=struct_folder, struct_file=struct_file
    )
    orientations = results[1, :]
    full_sites = cfg.get_full_sites(results)
    print(len(full_sites))

    for site in full_sites:
        orientation = orientations[site]
        plot_cube_from_site_orientation(site, orientation, lx, ly)

    if fig_file:
        bpy.ops.wm.save_as_mainfile(filepath=fig_file)

    return


# plot_cube(euler_angles_xyz=[1.57, 0, 0])
# plot_cube_from_site_orientation(1, 1, 10, 10)

# cube_struct_path = cfg.python_path / "examples/cubic/data/01_crystal/structures"
# plot_cubes_from_simulation_results(struct_folder=cube_struct_path)
# blend_file_path = "my_scene.blend"  # Replace with your desired path
# bpy.ops.wm.save_as_mainfile(filepath=blend_file_path)
