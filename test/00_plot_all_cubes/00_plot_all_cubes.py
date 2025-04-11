"""
Vincent Ouazan-Reboul, 2025

Plots all of the cube orientation to make sure I got them right.
Should be left to right, orientations 0 to 23, with orientation i corresponding to patch i
replacing patch 0
"""

import plotting.plot_cubic as pc
import bpy

cube, material = pc.load_cube(pc.path_to_numbered_cube)
for orientation in range(24):
    pc.place_cube_copy_from_site_orientation(
        cube, material, orientation * 2, orientation, 25, 25
    )

pc.save_blender_fig("3dFigures/all_orientations.blend")
