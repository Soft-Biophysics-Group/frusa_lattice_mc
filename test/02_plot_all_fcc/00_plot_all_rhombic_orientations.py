"""Vincent Ouazan-Reboul, 2025

Plots all of the rhombic dodecahedra orientation to make sure I got them right.
Should be left to right, orientations 0 to 23, with orientation i corresponding to patch i
replacing patch 0
"""

from plotting import BlenderPlot
from geometry import ParticleGeometry
import bpy


rp = ParticleGeometry.from_lattice_name("fcc")

lx = 8
ly = 12
lz = 1
blender_plot = BlenderPlot.from_lattice_name("fcc", lx, ly, lz)

print(
    f"Possible particle paths for fcc lattice: {blender_plot.particle_paths.keys()}"
)

# Note that this line is optional, as we load this by default

for i, rotation in enumerate(rp.orientation_rotations):
    print(rotation.as_euler("xyz", degrees=True))
    blender_plot.place_obj_copy_from_site_orientation(i * 2, i)

blender_plot.save("./test_fcc_orientations.blend")
