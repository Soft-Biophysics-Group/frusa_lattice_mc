import config as cfg
from geometry import ParticleGeometry, LatticeGeometry
import numpy as np
from scipy.spatial.transform import Rotation as R
from plotting import BlenderPlot
from plotting.fcc import path_to_numbered_rhombic

p = ParticleGeometry.from_lattice_name("fcc")
print(len(p.orientation_rotations))

print("Testing generated rotation against printed dodecahedron")
all_rots_euler = np.array(
    [rot.as_euler("xyz", degrees=True) for rot in p.orientation_rotations]
)
print(f"Rotations in Euler angle, xyz order: {all_rots_euler}")

# print("Testing generated rotated bonds")
# print(p.rotation_to_bond)

# Should be rotation associated to orientation 15
test_rot = R.from_euler('xyz', (90, 180, 0), degrees = True)
print("\nTesting a newly-generated rotation against one of the cubes'")
test_rot_index = p.identify_orientation(test_rot)
print(f"Rotation identified as index {test_rot_index}")

# Apply rotation 4 to orientation 2.
# Result should be orientation 18
test_rot = 4
test_orientation = 2
new_ori = p.apply_rotation(test_rot, test_orientation)
print(
    f"\nApplying rotation associated with orientation {test_rot} "
    f"to orientation {test_orientation} yields orientation {new_ori}"
)

# Let's print all rotations to be safe, and plot them if we want to

blender_plot = BlenderPlot.from_lattice_name("fcc", 8, 12, 1, path_to_numbered_rhombic)

for i, rotation in enumerate(p.orientation_rotations):
    print(rotation.as_euler("xyz", degrees = True))
    blender_plot.place_obj_copy_from_site_orientation(i * 2, i)

blender_plot.save("./test_rhombic_orientations.blend")
