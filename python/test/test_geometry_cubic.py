import config as cfg
from geometry.cubic import CubicGeometry
import numpy as np
from scipy.spatial.transform import Rotation as R
from plotting.plot_cubic import place_cube_copy_from_site_orientation, save_blender_fig

cp = CubicGeometry.particle

print("Testing generated rotation against printed cube")
all_rots_euler = np.array(
    [rot.as_euler("xyz", degrees=True) for rot in cp.orientation_rotations]
)
print(f"Rotations in Euler angle, xyz order: {all_rots_euler}")

print("Testing generated rotated bonds")
print(cp.rotation_to_bond)

# Should be rotation associated to orientation 15
test_rot = R.from_euler('xyz', (90, 180, 0), degrees = True)
print("\nTesting a newly-generated rotation against one of the cubes'")
test_rot_index = cp.identify_orientation(test_rot)
print(f"Rotation identified as index {test_rot_index}")

# Apply rotation 4 to orientation 2.
# Result should be orientation 18
test_rot = 4
test_orientation = 2
new_ori = cp.apply_rotation(test_rot, test_orientation)
print(
    f"\nApplying rotation associated with orientation {test_rot} "
    f"to orientation {test_orientation} yields orientation {new_ori}"
)

# Let's print all rotations to be safe, and plot them if we want to
plot_flag = True
for i, rotation in enumerate(cp.orientation_rotations):
    print(rotation.as_euler("xyz"))
    place_cube_copy_from_site_orientation(i*2, i, 8, 12)

save_blender_fig("./test_cubic_orientation.blend")
