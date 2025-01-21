import config as cfg
from cubic import CubicParticle, get_opposite_orientation
import numpy as np
from scipy.spatial.transform import Rotation as R

cp = CubicParticle()

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

# Getting opposite orientation to 18
# Result should be 6
test_orientation = 18
opp_orientation = get_opposite_orientation(test_orientation)
print(f"Opposite orientation to {test_orientation} is {opp_orientation}")

# Putting one (face1, face2) formatted contact into (orientation1, orientation2, 0) format
# with orientation1 as low as can be given face1 == in the reference orientation associated with
# the face
# Contact between 21 and 4 should map to orientations 20 and 19
face1 = 21
face2 = 4
orientation1, orientation2, _ = cp.face_face_to_ref_bond(face1, face2)
print(
    f"\nContact between {face1} and {face2} maps to reference orientations "
    f"{orientation1} and {orientation2} along bond 0"
)
all_contacts = cp.get_all_orientation_bond_contacts(face1, face2)
print(
    f"All of the corresponding (orientation, orientation, bond) combinations are: {all_contacts}"
)
