# pyright basic
"""
Vincent Ouazan-Reboul, 2025
Tools to implement geometry of hexagonal particles, in order to easily generate contact maps and plot
simulation results

TODOS:
- Find a better way to explain the docstrings
"""

# mathutils is provided by bpy
import numpy as np
from scipy.spatial.transform import Rotation as R
from geometry.particle_geometry import ParticleGeometry
from geometry.lattice_geometry import LatticeGeometry

# Globals
ID_ROT = R.identity()
C6Z = R.from_euler("z", 60, degrees=True)
C2Z = R.from_euler("z", 180, degrees=True)
ORIENTATION_0_VEC = np.array([[1, 0, 0]])
BONDS = [
    (1, 0, 0),
    (0, 1, 0),
    (0, -1, 0),
    (-1, 0, 0),
    (-1, -1, 0),
    (1, -1, 0),
]
BOND_ORIENTATIONS_POSITIVE = [ID_ROT, C6Z, C6Z**2]
BOND_ROTATIONS = [
    *BOND_ORIENTATIONS_POSITIVE,
    *[C2Z * rot for rot in BOND_ORIENTATIONS_POSITIVE],
]
SR32: float = np.sqrt(3) / 2
BASIS_VECTORS = np.array([[1, 0.5, 0.], [0, SR32, 0.]])

class TriangularParticle(ParticleGeometry):
    def __init__(self):
        super().__init__(
            ORIENTATION_0_VEC,
            BONDS,
            BOND_ROTATIONS,
            BOND_ORIENTATIONS_POSITIVE,
            rotations_around_face_0=[ID_ROT],
            opposite_face_rotation=C2Z,
        )


class TriangularLattice(LatticeGeometry):
    def __init__(self, lx: int = 1, ly: int = 1, lattice_spacing: float = 1.0):
        super().__init__(BASIS_VECTORS, BONDS, lx, ly, 1, lattice_spacing)


class TriangularGeometry_bak:
    def __init__(self):
        # We keep track of only x and y internal vectors, z is redundant
        self.orientation_0_vector = np.array([[1, 0]])
        self.n_orientations = 6
        # self.orientation_rotations = self.gen_orientation_rotations()
        self.bond_to_bond_index = {
            ( 1,  0): 0,
            ( 0,  1): 1,
            ( 0, -1): 2,
            (-1,  0): 3,
            (-1, -1): 4,
            (1,  -1): 5,
        }
        # self.rotation_to_bond = self.gen_bond_rotations()

    def get_equivalent_face_pairs(self, orientation1, orientation2):
        return [[orientation1, orientation2]]

    # def get_reverse_orientations_from_orientations(self, orientation_1, orientation_2):
    #     if orientation_1 != self.get_opposite_orientation(orientation_2):
    #         return []
    #     else:
    #         return [[
    #             self.get_opposite_orientation(orientation_2),
    #             self.get_opposite_orientation(orientation_1),
    #         ]]
    #
    #
    # def gen_orientation_rotations(self):
    #     # Rotation operations are defined as quaternions and compositions of quaternions
    #     orientations = [ID_ROT for i in range(self.n_orientations)]
    #     # Generate first rotation
    #     # Generate the faces/orientations along + directions
    #     for i, orientation in enumerate(BOND_ORIENTATIONS_POSITIVE):
    #         for j in range(4):
    #             orientation_index = i* 4 + j
    #             orientations[orientation_index] = (
    #                 C4XM_POWERS[j] * BOND_ORIENTATIONS_POSITIVE[i]
    #             )
    #             orientations[get_opposite_orientation(orientation_index)] = (
    #                 C2Z * orientations[orientation_index]
    #             )
    #     return R.concatenate(orientations)
    #
    # def face_face_to_ref_bond(self, face1, face2, bond_index):
    #     """
    #     Takes a contact between 2 faces through a given bond index and returns the corresponding
    #     orientations of the particles if they were to touch through the reference bond (1, 0,
    #     0). 
    #     """
    #     # Put the particles back in the original orientation
    #     orientation_1_in_ref = self.identify_orientation(
    #         ALL_BOND_ORIENTATIONS[bond_index].inv() * self.orientation_rotations[face1]
    #     )
    #     opposite_orientation_2_in_ref = self.identify_orientation(
    #         ALL_BOND_ORIENTATIONS[bond_index].inv() * self.orientation_rotations[face2]
    #     )
    #     orientation_2_in_ref = get_opposite_orientation(opposite_orientation_2_in_ref)
    #
    #     return orientation_1_in_ref, orientation_2_in_ref
    #
    # def gen_bond_rotations(self):
    #     """
    #     Generate a list associating which bond (1, 0, 0) maps to under every rotation.
    #     """
    #     all_rotated_bonds = []
    #     zero_bond = np.array([1, 0, 0])
    #     for rot in self.orientation_rotations:
    #         # Change type to round coefficients to 0, 1, or -1
    #         new_bond_arr = rot.apply(zero_bond).astype(int)
    #         new_bond = tuple(new_bond_arr)
    #         all_rotated_bonds.append(self.bond_to_bond_index[new_bond])
    #     return all_rotated_bonds
    #
    # def identify_orientation(self, rotation):
    #     """
    #     Identifies which orientation the rotation `rotation` puts the particle in.
    #     Returns -1 if the supplied rotation does not put the particle in one of its 24 possible
    #     orientations.
    #     """
    #     for i, rot in enumerate(self.orientation_rotations):
    #         if rot.approx_equal(rotation):
    #             return i
    #     return -1
    #
    # def apply_rotation(self, rotation_idx, orientation):
    #     """
    #     Applies rotation rotation_idx to particle orientation orientation, and returns the new
    #     orientation of the particle.
    #     """
    #     # Mind the order when composing rotations!
    #     new_rot = (
    #         self.orientation_rotations[rotation_idx]
    #         * self.orientation_rotations[orientation]
    #     )
    #     return self.identify_orientation(new_rot)
    #
    # def get_all_orientations_reproducing_ref_contact(
    #     self, ref_orientation_1, ref_orientation_2
    # ):
    #     """
    #     From the orientations of particle 1 at lattice site (0,0,0) and of particle 2 at site
    #     (1,0,0), fill list `all_orientation_bond` with all the [orientation_1, orientation_2]
    #     combinations reproducing the contact between the 2 particles still occupying the same
    #     spots.
    #     """
    #     all_orientation_sets = []
    #     for rot in C4XM_POWERS:
    #         orientation_1 = self.identify_orientation(
    #             rot * self.orientation_rotations[ref_orientation_1]
    #         )
    #         orientation_2 = self.identify_orientation(
    #             rot * self.orientation_rotations[ref_orientation_2]
    #         )
    #         all_orientation_sets.append([orientation_1, orientation_2])
    #     return all_orientation_sets
    #
    # def get_bond_index(self, bond_vector):
    #     bond = tuple(bond_vector)
    #     if bond not in self.bond_to_bond_index.keys():
    #         print("Invalid bond vector supplied!")
    #         return None
    #     bond_index = self.bond_to_bond_index[bond]
    #     return bond_index
    #
    # def get_all_orientation_bond_contacts(self, face1, face2):
    #     """
    #     For a face-face contact between `face1` and `face2`,
    #     generates all the (orientation, orientation) pairs reproducing this contact through bond
    #     0 as a 2-column array.
    #     """
    #     # A note on how this works: if face1 and face2 are in contact, then they also are in
    #     # contact through bond 0. So we can consider the case where orientation1 = face1 and
    #     # orientation2 = opposite(face2).
    #     # We don't need to consider bonds different from 0: the C++ code already does that!
    #     all_orientation_bond = []
    #     ref_orientation_1 = face1
    #     ref_orientation_2 = get_opposite_orientation(face2)
    #     # ref_orientation_1, ref_orientation_2 = self.face_face_to_ref_bond(
    #     #     face1, face2, bond
    #     # )
    #     # If face1 n
    #     # Apply all possible rotations to the contact we are considering
    #     all_orientation_bond.extend(
    #         self.get_all_orientations_reproducing_ref_contact(
    #             ref_orientation_1, ref_orientation_2
    #         )
    #     )
    #
    #     # We also need to take into account the contacts when 2 is to the left of 1
    #     # in the reference bond!
    #     if face1 != face2:
    #         ref_orientation_1_reverse = get_opposite_orientation(face1)
    #         ref_orientation_2_reverse = face2
    #         all_orientation_bond.extend(
    #             self.get_all_orientations_reproducing_ref_contact(
    #                 ref_orientation_2_reverse, ref_orientation_1_reverse
    #             )
    #         )
    #
    #     return np.vstack(all_orientation_bond)
    #
    def get_opposite_orientation(self, orientation):
        return (orientation + 3) % 6
