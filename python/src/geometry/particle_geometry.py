"""Generic class describing lattice particles.

Uses: easy generation of contact maps and plotting of simulation results.
"""

# mathutils is provided by bpy
import numpy as np
from scipy.spatial.transform import Rotation as R
from numpy.typing import NDArray


class ParticleGeometry:
    """Generic class describing an abstract lattice particle.

    Attributes:
        orientation_0_vectors: array of int. Orientation of the particle's first 2 internal axes in
            orientation 0, in lattice coordinates.
        n_orientations: int, number of particle orientations.
        face_0_permutation_rotations: list of rotations which permute the face / feature 0
            (occupying position 0 in orientation 0) with other faces. Should only cover the
            "positive-oriented" faces, i.e. the ones along positive bond directions (+x, +y,
            +z) for cubes for instance
        rotations_around_face_0: list of all possible rotations which do not change permute
            faces. Ex: quarter turn rotations around one face of a cube.
        orientation_rotations: list of the rotations defining all possible orientations of the
            particle.
        bond_rotations: list of rotations which permute face 0 with the different possible bonds
            linking a particle to its nearest neighbours
    """
    def __init__(
        self,
        orientation_0_vectors: NDArray[np.int_],
        face_0_permutation_rotations,
        rotations_around_face_0,
        opposite_face_rotation,
        bond_rotations
    ):
        """Creates abstract lattice particle object.

        Args:
            orientation_0_vectors: array of int. Orientation of the particle's first 2 internal
                axes in orientation 0, in lattice coordinates.
            face_0_permutation_rotations: list of rotations which permute the face / feature 0
                (occupying position 0 in orientation 0) with other faces. Should only cover the
                "positive-oriented" faces, i.e. the ones along positive bond directions (+x, +y,
                +z) for cubes for instance.
            rotations_around_face_0: list of all possible rotations which do not change permute
                faces. Ex: quarter turn rotations around one face of a cube.
            opposite_face_rotation: rotation which maps one face to its opposite on the
                particle.
            bond_rotations: list of rotations which permute face 0 with the different possible
                bonds linking a particle to its nearest neighbours
        """
        # We keep track of only x and y internal vectors, z is redundant
        self.orientation_0_vectors = orientation_0_vectors
        self.n_orientations = (
            len(face_0_permutation_rotations) * len(rotations_around_face_0) * 2
        )
        self.face_0_permutation_rotations = face_0_permutation_rotations
        self.rotations_around_face_0 = rotations_around_face_0
        self.orientation_rotations = self.gen_face_orientations(
            face_0_permutation_rotations,
            rotations_around_face_0,
            opposite_face_rotation,
        )
        self.bond_rotations = bond_rotations
        self.bond_permutations = self.gen_bond_permutations()

    def gen_face_orientations(
        self,
        face_0_permutation_rotations,
        rotations_around_face_0,
        opposite_face_rotation,
    ):
        """Generates all the possible orientations of the particle.
        """
        # Rotation operations are defined as quaternions and compositions of quaternions
        rotations = [R.identity() for i in range(self.n_orientations)]
        # Generate first rotation
        # Generate the faces/orientations along + directions
        for i, bond_rotation in enumerate(face_0_permutation_rotations):
            for j in range(len(rotations_around_face_0)):
                face_index = i * len(rotations_around_face_0) + j
                rotations[face_index] = (
                    rotations_around_face_0[j] * face_0_permutation_rotations[i]
                )
                rotations[self.get_opposite_face_index(face_index)] = (
                    opposite_face_rotation * rotations[face_index]
                )
        return R.concatenate(rotations)

    def gen_bond_permutations(self):
        """ Returns an array describing how the particle orientations are permuted when a given bond
        is considered.
        """
        all_bond_permutations = []
        for bond_rotation in self.bond_rotations:
            these_permutations = []
            for orientation_rotation in self.orientation_rotations:
                permuted_orientation_rotation = bond_rotation * orientation_rotation
                these_permutations.append(
                    self.identify_orientation(permuted_orientation_rotation)
                )
            all_bond_permutations.append(these_permutations)

        return all_bond_permutations

    def identify_orientation(self, rotation):
        """Identifies which orientation the rotation `rotation` puts the particle in, defined as
        the face which takes the place of face 0.
        Returns -1 if the supplied rotation does not put the particle in one of its 24 possible
        orientations.
        """
        for i, rot in enumerate(self.orientation_rotations):
            if rot.approx_equal(rotation):
                return i
        return -1

    def apply_rotation(self, rotation_idx, orientation):
        """Applies rotation rotation_idx to particle orientation orientation, and returns the new
        orientation of the particle.
        """
        # Mind the order when composing rotations!
        new_rot = (
            self.orientation_rotations[rotation_idx]
            * self.orientation_rotations[orientation]
        )
        return self.identify_orientation(new_rot)

    def get_equivalent_face_pairs(self, face_1, face_2):
        """
        From the faces in contact of particle 1 at lattice site (0,0,0) and of particle 2 at site
        (1,0,0), fill list `all_face_pairs` with all the [face_1, face_2]
        combinations reproducing the contact between the 2 particles still occupying the same
        spots.
        """
        all_face_pairs = []
        opposite_face_2 = self.get_opposite_face_index(face_2)
        for rot in self.rotations_around_face_0:
            orientation_1 = self.identify_orientation(
                rot * self.orientation_rotations[face_1]
            )
            orientation_2 = self.identify_orientation(
                rot * self.orientation_rotations[opposite_face_2]
            )
            all_face_pairs.append(
                [orientation_1, self.get_opposite_face_index(orientation_2)]
            )
        return all_face_pairs

    def get_faces_in_contact(self, orientation1, orientation2, bond):
        if orientation1 == -1:
            face_1 = -1
        else:
            face_1 = self.bond_permutations[bond][orientation1]
        if orientation2 == -1:
            face_2 = -1
        else:
            face_2 = self.get_opposite_face_index(
                self.bond_permutations[bond][orientation2]
            )
        return face_1, face_2

    def get_equivalent_orientations_from_orientations(
        self, orientation_1, orientation_2
    ):
        all_equivalent_orientation_pairs = []
        # ref_orientation_1, ref_orientation_2 = self.face_face_to_ref_bond(
        #     face1, face2, bond
        # )
        # If face1 n
        # Apply all possible rotations to the contact we are considering
        all_equivalent_orientation_pairs.extend(
            self.get_equivalent_face_pairs(orientation_1, orientation_2)
        )

        return np.vstack(all_equivalent_orientation_pairs)

    def get_reverse_orientations_from_orientations(self, orientation_1, orientation_2):
        all_opposite_orientation_pairs = []
        if orientation_1 != self.get_opposite_face_index(orientation_2):
            orientation_1_reverse = self.get_opposite_face_index(orientation_1)
            orientation_2_reverse = self.get_opposite_face_index(orientation_2)
            return self.get_equivalent_orientations_from_orientations(
                orientation_2_reverse, orientation_1_reverse
            )
        else:
            return []

    def get_opposite_face_index(self, orientation):
        return (orientation + (self.n_orientations // 2)) % self.n_orientations
