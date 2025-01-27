# Andrey Zelenskiy, Vincent Ouazan-Reboul, 2024
# Functions and classes to generate contact maps for input of `frusa_mc`.
# pyright: basic
from sys import exception
import numpy as np
from numpy.typing import NDArray
from geometry import geometry
from geometry.cubic import CubicGeometry
from geometry.triangular import TriangularGeometry

LATTICE_NAMES = ["chain", "triangular", "cubic"]

class ContactMapWrapper:
    """
    Wrapper class to easily manipulate a contact map object.
    The mirror image of what is contained in the C++ code: changes made to one need to be made
    to the other.
    The class contains a `contact_map` member which is a flattened interaction map used
    directly as input for the C++ code.
    Its goal is to make filling that array as convenient as possible.

    Note that the interaction matrices manipulated by that class assume that all contacts happen
    through bond 0.
    The typical workflow when designing an interaction matrix will be:
    1. Create a class instance using the constructor associated to the lattice you will be
    working with
    2. Get the empty single-species and species-pair interaction matrices using the
    `get_single_species_contact_matrix` and `get_two_species_contact_matrix` methods.
    3. Setting the coefficients of the interaction matrices by hand.
    4. Fill up the flattened `contact_map` member array using the `set_single_species_contacts`
    and `set_two_species_contacts` method
    5. Get the properly formatted couplings using the `get_formatted_couplings` method.

    If you only need to set a couple of faces, you can do so using the syntax:
    `cmap_wrapper[face_1, type_1, face_2, type_2] = energy`.
    The __setitem__ function takes care of setting all the coefficients redundant with the
    initial geometry as well: i.e. the ones corresponding to particle 2 on the left and 1 on the
    right, the ones which should be equal through rotational invariance, etc...

    ## Constructors:
    - `triangular(n_types)`: creates a class instance for a triangular lattice containing
      `n_types` different particle types.
    - `cubic(n_types)`: creates a class instance for a cubic lattice containing
      `n_types` different particle types.
    - `__init(n_types, n_orientations)__`: creates a class instance for `n_types`
      different particle types and `n_orientations` particle orientations.
      Usually not called directly, but used by class constructors.
    """

    # ----- CONSTRUCTORS -----
    def __init__(
        self,
        n_types=1,
        n_orientations=6,
        lattice_geometry=TriangularGeometry(),
        init_energy=0.0,
    ):
        self.n_types = n_types
        self.n_orientations = n_orientations
        self.n_states = self.n_types * self.n_orientations
        self.contact_map = np.zeros(self.n_states**2)
        self.contact_map += init_energy
        self.geometry = lattice_geometry

    # Different lattices for which we can use this class
    @classmethod
    def triangular(cls, n_types, init_energy = 0.0):
        return cls(n_types, 6, TriangularGeometry(), init_energy = init_energy)

    @classmethod
    def cubic(cls, n_types, init_energy = 0.0):
        return cls(n_types, 24, CubicGeometry(), init_energy = init_energy)

    # ----- GETTER FUNCTION FOR COEFFICIENTS IN FLATTENED ARRAY -----
    def get_one_face_coeff(self, orientation, type):
        """
        Returns a single hashed index describing the face of a particle, taking into account its
        type: faces 0 to n_orientation-1 belong to species 0, 1 to 2*n_orientations-1 to 1,
        etc...
        """
        return type * self.n_orientations + orientation

    def get_interaction_coeff(self, face1, type1, face2, type2):
        coeff_1 = self.get_one_face_coeff(face1, type1)
        coeff_2 = self.get_one_face_coeff(face2, type2)
        return coeff_1 + coeff_2 * self.n_types * self.n_orientations


    # ----- MORE CONVENIENT ACCESS TO COEFFICIENTS -----
    def __getitem__(self, faces_types):
        if len(faces_types) == 2:
            face1, face2 = faces_types
            type1, type2 = 0, 0
        else:
            (face1, type1, face2, type2) = faces_types
        int_coeff = self.get_interaction_coeff(face1, type1, face2, type2)

        return self.contact_map[int_coeff]

    def __setitem__(self, faces_types, value):
        if len(faces_types) == 2:
            face1, face2 = faces_types
            type1, type2 = 0, 0
        else:
            (face1, type1, face2, type2) = faces_types
        # If we are in 3D, some pairs of faces might be equivalent to the one we are looking at
        # through rotational invariance. If so, we must generate these.
        equiv_faces = self.get_equivalent_face_pairs(face1, face2)
        # All equivalent contacts with particle 1 at the center and through bond 0
        for this_face1, this_face2 in equiv_faces:
            int_coeff = self.get_interaction_coeff(this_face1, type1, this_face2, type2)
            self.contact_map[int_coeff] = value
            # We also need to fill the coefficients with particle 1 and particle 2 swapped to
            # respect symmetry of the interactions
            int_coeff_reversed = self.get_interaction_coeff(
                this_face2, type2, this_face1, type1
            )
            self.contact_map[int_coeff_reversed] = value
        # for (
        #     this_face_2,
        #     this_face_1,
        # ) in self.geometry.get_reverse_orientations_from_orientations(face1, face2):
        #     int_coeff = self.get_interaction_coeff(
        #         this_face2, type2, this_face1, type1
        #     )
        #     self.contact_map[int_coeff] = value

    # -----SETTING EQUIVALENT COEFFICIENTS -----
    def get_equivalent_face_pairs(self, face1, face2):
        """
        Depending on the lattice geometry, when setting a contact with an (orientation,
        orientation) pair, there will be some amount of other (orientation, orientation) pairs
        reproducing the same contact through bond 0.
        This function generates all of them for the implemented lattices (so far only cubic).
        TODO: Implement triangular lattice
        """

        return self.geometry.get_equivalent_face_pairs(
            face1, face2
        )

    # ----- GETTING 2D CONTACT MATRICES -----
    def get_two_species_contact_matrix(self, type1, type2):
        contact_map = np.zeros((self.n_orientations, self.n_orientations))
        for face1 in range(self.n_orientations):
            for face2 in range(self.n_orientations):
                contact_map[face1, face2] = self[
                    face1, type1, face2, type2
                ]
        return contact_map

    def get_single_species_contact_matrix(self, type):
        return self.get_two_species_contact_matrix(type, type)

    # -----  GOING FROM 2D MATRICES TO FLATTENED ARRAYS -----
    def set_two_species_contacts(self, type1, type2, contact_matrix):
        """Contact_matrix has to be a n_orientations by n_orientations np array"""
        if contact_matrix.shape != (self.n_orientations, self.n_orientations):
            print("Invalid matrix format! Stopping now")
            return

        for face1 in range(self.n_orientations):
            for face2 in range(self.n_orientations):
                self[face1, type1, face2, type2] = contact_matrix[
                    face1, face2
                ]
                # self[
                #     self.geometry.get_opposite_orientation(face2),
                #     type2,
                #     self.geometry.get_opposite_orientation(face1),
                #     type1,
                # ] = contact_matrix[face1, face2]
        return

    def set_single_species_contacts(self, type:int, contact_matrix):
        """
        Set couplings of particle species with index type, from numpy array contact_matrix.
        contact_matrix has to be symmetric, or upper/lower triangular, otherwise the
        interactions between particles are non-reciprocal, which is cool but out of the scope of
        this program.
        """
        if (contact_matrix == contact_matrix.T).all():
            print("Symmetric contact matrix")
            self.set_two_species_contacts(type, type, contact_matrix)
        elif (np.tril(contact_matrix) == contact_matrix).all() or (
            np.triu(contact_matrix) == contact_matrix
        ).all():
            print("Triangular contact matrix")
            self.set_two_species_contacts(type, type, contact_matrix + contact_matrix.T)
        else:
            print(
                "Contact matrix is neither diagonal, nor upper or lower triangular."
                "Aborting"
            )

        return

    # ----- GETTING THE FORMATTED COUPLINGS FOR JSON INPUT ----- 
    def get_formatted_couplings(self):
        """Returns the contact  matrix in a format"""
        return list(self.contact_map)


# ---------- GENERAL FUNCTIONS ----------
def flatten_couplings(coupling_arr) -> list[float]:
    sym_coupling_arr = (coupling_arr + coupling_arr.T) / 2
    return sym_coupling_arr.flatten().tolist()


# SINGLE-TYPE SYSTEMS
def block_pauli_x(n: int):
    """Returns the block Pauli x matrix of size n.
    Necessary for moving from face-based to orientation-based matrices"""
    if n % 2 != 0:
        print("Dimension should be even!")
        return np.zeros(n)
    halfdim = n // 2
    zero_block = np.zeros(halfdim)
    id_block = np.eye(halfdim)
    return np.block(
        [
            [zero_block, id_block],
            [id_block, zero_block],
        ]
    )


def face_to_or(face_mat):
    """Converts a face-based interaction matrix to an orientation-based one,
    which is the format taken as input by the C++ code"""
    return np.matmul(face_mat, block_pauli_x(face_mat.shape[0]))


def symmetrize(face_mat):
    return (face_mat + face_mat.T) / 2


# SEVERAL-TYPES SYSTEMS
def block_block_pauli_x(n_faces: int, n_types: int):
    return np.tile(block_pauli_x(n_faces), (n_types, n_types))


# ---------- CHAIN FUNCTIONS ----------
def chain_LEL_1type(one_to_one, two_to_two, one_to_two):
    """ 
        Returns orientation-based 1D interaction matrix.
        Parameters:
        - one_to_one: face 1 to face 1 interaction
        - two_to_two: face 2 to face 2 interaction
        - one_to_two: face 1 to face 2 interaction
    """
    face_mat = np.array([
        [one_to_one, one_to_two],
        [one_to_two, two_to_two]])
    return face_to_or(face_mat)


def chain_LEL_2types(mat_11, mat_21, mat_22):
    full_face_matrix = np.block([
        [mat_11,   mat_21],
        [mat_21.T, mat_22]
    ])
    return np.matmul(block_block_pauli_x(2, 2), full_face_matrix)


# ---------- TRIANGULAR LATTICE FUNCTIONS ----------

def get_camembert_cmap(e_crystal, e_defect, e_repel):
    cmap_wrapper = ContactMapWrapper.triangular(1)
    contact_map_matrix = cmap_wrapper.get_single_species_contact_matrix(0)
    contact_map_matrix += e_repel
    # Crystal contacts: opposite faces every time
    for face in range(6):
        contact_map_matrix[face, (face + 3)%6 ] = e_crystal
    # Line contacts
    contact_map_matrix[0, 2] = e_defect
    contact_map_matrix[1, 5] = e_defect
    # We also need the symmetric contacts!
    contact_map_matrix[2, 0] = e_defect
    contact_map_matrix[5, 1] = e_defect
    print(contact_map_matrix)

    cmap_wrapper.set_single_species_contacts(0, contact_map_matrix)
    print(cmap_wrapper.get_single_species_contact_matrix(0))

    return cmap_wrapper.get_formatted_couplings()

