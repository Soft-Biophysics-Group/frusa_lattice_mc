from sys import exception
import numpy as np

# ---------- GENERAL FUNCTIONS ----------


def flatten_couplings(coupling_arr) -> list:
    sym_coupling_arr = (coupling_arr + coupling_arr.T)/2
    return sym_coupling_arr.flatten().tolist()

# SINGLE-TYPE SYSTEMS

def block_pauli_x(n: int):
    """ Returns the block Pauli x matrix of size n.
    Necessary for moving from face-based to orientation-based matrices """
    if n % 2 != 0:
        print("Dimension should be even!")
        return np.zeros(n)
    halfdim = n//2
    zero_block = np.zeros(halfdim)
    id_block = np.eye(halfdim)
    return np.block([
        [ zero_block, id_block  ],
        [ id_block,   zero_block],
        ])


def face_to_or(face_mat):
    """ Converts a face-based interaction matrix to an orientation-based one,
    which is the format taken as input by the C++ code """
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


def uniform_interaction_tri(e):
    """
        Returns an interaction matrix for triangular lattices where all
        coefficients are the same.
        Parameters:
        - e: interaction energy
    """
    return e * np.ones((6, 6))

def is_symmetric(matrix):
    return np.isclose(matrix, (matrix + matrix.T) / 2).prod()


class ContactMapWrapper:
    """
    Wrapper class to easily manipulate a contact map object.
    The mirror image of what is contained in the C++ code.
    """

    def __init__(self, n_types, n_orientations):
        self.n_types = n_types
        self.n_orientations = n_orientations
        self.n_states = self.n_types * self.n_orientations
        self.contact_map = np.zeros(self.n_states**2)


    # Different lattices for which we can use this class
    @classmethod
    def triangular(cls, n_types):
        return cls(n_types, 6)

    def get_one_face_coeff(self, face, type):
        return type * self.n_orientations + face

    def get_interaction_coeff(self, face1, type1, face2, type2):
        coeff_1 = self.get_one_face_coeff(face1, type1)
        coeff_2 = self.get_one_face_coeff(face2, type2)
        return coeff_1 + coeff_2 * self.n_types * self.n_orientations

    def __getitem__(self, faces_types):
        (face1, type1, face2, type2) = faces_types
        int_coeff = self.get_interaction_coeff(face1, type1, face2, type2)

        return self.contact_map[int_coeff]

    def __setitem__(self, faces_types, value):
        (face1, type1, face2, type2) = faces_types
        int_coeff = self.get_interaction_coeff(face1, type1, face2, type2)

        self.contact_map[int_coeff] = value

    def get_two_species_contact_matrix(self, type1, type2):
        contact_map = np.zeros((self.n_orientations, self.n_orientations))
        for face1 in range(self.n_orientations):
            for face2 in range(self.n_orientations):
                contact_map[face1, face2] = self[face1, type1, face2, type2]
        return contact_map

    def get_single_species_contact_matrix(self, type):
        return self.get_two_species_contact_matrix(type, type)

    def set_two_species_contacts(self, type1, type2, contact_matrix, symmetrize=True):
        """Contact_matrix has to be a n_orientations by n_orientations np array"""
        if contact_matrix.shape != (self.n_orientations, self.n_orientations):
            print("Invalid matrix format! Stopping now")
            return
        # Ensure contact matrix is symmetric
        if not is_symmetric(contact_matrix) and not symmetrize:
            print("Contact energies are not symmetric! Aborting")
            return
        else:
            contact_matrix = (contact_matrix + contact_matrix.T) / 2

        for face1 in range(self.n_orientations):
            for face2 in range(self.n_orientations):
                self[face1, type1, face2, type2] = contact_matrix[face1, face2]
        return

    def set_single_species_contact(self, type, contact_matrix):
        self.set_two_species_contacts(type, type, contact_matrix)
        return

    def get_formatted_couplings(self):
        """Returns the contact  matrix in a format"""
        # Check the contacts are symmetric
        for type1 in range(self.n_types):
            for type2 in range(self.n_types):
                if not is_symmetric(self.get_two_species_contact_matrix(type1, type2)):
                    print(
                        f"Contact matrix between species {type1} and {type2} is not symmetric."
                        "Aborting"
                    )
                    return
        return list(self.contact_map)


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

    cmap_wrapper.set_single_species_contact(0, contact_map_matrix)
    print(cmap_wrapper.get_single_species_contact_matrix(0))

    return cmap_wrapper.get_formatted_couplings()

