import numpy as np

# ---------- GENERAL FUNCTIONS ----------


def flatten_couplings(coupling_arr) -> list:
    return coupling_arr.flatten().tolist()

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
