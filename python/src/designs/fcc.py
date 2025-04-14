"""Vincent Ouazan-Rebo, 2025
Storing the contact map designs for fcc lattice particles.
"""

from contact_utils import ContactMapWrapper

CRYSTAL_CONTACTS = [(2 * i, (2 * i + 12) % 24) for i in range(6)]

RED_MINUS_FACES = [10, 11, 12, 13, 16, 17, 20, 21]
ORTHOGONAL_RED_FACES = [2, 3, 6, 7, 14, 15, 18, 19]
HEDGEHOG_CONTACTS = [
    (face_1, face_2)
    for face_1 in RED_MINUS_FACES
    for face_2 in RED_MINUS_FACES
    if face_1 != face_2
] + [(face_1, face_2) for face_1 in RED_MINUS_FACES for face_2 in ORTHOGONAL_RED_FACES]


def set_contacts(
    cmap: ContactMapWrapper, contact_list: list[tuple[int, int]], value: float
):
    for contact in contact_list:
        cmap[*contact] = value
    return


def set_crystal_contacts(crystal_e: float, cmap: ContactMapWrapper):
    for contact in CRYSTAL_CONTACTS:
        cmap[*contact] = crystal_e
    return


def set_hedgehog_camembert_contacts(defect_e: float, cmap: ContactMapWrapper):
    for contact in HEDGEHOG_CONTACTS:
        cmap[*contact] = defect_e

    return


def gen_hedgehog_camembert_cmap(crystal_e: float, defect_e: float, mismatch_e: float):
    cmap = ContactMapWrapper.cubic(1, mismatch_e)
    set_crystal_contacts(crystal_e, cmap)
    set_hedgehog_camembert_contacts(defect_e, cmap)

    return cmap
