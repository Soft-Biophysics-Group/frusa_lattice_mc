"""Vincent Ouazan-Rebo, 2025
Storing the contact map designs for fcc lattice particles.
"""

from contact_utils import ContactMapWrapper

CRYSTAL_CONTACTS: list[tuple[int, int]] = [(2 * i, (2 * i + 12) % 24) for i in range(6)]

RED_MINUS_FACES = [10, 11, 12, 13, 16, 17, 20, 21]
ORTHOGONAL_RED_FACES = [2, 3, 6, 7, 14, 15, 18, 19]
HEDGEHOG_CONTACTS = [
    (face_1, face_2)
    for face_1 in RED_MINUS_FACES
    for face_2 in RED_MINUS_FACES
    if face_1 != face_2
] + [(face_1, face_2) for face_1 in RED_MINUS_FACES for face_2 in ORTHOGONAL_RED_FACES]

# BACK_FACES = list(range(12, 24))
# HEDGEHOG_CONTACTS_ALT = [
#     (face_1, face_2)
#     for face_1 in BACK_FACES
#     for face_2 in BACK_FACES
#     if face_1 != face_2
# ]
# HEDGEHOG_CONTACTS_ALT_STRICTER = [
#     (face_1, face_2) for face_1 in (12, 13) for face_2 in BACK_FACES
# ]
#
FACE_01_CRYSTAL_CONTACTS = [
    (4, 16),
    (6, 17),
    (4, 19),
    (6, 18),
    (8, 20),
    (8, 15),
    (2, 14),
    (0, 12),
    (0, 13),
    (10, 22),
    (22, 23),
    (2, 21),
    (10, 11),
]

SEED_CONTACTS = [
    (14, 19),
    (14, 16),
    (18, 20),
    (16, 21),
]

# BACK_ORTH_ALT = [
#     (12, 22),
#     (13, 22),
#     (12, 10),
#     (13, 10),
# ]
#
# OPP_ORTH_ALT = [
#     (10, 21),
#     (10, 14),
#     (14, 23),
#     (20, 22),
#     (10, 17),
#     (16, 22),
#     (18, 23),
#     (10, 18),
# ]

HOLE_3_CONTACTS = [
    (4, 23),
    (4, 10),
    (6, 11),
    (6, 22),
    (2, 18),
    (2, 17),
    (8, 16),
    (8, 19),
    (16, 16),
    (16, 19),
    (18, 18),
    (4, 13),
    (6, 12),
    (6, 13),
    (4, 12),
    (2, 13),
    (2, 12),
    (8, 12),
    (8, 13),
    (14, 16),
    (16, 21),
    (18, 20),
    (14, 19),
]

FRONT_FACES = (2, 4, 7 ,9)
OPP_FACES = (12, 13)

# Promising
OPP_FRONT_ALT = [
    (face_1, face_2)
    for face_1 in FRONT_FACES
    for face_2 in OPP_FACES
]

# OPP_FRONT_ALT = [
#     (4, 12),
#     (4, 13),
#     (6, 12),
#     (6, 13),
#     (2, 12),
#     (2, 13),
#     (8, 12),
#     (8, 13),
# ]
#

# Far-fetched
BACK_ORTH_ALT = [
    (18, 23),
    (10, 18),
    (10, 17),
    (16, 22),
    (20, 22),
    (14, 23),
    (10, 14),
    (10, 21),
]

# Promising
BUTT_BACK_ALT = [
    (12, 18),
    (13, 18),
    (12, 21),
    (13, 21),
    (12, 16),
    (13, 16),
    (12, 14),
    (13, 14),
]


# Promising, but maybe problematic?
OPP_ORTH_ALT = [
    (12, 23),
    (13, 23),
    (10, 12),
    (10, 13),
]

# Maybe problematic
BUTT_FRONT_ALT = [
    (4, 12),
    (4, 13),
    (2, 12),
    (2, 13),
    (6, 12),
    (6, 13),
    (8, 12),
    (8, 13),
]

# Maybe problematic
BACK_FRONT_ALT = [
    (4, 21),
    (6, 20),
    (6, 15),
    (4, 14),
    (2, 19),
    (2, 16),
    (8, 17),
    (8, 18),
]

# Promising, but maybe problematic
BUTT_BUTT_ALT = [
    (12, 12),
    (12, 13),
]

SEED_CORNER = [
    (10, 11),
    (10, 21),
    (10, 12),
    (10, 17),
    (12, 13),
    (12, 16),
    (12, 20),
    (16, 21),
    (16, 17),
    (20, 21),
]

CRYSTAL_INTS_CORNER = [
    (0, 21),
    (0, 12),
    (0, 17),
    (0, 11),
    (4, 20),
    (4, 13),
    (4, 16),
    (4, 10),
    (20, 22),
    (12, 23),
    (16, 22),
    (10, 22),
    (8, 20),
    (8, 13),
    (8, 16),
    (8, 10),
    (6, 18),
    (18, 18),
    (14, 19),
    (14, 14),
    (6, 15),
    (2, 14),
    (2, 7),
    (2, 2),
    (2, 19),
    (6, 6),
]


SEED_HOLE_INTS_CORNER = [
    (2, 10),
    (2, 11),
    (2, 12),
    (2, 13),
    (2, 16),
    (2, 17),
    (2, 20),
    (2, 21),
    (6, 10),
    (6, 11),
    (6, 12),
    (6, 13),
    (6, 16),
    (6, 17),
    (6, 20),
    (6, 21),
    (10, 14),
    (10, 19),
    (10, 14),
    (10, 15),
    (11, 19),
    (12, 14),
    (12 ,15),
    (12, 18),
    (12, 19),
    (14, 16),
    (14, 17),
    (14, 20),
    (14, 21),
    (16, 19),
    (17, 19),
    (18, 20),
    (18, 21),
]

BUTT_TO_ORTH_CORNER = [
    (14, 17),
    (12, 14),
    (14, 21),
    (10, 15),
    (19, 20),
    (12, 18),
    (16, 19),
    (10, 19),
    (2, 21),
    (2, 11),
    (2, 17),
    (2 , 12),
    (6, 17),
    (6, 11),
    (6, 21),
    (6, 12),
]

DEFECT_CORNER = SEED_CORNER + SEED_HOLE_INTS_CORNER


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
    cmap = ContactMapWrapper.from_lattice_name("fcc", 1, mismatch_e)
    set_crystal_contacts(crystal_e, cmap)
    set_hedgehog_camembert_contacts(defect_e, cmap)

    return cmap

def gen_hedgehog_camembert_cmap_alt(
    crystal_e: float, defect_e: float, mismatch_e: float
):
    cmap = ContactMapWrapper.from_lattice_name("fcc", 1, mismatch_e)
    cmap.set_contacts(CRYSTAL_CONTACTS, crystal_e)
    cmap.set_contacts(HEDGEHOG_CONTACTS_ALT, defect_e)

    return cmap


def gen_hedgehog_camembert_cmap_alt_stricter(
    crystal_e: float, defect_e: float, mismatch_e: float
):
    cmap = ContactMapWrapper.from_lattice_name("fcc", 1, mismatch_e)
    cmap.set_contacts(CRYSTAL_CONTACTS, crystal_e)
    cmap.set_contacts(HEDGEHOG_CONTACTS_ALT_STRICTER, defect_e)

    return cmap
