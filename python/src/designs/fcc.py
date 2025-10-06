"""Vincent Ouazan-Reboul, 2025
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

DEFECT_CORNER = SEED_CORNER + SEED_HOLE_INTS_CORNER
