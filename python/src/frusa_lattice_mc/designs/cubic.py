""" Vincent Ouazan-Reboul, 2025
Storing the contact map designs for cubic lattice particles.
"""

from ..contact_utils import ContactMapWrapper
from ..geometry.particle_geometry import CubicParticle

CRYSTAL_CONTACTS = CubicParticle().get_all_canonical_contacts(
    [(4 * i, (4 * i + 12) % 24) for i in range(3)]
)
CRYSTAL_CONTACTS_ONE_AXIS = CubicParticle().get_all_canonical_contacts(
    [
        (0, 12),
        (0, 13),
        (0, 14),
        (0, 15),
        (4, 6),
        (4, 11),
        (4, 16),
        (4, 21),
        (8, 8),
        (8, 19),
        (8, 20),
        (16, 18),
        (16, 23),
        (20, 20),
    ]
)
ALL_TRUE_CAMEMBERT = {
    "Rminus_flag": True,
    "Gminus_flag": True,
    "Bminus_flag": True,
    "Rplus_flag": True,
    "Gplus_flag": True,
    "Bplus_flag": True,
}
HEDGEHOG_CONTACTS = CubicParticle().get_all_canonical_contacts(
    [
        (12, 4),
        (12, 5),
        (12, 6),
        (12, 7),
        (12, 8),
        (12, 9),
        (12, 10),
        (12, 11),
        (12, 16),
        (12, 17),
        (12, 18),
        (12, 19),
        (12, 20),
        (12, 21),
        (12, 22),
        (12, 23),
    ]
)
# HEDGEHOG_CONTACTS = [
#     (4, 12),
#     (4, 13),
#     (4, 14),
#     (4, 15),
#     (8, 12),
#     (8, 13),
#     (8, 14),
#     (8, 15),
#     (12, 20),
#     (12, 21),
#     (12, 22),
#     (12, 23),
#     (16, 12),
#     (16, 13),
#     (16, 14),
#     (16, 15),
#     # These 4 should be removable w/o creating issues. I'm too much of a coward to try though!
#     (20, 12),
#     (20, 13),
#     (20, 14),
#     (20, 15),
# ]
EXTRA_HEDGEHOG_DEFECT_CONTACTS = CubicParticle().get_all_canonical_contacts(
    [(12, 12), (12, 13), (12, 14), (12, 15)]
)


def get_triple_vortex_camembert_contacts(
    contact_flags: dict[str, bool] = ALL_TRUE_CAMEMBERT,
) -> frozenset[tuple[int, int]]:
    contacts_Rplus = CubicParticle().get_all_canonical_contacts(
        [(4, 13), (8, 13), (12, 17), (12, 23)]
    )
    contacts_Gplus = CubicParticle().get_all_canonical_contacts(
        [(0, 17), (8, 16), (12, 17), (16, 20)]
    )
    contacts_Bplus = CubicParticle().get_all_canonical_contacts(
        [(0, 21), (4, 22), (12, 23), (16, 20)]
    )
    contacts_Rminus = CubicParticle().get_all_canonical_contacts(
        [(0, 9), (0, 17), (0, 21), (0, 7)]
    )
    contacts_Gminus = CubicParticle().get_all_canonical_contacts(
        [(4, 22), (4, 13), (4, 10), (0, 7)]
    )
    contacts_Bminus = CubicParticle().get_all_canonical_contacts(
        [(0, 9), (4, 10), (8, 13), (8, 16)]
    )

    camembert_contacts = set()
    if contact_flags["Rplus_flag"]:
        camembert_contacts |= contacts_Rplus
    if contact_flags["Gplus_flag"]:
        camembert_contacts |= contacts_Gplus
    if contact_flags["Bplus_flag"]:
        camembert_contacts |= contacts_Bplus
    if contact_flags["Rminus_flag"]:
        camembert_contacts |= contacts_Rminus
    if contact_flags["Gminus_flag"]:
        camembert_contacts |= contacts_Gminus
    if contact_flags["Bminus_flag"]:
        camembert_contacts |= contacts_Bminus

    return frozenset(camembert_contacts)
