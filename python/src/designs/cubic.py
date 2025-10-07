""" Vincent Ouazan-Rebo, 2025
Storing the contact map designs for cubic lattice particles.
"""

from contact_utils import ContactMapWrapper

CRYSTAL_CONTACTS = [(4 * i, (4 * i + 12) % 24) for i in range(3)]
CRYSTAL_CONTACTS_ONE_AXIS = [
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
ALL_TRUE_CAMEMBERT = {
    "Rminus_flag": True,
    "Gminus_flag": True,
    "Bminus_flag": True,
    "Rplus_flag": True,
    "Gplus_flag": True,
    "Bplus_flag": True,
}
HEDGEHOG_CONTACTS = [
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
EXTRA_HEDGEHOG_DEFECT_CONTACTS = [(12, 12), (12, 13), (12, 14), (12, 15)]


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

def set_crystal_one_axis_contacts(crystal_e: float, cmap: ContactMapWrapper):
    set_contacts(cmap, CRYSTAL_CONTACTS_ONE_AXIS, crystal_e)


def set_triple_vortex_camembert_contacts(
    defect_e: float,
    cmap: ContactMapWrapper,
    contact_flags: dict[str, bool] = ALL_TRUE_CAMEMBERT,
):
    contacts_Rplus = [(4, 13), (8, 13), (12, 17), (12, 23)]
    contacts_Gplus = [(0, 17), (8, 16), (12, 17), (16, 20)]
    contacts_Bplus = [(0, 21), (4, 22), (12, 23), (16, 20)]
    contacts_Rminus = [(0, 9), (0, 17), (0, 21), (0, 7)]
    contacts_Gminus = [(4, 22), (4, 13), (4, 10), (0, 7)]
    contacts_Bminus = [(0, 9), (4, 10), (8, 13), (8, 16)]

    camembert_contacts = []
    if contact_flags["Rplus_flag"]:
        camembert_contacts.extend(contacts_Rplus)
    if contact_flags["Gplus_flag"]:
        camembert_contacts.extend(contacts_Gplus)
    if contact_flags["Bplus_flag"]:
        camembert_contacts.extend(contacts_Bplus)
    if contact_flags["Rminus_flag"]:
        camembert_contacts.extend(contacts_Rminus)
    if contact_flags["Gminus_flag"]:
        camembert_contacts.extend(contacts_Gminus)
    if contact_flags["Bminus_flag"]:
        camembert_contacts.extend(contacts_Bminus)

    for contact in camembert_contacts:
        cmap[*contact] = defect_e
    return

def gen_triple_vortex_camembert_contacts(crystal_e:float, defect_e:float, mismatch_e:float):
    cmap = ContactMapWrapper.cubic(1, mismatch_e)
    set_crystal_contacts(crystal_e, cmap)
    set_triple_vortex_camembert_contacts(defect_e, cmap)

    return cmap


def set_hedgehog_camembert_contacts(defect_e: float, cmap: ContactMapWrapper):
    for contact in HEDGEHOG_CONTACTS:
        cmap[*contact] = defect_e

    return 

def set_hedgehog_camembert_contacts_w_extra(defect_e: float, cmap: ContactMapWrapper):
    for contact in HEDGEHOG_CONTACTS:
        cmap[*contact] = defect_e
    for contact in EXTRA_HEDGEHOG_DEFECT_CONTACTS:
        cmap[*contact] = defect_e

    return 

def gen_hedgehog_camembert_cmap(crystal_e:float, defect_e:float, mismatch_e:float):
    cmap = ContactMapWrapper.cubic(1, mismatch_e)
    set_crystal_contacts(crystal_e, cmap)
    set_hedgehog_camembert_contacts(defect_e, cmap)

    return cmap


def gen_hedgehog_camembert_cmap_crystal_one_axis(
    crystal_e: float, defect_e: float, mismatch_e: float
):
    cmap = ContactMapWrapper.cubic(1, mismatch_e)
    set_crystal_one_axis_contacts(crystal_e, cmap)
    set_hedgehog_camembert_contacts(defect_e, cmap)

    return cmap


def gen_hedgehog_camembert_cmap_extra_contacts(
    crystal_e: float, defect_e: float, mismatch_e: float
):
    cmap = ContactMapWrapper.cubic(1, mismatch_e)
    set_crystal_one_axis_contacts(crystal_e, cmap)
    set_hedgehog_camembert_contacts_w_extra(defect_e, cmap)

    return cmap
