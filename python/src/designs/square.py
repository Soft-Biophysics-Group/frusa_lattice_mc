""" Vincent Ouazan-Reboul, 2026
Storing the contact map designs for square lattice particles.
"""

from contact_utils import ContactMapWrapper

CRYSTAL_CONTACTS = [(i, i+2) for i in range(2)]
SQUARE_VORTEX_CONTACTS = [(0, 1), (0, 3)]


def set_crystal_contacts(crystal_e:float, cmap: ContactMapWrapper):
    for contact in CRYSTAL_CONTACTS:
        cmap[*contact] = crystal_e
    return

def set_vortex_contacts(defect_e:float, cmap:ContactMapWrapper):
    for contact in SQUARE_VORTEX_CONTACTS:
        cmap[*contact] = defect_e
    return

def gen_vortex_contacts(
    defect_e: float, crystal_e: float = -0.5, mismatch_e: float = 10.0
) -> ContactMapWrapper:
    cmap = ContactMapWrapper.from_lattice_name(
        "square", n_types=1, init_energy=mismatch_e
    )
    set_crystal_contacts(crystal_e, cmap)
    set_vortex_contacts(defect_e, cmap)

    return cmap
