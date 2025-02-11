""" Vincent Ouazan-Rebo, 2025
Storing the contact map designs for triangular lattice particles (i.e. hexagons).
"""

from contact_utils import ContactMapWrapper

CRYSTAL_CONTACTS = [(i, i+3) for i in range(3)]

def set_crystal_contacts(crystal_e:float, cmap: ContactMapWrapper):
    for contact in CRYSTAL_CONTACTS:
        cmap[*contact] = crystal_e
    return

def set_vortex_camembert_contacts(defect_e:float, cmap:ContactMapWrapper):
    camembert_contacts = [(0,4), (1,5)]
    for contact in camembert_contacts:
        cmap[*contact] = defect_e
    return

def gen_vortex_camembert_contacts(crystal_e:float, defect_e:float, mismatch_e:float):
    cmap = ContactMapWrapper.triangular(1, mismatch_e)
    set_crystal_contacts(crystal_e, cmap)
    set_vortex_camembert_contacts(defect_e, cmap)

    return cmap
