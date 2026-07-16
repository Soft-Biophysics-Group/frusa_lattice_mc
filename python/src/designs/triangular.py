""" Vincent Ouazan-Rebo, 2025
Storing the contact map designs for triangular lattice particles (i.e. hexagons).
"""

from contact_utils import ContactMapWrapper

CRYSTAL_CONTACTS = [(i, i+3) for i in range(3)]
VORTEX_CAMEMBERT_CONTACTS = [(0,4), (1,5)]

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

def get_vortex_camembert_colormap(
    crystal_color="cyan", defect_color="navy", mismatch_color="red", surface_color="orange"
):
    colors = {}
    # Contacts with enpty site
    for i in range(6):
        colors[i, -1] = surface_color
    # All the other contacts we get through a bit of a twisted method:
    # we create a contact map, which will have all the coefficients figured out for us, and then
    # assign the colors based on their value.
    cmap = gen_vortex_camembert_contacts(0, 1, 2)
    for i in range(6):
        for j in range(6):
            if cmap[i,j] == 0:
                colors[i,j] = crystal_color
            elif cmap[i, j] == 1:
                colors[i, j] = defect_color
            else:
                colors[i,j] = mismatch_color

    return colors
