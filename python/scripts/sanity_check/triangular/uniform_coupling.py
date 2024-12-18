from helper_fct import sanity_check_tri
from contact_utils import uniform_interaction_tri

coupling_energy = - 1.0
coupling_matrix = uniform_interaction_tri(coupling_energy)

sanity_check_tri(coupling_matrix)
