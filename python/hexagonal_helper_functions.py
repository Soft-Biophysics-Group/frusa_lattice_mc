import numpy as np
from lattice_helper_classes import *


def uniform_couplings_hexagonal(type_couplings) -> list:
    n_types = type_couplings.shape[0]
    params = system_parameters(n_types, 6, 6)
    couplings = np.zeros(params.n_total_terms, float)
    # Set surface tensions to 0
    # for i in range(params.n_surface_tension_terms):
    # couplings[i] = 0.0
    for i in range(params.n_surface_tension_terms, params.n_total_terms):
        state1, state2, edge = states_edge_from_contact_index(i, params)
        print(f"{i}: {state1}, {state2}, {edge}")
        type1 = site_parameters.from_state(params, state1).get_ptype()
        type2 = site_parameters.from_state(params, state2).get_ptype()
        print(f"\t{type1}, {type2}")
        couplings[i] = type_couplings[type1, type2]
    return list(couplings)


def camembert_couplings_hexagonal(
    crystal_int: float, line_int: float, surf_tension: float, repulsion: float
) -> list:
    params = system_parameters(1, 6, 6)
    couplings = params.init_couplings()
    # Set surface tensions
    for i in range(params.n_surface_tension_terms):
        couplings[i] = surf_tension
    # Listing all the contacts corresponding to a grain boundary
    defect_configs = [(1,6,1), (2,3,1), (3,4,1), (6,5,1)]
    for i in range(params.n_surface_tension_terms, params.n_total_terms):
        # Here orientation and state are the same because we have only 1 type
        # of particle
        orientation1, orientation2, edge = states_edge_from_contact_index(
            i, params)
        if orientation1 == orientation2:
            couplings[i] = crystal_int
        elif (np.abs(orientation1-orientation2) % 6 == 5 or
              np.abs(orientation1-orientation2) % 6 == 1):
            couplings[i] = line_int
        else:
            couplings[i] = repulsion
    return list(couplings)
