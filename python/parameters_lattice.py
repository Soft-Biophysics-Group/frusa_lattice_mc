import numpy as np
import os
from json_dump import *

# Define model parameters

model_params = {}

model_params["n_types"] = 2
model_params["n_orientations"] = 6

# Lattice dimensions
model_params["lx"] = 5
model_params["ly"] = 5
model_params["lz"] = 1

# Number of particles
model_params["n_particles"] = [10, 10]


def signed_mod(a, b):
    return (a % b + b) % b


class system_parameters:
    def __init__(self, n_types, n_orientations, n_edges):
        self.n_types = n_types
        self.n_orientations = n_orientations
        self.n_edges = n_edges
        self.n_states = self.n_orientations * self.n_types
        self.n_surface_tension_terms = self.n_edges * self.n_types
        self.n_coupling_terms = self.n_states ** 2 * self.n_edges // 2
        self.n_total_terms = self.n_surface_tension_terms + self.n_coupling_terms

    def get_conjugate_edge(edge: int):
        return signed_mod(edge + self.n_edges // 2, self.n_edges)


def get_contact_index_full_empty(edge, ptype, orientation, n_edges):
    return ptype * n_edges + signed_mod((edge-orientation), n_edges)


def hash_3integers(x1, x2, y, x_range, y_range):
    return (x1-1)*x_range*y_range + (x2-1)*y_range + y


def unhash_3integers(hashed_int, x_range, y_range):
    x1 = (hashed_int // x_range // y_range) + 1
    x2 = (hashed_int // y_range) % x_range + 1
    y = hashed_int % y_range + 1
    return x1, x2, y

def states_edge_from_contact_index(contact_int, params):
    hashed_int = contact_int + 1 - params.n_types * params.n_edges
    return unhash_3integers(hashed_int-1, params.n_states, params.n_edges//2)


def get_contact_index_full_full(state1, state2, edge1, params):
    if edge1 <= params.n_edges//2:
        return n_types * n_edges - 1 + hash_3integers(state1, state2, edge1,
                                                      params.n_state, params.n_edges//2)
    else:
        edge_2 = params.get_conjugate_edge(edge1)
        return n_types * n_edges - 1 + hash_3integers(state2, state1, edge2,
                                                      params.n_state, params.n_edges//2)


class site_parameters:
    def __init__(self, params: system_parameters, orientation: int, ptype: int):
        self.params = params
        self.orientation = orientation
        self.ptype = ptype
        self.is_empty = (orientation == 0)
        self.state = self.calc_state()

    @classmethod
    def from_state(cls, params: system_parameters, state: int):
        params = site_parameters(params, 0, 0)
        params.state = state
        params.orientations = params.get_orientation()
        params.ptype = params.get_ptype()
        return params

    def calc_state(self):
        if self.is_empty:
            return 0
        else:
            return self.orientation + self.ptype * self.params.n_orientations

    def get_orientation(self):
        return self.state % self.params.n_orientations

    def get_ptype(self):
        return (self.state-1) // self.params.n_orientations

    def get_contact_index(site2, edge: int):
        if self.isempty():
            edge2 = self.params.get_conjugate_edge(edge)
            return get_contact_index_full_empty(edge2, site2.ptype,
                                                site2.orientation,
                                                self.params.n_edges)
        elif site2.isempty():
            return get_contact_index_full_empty(edge, self.ptype,
                                                self.orientaiton,
                                                self.params.n_edges)
        else:
            return get_contact_index_full_full(self.state, site2.state, edge,
                                               params)


# Model parameters
def uniform_couplings_hexagonal(type_couplings) -> list:
    n_types = type_couplings.shape[0]
    params = system_parameters(n_types, 6, 6)
    couplings = np.zeros(params.n_total_terms, float)
    # Set surface tensions to 0
    # for i in range(params.n_surface_tension_terms):
        # couplings[i] = 0.0
    for i in range(params.n_surface_tension_terms, params.n_total_terms):
        state1, state2, edge = states_edge_from_contact_index(i,
                                                              params)
        print(f'{i}: {state1}, {state2}, {edge}')
        type1 = site_parameters.from_state(params, state1).get_ptype()
        type2 = site_parameters.from_state(params, state2).get_ptype()
        print(f'\t{type1}, {type2}')
        couplings[i] = type_couplings[type1, type2]
    return list(couplings)


n_types = 2
type_couplings = np.array([[-1.0, -0.5],
                           [-0.5,  0.5]])

model_params["couplings"] = uniform_couplings_hexagonal(type_couplings)

model_params["initialize_option"] = "random_fixed_particle_numbers"

# # swap_empty_full_enum,
  # swap_full_full_enum,
  # rotate_enum,
  # mutate_enum,
model_params["move_probas"] = [
    1/4, # swap_empty_full
    1/4,   # swap_full_full
    1/4, # rotate
    1/4    # mutate
    ]


make_json_file(model_params, "../input/model_params.json")

# Define mc parameters

mc_params = {}

# Number of MC steps used for equilibration
mc_params["mcs_eq"] = 100

# Number of MC steps used for averaging
mc_params["mcs_av"] = 1

# Type of cooling schedule
mc_params["cooling_schedule"] = "exponential"

# Initial annealing temperature
mc_params["Ti"] = 0

# Final annealing temperature
mc_params["Tf"] = -5

# Number of annealing steps
mc_params["Nt"] = 10

# Option to collect state checkpoints at the end of each temperature cycle
mc_params["checkpoint_option"] = False

# If checkpoint is True, we need to provide the output address for the
# checkpoint files
if mc_params["checkpoint_option"]:
    mc_params["checkpoint_address"] = ""

make_json_file(mc_params, "../input/mc_params.json")
