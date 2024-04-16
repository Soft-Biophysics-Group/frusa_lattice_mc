import numpy as np


def signed_mod(a, b):
    return (a % b + b) % b


class system_parameters:
    def __init__(self, n_types, n_orientations, n_edges):
        self.n_types = n_types
        self.n_orientations = n_orientations
        self.n_edges = n_edges
        self.n_states = self.n_orientations * self.n_types
        self.n_surface_tension_terms = self.n_edges * self.n_types
        self.n_coupling_terms = self.n_states**2 * self.n_edges // 2
        self.n_total_terms = self.n_surface_tension_terms + self.n_coupling_terms

    def get_conjugate_edge(self, edge: int):
        return signed_mod(edge + self.n_edges // 2, self.n_edges)

    def init_couplings(self):
        return np.zeros(self.n_total_terms, float)


def get_contact_index_full_empty(edge, ptype, orientation, n_edges):
    return ptype * n_edges + signed_mod((edge - orientation), n_edges)


def hash_3integers(x1, x2, y, x_range, y_range):
    return (x1 - 1) * x_range * y_range + (x2 - 1) * y_range + y


def unhash_3integers(hashed_int, x_range, y_range):
    x1 = (hashed_int // x_range // y_range) + 1
    x2 = (hashed_int // y_range) % x_range + 1
    y = hashed_int % y_range + 1
    return x1, x2, y


def states_edge_from_contact_index(contact_int, params):
    hashed_int = contact_int + 1 - params.n_types * params.n_edges
    return unhash_3integers(hashed_int - 1, params.n_states, params.n_edges // 2)


def get_contact_index_full_full(state1, state2, edge1, params):
    if edge1 <= params.n_edges // 2:
        return (
            params.n_types * params.n_edges
            - 1
            + hash_3integers(state1, state2, edge1, params.n_state, params.n_edges // 2)
        )
    else:
        edge2 = params.get_conjugate_edge(edge1)
        return (
            params.n_types * params.n_edges
            - 1
            + hash_3integers(state2, state1, edge2, params.n_state, params.n_edges // 2)
        )


class site_parameters:
    def __init__(self, params: system_parameters, orientation: int, ptype: int):
        self.params = params
        self.orientation = orientation
        self.ptype = ptype
        self.is_empty = orientation == 0
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
        return (self.state - 1) // self.params.n_orientations

    def rotate_by(self, n_rotations):
        self.orientation = (self.orientation +
                            n_rotations) % self.n_orientations

    def get_contact_index(self, site2, edge: int):
        if self.isempty():
            edge2 = self.params.get_conjugate_edge(edge)
            return get_contact_index_full_empty(
                edge2, site2.ptype, site2.orientation, self.params.n_edges
            )
        elif site2.isempty():
            return get_contact_index_full_empty(
                edge, self.ptype, self.orientaiton, self.params.n_edges
            )
        else:
            return get_contact_index_full_full(
                self.state, site2.state, edge, self.params
            )

    def put_to_reference_orientation(self, orientation_ref=1):
        n_rotations = 0
        while self.orientation != orientation_ref:
            self.rotate_by(1)
            n_rotations += 1
        return n_rotations

    def put_other_in_ref(self, site2):
        n_rotations = site2.put_to_reference_orientation()
        self.rotate_by(n_rotations)
