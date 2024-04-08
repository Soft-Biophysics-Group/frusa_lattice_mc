#include "lattice_particles_interactions.h"
#include "lattice_particles_geometry.h"
#include "vector_utils.h"
#include <iomanip>

namespace lattice_particles_space {
void initialize_interactions(state_struct &state,
                             interactions_struct &interactions,
                             model_parameters_struct &parameters) {
  interactions.couplings = parameters.couplings;
  interactions.energy = get_energy(state, interactions);
}

double get_contact_energy(site_state &center_site, site_state &neighbour_site,
                          int center_site_edge, int neighbour_site_edge,
                          ContactMap contact_map, int n_states) {
  int contact_index{array_space::hash_into_contact(
      center_site.get_state(), neighbour_site.get_state(),
      center_site_edge, neighbour_site_edge, n_states)};
  return contact_map[contact_index];
}

template <int N>
double get_site_energy(int site_index, state_struct &state,
                       ContactMap contact_map, Neighbours<N> neighbours) {
  site_state& center_site {state.lattice_sites[site_index]};
  double site_energy{0.0};
  get_neighbours(neighbours, site_index, state.lx, state.ly, state.lz);
  // NOTE that the loop index is also the contact edge for the site @
  // site_index
  // TODO Check that's actually the case...
  for (std::size_t contact_edge {0}; contact_edge < N; contact_edge++) {
    int neighbour_ind{neighbours[contact_edge]};
    site_state& neighbour_site{state.lattice_sites[neighbour_ind]};
    int neighbour_site_edge{get_bond_direction(site_index, neighbour_ind,
                                                state.lx, state.ly, state.lz)};
    site_energy += get_contact_energy(
        center_site, neighbour_site, static_cast<int>(contact_edge),
        neighbour_site_edge, contact_map, state.n_states);
      }
  return site_energy;
}

template <int N>
double get_energy(state_struct& state, interactions_struct<N> interactions) {
  double energy {0.0};
  for (int site_index{0}; site_index < state.n_sites; site_index++)
    energy += get_site_energy(site_index, state, interactions.couplings,
                              interactions.neighbours);
  // Divide by 2, otherwise we'll be double counting!
  return energy / 2;
}
} // namespace lattice_particles_space
