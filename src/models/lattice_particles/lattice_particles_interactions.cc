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

double get_site_energy(site_state &site, state_struct &state,
                       interactions_struct& interactions) {
  double site_energy{0.0};
  get_neighbours(interactions.neighbours, site.get_site_index(),
      state.lx, state.ly, state.lz);
  // NOTE that the loop index is also the contact edge for the site @
  // site_index
  // TODO Check that's actually the case...
  for (std::size_t contact_edge {0}; contact_edge < N; contact_edge++) {
    int neighbour_ind{interactions.neighbours[contact_edge]};
    site_state& neighbour_site{state.lattice_sites[neighbour_ind]};
    // TODO Check I'm looking at the right directions
    int neighbour_site_edge{get_bond_direction(neighbour_ind, site.get_site_index(),
                                                state.lx, state.ly, state.lz)};
    site_energy += get_contact_energy(
        site, neighbour_site, static_cast<int>(contact_edge),
        neighbour_site_edge, interactions.couplings, state.n_states);
      }
  return site_energy;
}

double get_energy(state_struct& state, interactions_struct interactions) {
  double energy {0.0};
  for (site_state& site : state.lattice_sites)
    energy += get_site_energy(site, state, interactions);
  // Divide by 2, otherwise we'll be double counting!
  return energy / 2;
}
} // namespace lattice_particles_space
