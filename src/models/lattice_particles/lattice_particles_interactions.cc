#include "lattice_particles_interactions.h"
#include "lattice_particles_geometry.h"
#include "vector_utils.h"
#include <array>
#include <iomanip>

namespace lattice_particles_space {

void initialize_interactions(interactions_struct &interactions,
                             const state_struct &state,
                             const model_parameters_struct &parameters) {
  interactions.couplings = parameters.couplings;
  interactions.energy = get_energy(state, interactions);
  interactions.n_neighbours =
      static_cast<int>(std::size(interactions.neighbours));
}

double get_contact_energy(state_struct &state, int site1, int site2,
                          interactions_struct &interactions) {
  int state1{state.lattice_sites.get_state(site1)};
  int state2{state.lattice_sites.get_state(site2)};
  int edge1{get_bond_direction(site1, site2, state)};
  int edge2{get_bond_direction(site2, site1, state)};

  return get_contact_energy(interactions.couplings, state1, state2, edge1,
                            edge2, state.n_states);
}

double get_contact_energy(ContactMap contact_map, int state1, int state2,
                          int edge1, int edge2, int n_states) {
  int contact_index{
      array_space::hash_into_contact(state1, state2, edge1, edge2, n_states)};
  return contact_map[static_cast<std::size_t>(contact_index)];
}

double get_site_energy(state_struct &state, interactions_struct &interactions,
                       int site_index) {
  double site_energy{0.0};
  get_neighbours(interactions.neighbours, site_index, state);
  // NOTE that the loop index is also the contact edge for the site @
  // site_index
  // TODO Check that's actually the case...
  for (std::size_t contact_edge{0};
       contact_edge < interactions.neighbours.size(); contact_edge++) {
    int neighbour_ind{interactions.neighbours[contact_edge]};
    // TODO Check I'm looking at the right directions
    site_energy +=
        get_contact_energy(state, site_index, neighbour_ind, interactions);
  }
  return site_energy;
}

double get_energy(state_struct &state, interactions_struct interactions) {
  double energy{0.0};
  for (int i{0}; i < state.n_sites; i++)
    energy += get_site_energy(i, state, interactions);
  // Divide by 2, otherwise we'll be double counting!
  return energy / 2;
}
} // namespace lattice_particles_space
