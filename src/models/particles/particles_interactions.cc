#include "particles_interactions.h"
#include <array>
#include <iomanip>
#include <typeinfo>

namespace particles_space {

void initialize_interactions(state_struct &state,
                             interactions_struct &interactions,
                             model_parameters_struct &parameters,
                             geometry_space::Geometry geometry
                             ) {
  interactions.couplings = parameters.couplings;
  interactions.energy = get_energy(state, interactions, geometry);
}

double get_contact_energy(state_struct &state, int site1, int site2,
                          interactions_struct &interactions,
                          geometry_space::Geometry geometry) {
  // std::cout << "Edge between neighbours is " << edge1 << '\n' ;

  // Contacts with empty sites count as 0 energy
  if (state.lattice_sites.is_empty(site1) or
      state.lattice_sites.is_empty(site2)) {
    return 0.0;
  }

  int site_1_orientation {state.lattice_sites.get_orientation(site1)};
  int site_2_orientation {state.lattice_sites.get_orientation(site2)};
  int site_1_type {state.lattice_sites.get_type(site1)};
  int site_2_type {state.lattice_sites.get_type(site2)};

  return geometry.get_interaction(site_1_orientation, site_1_type, site1,
                                  site_2_orientation, site_2_type, site2,
                                  interactions.couplings);
}

double get_site_energy(state_struct &state, interactions_struct &interactions,
                       geometry_space::Geometry geometry, int site_index) {
  // std::cout << "\nGetting energy of site " << site_index << '\n' ;
  double site_energy{0.0};
  if (state.lattice_sites.is_empty(site_index)) {
    return 0.0;
  } else {
    for (int bond{0}; bond < geometry.get_n_neighbours(); ++bond) {
      // std::cout << "Checking neighbour " << neighbour_site << '\n' ;
      int neighbour_site{geometry.get_neighbour(site_index, bond)};
      site_energy += get_contact_energy(state, site_index, neighbour_site,
                                        interactions, geometry);
    }
  }
  return site_energy;
}

double get_energy(state_struct &state, interactions_struct &interactions,
                  geometry_space::Geometry geometry) {
  double energy{0.0};
  for (int i{0}; i < state.n_sites; i++)
    energy += get_site_energy(state, interactions, geometry, i);
  // Divide by 2, otherwise we'll be double counting!
  // std::cout << "Energy is now " << energy / 2 ;
  return energy / 2;
}

std::ostream &operator<<(std::ostream &out, interactions_struct &interactions) {
  out << "Printing interactions structure\n";
  out << "Printing coupling matrix: ";
  array_space::print_vector(out, interactions.couplings);
  out << '\n';
  out << "System energy is: " << interactions.energy << '\n';
  //out << "Last verified neighbours were: ";
  //array_space::print_array<int, 6>(out, interactions.neighbours);
  return out;
}

void print_interactions(interactions_struct &interactions) {
  std::cout << interactions;
}

void print_energy(interactions_struct &interactions) {
  std::cout << interactions.energy;
}
} // namespace lattice_particles_space
