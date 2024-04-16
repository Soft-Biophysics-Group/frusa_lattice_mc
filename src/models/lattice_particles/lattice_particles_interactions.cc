#include "lattice_particles_interactions.h"
#include "lattice_particles_geometry.h"
#include "lattice_particles_state.h"
#include "vector_utils.h"
#include <array>
#include <iomanip>
#include <typeinfo>

namespace lattice_particles_space {

void initialize_interactions(state_struct &state,
                             interactions_struct &interactions,
                             model_parameters_struct &parameters) {
  interactions.couplings = parameters.couplings;
  interactions.n_edges =
      static_cast<int>(std::size(interactions.neighbours));
  interactions.energy = get_energy(state, interactions);
  interactions.neighbours = Neighbours();
}

double get_contact_energy(state_struct &state, int site1, int site2,
                          interactions_struct &interactions) {
  // Note the +1: get_bond_direction returns values starting at 0, and our
  // edges run from 1 to n_edges, as 0 is empty
  int edge1{get_bond_direction(site1, site2, state)+1};
  //std::cout << "Edge between neighbours is " << edge1 << '\n' ;

  // An empty-empty contact counts as zero energy
  if (state.lattice_sites.is_empty(site1) and
      state.lattice_sites.is_empty(site2)) {
    return 0.0;
  }

  int contact_index{
      get_contact_index(state, site1, site2, edge1, interactions)};

  //std::cout << "Index of corresponding contact is " << contact_index << '\n' ;

  return interactions.couplings[static_cast<std::size_t>(contact_index)];
}

int get_contact_index_full_empty(int edge, int type, int orientation,
                                 int n_edges) {
  //std::cout << "Full-empty contact\n"
            //<< "(type " << type << ", orientation " << orientation << ", edge "
            //<< edge << ")\n";
  return type * n_edges + (std::abs(edge - orientation) % n_edges);
}

int get_contact_index_full_full(int state1, int state2, int edge1, int n_states,
                                int n_edges, int n_types) {
  // When 2 sides touch, they do so along 2 edges: one in the set
  // {1,2,...,n_edges/2}, the other in the set {n_edges/2+1,...,n_edges}
  // We then can get away with storing only half of the coefficients by
  // only storing the ones in the first category.
  // Note we apply an offset of n_edges*n_orientations because the first
  // componentes of the interaction vector are reserved for surface tension
  // terms
  //std::cout << "Full-full contact\n"<<
    //"(state1: " << state1 << ", state2 " << state2 << ", edge1: " << edge1 << ")\n";
  if (edge1 <= n_edges / 2) {
    return n_types * n_edges - 1 +
           array_space::hash_3integers(state1, state2, edge1, n_states,
                                       n_edges/2);
  } else {
    return n_types * n_edges - 1 +
           array_space::hash_3integers(
               state2, state1, get_conjugate_edge(edge1), n_states, n_edges/2);
  }
}

int get_contact_index(state_struct &state, int site1, int site2, int edge1,
                      interactions_struct &interactions) {
  int contact_index{0};
  // The first n_edges*n_types coefficients of the LEL are for surface tension
  // terms.
  // Coefficients 0 to n_states-1 are for the surface tensions of edges of
  // particles of type 1, n_states to 2(n_states - 1) for the edges of type 2,
  // etc...
  if (state.lattice_sites.is_empty(site1)) {
    contact_index = get_contact_index_full_empty(
        get_conjugate_edge(edge1), state.lattice_sites.get_type(site2),
        state.lattice_sites.get_orientation(site2), interactions.n_edges);}
  else if (state.lattice_sites.is_empty(site2)) {
    contact_index = get_contact_index_full_empty(
        edge1, state.lattice_sites.get_type(site1),
        state.lattice_sites.get_orientation(site1), interactions.n_edges);}
  else {
    int state1{state.lattice_sites.get_state(site1)};
    int state2{state.lattice_sites.get_state(site2)};
    contact_index =
        get_contact_index_full_full(state1, state2, edge1, state.n_states,
                                    interactions.n_edges, state.n_types);
  }
  return contact_index;
}

double get_site_energy(state_struct &state, interactions_struct &interactions,
                       int site_index) {
  //std::cout << "\nGetting energy of site " << site_index << '\n' ;
  double site_energy{0.0};
  get_neighbours(interactions.neighbours, site_index, state);
  // NOTE that the loop index is also the contact edge for the site @
  // site_index
  // TODO Check that's actually the case...
  for (int neighbour_site : interactions.neighbours) {
    //std::cout << "Checking neighbour " << neighbour_site << '\n' ;
    // TODO Check I'm looking at the right directions
    site_energy +=
        get_contact_energy(state, site_index, neighbour_site, interactions);
  }
  return site_energy;
}

double get_energy(state_struct &state, interactions_struct& interactions) {
  double energy{0.0};
  for (int i{0}; i < state.n_sites; i++)
    energy += get_site_energy(state, interactions, i);
  // Divide by 2, otherwise we'll be double counting!
  //std::cout << "Energy is now " << energy / 2 ;
  return energy / 2;
}

std::ostream& operator<< (std::ostream& out, interactions_struct& interactions) {
  out << "Printing interactions structure\n" ;
  out << "Printing coupling matrix: ";
  array_space::print_vector(out, interactions.couplings);
  out << '\n' ;
  out << "System energy is: " << interactions.energy << '\n';
  out << "Last verified neighbours were: ";
  array_space::print_array<int, 6>(out, interactions.neighbours);
  out << "Stored number of site edges:" << n_edges << '\n';

  return out;
}

void print_interactions(state_struct &state, interactions_struct &interactions) {
  std::cout << interactions;
}


void print_energy(state_struct &state, interactions_struct &interactions) {
  std::cout << interactions.energy ;
}
} // namespace lattice_particles_space
