#include "lattice_particles_interactions.h"
#include <iomanip>

namespace lattice_particles_space {
void initialize_interactions(state_struct &state,
                             interactions_struct &interactions,
                             model_parameters_struct &parameters) {
  interactions.couplings = parameters.couplings;
  interactions.energy = get_energy(state, interactions);
}

template <int N>
double get_site_energy(int site_index, state_struct &state,
                       ContactMap coupling_matrix, Neighbours<N> neighbours) {
    // TODO Stopped here
    double site_energy{0.0};
    get_neighbours(neighbours, site_index, state.lx, state.ly, state.lz);
    for (site_state& neigh_site : neighbours) {
      
    }
}
  // Print the summary of the interactions characteristics
  void print_interactions(state_struct &state,
                          interactions_struct &interactions);

  void print_energy(state_struct &state,
                    interactions_struct &interactions);

  double get_energy(state_struct& state, ContactMap interactions);

  double get_site_energy(int site_index, state_struct &state,
                         ContactMap coupling_matrix);

  // Function to update energy after concentration transfer from donor to
  // the acceptor field. r_d             - lattice position of the donor r_a
  // - lattice position of the acceptor index_d         - donor field
  // component index_a         - acceptor field component dc              -
  // amount of concentration transferred
  double get_energy_change(int r_d, int r_a, int index_d, int index_a,
                           double dc, state_struct &state,
                           interactions_struct &interactions);

} // namespace lattice_particles_space
