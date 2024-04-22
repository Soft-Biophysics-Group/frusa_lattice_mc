#include "particles_update.h"
#include <ranges>
#include <stdexcept>

namespace particles_space {

void update_system(state_struct &state, interactions_struct &interactions,
                   model_parameters_struct &parameters,
                   geometry_space::Geometry geometry, double T) {
  // Pick the kind of move we'll be making
  for (int i{0}; i < state.n_sites; i++) {
    mc_moves chosen_move{pick_random_move(parameters)};
    //std::cout << chosen_move << '\n';
    switch (chosen_move) {
    case mc_moves::swap_empty_full_enum:
      interactions.energy +=
          attempt_swap_empty_full(state, parameters, interactions, geometry, T);
      break;
    case mc_moves::swap_full_full_enum:
      interactions.energy +=
          attempt_swap_full_full(state, parameters, interactions, geometry, T);
      break;
    case mc_moves::rotate_enum:
      interactions.energy +=
          attempt_rotate(state, parameters, interactions, geometry, T);
      break;
    case mc_moves::mutate_enum:
      interactions.energy +=
          attempt_mutate(state, parameters, interactions, geometry, T);
      break;
    default:
      throw std::runtime_error("Something went wrong in the move selection");
    }
  }
}

mc_moves pick_random_move(model_parameters_struct &parameters) {
  real_dist move_dist{0.0, 1.0};
  double sampled_real{move_dist(parameters.rng)};
  // Standard tower sampling algorithm
  std::size_t move_index{0};
  double cumulative_prob{parameters.move_probas[0]};
  while (cumulative_prob < sampled_real) {
    ++move_index;
    cumulative_prob += parameters.move_probas[move_index];
  }

  return static_cast<mc_moves>(move_index);
}

double attempt_swap_sites(int index1, int index2, state_struct &state,
                          model_parameters_struct &parameters,
                          interactions_struct &interactions,
                          geometry_space::Geometry geometry, double T) {
  double energy_change{0.0};
  bool sites_are_neighbours{geometry.are_neighbours(index1, index2)};
  //std::cout << "Sites are neighbours: " << sites_are_neighbours << '\n' ;
  // Initial energy, possibly including neighbour correction
  energy_change -= get_site_energy(state, interactions, geometry, index1) +
                   get_site_energy(state, interactions, geometry, index2);
  // If the sites are neighbours, we need to avoid double counting
  if (sites_are_neighbours) {
    energy_change +=
        get_contact_energy(state, index1, index2, interactions, geometry);
  }
  // Make move and calculate energy after
  swap_sites(state, index1, index2);
  energy_change += get_site_energy(state, interactions, geometry, index1) +
                   get_site_energy(state, interactions, geometry, index2);
  if (sites_are_neighbours) {
    energy_change -=
        get_contact_energy(state, index1, index2, interactions, geometry);
  }
  //std::cout << "Energy change is: " << energy_change << '\n' ;
  // Accept or reject move
  if (is_move_accepted(energy_change, T, parameters)) {
    //std::cout << "Move accepted!\n" ;
    return energy_change;
  } else {
    //std::cout << "Move rejected :(\n";
    swap_sites(state, index1, index2);
    return 0.0;
  }
}

double attempt_swap_empty_full(state_struct &state,
                               model_parameters_struct &parameters,
                               interactions_struct &interactions,
                               geometry_space::Geometry geometry, double T) {
  int full_site_index{state.full_empty_sites.get_random_full_site(parameters)};
  int empty_site_index{
      state.full_empty_sites.get_random_empty_site(parameters)};
  //std::cout << "Attempting swap of empty site " << empty_site_index
            //<< " with full site " << full_site_index << '\n';
  return attempt_swap_sites(full_site_index, empty_site_index, state,
                            parameters, interactions, geometry, T);
}

double attempt_swap_full_full(state_struct &state,
                              model_parameters_struct &parameters,
                              interactions_struct &interactions,
                              geometry_space::Geometry geometry, double T) {
  int site1{state.full_empty_sites.get_random_full_site(parameters)};
  int site2{state.full_empty_sites.get_random_full_site(parameters)};
  // Let's avoid swapping a site with itself: if we picked the same site twice,
  // try again
  while (site2 == site1)
    site2 = state.full_empty_sites.get_random_full_site(parameters);
  //std::cout << "Attempting swap of full sites " << site1 << " and " << site2
            //<< '\n';
  return attempt_swap_sites(site1, site2, state, parameters, interactions,
                            geometry, T);
}

double attempt_rotate(state_struct &state, model_parameters_struct &parameters,
                      interactions_struct &interactions,
                      geometry_space::Geometry geometry, double T) {
  int site_index{state.full_empty_sites.get_random_full_site(parameters)};
  double energy_change{-get_site_energy(state, interactions, geometry, site_index)};

  // Lower bound for orientation is 1, as orientation 0 corresponds
  // to an empty site
  int_dist rot_dist{0, state.n_orientations-1};
  int new_orientation{rot_dist(parameters.rng)};
  int old_orientation{state.lattice_sites.get_orientation(site_index)};
  // Let's avoid doing a rotation to the same orientation as before
  while (new_orientation == old_orientation)
    new_orientation = rot_dist(parameters.rng);
  //std::cout << "Attempting rotation of site " << site_index
            //<< " with orientation " << old_orientation << " to "
            //<< new_orientation << '\n';
  state.lattice_sites.set_orientation(site_index, new_orientation);
  energy_change += get_site_energy(state, interactions, geometry, site_index);
  //std::cout << "Energy change is: " << energy_change << '\n' ;
  if (is_move_accepted(energy_change, T, parameters)) {
    //std::cout << "Move accepted!\n";
    return energy_change;
  } else {
    //std::cout << "Move rejected!\n";
    state.lattice_sites.set_orientation(site_index, old_orientation);
    return 0.0;
  }
}

double attempt_mutate(state_struct &state, model_parameters_struct &parameters,
                      interactions_struct &interactions,
                      geometry_space::Geometry geometry, double T) {
  int site_index{state.full_empty_sites.get_random_full_site(parameters)};
  double energy_change{-get_site_energy(state, interactions, geometry, site_index)};

  // Lower bound for orientation is 1, as orientation 0 corresponds
  // to an empty site
  int_dist type_dist{0, state.n_types-1};
  int new_type{type_dist(parameters.rng)};
  int old_type{state.lattice_sites.get_type(site_index)};
  // Let's avoid doing a mutation to the same type as before
  while (new_type == old_type)
    new_type = type_dist(parameters.rng);
  //std::cout << "Attempting rotation of site " << site_index
            //<< " with orientation " << old_orientation << " to "
            //<< new_orientation << '\n';
  state.lattice_sites.set_type(site_index, new_type);
  energy_change += get_site_energy(state, interactions, geometry, site_index);
  //std::cout << "Energy change is: " << energy_change << '\n' ;
  if (is_move_accepted(energy_change, T, parameters)) {
    //std::cout << "Move accepted!\n";
    return energy_change;
  } else {
    //std::cout << "Move rejected!\n";
    state.lattice_sites.set_type(site_index, old_type);
    return 0.0;
  }
}

bool is_move_accepted(double delta_e, double T,
                      model_parameters_struct &parameters) {
  if (delta_e < 0) {
    return true;
  } else {
    real_dist proba_dist(0, 1);
    double boltzmann_factor{std::exp(-delta_e / T)};
    return boltzmann_factor > proba_dist(parameters.rng);
  }
}

} // namespace lattice_particles_space
