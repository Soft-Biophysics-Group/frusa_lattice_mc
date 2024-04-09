#include "lattice_particles_update.h"
#include "lattice_particles_geometry.h"
#include "lattice_particles_interactions.h"
#include "lattice_particles_parameters.h"
#include "lattice_particles_state.h"
#include <stdexcept>
#include <ranges>

namespace lattice_particles_space {

// TODO Figure out where to stick move_probas
template <int N>
void update_system(state_struct &state, interactions_struct<N> &interactions,
                   model_parameters_struct &parameters, double T, move_probas_vec move_probas) {
  // Pick the kind of move we'll be making
  // TODO Initialize move_probas somewhere
  mc_moves chosen_move { pick_random_move(move_probas, parameters) };
  switch (chosen_move) {
    case mc_moves::swap_empty_full:
        interactions.energy += swap_empty_full(state, parameters, interactions, T);
    case mc_moves::swap_full_full:
        interactions.energy += swap_full_full(state, parameters, interactions, T);
    case mc_moves::rotate:
        interactions.energy += rotate(state, parameters, interactions, T);
    case mc_moves::mutate:
        interactions.energy += mutate(state, parameters, interactions, T);
    default:
        throw std::runtime_error("Something went wrong in the move selection");
  }
}

move_probas_vec get_move_probas(std::string &mc_json_file) {
    std::ifstream mc_json_f{mc_json_file};
    if (!mc_json_f) {
        std::cerr << "Could not find MC paramereters file" << std::endl;
        exit(1);
    }
    json mc_json{json::parse(mc_json_file)};
    move_probas_vec move_probas{
        mc_json["move_probas"].template get<move_probas_vec>()};
    // Move probabilities have to sum to 1
    if ((std::accumulate(move_probas.begin(), move_probas.end(), 0.) != 1.0))
        throw std::runtime_error("Move probabilities do not sum to 1!");

    return move_probas;
}

mc_moves pick_random_move(move_probas_vec move_probas,
                          model_parameters_struct &parameters) {
  real_dist move_dist {0.0, 1.0};
  double sampled_real { move_dist(parameters.rng) };
  // Standard tower sampling algorithm
  int move_index {0};
  double cumulative_prob {move_probas[0]};
  while (cumulative_prob < sampled_real) {
    ++move_index;
    cumulative_prob += move_probas[move_index];
  }

  return static_cast<mc_moves>(move_index);
}

//int select_random_full_index(state_struct &state,
                               //model_parameters_struct &parameters) {
  //// Using the new C++ 20 tools for lightweight access of arrays based on
  //// condition
  //// See
  //// https://learn.microsoft.com/en-us/cpp/standard-library/ranges?view=msvc-170
  //// TODO Test this and make sure it works
  //auto full_sites {state.lattice_sites
      //| std::views::filter(
              //[](site_state& site) { return site.is_empty(); };
              //)
  //};
  //// Generate a random index with the right upper bound
  //int_dist full_site_dist { int_dist(0, std::size(full_sites)) };
  //int index_in_view { full_site_dist(parameters.rng) };

  //// TODO See if I need to plug a view adaptor at this point
  //return full_sites[index_in_view].get_site_index();
//}

//int select_random_empty_site_index(state_struct &state,
                               //model_parameters_struct &parameters) {
  //// Using the new C++ 20 tools for lightweight access of arrays based on
  //// condition
  //// See
  //// https://learn.microsoft.com/en-us/cpp/standard-library/ranges?view=msvc-170
  //// TODO Test this and make sure it works
  //auto empty_sites {state.lattice_sites
      //| std::views::filter(
              //[](site_state& site) { return !site.is_empty(); };
              //)
  //};
  //// Generate a random index with the right upper bound
  //int_dist empty_site_dist { int_dist(0, std::size(empty_sites)) };
  //int index_in_view { empty_site_dist(parameters.rng) };

  //// TODO See if I need to plug a view adaptor at this point
  //return empty_sites[index_in_view].get_site_index();
//}

template <int N>
double swap_sites(site_state &site1, site_state &site2, state_struct &state,
                  model_parameters_struct &parameters,
                  interactions_struct<N> &interactions, double T) {
  double energy_change {0.0};
  bool sites_are_neighbours{are_neighbours(
      site1.get_site_index(), site2.get_site_index(), interactions, state)};
  // Initial energy, possibly including neighbour correction
  energy_change -=
      get_site_energy(site1, state, interactions.couplings) +
      get_site_energy(site2, state, interactions.couplings);
  if (sites_are_neighbours)
    energy_change += neighbour_correction(site1, site2,
                                          interactions, state, parameters);
  // Make move and calculate energy after
  site1.swap_with(site2);
  energy_change +=
      get_site_energy(site1, state, interactions.couplings) +
      get_site_energy(site2, state, interactions.couplings);
  if (sites_are_neighbours)
    energy_change -=
        neighbour_correction(site1, site2, interactions, state, parameters);
  // Accept or reject move
  if (move_accepted(energy_change, T))
    return energy_change;
  else {
    site1.swap_with(site2);
    return 0.0;
  }
}

template <int N>
double swap_empty_full(state_struct &state, model_parameters_struct &parameters,
                       interactions_struct<N> &interactions, double T) {
  int empty_site_index{select_random_empty_index(state, parameters)};
  int full_site_index{select_random_full_index(state, parameters)};
  site_state& empty_site {state.lattice_sites[empty_site_index]};
  site_state& full_site {state.lattice_sites[full_site_index]};
  return swap_sites(empty_site, full_site, state, parameters, interactions, T);
}

template <int N>
double swap_full_full(state_struct &state, model_parameters_struct &parameters,
                       interactions_struct<N> &interactions, double T) {
  int site_1_index{select_random_full_index(state, parameters)};
  int site_2_index{select_random_full_index(state, parameters)};
  site_state& site_1 {state.lattice_sites[site_1_index]};
  site_state& site_2 {state.lattice_sites[site_2_index]};
  return swap_sites(site_1, site_2, state, parameters, interactions, T);
}

template <int N>
double neighbour_correction(int site_1_index, int site_2_index,
                            interactions_struct<N> &interactions,
                            state_struct &state,
                            model_parameters_struct& parameters) {
  // TODO ALl of this is unholy. I have to do better once things are up and
  // running
  int edge_1{get_bond_direction(site_1_index, site_2_index, parameters.lx,
                                parameters.ly, parameters.lz)};
  int edge_2{get_bond_direction(site_2_index, site_1_index, parameters.lx,
                                parameters.ly, parameters.lz)};
  return get_contact_energy(state.lattice_sites[site_2_index],
                            state.lattice_sites[site_1_index], edge_1, edge_2,
                            interactions.couplings, state.n_states);
}

template <int N>
bool are_neighbours(int site_1_index, int site_2_index,
                    interactions_struct<N> &interactions, state_struct &state) {
  get_neighbours(interactions.neighbours, site_1_index, state);
  // Let's look for site_2 in the neighbours of site_1
  // If it's not in there, then the sites are not neighbours
  auto neigh_index{std::find(interactions.neighbours.begin(),
                             interactions.neighbours.end(), site_2_index)};
  return neigh_index != interactions.neighbours.end();
}

template <int N>
double rotate(state_struct &state, model_parameters_struct &parameters,
              interactions_struct<N> &interactions, double T) {
  int site_index {select_random_full_index(state, parameters)};
  site_state& site { state.lattice_sites[site_index] };
  double energy_change { - get_site_energy(site, state, interactions) };

  // Lower bound for orientation is 1, as orientation 0 corresponds
  // to an empty site
  int_dist rot_dist { 1, state.n_orientations };
  int new_orientation { rot_dist(parameters.rng) };
  int old_orientation { site.get_orientation() };
  site.set_state(site.get_type(), new_orientation);
  energy_change += get_site_energy(site, state, interactions);
  if (move_accepted(energy_change, T))
    return energy_change;
  else {
    site.set_state(site.get_type(), old_orientation);
    return 0.0;
  }
}

template <int N>
double mutate(state_struct &state, model_parameters_struct &parameters,
              interactions_struct<N> &interactions, double T) {
  int site_index{select_random_full_index(state, parameters)};
  site_state &site{state.lattice_sites[site_index]};
  double energy_change{-get_site_energy(site, state, interactions)};

  // Lower bound for orientation is 1, as orientation 0 corresponds
  // to an empty site
  int_dist type_dist{1, state.n_types};
  int new_type{type_dist(parameters.rng)};
  int old_type{site.get_type()};
  site.set_state(new_type, site.get_orientation());
  energy_change += get_site_energy(site, state, interactions);
  if (move_accepted(energy_change, T))
    return energy_change;
  else {
    site.set_state(old_type, site.get_orientation());
    return 0.0;
  }
}

} // namespace lattice_particles_space
