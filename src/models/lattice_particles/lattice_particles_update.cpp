#include "lattice_particles_update.h"
#include "lattice_particles_state.h"
#include <stdexcept>
#include <ranges>

namespace lattice_particles_space {

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

site_state select_random_full(state_struct &state,
                               model_parameters_struct &parameters) {
  // Using the new C++ 20 tools for lightweight access of arrays based on
  // condition
  // See
  // https://learn.microsoft.com/en-us/cpp/standard-library/ranges?view=msvc-170
  // TODO Test this and make sure it works
  auto full_sites {state.lattice_sites
      | std::views::filter(
              [](site_state& site) { return site.is_empty(); };
              )
  };
  // Generate a random index with the right upper bound
  int_dist full_site_dist { int_dist(0, std::size(full_sites)) };
  int site_index { full_site_dist(parameters.rng) };

  // TODO See if I need to plug a view adaptor at this point
  return full_sites[site_index];
}

site_state select_random_empty(state_struct &state,
                               model_parameters_struct &parameters) {
  // Using the new C++ 20 tools for lightweight access of arrays based on
  // condition
  // See
  // https://learn.microsoft.com/en-us/cpp/standard-library/ranges?view=msvc-170
  // TODO Test this and make sure it works
  auto empty_sites {state.lattice_sites
      | std::views::filter(
              [](site_state& site) { return !site.is_empty(); };
              )
  };
  // Generate a random index with the right upper bound
  int_dist empty_site_dist { int_dist(0, std::size(empty_sites)) };
  int site_index { empty_site_dist(parameters.rng) };

  // TODO See if I need to plug a view adaptor at this point
  return empty_sites[site_index];
}

double swap_empty_full(state_struct &state, model_parameters_struct &parameters) {
  site_state& empty_site { select_random_empty(state, parameters) };
  site_state& full_site { select_random_full(state, parameters) };
  // Determine if the sites are neighbours
}
void swap_full_full(state_struct &state, model_parameters_struct &parameters);
void rotate(state_struct &state, model_parameters_struct &parameters);
void mutate(state_struct &state, model_parameters_struct &parameters);

} // namespace lattice_particles_space
