#include "lattice_particles_update.h"
#include <stdexcept>

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

int select_random_full(state_struct& state, model_parameters_struct &parameters) {
    // Using the new C++ 20 tools for lightweight access of arrays based on
    // condition
    // See https://learn.microsoft.com/en-us/cpp/standard-library/ranges?view=msvc-170
}

} // namespace lattice_particles_space
