#include "lattice_particles_parameters.h"
#include "lattice_particles_state.h"

#include <iostream>

int main() {
    lattice_particles_space::model_parameters_struct params(
        "input/model_params.json");
    std::cout << params;

    lattice_particles_space::state_struct state {};
    lattice_particles_space::initialize_state(state, params);
    lattice_particles_space::save_state(state, "state.txt");

    lattice_particles_space::model_parameters_struct params_load(
        "input/model_params_new.json");
    lattice_particles_space::state_struct state_new {};
    lattice_particles_space::initialize_state(state_new, params_load);
    lattice_particles_space::save_state(state_new, "state_new.txt");
    return 0;
}
