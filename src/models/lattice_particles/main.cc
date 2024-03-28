#include "lattice_particles_parameters.h"
#include "lattice_particles_state.h"

#include <iostream>

int main() {
    lattice_particles_space::model_parameters_struct params(
        "../../../input/model_params.json");
    std::cout << params;
    return 0;
}
