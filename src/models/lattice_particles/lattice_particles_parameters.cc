#include "lattice_particles_parameters.h"

namespace lattice_particles_space {

model_parameters_struct::model_parameters_struct(const std::string& input_file) {

    std::ifstream model_input_f {input_file};
    if (!model_input_f) {
        std::cerr << "Could not open JSON model parameters file" << std::endl;
        exit(1);
    }

    json json_model_params = json::parse(model_input_f);

    n_types        = json_model_params["n_types"].template get<int>();
    n_orientations = json_model_params["n_orientations"].template get<int>();
    lx             = json_model_params["lx"].template get<int>();
    ly             = json_model_params["ly"].template get<int>();
    lz             = json_model_params["lz"].template get<int>();
    n_particles    = json_model_params["n_particles"].template get<vec1i>();
    couplings      = json_model_params["couplings"].template get<vec1d>();
    initialize_option =
        json_model_params["initialize_option"].template get<std::string>();

    if (initialize_option == "from_file") {
        state_input =
            json_model_params["state_input"].template get<std::string>();
    }

    std::random_device dev;
    EngineType engine(dev());

    rng = engine;
}
} // namespace lattice_particles_space
