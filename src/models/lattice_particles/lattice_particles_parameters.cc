#include "lattice_particles_parameters.h"
#include "vector_utils.h"

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
    n_particles    = json_model_params["n_particles"].template get<ContactMap>();
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

std::ostream& operator<< (std::ostream& out, model_parameters_struct& params) {
    out << "Printing the parameters for this simulation\n\n";
    out << "Number of particle types:" << params.n_types << '\n';
    out << "Number of possible particle orientations:" << params.n_orientations << '\n';
    out << "Lattice dimensions:" << '(' << params.lx << ", " << params.ly
        << ", " << params.lz << ")\n";
    out << "Number of particles of each type:" << '[';
    array_space::print_vector(out, params.n_particles);
    out << "]\n";
    out << "Flattened couplings matrix: [";
    array_space::print_vector(out, params.couplings);
    out << "]\n";
    out << "Chosen initialize option: " << params.initialize_option << '\n';
    if (params.initialize_option == "from_file") 
        out << "Initialized from file " << params.state_input << '\n';

    out << "End parameter print \n\n";

    return out;
}
} // namespace lattice_particles_space
