#include "particles_parameters.h"
#include "vector_utils.h"

namespace particles_space
{

model_parameters_struct::model_parameters_struct(const std::string& input_file)
{
  std::ifstream model_input_f {input_file};
  if (!model_input_f) {
    std::cerr << "Could not open JSON model parameters file" << '\n';
    exit(1);
  }

  json json_model_params = json::parse(model_input_f);

  // Parse parameters from JSON one by one
  n_types = json_model_params["n_types"].template get<int>();
  n_particles = json_model_params["n_particles"].template get<vec1i>();
  couplings = json_model_params["couplings"].template get<vec1d>();
  std::random_device dev;
  EngineType engine(dev());
  rng = engine;
  initialize_option =
      json_model_params["initialize_option"].template get<std::string>();
  if (initialize_option == "from_file") {
    state_input = json_model_params["state_input"].template get<std::string>();
  }
  move_probas = get_move_probas(input_file);
  e_av_option = json_model_params["e_av_option"].template get<bool>();
  if (e_av_option) {
    e_av_output = json_model_params["e_av_output"].template get<std::string>();
  }
  state_av_option = json_model_params["state_av_option"].template get<bool>();
  if (state_av_option) {
    state_av_output =
        json_model_params["state_av_output"].template get<std::string>();
  }
  e_record_option = json_model_params["e_record_option"].template get<bool>();
  if (e_record_option) {
    e_record_output =
        json_model_params["e_record_output"].template get<std::string>();
  }
}

std::ostream &operator<<(std::ostream &out, model_parameters_struct &params) {
  out << "Printing the parameters for this simulation\n\n";
  array_space::print_vector(out, params.n_particles);
  out << "]\n";
  out << "Flattened couplings matrix: [";
  array_space::print_vector(out, params.couplings);
  out << "]\n";
  out << "Chosen initialize option: " << params.initialize_option << '\n';
  if (params.initialize_option == "from_file") {
    out << "Initialized from file " << params.state_input << '\n';
  }
  out << "Move probabilities: ";
  // TODO Remove that ugly 6
  for (auto prob : params.move_probas) {
    out << prob << ", ";
  }
  out << '\n';

  out << "End parameter print \n\n";

  return out;
}

move_probas_arr get_move_probas(const std::string& model_input_file)
{
  std::ifstream mc_json_f {model_input_file};
  if (!mc_json_f) {
    std::cerr << "Could not find MC paramereters file" << '\n';
    exit(1);
  }

  json mc_json {json::parse(mc_json_f)};
  // This line is making the program fail on the cluster
  std::map<std::string, double> move_map {
      mc_json["move_probas"].template get<std::map<std::string, double>>()};

  move_probas_arr move_probas {};
  move_probas.fill(0.0);
  for (std::size_t move {0}; move < mc_moves::n_enum_moves; move++) {
    const std::string& move_name {mc_moves_str[move]};
    // look for entry with the name of the move and assign the right
    // probability if it exists
    if (move_map.contains(move_name)) {
      move_probas[move] = move_map[move_name];
    }
  }
  // Move probabilities have to sum to 1
  if ((std::accumulate(move_probas.begin(), move_probas.end(), 0.) != 1.0)) {
    throw std::runtime_error("Move probabilities do not sum to 1!");
  }

  return move_probas;
}
}  // namespace particles_space
