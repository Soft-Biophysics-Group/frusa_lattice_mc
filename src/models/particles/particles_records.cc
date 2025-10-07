#include "particles_records.h"


namespace particles_space
{
void update_records(model_parameters_struct& parameters,
                    interactions_struct& interactions,
                    records_struct& records)
{
  if (parameters.e_record_option == true) {
    records.e_records.push_back(interactions.energy);
  }
}

void save_records(model_parameters_struct& parameters,
                  double T,
                  records_struct& records)
{
  if (parameters.e_record_option == true) {
    vec1d output_vec {T};
    for (double el : records.e_records) {
      output_vec.push_back(el);
    }
    std::string output_file =
        parameters.e_record_output + "e_record_T_" + std::to_string(T) + ".dat";
    int out_size {static_cast<int>(output_vec.size())};
    io_space::save_vector(output_vec, out_size, output_file);
    records.e_records.clear();
  }
}

}  // namespace particles_space
