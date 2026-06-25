#include "particles_acceptance.h"

#include "io_utils.h"

namespace particles_space
{

void record_acceptance(acceptance_struct& acceptance, double T)
{
  vec1d row {T};

  double total_attempted {0.0};
  double total_accepted {0.0};
  for (std::size_t move {0}; move < mc_moves::n_enum_moves; ++move) {
    double att {acceptance.attempted[move]};
    double acc {acceptance.accepted[move]};
    row.push_back(acc / att);
    total_attempted += att;
    total_accepted += acc;
  }

  acceptance.overall_acceptance_rate = total_accepted / total_attempted;
  row.push_back(acceptance.overall_acceptance_rate);

  acceptance.history.push_back(row);
}

void save_acceptance(acceptance_struct& acceptance,
                     model_parameters_struct& parameters)
{
  if (!parameters.acceptance_option)
    return;

  // Written directly rather than via save_vector to prepend a column header:
  // the table is wide and unused moves show up as nan, so labels are needed.
  std::ofstream out {parameters.acceptance_output};
  if (!out) {
    std::cerr << "Unable to open file " << parameters.acceptance_output << '\n';
    exit(1);
  }

  out << "# T";
  for (const std::string& name : mc_moves_str)
    out << ' ' << name;
  out << " overall\n";

  out << std::setprecision(8);
  for (const vec1d& row : acceptance.history) {
    for (std::size_t col {0}; col < row.size(); ++col)
      out << (col == 0 ? "" : " ") << row[col];
    out << '\n';
  }
}

}  // namespace particles_space
