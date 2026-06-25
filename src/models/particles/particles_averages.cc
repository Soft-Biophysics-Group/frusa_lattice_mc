#include "particles_averages.h"

#include "io_utils.h"

namespace particles_space
{

/*
 * Definitions required for the public routines of the model class
 */
void initialize_averages(averages_struct& averages,
                         model_parameters_struct& parameters)
{
  if (parameters.e_av_option == true) {
    averages.e_av = 0.0;
    averages.e2_av = 0.0;
    averages.e_series.clear();
  }
}

void update_averages(averages_struct& averages,
                     [[maybe_unused]] state_struct& state,
                     interactions_struct& interactions,
                     model_parameters_struct& parameters,
                     [[maybe_unused]] double T)
{
  if (parameters.e_av_option == true) {
    averages.e_av += interactions.energy;
    averages.e2_av += interactions.energy * interactions.energy;
    averages.e_series.push_back(interactions.energy);
  }
}


// Calculate autocorrelation time according to the method developed by Sokal
// Returns 1.0 (no autocorrelation) if not enough samples or negative variance
// See https://emcee.readthedocs.io/en/stable/tutorials/autocorr/
// and Sokal 1997, Monte-Carlo Methods in Statistical Physics
static double calc_autocorr_time(const vec1d& e_series)
{
  const std::size_t n_steps {e_series.size()};
  double factor_window {5.0};

  if (n_steps < static_cast<std::size_t>(factor_window))
    return 1.0;

  double mean_e {0.0};
  for (double e : e_series)
    mean_e += e;
  mean_e /= static_cast<double>(n_steps);

  double var_e {0.0};
  for (double e : e_series)
    var_e += (e - mean_e) * (e - mean_e);
  var_e /= static_cast<double>(n_steps);

  double autocorr_time {0.0};
  double sum {0.0};
  for (std::size_t t = 1; t < n_steps; ++t) {
    // Calculate the normalized correlation function
    double corr_fct_of_t {0.0};
    for (std::size_t i = 0; i < n_steps - t; ++i) {
      corr_fct_of_t += (e_series[i] - mean_e) * (e_series[i + t] - mean_e);
    }
    corr_fct_of_t /= static_cast<double>(n_steps - t) * var_e;
    sum += corr_fct_of_t;

    // See docstring link for derivation
    autocorr_time = 1 + 2 * sum;
    if (static_cast<double>(t) > factor_window * autocorr_time)
      break;
  }
  return autocorr_time;
}

void save_averages(averages_struct& averages,
                   [[maybe_unused]] state_struct& state,
                   model_parameters_struct& parameters,
                   double T,
                   int mcs_av)
{
  if (parameters.e_av_option == true) {
    double autocorr_time {calc_autocorr_time(averages.e_series)};
    averages.e_av /= mcs_av;
    averages.e2_av /= mcs_av;

    vec1d output_vec = {T, averages.e_av, averages.e2_av, autocorr_time};
    std::string output_file =
        parameters.e_av_output + "esf_av_T_" + std::to_string(T) + ".dat";
    io_space::save_vector(output_vec, 4, output_file);
  }
  averages.e_series.clear();
}

/*
 * End of the required definitions for the model class
 */
}  // namespace particles_space
