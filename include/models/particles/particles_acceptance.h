#ifndef PARTICLES_ACCEPTANCE_H
#define PARTICLES_ACCEPTANCE_H

#include <array>

#include "particles_parameters.h"
#include "vector_utils.h"

namespace particles_space
{

struct acceptance_struct
{
  std::array<double, mc_moves::n_enum_moves> attempted {};
  std::array<double, mc_moves::n_enum_moves> accepted {};

  double overall_acceptance_rate {0.0};

  // Set by each attempt_*, read back by update_system to attribute the outcome
  // to the move type (avoids mistaking an accepted delta_e == 0 for a reject).
  bool last_accepted {false};

  // One row per temperature: T, per-move rate (mc_moves order), overall rate.
  vec2d history {};

  // Clears the per-temperature counters but keeps history.
  void reset()
  {
    attempted.fill(0);
    accepted.fill(0);
    overall_acceptance_rate = 0;
    last_accepted = false;
  }
};

// Appends the rates measured at temperature T to history.
void record_acceptance(acceptance_struct& acceptance, double T);

// Writes history to parameters.acceptance_output; no-op unless
// parameters.acceptance_option is set.
void save_acceptance(acceptance_struct& acceptance,
                     model_parameters_struct& parameters);

}  // namespace particles_space

#endif
