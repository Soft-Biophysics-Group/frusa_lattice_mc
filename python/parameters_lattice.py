import numpy as np
import os
from json_dump import *
from lattice_helper_classes import *
from hexagonal_helper_functions import *

# Define model parameters

model_params = {}

model_params["n_types"] = 1
model_params["n_orientations"] = 6

# Lattice dimensions
model_params["lx"] = 5
model_params["ly"] = 5
model_params["lz"] = 1

# Number of particles
model_params["n_particles"] = [10, 10]

# Model parameters

crystal_int = 0.5
line_int = -8
surf_tension = 6
repulsion = 15

model_params["couplings"] = camembert_couplings_hexagonal(
    crystal_int, line_int, surf_tension, repulsion
)

model_params["initialize_option"] = "random_fixed_particle_numbers"

# # swap_empty_full_enum,
# swap_full_full_enum,
# rotate_enum,
# mutate_enum,
model_params["move_probas"] = [
    1 / 2,  # swap_empty_full
    0,      # swap_full_full
    1 / 2,  # rotate
    0,      # mutate
]


make_json_file(model_params, "../input/model_params.json")

# Define mc parameters

mc_params = {}

# Number of MC steps used for equilibration
mc_params["mcs_eq"] = 1000

# Number of MC steps used for averaging
mc_params["mcs_av"] = 1

# Type of cooling schedule
# options: linear, exponential
mc_params["cooling_schedule"] = "linear"

# Initial annealing temperature
mc_params["Ti"] = 6

# Final annealing temperature
mc_params["Tf"] = 1

# Number of annealing steps
mc_params["Nt"] = 100

# Option to collect state checkpoints at the end of each temperature cycle
mc_params["checkpoint_option"] = False

# If checkpoint is True, we need to provide the output address for the
# checkpoint files
if mc_params["checkpoint_option"]:
    mc_params["checkpoint_address"] = ""

make_json_file(mc_params, "../input/mc_params.json")
