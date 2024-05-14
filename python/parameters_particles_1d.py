#! /usr/bin/env python3
from pathlib import Path
import numpy as np
from json_dump import *
import contact_utils as cu

### Define model parameters

model_params = {}

# ---------- LATTICE OPTIONS ----------

# Options:
# "chain", "square", "triangular", "cubic", "bcc", "fcc"
model_params["lattice"] = "chain"

# Lattice dimensions
model_params["lx"] = 100
model_params["ly"] = 1  # Has to be 1 for chain
model_params["lz"] = 1  # Has to be 1 for square & triangular

# ---------- MODEL PARAMETERS ----------

# Number of particle types
model_params["n_types"] = 1

# Number of particles of each type
model_params["n_particles"] = [20]

# Couplings is its own beast. Should be gotten with the appropriate helper
# function.

# For chain
mat_11 = cu.chain_LEL_1type(0.0, 0.0, 0.0)
model_params["couplings"] = cu.flatten_couplings(mat_11)
# mat_21 = np.array([[0.0, 0.0], [0.0, 0.0]])
# mat_22 = cu.chain_LEL_1type(0.0, 0.0, 0.0)
# model_params["couplings"] = cu.chain_LEL_2types(mat_11, mat_21, mat_22)

# Initialization option
model_params["initialize_option"] = "random"

# If "initialize_option" is set to "from_file", we must specify the location of
# the input file
# model_params["state_input"] = "./..."

# Options for average collection
model_params["state_av_option"] = True
model_params["e_av_option"] = True

data_path = Path("../data")
averages_path = data_path/"average_state"
if model_params["state_av_option"]:
    averages_path.mkdir(parents=True,exist_ok=True)
    model_params["state_av_output"] = str(averages_path)

data_path = Path("../data")
energy_path = data_path/"energy_moments"
if model_params["e_av_option"]:
    energy_path.mkdir(parents=True, exist_ok=True)
    model_params["e_av_output"] = "./data/energy_moments/"

# Pick the probabilities of different moves
# Options:
""""
    "swap_empty_full",
    "swap_full_full",
    "rotate",
    "mutate",
    "rotate_and_swap_w_empty"
"""

moves_dict = {}
moves_dict["swap_empty_full"] = 1/3
moves_dict["rotate"] = 1/3
moves_dict["rotate_and_swap_w_empty"] = 1/3

model_params["move_probas"] = moves_dict

make_json_file(model_params,"../input/model_params.json")

### Define mc parameters

mc_params = {}

# Number of MC steps used for equilibration
mc_params["mcs_eq"] = 100

# Number of MC steps used for averaging
mc_params["mcs_av"] = 1

# Type of cooling schedule
mc_params["cooling_schedule"] = "exponential"

# Initial annealing temperature
mc_params["Ti"] = 0

# Final annealing temperature
mc_params["Tf"] = -5

# Number of annealing steps
mc_params["Nt"] = 10

# Option to collect state checkpoints at the end of each temperature cycle
mc_params["checkpoint_option"] = False

# If checkpoint is True, we need to provide the output address for the 
# checkpoint files
if mc_params["checkpoint_option"]:
    mc_params["checkpoint_address"] = ""

# Output location of the final state configuration (must end with "/")
mc_params["final_structure_address"] = "./"

make_json_file(mc_params,"../input/mc_params.json")
