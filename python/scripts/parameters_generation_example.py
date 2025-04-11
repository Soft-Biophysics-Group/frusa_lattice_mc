"""
Vincent Ouazan-Reboul, 2025

Example python file to generate a .json input file for the frusa_mc program.
Simulates a few particles with uniform interactions on a triangular lattice.
"""

from json_dump import *
import contact_utils as cu
import config as cfg
from pathlib import Path
### Key parameter: how we will name this run
run_name = "test_run"

### Define model parameters

model_params = {}

# ---------- LATTICE OPTIONS ----------

# Options:
# "chain", "square", "triangular", "cubic", "bcc", "fcc"
model_params["lattice_name"] = "triangular"

# Lattice dimensions
model_params["lx"] = 40
model_params["ly"] = 40  # Has to be 1 for chain
model_params["lz"] = 1  # Has to be 1 for square & triangular

# ---------- MODEL PARAMETERS ----------

# Number of particle types
model_params["n_types"] = 1

# Number of particles of each type
model_params["n_particles"] = [400]

# Couplings is its own beast. Should be gotten with the appropriate helper
# function.

# Parameters I got from the draft of Lara's paper
e_crystal = -18.7
e_defect = 1/0.55 * e_crystal
e_repel = 10

# I even made a specific function for these couplings!
model_params["couplings"] = cu.get_camembert_cmap(e_crystal, e_defect, e_repel)

# Initialization option
model_params["initialize_option"] = "random"

# If "initialize_option" is set to "from_file", we must specify the location of
# the input file
# model_params["state_input"] = str(cfg.structures_path/"final_structure.dat")

# Options for average and record collection

model_params["state_av_option"] = True
model_params["e_av_option"] = False
model_params["e_record_option"] = True

if model_params["state_av_option"]:
    state_path = Path("./data/"+run_name+"/average_state")
    state_path.mkdir(parents = True, exist_ok = True)
    model_params["state_av_output"] =  str(state_path.resolve()) + "/"

if model_params["e_av_option"]:
    energy_path = Path("./data/"+run_name+"/average_energy")
    energy_path.mkdir(parents = True, exist_ok = True)
    model_params["e_av_output"] = str(energy_path.resolve())  + "/"

if model_params["e_record_option"]:
    e_records_path = Path("./data/" + run_name + "/energy_records")
    e_records_path.mkdir(parents=True, exist_ok=True)
    model_params["e_record_output"] = str(e_records_path.resolve()) + "/"

# Pick the probabilities of different moves. Has to sum to 1.
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

make_json_file(model_params, cfg.input_path/"test_model_params.json")

### Define mc parameters

mc_params = {}

# Number of MC steps used for equilibration
# Lara used 12k steps per particle, in our program the number of steps is given per site
# so for a start 12k * 400 / 40^2
mc_params["mcs_eq"] = 3000

# Number of MC steps used for averaging
mc_params["mcs_av"] = 100

# Type of cooling schedule
# Options: "exponential", "linear", "inverse"
# if exponential chosen: specify log10(T) as initial and final temperatures
# if inverse chosen: specify 1 / T as initial and final temperatures, in decreasing order.
mc_params["cooling_schedule"] = "inverse"

# Initial annealing temperature
mc_params["Ti"] = 0

# Final annealing temperature
mc_params["Tf"] = 1

# Number of annealing steps
mc_params["Nt"] = 200

# Option to collect state checkpoints at the end of each temperature cycle
mc_params["checkpoint_option"] = True

# If checkpoint is True, we need to provide the output address for the 
# checkpoint files
structures_path = Path("./data/"+run_name+"/structures")
structures_path.mkdir(parents = True, exist_ok = True)

if mc_params["checkpoint_option"]:
    mc_params["checkpoint_address"] = str(structures_path.resolve())+"/"

# Output location of the final state configuration (must end with "/")
mc_params["final_structure_address"] = str(structures_path.resolve())+"/"

make_json_file(mc_params, cfg.input_path/"test_mc_params.json")
