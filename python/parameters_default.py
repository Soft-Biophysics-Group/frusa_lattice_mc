import numpy as np
from json_dump import * 

### Define model parameters

model_params = {}

# Number of particles
model_params["N"] = 20

# Energy gap
model_params["delta"] = 1

# Initialization option
model_params["initialize_option"] = "random"

# If "initialize_option" is set to "from_file", we must specify the location of
# the input file
# model_params["state_input"] = "./..."

make_json_file(model_params,"../input/model_params.json")

# Define mc parameters

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
