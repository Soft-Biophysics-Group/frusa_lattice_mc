import numpy as np
from io_tools import *

#Define model parameters

model_params = {}

#Number of fields components
model_params["ns"] = 2

#Lattice dimensions
model_params["Lx"] = 100
model_params["Ly"] = 1
model_params["Lz"] = 1

#Number of particles
model_params["Np"] = 20

#Model parameters for a chain system
k11 = -1
k12 = 0
k21 = 0
T_model = 0.1
model_params["couplings"] = [T_model,k11,k12,k21,k11]

# Initialization option
model_params["initialize_option"] = "random"

# If "initialize_option" is set to "from_file", we must specify the location of
# the input file
# model_params["state_input"] = "./..."

# Options for average collection
model_params["e_av_option"] = True

if model_params["e_av_option"]:
    make_dir("../data/mf_energy_moments")
    model_params["e_av_output"] = "./data/mf_energy_moments/"

make_json_file(model_params,"../input/model_params.json")

#Define mc parameters

mc_params = {}

#Number of MC steps used for equilibration
mc_params["mcs_eq"] = 100

#Number of MC steps used for averaging
mc_params["mcs_av"] = 1

#Type of cooling schedule
mc_params["cooling_schedule"] = "exponential"

#Initial annealing temperature
mc_params["Ti"] = 0

#Final annealing temperature
mc_params["Tf"] = -5

#Number of annealing steps
mc_params["Nt"] = 10

#Option to collect state checkpoints at the end of each temperature cycle
mc_params["checkpoint_option"] = False

#If checkpoint is True, we need to provide the output address for the 
#checkpoint files
if mc_params["checkpoint_option"]:
    mc_params["checkpoint_address"] = ""

# Output location of the final state configuration (must end with "/")
mc_params["final_structure_address"] = "./"

make_json_file(mc_params,"../input/mc_params.json")
