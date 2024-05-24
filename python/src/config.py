#! /usr/bin/env python3
# Vincent Ouazan-Reboul, 2024.05.14
# Set of utilities which the Python code will have to use frequently

from pathlib import Path

# Reliably make absolute paths to the right place.
# #ILovePathLib <3 <3 <3
parent_path = Path(__file__).parent.parent.parent.absolute().resolve()

input_path = parent_path/"input/"

data_path = parent_path/"data/"
averages_path = data_path/"average_state/"
energy_path = data_path/"energy_moments/"
structures_path = data_path/"structures/"

exec_path = parent_path/"build/app/frusa_mc"
