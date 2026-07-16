"""
Vincent Ouazan-Reboul
2025

For now, just one function: analyze the energy records of our MC simulations
"""

from pathlib import Path
import numpy as np
from numpy.typing import NDArray

# Plotting the energy trace to see when the cube actually formed
def get_temperature_energy_records(
    records_folder: str | Path = ""
):
    all_T = []
    all_e_records = []
    records_path = Path(records_folder)
    for f in records_path.iterdir():
        this_record = np.loadtxt(f)
        all_T.append(this_record[0])
        all_e_records.append(this_record[1:])
    # Sort by decreasing order of temperature
    sorting_inds = np.argsort(all_T)[::-1]
    all_T_sorted = np.array(all_T)[sorting_inds]
    all_T_sorted = np.repeat(all_T_sorted, len(all_e_records[0]))
    all_e_sorted = np.array(all_e_records)[sorting_inds]

    return all_T_sorted, np.hstack(all_e_sorted)

def get_average_energies_and_squared(averages_folder: str | Path = ""):
    all_T = []
    all_avg_e = []
    all_avg_sq_e = []
    avgs_path = Path(averages_folder)
    for f in sorted(avgs_path.iterdir()):
        this_avg = np.loadtxt(f)
        all_T.append(this_avg[0])
        all_avg_e.append(this_avg[1])
        all_avg_sq_e.append(this_avg[2])
    sorting_inds = np.argsort(all_T)[::-1]
    all_T_sorted = np.array(all_T)[sorting_inds]
    all_e_sorted = np.array(all_avg_e)[sorting_inds]
    all_sq_e_sorted = np.array(all_avg_sq_e)[sorting_inds]
    return all_T_sorted, all_e_sorted, all_sq_e_sorted

def get_average_energies(averages_folder: str | Path = ""):

    all_T = []
    all_avg_e = []
    avgs_path = Path(averages_folder)
    for f in avgs_path.iterdir():
        this_avg = np.loadtxt(f)
        all_T.append(this_avg[0])
        all_avg_e.append(this_avg[1])
    # Sort by decreasing order of temperature
    sorting_inds = np.argsort(all_T)[::-1]
    all_T_sorted = np.array(all_T)[sorting_inds]
    all_e_sorted = np.array(all_avg_e)[sorting_inds]

    return all_T_sorted, all_e_sorted

def get_average_squared_energies(averages_folder: str | Path = ""):
    all_T = []
    all_avg_sq_e = []
    avgs_path = Path(averages_folder)
    for f in avgs_path.iterdir():
        this_avg = np.loadtxt(f)
        all_T.append(this_avg[0])
        all_avg_sq_e.append(this_avg[2])
    # Sort by decreasing order of temperature
    sorting_inds = np.argsort(all_T)[::-1]
    all_T_sorted = np.array(all_T)[sorting_inds]
    all_e_sorted:NDArray[np.float64] = np.array(all_avg_sq_e)[sorting_inds]

    return all_T_sorted, all_e_sorted
