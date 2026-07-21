"""
Vincent Ouazan-Reboul
2025

For now, just one function: analyze the energy records of our MC simulations
"""

from pathlib import Path
import numpy as np
from numpy.typing import NDArray

class Averages:
    temperatures: NDArray[np.float64]
    energies: NDArray[np.float64]
    squared_energies: NDArray[np.float64]
    autocorr_times: NDArray[np.float64]

    def __init__(self, averages_folder: Path | str):
        all_T = []
        all_e_records = []
        all_e_sq_records = []
        all_autocorr_times = []

        averages_path = Path(averages_folder)
        for f in averages_path.glob("esf_av_T_*.dat"):
            this_record = np.loadtxt(f)
            all_T.append(this_record[0])
            all_e_records.append(this_record[1])
            all_e_sq_records.append(this_record[2])
            all_autocorr_times.append(this_record[3])

        # Sort by decreasing order of temperature
        sorting_inds = np.argsort(all_T)[::-1]
        self.temperatures = np.array(all_T)[sorting_inds]
        self.energies = np.array(all_e_records)[sorting_inds]
        self.squared_energies = np.array(all_e_sq_records)[sorting_inds]
        self.autocorr_times = np.array(all_autocorr_times)[sorting_inds]

    @property
    def heat_capacity(self):
        variance = self.squared_energies - self.energies**2
        return variance / self.temperatures**2
