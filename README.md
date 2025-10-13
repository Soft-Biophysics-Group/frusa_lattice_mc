# 3DSizeLimitedAssembly

Codes and executable to run defect engineered three-dimensional assemblies introduced in [arXiv:2504.13073](https://doi.org/10.48550/arXiv.2504.13073).

## Contents

- A notebook, `example_simulation.ipynb`, which shows the basics of running a simulation
- A `scripts` folder to generate the input files for reproducing the results of Fig. 4H of
the article
- In `src/` and `include/`, the C++ files which yield the executable `frusa_mc` used to
actually run
simulations
- In `python/`, additional python code consisting of wrappers around the C++ code, utilities
  to generate interaction maps, and 3D visualization tools.
- In `input`, input files for reproducing the results of Fig.4H are supplied. Please note that
the simulations took three days on a high-performance cluster.
- In `data`, simulation results
- In `scripts`, files to generate the input files shown in `input` and to visualize the
simulation results are supplied. These can be ran using `uv` or enabling the virtual
environment first.

## Dependencies

The following dependencies are needed to use the code contained in this repository:

- Python version 3.11
- CMake for compiling the C++ code
- Blender for visualization (version 3.4.3 was used, which can be obtained [here](https://download.blender.org/release/Blender4.3/))
- We strongly advise using the [uv](https://docs.astral.sh/uv/) python package manager to
install the Python dependencies of the program and control the Python version used to run
notebooks/scripts.

## Installation

Either download this code as a .zip, or clone it using:

`git clone https://github.com/Soft-Biophysics-Group/frusa_lattice_mc.git -b VortexAssembly3D`

With all the dependencies installed, run: `make all`.
If you use uv, then run: `uv sync`

## A note on virtual environments

As some of the Python packages used in this code are heavy to install and may be subject to
change with new version of Python / Blender, we heavily recommend using a Python virtual
environment to maximize code reproducibility.

If you use the uv package manager, you can easily run the code in its virtual environment
with the command `uv run file.py` for running Python scripts, and `uv run jupyter
notebook` or `uv run jupyter lab` for running Jupyter notebooks.

Otherwise, you can activate the virtual environment from this folder with the Unix command
`source .venv/bin/activate`, after which you can run `python` and `jupyter` commands as
usual.

## Usage
A jupyter notebook, `example_simulation.ipynb`, is supplied as a minimal working example of
the code.
To run it, if you use uv, run `uv run jupyter lab` and select the notebook from the Jupyter
interface.
If you don't, after activating the virtual environment as outlined above, run `python -m
jupyter lab` and select the notebook from the jupyter interface.

The code can also be run from its executable `frusa_mc`. It takes two json input files, whose
generation is outlined in `example_simulation.ipynb`. The syntax for running it is
`./frusa_mc -m {path_to_model_file.json} -M {path_to_MC_file.json}`. The code can also be ran
with the thin wrapper `run_simu` contained in the `python/src/config.py` file, which could
help with compatibility if you are using a Windows system.
The output of the program is typically contained in the `data/` folder, and specified in the
`model_file.json` and `mc_file.json`.

### Input
The program takes two json files as an input, one containing the parameters related to the
model (lattice properties, interaction map between particles...) and the other to the Monte
Carlo parameters (temperature profile, number of steps...).
See the example notebook for the contents and generation of these files.

### Output
The results include, depending on what you want to include:
- an `average_e/` folder, containing one file per temperature step, whose contents consist of:
```
[Temperature]
[Average total energy]
[Energy second moment]
```
- an `e_record/` folder, with one file per temperature step. The first line of each file is
the temperature at which samples are recorded. Each following line corresponds to the energy
after each Monte-Carlo step taken at this temperature. Note that, for a large number of
steps, this creates very large files!
- a `structures/` folder, with one file per temperature step. Each file corresponds to the
final lattice configuration at the end of a temperature step. The file consists of two
lines, the first being the type of particle on each lattice site (always 0 here), the other
the orientation of the particle on each site, with orientation -1 corresponding to an empty
site. At the very end of the simulation, a final `final_structure.dat` file is recorded.
