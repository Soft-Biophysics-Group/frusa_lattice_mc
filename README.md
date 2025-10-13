# `frusa_mc` - Generic Monte Carlo module for mean-field and particle models of frustrated self-assembly

## Publications using this code

- Lara Koehler, Markus Eder, Vincent Ouazan-Reboul, Christoph Karfusehr, Andrey Zelenskiy,
Pierre Ronceray, Friedrich C. Simmel, Martin Lenz: **Topological defect engineering enables
size and shape control in self-assembly**,
[arXiv:2504.13073](https://arxiv.org/abs/2504.13073).
  - [Link to the branch used for the publication](https://github.com/Soft-Biophysics-Group/frusa_lattice_mc/tree/VortexAssembly3D)

## Scope of this code

Simulated-annealing Metropolis Monte Carlo on Bravais lattices.

Used by Vincent for: 

- Camembert simulations

## Future

* Next, key:  Add a description of how the geometry works.
* Next: Change the geometry class to accept user input
* When time: perform small refactor of the particles.update library, which is particularly
  ugly and reuses code

### General

* Documentation (priority)
* Write an interface for running jobs
* Implement calculations of average energy, heat capacity, and correlation
  functions in `averages_utils`

### `default` model (`default_model_dev` branch)

* Write the `default_averages` functions to collect average occupation numbers

### `fields` model on branch `fields_dev`

* Debug Segmentation fault error in the fields library

# `frusa_mc`

## Dependencies

### Venv

The python part of the code is developed as a package.
To use it for the first time, go through the following steps in the `python` folder:

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

Every time you want to use the Python code, activate the virtual environment first with
`source .venv/bin/activate`


### JSON

The relevant header for `nlohmann::json` c++ library is located
in `include/thirdparty/`.
The python `json` library is a built-in module.

### CLI11

The relevant header for `CLI11` c++ library is located
in `include/thirdparty/`.
Used for parsing command-line arguments, mostly to choose input files.
(Source)[https://github.com/CLIUtils/CLI11]


## Building using cmake

To build `frusa_mc` executable, simply use

```
cmake <PATH_TO_PROJECT> -DMODEL_TYPE=<MODEL_NAME>
cmake --build .
```
If `-DMODEL_TYPE` flag is omitted, the `default` library is used.

### Existing `model` classes

Currently, the available `<MODEL_NAME>` options are:

* `default` - two-state system (__stable__)
* `fields_chain` - concentration fields on a 1D lattice
* `fields_square` - concentration fields on a 2d square lattice
* `fields_hexagonal` - concentration fields on a 2d hexagonal lattice
* particles` - discrete particles on a lattice of your choice

## Using the `lattice_particles` model

The `particles` model consists of 2 portions: the C++ part which does the
computational heavy lifting, and the Python part which can be used to generate the contact
maps in various geometries and visualise the simulation results.

In principle, going to the `python/examples` folder and looking at the notebooks should give
you enough informations to run simulations.

Note that, with the lattice_particles model, the frusa_mc executable takes 2 arguments, the
help for which you can see by calling `frusa_mc -h`:

```
Options:
  -h,--help                   Print this help message and exit
  -m,--model-params TEXT      JSON input file for model parameters
  -M,--mc-params TEXT         JSON input file for Monte-Carlo annealing parameters
```

## Generating contact map and plotting using the python utilities

The `python` folder contains useful codes and utilities for generating contact maps, and
handling the particle geometry in a more intuitive manner than the C++ code.

To view how to generate contact maps, run simulations, and view results in 2 and 3D, the best
is to look at the notebooks in `python/examples`.

One particulartly important file is `python/scripts/generate_cubic_permutations.py`, which is
used for generating the permutations of particle orientations when viewed under different
bonds for use in the C++ code.

### Creating a custom `model` class

* [Tutorial: writing your own model]
* [Default model README](src/models/default/README.md)
* [How to connect your library to the `model` class](src/models/README.md)
