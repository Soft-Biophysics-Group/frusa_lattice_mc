# `frusa_mc` - Generic Monte Carlo module for mean-field and particle models of frustrated self-assembly

## Scope of this code

Simulated-annealing Metropolis Monte Carlo on Bravais lattices.

Used by Vincent for: 

- Camembert simulations

## Important: compiling the code for the cluster

The cluster uses an x86_64 linux architecture, but the version of the C++ compiler on it does
not support some of the C++20 features used in this code.

To overcome this, I have compiled it on my M3 macbook, which involved setting up a
cross-compilation toolchain in `toolchain-x86_64-linux.cmake`.

Before you use it on macOS, make sure to install the x86 cross-compiler by running the
MacPorts command:
`sudo port install x86_64-elf-gcc x86_64-elf-binutils `.

Then, add the installed files to your PATH variable:
`export PATH=/opt/local/bin:$PATH` (has to be done for every shell session, or added to your
`.bashrc`/`.zshrc` file

Finally, run the commands:
`cmake -DMODEL_TYPE="lattice_particles" -DCMAKE_TOOLCHAIN_FILE=toolchain-x86_64-linux.cmake -B build`



## Notes during code annotation

* Only chain and triangular lattices are so far implemented!!!

## Future

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

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
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
* `lattice_particles` - discrete particles on a lattice of your choice

## Using the `lattice_particles` model

The `lattice_particles` model consists of 2 portions: the C++ part which does the
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
