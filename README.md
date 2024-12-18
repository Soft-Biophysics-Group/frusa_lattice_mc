# `frusa_mc` - Generic Monte Carlo module for mean-field and particle models of frustrated self-assembly

## Future

* Next: Change the geometry class to accept user input
* Next: Implement Python classes to generate the interaction matrix
* Next: Implement 1D lattice geometry and test the simulation on it first
* When time: perform small refactor of the particles.update library, which is particularly
  ugly and reuses code

### General

* Documentation (priority)
* Geometry class (priority)
* Write an interface for running jobs
* Implement calculations of average energy, heat capacity, and correlation
  functions in `averages_utils`

### `default` model (`default_model_dev` branch)

* Write the `default_averages` functions to collect average occupation numbers

### Testing

* V: I have for now added a `test/test_geometry.cc` file, the contents of which can be
  copy-pasted into `app/main.cc` in order to test the Geometry class we have written.

### `fields` model on branch `fields_dev`

* Debug Segmentation fault error in the fields library

### `lattice_particles` model on branch `lattice_mc_dev`

- Someday: make a lattice library where we compile the functions (e.g. geometry) used by both
  the fields and lattice particles libraries.
- Someday/maybe: Pre-calculate the exponentials of the interactions so that we can calculate
  the energy differences faster during MC moves
- Next: figure out where to declare the move_probas array in the code
- When I start the cubic lattice: start specializing at compile time like Andrey does
- After I'm done: regularize the number of edges to be a constexpr int defined at compile
  time, depending on the compilation option chosen

# `frusa_mc`

## Dependencies

### Venv

The python part of the code is developed as a package.
To use it for the first time, go through the following steps in the `python` folder:

```bash
python -m venv .venv
source .venv/bin/activate
pip install requirements.txt
```

Every time you want to use the Python code, activate the virtual environment first with
`source .venv/bin/activate`


### JSON

The relevant header for `nlohmann::json` c++ library is located
in `include/thirdparty/`.
The python `json` library is a built-in module.


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

### Creating a custom `model` class

* [Tutorial: writing your own model]
* [Default model README](src/models/default/README.md)
* [How to connect your library to the `model` class](src/models/README.md)
