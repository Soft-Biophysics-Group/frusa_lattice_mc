# `frusa_mc` - Generic Monte Carlo module for mean-field and particle models of frustrated self-assembly

## Future

### General

* Documentation (priority)
* Write an interface for running jobs
* Implement calculations of average energy, heat capacity, and correlation
  functions in `averages_utils`

### `default` model (`default_model_dev` branch)

* Write the `default_averages` functions to collect average occupation numbers

### Tutorial and `ising` model (`ising_dev` branch)

* Write a tutorial on ising model (priority)

### `fields` model on branch `fields_dev`

* Debug Segmentation fault error in the fields library (priority)

### `lattice_particles` model on branch `lattice_mc_dev`

- Someday: make a more elegant iteration solution to loop through the sites and occupation
  arrays of the state struct
- Maybe: change the "initialize_state" functions into full-blown constructors?
- Next: check the hashing give correct values, and that I don't need to subtract 1 from the
  n_states and n_orientations.
- Someday: make a lattice library where we compile the functions (e.g. geometry) used by both
  the fields and lattice particles libraries.
- Next: write a Python code to generate the flattened LEL from 21 interaction energies for
  hexagonal particles
- Next: Figure out where to create the array of neighbours used in geometry
- When I implement the moves: make sure I generate only hashed integers, rather than pairs or
  triplets, in order to keep things faaaaast (i.e. generate a full state rather than an
  orientation + a type)
- Someday/maybe: Pre-calculate the exponentials of the interactions so that we can calculate
  the energy differences faster during MC moves
- Next: the interactions structure is redundant for lattice particles and can be absorbed in
  the state struct (I only need to keep track  of the energy)
- Now: partial reformatting
  - Site should be a class, not a struct. Should store state as well as orientation and type
    in order to be more consistent. Should be updated using external methods
  - When I get the contacts between the Sites, 

# `frusa_mc`

## Dependencies

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
