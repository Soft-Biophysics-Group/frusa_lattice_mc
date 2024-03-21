# `frusa_mc` - Generic Monte Carlo module for mean-field and particle models of frustrated self-assembly

## Future

### General

* Documentation (priority)
* Write an interface for running jobs
* Implement calculations of average energy, heat capacity, and correlation 
  functios in `averages_utils`

### `default` model (`default_model_dev` branch)

* Write the `default_averages` functions to collect average occupation numbers

### Tutorial and `ising` model (`ising_dev` branch)

* Write a tutorial on ising model (priority)

### `fields` model on branch `fields_dev`

* Debug Segmentation fault error in the fields library (priority)

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
