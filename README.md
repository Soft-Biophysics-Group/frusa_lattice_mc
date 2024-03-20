# `frusa_mc` - Generic Monte Carlo module for mean-field and particle models of frustrated self-assembly

## Future

* Debug Segmentation fault error in the fields library (priority, WIP on branch `fields_dev`)
* Write a test model ising (priority, WIP on branch `ising_dev`)
* Documentation (priority)
* Write an interface for running jobs

## `frusa_mc`

### Dependencies

* __JSON__ - the relevant header for `nlohmann::json` c++ library is located\
int `include/thirdparty/`.

### Building using `cmake`

To build `frusa_mc` executable, simply use

```
cmake <PATH_TO_PROJECT> -DMODEL_TYPE=<MODEL_NAME>
cmake --build .
```

### Existing `model` classes

Currently, the available `<MODEL_NAME>` options are:

* `fields_chain` - concentration fields on a 1D lattice
* `fields_square` - concentration fields on a 2d square lattice 
* `fields_hexagonal` - concentration fields on a 2d hexagonal lattice

__WIP:__ currently, simulations of the fields models produce occasional\
Segmentation fault error. The debugging is done on branch `fields_dev`.

### Creating a custom `model` class
