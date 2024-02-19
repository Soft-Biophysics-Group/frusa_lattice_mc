# `frusa_mc` - Generic Monte Carlo module for mean-field and particle models of frustrated self-assembly

## Future

* Define separate model libraries with cmake selection option (replace\
  template `mc_routines` by a class with `model` argument)
* Implement 2D square field model
* Implement 2D hexagonal field model
* Set PATH variables at compile time
* Implement useful averages (energy, heat capacity, correlation function,\
  cluster size) 

## `frusa_mc`

### Dependencies

* __JSON__

### Building using `cmake`

To build `frusa_mc` executable, simply use

```
cmake <PATH_TO_PROJECT> 
cmake --build .
```

### Existing `model` classes

### Creating a custom `model` class
