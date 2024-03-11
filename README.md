# `frusa_mc` - Generic Monte Carlo module for mean-field and particle models of frustrated self-assembly

## Future

* Currently, donor and acceptor lists are updated incorrectly. There is no 
  feature to check if site is already in the list and, as a result, the 
  simulation gives Segmentation fault
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
