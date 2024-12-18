# Default model library -- two state system 

This library is used as a default option in case no `<MODEL_TYPE>` flag is 
provided.

## Description of the system
The system consists of $N$ particles, which may be in one of two states with 
corresponding energies $\pm\Delta/2$.
The parameters of the model are

* $N$ - number of particles
* $\Delta$ - energy gap
* $T$ - temperature

The model is exactly solvable with partition function equal to

$$Z = 2^N\cosh^N\left(\frac{\beta\Delta}{2}\right).$$

From this, we can calculate the free energy

$$-\frac{\beta F}{N} = \ln{2} +\ln\left[\cosh\left(\frac{\beta\Delta}{2}\right)\right],$$

average energy

$$\frac{U}{N} = -\frac{\Delta}{2}\tanh\left(\frac{\beta\Delta}{2}\right),$$

and heat capacity

$$\frac{C_v}{N} = \frac{\beta^2\Delta^2}{4}\mathrm{sech}^2\left(\frac{\beta\Delta}{2}\right).$$

If the two states are defined in terms of occupation numbers $n=0,1$, the 
average state has the form 

$$\frac{\langle n\rangle}{N} = \frac{1}{2}\left[1-\tanh\left(\frac{\beta\Delta}{2}\right)\right].$$

## File structure
The core of the default (as well as any other user-defined) model library 
consists of five basic function groups: 

### `default_parameters`

This function defines a structure to store all user-specified, 
__model-specific__ parameters.

### `default_state`

Here, we store generic functions that create and manipulate the state 
structure (c++ `struct`) - an object that contains all configurational 
information about the system. 
In the case of the `default` library, the state structure contains the size of 
the system, an array containing state occupation (i.e. whether a particle is in
state 0 or 1), and the average occupation of the current state of the system.

### `default_interactions`

Similarly to the `default_state`, this file contains functions that create and 
manipulate the interaction properties of the system, which are contained in a 
c++`struct` object.
The present structure contains the size of the gap $\Delta$ and the energy of 
the current state of the system.

### `default_update`

This file contains a function used by MC engine to update the state of the 
system with Metropolis acceptance probability.

### `default_average`

This file contains routines that manipulate the c++ structure containing 
thermodynamic averages collected during the simulation. 
These functions initialize, update, and save the averages to the respective
files.
In the case of the two-state system, the averages include the first and second
energy moments (the latter can be used to calculate the heat capacity), and the
average state (or average occupation number).
These two types of averages have corresponding boolean option variables, which
are passed in the parameters structure.

## Notes on `CMakeLists.txt`

The dependencies and paths to the header files are specified in the cmake 
file.
In general, the structure of the `CMakeLists.txt` will be the same for all 
user-defined model libraries, but there is enough flexibility to allow for
very fancy compile processes (e.g. see `CMakeLists.txt` in the `fields` module).
