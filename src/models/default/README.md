# Default model library -- two state system

This library is used as a default option in case no `<MODEL_TYPE>` flag is provided.

## Description of the system The system consists of $N$ particles, which may be in one of two
states with corresponding energies $\pm\Delta/2$. The parameters of the model are

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

$$\frac{C_v}{N} =
\frac{\beta^2\Delta^2}{4}\mathrm{sech}^2\left(\frac{\beta\Delta}{2}\right).$$

## File structure The core of the default (as well as any other user-defined) model library
consists of five basic function groups:

### `default_parameters`

This function defines a structure to store all user-specified, __model-specific__ parameters.

### `default_state`

Here, we store generic functions that create and manipulate the state structure (c++
`struct`) - an object that contains all configurational information about the system. In the
case of the `default` library, the state structure contains the size of the system, an array
containing state occupation (i.e. whether a particle is in state 0 or 1), and the average
occupation of the current state of the system.

### `default_interactions`

Similarly to the `default_state`, this file contains functions that create and manipulate the
interaction properties of the system, which are contained in a c++`struct` object. The
present structure contains the size of the gap $\Delta$ and the energy of the current state
of the system.

### `default_update`

This file contains a function used by MC engine to update the state of the system with
Metropolis acceptance probability.

### `default_average`

__This file is currently empty!__ Reserved for functions that calculate __model-specific__
averages. Generic average quantities, such as average energy and heat capacity should be
implemented in the `utils` library as `utils_averages`. In the case of the two-state system,
an example (to be implemented) of a model-specific average quantity is average occupation
number.

## Notes on `CMakeLists.txt`

The dependencies and paths to the header files are specified in the cmake file. In general,
the structure of the `CMakeLists.txt` will be the same for all user-defined model libraries,
but there is enough flexibility to allow for very fancy compile processes (e.g. see
`CMakeLists.txt` in the `fields` module).
