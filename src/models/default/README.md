## Default model library -- two state system 

This library is used as a default option in case no `<MODEL_TYPE>` flag is\
provided.

### Description of the system
The system consists of $N$ particles, which may be in one of two states with\
corresponding energies $\pm\Delta/2$.
The parameters of the model are

* $N$ - number of particles
* $\Delta$ - energy gap
* $T$ - temperature

The model is exactly solvable with partition function equal to

$$Z = 2^N\cosh^N\left(\frac{\beta\Delta}{2}\right).$$

From this, we can calculate the free energy

$$-\frac{1}{N}\beta F = \ln{2} +\ln\left[\cosh\left(\frac{\beta\Delta}{2}\right)\right],

average energy

$$U = \frac{Delta}{2}\tanh\left(\frac{\beta\Delta}{2}\right)$$
