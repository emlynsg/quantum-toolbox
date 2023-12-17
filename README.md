QuantumToolbox
==============

The QuantumToolbox is a library for time-dependent calculations in quantum mechanics, developed by Emlyn Graham at the Department of Nuclear Physics at the Australian National University. The toolbox is structured around several key components, each represented by a class:

* Grid: This class represents the grid on which calculations are performed. The grid is defined by its size, minimum and maximum values, and a scale factor.

* Wavefunction: This class represents the wave function of the quantum system. It can be initialized using several functions (including Gaussian or constant functions), and its energy can be boosted for directional propagation.

* Potential: This class represents the potential in the quantum system. It can be initialized using several common nuclear physics potentials (including Coulomb and Woods-Saxon).

* System: This class represents the quantum system as a whole, consisting of a wave function and a potential. The system can evolve over time, either as a single internal state or as multiple coupled states (so-called coupled channels).

* Plotter: This class is used to visualize the quantum system. It can animate the evolution of the system over time.

The toolbox also includes scripts for running specific tests or simulations, such as tunneling through a square barrier or a Gaussian barrier with an absorbing potential. These scripts import the QuantumToolbox library, set up a grid, wave function, and potential, create a system, and then use a plotter to visualize the system's evolution over time.
