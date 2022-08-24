QuantumToolbox
==============

A toolbox for time dependent calculations in quantum mechanics.

Developed by Emlyn Graham at the Department of Nuclear Physics at the Australian National University.


Build Notes - Ed Simpson
------------------------

For Ubuntu 22.04.

Packages:

sudo apt install cmake libgsl-dev libboost-dev libboost-filesystem-dev libboost-system-dev libboost-regex-dev libboost-iostreams-dev libeigen3-dev

To do - repository and git tree
-------------------------------

@ Remove gnuplot-iostream from the tree and use git version. With instructions for installing.

@ Identify what madplotlib is used for, if anything. If not, remove from tree.

@ Create generic cmake build script for projects (i.e. scripts).

@ QTPlotter appears to be no longer used - superseded by gnuplot. Remove from tree.

@ Fix boost header errors


To do - physics
---------------

@ Implement expansion of state in Coulomb wave functions.

@ Come up with some way of looping/storage.

@ Investigate and check what is in Extras.cpp, particularly the
integration routines.