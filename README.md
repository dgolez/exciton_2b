# exciton_2b
This program (written by Denis Golež) uses the Keldysh formalism to evaluate the response of an excitonic insulator after a photoexcitation. The motivation for this comes from trARPES measurements in Ta2NiSe5.

The model describes two-orbital spinless electrons in a one-dimensional lattice, coupled to phonons. The coupling to the external EM field is included via the Peierls substitution. The program solves the Dyson equation on the Kadanoff-Baym (L-shaped) contour using the NESSI library. This program uses paralellization over momentum points using MPI.

## Structure of programs
The main program(s) is (are) in `program`. In `step` there is the most important part - the stepping procedure for solving the Dyson equation, whose solution are the electron and the phonon Green's function. Here the construction of the self-energy takes places, being split in different contributions. Here these are the Hartree, Fock, 2nd Born, and Migdal for the electronic self-energy. Before starting the stepping protocol, the hybridization based on the density of states in the bath and its coupling to the system. In `lattice` the lattice is defined, including the free-particle Hamiltonian, the interaction vertex, the current vertex, electron-phonon coupling, and also the conversion from electric field to vector potential. The NESSI library libcntr is imported. `parameters.hpp` and `print.hpp` within `program` are helpers for defining the types of parameters and for generating ouput files.

