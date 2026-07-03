# exciton_2b
This program (written by Denis Golež) uses the Keldysh formalism to evaluate the response of an excitonic insulator after a photoexcitation. The motivation for this comes from trARPES measurements in Ta2NiSe5.

The model describes two-orbital spinless electrons in a one-dimensional lattice, coupled to phonons. The coupling to the external EM field is included via the Peierls substitution. The program solves the Dyson equation on the Kadanoff-Baym (L-shaped) contour using the NESSI library. This program uses paralellization over momentum points using MPI.

## Structure of programs
The main program(s) is (are) in `program`. In `step` there is the most important part - the stepping procedure for solving the Dyson equation, whose solution are the electron and the phonon Green's function. Here the construction of the self-energy takes places, being split in different contributions. Here these are the Hartree, Fock, 2nd Born, and Migdal for the electronic self-energy. Before starting the stepping protocol, the hybridization based on the density of states in the bath and its coupling to the system. In `lattice` the lattice is defined, including the free-particle Hamiltonian, the interaction vertex, the current vertex, electron-phonon coupling, and also the conversion from electric field to vector potential. The NESSI library libcntr is imported. `parameters.hpp` and `print.hpp` within `program` are helpers for defining the types of parameters and for generating ouput files.

Many files are not needed for this purpose. I put them to `other`.

## Meaning of parameters

-__nt: number of time steps
-__ntau: number of steps on imaginary branch
__beta: inverse temperature
__h:
__mu: chemical potential 
__den: occupation


__itermax:
__errmax: 
__iter_rtime:
__err_rtime:
__walltime_limit
__kt: probably integration order used in solving Dyson's equation, but I need to check
__print_k:

__suscep
__omp_for_vie2
__outputfrequency
__restart_time

__test
__update

__fieldD
__fieldP
__eta

__nk: number of k points
__v01: Seed for pair breaking field (present in hamiltonian in first 5 iterations for equilibrium)
__v01_time
__mazza
__mix
__phonontype

__U: U, V and xi together define the interaction vertex as given in PHYSICAL REVIEW B 94, 035121 (2016); see `exciton_photo.hpp`
__V
__xi

__g: electron-phonon coupling; see `exciton_photo.hpp`
__omega0
__tt
__delta
__E
__dE
__dipolRe
__dipolIm
__gamma
__epsilon
__om0
__migdal
__s
__amp
__ratio


__bath_low
__bath_high
__gC_bath
__gV_bath
__gBATH
__omegaBATH


Redundant parameters:
__use_rpa: use_rpa is set to true
__v02: does not appear anywhere
__a_selective: does not appear anywhere
__approximation: removed it, I will use fixed approximation
