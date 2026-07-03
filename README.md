# exciton_2b
This program (written by Denis Golež) uses the Keldysh formalism to evaluate the response of an excitonic insulator after a photoexcitation. The motivation for this comes from trARPES measurements in Ta2NiSe5.

The model describes two-orbital spinless electrons in a one-dimensional lattice, coupled to phonons. The coupling to the external EM field is included via the Peierls substitution. The program solves the Dyson equation on the Kadanoff-Baym (L-shaped) contour using the NESSI library. This program uses paralellization over momentum points using MPI.

## Structure of programs
The main program(s) is (are) in `program`. In `step` there is the most important part - the stepping procedure for solving the Dyson equation, whose solution are the electron and the phonon Green's function. Here the construction of the self-energy takes places, being split in different contributions. Here these are the Hartree, Fock, 2nd Born, and Migdal for the electronic self-energy. Before starting the stepping protocol, the hybridization based on the density of states in the bath and its coupling to the system. In `lattice` the lattice is defined, including the free-particle Hamiltonian, the interaction vertex, the current vertex, electron-phonon coupling, and also the conversion from electric field to vector potential. The NESSI library libcntr is imported. `parameters.hpp` and `print.hpp` within `program` are helpers for defining the types of parameters and for generating ouput files.

Many files are not needed for this purpose. I put them to `other`.

## Meaning of parameters

- nt: number of time steps
- ntau: number of steps on imaginary branch
- beta: inverse temperature
- h:
- mu: chemical potential 
- den: occupation


- itermax:
- errmax: 
- iter_rtime:
- err_rtime:
- walltime_limit
- kt: probably integration order used in solving Dyson's equation, but I need to check
- print_k:

- suscep
- omp_for_vie2
- outputfrequency
- restart_time

- test
- update

- fieldD
- fieldP
- eta

- nk: number of k points
- v01: Seed for pair breaking field (present in hamiltonian in first 5 iterations for equilibrium)
- v01_time
- mazza
- mix
- phonontype

- U: U, V and xi together define the interaction vertex as given in PHYSICAL REVIEW B 94, 035121 (2016); see `exciton_photo.hpp`. in practice here V is put to 0 and in this case the value of xi doesn't matter
- V
- xi

- g: electron-phonon coupling; see `exciton_photo.hpp`
- omega0
- tt
- delta
- E
- dE
- dipolRe
- dipolIm
- gamma
- epsilon
- om0
- migdal
- s
- amp
- ratio


- bath_low
- bath_high
- gC_bath
- gV_bath
- gBATH
- omegaBATH


Redundant parameters:
- use_rpa: use_rpa is set to true
- v02: does not appear anywhere
- a_selective: does not appear anywhere
- approximation: removed it, I will use fixed approximation
