to do:

- finish documentation of various parameters below
- write section on how to run program

# exciton_2b
This program (written by Denis Golež) uses the non-equilibrium Green function formalism to evaluate the response of an excitonic insulator after a photoexcitation. The motivation for this comes from trARPES measurements in Ta2NiSe5.

The model describes two-orbital spinless electrons in a one-dimensional lattice, coupled to phonons. The coupling to the external EM field is included via the Peierls substitution. The program solves the Dyson equation on the Kadanoff-Baym (L-shaped) contour using the NESSI library. This program uses paralellization over momentum points using MPI.

## Structure of programs
The main program(s) is (are) in `program`. In `step` there is the most important part - the stepping procedure for solving the Dyson equation, whose solution are the electron and the phonon Green's function. Here the construction of the self-energy takes places, being split in different contributions. Here these are the Hartree, Fock, 2nd Born, and Migdal for the electronic self-energy. Before starting the stepping protocol, the hybridization based on the density of states in the bath and its coupling to the system. In `lattice` the lattice is defined, including the free-particle Hamiltonian, the interaction vertex, the current vertex, electron-phonon coupling, and also the conversion from electric field to vector potential. The NESSI library libcntr is imported. `parameters.hpp` and `print.hpp` within `program` are helpers for defining the types of parameters and for generating ouput files.

Many files are not needed for this purpose. I put them to `other`.

## Meaning of parameters

- nt: number of time steps
- ntau: number of steps on imaginary branch
- beta: inverse temperature
- h: time step
- mu: chemical potential 
- den: occupation

- itermax: max iterations of self-consistency loop for equilibrium G
- errmax: tolerance in self-consistency loop for equilibrium G
- iter_rtime: max iterations of self-consistency loop for non-equilibrium G
- err_rtime: tolerance in self-consistency loop for non-equilibrium G
- walltime_limit: sets the maximum elapsed real-world runtime, in seconds; the timer starts just after MPI initialization
- kt: probably integration order used in solving Dyson's equation, but I need to check
- print_k: if print_k=1, k-resolved Green's function is saved for every k. it saves ret, les, tv
- outputfrequency: time frequency of saving in case print_k=1; equilibrium is always included (timestep 0)

- omp_for_vie2: omp stands for OpenMP (library for parallel programming). but it is put to 0?
- restart_time: currently program works only for restart_time=-1, meaning start is with search for equilibrium solution. restart_time!=-1 is meant for the case when you give input state, however, for this case code has not been implemented fully

- test: if true (1), then for steps<=10 values of propagators are saved in all steps of obtaining self-consistency equation
- update: if true (1), then Hartree and Fock are updated in each time step, calculated from rho(t). if false, equilibrium values are used instead 

- nk: number of k points
- v01: Seed for pair breaking field (present in hamiltonian in first 5 iterations for equilibrium)
- v01_time: amplitude of Gaussian for pair breaking field in time simulation. this should be non-zero to obtain excitonic susceptibility
- mix: blends the previous iteration with the newly calculated result; for G and density matrix rho
- phonontype: takes values 0,1,or 3. means that the phonon mode couples to \Psi^\dag \sigma_{phonontype} \Psi. for excitonic case, we take phonontype=1

- mazza: interorbital hopping. beside dipole term, mazza determines off-diagonal elements of hk. its value is put to 0.116 (rescaling is done then by 0.3, which is Ni hopping)
- delta: orbital energies (diagonal elements of hk) are shifted by \pm delta/2 (and also equilibrium Hartree shift is subtracted!!). delta either set to constant or a time sequence of delta(t) is supplied

these parameters govern coupling to electromagnetic field
- E: either fixed (put to 0.0 for equilibrium) or give file with time series E(t); from that A(t) is calculated (see `exciton_2b.hpp`)
- dipolRe: real part of dipole moment, which multiplies E(t) to determine electromagnetic coupling. set to 1.0
- dipolIm: imaginary part -||-. set to 0.0
- fieldD: multiplicative factor in dipole term; multiplies E(t). set to 1.0
- fieldP: multiplicative factor in Peierls substitution; multiplies A(t) in k-shift. set to 1.0

- U, V and xi: they define the interaction vertex as given in PHYSICAL REVIEW B 94, 035121 (2016) (`exciton_2b.hpp`). in practice here V is put to 0 and in this case the value of xi doesn't matter, so that U is the only meaningful parameter here

- g: electron-phonon coupling (see `exciton_2b.hpp`)
- omega0: frequency of phonon mode; only works for omega0 a constant, do not suppy time series because code was not implemented for this
- migdal: if false (0), then the dynamic part of Sigma_ph ~ g D_0 g G, where g is el-ph coupling and D_0 is bare phonon propagator. if true (1), then dressed phonon propagator is solved via D=D_0 + D_0 Pi_ph D, where Pi_ph ~ gGGg. then, Sigma_ph ~ gDgG. migdal=true provides inclusions of fluctuations beyond classical el-ph coupling

- gBATH: coupling of electron to bath phonon
- omegaBATH: frequency of bath phonon mode

parameters not used in current simulation but kept for possible future use
- eta: artificial broadening of single-particle spectral functions used e.g. in calculating Seebeck
- tt: time series of hopping. not used in current simulation, but kept for possible future use

Redundant parameters, so I removed them:
- use_rpa: use_rpa was set to true
- v02: does not appear anywhere
- a_selective: does not appear anywhere
- approximation: removed it, I will use fixed approximation
- susucep: does not appear anywhere
- dE: not used anywhere
- gamma: relaxation rate in von Neumann equation, used only if mean field approximation is done. not needed here
- epsilon: seed for U_k e.g. to measure susceptibility if mean field approximation is done. not needed here
- om0, s, emp: they determine the nonthermal occupation perturbation if mean field approximation is done. not needed here
- ratio: used in vkFREE, which is not needed her
- bath_low, bath_high, gC_bath, gV_bath: they determine fermionic bath, which is however not used here (instead, a phononic bath is used)
  
## How one runs the program
build the executable with ./build.sh
run program by giving input file and output prefix/suffix: e.g. `.....cbuild/program/exciton_2b.ex input.txt out` 

## Possible TODO's i.e. ways to add functionalities
- make omega0 as time series work
- make tt as time series work (hkfree in `exciton_2b.hpp` needs to be fixed)
- transport coefficients L11,L12 and Seebeck are not properly calculated (fix get_seebeck in `step_impl_optical.hpp`)
- fix how get_seebeck and get_optical0 use vk (vkFREE no longer exists) if ever to be called
- include fermionic bath, add possibility that either phononic or fermionic bath is used. then parameters bath_low, bath_high, gC_bath, gV_bath would have to be included again
