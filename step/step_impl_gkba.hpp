#pragma once

#include "step.hpp"


template <class LATTICE,class PHONON>
mpi_lattice_step_GKBA<LATTICE,PHONON>::mpi_lattice_step_GKBA(int nt,int ntau,int size,double beta,double h,bool use_omp_for_vie2,double mu,double epsilon,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &g,CFUNC &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,int nk,int kt):
	base_type::mpi_lattice_step(nt,ntau,size,beta,h,use_omp_for_vie2,mu,epsilon,tt,U,V,g,omega0,xi,delta,A,dA,v01,nk,kt){
}
