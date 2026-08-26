#pragma once

#include <cmath>
#include <cassert>
#include <iostream>
#include <complex>
#include <vector>
#include "cntr/cntr.hpp"
#include "cntr/hdf5/hdf5_interface.hpp"
#include "cntr/hdf5/hdf5_interface_cntr.hpp"
//#include <boost/numeric/odeint.hpp>

#define ASSERT_0 1
#define ASSERT_1 1
#define FLEX_CAN_USE_OMP 1
#define FLEX_USE_MPI 1

#define CPLX std::complex<double>
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>
#define GREEN_TVIEW cntr::herm_matrix_timestep_view<double>
#define CFUNC cntr::function<double>



template <class LATTICE> class kpoint_density{
public:
	kpoint_density(void);
	kpoint_density(int nt,int ntau,int size,double beta,double h,double kk,LATTICE &latt,double mu,double den,double mix=0.0,int phonontype=0);
	//////////////////////////////////////////////////////////////////////////////////////////
	// OMP parallelization: (parellel solution of RPA integral equations)
	// must be switched on explicitly !!
	
	void init_rho_free(LATTICE &latt);
	void use_omp(bool onoff);
	void init_G_mat_nointeraction(LATTICE &latt);
	void set_hk(int tstp,int iter,LATTICE &latt);
	void set_vertex(int tstp,LATTICE &latt);
	void step_dyson(int tstp,int iter,LATTICE &latt);
	double step_dyson_with_error(int tstp,int iter,LATTICE &latt);
	double heaviside(double ene);
	void write_to_hdf5(hid_t group_id);
	void write_to_hdf5(const char *filename);
	void write_to_hdf5_density(hid_t group_id,int dt,int tid);
	double beta_;
	double h_;
	int nt_;
	int ntau_;
	int nrpa_; // Dimension of the Green's function
	CFUNC rho_,rho_diag_,SHartree_,SFock_,hk_,hkeff_,vertex_;
	CFUNC eigen_vec_,hkeff_eigen_;
	CFUNC rho_loc_;
	cdmatrix rho_eq_;
	double mu_,den_;
	double kk_;
	bool use_omp_;
	double mix_;
	int phonontype_;
};


template <class LATTICE> class kpoint_green{
public:
	double beta_;
	double h_;
	int nt_;
	int ntau_;
	int nrpa_;
	int migdal_;
	double mu_;
	double kk_;
	double mix_;
	bool use_omp_;
	CFUNC gph_;
	GREEN G_,G0_,G0Sigma_,SigmaG0_,Sigma_,chi_,P_,VP_,PV_,W_,D0_,D_,Pph_,D0Pph_,PphD0_;
	
	kpoint_green(void);
    kpoint_green(int nt,int ntau,int size,double beta,double h,double kk,LATTICE &latt,double mu,CFUNC &omega,CFUNC &g,double mix,int migdal=0);
	void use_omp(bool onoff);
	void get_Density_matrix(int tstp,kpoint_density<LATTICE> &density);
	void step_chi(int tstp,int kt,LATTICE &latt,kpoint_density<LATTICE> & density);
	void step_W(int tstp,int kt,LATTICE &latt,kpoint_density<LATTICE> & density);
	void step_W2b(int tstp,int kt,LATTICE &latt,kpoint_density<LATTICE> & density);
	void step_D(int tstp,int kt,LATTICE &latt,kpoint_density<LATTICE> & density);
	void step_dyson(int tstp,int iter,int kt,LATTICE &latt,kpoint_density<LATTICE> &density);
	double step_dyson_with_error(int tstp,int iter,int kt,LATTICE &latt,kpoint_density<LATTICE> &density);
	void step_dyson_integral(int tstp,int iter,int kt,LATTICE &latt,kpoint_density<LATTICE> &density);
	double step_dyson_with_error_integral(int tstp,int iter,int kt,LATTICE &latt,kpoint_density<LATTICE> &density);
	void init_G_mat_nointeraction(LATTICE &latt,kpoint_density<LATTICE> & density,int kt);
	// Reading/printing
	// void read_from_file(int nt1,const char *folder,const char *suffix,LATTICE &latt);
	void read_from_hdf5(int tstp,hid_t group_id,LATTICE &latt,kpoint_density<LATTICE> &density);
	void read_from_hdf5(int nt1,const char *filename,LATTICE &latt,kpoint_density<LATTICE> &density);
	// void print_to_file(const char *folder,const char *suffix);
	void write_to_hdf5(hid_t group_id,kpoint_density<LATTICE> &density);
	void write_to_hdf5(const char *filename,kpoint_density<LATTICE> &density);
	void write_to_hdf5_slices(hid_t group_id,int dt,int tid);
	void write_to_hdf5_slices(const char *filename,int dt,int tid);
	void append_to_hdf5(hid_t group_id,int tstp,bool write_parameters);
};
