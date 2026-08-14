#pragma once

#include <cmath>
#include <cassert>
#include <iostream>
#include <complex>
#include <vector>
#include "cntr/cntr.hpp"
#include "cntr/hdf5/hdf5_interface.hpp"
#include "cntr/hdf5/hdf5_interface_cntr.hpp"
#include "parameters.hpp"

#define ASSERT_0 1
#define ASSERT_1 1
#define FLEX_CAN_USE_OMP 1
#define FLEX_USE_MPI 1

#define CPLX std::complex<double>
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>
#define GREEN_TVIEW cntr::herm_matrix_timestep_view<double>
#define CFUNC cntr::function<double>

class phonon
{
  CFUNC omega0_;
  CFUNC g_;
  CFUNC density_; // Local single particle density matrix
  int phonontype_;
public:
  phonon(void);
  phonon(CFUNC &omega0,CFUNC &g,CFUNC &density,int phonontype=1);
  void eq(CFUNC &X,CFUNC &Pi);
  void step(int tstp,CFUNC &X,CFUNC &Pi,int kt,double h);
  double dx1dt(int tstp,CFUNC &Pi);
  double dx2dt(int tstp,CFUNC &X);
  void hartree(cdmatrix &inter,cdmatrix &X);
  void rpa(int tstp,GREEN_TSTP &SigmaPh,cntr::herm_matrix_timestep_view<double> & G,cntr::herm_matrix_timestep_view<double> & D);
  void matrixele(int i1,int i2,int &j1, int &j2,int nrpa);
};

template <class LATTICE> class approx{
	public:
	LATTICE latt_;
	phonon phonon_;
	int nk_;   // taken from lattice
	std::vector<kpoint_density<LATTICE> >  density_k_;
	std::vector<CFUNC> vertex_;
	cntr::distributed_array<double> convergence_error_;
	double beta_,epsilon_;
	double h_;
	std::vector<cdmatrix> rk_,sh_,sf_,s2rk_;
	CFUNC rho_loc_,rho_sym_,X_,Pi_,order_;
	std::vector<int> tid_map_;
	int ntasks_,tid_,tid_root_;
	int nt_;
	int ntau_;
	int norb_;
	int nrpa_;
	int kt_;
	// Phonons
	// boost::numeric::odeint::runge_kutta4< std::vector<double> > stepper_;

	virtual double step(int tstp,int iter,int kt,double om0,double s,double amp)=0;
	// virtual double density_matrix(int tstp,cdmatrix density)=0;
	virtual double get_ekin(int tstp)=0;
	virtual void   get_ekin(int tstp,cdmatrix &ekin)=0;
	virtual double get_curr(int tstp)=0;
	virtual double get_dAk(int tstp)=0;
	virtual double get_eneHF(int tstp)=0;
    //virtual double get_optical0_eq(double om,int mu,int nu)=0;
	virtual void read_from_file_hdf5(int tstp,const char *filename_prefix,int kt)=0;
	virtual void print_to_file_hdf5(const char *filename_prefix,int print_k)=0;
	virtual void print_to_file_hdf5_slice(const char *filename_prefix,int outputfrequency, int print_k)=0;
	virtual ~approx(void)=default;
};

template <class LATTICE> class mpi_lattice_step_optical: public approx<LATTICE>{
public:
	using approx<LATTICE>::latt_;
	using approx<LATTICE>::nk_;
	using approx<LATTICE>::density_k_;
	using approx<LATTICE>::vertex_;
	using approx<LATTICE>::convergence_error_;
	using approx<LATTICE>::beta_;
	using approx<LATTICE>::epsilon_;
	using approx<LATTICE>::h_;
	using approx<LATTICE>::rk_;
	using approx<LATTICE>::sh_;
	using approx<LATTICE>::sf_;
	using approx<LATTICE>::s2rk_;
	using approx<LATTICE>::rho_loc_;
	using approx<LATTICE>::rho_sym_;
	using approx<LATTICE>::order_;
	using approx<LATTICE>::X_;
	using approx<LATTICE>::Pi_;
	using approx<LATTICE>::tid_map_;
	using approx<LATTICE>::ntasks_;
	using approx<LATTICE>::tid_;
	using approx<LATTICE>::tid_root_;
	using approx<LATTICE>::nt_;
	using approx<LATTICE>::ntau_;
	using approx<LATTICE>::norb_;
	using approx<LATTICE>::nrpa_;
	// using approx<LATTICE>::stepper_;
	using approx<LATTICE>::kt_;
	

	bool update_;
	double eta_; //Broadening
	double mazza_;
	parameters param_;
	// int migdal_;
	/////////////////////////////////////////////////////////
	// before this, latt_ must be initialized !!
  	mpi_lattice_step_optical(parameters &param);
	// void init(int nt,int ntau,int size,double beta,double h,bool use_omp_for_vie2,double mu,double epsilon);
	void get_Sigma_Hartree_electronic(int tstp,cdmatrix &S);
	void get_Sigma_phonon_mean_field(int tstp,cdmatrix &S);
	void get_Sigma_Hartree(int tstp,int kk,CFUNC &S);
	void get_Sigma_Fock(int tstp,int kk,CFUNC &S);
	void set_local(int tstp);
	void set_sym(int tstp);
	void set_order(int tstp);
	void extrapolate_rho(int tstp);
	// void read_from_file(int nt1,const char *folder);
	// void read_from_file_hdf5(int tstp,const char *filename_prefix,int kt);
	///////////////////////////////////////////////////////////
	// double step(int tstp,int iter,int kt);
	void step_boson(int tstp);
        double get_dip(int tstp);
	double get_curr_dip(int tstp);
  double get_curr_peierls(int tstp);
  std::complex<double> get_optical0(CFUNC &optics,double domega,int nomega,double om0,double s,double amp);
  std::complex<double> get_chi0(double omega,int mu,int nu,double om0,double s,double amp);
  std::complex<double> get_seebeck(CFUNC &seebeck,CFUNC &omega,double domega,int nomega);
  double get_seebeck_boltzmann(double beta);
  double dfermi(double omega,double beta);
  std::complex<double> get_dos(CFUNC &dos,double domega,int nomega);
  double gauss(double om,double om0,double s);

	virtual double step(int tstp,int iter,int kt,double om0,double s,double amp) override;
	virtual void read_from_file_hdf5(int tstp,const char *filename_prefix,int kt) override;
	virtual void print_to_file_hdf5(const char *filename_prefix,int print_k);
	virtual void print_to_file_hdf5_slice(const char *filename_prefix,int outputfrequency,int print_k);
	virtual double get_ekin(int tstp) override;
	virtual void get_ekin(int tstp,cdmatrix &ekin) override;
	virtual double get_curr(int tstp) override;
	virtual double get_dAk(int tstp) override;
	virtual double get_eneHF(int tstp) override;
	
};

template <class LATTICE> class mpi_lattice_step_2b_optical: public mpi_lattice_step_optical<LATTICE>{
public:
	typedef mpi_lattice_step_optical<LATTICE> base_type;
	using approx<LATTICE>::latt_;
	using approx<LATTICE>::nk_;
	using approx<LATTICE>::vertex_;
	using approx<LATTICE>::convergence_error_;
	using approx<LATTICE>::beta_;
	using approx<LATTICE>::epsilon_;
	using approx<LATTICE>::h_;
	using approx<LATTICE>::rk_;
	using approx<LATTICE>::sh_;
	using approx<LATTICE>::sf_;
	using approx<LATTICE>::s2rk_;
	using approx<LATTICE>::X_;
	using approx<LATTICE>::Pi_;
	using approx<LATTICE>::tid_map_;
	using approx<LATTICE>::ntasks_;
	using approx<LATTICE>::tid_;
	using approx<LATTICE>::tid_root_;
	using approx<LATTICE>::nt_;
	using approx<LATTICE>::ntau_;
	using approx<LATTICE>::norb_;
	using approx<LATTICE>::nrpa_;
	// using approx<LATTICE>::stepper_;
	using approx<LATTICE>::kt_;

	std::vector<kpoint_green<LATTICE> > green_k_;
	cntr::distributed_timestep_array<double>  gk_all_timesteps_;  //Timestep of all k dependent Green's functions in array
	cntr::distributed_timestep_array<double>  wk_all_timesteps_;
	cntr::distributed_timestep_array<double>  dk_all_timesteps_;
	GREEN Gloc_,Wloc_,Dloc_,Bath_,D0loc_,D0hol_;
	CFUNC g_;
	/// kk-dependent observables, to be collected at rank tid_root_
	bool test_;
	bool use_omp_for_vie2_;
	int migdal_,phonontype_;
	/////////////////////////////////////////////////////////
	// before this, latt_ must be initialized !!
  	mpi_lattice_step_2b_optical(parameters &param);
	void init_G_mat_nointeraction(void);
	void set_rpa(bool onoff);
	void extrapolate_timestep_G(int tstp,int kt);
	void extrapolate_timestep_D(int tstp,int kt);
	void extrapolate_timestep_W(int tstp,int kt);
	void extrapolate_timestep_Sigma(int tstp,int kt);
	void extrapolate_timestep_Pi(int tstp,int kt);
	void symmetrise_mat(GREEN &G);
	void step_boson(int tstp);
	void gather_gk_timestep(int tstp);
	void set_density_k(int tstp);
	void gather_Wk_timestep(int tstp);
	void gather_Dk_timestep(int tstp);
	void get_PP_Bubble(int tstp,int qq,GREEN &P);
	void get_phonon_Bubble(int tstp,int qq,GREEN &P);
	void get_Sigma_2b(int tstp,int kk,GREEN &S,int iter);
	void get_Gloc(int tstp);
	void get_Wloc(int tstp);
	void get_Dloc(int tstp);
	void print_diagonal_occupations(const char *file_prefix,int t1,int t2,int kt,double h,int size);
	void gather_kk_observables(int tstp, int kt);
	double get_eneRPA(int tstp,int kt);

	virtual double step(int tstp,int iter,int kt,double om0,double s,double amp) override;
	virtual void read_from_file_hdf5(int tstp,const char *filename_prefix,int kt) override;
	virtual void print_to_file_hdf5(const char *filename_prefix,int print_k);
	virtual void print_to_file_hdf5_slice(const char *filename_prefix,int outputfrequency,int print_k);
};
