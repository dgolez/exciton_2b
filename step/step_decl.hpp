#pragma once

#include <cmath>
#include <cassert>
#include <iostream>
#include <complex>
#include <vector>
#include "cntr/cntr.hpp"
#include "cntr/hdf5/hdf5_interface.hpp"
#include "cntr/hdf5/hdf5_interface_cntr.hpp"
#include "../program/parameters.hpp"
// #include <boost/numeric/odeint.hpp>

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
  void step(int tstp,CFUNC &X,CFUNC &Pi,int kt,double dt);
  double dx1dt(int tstp,CFUNC &Pi);
  double dx2dt(int tstp,CFUNC &X);
  void hartree(cdmatrix &inter,cdmatrix &X);
  void rpa(int tstp,GREEN_TSTP &SigmaPh,cntr::herm_matrix_timestep_view<double> & G,cntr::herm_matrix_timestep_view<double> & D);
  void matrixele(int i1,int i2,int &j1, int &j2,int nrpa);
};


// class phonon_hol
// {
// 	CFUNC omega0_;
// 	CFUNC g_;
// 	CFUNC density_; // Local single particle density matrix
// public:
//   phonon_hol(void);
//   phonon_hol(CFUNC &omega0,CFUNC &g,CFUNC &density);
//   void eq(CFUNC &X,CFUNC &Pi);
//   void step(int tstp,CFUNC &X,CFUNC &Pi,int kt,double dt);
//   double dx1dt(int tstp,CFUNC &Pi);
//   double dx2dt(int tstp,CFUNC &X);
//   void rpa(int tstp,int nrpa, GREEN_TSTP &SigmaPh,cntr::herm_matrix_timestep_view<double> & G,cntr::herm_matrix_timestep_view<double> & D);
//   void hartree(cdmatrix &inter,cdmatrix X);
// };


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
	virtual void print_to_file_hdf5_slice(const char *filename_prefix,int dt, int print_k)=0;
	virtual ~approx(void)=default;
};


template <class LATTICE> class mpi_lattice_step: public approx<LATTICE>{
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

	/////////////////////////////////////////////////////////
	// before this, latt_ must be initialized !!
  mpi_lattice_step(int nt,int ntau,int size,double beta,double h,bool use_omp_for_vie2,double mu,double epsilon,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &g,CFUNC &omega0,CFUNC &gBATH,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,int nk,int kt);
	// void init(int nt,int ntau,int size,double beta,double h,bool use_omp_for_vie2,double mu,double epsilon);
	void get_Sigma_Hartree(int tstp,int kk,CFUNC &S);
	void get_Sigma_Fock(int tstp,int kk,CFUNC &S);
	void set_local(int tstp);
	void set_sym(int tstp);
	void extrapolate_rho(int tstp);
	// void read_from_file(int nt1,const char *folder);
	// void read_from_file_hdf5(int tstp,const char *filename_prefix,int kt);
	///////////////////////////////////////////////////////////
	// double step(int tstp,int iter,int kt);
	void step_boson(int tstp);

	virtual double step(int tstp,int iter,int kt,double om0,double s,double amp) override;
	virtual void read_from_file_hdf5(int tstp,const char *filename_prefix,int kt) override;
	virtual void print_to_file_hdf5(const char *filename_prefix,int print_k);
	virtual void print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k);
	virtual double get_ekin(int tstp) override;
	virtual void get_ekin(int tstp,cdmatrix &ekin) override;
	virtual double get_curr(int tstp) override;
	virtual double get_dAk(int tstp) override;
	virtual double get_eneHF(int tstp) override;
	
};

template <class LATTICE> class mpi_lattice_deb_step: public approx<LATTICE>{
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

	/////////////////////////////////////////////////////////
	// before this, latt_ must be initialized !!
  mpi_lattice_deb_step(parameters &param);
	// void init(int nt,int ntau,int size,double beta,double h,bool use_omp_for_vie2,double mu,double epsilon);
	void get_Sigma_Hartree(int tstp,int kk,CFUNC &S);
	void get_Sigma_Fock(int tstp,int kk,CFUNC &S);
	void set_local(int tstp);
	void extrapolate_rho(int tstp);
	// void read_from_file(int nt1,const char *folder);
	// void read_from_file_hdf5(int tstp,const char *filename_prefix,int kt);
	///////////////////////////////////////////////////////////
	// double step(int tstp,int iter,int kt);
	void step_boson(int tstp);

	virtual double step(int tstp,int iter,int kt,double om0,double s,double amp) override;
	virtual void read_from_file_hdf5(int tstp,const char *filename_prefix,int kt) override;
	virtual void print_to_file_hdf5(const char *filename_prefix,int print_k);
	virtual void print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k);
	virtual double get_ekin(int tstp) override;
	virtual void get_ekin(int tstp,cdmatrix &ekin) override;
	virtual double get_curr(int tstp) override;
	virtual double get_dAk(int tstp) override;
	virtual double get_eneHF(int tstp) override;
	
};


template <class LATTICE> class mpi_lattice_step_cdw: public approx<LATTICE>{
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

	/////////////////////////////////////////////////////////
	// before this, latt_ must be initialized !!
  	mpi_lattice_step_cdw(int nt,int ntau,int size,double beta,double h,bool use_omp_for_vie2,double mu,double epsilon,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &g,CFUNC &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,int nk,int kt);
	// void init(int nt,int ntau,int size,double beta,double h,bool use_omp_for_vie2,double mu,double epsilon);
	void get_Sigma_Hartree(int tstp,int kk,CFUNC &S);
	void get_Sigma_Fock(int tstp,int kk,CFUNC &S);
	void set_local(int tstp);
	void extrapolate_rho(int tstp);
	// void read_from_file(int nt1,const char *folder);
	// void read_from_file_hdf5(int tstp,const char *filename_prefix,int kt);
	///////////////////////////////////////////////////////////
	// double step(int tstp,int iter,int kt);
	void step_boson(int tstp);

	virtual double step(int tstp,int iter,int kt,double om0,double s,double amp) override;
	virtual void read_from_file_hdf5(int tstp,const char *filename_prefix,int kt) override;
	virtual void print_to_file_hdf5(const char *filename_prefix,int print_k);
	virtual void print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k);
	virtual double get_ekin(int tstp) override;
	virtual void get_ekin(int tstp,cdmatrix &ekin) override;
	virtual double get_curr(int tstp) override;
	virtual double get_dAk(int tstp) override;
	virtual double get_eneHF(int tstp) override;
	
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
	int suscep_;
	double eta_; //Broadening
	double mazza_;
	parameters param_;
	// int migdal_;
	/////////////////////////////////////////////////////////
	// before this, latt_ must be initialized !!
  	mpi_lattice_step_optical(parameters &param);
	// void init(int nt,int ntau,int size,double beta,double h,bool use_omp_for_vie2,double mu,double epsilon);
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
	virtual void print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k);
	virtual double get_ekin(int tstp) override;
	virtual void get_ekin(int tstp,cdmatrix &ekin) override;
	virtual double get_curr(int tstp) override;
	virtual double get_dAk(int tstp) override;
	virtual double get_eneHF(int tstp) override;
	
};



template <class LATTICE> class mpi_lattice_step_RPA: public mpi_lattice_step<LATTICE>{
public:
	typedef mpi_lattice_step<LATTICE> base_type;
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
	GREEN Gloc_,Wloc_,Dloc_;
	/// kk-dependent observables, to be collected at rank tid_root_
	bool test_;
	bool use_omp_for_vie2_;
	/////////////////////////////////////////////////////////
	// before this, latt_ must be initialized !!
	mpi_lattice_step_RPA(int nt,int ntau,int size,double beta,double h,bool use_omp_for_vie2,double mu,double epsilon,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &g,CFUNC &omega0,CFUNC &gBATH,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,int nk,int kt,bool test,double mix);
	void init_G_mat_nointeraction(void);
	void set_rpa(bool onoff);
	void extrapolate_timestep_G(int tstp,int kt);
	void extrapolate_timestep_W(int tstp,int kt);
	void extrapolate_timestep_Sigma(int tstp,int kt);
	void extrapolate_timestep_Pi(int tstp,int kt);
	void symmetrise_mat(GREEN &G);
	void step_boson(int tstp);
	void gather_gk_timestep(int tstp);
	void set_density_k(int tstp);
	void gather_Wk_timestep(int tstp);
	void gather_Dk_timestep(int tstp);
	void get_Polarization_Bubble(int tstp,int qq,GREEN &P);
	void get_Sigma_RPA(int tstp,int kk,GREEN &S,int iter);
	void get_Gloc(int tstp);
	void get_Wloc(int tstp);
	void get_Dloc(int tstp);
	void print_diagonal_occupations(const char *file_prefix,int t1,int t2,int kt,double h,int size);
	void gather_kk_observables(int tstp, int kt);
	double get_eneRPA(int tstp,int kt);

	virtual double step(int tstp,int iter,int kt,double om0,double s,double amp) override;
	virtual void read_from_file_hdf5(int tstp,const char *filename_prefix,int kt) override;
	virtual void print_to_file_hdf5(const char *filename_prefix,int print_k);
	virtual void print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k);
};



template <class LATTICE> class mpi_lattice_step_RPA_cdw: public mpi_lattice_step_cdw<LATTICE>{
public:
	typedef mpi_lattice_step_cdw<LATTICE> base_type;
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
	GREEN Gloc_,Wloc_;
	/// kk-dependent observables, to be collected at rank tid_root_
	bool test_;
	bool use_omp_for_vie2_;
	/////////////////////////////////////////////////////////
	// before this, latt_ must be initialized !!
	mpi_lattice_step_RPA_cdw(int nt,int ntau,int size,double beta,double h,bool use_omp_for_vie2,double mu,double epsilon,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &g,CFUNC &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,int nk,int kt,bool test,double mix);
	void init_G_mat_nointeraction(void);
	void set_rpa(bool onoff);
	void extrapolate_timestep_G(int tstp,int kt);
	void extrapolate_timestep_W(int tstp,int kt);
	void extrapolate_timestep_Sigma(int tstp,int kt);
	void extrapolate_timestep_Pi(int tstp,int kt);
	void symmetrise_mat(GREEN &G);
	void step_boson(int tstp);
	void gather_gk_timestep(int tstp);
	void set_density_k(int tstp);
	void gather_Wk_timestep(int tstp);
	void gather_Dk_timestep(int tstp);
	void get_Polarization_Bubble(int tstp,int qq,GREEN &P);
	void get_Sigma_Hartree_cdw(int n,int k,CFUNC &S);
	void get_Sigma_Fock_cdw(int n,int k,CFUNC &S);
	void get_Sigma_RPA_cdw(int tstp,int kk,GREEN &S,int iter);
	void get_Gloc(int tstp);
	void get_Wloc(int tstp);
	void hkkin(cdmatrix &hkmatrix,int tstp,double kk);
	void print_diagonal_occupations(const char *file_prefix,int t1,int t2,int kt,double h,int size);
	void gather_kk_observables(int tstp, int kt);
	double get_eneRPA(int tstp,int kt);

	virtual double step(int tstp,int iter,int kt,double om0,double s,double amp) override;
	virtual void read_from_file_hdf5(int tstp,const char *filename_prefix,int kt) override;
	virtual void print_to_file_hdf5(const char *filename_prefix,int print_k);
	virtual void print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k);
};


template <class LATTICE> class mpi_lattice_step_2b_cdw: public mpi_lattice_step_cdw<LATTICE>{
public:
	typedef mpi_lattice_step_cdw<LATTICE> base_type;
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
	GREEN Gloc_,Wloc_;
	/// kk-dependent observables, to be collected at rank tid_root_
	bool test_;
	bool use_omp_for_vie2_;
	/////////////////////////////////////////////////////////
	// before this, latt_ must be initialized !!
	mpi_lattice_step_2b_cdw(int nt,int ntau,int size,double beta,double h,bool use_omp_for_vie2,double mu,double epsilon,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &g,CFUNC &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,int nk,int kt,bool test,double mix);
	void init_G_mat_nointeraction(void);
	void set_rpa(bool onoff);
	void extrapolate_timestep_G(int tstp,int kt);
	void extrapolate_timestep_W(int tstp,int kt);
	void extrapolate_timestep_Sigma(int tstp,int kt);
	void extrapolate_timestep_Pi(int tstp,int kt);
	void symmetrise_mat(GREEN &G);
	void step_boson(int tstp);
	void gather_gk_timestep(int tstp);
	void set_density_k(int tstp);
	void gather_Wk_timestep(int tstp);
	void gather_Dk_timestep(int tstp);
	void get_Polarization_Bubble(int tstp,int qq,GREEN &P);
	void get_Sigma_Hartree_cdw(int n,int k,CFUNC &S);
	void get_Sigma_Fock_cdw(int n,int k,CFUNC &S);
	void get_Sigma_RPA_cdw(int tstp,int kk,GREEN &S,int iter);
	void get_Gloc(int tstp);
	void get_Wloc(int tstp);
	void hkkin(cdmatrix &hkmatrix,int tstp,double kk);
	void print_diagonal_occupations(const char *file_prefix,int t1,int t2,int kt,double h,int size);
	void gather_kk_observables(int tstp, int kt);
	double get_eneRPA(int tstp,int kt);

	virtual double step(int tstp,int iter,int kt,double om0,double s,double amp) override;
	virtual void read_from_file_hdf5(int tstp,const char *filename_prefix,int kt) override;
	virtual void print_to_file_hdf5(const char *filename_prefix,int print_k);
	virtual void print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k);
};

template <class LATTICE> class mpi_lattice_step_2b: public mpi_lattice_step<LATTICE>{
public:
	typedef mpi_lattice_step<LATTICE> base_type;
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
	GREEN Gloc_,Wloc_,Dloc_;
	/// kk-dependent observables, to be collected at rank tid_root_
	bool test_;
	bool use_omp_for_vie2_;
	/////////////////////////////////////////////////////////
	// before this, latt_ must be initialized !!
	mpi_lattice_step_2b(int nt,int ntau,int size,double beta,double h,bool use_omp_for_vie2,double mu,double epsilon,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &g,CFUNC &omega0,CFUNC &gBATH,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,int nk,int kt,bool test,double mix);
	void init_G_mat_nointeraction(void);
	void set_rpa(bool onoff);
	void extrapolate_timestep_G(int tstp,int kt);
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
	void get_Sigma_2b(int tstp,int kk,GREEN &S,int iter);
	void get_Gloc(int tstp);
	void get_Wloc(int tstp);
	void get_Dloc(int tstp);
	void print_diagonal_occupations(const char *file_prefix,int t1,int t2,int kt,double h,int size);
	void gather_kk_observables(int tstp, int kt);
	double get_eneRPA(int tstp,int kt);

	// virtual double get_ekin(int tstp) override;
	// virtual void get_ekin(int tstp,cdmatrix &ekin) override;
	virtual double step(int tstp,int iter,int kt,double om0,double s,double amp) override;
	virtual void read_from_file_hdf5(int tstp,const char *filename_prefix,int kt) override;
	virtual void print_to_file_hdf5(const char *filename_prefix,int print_k);
	virtual void print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k);
};

template <class LATTICE> class mpi_lattice_step_2b_deb: public mpi_lattice_deb_step<LATTICE>{
public:
	typedef mpi_lattice_deb_step<LATTICE> base_type;
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
	GREEN Gloc_,Wloc_,Dloc_,D0loc_,D0hol_,Bath_;
	CFUNC g_;
	int migdal_,phonontype_;
	/// kk-dependent observables, to be collected at rank tid_root_
	bool test_;
	bool use_omp_for_vie2_;
	/////////////////////////////////////////////////////////
	// before this, latt_ must be initialized !!
	mpi_lattice_step_2b_deb(parameters &param);
	void init_G_mat_nointeraction(void);
	void set_rpa(bool onoff);
	void extrapolate_timestep_G(int tstp,int kt);
	void extrapolate_timestep_W(int tstp,int kt);
	void extrapolate_timestep_Sigma(int tstp,int kt);
	void extrapolate_timestep_Pi(int tstp,int kt);
	void extrapolate_timestep_D(int tstp,int kt);
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

	// virtual double get_ekin(int tstp) override;
	// virtual void get_ekin(int tstp,cdmatrix &ekin) override;
	virtual double step(int tstp,int iter,int kt,double om0,double s,double amp) override;
	virtual void read_from_file_hdf5(int tstp,const char *filename_prefix,int kt) override;
	virtual void print_to_file_hdf5(const char *filename_prefix,int print_k);
	virtual void print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k);
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
	virtual void print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k);
};

template <class LATTICE> class mpi_lattice_step_hubbard: public approx<LATTICE>{
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
	int suscep_;
	double eta_; //Broadening
	double mazza_;
	// int migdal_;
	/////////////////////////////////////////////////////////
	// before this, latt_ must be initialized !!
  mpi_lattice_step_hubbard(parameters_hub &param);
  mpi_lattice_step_hubbard(parameters_tns &param);
	// void init(int nt,int ntau,int size,double beta,double h,bool use_omp_for_vie2,double mu,double epsilon);
	void get_Sigma_Hartree(int tstp,int kk,CFUNC &S);
	void get_Sigma_Fock(int tstp,int kk,CFUNC &S);
	void set_local(int tstp);
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
	virtual void print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k);
	virtual double get_ekin(int tstp) override;
	virtual void get_ekin(int tstp,cdmatrix &ekin) override;
	virtual double get_curr(int tstp) override;
	virtual double get_dAk(int tstp) override;
	virtual double get_eneHF(int tstp) override;
	
};


template <class LATTICE> class mpi_lattice_step_hubbard_RPA: public mpi_lattice_step_hubbard<LATTICE>{
public:
	typedef mpi_lattice_step_hubbard<LATTICE> base_type;
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
	GREEN Gloc_,Wloc_,Dloc_,Bath_,D0loc_;
	CFUNC g_;
	/// kk-dependent observables, to be collected at rank tid_root_
	bool test_;
	bool use_omp_for_vie2_;
	int migdal_,phonontype_;
	/////////////////////////////////////////////////////////
	// before this, latt_ must be initialized !!
  	mpi_lattice_step_hubbard_RPA(parameters_hub &param);
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
	void get_Bubble(int tstp,int qq,GREEN &P);
	void get_phonon_Bubble(int tstp,int qq,GREEN &P);
	void get_Sigma_rpa(int tstp,int kk,GREEN &S,int iter);
	void get_Gloc(int tstp);
	void get_Wloc(int tstp);
	void get_Dloc(int tstp);
	void print_diagonal_occupations(const char *file_prefix,int t1,int t2,int kt,double h,int size);
	void gather_kk_observables(int tstp, int kt);
	double get_eneRPA(int tstp,int kt);

	virtual double step(int tstp,int iter,int kt,double om0,double s,double amp) override;
	virtual void read_from_file_hdf5(int tstp,const char *filename_prefix,int kt) override;
	virtual void print_to_file_hdf5(const char *filename_prefix,int print_k);
	virtual void print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k);
};

template <class LATTICE> class mpi_lattice_step_hubbard_2b: public mpi_lattice_step_hubbard<LATTICE>{
public:
	typedef mpi_lattice_step_hubbard<LATTICE> base_type;
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
	GREEN Gloc_,Wloc_,Dloc_,Bath_,D0loc_;
	CFUNC g_;
	/// kk-dependent observables, to be collected at rank tid_root_
	bool test_;
	bool use_omp_for_vie2_;
	int migdal_,phonontype_;
	/////////////////////////////////////////////////////////
	// before this, latt_ must be initialized !!
  	mpi_lattice_step_hubbard_2b(parameters_hub &param);
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
	virtual void print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k);
};


template <class LATTICE> class mpi_lattice_step_GKBA: public mpi_lattice_step<LATTICE>{
public:
	typedef mpi_lattice_step<LATTICE> base_type;
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

	bool use_omp_for_vie2_;
	/////////////////////////////////////////////////////////
	// before this, latt_ must be initialized !!
	mpi_lattice_step_GKBA(int nt,int ntau,int size,double beta,double h,bool use_omp_for_vie2,double mu,double epsilon,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &g,CFUNC &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,int nk,int kt);

	virtual double step(int tstp,int iter,int kt,double om0,double s,double amp) override{};
	virtual void read_from_file_hdf5(int tstp,const char *filename_prefix,int kt) override;
	virtual void print_to_file_hdf5(const char *filename_prefix,int print_k);
	virtual void print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k);
	virtual double get_eneHF(int tstp) override{}
};
