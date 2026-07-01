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

// TODO: Make also here inherited class structure !!

class lattice_1d_2b_nofield{
/*///////////////////////////////////////////////////////////////////////////

1d-lattice, 2 bands and interband-transitions


setup of kpoints: 

** The model stays inversion symmetric (k <-> k-1); (there is no electric 
field, excitation via interband dipole transitons);
   
** each point in the BZ has one representatives under the point-group symmetry;
   - kpoints_[0...nk_-1] is the list of all representatives (we choose kk<=0)
   - later, Green functions are only stored for representatives


G=(nk-1),  (X]: point X included, [X: point X excluded,).

       kk             -PI   0     PI
 	
index k (full BZ)      0    G    [2G
index k (stored)       0    G ]

mappings: 
BZ -> stored:
representative_kk  (kbz<=G ? kbz : 2G-kbz) 

///////////////////////////////////////////////////////////////////////////*/
public:
	int nt_;    
	int nk_; // the RBZ  
	int G_; // = index of Gamma point in BZ and stored sector
	int nk_bz_;

	std::vector<double>  kpoints_bz_;  // 0...nk_bz_-1
	std::vector<double>  kpoints_;  // 0...nk_bz_-1
	std::vector<double>  kweight_bz_;  // 0...nk_bz_-1, normalized to 1
	//////////////////////////////////////////////////////////////////////
	// MODEL PARAMETERS
	int nrpa_;
	double mu_;
	std::vector<double>   tt_; // time-dependent hopping
	std::vector<double>   U_; 
	std::vector<double>   V_;
	CFUNC   g_,omega0_,gBATH_; 
	double xi_; // interaction falls off as V_l = V_ * (1/xi)^{|l|-1} for l>0
	std::vector<double>   delta_; // band splitting 
	std::vector<cdouble>  A_,dA_; //
	std::vector<cdouble>  v01_; // exciton pair field

	lattice_1d_2b_nofield(void);
	lattice_1d_2b_nofield(int nk,int nt,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &g,CFUNC &omega0,CFUNC &gBATH,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,double mu,int size,int kt,double h);
	int representative_kk(int kbz);
	int add_kpoints(int k1,int s1,int k2,int s2);
	void init_kk(int nk);
	double V(int tstp,double qq,int j1,int j2);
	std::complex<double> Velph(int tstp,std::complex<double> &X,int j1,int j2);
	void V(cdmatrix &Vmatrix,int tstp,double qq);
	void Velph(cdmatrix &Vmatrix,int tstp,std::complex<double> &X);
	// the dispersion (particle-hole symmetric for A=0)
	void hk(cdmatrix &hkmatrix,int tstp,double kk,int iter);
	void hkfree(cdmatrix &hkmatrix,int tstp,double kk);

	//Velocity is still 1d vector
	void vk(cdmatrix &vkmatrix,int tstp,double kk);

	// Dipolar interaction
	void Ak(cdmatrix &hkmatrix,int tstp,double kk);

	// dAk interaction
	void dAk(cdmatrix &hkmatrix,int tstp,double kk);
	
	////////////////////////////////////////////////////////
	
	// void init(int nk,int nt,double tt,double U,double V,double g,double omega0,double xi,double delta);
	// void init(int nk,int nt,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,std::vector<double> &g,std::vector<double> &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,double mu);
};


class lattice_1d_2b_nofield_deb{
/*///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////*/
public:
	int nt_;    
	int nk_; // the RBZ  
	int G_; // = index of Gamma point in BZ and stored sector
	int nk_bz_;

	std::vector<double>  kpoints_bz_;  // 0...nk_bz_-1
	std::vector<double>  kpoints_;  // 0...nk_bz_-1
	std::vector<double>  kweight_bz_;  // 0...nk_bz_-1, normalized to 1
	//////////////////////////////////////////////////////////////////////
	// MODEL PARAMETERS
	int nrpa_;
	double mu_,mazza_;
	std::vector<double>   tt_; // time-dependent hopping
	std::vector<double>   U_; 
	std::vector<double>   V_;
	CFUNC   g_,omega0_,gBATH_; 
	double xi_; // interaction falls off as V_l = V_ * (1/xi)^{|l|-1} for l>0
	std::vector<double>   delta_; // band splitting 
	std::vector<double>  E_,dE_,A_;
	std::vector<cdouble>  v01_; // exciton pair field

	lattice_1d_2b_nofield_deb(void);
	lattice_1d_2b_nofield_deb(parameters &param);
	void efield_to_afield(int nt,double h,std::vector<double> &efield,std::vector<double> &afield,int kt);
	int representative_kk(int kbz);
	int add_kpoints(int k1,int s1,int k2,int s2);
	void init_kk(int nk);
	double V(int tstp,double qq,int j1,int j2);
	std::complex<double> Velph(int tstp,std::complex<double> &X,int j1,int j2);
	void V(cdmatrix &Vmatrix,int tstp,double qq);
	void Velph(cdmatrix &Vmatrix,int tstp,std::complex<double> &X);
	// the dispersion (particle-hole symmetric for A=0)
	void hk(cdmatrix &hkmatrix,int tstp,double kk,int iter);
	void hkfree(cdmatrix &hkmatrix,int tstp,double kk);

	//Velocity is still 1d vector
	void vk(cdmatrix &vkmatrix,int tstp,double kk);

	// Dipolar interaction
	void Ak(cdmatrix &hkmatrix,int tstp,double kk);

	// dAk interaction
	void dAk(cdmatrix &hkmatrix,int tstp,double kk);
	
	////////////////////////////////////////////////////////
	
	// void init(int nk,int nt,double tt,double U,double V,double g,double omega0,double xi,double delta);
	// void init(int nk,int nt,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,std::vector<double> &g,std::vector<double> &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,double mu);
};



class lattice_1d_rpa_cdw_nofield{
/*///////////////////////////////////////////////////////////////////////////

1d-lattice, 2 bands and interband-transitions


setup of kpoints: 

** The model stays inversion symmetric (k <-> k-1); (there is no electric 
field, excitation via interband dipole transitons);
   
** each point in the BZ has one representatives under the point-group symmetry;
   - kpoints_[0...nk_-1] is the list of all representatives (we choose kk<=0)
   - later, Green functions are only stored for representatives


G=(nk-1),  (X]: point X included, [X: point X excluded,).

       kk             -PI   0     PI
 	
index k (full BZ)      0    G    [2G
index k (stored)       0    G ]

mappings: 
BZ -> stored:
representative_kk  (kbz<=G ? kbz : 2G-kbz) 

///////////////////////////////////////////////////////////////////////////*/
public:
	int nt_;    
	int nk_; // the RBZ  
	int G_; // = index of Gamma point in BZ and stored sector
	int nk_rbz_,nk_fbz_,nk_bz_;

	std::vector<double>  kpoints_rbz_,kpoints_fbz_;  // 0...nk_bz_-1
	std::vector<double>  kpoints_;  // 0...nk_bz_-1
	std::vector<double>  kweight_bz_;  // 0...nk_bz_-1, normalized to 1
	//////////////////////////////////////////////////////////////////////
	// MODEL PARAMETERS
	int nrpa_,norb_;
	double mu_;
	std::vector<double>   tt_; // time-dependent hopping
	std::vector<double>   U_; 
	std::vector<double>   V_;
	CFUNC   g_,omega0_; 
	double xi_,Q_; // interaction falls off as V_l = V_ * (1/xi)^{|l|-1} for l>0
	std::vector<double>   delta_; // band splitting 
	std::vector<cdouble>  A_,dA_; //
	std::vector<cdouble>  v01_; // exciton pair field

	lattice_1d_rpa_cdw_nofield(void);
	lattice_1d_rpa_cdw_nofield(int nk,int nt,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &g,CFUNC &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,double mu,int size,int kt,double h);
	int representative_kk(int krbz);
	int add_kpoints_fbz(int s1,int k1,int s2,int k2);
	int add_kpoints(int s1,int k1,int s2,int k2);
	int add_kpoints_rbz(int &k12, int &gamma,int s1,int k1,int s2,int k2);
	void decompose_kfbz(int kfbz,int &krbz,int &gamma);
	int k_fbz(int krbz);
	void init_kk(int nk);
	double V(int tstp,double qq,int j1,int j2);
	std::complex<double> Velph(int tstp,std::complex<double> &X,int j1,int j2);
	void V(cdmatrix &Vmatrix,int tstp,double qq);
	void Velph(cdmatrix &Vmatrix,int tstp,std::complex<double> &X);
	// the dispersion (particle-hole symmetric for A=0)
	void hk(cdmatrix &hkmatrix,int tstp,double kk,int iter);
	void hkfree(cdmatrix &hkmatrix,int tstp,double kk);
	void hkkin(cdmatrix &hkmatrix,int tstp,double kk);

	int idx(int orb,int sl);
	int idx_a(int j);
	int idx_c(int j);
	//Velocity is still 1d vector
	void vk(cdmatrix &vkmatrix,int tstp,double kk);

	// Dipolar interaction
	void Ak(cdmatrix &hkmatrix,int tstp,double kk);

	// dAk interaction
	void dAk(cdmatrix &hkmatrix,int tstp,double kk);
	
	////////////////////////////////////////////////////////
	
	// void init(int nk,int nt,double tt,double U,double V,double g,double omega0,double xi,double delta);
	// void init(int nk,int nt,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,std::vector<double> &g,std::vector<double> &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,double mu);
};



class lattice_1d_2b_optical_nofield{
/*///////////////////////////////////////////////////////////////////////////

1d-lattice, 2 bands and interband-transitions


setup of kpoints: 

** The model stays inversion symmetric (k <-> k-1); (there is no electric 
field, excitation via interband dipole transitons);
   
** each point in the BZ has one representatives under the point-group symmetry;
   - kpoints_[0...nk_-1] is the list of all representatives (we choose kk<=0)
   - later, Green functions are only stored for representatives


G=(nk-1),  (X]: point X included, [X: point X excluded,).

       kk             -PI   0     PI
 	
index k (full BZ)      0    G    [2G
index k (stored)       0    G ]

mappings: 
BZ -> stored:
representative_kk  (kbz<=G ? kbz : 2G-kbz) 

///////////////////////////////////////////////////////////////////////////*/
public:
	int nt_;    
	int nk_; // the RBZ  
	int G_; // = index of Gamma point in BZ and stored sector
	int nk_bz_;

	std::vector<double>  kpoints_bz_;  // 0...nk_bz_-1
	std::vector<double>  kpoints_;  // 0...nk_bz_-1
	std::vector<double>  kweight_bz_;  // 0...nk_bz_-1, normalized to 1
	//////////////////////////////////////////////////////////////////////
	// MODEL PARAMETERS
	parameters params_;
	int nrpa_;
	double mu_,mazza_;

	std::vector<double>   tt_; // time-dependent hopping
	std::vector<double>   U_; 
	std::vector<double>   V_;
	CFUNC   g_,omega0_; 
	double fieldP_,fieldD_; // We can choose to have only one
	double xi_; // interaction falls off as V_l = V_ * (1/xi)^{|l|-1} for l>0
	std::vector<double>   delta_; // band splitting 
	std::vector<double>  E_,dE_,A_; //
	std::vector<cdouble>  v01_; // exciton pair field
	std::complex<double> dipol_;
	double ratio_;

	lattice_1d_2b_optical_nofield(void);
	lattice_1d_2b_optical_nofield(parameters &param);
	void efield_to_afield(int nt,double h,std::vector<double> &efield,std::vector<double> &afield,int kt);
	int representative_kk(int kbz);
	int add_kpoints(int k1,int s1,int k2,int s2);
	void init_kk(int nk);
	double V(int tstp,double qq,int j1,int j2);
	std::complex<double> Velph(int tstp,std::complex<double> &X,int j1,int j2);
	void V(cdmatrix &Vmatrix,int tstp,double qq);
	void Velph(cdmatrix &Vmatrix,int tstp,std::complex<double> &X);
	// the dispersion (particle-hole symmetric for A=0)
	void hk(cdmatrix &hkmatrix,int tstp,double kk,int iter);
	void hkfree(cdmatrix &hkmatrix,int tstp,double kk);
	//Velocity is still 1d vector
	void vk(cdmatrix &vkmatrix,int tstp,double kk);
        void vkFULL(cdmatrix &vkmatrix,int tstp,double kk); //Peierls+Dipolar
	// Dipolar interaction
	void Ak(cdmatrix &hkmatrix,int tstp,double kk);

	// dAk interaction
	void dAk(cdmatrix &hkmatrix,int tstp,double kk);
	
	////////////////////////////////////////////////////////
	
	// void init(int nk,int nt,double tt,double U,double V,double g,double omega0,double xi,double delta);
	// void init(int nk,int nt,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,std::vector<double> &g,std::vector<double> &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,double mu);
};


class lattice_1d_2b_optical_nofield_abinitio{
/*///////////////////////////////////////////////////////////////////////////

1d-lattice, 2 bands and interband-transitions


setup of kpoints: 

** The model stays inversion symmetric (k <-> k-1); (there is no electric 
field, excitation via interband dipole transitons);
   
** each point in the BZ has one representatives under the point-group symmetry;
   - kpoints_[0...nk_-1] is the list of all representatives (we choose kk<=0)
   - later, Green functions are only stored for representatives


G=(nk-1),  (X]: point X included, [X: point X excluded,).

       kk             -PI   0     PI
 	
index k (full BZ)      0    G    [2G
index k (stored)       0    G ]

mappings: 
BZ -> stored:
representative_kk  (kbz<=G ? kbz : 2G-kbz) 

///////////////////////////////////////////////////////////////////////////*/
public:
	int nt_;    
	int nk_; // the RBZ  
	int G_; // = index of Gamma point in BZ and stored sector
	int nk_bz_;

	std::vector<double>  kpoints_bz_;  // 0...nk_bz_-1
	std::vector<double>  kpoints_;  // 0...nk_bz_-1
	std::vector<double>  kweight_bz_;  // 0...nk_bz_-1, normalized to 1
	//////////////////////////////////////////////////////////////////////
	// MODEL PARAMETERS
	parameters params_;
	int nrpa_;
	double mu_,mazza_;

	std::vector<double>   tt_; // time-dependent hopping
	std::vector<double>   U_; 
	std::vector<double>   V_;
	CFUNC   g_,omega0_; 
	CFUNC   gBATH_,omegaBATH_;
	double fieldP_,fieldD_; // We can choose to have only one
	double xi_; // interaction falls off as V_l = V_ * (1/xi)^{|l|-1} for l>0
	std::vector<double>   delta_; // band splitting 
	std::vector<double>  E_,dE_,A_; //
	std::vector<cdouble>  v01_; // exciton pair field
	std::complex<double> dipol_;
	cdmatrix rho_eq_;
	double ratio_;

	lattice_1d_2b_optical_nofield_abinitio(void);
	lattice_1d_2b_optical_nofield_abinitio(parameters &param);
	void efield_to_afield(int nt,double h,std::vector<double> &efield,std::vector<double> &afield,int kt);
	int representative_kk(int kbz);
	int add_kpoints(int k1,int s1,int k2,int s2);
	void init_kk(int nk);
	double V(int tstp,double qq,int j1,int j2);
	std::complex<double> Velph(int tstp,std::complex<double> &X,int j1,int j2);
	void V(cdmatrix &Vmatrix,int tstp,double qq);
	void Velph(cdmatrix &Vmatrix,int tstp,std::complex<double> &X);
	// the dispersion (particle-hole symmetric for A=0)
	void hk(cdmatrix &hkmatrix,int tstp,double kk,int iter);
	void hkfree(cdmatrix &hkmatrix,int tstp,double kk);
	//Velocity is still 1d vector
	void vk(cdmatrix &vkmatrix,int tstp,double kk);
        void vkFULL(cdmatrix &vkmatrix,int tstp,double kk); //Peierls+Dipolar
	// Dipolar interaction
	void Ak(cdmatrix &hkmatrix,int tstp,double kk);

	// dAk interaction
	void dAk(cdmatrix &hkmatrix,int tstp,double kk);
	
	////////////////////////////////////////////////////////
	
	// void init(int nk,int nt,double tt,double U,double V,double g,double omega0,double xi,double delta);
	// void init(int nk,int nt,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,std::vector<double> &g,std::vector<double> &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,double mu);
};

class lattice_1d_hubbard{
/*///////////////////////////////////////////////////////////////////////////

1d-lattice, 2 bands and Hubbard-interaction


setup of kpoints: 

** The model stays inversion symmetric (k <-> k-1); (there is no electric 
field, excitation via interband dipole transitons);
   
** each point in the BZ has one representatives under the point-group symmetry;
   - kpoints_[0...nk_-1] is the list of all representatives (we choose kk<=0)
   - later, Green functions are only stored for representatives


G=(nk-1),  (X]: point X included, [X: point X excluded,).

       kk             -PI   0     PI
 	
index k (full BZ)      0    G    [2G
index k (stored)       0    G ]

mappings: 
BZ -> stored:
representative_kk  (kbz<=G ? kbz : 2G-kbz) 

///////////////////////////////////////////////////////////////////////////*/
public:
	int nt_;    
	int nk_; // the RBZ  
	int G_; // = index of Gamma point in BZ and stored sector
	int nk_bz_;

	std::vector<double>  kpoints_bz_;  // 0...nk_bz_-1
	std::vector<double>  kpoints_;  // 0...nk_bz_-1
	std::vector<double>  kweight_bz_;  // 0...nk_bz_-1, normalized to 1
	//////////////////////////////////////////////////////////////////////
	// MODEL PARAMETERS
	parameters params_;
	int nrpa_;
	double mu_,mazza_;

	std::vector<double>   tt_; // time-dependent hopping
	std::vector<double>   U_; 
	std::vector<double>   V_;
	CFUNC   g_,omega0_; 
	double fieldP_,fieldD_; // We can choose to have only one
	double xi_; // interaction falls off as V_l = V_ * (1/xi)^{|l|-1} for l>0
	std::vector<double>   delta_; // band splitting 
	std::vector<double>  E_,dE_,A_; //
	std::vector<cdouble>  v01_; // exciton pair field
	std::complex<double> dipol_;
	cdmatrix rho_eq_;
	double ratio_;

	lattice_1d_hubbard(void);
	lattice_1d_hubbard(parameters_hub &param);
	void efield_to_afield(int nt,double h,std::vector<double> &efield,std::vector<double> &afield,int kt);
	int representative_kk(int kbz);
	int add_kpoints(int k1,int s1,int k2,int s2);
	void init_kk(int nk);
	double V(int tstp,double qq,int j1,int j2);
	std::complex<double> Velph(int tstp,std::complex<double> &X,int j1,int j2);
	void V(cdmatrix &Vmatrix,int tstp,double qq);
	void Velph(cdmatrix &Vmatrix,int tstp,std::complex<double> &X);
	// the dispersion (particle-hole symmetric for A=0)
	void hk(cdmatrix &hkmatrix,int tstp,double kk,int iter);
	void hkfree(cdmatrix &hkmatrix,int tstp,double kk);
	//Velocity is still 1d vector
	void vk(cdmatrix &vkmatrix,int tstp,double kk);
        void vkFULL(cdmatrix &vkmatrix,int tstp,double kk); //Peierls+Dipolar
	// Dipolar interaction
	void Ak(cdmatrix &hkmatrix,int tstp,double kk);

	// dAk interaction
	void dAk(cdmatrix &hkmatrix,int tstp,double kk);
	
	////////////////////////////////////////////////////////
	
	// void init(int nk,int nt,double tt,double U,double V,double g,double omega0,double xi,double delta);
	// void init(int nk,int nt,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,std::vector<double> &g,std::vector<double> &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,double mu);
};


class lattice_1d_1b_nofield{
/*///////////////////////////////////////////////////////////////////////////

1d-lattice, 2 bands and interband-transitions


setup of kpoints: 

** The model stays inversion symmetric (k <-> k-1); (there is no electric 
field, excitation via interband dipole transitons);
   
** each point in the BZ has one representatives under the point-group symmetry;
   - kpoints_[0...nk_-1] is the list of all representatives (we choose kk<=0)
   - later, Green functions are only stored for representatives


G=(nk-1),  (X]: point X included, [X: point X excluded,).

       kk             -PI   0     PI
 	
index k (full BZ)      0    G    [2G
index k (stored)       0    G ]

mappings: 
BZ -> stored:
representative_kk  (kbz<=G ? kbz : 2G-kbz) 

///////////////////////////////////////////////////////////////////////////*/
public:
	int nt_;    
	int nk_; // the RBZ  
	int G_; // = index of Gamma point in BZ and stored sector
	int nk_bz_;

	std::vector<double>  kpoints_bz_;  // 0...nk_bz_-1
	std::vector<double>  kpoints_;  // 0...nk_bz_-1
	std::vector<double>  kweight_bz_;  // 0...nk_bz_-1, normalized to 1
	//////////////////////////////////////////////////////////////////////
	// MODEL PARAMETERS
	int nrpa_;
	double mu_;
	std::vector<double>   tt_; // time-dependent hopping
	std::vector<double>   U_; 
	std::vector<double>   V_;
	CFUNC  g_,omega0_; 
	double xi_; // interaction falls off as V_l = V_ * (1/xi)^{|l|-1} for l>0
	std::vector<double>   delta_; // band splitting 
	std::vector<double>  A_,dA_; // Peierls phase
	std::vector<double>  E_; // Electric field
	std::vector<cdouble>  v01_; // exciton pair field

	lattice_1d_1b_nofield(void);
	lattice_1d_1b_nofield(int nk,int nt,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &g,CFUNC &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,double mu,int size,int kt,double h);
	int representative_kk(int kbz);
	int add_kpoints(int k1,int s1,int k2,int s2);
	void init_kk(int nk);
	double V(int tstp,double qq,int j1,int j2);
	std::complex<double> Velph(int tstp,std::complex<double> &X,int j1,int j2);
	void V(cdmatrix &Vmatrix,int tstp,double qq);
	void Velph(cdmatrix &Vmatrix,int tstp,std::complex<double> &X);
	// the dispersion (particle-hole symmetric for A=0)
	void hk(cdmatrix &hkmatrix,int tstp,double kk,int iter);
	void hkfree(cdmatrix &hkmatrix,int tstp,double kk);
	void efield_to_afield(int nt,double h,std::vector<double> &efield,std::vector<double> &afield,int kt);
	//Velocity is still 1d vector
	void vk(cdmatrix &vkmatrix,int tstp,double kk);

	// Dipolar interaction
	void Ak(cdmatrix &hkmatrix,int tstp,double kk);

	// dAk interaction
	void dAk(cdmatrix &hkmatrix,int tstp,double kk);
	
	////////////////////////////////////////////////////////
	
	// void init(int nk,int nt,double tt,double U,double V,double g,double omega0,double xi,double delta);
	// void init(int nk,int nt,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,std::vector<double> &g,std::vector<double> &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,double mu);
};






// class lattice_2d_tns{
// /*///////////////////////////////////////////////////////////////////////////

// 2d-lattice, 6 bands and Hubbard-interaction

// ///////////////////////////////////////////////////////////////////////////*/
// public:
// 	int nt_;    
// 	int nk_;
// 	int G_; // = index of Gamma point in BZ and stored sector
// 	int nkx_;  // number of k-points along x-direction

// 	std::vector<dvector>  kpoints_;  // 0...nk_bz_-1
// 	std::vector<double>  kweight_;  // 0...nk_bz_-1, normalized to 1
// 	//////////////////////////////////////////////////////////////////////
// 	// MODEL PARAMETERS
// 	parameters_tns params_;
// 	int nrpa_;
// 	double mu_,mazza_;

// 	std::vector<double>   tt_; // time-dependent hopping
// 	std::vector<double>   U_; 
// 	std::vector<double>   V_;
// 	CFUNC   g_,omega0_; 
// 	// double fieldP_,fieldD_; // We can choose to have only one
// 	double xi_; // interaction falls off as V_l = V_ * (1/xi)^{|l|-1} for l>0
// 	std::vector<double>   delta_; // band splitting 
// 	std::vector<double>  Ex_,Ax_,Ey_,Ay_; //
// 	std::vector<cdouble>  v01_; // exciton pair field
// 	std::complex<double> dipol_;
// 	cdmatrix rho_eq_;
// 	double ratio_;

// 	lattice_2d_tns(void);
// 	lattice_2d_tns(parameters_tns &param);
// 	void efield_to_afield(int nt,double h,std::vector<double> &E,std::vector<double> &A,int kt);
// 	int add_kpoints(int k1,int s1,int k2,int s2);
// 	void init_kk(int nk);
// 	std::complex<double> V(int tstp,dvector &qq,int j1,int j2);
// 	// std::complex<double> Velph(int tstp,std::complex<double> &X,int j1,int j2);
// 	void V(cdmatrix &Vmatrix,int tstp,dvector &qq);
// 	void Velph(cdmatrix &Vmatrix,int tstp,std::complex<double> &X);
// 	// the dispersion (particle-hole symmetric for A=0)
// 	void hk(cdmatrix &hkmatrix,int tstp,dvector &kk,int iter);
// 	void hkfree(cdmatrix &hkmatrix,int tstp,dvector &kk);
// 	//Velocity is still 1d vector
// 	void vk(cdmatrix &vkmatrix,int tstp,double kk);
//         void vkFULL(cdmatrix &vkmatrix,int tstp,double kk); //Peierls+Dipolar
// 	// Dipolar interaction
// 	void Ak(cdmatrix &hkmatrix,int tstp,double kk);
// 	// dAk interaction
// 	void dAk(cdmatrix &hkmatrix,int tstp,double kk);
// };
