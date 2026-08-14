#pragma once

#include "inclusions.hpp"
#include "parameters.hpp"
#include <iostream>

lattice_1d_2b_optical_nofield_abinitio::lattice_1d_2b_optical_nofield_abinitio(void){
}

lattice_1d_2b_optical_nofield_abinitio::lattice_1d_2b_optical_nofield_abinitio(parameters &param):
	nt_(param.nt),nrpa_(param.size),xi_(param.xi),mu_(0.0)
{
	int nt=param.nt;
	assert(-1<=nt);
  	assert(1.0<=param.xi);
  	assert(-1<=param.nt);
  	assert(nt+2==(int)param.U.size());

  	tt_.resize(nt+2);
  	U_.resize(nt+2);
  	V_.resize(nt+2);
  	omega0_=param.omega0;
  	g_=param.g;
  	mazza_=param.mazza;
  	A_.resize(nt+2);
  	E_.resize(nt+2);
	dE_.resize(nt+2);
	v01_.resize(nt+2);
	delta_.resize(nt+2);
	dipol_=param.dipol;
	mu_=param.mu;
	fieldD_=param.fieldD;
	fieldP_=param.fieldP;
	ratio_=param.ratio;
	rho_eq_.resize(nrpa_,nrpa_);
	
	v01_[0]=param.v01; //Artificial symmetry breaking field only for equilibrium
	for(int tstp=-1;tstp<=nt;tstp++){
		tt_[tstp+1]=param.tt[tstp+1];
		U_[tstp+1]=param.U[tstp+1];
		V_[tstp+1]=param.V[tstp+1];
		E_[tstp+1]=param.E[tstp+1];
		dE_[tstp+1]=param.dE[tstp+1];
		v01_[tstp+1]=0.0;
		delta_[tstp+1]=param.delta[tstp+1];
		if(tstp>0){
    		// return Eseed*np.exp(-4.6*np.power(t-tProbe,2)/np.power(sProbe,2))
			double sprobe=10*param.h;
			double tprobe=5*sprobe;
			double t=(tstp+1)*param.h;
			v01_[tstp+1]=param.v01_time*exp(-4.6*(t-tprobe)*(t-tprobe)/(sprobe*sprobe));
		}
	}

	std::ofstream out;
	out.open("v01.out");
	for(int tstp=-1;tstp<=nt;tstp++){
		double t=(tstp+1)*param.h;
		out << t << " " << v01_[tstp+1] << std::endl;
	}
	efield_to_afield(nt,param.h,param.E,A_,param.kt);

	
	init_kk(param.nk);
	// TODO: check Ak, dAk kako se kliceta v programu!
}

void lattice_1d_2b_optical_nofield_abinitio::efield_to_afield(int nt,double h,std::vector<double> &efield,std::vector<double> &afield,int kt){
	int kt1=(nt>=kt ? kt : nt),n,n1,tstp;
	double At;
	assert((int)efield.size()>=nt+2);
	afield.resize(nt+2);
	afield[0]=0.0; // tstp=-1
	for(tstp=0;tstp<=nt;tstp++){
		At=0.0;
		n1=(tstp<kt1 ? kt1 : tstp);
		for(n=0;n<=n1;n++){
			At += integration::I<double>(kt1).gregory_weights(tstp,n)*efield[n+1];
		}
		afield[tstp+1]=At*(-h);
	}
}



int lattice_1d_2b_optical_nofield_abinitio::add_kpoints(int k1,int s1,int k2,int s2){
  int k12=s1*(k1-G_)+s2*(k2-G_)+G_;
  while (k12<0){k12 += nk_;}
  while (k12>=nk_){k12 -= nk_;}
  return k12;
}

void lattice_1d_2b_optical_nofield_abinitio::init_kk(int nk){
	double dk,kw;
	assert(nk>1);
	nk_=nk;
	G_=nk_/2;
	dk=2.0*PI/nk_; // kk= (index-G_)*dk 
	kw=1.0/nk_;
	kpoints_.resize(nk_);
	kweight_bz_.resize(nk_);
	for(int i1=0;i1<nk_;i1++) kpoints_[i1]=(i1-G_)*dk;
	for(int i1=0;i1<nk_;i1++) kweight_bz_[i1]=kw;

	//for(int i1=0;i1<nk_;i1++){
	//  std::cout << "Kpoints: " << kpoints_[i1] << " +  " << kpoints_[3] << " = " << kpoints_[add_kpoints(i1,1,3,1)]  <<  std::endl;
	//}
}

// interaction vertex V_{a1;a2}(q) = V_{a1,a2}(q) 
// look at PHYSICAL REVIEW B 94, 035121 (2016) for context !
double lattice_1d_2b_optical_nofield_abinitio::V(int tstp,double qq,int j1,int j2){
	if(j1!=j2){
		double xi2=xi_*xi_;
		double cq=cos(qq);
		return U_[tstp+1]+2.0*V_[tstp+1]*(xi2*cq-xi_)/(xi2-2.0*cq*xi_+1.0);
	}else{
		return 0.0;
	}
}
	
std::complex<double> lattice_1d_2b_optical_nofield_abinitio::Velph(int tstp,std::complex<double> &X,int j1,int j2){
	if(j1!=j2){
		return g_[tstp+1]*X;
	}else{
		return 0.0;
	}
}
	
void lattice_1d_2b_optical_nofield_abinitio::V(cdmatrix &Vmatrix,int tstp,double qq){
	// not efficient, but used only for initialization
	Vmatrix.resize(nrpa_,nrpa_);
	for(int j1=0;j1<nrpa_;j1++){
		for(int j2=0;j2<nrpa_;j2++){
			Vmatrix(j1,j2)=V(tstp,qq,j1,j2);
		}
	}
}
	
void lattice_1d_2b_optical_nofield_abinitio::Velph(cdmatrix &Vmatrix,int tstp,std::complex<double> &X){
		// not efficient, but used only for initialization
	Vmatrix.resize(nrpa_,nrpa_);
	for(int j1=0;j1<nrpa_;j1++){
		for(int j2=0;j2<nrpa_;j2++){
			Vmatrix(j1,j2)=Velph(tstp,X,j1,j2);
		}
	}
}

// the dispersion (particle-hole symmetric for A=0)
void lattice_1d_2b_optical_nofield_abinitio::hk(cdmatrix &hkmatrix,int tstp,double kk,int iter){
	// All terms normalized to the hopping of the valence band [0.3]
	double kkshift=kk+A_[tstp+1]*fieldP_;
	double tTa=-0.72/0.3;
	double tNi=1.0;
	double epsTa=1.35/0.3;
	double epsNi=-0.36/0.3;
	double epskTa=2.0*tTa*cos(kkshift);
	double epskNi=2.0*tNi*cos(kkshift);
	// Off-diagonal components (assumption of a/2 distance between atoms)
	cdouble ak=dipol_*E_[tstp+1]*fieldD_;
	// std::cout << "Ak " << dipol_*E_[tstp+1]*fieldD_ << std::endl;

	double vq0_01=V(tstp,0.0,0,1);
	double d0=-vq0_01*std::real(rho_eq_(1,1))-0.5*delta_[tstp+1];
	double d1=-vq0_01*std::real(rho_eq_(0,0))+0.5*delta_[tstp+1];
	// double d0=-0.5*delta_[tstp+1];
	// double d1=+0.5*delta_[tstp+1];
	// std::cout << "hk " <<  std::real(rho_eq_(0,0)) << " " << std::real(rho_eq_(1,1)) << std::endl;
	std::complex<double> I(0.0,1.0);
	// std::complex<double> maz=2.0*sqrt(2)*I*((mazza_/0.3)*sin(kk/2.0));
	std::complex<double> maz=(1.0-exp(I*kk))*(mazza_/0.3)*sqrt(2);
	// 2.0*sqrt(2)*I*((mazza_/0.3)*sin(kk/2.0));

	hkmatrix.resize(nrpa_,nrpa_);
	hkmatrix.setZero();
	hkmatrix(0,0)=epskTa+epsTa+d0;
	hkmatrix(1,1)=epskNi+epsNi+d1;
	hkmatrix(0,1)=ak+maz;
	hkmatrix(1,0)=conj(ak+maz);
	// std::cout << "off diagonal " << ak+maz << std::endl;

	if(tstp==-1 and 0<iter and iter<10){
		hkmatrix(0,1)=v01_[tstp+1]*cos(kk/2.0);
		hkmatrix(1,0)=conj(v01_[tstp+1]*cos(kk/2.0));
	}
}
	
void lattice_1d_2b_optical_nofield_abinitio::hkfree(cdmatrix &hkmatrix,int tstp,double kk){
	// All terms normalized to the hopping of the valence band [0.36]
	double tTa=-0.72/0.3;
	double tNi=1.0;
	double epsTa=-0.36/0.3;
	double epsNi=1.35/0.3;
	double epskTa=2.0*tTa*cos(kk);
	double epskNi=2.0*tNi*cos(kk);
	double d0=-0.5*delta_[tstp+1];
	double d1=0.5*delta_[tstp+1];

	// Off-diagonal components (assumption of a/2 distance between atoms)
	cdouble ak=dipol_*E_[tstp+1]*fieldD_;
	std::complex<double> I(0.0,1.0);
	std::complex<double> maz=2.0*sqrt(2)*I*(mazza_/0.3)*sin(kk/2.0);
	
	hkmatrix.resize(nrpa_,nrpa_);
	hkmatrix.setZero();
	hkmatrix(0,0)=epskTa+epsTa+d0;
	hkmatrix(1,1)=epskNi+epsNi+d1;
	hkmatrix(0,1)=maz;
	hkmatrix(1,0)=conj(maz);

}

//Velocity is still 1d vector
void lattice_1d_2b_optical_nofield_abinitio::vk(cdmatrix &vkmatrix,int tstp,double kk){
	// Diagonal terms
	double kkshift=kk+A_[tstp+1];
	double tTa=-0.72/0.3;
	double tNi=1.0;

	double vkTa=-2.0*tTa*sin(kkshift);
	double vkNi=-2.0*tNi*sin(kkshift);
	std::complex<double> I(0.0,1.0);
	std::complex<double> vkMaz=2.0*sqrt(2)*I*(mazza_/0.3);

	// double vk=2.0*tt_[tstp+1]*sin(kk);
	vkmatrix.resize(nrpa_,nrpa_);
	vkmatrix.setZero();
	vkmatrix(0,0)=vkTa;
	vkmatrix(1,1)=vkNi;
	vkmatrix(0,1)=vkMaz;
	vkmatrix(1,0)=conj(vkMaz);
}

//Velocity is still 1d vector
void lattice_1d_2b_optical_nofield_abinitio::vkFULL(cdmatrix &vkmatrix,int tstp,double kk){
	// Diagonal terms
	double kkshift=kk+A_[tstp+1];
	double vk=2.0*tt_[tstp+1]*sin(kkshift);
	cdouble vD=8.0*std::complex<double>(1.0,0.0)*tt_[tstp+1]*dipol_*cos(kk);
	// double vk=2.0*tt_[tstp+1]*sin(kk);
	vkmatrix.resize(nrpa_,nrpa_);
	vkmatrix.setZero();
	vkmatrix(0,0)=-vk;
	vkmatrix(1,1)=vk*ratio_;
	vkmatrix(0,1)=vD;
	vkmatrix(1,0)=std::conj(vD);
}

// Dipolar interaction
void lattice_1d_2b_optical_nofield_abinitio::Ak(cdmatrix &hkmatrix,int tstp,double kk){
        double kkshift=kk;
	cdouble ak=dipol_;
	// cdouble ak=I*dipol_*E_[tstp+1]*std::cos(kkshift/2.0);

	hkmatrix.resize(nrpa_,nrpa_);
	hkmatrix.setZero();
	hkmatrix(0,1)=conj(ak);
	hkmatrix(1,0)=ak;
}

// dAk interaction
void lattice_1d_2b_optical_nofield_abinitio::dAk(cdmatrix &hkmatrix,int tstp,double kk){
	// cdouble ak=dA_[tstp+1];
	// double kkshift=kk+A_[tstp+1];
	// cdouble ak=s*dipol_*dE_[tstp+1]*std::cos(kkshift/2.0);
	// hkmatrix.resize(nrpa_,nrpa_);
	// hkmatrix.setZero();
	// hkmatrix(0,1)=conj(ak);
	// hkmatrix(1,0)=ak;
}	
