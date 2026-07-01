#pragma once

#include "../step/step.hpp"

lattice_1d_2b_nofield_deb::lattice_1d_2b_nofield_deb(void){
}

lattice_1d_2b_nofield_deb::lattice_1d_2b_nofield_deb(parameters &param):
	nt_(param.nt),nrpa_(param.size),xi_(param.xi),mu_(param.mu)
{
  assert(-1<=param.nt);
  assert(1.0<=param.xi);
  assert(-1<=param.nt);
  assert(param.nt+2==(int)param.U.size());

	tt_.resize(param.nt+2);
	U_.resize(param.nt+2);
	V_.resize(param.nt+2);
	omega0_=param.omega0;
	g_=param.g;
	E_.resize(param.nt+2);
	dE_.resize(param.nt+2);
	v01_.resize(param.nt+2);
	delta_.resize(param.nt+2);

	mu_=param.mu;
	mazza_=param.mazza;
	
	for(int tstp=-1;tstp<=param.nt;tstp++){
		tt_[tstp+1]=param.tt[tstp+1];
		U_[tstp+1]=param.U[tstp+1];
		V_[tstp+1]=param.V[tstp+1];
		E_[tstp+1]=param.E[tstp+1];
		dE_[tstp+1]=param.dE[tstp+1];
		v01_[tstp+1]=0.0;
		delta_[tstp+1]=param.delta[tstp+1];
	}
	for(int k=0;k<param.nt;k++){
		std::cout << "Lattice: " << param.tt[k] << " " << param.U[k] << " " << param.V[k] << " " << param.E[k] << " " << param.dE[k] << " " << param.delta[k] << std::endl;
	}
	efield_to_afield(param.nt,param.h,param.E,A_,param.kt);

	v01_[0]=param.v01; //Artificial symmetry breaking field only for equilibrium
	init_kk(param.nk);
}

void lattice_1d_2b_nofield_deb::efield_to_afield(int nt,double h,std::vector<double> &efield,std::vector<double> &afield,int kt){
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

int lattice_1d_2b_nofield_deb::add_kpoints(int k1,int s1,int k2,int s2){
  int k12=s1*(k1-G_)+s2*(k2-G_)+G_;
  while (k12<0){k12 += nk_;}
  while (k12>=nk_){k12 -= nk_;}
  return k12;
}

void lattice_1d_2b_nofield_deb::init_kk(int nk){
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

	// for(int i1=0;i1<nk_;i1++){
	//  std::cout << "Kpoints: " << kpoints_[i1] << " " << kw  <<  std::endl;
	// }
}

// interaction vertex V_{a1;a2}(q) = V_{a1,a2}(q) 
double lattice_1d_2b_nofield_deb::V(int tstp,double qq,int j1,int j2){
	if(j1!=j2){
		double xi2=xi_*xi_;
		double cq=cos(qq);
		return U_[tstp+1]+2.0*V_[tstp+1]*(xi2*cq-xi_)/(xi2-2.0*cq*xi_+1.0);
	}else{
		return 0.0;
	}
}
	
std::complex<double> lattice_1d_2b_nofield_deb::Velph(int tstp,std::complex<double> &X,int j1,int j2){
	if(j1!=j2){
		return g_[tstp+1]*X;
	}else{
		return 0.0;
	}
}
	
void lattice_1d_2b_nofield_deb::V(cdmatrix &Vmatrix,int tstp,double qq){
	// not efficient, but used only for initialization
	Vmatrix.resize(nrpa_,nrpa_);
	for(int j1=0;j1<nrpa_;j1++){
		for(int j2=0;j2<nrpa_;j2++){
			Vmatrix(j1,j2)=V(tstp,qq,j1,j2);
		}
	}
}
	
void lattice_1d_2b_nofield_deb::Velph(cdmatrix &Vmatrix,int tstp,std::complex<double> &X){
		// not efficient, but used only for initialization
	Vmatrix.resize(nrpa_,nrpa_);
	for(int j1=0;j1<nrpa_;j1++){
		for(int j2=0;j2<nrpa_;j2++){
			Vmatrix(j1,j2)=Velph(tstp,X,j1,j2);
		}
	}
}

// the dispersion (particle-hole symmetric for A=0)
void lattice_1d_2b_nofield_deb::hk(cdmatrix &hkmatrix,int tstp,double kk,int iter){
	// All terms normalized to the hopping of the valence band [0.3]
	double fieldP=1.0;double fieldD=1.0;double dipol=1.0;double mazza=0.116;

	double vq0_00=V(tstp,0.0,0,0);
	double vq0_01=V(tstp,0.0,0,1);

	// TODO: if you want proper gauge
	double kkshift=kk+A_[tstp+1]*fieldP;
	double tTa=-0.72/0.3;
	double tNi=1.0;
	double epsTa=-2.0*tTa*cos(0.0)+0.5*delta_[tstp+1]-0.5*(vq0_00+vq0_01);
	double epsNi=-2.0*tNi*cos(0.0)-0.5*delta_[tstp+1]-0.5*(vq0_00+vq0_01);
	double epskTa=2.0*tTa*cos(kkshift);
	double epskNi=2.0*tNi*cos(kkshift);
	// Off-diagonal components (assumption of a/2 distance between atoms)
	cdouble ak=dipol*E_[tstp+1]*fieldD;

	std::complex<double> I(0.0,1.0);
	std::complex<double> maz=2.0*sqrt(2)*I*((mazza/0.3)*exp(I*kk/2.0)*sin(kk/2.0));

	hkmatrix.resize(nrpa_,nrpa_);
	hkmatrix.setZero();
	hkmatrix(0,0)=epskTa+epsTa;
	hkmatrix(1,1)=epskNi+epsNi;
	hkmatrix(0,1)=ak+maz;
	hkmatrix(1,0)=conj(ak+maz);
	// std::cout << " kkshift " << kk << " " << kkshift << " " << ak<< std::endl;

	if(tstp==-1 and 0<iter and iter<10){
		hkmatrix(0,1)=v01_[tstp+1]+maz;
		hkmatrix(1,0)=conj(v01_[tstp+1]+maz);
	}
}
	
void lattice_1d_2b_nofield_deb::hkfree(cdmatrix &hkmatrix,int tstp,double kk){
	// All terms normalized to the hopping of the valence band [0.36]
	double fieldP=1.0;double fieldD=1.0;double dipol=1.0;double mazza=0.116;

	double tTa=-0.72/0.3;
	double tNi=1.0;
	double epsTa=-2.0*tTa*cos(0.0)+0.5*delta_[tstp+1];
	double epsNi=-2.0*tNi*cos(0.0)-0.5*delta_[tstp+1];
	double epskTa=2.0*tTa*cos(kk);
	double epskNi=2.0*tNi*cos(kk);
	
	// Off-diagonal components (assumption of a/2 distance between atoms)
	std::complex<double> I(0.0,1.0);
	std::complex<double> maz=2.0*sqrt(2)*I*(mazza/0.3)*exp(I*kk/2.0)*sin(kk/2.0);
	
	hkmatrix.resize(nrpa_,nrpa_);
	hkmatrix.setZero();
	hkmatrix(0,0)=epskTa+epsTa;
	hkmatrix(1,1)=epskNi+epsNi;
	hkmatrix(0,1)=maz;
	hkmatrix(1,0)=conj(maz);

}

//Velocity is still 1d vector
void lattice_1d_2b_nofield_deb::vk(cdmatrix &vkmatrix,int tstp,double kk){
	double fieldP=1.0;double fieldD=1.0;double dipol=1.0;double mazza=0.116;
	// Diagonal terms
	double kkshift=kk+A_[tstp+1];
	// double kkshift=kk;
	double tTa=-0.72/0.3;
	double tNi=1.0;

	double vkTa=-2.0*tTa*sin(kkshift);
	double vkNi=-2.0*tNi*sin(kkshift);
	std::complex<double> I(0.0,1.0);
	std::complex<double> vkMaz=2.0*sqrt(2)*I*(mazza/0.3)*exp((I*kk)/2.0);

	// double vk=2.0*tt_[tstp+1]*sin(kk);
	vkmatrix.resize(nrpa_,nrpa_);
	vkmatrix.setZero();
	vkmatrix(0,0)=vkTa;
	vkmatrix(1,1)=vkNi;
	vkmatrix(0,1)=vkMaz;
	vkmatrix(1,0)=conj(vkMaz);
}

// Dipolar interaction
void lattice_1d_2b_nofield_deb::Ak(cdmatrix &hkmatrix,int tstp,double kk){
	// cdouble ak=A_[tstp+1];
	// hkmatrix.resize(nrpa_,nrpa_);
	// hkmatrix.setZero();
	// cdouble ak=dipol_;
	// hkmatrix(0,1)=conj(ak);
	// hkmatrix(1,0)=ak;
}

// dAk interaction
void lattice_1d_2b_nofield_deb::dAk(cdmatrix &hkmatrix,int tstp,double kk){
	// cdouble ak=dA_[tstp+1];
	// hkmatrix.resize(nrpa_,nrpa_);
	// hkmatrix.setZero();
	// hkmatrix(0,1)=conj(ak);
	// hkmatrix(1,0)=ak;
}	
