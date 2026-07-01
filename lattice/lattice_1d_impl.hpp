#pragma once

#include "../step/step.hpp"
lattice_1d_nofield::lattice_1d_nofield(void){
}

lattice_1d_nofield::lattice_1d_nofield(int nk,int nt,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &g,CFUNC &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,double mu,int size):
	nt_(nt),nrpa_(size),xi_(xi),mu_(0.0)
{
	CNTR_ASSERT_LESEQ(ASSERT_0,-1,nt,__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ(ASSERT_0,1.00001,xi,__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ(ASSERT_0,-1,nt,__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ(ASSERT_0,1.0,xi,__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ(ASSERT_0,nt+2,(int)U.size(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ(ASSERT_0,nt+2,(int)V.size(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ(ASSERT_0,nt+2,(int)tt.size(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ(ASSERT_0,nt+2,(int)A.size(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ(ASSERT_0,nt+2,(int)dA.size(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ(ASSERT_0,nt+2,(int)delta.size(),__PRETTY_FUNCTION__)

	tt_.resize(nt+2);
	U_.resize(nt+2);
	V_.resize(nt+2);
	omega0_=omega0;
	g_=g;
	A_.resize(nt+2);
	dA_.resize(nt+2);
	v01_.resize(nt+2);
	delta_.resize(nt+2);

	mu_=mu;
	
	for(int tstp=-1;tstp<=nt;tstp++){
		tt_[tstp+1]=tt[tstp+1];
		U_[tstp+1]=U[tstp+1];
		V_[tstp+1]=V[tstp+1];
		A_[tstp+1]=A[tstp+1];
		dA_[tstp+1]=dA[tstp+1];
		v01_[tstp+1]=0.0;
		delta_[tstp+1]=delta[tstp+1];
	}
	v01_[0]=v01; //Artificial symmetry breaking field only for equilibrium
	init_kk(nk);
}

int lattice_1d_nofield::representative_kk(int kbz){
	CNTR_ASSERT_LESEQ_3(ASSERT_0,0,kbz,nk_bz_-1,__PRETTY_FUNCTION__)	
	return ( kbz <= G_ ? kbz  : 2*G_-kbz);
}

int lattice_1d_nofield::add_kpoints(int k1,int s1,int k2,int s2){
	int k12=s1*(k1-G_)+s2*(k2-G_)+G_;
	while (k12<0){k12 += nk_bz_;}
	while (k12>=nk_bz_){k12 -= nk_bz_;}
	return k12;
}

void lattice_1d_nofield::init_kk(int nk){
	double dk,kw;
	CNTR_ASSERT(ASSERT_0,nk>1,__PRETTY_FUNCTION__)
	nk_=nk;
	G_=nk_-1;
	nk_bz_=2*G_;
	dk=2.0*PI/nk_bz_; // kk= (index-G_)*dk 
	kw=1.0/nk_bz_;
	kpoints_.resize(nk_);
	kweight_bz_.resize(nk_bz_);	
	kpoints_bz_.resize(nk_bz_);
	for(int i1=0;i1<nk_;i1++) kpoints_[i1]=(i1-G_)*dk;
	for(int i1=0;i1<nk_bz_;i1++) kpoints_bz_[i1]=(i1-G_)*dk;
	for(int i1=0;i1<nk_bz_;i1++) kweight_bz_[i1]=kw;
}

// interaction vertex V_{a1;a2}(q) = V_{a1,a2}(q) 
double lattice_1d_nofield::V(int tstp,double qq,int j1,int j2){
	if(j1!=j2){
		double xi2=xi_*xi_;
		double cq=cos(qq);
		return U_[tstp+1]+2.0*V_[tstp+1]*(xi2*cq-xi_)/(xi2-2.0*cq*xi_+1.0);
	}else{
		return 0.0;
	}
}
	
std::complex<double> lattice_1d_nofield::Velph(int tstp,std::complex<double> &X,int j1,int j2){
	if(j1!=j2){
		return g_[tstp+1]*X;
	}else{
		return 0.0;
	}
}
	
void lattice_1d_nofield::V(cdmatrix &Vmatrix,int tstp,double qq){
	// not efficient, but used only for initialization
	Vmatrix.resize(nrpa_,nrpa_);
	for(int j1=0;j1<nrpa_;j1++){
		for(int j2=0;j2<nrpa_;j2++){
			Vmatrix(j1,j2)=V(tstp,qq,j1,j2);
		}
	}
}
	
void lattice_1d_nofield::Velph(cdmatrix &Vmatrix,int tstp,std::complex<double> &X){
		// not efficient, but used only for initialization
	Vmatrix.resize(nrpa_,nrpa_);
	for(int j1=0;j1<nrpa_;j1++){
		for(int j2=0;j2<nrpa_;j2++){
			Vmatrix(j1,j2)=Velph(tstp,X,j1,j2);
		}
	}
}

// the dispersion (particle-hole symmetric for A=0)
void lattice_1d_nofield::hk(cdmatrix &hkmatrix,int tstp,double kk,int iter){
	double epsk=-2.0*tt_[tstp+1]*cos(kk);
	double epsk1=2.0*tt_[tstp+1]*cos(kk);
	double vq0_00=V(tstp,0.0,0,0);
	double vq0_01=V(tstp,0.0,0,1);
	double dbar0=-0.5*(vq0_00+vq0_01)-2.0*tt_[tstp+1]*cos(0.0);
	double dbar1=-0.5*(vq0_00+vq0_01)+2.0*tt_[tstp+1]*cos(0.0);
	double d0=dbar0-0.5*delta_[tstp+1];
	double d1=dbar1+0.5*delta_[tstp+1];
	cdouble ak=A_[tstp+1];
	hkmatrix.resize(nrpa_,nrpa_);
	hkmatrix.setZero();
	hkmatrix(0,0)=-epsk+d0;
	hkmatrix(1,1)=epsk+d1;
	hkmatrix(0,1)=ak;
	hkmatrix(1,0)=ak;

	if(tstp==-1 and 0<iter and iter<10){
		hkmatrix(0,1)=v01_[tstp+1];
		hkmatrix(1,0)=v01_[tstp+1];
	}
	// std::cout << "Herm " << hkmatrix << std::endl;
}
	
void lattice_1d_nofield::hkfree(cdmatrix &hkmatrix,int tstp,double kk){
	double epsk=-2.0*tt_[tstp+1]*cos(kk);
	double dbar=-2.0*tt_[tstp+1]*cos(0.0);
	double d0=dbar-0.5*delta_[tstp+1];
	double d1=-dbar+0.5*delta_[tstp+1];
	cdouble ak=A_[tstp+1];
	hkmatrix.resize(nrpa_,nrpa_);
	hkmatrix.setZero();
	hkmatrix(0,0)=-epsk+d0;
	hkmatrix(1,1)=epsk+d1;
	hkmatrix(0,1)=0;
	hkmatrix(1,0)=0;
}

//Velocity is still 1d vector
void lattice_1d_nofield::vk(cdmatrix &vkmatrix,int tstp,double kk){
	double vk=2.0*tt_[tstp+1]*sin(kk);
	vkmatrix.resize(nrpa_,nrpa_);
	vkmatrix.setZero();
	vkmatrix(0,0)=vk;
	vkmatrix(1,1)=-vk;
	vkmatrix(0,1)=vk;
	vkmatrix(1,0)=-vk;
}

// Dipolar interaction
void lattice_1d_nofield::Ak(cdmatrix &hkmatrix,int tstp,double kk){
	cdouble ak=A_[tstp+1];
	hkmatrix.resize(nrpa_,nrpa_);
	hkmatrix.setZero();
	hkmatrix(0,1)=conj(ak);
	hkmatrix(1,0)=ak;
}

// dAk interaction
void lattice_1d_nofield::dAk(cdmatrix &hkmatrix,int tstp,double kk){
	cdouble ak=dA_[tstp+1];
	hkmatrix.resize(nrpa_,nrpa_);
	hkmatrix.setZero();
	hkmatrix(0,1)=conj(ak);
	hkmatrix(1,0)=ak;
}	
