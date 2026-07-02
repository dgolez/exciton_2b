#pragma once

#include "../step/step.hpp"
lattice_1d_rpa_cdw_nofield::lattice_1d_rpa_cdw_nofield(void){
}

lattice_1d_rpa_cdw_nofield::lattice_1d_rpa_cdw_nofield(int nk,int nt,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &g,CFUNC &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,double mu,int size,int kt,double h):
	nt_(nt),nrpa_(size),xi_(xi),mu_(0.0)
{
  assert(-1<=nt);
  assert(1.0<=xi);
  assert(-1<=nt);
  assert(nt+2==(int)U.size());

	tt_.resize(nt+2);
	U_.resize(nt+2);
	V_.resize(nt+2);
	omega0_=omega0;
	g_=g;
	A_.resize(nt+2);
	dA_.resize(nt+2);
	v01_.resize(nt+2);
	delta_.resize(nt+2);
	norb_=2;
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

int lattice_1d_rpa_cdw_nofield::representative_kk(int krbz){
		// CNTR_ASSERT_LESEQ_3(ASSERT_0,0,krbz,nk_rbz_-1,__PRETTY_FUNCTION__)	
		return ( krbz <= G_ ? krbz  : 2*G_-krbz);
}

int lattice_1d_rpa_cdw_nofield::k_fbz(int krbz){
		// CNTR_ASSERT_LESEQ_3(ASSERT_0,0,krbz,nk_rbz_-1,__PRETTY_FUNCTION__)	
		return krbz+G_;
}

void lattice_1d_rpa_cdw_nofield::decompose_kfbz(int kfbz,int &krbz,int &gamma){
		// CNTR_ASSERT_LESEQ_3(ASSERT_0,0,kfbz,nk_fbz_-1,__PRETTY_FUNCTION__)	
		if(kfbz<G_){krbz=kfbz+G_;gamma=1;}
		else if(kfbz<3*G_){krbz=kfbz-G_;gamma=0;}
		else{krbz=kfbz-3*G_;gamma=1;}
}



int lattice_1d_rpa_cdw_nofield::add_kpoints_fbz(int s1,int k1,int s2,int k2){
		int k12=s1*(k1-2*G_)+s2*(k2-2*G_)+2*G_;
		while (k12<0){k12 += nk_fbz_;}
		while (k12>=nk_fbz_){k12 -= nk_fbz_;}
		return k12;
}

int lattice_1d_rpa_cdw_nofield::add_kpoints(int s1,int k1,int s2,int k2){
		int k12=s1*(k1-2*G_)+s2*(k2-2*G_)+2*G_;
		while (k12<0){k12 += nk_fbz_;}
		while (k12>=nk_fbz_){k12 -= nk_fbz_;}
		return k12;
}

// add and project back to RBZ
int lattice_1d_rpa_cdw_nofield::add_kpoints_rbz(int &k12, int &gamma,int s1,int k1,int s2,int k2){
		int k3=add_kpoints_fbz(s1,k_fbz(k1),s2,k_fbz(k2));
		decompose_kfbz(k3,k12,gamma);
}



void lattice_1d_rpa_cdw_nofield::init_kk(int nk){
	double dk,kw;
	assert(nk>1);
	nk_=nk;
	G_=nk_-1;
	nk_rbz_=2*G_;
	nk_fbz_=4*G_;
	dk=PI/nk_rbz_; // kk= (index-G_)*dk
	kw=1.0/nk_rbz_;
	kpoints_.resize(nk_);
	kpoints_rbz_.resize(nk_rbz_);
	kweight_bz_.resize(nk_rbz_); //Actually weights over reduced bz 
	kpoints_fbz_.resize(nk_fbz_);
	for(int i1=0;i1<nk_;i1++) kpoints_[i1]=(i1-G_)*dk;
	for(int  i1=0;i1<nk_rbz_;i1++) kpoints_rbz_[i1]=(i1-G_)*dk;
	for(int i1=0;i1<nk_rbz_;i1++) kweight_bz_[i1]=kw;
	for(int i1=0;i1<nk_fbz_;i1++) kpoints_fbz_[i1]=(i1-2*G_)*dk;
	Q_=PI;
}

int lattice_1d_rpa_cdw_nofield::idx(int orb,int sl){return (orb%2)*2+(sl%2);}
int lattice_1d_rpa_cdw_nofield::idx_a(int j){ return (j%4)/2;}
int lattice_1d_rpa_cdw_nofield::idx_c(int j){ return  (j%4)%2;}



double lattice_1d_rpa_cdw_nofield::V(int tstp,double qq,int j1,int j2){
	int a1=idx_a(j1);
	int c1=idx_c(j1);
	int a2=idx_a(j2);
	int c2=idx_c(j2);
	if(c1==c2){
		double xi2=xi_*xi_;
		double cq=cos(qq+c1*Q_);
		return U_[tstp+1]*(a1!=a2 ? 1 : 0)+2.0*V_[tstp+1]*(xi2*cq-xi_)/(xi2-2.0*cq*xi_+1.0);
	}else{
		return 0.0;
	}
}

void lattice_1d_rpa_cdw_nofield::V(cdmatrix &Vmatrix,int tstp,double qq){
	// not efficient, but used only for initialization
	Vmatrix.resize(nrpa_,nrpa_);
	for(int j1=0;j1<nrpa_;j1++){
		for(int j2=0;j2<nrpa_;j2++){
			Vmatrix(j1,j2)=V(tstp,qq,j1,j2);
		}
	}
}


	
std::complex<double> lattice_1d_rpa_cdw_nofield::Velph(int tstp,std::complex<double> &X,int j1,int j2){
	if(j1!=j2){
		return g_[tstp+1]*X;
	}else{
		return 0.0;
	}
}
	
void lattice_1d_rpa_cdw_nofield::Velph(cdmatrix &Vmatrix,int tstp,std::complex<double> &X){
		// not efficient, but used only for initialization
	Vmatrix.resize(nrpa_,nrpa_);
	for(int j1=0;j1<nrpa_;j1++){
		for(int j2=0;j2<nrpa_;j2++){
			Vmatrix(j1,j2)=Velph(tstp,X,j1,j2);
		}
	}
}

// the dispersion (particle-hole symmetric for A=0)

void lattice_1d_rpa_cdw_nofield::hk(cdmatrix &hkmatrix,int tstp,double kk,int iter=-1){
		double epsk=-2.0*tt_[tstp+1]*cos(kk);
		double vq0_00=V(tstp,0.0,idx(0,0),idx(0,0));
		double vq0_01=V(tstp,0.0,idx(0,0),idx(1,0));
		double dbar=-0.5*(vq0_00+vq0_01);
		double d0=dbar-0.5*delta_[tstp+1];
		double d1=dbar+0.5*delta_[tstp+1];
		// std::cout << "param " <<  kk << " " << epsk << " " << vq0_00 << " " << vq0_01 << " " << dbar << " " <<  d0 << " " << d1 << std::endl;
		cdouble ak=A_[tstp+1];
		hkmatrix.resize(nrpa_,nrpa_);
		hkmatrix.setZero();
		hkmatrix(0,0)=epsk+d0;
		hkmatrix(1,1)=-epsk+d0; // eps(k+Q)=-eps(k)
		hkmatrix(2,2)=epsk+d1;
		hkmatrix(3,3)=-epsk+d1;
		
		hkmatrix(0,2)=conj(ak);
		hkmatrix(2,0)=ak;
		hkmatrix(1,3)=conj(ak);
		hkmatrix(3,1)=ak;

		if(tstp==-1 and 0<iter and iter<15){
			hkmatrix(3,0)=v01_[tstp+1];
			hkmatrix(2,1)=v01_[tstp+1];
			hkmatrix(1,2)=v01_[tstp+1];
			hkmatrix(0,3)=v01_[tstp+1];
		}
}

void lattice_1d_rpa_cdw_nofield::hkfree(cdmatrix &hkmatrix,int tstp,double kk){
		double epsk=-2.0*tt_[tstp+1]*cos(kk);
		double d0=-0.5*delta_[tstp+1];
		double d1=+0.5*delta_[tstp+1];
		cdouble ak=A_[tstp+1];
		hkmatrix.resize(nrpa_,nrpa_);
		hkmatrix.setZero();
		hkmatrix(0,0)=epsk+d0;
		hkmatrix(1,1)=-epsk+d0; // eps(k+Q)=-eps(k)
		hkmatrix(2,2)=epsk+d1;
		hkmatrix(3,3)=-epsk+d1;
		hkmatrix(0,2)=conj(ak);
		hkmatrix(2,0)=ak;
		hkmatrix(1,3)=conj(ak);
		hkmatrix(3,1)=ak;
	}

void lattice_1d_rpa_cdw_nofield::hkkin(cdmatrix &hkmatrix,int tstp,double kk){
		double epsk=-2.0*tt_[tstp+1]*cos(kk);
		double d0=-0.5*delta_[tstp+1];
		double d1=+0.5*delta_[tstp+1];
		cdouble ak=A_[tstp+1];
		hkmatrix.resize(nrpa_,nrpa_);
		hkmatrix.setZero();
		hkmatrix(0,0)=epsk+d0;
		hkmatrix(1,1)=-epsk+d0; // eps(k+Q)=-eps(k)
		hkmatrix(2,2)=epsk+d1;
		hkmatrix(3,3)=-epsk+d1;
		hkmatrix(0,2)=conj(ak);
		hkmatrix(2,0)=ak;
		hkmatrix(1,3)=conj(ak);
		hkmatrix(3,1)=ak;
}
	

//Velocity is still 1d vector
void lattice_1d_rpa_cdw_nofield::vk(cdmatrix &vkmatrix,int tstp,double kk){
		double vk=2.0*tt_[tstp+1]*sin(kk);
		vkmatrix.resize(nrpa_,nrpa_);
		vkmatrix.setZero();
		vkmatrix(0,0)=vk;
		vkmatrix(1,1)=-vk; // eps(k+Q)=-eps(k)
		vkmatrix(2,2)=vk;
		vkmatrix(3,3)=-vk;
	}

// Dipolar interaction
void lattice_1d_rpa_cdw_nofield::Ak(cdmatrix &hkmatrix,int tstp,double kk){
		cdouble ak=A_[tstp+1];
		hkmatrix.resize(nrpa_,nrpa_);
		hkmatrix.setZero();
		hkmatrix(0,2)=conj(ak);
		hkmatrix(2,0)=ak;
		hkmatrix(1,3)=conj(ak);
		hkmatrix(3,1)=ak;
	}

// dAk interaction
void lattice_1d_rpa_cdw_nofield::dAk(cdmatrix &hkmatrix,int tstp,double kk){
		cdouble ak=dA_[tstp+1];
		hkmatrix.resize(nrpa_,nrpa_);
		hkmatrix.setZero();
		hkmatrix(0,2)=conj(ak);
		hkmatrix(2,0)=ak;
		hkmatrix(1,3)=conj(ak);
		hkmatrix(3,1)=ak;
	}	
