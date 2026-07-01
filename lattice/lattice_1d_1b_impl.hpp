#pragma once

#include "../step/step.hpp"
lattice_1d_1b_nofield::lattice_1d_1b_nofield(void){
}

lattice_1d_1b_nofield::lattice_1d_1b_nofield(int nk,int nt,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &g,CFUNC &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,double mu,int size,int kt,double h):
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
  E_.resize(nt+2);
  dA_.resize(nt+2);
  v01_.resize(nt+2);
  delta_.resize(nt+2);
  mu_=mu;
  
  for(int tstp=-1;tstp<=nt;tstp++){
    tt_[tstp+1]=tt[tstp+1];
    U_[tstp+1]=U[tstp+1];
    V_[tstp+1]=V[tstp+1];
    E_[tstp+1]=A[tstp+1]; // A is used as the electric field - fix if not nice
    dA_[tstp+1]=dA[tstp+1];
    v01_[tstp+1]=0.0;
    delta_[tstp+1]=delta[tstp+1];
  }
  efield_to_afield(nt,h,E_,A_,kt);
  v01_[0]=v01; //Artificial symmetry breaking field only for equilibrium
  init_kk(nk);
}

// A(t)=-int_0^t ds E(s) afield[0]=0, afield[i]=A(h*(tstp-1)], i=1...nt, same or efield
void lattice_1d_1b_nofield::efield_to_afield(int nt,double h,std::vector<double> &efield,std::vector<double> &afield,int kt){
  int kt1=(nt>=kt ? kt : nt),n,n1,tstp;
  double At;
  assert((int)efield.size()==nt+2);
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


int lattice_1d_1b_nofield::representative_kk(int kbz){
  assert(kbz<=nk_bz_-1);	
  return ( kbz <= G_ ? kbz  : 2*G_-kbz);
}

int lattice_1d_1b_nofield::add_kpoints(int k1,int s1,int k2,int s2){
  int k12=s1*(k1-G_)+s2*(k2-G_)+G_;
  while (k12<0){k12 += nk_bz_;}
  while (k12>=nk_bz_){k12 -= nk_bz_;}
  return k12;
}

void lattice_1d_1b_nofield::init_kk(int nk){
  double dk,kw;
  assert(nk>1);
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
double lattice_1d_1b_nofield::V(int tstp,double qq,int j1,int j2){
  if(j1!=j2){
    double xi2=xi_*xi_;
    double cq=cos(qq);
    return U_[tstp+1]+2.0*V_[tstp+1]*(xi2*cq-xi_)/(xi2-2.0*cq*xi_+1.0);
  }else{
    return 0.0;
  }
}
	
std::complex<double> lattice_1d_1b_nofield::Velph(int tstp,std::complex<double> &X,int j1,int j2){
  if(j1!=j2){
    return g_[tstp+1]*X;
  }else{
    return 0.0;
  }
}
	
void lattice_1d_1b_nofield::V(cdmatrix &Vmatrix,int tstp,double qq){
	// not efficient, but used only for initialization
  Vmatrix.resize(nrpa_,nrpa_);
  for(int j1=0;j1<nrpa_;j1++){
    for(int j2=0;j2<nrpa_;j2++){
      Vmatrix(j1,j2)=V(tstp,qq,j1,j2);
    }
  }
}
	
void lattice_1d_1b_nofield::Velph(cdmatrix &Vmatrix,int tstp,std::complex<double> &X){
  // not efficient, but used only for initialization
  Vmatrix.resize(nrpa_,nrpa_);
  for(int j1=0;j1<nrpa_;j1++){
    for(int j2=0;j2<nrpa_;j2++){
      Vmatrix(j1,j2)=Velph(tstp,X,j1,j2);
    }
  }
}

// the dispersion (particle-hole symmetric for A=0)
void lattice_1d_1b_nofield::hk(cdmatrix &hkmatrix,int tstp,double kk,int iter){
  double epsk=-2.0*tt_[tstp+1]*cos(kk-A_[tstp+1]);
  hkmatrix.resize(nrpa_,nrpa_);
  hkmatrix.setZero();
  hkmatrix(0,0)=epsk;
}
	
void lattice_1d_1b_nofield::hkfree(cdmatrix &hkmatrix,int tstp,double kk){
	double epsk=-2.0*tt_[tstp+1]*cos(kk);
	hkmatrix.resize(nrpa_,nrpa_);
	hkmatrix.setZero();
	hkmatrix(0,0)=epsk;
}

//Velocity is still 1d vector
void lattice_1d_1b_nofield::vk(cdmatrix &vkmatrix,int tstp,double kk){
	double vk=2.0*tt_[tstp+1]*sin(kk-A_[tstp+1]);
	vkmatrix.resize(nrpa_,nrpa_);
	vkmatrix.setZero();
	vkmatrix(0,0)=vk;
}

// Peierls phase
void lattice_1d_1b_nofield::Ak(cdmatrix &hkmatrix,int tstp,double kk){
	cdouble ak=A_[tstp+1];
	hkmatrix.resize(nrpa_,nrpa_);
	hkmatrix.setZero();
	// hkmatrix(0,1)=conj(ak);
	// hkmatrix(1,0)=ak;
}

// dAk interaction
void lattice_1d_1b_nofield::dAk(cdmatrix &hkmatrix,int tstp,double kk){
	cdouble ak=dA_[tstp+1];
	hkmatrix.resize(nrpa_,nrpa_);
	hkmatrix.setZero();
	// hkmatrix(0,1)=conj(ak);
	// hkmatrix(1,0)=ak;
}	
