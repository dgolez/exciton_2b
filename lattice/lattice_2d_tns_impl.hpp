#pragma once

#include "../step/step.hpp"
#include "../program/parameters.hpp"

lattice_2d_tns::lattice_2d_tns(void){
}

lattice_2d_tns::lattice_2d_tns(parameters_tns &param):
	nt_(param.nt),nrpa_(param.size)
{
	int nt=param.nt;
	assert(-1<=nt);
  assert(-1<=param.nt);
  assert(nt+2==(int)param.U.size());

  tt_.resize(nt+2);
  U_.resize(nt+2);
  V_.resize(nt+2);
  Ax_.resize(nt+2);
  Ay_.resize(nt+2);
  Ex_.resize(nt+2);
  Ey_.resize(nt+2);
  nkx_=param.nk;

	v01_.resize(nt+2);
	dipol_=param.dipol;
	mu_=param.mu;
	ratio_=param.ratio;
	rho_eq_.resize(nrpa_,nrpa_);
	
	for(int tstp=-1;tstp<=nt;tstp++){
		tt_[tstp+1]=param.tt[tstp+1];
		U_[tstp+1]=param.U[tstp+1];
		V_[tstp+1]=param.V[tstp+1];
		Ex_[tstp+1]=param.Ex[tstp+1];
		Ey_[tstp+1]=param.Ey[tstp+1];
		Ex_[tstp+1]=param.Ex[tstp+1];
		Ey_[tstp+1]=param.Ey[tstp+1];
		v01_[tstp+1]=0.0;
	}
	efield_to_afield(nt,param.h,param.Ex,Ax_,param.kt);
	efield_to_afield(nt,param.h,param.Ey,Ay_,param.kt);

	v01_[0]=param.v01; //Artificial symmetry breaking field only for equilibrium
	
	init_kk(nkx_);
}

void lattice_2d_tns::efield_to_afield(int nt,double h,std::vector<double> &E,std::vector<double> &A,int kt){
	int kt1=(nt>=kt ? kt : nt),n,n1,tstp;
	double At;
	assert((int)E.size()>=nt+2);
	A.resize(nt+2);
	A[0]=0.0; // tstp=-1
	for(tstp=0;tstp<=nt;tstp++){
		At=0.0;
		n1=(tstp<kt1 ? kt1 : tstp);
		for(n=0;n<=n1;n++){
			At += integration::I<double>(kt1).gregory_weights(tstp,n)*E[n+1];
		}
		A[tstp+1]=At*(-h);
	}
}

//Sum vectors in bz
int lattice_2d_tns::add_kpoints(int k1,int s1,int k2,int s2){
  div_t divk1, divk2;
  int nk;

  divk1=div(k1,nkx_);
  divk2=div(k2,nkx_);

  int k12x=s1*(divk1.rem-G_)+s2*(divk2.rem-G_)+G_;
  int k12y=s1*(divk1.quot-G_)+s2*(divk2.quot-G_)+G_;
  int k12;
  while (k12x<0){k12x += nkx_;}
  while (k12y<0){k12y += nkx_;}
  while (k12x>=nk){k12x -= nkx_;}
  while (k12y>=nk){k12y -= nkx_;}
  k12 = k12x+nkx_*k12y;
  return k12;
}

void lattice_2d_tns::init_kk(int nkx){
  int ix=0;
  double dk,kw;
  dvector kk_;
  assert(nkx>1);
  nk_=nkx*nkx;	// square lattice
  //nk_=nk;
  G_=nkx/2;
  dk=2*PI/nkx; // kk= (index-G_)*dk 
  kw=1.0/nk_;
  kk_=dvector(2);
  kpoints_.resize(nk_);
  kweight_.resize(nk_);	
  for(int i1=0;i1<nkx;i1++){
    for(int i2=0;i2<nkx;i2++){
      kk_(0)=(i1-G_)*dk;
      kk_(1)=(i2-G_)*dk;
//      kk_(0)=(i1-G_)*dk/(std::sqrt(2))+(i2-G_)*dk/(std::sqrt(2));
//      kk_(1)=(i1-G_)*dk/(std::sqrt(2))-(i2-G_)*dk/(std::sqrt(2));
      kpoints_[ix]=kk_;
//	  std::cout<< "kx=" << kk_(0) << " ky="<<kk_(1) << std::endl;
	  ix++;
    }
  }
  for(int i1=0;i1<nk_;i1++) kweight_[i1]=kw;
}


// interaction vertex V_{a1;a2}(q) = V_{a1,a2}(q) 
std::complex<double> lattice_2d_tns::V(int tstp,dvector &kk,int j1,int j2){
	std::complex<double> I(0.0,1.0);

	if(j1==j2){
		return U_[tstp+1]/2.0;
	}else if(j1==0 && j2==4){ // V_1 entries
		return V_[tstp]*(1.0+exp(I*kk));
	}else if(j1==1 && j2==4){
		return V_[tstp]*(1.0+exp(I*kk));
	}else if(j1==4 && j2==0){
		return V_[tstp]*(1.0+exp(-I*kk));
	}else if(j1==4 && j2==1){
		return V_[tstp]*(1.0+exp(-I*kk));
	}else if(j1==2 && j2==5){ // V_2 entries
		return V_[tstp]*(1.0+exp(-I*kk));
	}else if(j1==3 && j2==5){
		return V_[tstp]*(1.0+exp(-I*kk));
	}else if(j1==5 && j2==2){ // V_2 entries
		return V_[tstp]*(1.0+exp(I*kk));
	}else if(j1==5 && j2==3){
		return V_[tstp]*(1.0+exp(I*kk));
	}else{
		return 0;
	}
}
	
void lattice_2d_tns::V(cdmatrix &Vmatrix,int tstp,dvector &qq){
	// not efficient, but used only for initialization
	Vmatrix.resize(nrpa_,nrpa_);
	for(int j1=0;j1<nrpa_;j1++){
		for(int j2=0;j2<nrpa_;j2++){
			Vmatrix(j1,j2)=V(tstp,qq,j1,j2);
		}
	}
}

// TODO
void lattice_2d_tns::Velph(cdmatrix &Vmatrix,int tstp,std::complex<double> &X){
	// Vmatrix.resize(nrpa_,nrpa_);
	// for(int j1=0;j1<nrpa_;j1++){
	// 	for(int j2=0;j2<nrpa_;j2++){
	// 		Vmatrix(j1,j2)=Velph(tstp,X,j1,j2);
	// 	}
	// }
}


void lattice_2d_tns::hk(cdmatrix &hkmatrix,int tstp,dvector &kk,int iter){
	hkmatrix.resize(nrpa_,nrpa_);
	hkmatrix.setZero();
	
	hkfree(hk,tstp,kk);

	// TODO: Subtract hartree (double counting)

	// Set seed
	if(tstp==-1 and 0<iter and iter<10){
		hkmatrix(0,1)=v01_[tstp+1];
		hkmatrix(1,0)=conj(v01_[tstp+1]);
	}
}
	
void lattice_2d_tns::hkfree(cdmatrix &hkmatrix,int tstp,dvector &kk){
	double kkshiftx=kk(0)+Ax_[tstp+1];
	double kkshifty=kk(1)+Ay_[tstp+1];
	double norm=0.3; // All energies are normalized to the hopping of the valence band [0.3]
	double tNi=1.0;
	double tTa1=-0.72/norm;
	double epsTa=1.35/norm;
	double epsNi=-0.36/0.3;
	double t15=0.035/norm;
	double t45=0.04/norm;
	double t65=0.03/norm;
	double t23=0.02/norm;


	// Determine dispersions
	std::complex<double> I(0.0,1.0);

	std::complex<double> epskTa=2.0*tTa*cos(kkshiftx)+epsTa;
	std::complex<double> epskNi=2.0*tNi*cos(kkshiftx)+epsNi;
	std::complex<double> epsk15=2*I*t15*exp(-I*kkshiftx/2.0)*sin(kkshiftx/2.0);
	std::complex<double> epsk25=epsk15;
	std::complex<double> eps36=-eps25;
	std::complex<double> eps46=-eps25;
	std::complex<double> eps45=2*I*t45*exp(I*kkshifty)*sin(kkshiftx);
	std::complex<double> eps35=2*I*t45*sin(kkshiftx);
	std::complex<double> eps26=-eps35;
	std::complex<double> eps16=-2*I*t45*exp(-I*kkshifty)*sin(kkshiftx);
	std::complex<double> eps56=t65*(1+exp(I*kkshifty))*(1+exp(kkshiftx));
	std::complex<double> eps23=t23*(1+exp(kkshiftx));

	hkmatrix.resize(nrpa_,nrpa_);
	hkmatrix.setZero();

	// diagonal
	hkmatrix(0,0)=epskTa;
	hkmatrix(1,1)=epskTa;
	hkmatrix(2,2)=epskTa;
	hkmatrix(3,3)=epskTa;
	hkmatrix(4,4)=epskNi;
	hkmatrix(5,5)=epskNi;

	// 
	hkmatrix(0,4)=epsk15;
	hkmatrix(4,0)=conj(epsk15);
	hkmatrix(1,4)=epsk15;
	hkmatrix(4,1)=conj(epsk15);
	hkmatrix(2,5)=-epsk15;
	hkmatrix(5,2)=-conj(epsk15);
	hkmatrix(3,5)=-epsk15;
	hkmatrix(5,3)=-conj(epsk15);

	hkmatrix(3,4)=epsk45;
	hkmatrix(4,3)=conj(epsk45);

	hkmatrix(2,4)=epsk35;
	hkmatrix(4,2)=conj(epsk35);

	hkmatrix(0,5)=epsk16;
	hkmatrix(5,0)=conj(epsk16);

	hkmatrix(4,5)=epsk56;
	hkmatrix(5,4)=conj(epsk56);

	hkmatrix(1,2)=epsk23;
	hkmatrix(2,1)=conj(epsk23);

}

//Velocity is still 1d vector
void lattice_2d_tns::vk(cdmatrix &vkmatrix,int tstp,double kk){
	// Diagonal terms
	double kkshift=kk+A_[tstp+1];
	double tTa=-0.72/0.3;
	double tNi=1.0;

	double vkTa=-2.0*tTa*sin(kkshift);
	double vkNi=-2.0*tNi*sin(kkshift);
	std::complex<double> I(0.0,1.0);
	std::complex<double> vkMaz=2.0*sqrt(2)*I*(mazza_/0.3)*exp((I*kk)/2.0);

	// double vk=2.0*tt_[tstp+1]*sin(kk);
	vkmatrix.resize(nrpa_,nrpa_);
	vkmatrix.setZero();
	vkmatrix(0,0)=vkTa;
	vkmatrix(1,1)=vkNi;
	vkmatrix(0,1)=vkMaz;
	vkmatrix(1,0)=conj(vkMaz);
}

// the dispersion (particle-hole symmetric for A=0)
// void lattice_2d_tns::hk(cdmatrix &hkmatrix,int tstp,double kk,int iter){
// 	// All terms normalized to the hopping of the valence band [0.3]
// 	double kkshift=kk+A_[tstp+1]*fieldP_;
// 	double tTa=-1.0;
// 	double tNi=1.0;
// 	double epsTa=1.0;
// 	double epsNi=-1.0;
// 	double epskTa=2.0*tTa*cos(kkshift);
// 	double epskNi=2.0*tNi*cos(kkshift);
// 	// Off-diagonal components (assumption of a/2 distance between atoms)
// 	cdouble ak=dipol_*E_[tstp+1]*fieldD_;

// 	// double vq0_01=V(tstp,0.0,0,1);
// 	double d0=-V(tstp,0.0,0,0)*std::real(rho_eq_(0,0))-V(tstp,0.0,0,1)*std::real(rho_eq_(1,1))-0.5*delta_[tstp+1];
// 	double d1=-V(tstp,0.0,1,1)*std::real(rho_eq_(1,1))-V(tstp,0.0,0,1)*std::real(rho_eq_(0,0))+0.5*delta_[tstp+1];
// 	std::complex<double> I(0.0,1.0);
// 	std::complex<double> maz=2.0*sqrt(2)*I*((mazza_/0.3)*exp(I*kk/2.0)*sin(kk/2.0));

// 	hkmatrix.resize(nrpa_,nrpa_);
// 	hkmatrix.setZero();
// 	hkmatrix(0,0)=epskTa+epsTa+d0;
// 	hkmatrix(1,1)=epskNi+epsNi+d1;
// 	hkmatrix(0,1)=ak+maz;
// 	hkmatrix(1,0)=conj(ak+maz);

// 	if(tstp==-1 and 0<iter and iter<10){
// 		hkmatrix(0,1)=v01_[tstp+1]+maz;
// 		hkmatrix(1,0)=conj(v01_[tstp+1]+maz);
// 	}
// }
	
// void lattice_2d_tns::hkfree(cdmatrix &hkmatrix,int tstp,double kk){
// 	// All terms normalized to the hopping of the valence band [0.36]
// 	double tTa=-1.0;
// 	double tNi=1.0;
// 	double epsTa=1.0;
// 	double epsNi=-1.0;
// 	double epskTa=2.0*tTa*cos(kk);
// 	double epskNi=2.0*tNi*cos(kk);
// 	double d0=-0.5*delta_[tstp+1];
// 	double d1=0.5*delta_[tstp+1];

// 	// Off-diagonal components (assumption of a/2 distance between atoms)
// 	cdouble ak=dipol_*E_[tstp+1]*fieldD_;
// 	std::complex<double> I(0.0,1.0);
// 	std::complex<double> maz=2.0*sqrt(2)*I*(mazza_/0.3)*exp(I*kk/2.0)*sin(kk/2.0);
	
// 	hkmatrix.resize(nrpa_,nrpa_);
// 	hkmatrix.setZero();
// 	hkmatrix(0,0)=epskTa+epsTa+d0;
// 	hkmatrix(1,1)=epskNi+epsNi+d1;
// 	hkmatrix(0,1)=maz;
// 	hkmatrix(1,0)=conj(maz);

// }

// //Velocity is still 1d vector
// void lattice_2d_tns::vk(cdmatrix &vkmatrix,int tstp,double kk){
// 	// Diagonal terms
// 	double kkshift=kk+A_[tstp+1];
// 	double tTa=-0.72/0.3;
// 	double tNi=1.0;

// 	double vkTa=-2.0*tTa*sin(kkshift);
// 	double vkNi=-2.0*tNi*sin(kkshift);
// 	std::complex<double> I(0.0,1.0);
// 	std::complex<double> vkMaz=2.0*sqrt(2)*I*(mazza_/0.3)*exp((I*kk)/2.0);

// 	// double vk=2.0*tt_[tstp+1]*sin(kk);
// 	vkmatrix.resize(nrpa_,nrpa_);
// 	vkmatrix.setZero();
// 	vkmatrix(0,0)=vkTa;
// 	vkmatrix(1,1)=vkNi;
// 	vkmatrix(0,1)=vkMaz;
// 	vkmatrix(1,0)=conj(vkMaz);
// }

//Velocity is still 1d vector
void lattice_2d_tns::vkFULL(cdmatrix &vkmatrix,int tstp,double kk){
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
void lattice_2d_tns::Ak(cdmatrix &hkmatrix,int tstp,double kk){
        double kkshift=kk;
	cdouble ak=dipol_;
	// cdouble ak=I*dipol_*E_[tstp+1]*std::cos(kkshift/2.0);

	hkmatrix.resize(nrpa_,nrpa_);
	hkmatrix.setZero();
	hkmatrix(0,1)=conj(ak);
	hkmatrix(1,0)=ak;
}

// dAk interaction
void lattice_2d_tns::dAk(cdmatrix &hkmatrix,int tstp,double kk){
	// cdouble ak=dA_[tstp+1];
	// double kkshift=kk+A_[tstp+1];
	// cdouble ak=s*dipol_*dE_[tstp+1]*std::cos(kkshift/2.0);
	// hkmatrix.resize(nrpa_,nrpa_);
	// hkmatrix.setZero();
	// hkmatrix(0,1)=conj(ak);
	// hkmatrix(1,0)=ak;
}	
