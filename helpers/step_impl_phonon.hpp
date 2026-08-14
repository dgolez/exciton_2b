#pragma once

#include "inclusions.hpp"

// phonontype = 0: Holstein-like coupling to total density; produces a diagonal common energy shift.
// phonontype = 1: coupling to real interband coherence / hybridization.
// phonontype = 3: coupling to the orbital-density difference \(n_1-n_0\); produces opposite diagonal shifts, like a phonon-modulated band splitting.

// DIPOLAR COUPLING
phonon::phonon(void){
}

phonon::phonon(CFUNC &omega0,CFUNC &g,CFUNC &density,int phonontype){
	omega0_=omega0;
	g_=g;
	density_=density;
  phonontype_=phonontype;
}

void phonon::step(int tstp,CFUNC &X,CFUNC &Pi,int kt,double h){
  // Fix first 5 timesteps
  // std::cout << "phonon step " << tstp << " " << kt  << std::endl;
  if(0<=tstp && tstp <=kt){
    cdmatrix tmpX,tmpPi;
    X.get_value(-1,tmpX);
    Pi.get_value(-1,tmpPi);
    // std::cout << "phonon step " << tstp << " " << tmpX << " " <<tmpPi << std::endl;
    X.set_value(tstp,tmpX);
    Pi.set_value(tstp,tmpPi);
    // std::cout << "phonon step " << tstp << " " << tmpX << " " <<tmpPi << std::endl;
  }else{
    // Adams-Moulton method
    cdmatrix  tmpPi(1,1),tmpX(1,1),tmp(1,1);
    cntr::extrapolate_timestep(tstp-1,X,integration::I<double>(kt));
    cntr::extrapolate_timestep(tstp-1,Pi,integration::I<double>(kt));
    X.get_value(tstp-1,tmp);
    tmpX(0,0)=tmp(0,0)+h*(251.0*dx1dt(tstp,Pi)+646.0*dx1dt(tstp-1,Pi)-264.0*dx1dt(tstp-2,Pi)+106.0*dx1dt(tstp-3,Pi)-19.0*dx1dt(tstp-4,Pi))/720.0;
    X.set_value(tstp,tmpX);
    Pi.get_value(tstp-1,tmp);
    tmpPi(0,0)=tmp(0,0)+h*(251.0*dx2dt(tstp,X)+646.0*dx2dt(tstp-1,X)-264.0*dx2dt(tstp-2,X)+106.0*dx2dt(tstp-3,X)-19.0*dx2dt(tstp-4,X))/720.0;
    Pi.set_value(tstp,tmpPi);
  }
}

// phonon EoM \partial P(t) = omega0*P
double phonon::dx1dt(int tstp,CFUNC &Pi){
  cdmatrix tmpom(1,1),tmpPi(1,1);
  omega0_.get_value(tstp,tmpom);
  Pi.get_value(tstp,tmpPi);
  double tmp=std::real(tmpPi(0,0)*tmpom(0,0));
  return tmp;
}

// phonon EoM \partial X(t) = -omega0*X + other, where other depends on phonontype
double phonon::dx2dt(int tstp,CFUNC &X){
  cdmatrix tmpom(1,1),tmpX(1,1),tmpG(1,1);
  int size=density_.size1_;
  cdmatrix tmprho(size,size);
  omega0_.get_value(tstp,tmpom);
  g_.get_value(tstp,tmpG);
  density_.get_value(tstp,tmprho);
  X.get_value(tstp,tmpX);
  double tmp=0.0;
  if(phonontype_==0){
    tmp=std::real(-tmpom(0,0)*tmpX(0,0));
  }else if(phonontype_==1){
    tmp=std::real(-tmpom(0,0)*tmpX(0,0)-4.0*tmpG(0,0)*tmprho(0,1).real());
  }else if(phonontype_==3){
    tmp=std::real(-tmpom(0,0)*tmpX(0,0)-2.0*tmpG(0,0)*(tmprho(1,1).real()-tmprho(0,0).real()));
  }
  return tmp;
}

void phonon::eq(CFUNC &X,CFUNC &Pi){
  int size=density_.size1_;
  cdmatrix tmpom(1,1),tmpX(1,1),tmpG(1,1);
  cdmatrix tmprho(size,size);
  cdmatrix Xtmp(1,1),Ptmp(1,1);
  omega0_.get_value(-1,tmpom);
  // std::cout << "eq 2" << std::endl;
  density_.get_value(-1,tmprho);
  // std::cout << "Phonon" << std::endl;
  g_.get_value(-1,tmpG);
  // std::cout << "Phonon " << X.size1_  <<tmpG  << std::endl;
  if(phonontype_==0){
    Xtmp(0,0)=0.0; //Compensated by phonon distortion
    Ptmp(0,0)=0.0;
  }else if(phonontype_==1){

    Xtmp(0,0)=-4.0*tmpG(0,0)*tmprho(0,1).real()/(tmpom(0,0));
    Ptmp(0,0)=0.0;

  }else if(phonontype_==3){
    Xtmp(0,0)=-2.0*tmpG(0,0)*(tmprho(1,1).real()-tmprho(0,0).real())/(tmpom(0,0));
    Ptmp(0,0)=0.0;
  }
  X.set_value(-1,Xtmp);
  Pi.set_value(-1,Ptmp);
  // std::cout << "Phonon " << Xtmp << " " << Ptmp << std::endl;
}

void phonon::hartree(cdmatrix &inter,cdmatrix &X){
  // std::cout << " Hartree Phonon inside   " << std::endl;
  cdmatrix tmpG(1,1);
  g_.get_value(-1,tmpG);
  if(phonontype_==0){
    inter(0,0)=2.0*tmpG(0,0)*X(0,0);
    inter(0,1)=0.0;
    inter(1,0)=0.0;
    inter(1,1)=2.0*tmpG(0,0)*X(0,0); 
  }else if(phonontype_==1){
    inter(0,0)=0.0;
    // inter(0,1)=2.0*tmpG(0,0)*X(0,0);
    // inter(1,0)=2.0*tmpG(0,0)*X(0,0);
    // TODO: check if 2 or without
    inter(0,1)=tmpG(0,0)*X(0,0);
    inter(1,0)=tmpG(0,0)*X(0,0);
    inter(1,1)=0.0;
    g_.get_value(-1,tmpG);
  }else if(phonontype_==3){
    inter(0,0)=-2.0*tmpG(0,0)*X(0,0);
    inter(0,1)=0.0;
    inter(1,0)=0.0;
    inter(1,1)=2.0*tmpG(0,0)*X(0,0);
  }
  // std::cout << "Hartree Phonon inside  "  <<tmpG  << " " << X << std::endl;
}

void phonon::rpa(int tstp,GREEN_TSTP &SigmaPh,cntr::herm_matrix_timestep_view<double> & G,cntr::herm_matrix_timestep_view<double> & D){
	int j1,j2;
  GREEN_TSTP tmp(tstp,SigmaPh.ntau_,1);
  if(phonontype_==0){ // Holstein
	  cntr::Bubble2(tstp,SigmaPh,0,0,G,G,0,0,D,D,0,0);
    cntr::Bubble2(tstp,SigmaPh,0,1,G,G,0,1,D,D,0,0);
    cntr::Bubble2(tstp,SigmaPh,1,0,G,G,1,0,D,D,0,0);
    cntr::Bubble2(tstp,SigmaPh,1,1,G,G,1,1,D,D,0,0);
  }else if(phonontype_==1){ // Coupling to c^{\dagger} c - real
    cntr::Bubble2(tstp,SigmaPh,0,0,G,G,1,1,D,D,0,0);
    cntr::Bubble2(tstp,SigmaPh,1,1,G,G,0,0,D,D,0,0);
    cntr::Bubble2(tstp,SigmaPh,0,1,G,G,1,0,D,D,0,0);
    cntr::Bubble2(tstp,SigmaPh,1,0,G,G,0,1,D,D,0,0);
  }else if(phonontype_==3){ // Coupling to n_1-n_0
    cntr::Bubble2(tstp,SigmaPh,0,0,G,G,0,0,D,D,0,0);
    
    cntr::Bubble2(tstp,tmp,0,0,G,G,0,1,D,D,0,0);
    tmp.smul(tstp,-1.0);
    SigmaPh.set_matrixelement(0,1,tmp,0,0);
    
    cntr::Bubble2(tstp,tmp,0,0,G,G,1,0,D,D,0,0);
    tmp.smul(tstp,-1.0);
    SigmaPh.set_matrixelement(1,0,tmp,0,0);
    
    cntr::Bubble2(tstp,SigmaPh,1,1,G,G,1,1,D,D,0,0);
  }
}

// Do this more generally
void phonon::matrixele(int i1,int i2,int &j1, int &j2,int nrpa){
	// j1=(i1+1)%nrpa;
	// j2=(i2+1)%nrpa;
	// std::cout << "Inside " << i1 <<" " << i2 <<" " << j1 << " " <<j2 <<std::endl;
}