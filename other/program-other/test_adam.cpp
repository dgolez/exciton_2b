#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include <memory>
#include <iostream>
#include <fstream>
#define CNTR_USE_OMP
#define CNTR_USE_MPI
#include "cntr/cntr.hpp"
#include "cntr/utils/read_inputfile.hpp"
#include "../step/step.hpp"
#include "print.hpp"
#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision

using namespace std; 


double dx1dt(int tstp,CFUNC &Pi,double omega0){
  cdmatrix tmpPi(1,1);
  Pi.get_value(tstp,tmpPi);
  return std::real(tmpPi(0.0)*omega0);
}

double dx2dt(int tstp,CFUNC &X,CFUNC &rho,double omega0,double g){
  cdmatrix tmpom(1,1),tmpX(1,1),tmpG(1,1);
  cdmatrix tmprho(1,1);
  rho.get_value(tstp,tmprho);
  X.get_value(tstp,tmpX);
  return std::real(-omega0*tmpX(0,0)-4.0*g*tmprho(0,0).real());
}

int main(int argc,char *argv[]){
  int kt=5;
  double omega0=1.0;
  double tmax=20;
  double g=0.5;
  double dt=0.01;
  //find_param(argv[1],"__dt=",dt);
  int nt=int(tmax/dt);
  CFUNC X(nt,1),Pi(nt,1),rho(nt,1),Xext(nt,1);

  ofstream myfile;
  myfile.open ("phonon01.dat");

  ofstream rho1;
  rho1.open ("rho.dat");
  
  //The driving "force"
  for(int tstp=-1;tstp<=nt;tstp++){
    cdmatrix tmprho(1,1);
    if(tstp<=kt){
      tmprho(0,0)=0.0;
    }else{
      tmprho(0,0)=(tstp-5)*dt;
    }
    rho1 << tstp*dt << " " << tmprho(0,0).real() << std::endl;
    rho.set_value(tstp,tmprho);
  }
  rho1.close();

  for(int tstp=-1;tstp<=kt;tstp++){
    cdmatrix tmpX(1,1),tmpPi(1,1),tmprho(1,1);
    rho.get_value(tstp,tmprho);
    tmpX(0,0)=-4.0*g*tmprho(0,0)/omega0;
    tmpPi(0,0)=0.0;
    // tmpX(0,0)=-4.0*g*tmprho(0,0)/omega0;
    // tmpX(0,0)=cos(tstp*dt*omega0);
    // tmpPi(0,0)=-omega0*sin(tstp*dt*omega0);
    X.set_value(tstp,tmpX);
    Pi.set_value(tstp,tmpPi);
  }

  for(int tstp=kt+1;tstp<=nt;tstp++){
    for(int iter=0;iter<2;iter++){
      //Adams-Moulton method
      cdmatrix  tmpPi(1,1),tmpX(1,1),tmp(1,1);
      if(iter==0){
	       // cntr::extrapolate_timestep(tstp-1,X,integration::I<double>(5));
	       // cntr::extrapolate_timestep(tstp-1,Pi,integration::I<double>(5));
  	       X.get_value(tstp-1,tmp);
  	       tmpX(0,0)=tmp(0,0)+dt*(23.0*dx1dt(tstp-1,Pi,omega0)/12.0-4.0*dx1dt(tstp-2,Pi,omega0)/3.0+5.0*dx1dt(tstp-3,Pi,omega0)/12.0);
  	       X.set_value(tstp,tmpX);
  	       Pi.get_value(tstp-1,tmp);
  	       tmpPi(0,0)=tmp(0,0)+dt*(23.0*dx2dt(tstp-1,X,rho,omega0,g)/12.0-4.0*dx2dt(tstp-2,X,rho,omega0,g)/3.0+5.0*dx2dt(tstp-3,X,rho,omega0,g)/12.0);
  	       Pi.set_value(tstp,tmpPi);
      }
      X.get_value(tstp-1,tmp);
      tmpX(0,0)=tmp(0,0)+dt*(251.0*dx1dt(tstp,Pi,omega0)+646.0*dx1dt(tstp-1,Pi,omega0)-264.0*dx1dt(tstp-2,Pi,omega0)+106.0*dx1dt(tstp-3,Pi,omega0)-19.0*dx1dt(tstp-4,Pi,omega0))/720.0;
      // tmpX(0,0)=tmp(0,0)+dt*dx1dt(tstp,Pi,omega0);
      X.set_value(tstp,tmpX);
      Pi.get_value(tstp-1,tmp);
      tmpPi(0,0)=tmp(0,0)+dt*(251.0*dx2dt(tstp,X,rho,omega0,g)+646.0*dx2dt(tstp-1,X,rho,omega0,g)-264.0*dx2dt(tstp-2,X,rho,omega0,g)+106.0*dx2dt(tstp-3,X,rho,omega0,g)-19.0*dx2dt(tstp-4,X,rho,omega0,g))/720.0;
      // tmpPi(0,0)=tmp(0,0)+dt*dx2dt(tstp,X,rho,omega0,g);
      Pi.set_value(tstp,tmpPi);
    }
  }
  
  for(int tstp=0;tstp<=nt;tstp++){
    cdmatrix  tmpPi(1,1),tmpX(1,1);
    X.get_value(tstp,tmpX);
    Pi.get_value(tstp,tmpPi);
    myfile << std::setprecision(16) << tstp*dt << " " <<   tmpX(0,0).real() << " " << tmpPi(0,0).real() << std::endl;
  }
  myfile.close();
  
  X.write_to_hdf5("X.h5","X");
  Pi.write_to_hdf5("Pi.h5","Pi");
  
  return 0;
}




