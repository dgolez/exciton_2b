#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <chrono>
#include <ctime>

#include "cntr/cntr.hpp"
#include <Eigen/Eigenvalues>

using namespace std;
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double> 
#define CFUNCTION cntr::function<double>

/*///////////////////////////////////////////////////////////////////////////////////////

Test of the downfolding on the three level system

///////////////////////////////////////////////////////////////////////////////////////*/ 

// Embeded self energy for subspace S elements (i,j)
// from the total space T, where 
// the rest of the space R is integrated out 
// Assumption that the off-diagonal terms are not retarded

void downfold(int tstp,int i,int j,GREEN &R,GREEN &S,CFUNCTION &epsilon){
	int nt=R.nt();
	int ntau=R.ntau();
	int ssize=S.size1();
	int Rsize=R.size1();

	GREEN_TSTP etatmp(tstp,ntau,ssize);
	cntr::function<double> epsR(nt,ssize),epsL(nt,ssize);
	etatmp.clear();
	for(int k=0;k<Rsize;k++){
		for(int l=0;l<Rsize;l++){
			GREEN_TSTP tmp(tstp,ntau,1);
			tmp.set_matrixelement(0,0,R,k,l);
			epsR.set_matrixelement(0,0,epsilon,ssize+l,j);
			epsL.set_matrixelement(0,0,epsilon,i,ssize+k);
			tmp.right_multiply(epsR,1.0);
			tmp.left_multiply(epsL,1.0);
			etatmp.incr(tmp,1.0);
		}
	}
	S.set_matrixelement(tstp,i,j,etatmp);
}


int main(int argc, char *argv[]) {
	int size=2,ssize=1;
	int nt=500;
	int ntau=300;
	double beta=1.0,eps=1e-6;
	double h=0.03;
	int kt=5;

	const std::complex<double> i(0, 1);
	const double pi = std::acos(-1);

	GREEN A=GREEN(nt,ntau,size,-1);
	GREEN B=GREEN(nt,ntau,size,-1);

	CFUNCTION epsilon(nt,size,size),epsilon11(nt,size-ssize,size-ssize),epsilon00(nt,ssize,ssize);
	CFUNCTION U(nt,size,size),Ut(nt,size,size);
	
	
		
	cdmatrix a(size,size);
	a(0,0)=1.0;
	a(1,1)=2.0;
	a(0,1)=std::exp(i*pi);
	a(1,0)=std::exp(-i*pi);
	epsilon.set_constant(a);
	cntr::green_from_H(A,0.0,epsilon,beta,h,5,4,"true");

	A.write_to_hdf5("A.h5","G");
	// And eigen-basis
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> est(a);
	std::cout << "Eigenval "<<  est.eigenvalues()[0] << " " << est.eigenvalues()[1]  << std::endl; 

	cdmatrix b(size,size);
	b=est.eigenvectors().transpose()*a*est.eigenvectors();
	std::cout << "Check " << b << std::endl;
	U.set_constant(est.eigenvectors());
	Ut.set_constant(est.eigenvectors().transpose());

	for(int tstp=-1;tstp<nt;tstp++){
		A.right_multiply(tstp,U,1.0);
		A.left_multiply(tstp,Ut,1.0);
	}

	A.write_to_hdf5("A1.h5","G");
}

