#pragma once

// parameter classes also for hubbard and tns programs, but not used in exciton_photo program

class parameters{
	public:
		parameters(){};
		parameters(int nt,int kt,int nk,int ntau,int size,double beta,double h,double mu,double den, double epsilon, double xi,std::vector<double> &delta,double v01,double v01_time,int suscep,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &omega0,CFUNC &g,std::vector<double> &E,std::vector<double> &dE,double dipolRe,double dipolIm,double ratio,double fieldP,double fieldD,double mazza,bool update,double eta,double gamma,double mix,bool test,int omp_for_vie2,int phonontype,double bath_low,double bath_high,std::vector<double> &gC_bath,std::vector<double> &gV_bath,int migdal,std::vector<double> &gBATH,std::vector<double> &omegaBATH);

	int nt,ntau,size,nk,kt,suscep,omp_for_vie2,phonontype,migdal;
	double beta,h,mu,den,epsilon,xi,v01,v01_time,ratio,fieldP,fieldD,mazza,eta,gamma,mix;
	double bath_low,bath_high;
	std::vector<double> gC_bath,gV_bath;
	std::vector<double> gBATH,omegaBATH;
 	std::complex<double> dipol;
	std::vector<double> tt,U,V,E,dE,delta;
	CFUNC omega0,g;
	bool update,test;

};

class parameters_hub{
	public:
		parameters_hub(){};
		parameters_hub(int nt,int kt,int nk,int ntau,int size,double beta,double h,double mu,double den, double epsilon, double xi,std::vector<double> &delta,double v01,int suscep,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &omega0,CFUNC &g,std::vector<double> &E,std::vector<double> &dE,double dipolRe,double dipolIm,double ratio,double fieldP,double fieldD,double mazza,bool update,double eta,double gamma,double mix,bool test,int omp_for_vie2,int phonontype,double bath_low,double bath_high,std::vector<double> &gC_bath,std::vector<double> &gV_bath,int migdal);

	int nt,ntau,size,nk,kt,suscep,omp_for_vie2,phonontype,migdal;
	double beta,h,mu,den,epsilon,xi,v01,ratio,fieldP,fieldD,mazza,eta,gamma,mix;
	double bath_low,bath_high;
	std::vector<double> gC_bath,gV_bath;
	std::complex<double> dipol;
	std::vector<double> tt,U,V,E,dE,delta;
	CFUNC omega0,g;
	bool update,test;

};

class parameters_tns{
	public:
		parameters_tns(){};
		parameters_tns(int nt,int kt,int nk,int ntau,int size,double beta,double h,double mu,double den, double epsilon,double v01,std::vector<double> &U,std::vector<double> &V,std::vector<double> &Ex,std::vector<double> &Ey,std::complex<double> &dipol,bool update,double eta,double gamma,double mix,int omp_for_vie2);

	int nt,ntau,size,nk,kt,omp_for_vie2;
	double beta,h,mu,den,epsilon,v01,ratio,fieldP,fieldD,eta,gamma,mix;
	std::complex<double> dipol;
	std::vector<double> tt,U,V,Ex,Ey;
	bool update;
};