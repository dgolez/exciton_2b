#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include <memory>
#ifndef CNTR_USE_OMP
#define CNTR_USE_OMP
#endif
#include "cntr/cntr.hpp"
#include "cntr/utils/read_inputfile.hpp"
#include "../helpers/inclusions.hpp"
#include "../helpers/print.hpp"
#include "../helpers/parameters.hpp"

using namespace std; 
#if CNTR_USE_MPI==1
#define USE_MPI 1
#endif

#include <time.h>
#include <sys/time.h>

parameters::parameters(int nt,int kt,int nk,int ntau,int size,double beta,double h,double mu,double den, double xi,std::vector<double> &delta,double v01,double v01_time,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &omega0,CFUNC &g,std::vector<double> &E,double dipolRe,double dipolIm,double ratio,double fieldP,double fieldD,double mazza,bool update,double eta,double mix,bool test,int omp_for_vie2,int phonontype,double bath_low,double bath_high,std::vector<double> &gC_bath,std::vector<double> &gV_bath,int migdal,std::vector<double> &gBATH,std::vector<double> &omegaBATH): 
		nt(nt), ntau(ntau), size(size), nk(nk), kt(kt), beta(beta), h(h), mu(mu), den(den),xi(xi),delta(delta),v01(v01),v01_time(v01_time),
		tt(tt),U(U),V(V),omega0(omega0),g(g),E(E),dipol(dipolRe,dipolIm),ratio(ratio),fieldP(fieldP),fieldD(fieldD),mazza(mazza),update(update),eta(eta),mix(mix),test(test),omp_for_vie2(omp_for_vie2),phonontype(phonontype),bath_low(bath_low),bath_high(bath_high),gC_bath(gC_bath),gV_bath(gV_bath),migdal(migdal),gBATH(gBATH),omegaBATH(omegaBATH){};

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

void spy(int tstp,const char *file_prefix){
	char filename[1000];
	FILE *out;
	sprintf(filename,"%s_status.out",file_prefix);
	out=std::fopen(filename,"w");
	if(out==0){
		std::cout << "cannot create file " << filename << std::endl;
		abort();
	}
	fprintf(out,"finished step %d",tstp);
	std::fclose(out);
}

int main(int argc,char *argv[]){
	// std::cout << "tu smo 1 " << std::endl;
	int tid,ntasks,ntau,nt,itermax,iter_rtime,omp_for_vie2,kt,nk,tstp=-2,iter,print_k,outputfrequency,size=2;
	int restart_time,cont,use_rpa=true;
	bool test,update;
	double mu,beta,h,errmax,errmax_rtime,err,v01,v01_time,mazza,mix,den;
	double wtime_limit,wtime_start,wtime_end,xi,dipolRe,dipolIm,eta;
	std::vector<double> U,V,tt,delta,E,gVec,omega0Vec,gBATH,omegaBATH;
	approx<lattice_1d_2b_optical_nofield_abinitio> *lattice;
	int migdal;
	double om0,s,amp,fieldD,fieldP,ratio,phonontype;
	//Bath parameters
	double bath_low,bath_high;
	std::vector<double> gC_bath,gV_bath;
	
	// (I) MPI initialization ... uses MPI over kpoints with one omp task per rank
	{
		#if USE_MPI==1
			MPI_Init(&argc,&argv);
			MPI_Comm_size(MPI_COMM_WORLD,&ntasks);
			MPI_Comm_rank(MPI_COMM_WORLD, &tid);
		#else
			ntasks=1;
			tid=0;
		#endif
		wtime_start=get_wall_time();
	}
	
	////////////////////////////////////////////////////////////////////////////////////	
	// (II) READ GENERAL INPUT
	{
		// if(argc<3) throw("COMMAND LINE ARGUMENT MISSING");
		// argv[1]: input-file (including full path)
		// argv[2]: output-file prefix (including full path): output-files are called argv[2]_some_observable etc.
		// argv[3]: needed only for restart at timestep >= -1
		find_param(argv[1],"__nt=",nt);
		find_param(argv[1],"__ntau=",ntau);
		find_param(argv[1],"__beta=",beta);
		find_param(argv[1],"__h=",h);    
		find_param(argv[1],"__mu=",mu);
		find_param(argv[1],"__den=",den);
		find_param(argv[1],"__itermax=",itermax);
		find_param(argv[1],"__errmax=",errmax);
		find_param(argv[1],"__iter_rtime=",iter_rtime);
		find_param(argv[1],"__err_rtime=",errmax_rtime);
		find_param(argv[1],"__kt=",kt);
		find_param(argv[1],"__print_k=",print_k); //Write k-dependent quantities to output
		find_param(argv[1],"__omp_for_vie2=",omp_for_vie2);
		find_param(argv[1],"__outputfrequency=",outputfrequency);
		find_param(argv[1],"__restart_time=",restart_time);
		find_param(argv[1],"__walltime_limit=",wtime_limit);
		find_param(argv[1],"__test=",test);
		find_param(argv[1],"__update=",update);
		find_param(argv[1],"__fieldD=",fieldD);
		find_param(argv[1],"__fieldP=",fieldP);
		find_param(argv[1],"__eta=",eta);
	}
	//////////////////////////////////////////////////////////////////////////////////
	// init the lattice and the physical input parameters
	{
		find_param(argv[1],"__nk=",nk);
		find_param(argv[1],"__xi=",xi);
		find_param(argv[1],"__v01=",v01);	//Seed for pair breaking field (present in hamiltonian in first 5 iterations for equilibrium)
		find_param(argv[1],"__v01_time=",v01_time); // Source field for evaluating excitonic susceptibility
		find_param(argv[1],"__mazza=",mazza);
		find_param(argv[1],"__mix=",mix);
		find_param(argv[1],"__phonontype=",phonontype);
		find_param_tvector(argv[1],"__U=",U,nt);
		find_param_tvector(argv[1],"__V=",V,nt);
		find_param_tvector(argv[1],"__g=",gVec,nt);
		find_param_tvector(argv[1],"__omega0=",omega0Vec,nt);
		find_param_tvector(argv[1],"__tt=",tt,nt);
		find_param_tvector(argv[1],"__delta=",delta,nt);
		find_param_tvector(argv[1],"__E=",E,nt);
		find_param(argv[1],"__dipolRe=",dipolRe);
		find_param(argv[1],"__dipolIm=",dipolIm);
		find_param(argv[1],"__om0=",om0);
		find_param(argv[1],"__migdal=",migdal);
		find_param(argv[1],"__s=",s);
		find_param(argv[1],"__amp=",amp);
		find_param(argv[1],"__ratio=",ratio); //Ratio  of masses

		// Add bath fermions
      	find_param(argv[1], "__bath_low=", bath_low);
      	find_param(argv[1], "__bath_high=", bath_high);
      	find_param_tvector(argv[1],"__gC_bath=",gC_bath,nt);
      	find_param_tvector(argv[1],"__gV_bath=",gV_bath,nt);
      	// Add bath Holstein like phonons
      	find_param_tvector(argv[1],"__gBATH=",gBATH,nt);
		find_param_tvector(argv[1],"__omegaBATH=",omegaBATH,nt);

		 // Change in the initial conditions to obtain exciton response function
		// NOTE: whenever parameters in lattice_rpa.latt_ are reset,
		// this is automatically recognized by lattice_rpa
	}
	std::cout << "tu smo2 "<< std::endl;
	if(nt<=10 && test){
		test=test;
	}else{
		test=false;
	}
	//Set cfunc for g, omega0
	CFUNC g(nt,1),omega0(nt,1);
	{
	  cdmatrix tmpG(1,1),tmpOm(1,1);
	  for(int tstp=-1;tstp<=nt;tstp++){
	  	// gX
	    tmpG(0,0)=gVec[tstp+1];
	    g.set_value(tstp,tmpG);
	    // omega
	    tmpOm(0,0)=omega0Vec[tstp+1];
	    omega0.set_value(tstp,tmpOm);
	  }
	}

	std::cout << den << std::endl;
	parameters param(nt,kt,nk,ntau,size,beta,h,mu,den,xi,delta,v01,v01_time,tt,U,V,omega0,g,E,dipolRe,dipolIm,ratio,fieldP,fieldD,mazza,update,eta,mix,test,(omp_for_vie2==1 ? true : false),phonontype,bath_low,bath_high,gC_bath,gV_bath,migdal,gBATH,omegaBATH);

	lattice=new mpi_lattice_step_2b_optical<lattice_1d_2b_optical_nofield_abinitio>(param);

	if(restart_time!=-1){
		if(argc<4) throw("INPUT FILE PREFIX (including path) NEEDED FOR RESTART");
		if(restart_time>nt) throw("RESTART TIME > NT ");
		lattice->read_from_file_hdf5(restart_time-1,argv[3],kt);
		// TODO I you want restart
		// for(tstp=-1;tstp<restart_time;tstp++) lattice_rpa.gather_kk_observables(tstp,kt);
		// TODO I you want restart
	}else if(restart_time==-1 && argc==4){
		if(argc<4) throw("INPUT FILE PREFIX (including path) NEEDED FOR RESTART");
		std::cout << "Seeding the equilibrium solution" << std::endl;
		lattice->read_from_file_hdf5(-1,argv[3],kt);
		// lattice->gather_kk_observables(-1,kt);
	}
	std::ofstream out;
	out.open("disper.out");
	for(int k=0;k<lattice->latt_.nk_;k++){
	  Eigen::MatrixXcd vtmp;
	  out << lattice->latt_.kpoints_[k];
	  lattice->latt_.hk(vtmp,-1,lattice->latt_.kpoints_[k],1);
	  out << " " << vtmp(0,0).real() << " " << vtmp(1,1).real()  <<  " ";
	  lattice->latt_.hkfree(vtmp,-1,lattice->latt_.kpoints_[k]);
	  out << " " << vtmp(0,0).real() << " " << vtmp(1,1).real() << std::endl;
	}
	out.close();
	if(tid==0) spy(tstp,argv[2]);
	//////////////////////////////////////////////////////////////////////////////////
	// initialization: G is initially zero, thus also Sigma=0 and the first step 
	// should produce non-interacting Green's functions
	for(tstp=restart_time;tstp<=nt;tstp++){
		int itermax1;
		double errmax1;
		int kt1;
		if(tstp==-1){
			kt1=kt;
			errmax1=errmax;
			itermax1=itermax;
		}else{
			kt1=(tstp>=kt ? kt : tstp);
			errmax1=errmax_rtime;
			itermax1=iter_rtime;			
		}

		// if(tstp>0) lattice->
		for(iter=1;iter<=itermax1;iter++){
		  err=lattice->step(tstp,iter,kt,om0,s,amp); /// should work at time zero!
		  if(tid==0){
		  	cdmatrix tmp(2,2),tmpsym(2,2),ord(1,1);
		  	lattice->rho_loc_.get_value(tstp,tmp);
		  	lattice->rho_sym_.get_value(tstp,tmpsym);
			lattice->order_.get_value(tstp,ord);
		  	cout << "tstp= " << tstp << " iter:  " << iter << " err: " << err << " - " << tmp << " " << tmp.trace() << " order_cos: " << tmpsym(0,1) << " order: " << ord(0,0)  << " "   << endl;
		  }
		  // Set up the local density matrix for the Hartree shift
		  cdmatrix tmp;
		  lattice->rho_loc_.get_value(-1,tmp);
		  lattice->latt_.rho_eq_=tmp;
		  if(err<errmax1 && iter>15 ){
		    if(tstp==-1){
		      out.open("disper_eq.out");
		      for(int k=0;k<lattice->latt_.nk_;k++){
				Eigen::MatrixXcd vtmp;
				out << lattice->latt_.kpoints_[k];
				lattice->density_k_[k].hkeff_.get_value(-1,vtmp);
				out << " " << vtmp(0,0).real() << " " << vtmp(0,1).real() << " " << vtmp(1,0).real() << " " << vtmp(1,1).real()  <<  " ";
				lattice->density_k_[k].hkeff_eigen_.get_value(-1,vtmp);
				out << " " << vtmp(0,0).real() << " " << vtmp(1,1).real()  <<  " ";
				lattice->latt_.hkfree(vtmp,-1,lattice->latt_.kpoints_[k]);
				out << " " << vtmp(0,0).real() << " " << vtmp(0,1).real() << " " << vtmp(1,0).real() << " " << vtmp(1,1).real()  << std::endl;
		      }
		      out.close();

		    }
		    break;
		  }
		}
		if(tstp==-1 and err>errmax1){
			std::cout << "Matsubara didn't converged" << std::endl;
			abort();
		}
		//lattice->get_optical0_eq(0.0,0,0);
		//std::cout << "Matsubara converged" << std::endl;
		/// timestep converged: collect observables
		//lattice_rpa.gather_kk_observables(tstp,kt);
		//std::cout << "obtained kk" << std::endl;
		// lattice_rpa.get_Gloc(tstp);
		//std::cout << "obtained gloc" << std::endl;
		// if the walltime exceeds the limit: break and save the result
		if(tid==0){
		  wtime_end=get_wall_time();	
		  cont=( wtime_end-wtime_start>wtime_limit ? 0 : 1);
		  std::cout << "finished tstp= " << tstp << " after " << wtime_end-wtime_start << " seconds" << std::endl;
		}
		MPI_Bcast(&cont, 1, MPI_INTEGER,0, MPI_COMM_WORLD);
		if(!cont){
		  if(tid==0){
		    std::cout << "ended calculations after " << wtime_end-wtime_start << " seconds" << " at tstp= " << tstp << std::endl;
		  }
		  tstp++;
		  break;
		}
		if(tstp==-1 || tstp%outputfrequency==0) if(tid==0) spy(tstp,argv[2]);
	}
	// lattice->print_to_file_hdf5_slice(argv[2],outputfrequency);
	// lattice_rpa.print_diagonal_occupations(argv[2],-1,tstp-1,kt,h,size);
	//Evalutate the order parameter


	char filename[1000];
	std::cout << "output " << argv[2]<< std::endl;
	sprintf(filename,"%s_full",argv[2]);
	std::cout << "tu smo 1 " << tid  << std::endl;
	lattice->print_to_file_hdf5(filename,print_k);
	std::cout << "tu smo 2 " << tid  << std::endl;
	lattice->print_to_file_hdf5_slice(argv[2],outputfrequency,print_k);
	////////////////////////////////////////////////////////////////////////////////////
	// .....
	///////////////////////////////////////////////////////////////////////////////////
	
	MPI_Finalize();	
	return 0;
}
