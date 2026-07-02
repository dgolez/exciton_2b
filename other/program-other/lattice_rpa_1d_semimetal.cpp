#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include <memory>
#define CNTR_USE_OMP
#define CNTR_USE_MPI
#include "cntr/cntr.hpp"
#include "cntr/utils/read_inputfile.hpp"
#include "../step/step.hpp"
#include "print.hpp"

using namespace std; 
#if CNTR_USE_MPI==1
#define USE_MPI 1
#endif


#include <time.h>
#include <sys/time.h>
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
	int tid,ntasks,ntau,nt,itermax,iter_rtime,omp_for_vie2,kt,nk,tstp=-2,iter,print_k,outputfrequency,size=4;
	int restart_time,cont;
	bool test;
	double mu,beta,h,errmax,errmax_rtime,err,v01,epsilon,mix;
	double wtime_limit,wtime_start,wtime_end,xi;
	std::vector<double> U,V,tt,delta,Ak,dAk,gVec,omega0Vec;
	approx<lattice_1d_rpa_cdw_nofield> *lattice;
	int approximation;
	// (I) MPI initialization ... uses MPI over kpoints with one omp task per rank
	{
		#if USE_MPI==1
			MPI_Init(&argc,&argv);
			MPI_Comm_size (MPI_COMM_WORLD,
&ntasks);
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
		find_param(argv[1],"__approximation=",approximation); // 0 - mean field, 1 - gw, 2 - gkba
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
	}
	//////////////////////////////////////////////////////////////////////////////////
	// init the lattice and the physical input parameters
	{
		find_param(argv[1],"__nk=",nk);
		find_param(argv[1],"__xi=",xi);
		find_param(argv[1],"__v01=",v01);	//Seed for pair breaking field (present in hamiltonian in first 5 iterations for equilibrium)
		find_param(argv[1],"__mix=",mix);
		find_param_tvector(argv[1],"__U=",U,nt);
		find_param_tvector(argv[1],"__V=",V,nt);
		find_param_tvector(argv[1],"__g=",gVec,nt);
		find_param_tvector(argv[1],"__omega0=",omega0Vec,nt);
		find_param_tvector(argv[1],"__tt=",tt,nt);
		find_param_tvector(argv[1],"__delta=",delta,nt);
		find_param_tvector(argv[1],"__Ak=",Ak,nt);
		find_param_tvector(argv[1],"__dAk=",dAk,nt);
		find_param(argv[1],"__epsilon=",epsilon); // Change in the initial conditions to obtain exciton response function
		// NOTE: whenever parameters in lattice_rpa.latt_ are reset,
		// this is automatically recognized by lattice_rpa
	}
	if(nt<=10 && test){
		test=test;
	}else{
		test=false;
	}
	//Set cfunc for g, omega0
	CFUNC g(nt,1),omega0(nt,1);
	{
	  for(int tstp=-1;tstp<=nt;tstp++){
	    cdmatrix tmpG(1,1),tmpOm(1,1);
	    tmpG(0,0)=gVec[tstp+1];
	    tmpOm(0,0)=omega0Vec[tstp+1];
	    g.set_value(tstp,tmpG);
	    omega0.set_value(tstp,tmpOm);
	  }
	}
	switch(approximation){
		case 0 : lattice=new mpi_lattice_step_cdw<lattice_1d_rpa_cdw_nofield>(nt,ntau,size,beta,h,(omp_for_vie2==1 ? true : false),mu,epsilon,tt,U,V,g,omega0,xi,delta,Ak,dAk,v01,nk,kt);break;
		// case 1 : lattice=new mpi_lattice_step_GKBA<lattice_1d_2b_nofield>(nt,ntau,size,beta,h,(omp_for_vie2==1 ? true : false),mu,epsilon,tt,U,V,g,omega0,xi,delta,Ak,dAk,v01,nk,kt);break;
		case 2 : lattice=new mpi_lattice_step_RPA_cdw<lattice_1d_rpa_cdw_nofield>(nt,ntau,size,beta,h,(omp_for_vie2==1 ? true : false),mu,epsilon,tt,U,V,g,omega0,xi,delta,Ak,dAk,v01,nk,kt,test,mix);break;
	case 3 : lattice=new mpi_lattice_step_2b_cdw<lattice_1d_rpa_cdw_nofield>(nt,ntau,size,beta,h,(omp_for_vie2==1 ? true : false),mu,epsilon,tt,U,V,g,omega0,xi,delta,Ak,dAk,v01,nk,kt,test,mix);break;
	}
	// TODO
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
	std::ofstream out,outeig,outfree;
	out.open("disper.out");
	for(int k=0;k<lattice->latt_.nk_rbz_;k++){
		Eigen::MatrixXcd vtmp;
		out << lattice->latt_.kpoints_rbz_[k];
		lattice->latt_.hk(vtmp,-1,lattice->latt_.kpoints_rbz_[k],1);
		out << " " << vtmp(0,0).real() << " " << vtmp(1,1).real()  <<  " ";
		lattice->latt_.hkfree(vtmp,-1,lattice->latt_.kpoints_rbz_[k]);
		out << " " << vtmp(0,0).real() << " " << vtmp(1,1).real() << std::endl;
	}
	out.close();
	if(tid==0) spy(tstp,argv[2]);
	//////////////////////////////////////////////////////////////////////////////////
	// initialization: G is initially zero, thus also Simga=0 and the first step 
	// should produce non-interactng Greens functions
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
		  err=lattice->step(tstp,iter,kt,0.1,0.1,0.0); /// should work at time zero!
		  if(tid==0){
		  	cdmatrix tmp;
		  	lattice->rho_loc_.get_value(tstp,tmp);
		  	cout << "tstp= " << tstp << " iter:  " << iter << " err: " << err << std::endl;
		  	std::cout<< " Den mat: " << tmp << std::endl;
		  	std::cout << " Trace " << tmp.trace() << endl;
		  }
		  if(err<errmax1 && iter>2 ){
		  // 	if(tstp==-1){
		  //   	out.open("disper_eff.out");
		  //   	outeig.open("disper_eig.out");
		  //   	outfree.open("disper_free.out");

		  //     	for(int k=0;k<lattice->latt_.nk_rbz_;k++){
			 //    Eigen::MatrixXcd vtmp;
				// out << lattice->latt_.kpoints_[k];
				// lattice->density_k_[k].hkeff_.get_value(-1,vtmp);
				// for(int i=0;i<4;i++){
				// 	for(int j=0;j<4;j++){
				// 		out << " " << vtmp(i,j).real() << " ";
				// 	}
				// }
				// lattice->density_k_[k].hkeff_eigen_.get_value(-1,vtmp);
				// outeig << lattice->latt_.kpoints_[k];
				// for(int i=0;i<4;i++){
				// 		outeig << " " << vtmp(i,i).real() << " ";
				// }
				
				// lattice->latt_.hkfree(vtmp,-1,lattice->latt_.kpoints_[k]);
				// outfree << lattice->latt_.kpoints_[k];
				// for(int i=0;i<4;i++){
				// 	for(int j=0;j<4;j++){
				// 		outfree << " " << vtmp(i,j).real() << " ";
				// 	}
				// }
				// out << std::endl;
				// outeig << std::endl;
				// outfree << std::endl;
		  //     }
		  //     out.close();
		  //     outfree.close();
		  //     outeig.close();
		  //   }
		    break;
		  } 
		}
		if(tstp==-1 and err>errmax1){
			std::cout << "Matsubara didn't converged" << std::endl;
			abort();
		}
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
	// std::cout << "Print " << print_k << " " << approximation << std::endl;
	char filename[1000];
	sprintf(filename,"%s_full",argv[2]);
	lattice->print_to_file_hdf5(filename,print_k);
	
	lattice->print_to_file_hdf5_slice(argv[2],outputfrequency,print_k);
	////////////////////////////////////////////////////////////////////////////////////
	// .....
	///////////////////////////////////////////////////////////////////////////////////
	
	MPI_Finalize();	
	return 0;
}




