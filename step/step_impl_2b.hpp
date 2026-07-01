#pragma once

#include "step.hpp"

template <class LATTICE>
mpi_lattice_step_2b<LATTICE>::mpi_lattice_step_2b(int nt,int ntau,int size,double beta,double h,bool use_omp_for_vie2,double mu,double epsilon,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &g,CFUNC &omega0,CFUNC &gBATH,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,int nk,int kt,bool test,double mix):
	base_type::mpi_lattice_step(nt,ntau,size,beta,h,use_omp_for_vie2,mu,epsilon,tt,U,V,g,omega0,gBATH,xi,delta,A,dA,v01,nk,kt){
	// This is a different type that the base one(includes the info about the RPA )
	// One can use this to reduce the amount of this-> for repeated variables
	// std::cout << "step 1 " << std::endl;
	int & nk_=this->nk_;
	int & nt_=this->nt_;
	int & ntau_=this->ntau_;
	int & nrpa_=this->nrpa_;
	test_=test;
	LATTICE & latt_=this->latt_;
	// std::cout << "step 2 " << std::endl;
	green_k_.resize(nk_);
	assert(nt==latt_.nt_);
	assert(0<=ntau_);
	// std::cout << "step 3 " << std::endl;
	// base_type::init(nt,ntau,size,beta,h,use_omp_for_vie2,mu,epsilon);
	gk_all_timesteps_=cntr::distributed_timestep_array<double>(nk_,nt_,ntau_,nrpa_,-1,true);
	wk_all_timesteps_=cntr::distributed_timestep_array<double>(nk_,nt_,ntau_,nrpa_,+1,true);
	dk_all_timesteps_=cntr::distributed_timestep_array<double>(nk_,nt_,ntau_,1,+1,true);
	this->convergence_error_=cntr::distributed_array<double>(nk_,1,true);
	this->convergence_error_.reset_blocksize(1);
	// Times 4 due to four variables stored in the k_observables
	// this->k_observables_=cntr::distributed_array<cdouble>(nk_,4*nrpa_*nrpa_,true);
	// this->k_observables_.reset_blocksize(4*nrpa_*nrpa_);
	// Why is there 4 ?
	// std::cout << "step 4 " << std::endl;
	/////////////////////////////////////////////////
	// MPI Setup;  gk_..., wk_..., and *this will have the same mpi-layout
	this->tid_=gk_all_timesteps_.tid(); 
	this->ntasks_=gk_all_timesteps_.ntasks();
	this->tid_map_=gk_all_timesteps_.data().tid_map();
	this->tid_root_=0;
	if(this->tid_==this->tid_root_){
		Gloc_=GREEN(nt_,ntau_,nrpa_,-1);
		Wloc_=GREEN(nt_,ntau_,nrpa_,1);
		Dloc_=GREEN(nt_,ntau_,1,1);
	} 
	// using OMP in the solution of vie2 at each kpoint?
	use_omp_for_vie2_=use_omp_for_vie2;
	/////////////////////////////////////////////////
	
    // std::cout << "step 5d " << std::endl;
	for(int k=0;k<nk_;k++){
		// std::cout << "step 5e " << std::endl;
	  if(this->tid_map_[k]==this->tid_){
	  	// std::cout << "step 5f " << std::endl;
		this->green_k_[k]=kpoint_green<LATTICE>(nt,ntau,size,beta,h,latt_.kpoints_[k],latt_,mu,this->epsilon_,omega0,g,mix);
		// std::cout << "step 5g " << std::endl;
		this->green_k_[k].use_omp(use_omp_for_vie2_);
		// std::cout << "step 5h " << std::endl;
	  }
	}
	// std::cout << "step 6 " << std::endl;
	init_G_mat_nointeraction();
	// std::cout << "step 7 " << std::endl;
	// gather_kk_observables(-1,kt);
	// std::cout << "step 8 " << std::endl;
}



template <class LATTICE>
void mpi_lattice_step_2b<LATTICE>::init_G_mat_nointeraction(void){
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			// std::cout << "Init " << k << std::endl;
			this->green_k_[k].init_G_mat_nointeraction(this->latt_,this->density_k_[k],kt_);
		}
	}
}

template <class LATTICE>
void mpi_lattice_step_2b<LATTICE>::extrapolate_timestep_G(int tstp,int kt){
  assert(tstp<=this->nt_-1);
	// extrapolate from tstp to tstp+1
	if(tstp<-1) return;
	if(tstp==-1){
		for(int k=0;k<this->nk_;k++){
		  if(this->tid_map_[k]==this->tid_) cntr::set_t0_from_mat(this->green_k_[k].G_);
		}
	}
	for(int k=0;k<this->nk_;k++){
	  if(this->tid_map_[k]==this->tid_) cntr::extrapolate_timestep(tstp,this->green_k_[k].G_,integration::I<double>(kt));
	}
}

template <class LATTICE>
void mpi_lattice_step_2b<LATTICE>::extrapolate_timestep_W(int tstp,int kt){
  assert(tstp<=this->nt_-1);
	// extrapolate from tstp to tstp+1
	if(tstp<-1) return;
	if(tstp==-1){
		for(int k=0;k<this->nk_;k++){
		  if(this->tid_map_[k]==this->tid_) cntr::set_t0_from_mat(this->green_k_[k].W_);
		}
	}
	for(int k=0;k<this->nk_;k++){
	  if(this->tid_map_[k]==this->tid_) cntr::extrapolate_timestep(tstp,this->green_k_[k].W_,integration::I<double>(kt));
	}
}

template <class LATTICE>
void mpi_lattice_step_2b<LATTICE>::symmetrise_mat(GREEN &G){
	cdmatrix M,M1;
	for(int i=0;i<this->ntau_;i++){
		G.get_mat(i,M);
		M1=(M+M.adjoint())/2.0;
		G.set_mat(i,M1);
		}
}

// TODO: check if it is better to gather density from extrapolated Gk or from rho. I'd guess rho
// TODO: CLEAN UP declaration of same variables in the base and this class
template <class LATTICE>
double mpi_lattice_step_2b<LATTICE>::step(int tstp,int iter,int kt,double om0,double s,double amp){
  assert(tstp<=this->nt_);
	double err;
	int n;
	int n1=(tstp==-1 || tstp>kt ? tstp : 0);
	int n2=(tstp==-1 || tstp>kt ? tstp : kt);
	//std::cout << "STEP: " << tstp << ", n1: " << n1 << ", n2: " << n2 << " " << kt  << " " <<  (tstp==-1 || tstp>kt) << std::endl;
	if(tstp>0 && iter==1){
		this->extrapolate_rho(tstp);
		extrapolate_timestep_G(tstp-1,(kt<tstp-1 ? kt : tstp-1));
		extrapolate_timestep_W(tstp-1,(kt<tstp-1 ? kt : tstp-1));
	}
	// std::cout << "Tu smo 1 " << std::endl;
	for(n=n1;n<=n2;n++){
	  //std::cout << "Loop Pol " << tstp << " " <<  n << " " << n1 << " " << " " << n2 << std::endl;
		//if(n1!=n2)
		gather_gk_timestep(n); //Gather timesteps of all k dependent green's functions
		// set_density_k(n);
		// std::cout << "Tu smo 2 " << std::endl;
		for(int k=0;k<this->nk_;k++){
		  if(this->tid_map_[k]==this->tid_){
		    if(test_){
		      std::ostringstream nameG;
		      nameG << "Gk0_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
		      green_k_[k].G_.write_to_hdf5(nameG.str().c_str(),"G");	
		    }
		    // std::cout << "Tu smo 3 " << std::endl;
		    // get_Polarization_Bubble(n,k,green_k_[k].P_);
		    // The full 2b evaluates a different "bubble" than GW 
		    get_PP_Bubble(n,k,green_k_[k].P_);
		    // std::cout << "Tu smo 4 " << std::endl;
		    //cdmatrix tmp;
		    //green_k_[k].P_.get_mat(0,tmp);
		    //std::cout << "Polar " << this->density_k_[k].kk_ << " " << tmp << std::endl; 
		  }
		}
	}
	// std::cout << "Tu smo 5 " << std::endl;
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			green_k_[k].step_W2b(tstp,kt,this->latt_,this->density_k_[k]);
			// std::cout << "Tu smo 5a " << std::endl;
			green_k_[k].step_D(tstp,kt,this->latt_,this->density_k_[k]);
			// std::cout << "Tu smo 5b " << std::endl;
			if(test_){
			  cdmatrix tmp;
			  green_k_[k].W_.get_mat(0,tmp);
			  //std::cout << "W " << this->density_k_[k].kk_ << " " << tmp << std::endl;
				for(int n=n1;n<=n2;n++){
					std::ostringstream nameP;
					nameP << "P_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
					green_k_[k].P_.write_to_hdf5(nameP.str().c_str(),"G");

					std::ostringstream nameW;
					nameW << "W_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
					green_k_[k].W_.write_to_hdf5(nameW.str().c_str(),"G");
				}	
			}	
		} 
	}
	// std::cout << "Tu smo 6 " << std::endl;
	// Set local
	cdmatrix rtmp;
	this->set_local(tstp);
	// std::cout << "Tu smo 6a " << std::endl;
	this->phonon_=phonon(latt_.omega0_,latt_.g_,this->rho_loc_);
	// std::cout << "Tu smo 6b " << std::endl;
	for(n=n1;n<=n2;n++){

		cdmatrix tmp;
		// std::cout << "Loop Sigma " << tstp << " " <<  n << " " << n1 << " " << " " << n2 << std::endl;
		// std::cout << "Tu smo 6c " << std::endl;
		gather_Wk_timestep(n);
		// std::cout << "Tu smo 6d " << std::endl;
		gather_Dk_timestep(n);
		// std::cout << "Tu smo 6e " << std::endl;
		// gather_kk_observables(n,kt);
		if(test_){
			for(int k=0;k<this->latt_.nk_;k++){
				if(this->tid_map_[k]==this->tid_){
				
					std::ostringstream nameW;
					nameW << "Wp_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
					wk_all_timesteps_.G()[k].write_to_hdf5(nameW.str().c_str(),"G");

					std::ostringstream nameD;
					nameD << "D_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
					dk_all_timesteps_.G()[k].write_to_hdf5(nameD.str().c_str(),"G");
				}
			}
		}
		// std::cout << "Tu smo 7 " << std::endl;
		gather_gk_timestep(n);
		// std::cout << "Tu smo 8 " << std::endl;
		base_type::step_boson(n);
		// std::cout << "Tu smo 9 " << std::endl;
		for(int k=0;k<this->nk_;k++){
			if(this->tid_map_[k]==this->tid_){
				
				base_type::get_Sigma_Hartree(n,k,this->density_k_[k].SHartree_);
				
				this->density_k_[k].SHartree_.get_value(n,tmp);
				std::cout << "Hartree " << tmp << std::endl; 

				base_type::get_Sigma_Fock(n,k,this->density_k_[k].SFock_);
				this->density_k_[k].SFock_.get_value(n,tmp);
				std::cout << "Fock " << tmp << std::endl;
				
				get_Sigma_2b(n,k,green_k_[k].Sigma_,iter);
				green_k_[k].Sigma_.get_mat(0,tmp);
				std::cout <<  "SigmaRPA "   << n << " " << this->density_k_[k].kk_ << " " << tmp << std::endl;
				if(test_){
					std::ostringstream nameW;
					nameW << "Sigma_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
					green_k_[k].Sigma_.write_to_hdf5(nameW.str().c_str(),"G");	
				}
			}
		}
		
		// std::cout << "step 6e" <<std::endl;
	}
	// std::cout << "Tu smo 10 " << std::endl;
	
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			if(test_){
				std::ostringstream nameG;
				nameG << "Gkp_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
				green_k_[k].G_.write_to_hdf5(nameG.str().c_str(),"G");	
			}
			this->convergence_error_.block(k)[0]=green_k_[k].step_dyson_with_error_integral(tstp,iter,kt,this->latt_,this->density_k_[k]);
			if(test_){
				std::ostringstream nameW;
				nameW << "Gk_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
				green_k_[k].G_.write_to_hdf5(nameW.str().c_str(),"G");
			}
		}
	}
	// std::cout << "Tu smo 11 " << std::endl;
	//Symmetrize  the solutions
	// if(tstp==-1){
	// 	for(int k=0;k<this->nk_;k++){
	// 		if(this->tid_map_[k]==this->tid_){	
	// 			symmetrise_mat(green_k_[k].G_);
	// 			symmetrise_mat(green_k_[k].Sigma_);
	// 		}
	// 	}
	// }
	
	this->convergence_error_.mpi_bcast_all();
	err=0.0;
	for(int k=0;k<this->nk_;k++) err+=this->convergence_error_.block(k)[0];
	//Symmetrize  the solutions
	if(tstp==-1){
		for(int k=0;k<this->nk_;k++){
			if(this->tid_map_[k]==this->tid_){
				cdmatrix tmp;
				green_k_[k].G_.set_mat_herm();
				green_k_[k].Sigma_.set_mat_herm();
				green_k_[k].chi_.set_mat_herm();
				green_k_[k].P_.set_mat_herm();
				green_k_[k].W_.set_mat_herm();

				// symmetrise_mat(green_k_[k].G_);
				// symmetrise_mat(green_k_[k].Sigma_);
				// symmetrise_mat(green_k_[k].chi_);
				// symmetrise_mat(green_k_[k].P_);
				// ADD symmetrization of W !!
				// Change the order of step maybe it helps !!
				// Change the evaluation of W to the old version, maybe something tricky is happening !!
			}
		}
	}
	// std::cout << "Tu smo 12" << std::endl;
	// Later remove since it's doubling of the communication
	for(n=n1;n<=n2;n++){
		gather_gk_timestep(n); //Gather timesteps of all k dependent green's functions
		set_density_k(n);
		base_type::set_local(n);
		get_Gloc(n);
		get_Wloc(n);
		get_Dloc(n);
	}
	// std::cout << "Tu smo 13 " << std::endl;
	return err;
}


template <class LATTICE>
void mpi_lattice_step_2b<LATTICE>::gather_gk_timestep(int tstp){
	gk_all_timesteps_.reset_tstp(tstp);
	// read Gk_ to timestep
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
		  gk_all_timesteps_.G()[k].get_data(green_k_[k].G_);
		} 
	}
	// distribute to all nodes
	gk_all_timesteps_.mpi_bcast_all();
}

template <class LATTICE>
void mpi_lattice_step_2b<LATTICE>::set_density_k(int tstp){
	cdmatrix tmp(this->nrpa_,this->nrpa_);
	for(int k=0;k<this->nk_;k++){
	  gk_all_timesteps_.G()[this->latt_.representative_kk(k)].density_matrix(tstp,tmp);
		this->density_k_[this->latt_.representative_kk(k)].rho_.set_value(tstp,tmp);
	}
}

template <class LATTICE>
void mpi_lattice_step_2b<LATTICE>::gather_Wk_timestep(int tstp){
	wk_all_timesteps_.reset_tstp(tstp);
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
		  wk_all_timesteps_.G()[k].get_data(green_k_[k].W_);
		}
	}
	// distribute to all nodes
	wk_all_timesteps_.mpi_bcast_all();
}

template <class LATTICE>
void mpi_lattice_step_2b<LATTICE>::gather_Dk_timestep(int tstp){
	dk_all_timesteps_.reset_tstp(tstp);
	// std::cout << "Gather dk timestep 1" <<std::endl;
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			// std::cout << "Gather dk timestep 2 " <<std::endl;
			// cdmatrix tmp;
			// green_k_[k].D_.get_mat(0,tmp);
			// std::cout << "gather S " << tmp << std::endl;
		  	dk_all_timesteps_.G()[k].get_data(green_k_[k].D_);
		  	// dk_all_timesteps_.G()[k].get_mat(0,tmp);
		  	// std::cout << "gather S2 " << tmp << std::endl;
			// std::cout << "Gather dk timestep 3 " <<std::endl;
		}
	}
	// std::cout << "Gather dk timestep 4 " <<std::endl;
	// distribute to all nodes
	dk_all_timesteps_.mpi_bcast_all();
	// std::cout << "Gather dk timestep 5 " <<std::endl;
}


//get P_pp_q (a1,a2) = sum_k ii*G_k(a1,a2_;t,t')G_{q-k}(a1_,a2;t,t')
template <class LATTICE>
void mpi_lattice_step_2b<LATTICE>::get_PP_Bubble(int tstp,int qq,GREEN &P){
  assert(tstp==gk_all_timesteps_.tstp());	
	// compute P_{a,a'}(q;t,t') on a timestep
	int qk,i1,i2,qk1,kk1;
	GREEN_TSTP ptmp(tstp,this->ntau_,this->nrpa_,+1);
	P.set_timestep_zero(tstp);
	for(int kk=0;kk<this->latt_.nk_bz_;kk++){
		double wk=this->latt_.kweight_bz_[kk]; // factor 2: ksum normalized for RBZ
		int gammakq;
		// get bubble for 2b approximation
		// 
		//kq=latt_.add_kpoints(1,kk,1,qq);
		qk=this->latt_.add_kpoints(kk,-1,qq,1);
		//latt_.add_kpoints(kq, gammakq,1,kk,1,qq);
		qk1=this->latt_.representative_kk(qk);
		kk1=this->latt_.representative_kk(kk);
		ptmp.clear();

		//Note:BUBBLE2 IS C_{c1,c2}(t1,t2) = ii * A_{a1,a2}(t1,t2) * B_{b1,b2}(t1,t2) 
		cntr::Bubble2(tstp,ptmp,0,0,gk_all_timesteps_.G()[kk1],gk_all_timesteps_.G()[kk1],0,1,gk_all_timesteps_.G()[qk1],gk_all_timesteps_.G()[qk1],1,0);
		cntr::Bubble2(tstp,ptmp,0,1,gk_all_timesteps_.G()[kk1],gk_all_timesteps_.G()[kk1],0,0,gk_all_timesteps_.G()[qk1],gk_all_timesteps_.G()[qk1],1,1);
		cntr::Bubble2(tstp,ptmp,1,0,gk_all_timesteps_.G()[kk1],gk_all_timesteps_.G()[kk1],1,1,gk_all_timesteps_.G()[qk1],gk_all_timesteps_.G()[qk1],0,0);
		cntr::Bubble2(tstp,ptmp,1,1,gk_all_timesteps_.G()[kk1],gk_all_timesteps_.G()[kk1],1,0,gk_all_timesteps_.G()[qk1],gk_all_timesteps_.G()[qk1],0,1);

		
		// increment P += -wt*Ptmp
		P.incr_timestep(tstp,ptmp,wk);
		//cdmatrix out;
		//P.get_mat(0,out);
		//std::cout << "bubble 3 " << kk << " " << out  << std::endl;
		//std::cout << " -------------------- " << std::endl;
	}
}

// \Sig_{a1,a2}(k,t,t')=-i\sum_q [gamma_pp_{a1,a2_}(q+k,t,t')- gamma_pp_{a1_,a2_}(q+k,t,t')]G_{a2_,a1_}(q,t',t)
// where gamma_pp_{a1,a2}(q,t,t')=W_{a1,a2}(q,t,t')

template <class LATTICE>
void mpi_lattice_step_2b<LATTICE>::get_Sigma_2b(int tstp,int kk,GREEN &S,int iter){
  assert(tstp==gk_all_timesteps_.tstp());
  assert(tstp==wk_all_timesteps_.tstp());
	GREEN_TSTP stmp(tstp,this->ntau_,this->nrpa_,-1);
	GREEN_TSTP PHtmp(tstp,this->ntau_,this->nrpa_,-1);
 	S.set_timestep_zero(tstp);

	for(int q=0;q<this->latt_.nk_bz_;q++){
		double wk=this->latt_.kweight_bz_[q]; // factor 2: ksum normalized for RBZ
		int kq=this->latt_.add_kpoints(kk,1,q,1);
		int kq1=this->latt_.representative_kk(kq);
		int q1=this->latt_.representative_kk(q);
		stmp.clear();

		// Electronic 2nd born
		cntr::Bubble1(tstp,stmp,0,0,wk_all_timesteps_.G()[kq1],wk_all_timesteps_.G()[kq1],0,1,gk_all_timesteps_.G()[q1],gk_all_timesteps_.G()[q1],1,1);
		cntr::Bubble1(tstp,stmp,1,0,wk_all_timesteps_.G()[kq1],wk_all_timesteps_.G()[kq1],1,1,gk_all_timesteps_.G()[q1],gk_all_timesteps_.G()[q1],0,1);
		cntr::Bubble1(tstp,stmp,0,1,wk_all_timesteps_.G()[kq1],wk_all_timesteps_.G()[kq1],0,0,gk_all_timesteps_.G()[q1],gk_all_timesteps_.G()[q1],1,0);
		cntr::Bubble1(tstp,stmp,1,1,wk_all_timesteps_.G()[kq1],wk_all_timesteps_.G()[kq1],1,0,gk_all_timesteps_.G()[q1],gk_all_timesteps_.G()[q1],0,0);

		S.incr_timestep(tstp,stmp,-wk);

		cntr::Bubble1(tstp,stmp,0,0,wk_all_timesteps_.G()[kq1],wk_all_timesteps_.G()[kq1],1,1,gk_all_timesteps_.G()[q1],gk_all_timesteps_.G()[q1],1,1);
		cntr::Bubble1(tstp,stmp,1,0,wk_all_timesteps_.G()[kq1],wk_all_timesteps_.G()[kq1],0,1,gk_all_timesteps_.G()[q1],gk_all_timesteps_.G()[q1],0,1);
		cntr::Bubble1(tstp,stmp,0,1,wk_all_timesteps_.G()[kq1],wk_all_timesteps_.G()[kq1],1,0,gk_all_timesteps_.G()[q1],gk_all_timesteps_.G()[q1],1,0);
		cntr::Bubble1(tstp,stmp,1,1,wk_all_timesteps_.G()[kq1],wk_all_timesteps_.G()[kq1],0,0,gk_all_timesteps_.G()[q1],gk_all_timesteps_.G()[q1],0,0);

		S.incr_timestep(tstp,stmp,wk);
		
		// cdmatrix tmp;
		// PHtmp.get_mat(ntau_,tmp);
		// std::cout << "Stmp " << tmp << std::endl;

		// TODO ADD phonon self energy if needed
		// cdmatrix tmp;

		// Add Holstein like bath
		this->phonon_.rpa(tstp,PHtmp,gk_all_timesteps_.G()[kq1],dk_all_timesteps_.G()[q1]);
		// cntr::Bubble2(tstp,PHtmp,0,0,gk_all_timesteps_.G()[kq1],gk_all_timesteps_.G()[kq1],0,0,dk_all_timesteps_.G()[q1],dk_all_timesteps_.G()[q1],0,0);
		// cntr::Bubble2(tstp,PHtmp,1,1,gk_all_timesteps_.G()[kq1],gk_all_timesteps_.G()[kq1],1,1,dk_all_timesteps_.G()[q1],dk_all_timesteps_.G()[q1],0,0);
		
		// cdmatrix tmp;
		// PHtmp.get_mat(ntau_/2,tmp);
		// std::cout << "PHtmp " << tmp << std::endl;
		S.incr_timestep(tstp,PHtmp,wk);
	}
	// std::cout << "Sigma RPA 6 " <<std::endl;
}

// template <class LATTICE>
// double mpi_lattice_step<LATTICE>::extrapolate(int tstp,int kt);

template <class LATTICE>
void mpi_lattice_step_2b<LATTICE>::get_Gloc(int tstp){
	// collect the local Greenfunction at rank with id tid
	GREEN_TSTP gtmp(tstp,this->ntau_,this->nrpa_,-1);
	for(int kk=0;kk<this->latt_.nk_bz_;kk++){
		double wt=this->latt_.kweight_bz_[kk];
		int kk1=this->latt_.representative_kk(kk);
		if(this->tid_map_[kk1]==this->tid_){
			GREEN_TVIEW tview(tstp,gtmp);
			tview.incr_timestep(green_k_[kk1].G_,CPLX(wt,0.0));
		}
	}
	// add all to task tid
	#if FLEX_USE_MPI==1
		gtmp.Reduce_timestep(this->tid_root_);
	#endif
	if(this->tid_==this->tid_root_){
		cdmatrix tmp;
		Gloc_.set_timestep(tstp,gtmp);
		// Gloc.density_matrix(tstp,tmp);
		// rho_loc_.set_value(tstp,tmp);	
	} 
}

template <class LATTICE>
void mpi_lattice_step_2b<LATTICE>::get_Wloc(int tstp){
	// collect the local Greenfunction at rank with id tid
	GREEN_TSTP gtmp(tstp,this->ntau_,this->nrpa_,1);
	for(int kk=0;kk<this->latt_.nk_bz_;kk++){
		double wt=this->latt_.kweight_bz_[kk];
		int kk1=this->latt_.representative_kk(kk);
		if(this->tid_map_[kk1]==this->tid_){
			GREEN_TVIEW tview(tstp,gtmp);
			tview.incr_timestep(green_k_[kk1].W_,CPLX(wt,0.0));
		}
	}
	// add all to task tid
	#if FLEX_USE_MPI==1
		gtmp.Reduce_timestep(this->tid_root_);
	#endif
	if(this->tid_==this->tid_root_){
		cdmatrix tmp;
		Wloc_.set_timestep(tstp,gtmp);
		// Gloc.density_matrix(tstp,tmp);
		// rho_loc_.set_value(tstp,tmp);	
	} 
}


template <class LATTICE>
void mpi_lattice_step_2b<LATTICE>::get_Dloc(int tstp){
	// collect the local Greenfunction at rank with id tid
	// std::cout << "Dloc 1" << std::endl;
	GREEN_TSTP gtmp(tstp,ntau_,1);
	// std::cout << "Dloc 2" << std::endl;
	for(int kk=0;kk<this->latt_.nk_bz_;kk++){
		double wt=this->latt_.kweight_bz_[kk];
		// std::cout << "Dloc 3" << std::endl;
		int kk1=this->latt_.representative_kk(kk);
		
		if(this->tid_map_[kk1]==this->tid_){
			GREEN_TVIEW tview(tstp,gtmp);
			// std::cout << "Dloc 4 " <<  gtmp.tstp_ << " " <<  green_k_[kk1].D_.ntau() << " " << green_k_[kk1].D_.size1() << " " << green_k_[kk1].D_.size2() << std::endl;
			// std::cout << "Dloc 5 " <<  tview.tstp_ << " " <<  tview.ntau_ << " " << tview.size1_ << " " << tview.size2_ << std::endl;
			tview.incr_timestep(green_k_[kk1].D_,CPLX(wt,0.0));
		}
		// std::cout << "Dloc 6" << std::endl;
	}
	// add all to task tid
	#if FLEX_USE_MPI==1
		gtmp.Reduce_timestep(this->tid_root_);
	#endif
	if(this->tid_==this->tid_root_){
		cdmatrix tmp;
		Dloc_.set_timestep(tstp,gtmp);
		// Gloc.density_matrix(tstp,tmp);
		// rho_loc_.set_value(tstp,tmp);	
	} 
}


template <class LATTICE>
double mpi_lattice_step_2b<LATTICE>::get_eneRPA(int tstp,int kt){
	cdmatrix tmp(this->nrpa_,this->nrpa_);
	double eRPAtmp=0.0,eRPA=0.0;
	
	for(int kk=0;kk<this->latt_.nk_;kk++){
		double wt=this->latt_.kweight_bz_[kk];
		if(this->tid_map_[kk]==this->tid_){
			cntr::convolution_density_matrix(tstp,tmp,green_k_[kk].Sigma_,green_k_[kk].G_,integration::I<double>(kt),this->beta_,this->h_);
			eRPAtmp+=std::real(tmp.trace())*wt;
		}
	}
	std::cout << " ene  " <<eRPAtmp << std::endl;
	eRPA=eRPAtmp;
	// add all to task tid
	#if FLEX_USE_MPI==1
		eRPA=0.0;
		MPI_Allreduce(&eRPAtmp,&eRPA,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		// MPI::COMM_WORLD.Allreduce(&eRPAtmp,&eRPA,1,MPI::DOUBLE,MPI_SUM);
		std::cout << "Fin ene " << tstp << " "  << eRPAtmp << " " << eRPA  << std::endl;
	#endif
	return eRPA;
}