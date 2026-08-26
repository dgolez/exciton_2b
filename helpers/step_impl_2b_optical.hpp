#pragma once

#include "inclusions.hpp"

template <class LATTICE>
mpi_lattice_step_2b_optical<LATTICE>::mpi_lattice_step_2b_optical(parameters &param):
	base_type::mpi_lattice_step_optical(param){
	// This is a different type that the base one(includes the info about the RPA )
	// One can use this to reduce the amount of this-> for repeated variables
	int & nk_=this->nk_;
	int & nt_=this->nt_;
	int & ntau_=this->ntau_;
	int & nrpa_=this->nrpa_;
	test_=param.test;
	LATTICE & latt_=this->latt_;
	use_omp_for_vie2_=param.omp_for_vie2;
	// std::cout << "step 2 " << std::endl;
	green_k_.resize(nk_);
	assert(param.nt==latt_.nt_);
	assert(0<=ntau_);
	// std::cout << "step 3 " << std::endl;
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
		D0loc_=GREEN(nt_,ntau_,1,1);
	}
	// Electron-phonon coupling
	g_=param.g;
	migdal_=param.migdal;
	phonontype_=param.phonontype;

/*  Archived fermionic-reservoir bath construction.
	
	It builds the diagonal hybridization self-energy
	Bath_a(t,t') = g_a(t) g^0_{bath,a}(t,t') g_a(t')
	for a conduction reservoir in [bath_low, bath_high] and a valence
	reservoir in [-bath_high, -bath_low].  It was disabled because Bath_
	was never inserted into the Dyson equation.  To reactivate it, restore
	Bath_ in step_decl.hpp, restore bath_low/bath_high/gC_bath/gV_bath in
	parameters and input parsing, and add Bath_ explicitly to the electronic
	self-energy used by kpoint_green::step_dyson.
	Bath_=GREEN(nt_,ntau_,2,-1);
	GREEN Bathtmp(nt_,ntau_,1,-1);
	CFUNC gtmp(param.nt,1);
	cntr::smooth_box dos1(param.bath_low,param.bath_high,30);
    green_equilibrium(Bathtmp,dos1,param.beta,param.h,param.bath_high*2.0,30,30);
    for(int tstp=-1;tstp<=nt_;tstp++){
    	cdmatrix tmp(1,1);
    	tmp(0,0)=param.gC_bath[tstp+1];
    	gtmp.set_value(tstp,tmp);
    }
    for(int tstp=-1;tstp<=nt_;tstp++){
    	Bathtmp.right_multiply(tstp,gtmp);
    	Bathtmp.left_multiply(tstp,gtmp);
    }
    Bath_.set_matrixelement(1,1,Bathtmp,0,0);

    Bathtmp.clear();
    gtmp.set_zero();
    cntr::smooth_box dos2(-param.bath_high,-param.bath_low,30);
    green_equilibrium(Bathtmp,dos2,param.beta,param.h,-param.bath_high*2.0,30,30);
    for(int tstp=-1;tstp<=nt_;tstp++){
    	cdmatrix tmp(1,1);
    	tmp(0,0)=param.gV_bath[tstp+1];
    	gtmp.set_value(tstp,tmp);
    }
    for(int tstp=-1;tstp<=nt_;tstp++){
    	Bathtmp.right_multiply(tstp,gtmp);
    	Bathtmp.left_multiply(tstp,gtmp);
	}
	Bath_.set_matrixelement(0,0,Bathtmp,0,0);
*/

	// Set holstein like bath
	// std::cout << "step 5" << param.omegaBATH[0] <<  std::endl;
	D0hol_=GREEN(nt_,ntau_,1,1);
	cdmatrix tmpOmega(1,1);
	tmpOmega(0,0)=param.omegaBATH[0];
	cntr::green_single_pole_XX(D0hol_,param.omegaBATH[0],param.beta,param.h);
	CFUNC g2tmp(param.nt,1);
	// std::cout << "step 6" << std::endl;
	for(int tstp=-1;tstp<=nt_;tstp++){
    	cdmatrix tmp(1,1);
    	tmp(0,0)=param.gBATH[tstp+1];
    	g2tmp.set_value(tstp,tmp);
    	// std::cout << "coupl " << tstp << " " << tmp << std::endl;
    }
    // std::cout << "step 7" << std::endl;
    for(int tstp=-1;tstp<=nt_;tstp++){
    	D0hol_.right_multiply(tstp,g2tmp);
    	D0hol_.left_multiply(tstp,g2tmp);
	}

	/////////////////////////////////////////////////
	for(int k=0;k<nk_;k++){
	  if(this->tid_map_[k]==this->tid_){
		this->green_k_[k]=kpoint_green<LATTICE>(nt_,ntau_,param.size,param.beta,param.h,latt_.kpoints_[k],latt_,param.mu,param.omega0,param.g,param.mix,param.migdal);
		this->green_k_[k].use_omp(use_omp_for_vie2_);
	  }
	}
	// std::cout << "step 6 " << std::endl;
	init_G_mat_nointeraction();
	
}



template <class LATTICE>
void mpi_lattice_step_2b_optical<LATTICE>::init_G_mat_nointeraction(void){
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			// std::cout << "Init " << k << std::endl;
			this->green_k_[k].init_G_mat_nointeraction(this->latt_,this->density_k_[k],kt_);
		}
	}
}

template <class LATTICE>
void mpi_lattice_step_2b_optical<LATTICE>::extrapolate_timestep_G(int tstp,int kt){
  // assert(tstp<=this->nt_-1);
	// extrapolate from tstp to tstp+1
	if(tstp<-1) return;
	if(tstp==-1){
		for(int k=0;k<this->nk_;k++){
		  if(this->tid_map_[k]==this->tid_) cntr::set_t0_from_mat(this->green_k_[k].G_);
		}
	}else{
		for(int k=0;k<this->nk_;k++){
	  	if(this->tid_map_[k]==this->tid_) cntr::extrapolate_timestep(tstp,this->green_k_[k].G_,integration::I<double>(kt));
		}	
	}
}

template <class LATTICE>
void mpi_lattice_step_2b_optical<LATTICE>::extrapolate_timestep_W(int tstp,int kt){
  // assert(tstp<=this->nt_-1);
	// extrapolate from tstp to tstp+1
	if(tstp<-1) return;
	if(tstp==-1){
		for(int k=0;k<this->nk_;k++){
		  if(this->tid_map_[k]==this->tid_) cntr::set_t0_from_mat(this->green_k_[k].W_);
		}
	}else{
		for(int k=0;k<this->nk_;k++){
	  	if(this->tid_map_[k]==this->tid_) cntr::extrapolate_timestep(tstp,this->green_k_[k].W_,integration::I<double>(kt));
		}	
	}
}

template <class LATTICE>
void mpi_lattice_step_2b_optical<LATTICE>::extrapolate_timestep_D(int tstp,int kt){
  // assert(tstp<=this->nt_-1);
	// extrapolate from tstp to tstp+1
	if(tstp<-1) return;
	if(tstp==-1){
		for(int k=0;k<this->nk_;k++){
		  if(this->tid_map_[k]==this->tid_) cntr::set_t0_from_mat(this->green_k_[k].D_);
		}
	}else{
		for(int k=0;k<this->nk_;k++){
	  	if(this->tid_map_[k]==this->tid_) cntr::extrapolate_timestep(tstp,this->green_k_[k].D_,integration::I<double>(kt));
		}	
	}
}

template <class LATTICE>
void mpi_lattice_step_2b_optical<LATTICE>::symmetrise_mat(GREEN &G){
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
double mpi_lattice_step_2b_optical<LATTICE>::step(int tstp,int iter,int kt){
	// assert(tstp<=this->nt_-1);
	double err;
	int n;
	// if tstp equals -1 or is bigger than kt, both variables get the value of tstp
	// if not, n1 gets 0 and n2 gets kt
	int n1=(tstp==-1 || tstp>kt ? tstp : 0);
	int n2=(tstp==-1 || tstp>kt ? tstp : kt);
	if(tstp>0 && iter==1){
		this->extrapolate_rho(tstp);
		extrapolate_timestep_G(tstp-1,(kt<tstp-1 ? kt : tstp-1));
		extrapolate_timestep_W(tstp-1,(kt<tstp-1 ? kt : tstp-1));
		if(migdal_==1) extrapolate_timestep_D(tstp-1,(kt<tstp-1 ? kt : tstp-1));
	}
	for(n=n1;n<=n2;n++){
		//if(n1!=n2)
		gather_gk_timestep(n); //Gather timesteps of all k dependent green's functions
		// set_density_k(n);
		for(int k=0;k<this->nk_;k++){
		  if(this->tid_map_[k]==this->tid_){
		    if(test_){
		      std::ostringstream nameG;
		      nameG << "Gk0_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
		      green_k_[k].G_.write_to_hdf5(nameG.str().c_str(),"G");	
		    }
		    // get_Polarization_Bubble(n,k,green_k_[k].P_);
		    // The full 2b evaluates a different "bubble" than GW 
		    get_PP_Bubble(n,k,green_k_[k].P_);
		    if(migdal_){
		    	get_phonon_Bubble(n,k,green_k_[k].Pph_);
		    }
		    // cdmatrix tmp;
		    // green_k_[k].P_.get_mat(0,tmp);
		    // std::cout << "Polar " << this->density_k_[k].kk_ << " " << tmp << std::endl; 
		  }
		}
	}
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			green_k_[k].step_W2b(tstp,kt,this->latt_,this->density_k_[k]);
			// spremeni tu in coupling in print!!!
			green_k_[k].step_D(tstp,kt,this->latt_,this->density_k_[k]);
			if(test_){
			  cdmatrix tmp;
			  green_k_[k].W_.get_mat(0,tmp);
			  std::cout << "W " << this->density_k_[k].kk_ << " " << tmp << std::endl;
				for(int n=n1;n<=n2;n++){
					std::ostringstream nameP;
					nameP << "P_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
					green_k_[k].P_.write_to_hdf5(nameP.str().c_str(),"G");

					std::ostringstream nameW;
					nameW << "W_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
					green_k_[k].W_.write_to_hdf5(nameW.str().c_str(),"G");

					std::ostringstream nameDk;
					nameDk << "Dk_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
					green_k_[k].D_.write_to_hdf5(nameDk.str().c_str(),"G");

					std::ostringstream nameD0;
					nameD0 << "D0_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
					green_k_[k].D0_.write_to_hdf5(nameD0.str().c_str(),"G");

					std::ostringstream namePph;
					namePph << "Pph_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
					green_k_[k].Pph_.write_to_hdf5(namePph.str().c_str(),"G");

					std::ostringstream nameD0Pph;
					nameD0Pph << "D0Pph_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
					green_k_[k].D0Pph_.write_to_hdf5(nameD0Pph.str().c_str(),"G");

					std::ostringstream namePphD0;
					namePphD0 << "PphD0_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
					green_k_[k].PphD0_.write_to_hdf5(namePphD0.str().c_str(),"G");
				}	
			}	
		} 
	}
	
	// Set local
	cdmatrix rtmp;
	this->set_local(tstp);
	this->set_sym(tstp);
	// std::cout << "set order 2 " << std::endl;
	this->set_order(tstp);
	this->phonon_=phonon(latt_.omega0_,latt_.g_,this->rho_loc_,this->param_.phonontype);
	for(n=n1;n<=n2;n++){

		cdmatrix tmp;
		// std::cout << "Loop Sigma " << tstp << " " <<  n << " " << n1 << " " << " " << n2 << std::endl;
		gather_Wk_timestep(n);
		gather_Dk_timestep(n);
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
		gather_gk_timestep(n);
		base_type::step_boson(n);
		for(int k=0;k<this->nk_;k++){
			if(this->tid_map_[k]==this->tid_){
				// std::cout <<"Before h" << std::endl;
				base_type::get_Sigma_Hartree(n,this->density_k_[k].SHartree_);
				this->density_k_[k].SHartree_.get_value(n,tmp);
				// std::cout <<  "SigmaH "   << n << " " << this->density_k_[k].kk_ << " " << tmp << std::endl;
				base_type::get_Sigma_Fock(n,k,this->density_k_[k].SFock_);
				this->density_k_[k].SFock_.get_value(n,tmp);
				// std::cout <<  "SigmaF "   << n << " " << this->density_k_[k].kk_ << " " << tmp << std::endl;
				get_Sigma_2b(n,k,green_k_[k].Sigma_,iter);
				green_k_[k].Sigma_.get_mat(0,tmp);
				// std::cout <<  "SigmaRPA "   << n << " " << this->density_k_[k].kk_ << " " << tmp << std::endl;
				if(test_){
					std::ostringstream nameW;
					nameW << "Sigma_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
					green_k_[k].Sigma_.write_to_hdf5(nameW.str().c_str(),"G");	
				}
			}
		}
		// std::cout << "step 6e" <<std::endl;
	}
	// std::cout << "STEP12: " << tstp  << std::endl;
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
	// Later remove since it's doubling of the communication
	for(n=n1;n<=n2;n++){
		gather_gk_timestep(n); //Gather timesteps of all k dependent green's functions
		set_density_k(n);
		base_type::set_local(n);
		base_type::set_sym(n);
		base_type::set_order(n);
		get_Gloc(n);
		get_Wloc(n);
		get_Dloc(n);
	}
	return err;
}


template <class LATTICE>
void mpi_lattice_step_2b_optical<LATTICE>::gather_gk_timestep(int tstp){
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
void mpi_lattice_step_2b_optical<LATTICE>::set_density_k(int tstp){
	cdmatrix tmp(this->nrpa_,this->nrpa_);
	for(int k=0;k<this->nk_;k++){
	  gk_all_timesteps_.G()[k].density_matrix(tstp,tmp);
		this->density_k_[k].rho_.set_value(tstp,tmp);
	}
}

template <class LATTICE>
void mpi_lattice_step_2b_optical<LATTICE>::gather_Wk_timestep(int tstp){
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
void mpi_lattice_step_2b_optical<LATTICE>::gather_Dk_timestep(int tstp){
	dk_all_timesteps_.reset_tstp(tstp);
	// std::cout << "Gather dk timestep 1" <<std::endl;
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			// std::cout << "Gather dk timestep 2 " <<std::endl;
		  dk_all_timesteps_.G()[k].get_data(green_k_[k].D_);
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
void mpi_lattice_step_2b_optical<LATTICE>::get_PP_Bubble(int tstp,int qq,GREEN &P){
  assert(tstp==gk_all_timesteps_.tstp());	
	// compute P_{a,a'}(q;t,t') on a timestep
	int qk,i1,i2,qk1,kk1;
	GREEN_TSTP ptmp(tstp,this->ntau_,this->nrpa_,+1);
	P.set_timestep_zero(tstp);
	for(int kk=0;kk<this->latt_.nk_;kk++){
		double wk=this->latt_.kweight_bz_[kk]; // factor 2: ksum normalized for RBZ
		int gammakq;
		// get bubble for 2b approximation
		// 
		//kq=latt_.add_kpoints(1,kk,1,qq);
		qk=this->latt_.add_kpoints(kk,-1,qq,1);
		//latt_.add_kpoints(kq, gammakq,1,kk,1,qq);
		//qk1=this->latt_.representative_kk(qk);
		//kk1=this->latt_.representative_kk(kk);
		ptmp.clear();

		//Note:BUBBLE2 IS C_{c1,c2}(t1,t2) = ii * A_{a1,a2}(t1,t2) * B_{b1,b2}(t1,t2) 
		cntr::Bubble2(tstp,ptmp,0,0,gk_all_timesteps_.G()[kk],gk_all_timesteps_.G()[kk],0,1,gk_all_timesteps_.G()[qk],gk_all_timesteps_.G()[qk],1,0);
		cntr::Bubble2(tstp,ptmp,0,1,gk_all_timesteps_.G()[kk],gk_all_timesteps_.G()[kk],0,0,gk_all_timesteps_.G()[qk],gk_all_timesteps_.G()[qk],1,1);
		cntr::Bubble2(tstp,ptmp,1,0,gk_all_timesteps_.G()[kk],gk_all_timesteps_.G()[kk],1,1,gk_all_timesteps_.G()[qk],gk_all_timesteps_.G()[qk],0,0);
		cntr::Bubble2(tstp,ptmp,1,1,gk_all_timesteps_.G()[kk],gk_all_timesteps_.G()[kk],1,0,gk_all_timesteps_.G()[qk],gk_all_timesteps_.G()[qk],0,1);

		
		//The sign is different for the pp channel and ph channel
		P.incr_timestep(tstp,ptmp,wk);
		//cdmatrix out;
		//P.get_mat(0,out);
		//std::cout << "bubble 3 " << kk << " " << out  << std::endl;
		//std::cout << " -------------------- " << std::endl;
	}
}

//get P_ph_q (a1,a2) = sum_k ii*G_k(a1,a2_;t,t')G_{q-k}(a1_,a2;t,t')

template <class LATTICE>
void mpi_lattice_step_2b_optical<LATTICE>::get_phonon_Bubble(int tstp,int qq,GREEN &P){
  assert(tstp==gk_all_timesteps_.tstp());	
	// compute P_{a,a'}(q;t,t') on a timestep
	int qk,i1,i2,qk1,kk1;
	GREEN_TSTP ptmp(tstp,this->ntau_,1,+1);
	P.set_timestep_zero(tstp);
	for(int kk=0;kk<this->latt_.nk_;kk++){
		double wk=this->latt_.kweight_bz_[kk]; // factor 2: ksum normalized for RBZ
		int gammakq;
		// get bubble for 2b approximation
		// 
		//kq=latt_.add_kpoints(1,kk,1,qq);
		qk=this->latt_.add_kpoints(kk,-1,qq,1);
		//latt_.add_kpoints(kq, gammakq,1,kk,1,qq);
		//qk1=this->latt_.representative_kk(qk);
		//kk1=this->latt_.representative_kk(kk);
		ptmp.clear();

		//Note:BUBBLE2 IS C_{c1,c2}(t1,t2) = ii * A_{a1,a2}(t1,t2) * B_{b1,b2}(t1,t2) 
		if(phonontype_==0){ // Holstein
			cntr::Bubble2(tstp,ptmp,0,0,gk_all_timesteps_.G()[kk],gk_all_timesteps_.G()[kk],0,0,gk_all_timesteps_.G()[qk],gk_all_timesteps_.G()[qk],0,0);
			P.incr_timestep(tstp,ptmp,wk);
			cntr::Bubble2(tstp,ptmp,0,0,gk_all_timesteps_.G()[kk],gk_all_timesteps_.G()[kk],0,1,gk_all_timesteps_.G()[qk],gk_all_timesteps_.G()[qk],1,0);
			P.incr_timestep(tstp,ptmp,wk);
			cntr::Bubble2(tstp,ptmp,0,0,gk_all_timesteps_.G()[kk],gk_all_timesteps_.G()[kk],1,1,gk_all_timesteps_.G()[qk],gk_all_timesteps_.G()[qk],1,1);
			P.incr_timestep(tstp,ptmp,wk);
			cntr::Bubble2(tstp,ptmp,0,0,gk_all_timesteps_.G()[kk],gk_all_timesteps_.G()[kk],1,0,gk_all_timesteps_.G()[qk],gk_all_timesteps_.G()[qk],0,1);
			P.incr_timestep(tstp,ptmp,wk);
		}else if(phonontype_==1){ // Coupling to c^{\dagger} c - real
			cntr::Bubble2(tstp,ptmp,0,0,gk_all_timesteps_.G()[kk],gk_all_timesteps_.G()[kk],0,0,gk_all_timesteps_.G()[qk],gk_all_timesteps_.G()[qk],1,1);
			P.incr_timestep(tstp,ptmp,wk);
			cntr::Bubble2(tstp,ptmp,0,0,gk_all_timesteps_.G()[kk],gk_all_timesteps_.G()[kk],0,1,gk_all_timesteps_.G()[qk],gk_all_timesteps_.G()[qk],1,0);
			P.incr_timestep(tstp,ptmp,wk);
			cntr::Bubble2(tstp,ptmp,0,0,gk_all_timesteps_.G()[kk],gk_all_timesteps_.G()[kk],1,1,gk_all_timesteps_.G()[qk],gk_all_timesteps_.G()[qk],0,0);
			P.incr_timestep(tstp,ptmp,wk);
			cntr::Bubble2(tstp,ptmp,0,0,gk_all_timesteps_.G()[kk],gk_all_timesteps_.G()[kk],1,0,gk_all_timesteps_.G()[qk],gk_all_timesteps_.G()[qk],0,1);
			P.incr_timestep(tstp,ptmp,wk);
		}else if(phonontype_==3){ // Coupling to n_1-n_0
			cntr::Bubble2(tstp,ptmp,0,0,gk_all_timesteps_.G()[kk],gk_all_timesteps_.G()[kk],0,0,gk_all_timesteps_.G()[qk],gk_all_timesteps_.G()[qk],0,0);
			P.incr_timestep(tstp,ptmp,wk);
			cntr::Bubble2(tstp,ptmp,0,0,gk_all_timesteps_.G()[kk],gk_all_timesteps_.G()[kk],0,1,gk_all_timesteps_.G()[qk],gk_all_timesteps_.G()[qk],1,0);
			P.incr_timestep(tstp,ptmp,-wk);
			cntr::Bubble2(tstp,ptmp,0,0,gk_all_timesteps_.G()[kk],gk_all_timesteps_.G()[kk],1,0,gk_all_timesteps_.G()[qk],gk_all_timesteps_.G()[qk],0,1);
			P.incr_timestep(tstp,ptmp,-wk);
			cntr::Bubble2(tstp,ptmp,0,0,gk_all_timesteps_.G()[kk],gk_all_timesteps_.G()[kk],1,1,gk_all_timesteps_.G()[qk],gk_all_timesteps_.G()[qk],1,1);
			P.incr_timestep(tstp,ptmp,wk);
		}
		
		// increment P += -wt*Ptmp
		
		//cdmatrix out;
		//P.get_mat(0,out);
		//std::cout << "bubble 3 " << kk << " " << out  << std::endl;
		//std::cout << " -------------------- " << std::endl;
	}
	P.right_multiply(tstp,g_);
	P.left_multiply(tstp,g_);
}

// \Sig_{a1,a2}(k,t,t')=-i\sum_q [gamma_pp_{a1,a2_}(q+k,t,t')- gamma_pp_{a1_,a2_}(q+k,t,t')]G_{a2_,a1_}(q,t',t)
// where gamma_pp_{a1,a2}(q,t,t')=W_{a1,a2}(q,t,t')

template <class LATTICE>
void mpi_lattice_step_2b_optical<LATTICE>::get_Sigma_2b(int tstp,int kk,GREEN &S,int iter){
  assert(tstp==gk_all_timesteps_.tstp());
  assert(tstp==wk_all_timesteps_.tstp());
	GREEN_TSTP stmp(tstp,this->ntau_,this->nrpa_,-1);
	GREEN_TSTP PHtmp(tstp,this->ntau_,this->nrpa_,-1);
	GREEN_TSTP HOLtmp(tstp,this->ntau_,this->nrpa_,-1);
 	S.set_timestep_zero(tstp);

	for(int q=0;q<this->latt_.nk_;q++){
		double wk=this->latt_.kweight_bz_[q]; // factor 2: ksum normalized for RBZ
		int kq=this->latt_.add_kpoints(kk,1,q,1);
		//int kq1=this->latt_.representative_kk(kq);
		//int q1=this->latt_.representative_kk(q);
		stmp.clear();
		cntr::Bubble1(tstp,stmp,0,0,wk_all_timesteps_.G()[kq],wk_all_timesteps_.G()[kq],0,1,gk_all_timesteps_.G()[q],gk_all_timesteps_.G()[q],1,1);
		cntr::Bubble1(tstp,stmp,1,0,wk_all_timesteps_.G()[kq],wk_all_timesteps_.G()[kq],1,1,gk_all_timesteps_.G()[q],gk_all_timesteps_.G()[q],0,1);
		cntr::Bubble1(tstp,stmp,0,1,wk_all_timesteps_.G()[kq],wk_all_timesteps_.G()[kq],0,0,gk_all_timesteps_.G()[q],gk_all_timesteps_.G()[q],1,0);
		cntr::Bubble1(tstp,stmp,1,1,wk_all_timesteps_.G()[kq],wk_all_timesteps_.G()[kq],1,0,gk_all_timesteps_.G()[q],gk_all_timesteps_.G()[q],0,0);

		S.incr_timestep(tstp,stmp,-wk);

		cntr::Bubble1(tstp,stmp,0,0,wk_all_timesteps_.G()[kq],wk_all_timesteps_.G()[kq],1,1,gk_all_timesteps_.G()[q],gk_all_timesteps_.G()[q],1,1);
		cntr::Bubble1(tstp,stmp,1,0,wk_all_timesteps_.G()[kq],wk_all_timesteps_.G()[kq],0,1,gk_all_timesteps_.G()[q],gk_all_timesteps_.G()[q],0,1);
		cntr::Bubble1(tstp,stmp,0,1,wk_all_timesteps_.G()[kq],wk_all_timesteps_.G()[kq],1,0,gk_all_timesteps_.G()[q],gk_all_timesteps_.G()[q],1,0);
		cntr::Bubble1(tstp,stmp,1,1,wk_all_timesteps_.G()[kq],wk_all_timesteps_.G()[kq],0,0,gk_all_timesteps_.G()[q],gk_all_timesteps_.G()[q],0,0);

		S.incr_timestep(tstp,stmp,wk);

		// Add phonon self energy
		PHtmp.clear();
		this->phonon_.rpa(tstp,PHtmp,gk_all_timesteps_.G()[kq],dk_all_timesteps_.G()[q]);
		S.incr_timestep(tstp,PHtmp,wk);

		// // Add holstein like bath
		HOLtmp.clear();
		cntr::Bubble2(tstp,HOLtmp,0,0,gk_all_timesteps_.G()[q],gk_all_timesteps_.G()[q],0,0,D0hol_,D0hol_,0,0);
   	cntr::Bubble2(tstp,HOLtmp,0,1,gk_all_timesteps_.G()[q],gk_all_timesteps_.G()[q],0,1,D0hol_,D0hol_,0,0);
    cntr::Bubble2(tstp,HOLtmp,1,0,gk_all_timesteps_.G()[q],gk_all_timesteps_.G()[q],1,0,D0hol_,D0hol_,0,0);
    cntr::Bubble2(tstp,HOLtmp,1,1,gk_all_timesteps_.G()[q],gk_all_timesteps_.G()[q],1,1,D0hol_,D0hol_,0,0);
		S.incr_timestep(tstp,HOLtmp,wk);
	}
}

// template <class LATTICE>
// double mpi_lattice_step<LATTICE>::extrapolate(int tstp,int kt);

template <class LATTICE>
void mpi_lattice_step_2b_optical<LATTICE>::get_Gloc(int tstp){
	// collect the local Greenfunction at rank with id tid
	GREEN_TSTP gtmp(tstp,this->ntau_,this->nrpa_,-1);
	for(int kk=0;kk<this->latt_.nk_;kk++){
		double wt=this->latt_.kweight_bz_[kk];
		//int kk1=this->latt_.representative_kk(kk);
		if(this->tid_map_[kk]==this->tid_){
			GREEN_TVIEW tview(tstp,gtmp);
			tview.incr_timestep(green_k_[kk].G_,CPLX(wt,0.0));
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
void mpi_lattice_step_2b_optical<LATTICE>::get_Wloc(int tstp){
	// collect the local Greenfunction at rank with id tid
	GREEN_TSTP gtmp(tstp,this->ntau_,this->nrpa_,1);
	for(int kk=0;kk<this->latt_.nk_;kk++){
		double wt=this->latt_.kweight_bz_[kk];
		//int kk1=this->latt_.representative_kk(kk);
		if(this->tid_map_[kk]==this->tid_){
			GREEN_TVIEW tview(tstp,gtmp);
			tview.incr_timestep(green_k_[kk].W_,CPLX(wt,0.0));
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
void mpi_lattice_step_2b_optical<LATTICE>::get_Dloc(int tstp){
	// collect the local Greenfunction at rank with id tid
	GREEN_TSTP gtmp(tstp,this->ntau_,1,1),g0tmp(tstp,this->ntau_,1,1);
	for(int kk=0;kk<this->latt_.nk_;kk++){
		double wt=this->latt_.kweight_bz_[kk];
		if(this->tid_map_[kk]==this->tid_){
			GREEN_TVIEW tview(tstp,gtmp),tview0(tstp,g0tmp);
			tview.incr_timestep(green_k_[kk].D_,CPLX(wt,0.0));
			tview0.incr_timestep(green_k_[kk].D0_,CPLX(wt,0.0));
		}
	}
	// add all to task tid
	#if FLEX_USE_MPI==1
		gtmp.Reduce_timestep(this->tid_root_);
		g0tmp.Reduce_timestep(this->tid_root_);
	#endif
	if(this->tid_==this->tid_root_){
		Dloc_.set_timestep(tstp,gtmp);
		D0loc_.set_timestep(tstp,g0tmp);
	}
}



template <class LATTICE>
double mpi_lattice_step_2b_optical<LATTICE>::get_eneRPA(int tstp,int kt){
	cdmatrix tmp(this->nrpa_,this->nrpa_);
	double eRPAtmp=0.0,eRPA=0.0;
	
	for(int kk=0;kk<this->latt_.nk_;kk++){
		double wt=this->latt_.kweight_bz_[kk];
		if(this->tid_map_[kk]==this->tid_){
			cntr::convolution_density_matrix(tstp,tmp,green_k_[kk].Sigma_,green_k_[kk].G_,integration::I<double>(kt),this->beta_,this->h_);
			eRPAtmp+=std::real(tmp.trace())*wt;
		}
	}
	eRPA=eRPAtmp;
	// add all to task tid
	#if FLEX_USE_MPI==1
		// std::cout << "out " << FLEX_USE_MPI << std::endl;
		eRPA=0.0;
		MPI_Allreduce(&eRPAtmp,&eRPA,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		// std::cout << "ene " << tstp << " "  << eRPAtmp << " " << eRPA  << std::endl;
	#endif
	return eRPA;
}
