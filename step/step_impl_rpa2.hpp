#pragma once

#include "step.hpp"

template <class LATTICE,class PHONON>
mpi_lattice_step_RPA2<LATTICE,PHONON>::mpi_lattice_step_RPA2(int nt,int ntau,int size,double beta,double h,bool use_omp_for_vie2,double mu,double epsilon,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,std::vector<double> &g,std::vector<double> &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,int nk,int kt,bool test):
	base_type::mpi_lattice_step(nt,ntau,size,beta,h,use_omp_for_vie2,mu,epsilon,tt,U,V,g,omega0,xi,delta,A,dA,v01,nk,kt){
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
	CNTR_ASSERT_EQ(ASSERT_0,nt,latt_.nt_,__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ(ASSERT_0,0,ntau,__PRETTY_FUNCTION__)
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
	this->tid_=gk_all_timesteps_.tid_; 
	this->ntasks_=gk_all_timesteps_.ntasks_;
	this->tid_map_=gk_all_timesteps_.data_.tid_map_;
	this->tid_root_=0;
	if(this->tid_==this->tid_root_){
		Gloc_=GREEN(nt_,ntau_,nrpa_,-1);
		Wloc_=GREEN(nt_,ntau_,nrpa_,1);
	} 
	// using OMP in the solution of vie2 at each kpoint?
	use_omp_for_vie2_=use_omp_for_vie2;
	/////////////////////////////////////////////////
	// std::cout << "step 5 " << std::endl;
	cntr::function<double> gfunc(nt);
    for(int tstp=-1;tstp<=nt;tstp++){
      cdmatrix tmp(1,1);
      tmp(0,0)=g[tstp+1];
      gfunc.set_value(tstp,tmp);
    }
	for(int k=0;k<nk_;k++){
	  if(this->tid_map_[k]==this->tid_){
		this->green_k_[k]=kpoint_green<LATTICE>(nt,ntau,size,beta,h,latt_.kpoints_[k],latt_,mu,this->epsilon_,omega0,gfunc);
		this->green_k_[k].use_omp(use_omp_for_vie2_);
	  }
	}
	// std::cout << "step 6 " << std::endl;
	init_G_mat_nointeraction();
	// std::cout << "step 7 " << std::endl;
	// gather_kk_observables(-1,kt);
	// std::cout << "step 8 " << std::endl;
}



template <class LATTICE,class PHONON>
void mpi_lattice_step_RPA2<LATTICE,PHONON>::init_G_mat_nointeraction(void){
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			// std::cout << "Init " << k << std::endl;
			this->green_k_[k].init_G_mat_nointeraction(this->latt_,this->density_k_[k],kt_);
		}
	}
}

template <class LATTICE,class PHONON>
void mpi_lattice_step_RPA2<LATTICE,PHONON>::extrapolate_timestep_G(int tstp,int kt){
	CNTR_ASSERT_LESEQ(ASSERT_0,tstp,this->nt_-1,__PRETTY_FUNCTION__)
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

template <class LATTICE,class PHONON>
void mpi_lattice_step_RPA2<LATTICE,PHONON>::extrapolate_timestep_W(int tstp,int kt){
	CNTR_ASSERT_LESEQ(ASSERT_0,tstp,this->nt_-1,__PRETTY_FUNCTION__)
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

template <class LATTICE,class PHONON>
void mpi_lattice_step_RPA2<LATTICE,PHONON>::symmetrise_mat(GREEN &G){
	cdmatrix M,M1;
	for(int i=0;i<this->ntau_;i++){
		G.get_mat(i,M);
		M1=(M+M.adjoint())/2.0;
		G.set_mat(i,M1);
		}
}

// TODO: check if it is better to gather density from extrapolated Gk or from rho. I'd guess rho
// TODO: CLEAN UP declaration of same variables in the base and this class
template <class LATTICE,class PHONON>
double mpi_lattice_step_RPA2<LATTICE,PHONON>::step(int tstp,int iter,int kt){
	CNTR_ASSERT_LESEQ(ASSERT_0,tstp,this->nt_,__PRETTY_FUNCTION__)
	double err;
	int n;
	int n1=(tstp==-1 || tstp>kt ? tstp : 0);
	int n2=(tstp==-1 || tstp>kt ? tstp : kt);
	// std::cout << "Step 1  " << std::endl;
	// std::cout << "STEP: " << tstp << ", n1: " << n1 << ", n2: " << n2 << " " << kt  << " " <<  (tstp==-1 || tstp>kt) << std::endl;
	if(tstp>0 && iter==1){
		this->extrapolate_rho(tstp);
		extrapolate_timestep_G(tstp-1,(kt<tstp-1 ? kt : tstp-1));
		extrapolate_timestep_W(tstp-1,(kt<tstp-1 ? kt : tstp-1));
	}
	// std::cout << "Step 2  " << std::endl;
	// Set local
	cdmatrix rtmp;
	this->set_local(tstp);
	this->rho_loc_.get_value(tstp,rtmp);
	this->phonon_=PHONON(latt_.omega0_[tstp+1],latt_.g_[tstp+1],rtmp);
	// std::cout << "Step 3  " << std::endl;
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
					wk_all_timesteps_.G_[k].write_to_hdf5(nameW.str().c_str(),"G");

					std::ostringstream nameD;
					nameD << "D_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
					dk_all_timesteps_.G_[k].write_to_hdf5(nameD.str().c_str(),"G");
				}
			}
		}
		gather_gk_timestep(n);
		base_type::step_boson(n);
		// std::cout << "step 3 " << n << std::endl;
		for(int k=0;k<this->nk_;k++){
			if(this->tid_map_[k]==this->tid_){
				cdmatrix tmp;
				// std::cout << "step 4 " << n << std::endl;
				base_type::get_Sigma_Hartree(n,k,this->density_k_[k].SHartree_);
				// this->density_k_[k].SHartree_.get_value(n,tmp);
				// std::cout <<  "Sigma H "   << n << " " << this->density_k_[k].kk_ << " " << tmp  << std::endl;				
				// std::cout << "step 5 " << n << std::endl;
				base_type::get_Sigma_Fock(n,k,this->density_k_[k].SFock_);
				// this->density_k_[k].SFock_.get_value(n,tmp);
				// std::cout <<  "Sigma F "   << n << " " << this->density_k_[k].kk_ << " " << tmp << std::endl;
				// std::cout << "step 6 " << n << std::endl;
				get_Sigma_RPA(n,k,green_k_[k].Sigma_,iter);

				if(test_){
					std::ostringstream nameW;
					nameW << "Sigma_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
					green_k_[k].Sigma_.write_to_hdf5(nameW.str().c_str(),"G");	
				}
			}
			// std::cout << "step 6c" <<std::endl;
		}
		// std::cout << "step 6d" <<std::endl;
		
		// std::cout << "step 6e" <<std::endl;
	}
	// std::cout << "step 7 " << n << std::endl;
	
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
	// std::cout << "step 8 " << n << std::endl;
	for(n=n1;n<=n2;n++){
		// std::cout << "Loop Pol " << tstp << " " <<  n << " " << n1 << " " << " " << n2 << std::endl;
		if(n1!=n2) gather_gk_timestep(n); //Gather timesteps of all k dependent green's functions
		// set_density_k(n);
		for(int k=0;k<this->nk_;k++){
			if(this->tid_map_[k]==this->tid_) get_Polarization_Bubble(n,k,green_k_[k].P_);
		}
	}
	// std::cout << "step 9 " << n << std::endl;
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			green_k_[k].step_W(tstp,kt,this->latt_,this->density_k_[k]);

			if(test_){
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
	// std::cout << "step 10 " << n << std::endl;
	//Symmetrize  the solutions
	// if(tstp==-1){
	// 	for(int k=0;k<this->nk_;k++){
	// 		if(this->tid_map_[k]==this->tid_){	
	// 			symmetrise_mat(green_k_[k].G_);
	// 			symmetrise_mat(green_k_[k].Sigma_);
	// 		}
	// 	}
	// }
	

	// std::cout << "step 11 " << n << std::endl;
	this->convergence_error_.mpi_bcast_all();
	err=0.0;
	for(int k=0;k<this->nk_;k++) err+=this->convergence_error_.block(k)[0];
	//Symmetrize  the solutions
	if(tstp==-1){
		for(int k=0;k<this->nk_;k++){
			if(this->tid_map_[k]==this->tid_){
				cdmatrix tmp;	
				symmetrise_mat(green_k_[k].G_);
				symmetrise_mat(green_k_[k].Sigma_);
				symmetrise_mat(green_k_[k].chi_);
				symmetrise_mat(green_k_[k].P_);
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
		get_Gloc(n);
		get_Wloc(n);
	}
	return err;
}


template <class LATTICE,class PHONON>
void mpi_lattice_step_RPA2<LATTICE,PHONON>::gather_gk_timestep(int tstp){
	gk_all_timesteps_.reset_tstp(tstp);
	// read Gk_ to timestep
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			gk_all_timesteps_.G_[k].get_data(green_k_[k].G_);
		} 
	}
	// distribute to all nodes
	gk_all_timesteps_.mpi_bcast_all();
}

template <class LATTICE,class PHONON>
void mpi_lattice_step_RPA2<LATTICE,PHONON>::set_density_k(int tstp){
	cdmatrix tmp(this->nrpa_,this->nrpa_);
	for(int k=0;k<this->nk_;k++){
		gk_all_timesteps_.G_[this->latt_.representative_kk(k)].density_matrix(tstp,tmp);
		this->density_k_[this->latt_.representative_kk(k)].rho_.set_value(tstp,tmp);
	}
}

template <class LATTICE,class PHONON>
void mpi_lattice_step_RPA2<LATTICE,PHONON>::gather_Wk_timestep(int tstp){
	wk_all_timesteps_.reset_tstp(tstp);
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			wk_all_timesteps_.G_[k].get_data(green_k_[k].W_);
		}
	}
	// distribute to all nodes
	wk_all_timesteps_.mpi_bcast_all();
}

template <class LATTICE,class PHONON>
void mpi_lattice_step_RPA2<LATTICE,PHONON>::gather_Dk_timestep(int tstp){
	dk_all_timesteps_.reset_tstp(tstp);
	// std::cout << "Gather dk timestep 1" <<std::endl;
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			// std::cout << "Gather dk timestep 2 " <<std::endl;
			dk_all_timesteps_.G_[k].get_data(green_k_[k].D_);
			// std::cout << "Gather dk timestep 3 " <<std::endl;
		}
	}
	// std::cout << "Gather dk timestep 4 " <<std::endl;
	// distribute to all nodes
	dk_all_timesteps_.mpi_bcast_all();
	// std::cout << "Gather dk timestep 5 " <<std::endl;
}


template <class LATTICE,class PHONON>
void mpi_lattice_step_RPA2<LATTICE,PHONON>::get_Polarization_Bubble(int tstp,int qq,GREEN &P){
	CNTR_ASSERT_EQ(ASSERT_0,tstp,gk_all_timesteps_.tstp_,__PRETTY_FUNCTION__)		
	// compute P_{a,a'}(q;t,t') on a timestep
	int kq,i1,i2,kq1,kk1;
	GREEN_TSTP ptmp(tstp,this->ntau_,this->nrpa_,+1);
	P.set_timestep_zero(tstp);
	for(int kk=0;kk<this->latt_.nk_bz_;kk++){
		double wk=this->latt_.kweight_bz_[kk]; // factor 2: ksum normalized for RBZ
		int gammakq;
		// get bubble 
		// -ii*sum_{j1,j2} G_{(a1),(a2)}(k-q;t,t') G_{(a2),(a1)}(k;t',t)
		//kq=latt_.add_kpoints(1,kk,1,qq);
		kq=this->latt_.add_kpoints(1,kk,1,qq);
		// latt_.add_kpoints(kq, gammakq,1,kk,1,qq);
		kq1=this->latt_.representative_kk(kq);
		kk1=this->latt_.representative_kk(kk);
		ptmp.clear();
		for(i1=0;i1<this->nrpa_;i1++){
			for(i2=0;i2<this->nrpa_;i2++){
				cntr::Bubble1(tstp,ptmp,i1,i2,gk_all_timesteps_.G_[kq1],gk_all_timesteps_.G_[kq1],i1,i2,gk_all_timesteps_.G_[kk1],gk_all_timesteps_.G_[kk1],i2,i1);										  
			}
		}
		// increment P += -wt*Ptmp
		P.incr_timestep(tstp,ptmp,-wk);
	}
}


template <class LATTICE,class PHONON>
void mpi_lattice_step_RPA2<LATTICE,PHONON>::get_Sigma_RPA(int tstp,int kk,GREEN &S,int iter){
	CNTR_ASSERT_EQ(ASSERT_0,tstp,gk_all_timesteps_.tstp_,__PRETTY_FUNCTION__)
	CNTR_ASSERT_EQ(ASSERT_0,tstp,wk_all_timesteps_.tstp_,__PRETTY_FUNCTION__)
	GREEN_TSTP stmp(tstp,this->ntau_,this->nrpa_,-1),PHtmp(tstp,this->ntau_,this->nrpa_,-1);
	// GREEN_TSTP stmp1(tstp,this->ntau_,1,-1),stmp2(tstp,this->ntau_,1,-1);
	// std::cout << "Sigma RPA 1" << this->nrpa_ <<std::endl;
 	S.set_timestep_zero(tstp);
	for(int q=0;q<this->latt_.nk_bz_;q++){
		// std::cout << "Sigma RPA 1a " << this->nrpa_ <<std::endl;
		double wk=this->latt_.kweight_bz_[q]; // factor 2: ksum normalized for RBZ
		// std::cout << "Sigma RPA 1b " << this->nrpa_ <<std::endl;
		int i1,i2,ct1,ct2;
		cdmatrix rtmp,vtmp;
		// std::cout << "Sigma RPA 1c " << this->nrpa_ <<std::endl;
		int kq=this->latt_.add_kpoints(0,kk,-1,q);
		int kq1,q1;
		// std::cout << "Sigma RPA 2 " << q << this->nrpa_ <<std::endl;
		// latt_.add_kpoints(kq,gammakq,1,kk,-1,q);
		kq1=this->latt_.representative_kk(kq);
		q1=this->latt_.representative_kk(q);
		stmp.clear();
		PHtmp.clear();
		// std::cout << "Sigma RPA 3 " <<  q << this->nrpa_ <<std::endl;
		for(i1=0;i1<this->nrpa_;i1++){
			for(i2=0;i2<this->nrpa_;i2++){
				cntr::Bubble2(tstp,stmp,i1,i2,gk_all_timesteps_.G_[kq1],gk_all_timesteps_.G_[kq1],i1,i2,wk_all_timesteps_.G_[q1],wk_all_timesteps_.G_[q1],i1,i2);
				// !! NAPAKA El-ph scattering ti doloci ta clen !! 
				// cntr::Bubble2(tstp,PHtmp,i1,i2,gk_all_timesteps_.G_[kq1],gk_all_timesteps_.G_[kq1],i1,i2,dk_all_timesteps_.G_[q1],dk_all_timesteps_.G_[q1],i1,i2);
			}
		}
		// TODO ADD
		cdmatrix tmp;
		this->phonon_.rpa(tstp,this->nrpa_,PHtmp,gk_all_timesteps_.G_[kq1],dk_all_timesteps_.G_[q1]);
		PHtmp.get_mat(ntau_,tmp);
		// std::cout << "Phtmp " << tmp << std::endl;
		if(test_){
			std::ostringstream nameW;
			nameW << "SigmaPh_" << tstp <<"_" << iter << "_" << q<< "_" << kk <<  ".h5";
			PHtmp.write_to_hdf5(nameW.str().c_str(),"G");	
		}

		// std::cout << "Phtmp  " << tmp  <<std::endl;
		S.incr_timestep(tstp,stmp,wk);
		// TODO ADD
		S.incr_timestep(tstp,PHtmp,wk);
		// TODO ADD
		// std::cout << "Sigma RPA 5 " << q << " " << this->nrpa_ << " " << stmp.size1()  <<" " << S.size1() << " " << S.size2() <<std::endl;
	}
	// std::cout << "Sigma RPA 6 " <<std::endl;
}

// template <class LATTICE,class PHONON>
// double mpi_lattice_step<LATTICE>::extrapolate(int tstp,int kt);

template <class LATTICE,class PHONON>
void mpi_lattice_step_RPA2<LATTICE,PHONON>::get_Gloc(int tstp){
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
		gtmp.MPI_Reduce(this->tid_root_);
	#endif
	if(this->tid_==this->tid_root_){
		cdmatrix tmp;
		Gloc_.set_timestep(tstp,gtmp);
		// Gloc.density_matrix(tstp,tmp);
		// rho_loc_.set_value(tstp,tmp);	
	} 
}

template <class LATTICE,class PHONON>
void mpi_lattice_step_RPA2<LATTICE,PHONON>::get_Wloc(int tstp){
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
		gtmp.MPI_Reduce(this->tid_root_);
	#endif
	if(this->tid_==this->tid_root_){
		cdmatrix tmp;
		Wloc_.set_timestep(tstp,gtmp);
		// Gloc.density_matrix(tstp,tmp);
		// rho_loc_.set_value(tstp,tmp);	
	} 
}



template <class LATTICE,class PHONON>
double mpi_lattice_step_RPA2<LATTICE,PHONON>::get_eneRPA(int tstp,int kt){
	cdmatrix tmp,test;
	tmp.resize(this->nrpa_,this->nrpa_);
	int n2=this->nrpa_*this->nrpa_;
	cdouble * tmp1= new cdouble [n2];
	double eRPAtmp=0.0,eRPA=0.0;
	GREEN SigmaG=GREEN(this->nt_,this->ntau_,this->nrpa_,-1);

	for(int k=0;k<this->latt_.nk_bz_;k++){
		double wt=this->latt_.kweight_bz_[k]; // factor 2: ksum normalized for RBZ
		int kk1=this->latt_.representative_kk(k);
		if(this->tid_map_[kk1]==this->tid_){
			//Evaluate the trace of matrix product
			cntr::convolution_density_matrix(tstp,tmp1,green_k_[kk1].Sigma_,green_k_[kk1].G_,integration::I<double>(kt),this->beta_,this->h_);
			for(int i=0;i<n2;i++){
				tmp(i/this->nrpa_,i%this->nrpa_)=tmp1[i];
			}
			eRPAtmp+=std::real(tmp.trace())*wt;
		}
	}
	eRPA=eRPAtmp;
	// add all to task tid
	#if FLEX_USE_MPI==1
		eRPA=0.0;
		MPI::COMM_WORLD.Allreduce(&eRPAtmp,&eRPA,1,MPI::DOUBLE,MPI_SUM);
	#endif
	return eRPA;
}
