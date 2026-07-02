#pragma once

#include "step.hpp"

template <class LATTICE>
mpi_lattice_step_2b_cdw<LATTICE>::mpi_lattice_step_2b_cdw(int nt,int ntau,int size,double beta,double h,bool use_omp_for_vie2,double mu,double epsilon,std::vector<double> &tt,std::vector<double> &U,std::vector<double> &V,CFUNC &g,CFUNC &omega0,double xi,std::vector<double> &delta,std::vector<double> &A,std::vector<double> &dA,double v01,int nk,int kt,bool test,double mix):
  base_type::mpi_lattice_step_cdw(nt,ntau,size,beta,h,use_omp_for_vie2,mu,epsilon,tt,U,V,g,omega0,xi,delta,A,dA,v01,nk,kt){
        // This is a different type that the base one(includes the info about the RPA )
	// One can use this to reduce the amount of this-> for repeated variables
	// std::cout << "step 1 " << std::endl;
	int & nk_=this->nk_;
	int & nt_=this->nt_;
	int & ntau_=this->ntau_;
	int & nrpa_=this->nrpa_;
	test_=test;
	LATTICE & latt_=this->latt_;
	green_k_.resize(nk_);
	assert(nt==latt_.nt_);
	assert(0<=ntau_);
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
	/////////////////////////////////////////////////
	// MPI Setup;  gk_..., wk_..., and *this will have the same mpi-layout
	this->tid_=gk_all_timesteps_.tid(); 
	this->ntasks_=gk_all_timesteps_.ntasks();
	this->tid_map_=gk_all_timesteps_.data().tid_map();
	this->tid_root_=0;
	if(this->tid_==this->tid_root_){
		Gloc_=GREEN(nt_,ntau_,nrpa_,-1);
		Wloc_=GREEN(nt_,ntau_,nrpa_,1);
	} 
	// using OMP in the solution of vie2 at each kpoint?
	use_omp_for_vie2_=use_omp_for_vie2;
	/////////////////////////////////////////////////
	for(int k=0;k<nk_;k++){
	  if(this->tid_map_[k]==this->tid_){
		this->green_k_[k]=kpoint_green<LATTICE>(nt,ntau,size,beta,h,latt_.kpoints_[k],latt_,mu,this->epsilon_,omega0,g,mix);
		this->green_k_[k].use_omp(use_omp_for_vie2_);
	  }
	}
	init_G_mat_nointeraction();
	gather_gk_timestep(-1);
	set_density_k(-1);
}



template <class LATTICE>
void mpi_lattice_step_2b_cdw<LATTICE>::init_G_mat_nointeraction(void){
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			// std::cout << "Init " << k << std::endl;
			this->green_k_[k].init_G_mat_nointeraction(this->latt_,this->density_k_[k],kt_);
		}
	}
}

template <class LATTICE>
void mpi_lattice_step_2b_cdw<LATTICE>::extrapolate_timestep_G(int tstp,int kt){
  assert(tstp<=this->nt_);
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
void mpi_lattice_step_2b_cdw<LATTICE>::extrapolate_timestep_W(int tstp,int kt){
  assert(tstp<=this->nt_);
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
void mpi_lattice_step_2b_cdw<LATTICE>::symmetrise_mat(GREEN &G){
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
double mpi_lattice_step_2b_cdw<LATTICE>::step(int tstp,int iter,int kt,double om0,double s,double amp){
  assert(tstp<=this->nt_);
	double err;
	int n;
	int n1=(tstp==-1 || tstp>kt ? tstp : 0);
	int n2=(tstp==-1 || tstp>kt ? tstp : kt);
	if(tstp>0 && iter==1){
		this->extrapolate_rho(tstp);
		extrapolate_timestep_G(tstp-1,(kt<tstp-1 ? kt : tstp-1));
		extrapolate_timestep_W(tstp-1,(kt<tstp-1 ? kt : tstp-1));
	}
	// Set local
	cdmatrix rtmp;
	this->set_local(tstp);
	this->phonon_=phonon(latt_.omega0_,latt_.g_,this->rho_loc_);
	for(n=n1;n<=n2;n++){
		cdmatrix tmp;
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
				cdmatrix tmp;
				base_type::get_Sigma_Hartree(n,k,this->density_k_[k].SHartree_);
				// this->density_k_[k].SHartree_.get_value(n,tmp);
				// std::cout <<  "Sigma H "   << n << " " << this->density_k_[k].kk_ << " " << tmp  << std::endl;			    
				base_type::get_Sigma_Fock(n,k,this->density_k_[k].SFock_);
				// this->density_k_[k].SFock_.get_value(n,tmp);
				// std::cout <<  "Sigma F "   << n << " " << this->density_k_[k].kk_ << " " << tmp << std::endl;
				get_Sigma_RPA_cdw(n,k,green_k_[k].Sigma_,iter);

				if(test_){
					std::ostringstream nameW;
					nameW << "Sigma_" << tstp <<"_" << n <<"_" << iter << "_" << k <<  ".h5";
					green_k_[k].Sigma_.write_to_hdf5(nameW.str().c_str(),"G");	
				}
			}
		}
	}
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
	for(n=n1;n<=n2;n++){
		// std::cout << "Loop Pol " << tstp << " " <<  n << " " << n1 << " " << " " << n2 << std::endl;
		if(n1!=n2) gather_gk_timestep(n); //Gather timesteps of all k dependent green's functions
		// set_density_k(n);
		for(int k=0;k<this->nk_;k++){
			if(this->tid_map_[k]==this->tid_) get_Polarization_Bubble(n,k,green_k_[k].P_);
		}
	}
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			green_k_[k].step_W2b(tstp,kt,this->latt_,this->density_k_[k]);

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


template <class LATTICE>
void mpi_lattice_step_2b_cdw<LATTICE>::gather_gk_timestep(int tstp){
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
void mpi_lattice_step_2b_cdw<LATTICE>::set_density_k(int tstp){
	cdmatrix tmp(this->nrpa_,this->nrpa_);
	for(int k=0;k<this->nk_;k++){
	  gk_all_timesteps_.G()[this->latt_.representative_kk(k)].density_matrix(tstp,tmp);
	  this->density_k_[this->latt_.representative_kk(k)].rho_.set_value(tstp,tmp);
	}
}

template <class LATTICE>
void mpi_lattice_step_2b_cdw<LATTICE>::gather_Wk_timestep(int tstp){
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
void mpi_lattice_step_2b_cdw<LATTICE>::gather_Dk_timestep(int tstp){
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


template <class LATTICE>
void mpi_lattice_step_2b_cdw<LATTICE>::get_Sigma_RPA_cdw(int tstp,int kk,GREEN &S,int iter){
  assert(tstp==gk_all_timesteps_.tstp());
  assert(tstp==wk_all_timesteps_.tstp());
	GREEN_TSTP stmp(tstp,this->ntau_,this->nrpa_,-1),PHtmp(tstp,this->ntau_,this->nrpa_,-1);
	GREEN_TSTP stmp1(tstp,this->ntau_,this->nrpa_,-1),stmp2(tstp,this->ntau_,this->nrpa_,-1);
 	S.set_timestep_zero(tstp);

 	for(int q=0;q<latt_.nk_rbz_;q++){
		double wk=0.5*latt_.kweight_bz_[q]; // factor 2: ksum normalized for RBZ
		int i1,i2,ct1,ct2;
		cdmatrix rtmp,vtmp;
		//int kq=latt_.add_kpoints(0,kk,-1,q);
		int kq,gammakq,kq1,q1;
		latt_.add_kpoints_rbz(kq,gammakq,1,kk,-1,q);
		kq1=latt_.representative_kk(kq);
		q1=latt_.representative_kk(q);
		for(i1=0;i1<nrpa_;i1++){
			for(i2=0;i2<nrpa_;i2++){
				int a1=latt_.idx_a(i1);
				int c1=latt_.idx_c(i1);
				int a2=latt_.idx_a(i2);
				int c2=latt_.idx_c(i2);
				stmp1.clear();
				for(ct1=0;ct1<2;ct1++){
					for(ct2=0;ct2<2;ct2++){
						int m1=latt_.idx(a1,(c1+ct1+gammakq)%2);
						int m2=latt_.idx(a2,(c2+ct2+gammakq)%2);
						int l1=latt_.idx(a1,ct1);
						int l2=latt_.idx(a2,ct2);
						cntr::Bubble2(tstp,stmp2,0,0,gk_all_timesteps_.G()[kq1],gk_all_timesteps_.G()[kq1],m1,m2,wk_all_timesteps_.G()[q1],wk_all_timesteps_.G()[q1],l1,l2);
					    stmp1.incr(stmp2,1.0);
					}
				}
		  		stmp.set_matrixelement(i1,i2,stmp1,0,0);
			}
		}
		S.incr_timestep(tstp,stmp,wk);
	}
}


template <class LATTICE>
void mpi_lattice_step_2b_cdw<LATTICE>::get_Polarization_Bubble(int tstp,int qq,GREEN &P){
   assert(tstp==gk_all_timesteps_.tstp());	
	// compute P_{a,a'}(q;t,t') on a timestep
	int kq,i1,i2,kq1,kk1;
	GREEN_TSTP ptmp(tstp,this->ntau_,this->nrpa_,+1),ptmp1(tstp,this->ntau_,this->nrpa_,+1),ptmp2(tstp,this->ntau_,this->nrpa_,+1);
	P.set_timestep_zero(tstp);
	for(int kk=0;kk<this->latt_.nk_rbz_;kk++){
		double wk=0.5*this->latt_.kweight_bz_[kk]; // factor 2: ksum normalized for RBZ
		int gammakq;
		// get bubble 
		// -ii*sum_{j1,j2} G_{(a1),(a2)}(k-q;t,t') G_{(a2),(a1)}(k;t',t)
		//kq=latt_.add_kpoints(1,kk,1,qq);
		this->latt_.add_kpoints_rbz(kq, gammakq,1,kk,1,qq);
		kq1=this->latt_.representative_kk(kq);
		kk1=this->latt_.representative_kk(kk);

		for(i1=0;i1<nrpa_;i1++){
			for(i2=0;i2<nrpa_;i2++){
				int a1=latt_.idx_a(i1);
				int c1=latt_.idx_c(i1);
				int a2=latt_.idx_a(i2);
				int c2=latt_.idx_c(i2);
				ptmp1.clear();
				for(int j1=0;j1<2;j1++){
					for(int j2=0;j2<2;j2++){
						int m1=latt_.idx(a1,(c1+j1+gammakq)%2);
						int m2=latt_.idx(a2,(c2+j2+gammakq)%2);
						int l1=latt_.idx(a1,j1);
						int l2=latt_.idx(a2,j2);
						cntr::Bubble1(tstp,ptmp2,0,0,gk_all_timesteps_.G()[kq1],gk_all_timesteps_.G()[kq1],m1,m2,
										  gk_all_timesteps_.G()[kk1],gk_all_timesteps_.G()[kk1],l1,l2);
						ptmp1.incr(ptmp2,1.0);
					}
				}
				ptmp.set_matrixelement(i1,i2,ptmp1,0,0);
			}
		}
		P.incr_timestep(tstp,ptmp,-wk);
	}
}

template <class LATTICE>
void mpi_lattice_step_2b_cdw<LATTICE>::get_Gloc(int tstp){
	// collect the local Greenfunction at rank with id tid
	GREEN_TSTP gtmp(tstp,this->ntau_,this->nrpa_,-1);
	for(int kk=0;kk<this->latt_.nk_rbz_;kk++){
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
void mpi_lattice_step_2b_cdw<LATTICE>::get_Wloc(int tstp){
	// collect the local Greenfunction at rank with id tid
	GREEN_TSTP gtmp(tstp,this->ntau_,this->nrpa_,1);
	for(int kk=0;kk<this->latt_.nk_rbz_;kk++){
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
double mpi_lattice_step_2b_cdw<LATTICE>::get_eneRPA(int tstp,int kt){
	cdmatrix tmp(this->nrpa_,this->nrpa_);
	int n2=this->nrpa_*this->nrpa_;
	cdouble * tmp1= new cdouble [n2];
	double eRPAtmp=0.0,eRPA=0.0;
	GREEN SigmaG=GREEN(this->nt_,this->ntau_,this->nrpa_,-1);

	for(int k=0;k<this->latt_.nk_rbz_;k++){
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
		MPI_Allreduce(&eRPAtmp,&eRPA,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		// MPI::COMM_WORLD.Allreduce(&eRPAtmp,&eRPA,1,MPI::DOUBLE,MPI_SUM);
	#endif
	return eRPA;
}
