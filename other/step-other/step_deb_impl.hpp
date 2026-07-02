#pragma once

#include "step.hpp"
#include "../program/parameters.hpp"

template <class LATTICE>
mpi_lattice_deb_step<LATTICE>::mpi_lattice_deb_step(parameters &param){
	epsilon_=param.epsilon;
	nt_=param.nt;
	ntau_=param.ntau;
	beta_=param.beta;
	h_=param.h;
	norb_=param.size;
	nrpa_=param.size;
	nk_=param.nk;
	kt_=param.kt;
	// std::cout << "Init MF"<< std::endl;
	latt_=LATTICE(param);
	density_k_.resize(nk_);
	vertex_.resize(nk_);
	rho_loc_=CFUNC(nt_,nrpa_);
	rho_loc_.set_zero();
	// Phonon
	X_=CFUNC(nt_,1);
	Pi_=CFUNC(nt_,1);
	X_.set_zero();
	Pi_.set_zero();
	MPI_Comm_rank(MPI_COMM_WORLD,&this->tid_);
        //this->tid_=MPI_COMM_WORLD.Get_rank();                                                                                                                                                          
        MPI_Comm_size(MPI_COMM_WORLD,&this->ntasks_);
        //this->ntasks_=MPI_COMM_WORLD.Get_size();  
	//this->tid_=MPI::COMM_WORLD.Get_rank();
	//this->ntasks_=MPI::COMM_WORLD.Get_size();
	this->tid_root_=0;

	// Just define tid_map
	cntr::distributed_timestep_array<double>  gk_all_timesteps_=cntr::distributed_timestep_array<double>(nk_,nt_,ntau_,nrpa_,-1,true);
	this->tid_=gk_all_timesteps_.tid(); 
	this->ntasks_=gk_all_timesteps_.ntasks();
	this->tid_map_=gk_all_timesteps_.data().tid_map();
	this->tid_root_=0;
	// std::cout << "Init MF 1"<< std::endl;
	/////////////////////////////////////////////////
	double gamma=0.0;
	double mix=0.0;
	int phonontype=1;

	for(int k=0;k<nk_;k++){
		density_k_[k]=kpoint_density<LATTICE>(nt_,ntau_,norb_,beta_,h_,latt_.kpoints_[k],latt_,param.mu,param.epsilon,param.gamma,param.mix,param.phonontype);
	  	// set the vertex at *every* kk
	  	vertex_[k]=CFUNC(nt_,nrpa_);
	  	for(int tstp=-1;tstp<=nt_;tstp++){
			cdmatrix vtmp;
			latt_.V(vtmp,tstp,latt_.kpoints_[k]);
			vertex_[k].set_value(tstp,vtmp);
	  	}
	}
	// std::cout << "Init MF 2"<< std::endl;
	for(int k=0;k<nk_;k++) density_k_[k].init_rho_free(latt_);
	// std::cout << "Init MF 3"<< std::endl;
}

template <class LATTICE>
double mpi_lattice_deb_step<LATTICE>::step(int tstp,int iter,int kt,double om0,double s,double amp){
  assert(tstp<=nt_);
  int n;
  int n1=(tstp==-1 || tstp>kt ? tstp : 0);
  int n2=(tstp==-1 || tstp>kt ? tstp : kt);
  cdmatrix rtmp;
  // O-th order extrapolation TODO: make it higher order
  // extrapolate_rho(tstp);
  // set_local(tstp);
  // this->phonon_=PHONON(latt_.omega0_,latt_.g_,this->rho_loc_);
  
  // gather_gk_timestep(n);
  // step_boson(tstp);

  for(int k=0;k<nk_;k++){
    cdmatrix tmp;
    get_Sigma_Hartree(tstp,k,density_k_[k].SHartree_);
    // density_k_[k].SHartree_.get_value(tstp,tmp);
    get_Sigma_Fock(tstp,k,density_k_[k].SFock_);
    // density_k_[k].SFock_.get_value(tstp,tmp);
  }
  double err=0.0;
  for(int k=0;k<nk_;k++) err += density_k_[k].step_dyson_with_error(tstp,iter,latt_,om0,s,amp);
  //Symmetrize  the solutions
  // std::cout <<"Before sym " << err <<std::endl;
  for(int k=0;k<nk_;k++){
    cdmatrix rtmp;
    density_k_[k].rho_.get_value(tstp,rtmp);
    rtmp=(rtmp+rtmp.adjoint())*0.5;
    density_k_[k].rho_.set_value(tstp,rtmp);
  }
  set_local(tstp);
  return err;		
}


template <class LATTICE>
void mpi_lattice_deb_step<LATTICE>::extrapolate_rho(int tstp){
	if(tstp<-1) return;
	if(tstp==0){
		for(int k=0;k<nk_;k++){
			cdmatrix tmp;
			density_k_[k].rho_.get_value(-1,tmp);
			density_k_[k].rho_.set_value(0,tmp);
		}
	}
	// std::cout << "Extrap 5 " << std::endl;
	if(tstp>0){
		for(int k=0;k<nk_;k++){
				cdmatrix tmp;
				density_k_[k].rho_.get_value(tstp-1,tmp);
				density_k_[k].rho_.set_value(tstp,tmp);
		}
	}
}

template <class LATTICE>
void mpi_lattice_deb_step<LATTICE>::step_boson(int tstp){
	// Set phonons
	cdmatrix Xtmp(1,1),Ptmp(1,1),rtmp(nrpa_,nrpa_);
	// rho_loc_.get_value(tstp,rtmp);
	// std::cout << "Bos 1"<<std::endl;
	if(tstp==-1){
		// std::cout << "Bos -1 a" << latt_.omega0_[tstp+1] << " " <<  latt_.g_[tstp+1] << " " << rtmp <<std::endl;
		// PHONON ho(latt_.omega0_[tstp+1],latt_.g_[tstp+1],rtmp);
		// std::cout << "Bos -1 b"<<std::endl;
		this->phonon_.eq(X_,Pi_);
		// std::cout << "Bos -1 c"<<std::endl;
		// Xtmp(0,0)=-4.0*latt_.g_[tstp+1]*rtmp(0,1).real()/latt_.omega0_[tstp+1];
		// X_.set_value(tstp,Xtmp);
		// std::cout << "X distor: " <<tstp << " " <<  Xtmp  <<  std::endl;
		// Ptmp(0,0)=0.0;
		// Pi_.set_value(tstp,Ptmp);
	}else if(tstp==0){
		// std::cout << "Bos 0a"<<std::endl;
		X_.get_value(-1,Xtmp);
		// std::cout << "Bos 0b"<<std::endl;
		Pi_.get_value(-1,Ptmp);
		// std::cout << "Bos 0c"<<std::endl;
		X_.set_value(0,Xtmp);
		// std::cout << "Bos 0d"<<std::endl;
		Pi_.set_value(0,Ptmp);
		// std::cout << "Bos 0e"<<std::endl;
	}else{
		// Phonons
		this->phonon_.step(tstp,X_,Pi_,kt_,h_);
		
	}
}

template <class LATTICE>
void mpi_lattice_deb_step<LATTICE>::get_Sigma_Hartree(int tstp,int kk,CFUNC &S){
	// kk-independent, nut anyway
	cdmatrix stmp(nrpa_,nrpa_),v0;
	stmp.setZero();
	//// the vertex V_{a1,a2} = vertex(q=0)_{a1,a2}
	vertex_[latt_.G_].get_value(tstp,v0);
	for(int k=0;k<latt_.nk_;k++){
		double wk=latt_.kweight_bz_[k]; 
		int i1,i2,c12,at,ct;
		cdmatrix rtmp;
		density_k_[k].rho_.get_value(tstp,rtmp);
		// std::cout << "Hartree " << kk << " " << rtmp << 
		for(int i1=0;i1<nrpa_;i1++){
			for(int itmp=0;itmp<nrpa_;itmp++){
				stmp(i1,i1)+=v0(i1,itmp)*rtmp(itmp,itmp)*wk;
				// std::cout << "Inside hartree " << k << " " << rtmp << " " << wk <<  " " << v0 << " " << rtmp <<  " " << stmp << std::endl;
			}
		}
		// stmp*=wk;
	}

	// Electron-phonon interaction
	cdmatrix tmp,tmpP,Vtmp(nrpa_,nrpa_),Xtmp(1,1);
	X_.get_value(tstp,Xtmp);
	this->phonon_.hartree(Vtmp,Xtmp);
	
	stmp+=Vtmp;
	// std::cout << "Inside hartree2 " << Vtmp << " " << Xtmp << " " << stmp << std::endl;
	S.set_value(tstp,stmp);
}

template <class LATTICE>
void mpi_lattice_deb_step<LATTICE>::set_local(int tstp){
	cdmatrix loc(nrpa_,nrpa_),tmp(nrpa_,nrpa_),rtmp(nrpa_,nrpa_);
	// Local density matrix
	loc.setZero();
	for(int k=0;k<latt_.nk_;k++){
		double wk=latt_.kweight_bz_[k];
		density_k_[k].rho_.get_value(tstp,rtmp);
		loc+=rtmp*wk;
	}
	rho_loc_.set_value(tstp,loc);
}


template <class LATTICE>
void mpi_lattice_deb_step<LATTICE>::get_Sigma_Fock(int tstp,int kk,CFUNC &S){
	cdmatrix stmp(nrpa_,nrpa_),inter_tmp(nrpa_,nrpa_);
	cdmatrix Xtmp(1,1);
	stmp.setZero();
	// Electrons
	for(int q=0;q<latt_.nk_;q++){
		double wk=latt_.kweight_bz_[q];
		int ct;
		cdmatrix rtmp,vtmp;
		int kq,gammakq;
		kq=latt_.add_kpoints(kk,1,q,-1);
		// gk_all_timesteps_.G_[latt_.representative_kk(kq)].density_matrix(tstp,rtmp);
		// density_k_[latt_.representative_kk(kq)].rho_.get_value(tstp,rtmp);
		density_k_[kq].rho_.get_value(tstp,rtmp);
		// kk_functions_[latt_.representative_kk(kq)].rho_.get_value(tstp,rtmp);		
		vertex_[q].get_value(tstp,vtmp);
		// vertex_[latt_.representative_kk(q)].get_value(tstp,vtmp);
		// std::cout <<" Inside S fock " << this->tid_ << " " << kk << " " << latt_.representative_kk(q)  << " " << vtmp << std::endl;
		// std::cout <<  " " << rtmp << " " << wk << std::endl;
		for(int i1=0;i1<nrpa_;i1++){
			for(int i2=0;i2<nrpa_;i2++){
				stmp(i1,i2)-= vtmp(i1,i2)*rtmp(i1,i2)*wk;
				// std::cout << "Inside Fock " << vtmp(i1,i2)<<std::endl;
				// std::cout << rtmp(i1,i2) << " " <<wk << std::endl;
			}
		}
	}
	inter_tmp=stmp;
	// std::cout << "Fock " << kk << " " << stmp << std::endl;
	// std::cout << "After " << tstp << stmp-inter_tmp  << std::endl;
	// std::cout << "After " << tstp << stmp  << std::endl;
	S.set_value(tstp,stmp);
}

template <class LATTICE>
double mpi_lattice_deb_step<LATTICE>::get_ekin(int tstp){
	cdmatrix rtmp,hktmp;
	double ekintmp=0.0,ekin=0.0;
	if(tid_==tid_root_){
		for(int k=0;k<latt_.nk_;k++){
			double wt=latt_.kweight_bz_[k];
			// int kk1=latt_.representative_kk(k);
			density_k_[k].rho_.get_value(tstp,rtmp);
			// density_k_[latt_.representative_kk(kk1)].rho_.get_value(tstp,rtmp);
			latt_.hk(hktmp,tstp,density_k_[k].kk_,12);
			// latt_.hkfree(hktmp,tstp,density_k_[latt_.representative_kk(kk1)].kk_);
			ekintmp+=std::real((hktmp*rtmp).trace())*wt;
		}
		ekin=ekintmp;
	}
	#if FLEX_USE_MPI==1
	    MPI_Bcast(&ekin,1,MPI_DOUBLE,tid_root_,MPI_COMM_WORLD);
	#endif
	return ekin;
}

template <class LATTICE>
void mpi_lattice_deb_step<LATTICE>::get_ekin(int tstp,cdmatrix &ekin){
	cdmatrix rtmp,hktmp,vec,tmp1,tmp2,tmp3,hkeff;
	ekin.setZero();
	if(tid_==tid_root_){
		for(int k=0;k<latt_.nk_;k++){
			double wt=latt_.kweight_bz_[k];
			density_k_[k].rho_.get_value(tstp,rtmp);
			latt_.hk(hktmp,tstp,density_k_[k].kk_,12);
			density_k_[k].eigen_vec_.get_value(tstp,vec);

			ekin+=(vec.adjoint()*hktmp*rtmp*vec)*wt;
		}
	}
	#if FLEX_USE_MPI==1
	    MPI_Bcast(ekin.data(),ekin.size(),MPI_COMPLEX,tid_root_,MPI_COMM_WORLD);
	#endif
}

template <class LATTICE>
double mpi_lattice_deb_step<LATTICE>::get_curr(int tstp){
	cdmatrix rtmp,hktmp;
	double vktmp=0.0,vk=0.0;

	if(tid_==tid_root_){
		for(int k=0;k<latt_.nk_;k++){
			double wt=latt_.kweight_bz_[k]; // factor 2: ksum normalized for RBZ
			density_k_[k].rho_.get_value(tstp,rtmp);
			latt_.vk(hktmp,tstp,latt_.kpoints_[k]);
			vktmp+=std::real((hktmp*rtmp).trace())*wt;

			// double sign=(k<=this->latt_.G_ ? -1.0 : 1.0 ); // Parity sym has -1
			// int kk1=latt_.representative_kk(k);
			// density_k_[latt_.representative_kk(kk1)].rho_.get_value(tstp,rtmp);
			// latt_.vk(hktmp,tstp,density_k_[latt_.representative_kk(kk1)].kk_);
			// vktmp+=std::real((hktmp*rtmp).trace())*wt*sign;
		}
		vk=vktmp;
	}
	// add all to task tid
	#if FLEX_USE_MPI==1
	    MPI_Bcast(&vktmp,1,MPI_DOUBLE,tid_root_,MPI_COMM_WORLD);
	#endif
	return vk;
}

template <class LATTICE>
double mpi_lattice_deb_step<LATTICE>::get_dAk(int tstp){
	// std::cout << "daK 0" << std::endl;
	cdmatrix rtmp,hktmp;
	double vktmp=0.0,vk=0.0;
	// std::cout << "daK 1" << std::endl;
	if(tid_==tid_root_){
				for(int k=0;k<latt_.nk_;k++){
					double wt=latt_.kweight_bz_[k]; // factor 2: ksum normalized for RBZ
					density_k_[k].rho_.get_value(tstp,rtmp);
					// green_k_[kk1].G_.density_matrix(tstp,rtmp);
					latt_.dAk(hktmp,tstp,density_k_[k].kk_);
					vktmp+=std::real((hktmp*rtmp).trace())*wt;
				}
		vk=vktmp;
	}
	// add all to task tid
	#if FLEX_USE_MPI==1
	    MPI_Bcast(&vktmp,1,MPI_DOUBLE,tid_root_,MPI_COMM_WORLD);
	#endif
	return vk;
}

template <class LATTICE>
double mpi_lattice_deb_step<LATTICE>::get_eneHF(int tstp){
	cdmatrix tmpH,tmpF,rtmp;
	double eHFtmp=0.0,eHF=0.0;
	for(int k=0;k<latt_.nk_;k++){
	  double wt=latt_.kweight_bz_[k]; // factor 2: ksum normalized for RBZ
	  
	    density_k_[k].rho_.get_value(tstp,rtmp);
            density_k_[k].SHartree_.get_value(tstp,tmpH);
            density_k_[k].SFock_.get_value(tstp,tmpF);
            eHFtmp+=std::real(((tmpH+tmpF)*rtmp).trace())*wt;
	  
	}
	eHF=eHFtmp;
	// add all to task tid
	#if FLEX_USE_MPI==1
	    eHF=0.0;
	    MPI_Allreduce(&eHFtmp,&eHF,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	#endif
	return eHF;
}
