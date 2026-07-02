#pragma once

#include "inclusions.hpp"
#include <chrono> 
#include "parameters.hpp"
using namespace std::chrono;

template <class LATTICE>
mpi_lattice_step_optical<LATTICE>::mpi_lattice_step_optical(parameters &param){
	param_=param;
	epsilon_=param.epsilon;
	nt_=param.nt;
	ntau_=param.ntau;
	beta_=param.beta;
	h_=param.h;
	norb_=param.size;
	nrpa_=param.size;
	nk_=param.nk;
	kt_=param.kt;
	update_=param.update;
	suscep_=param.suscep;
	eta_=param.eta;
	latt_=LATTICE(param);
	density_k_.resize(nk_);
	vertex_.resize(nk_);
	rho_loc_=CFUNC(nt_,nrpa_);
	rho_loc_.set_zero();
	rho_sym_=CFUNC(nt_,nrpa_);
	rho_sym_.set_zero();
	// std::cout << "set order " << std::endl;
	order_=CFUNC(nt_,1);
	order_.set_zero();
	// Phonon
	X_=CFUNC(nt_,1);
	Pi_=CFUNC(nt_,1);
	X_.set_zero();
	Pi_.set_zero();
	MPI_Comm_rank(MPI_COMM_WORLD,&this->tid_);
	MPI_Comm_size(MPI_COMM_WORLD,&this->ntasks_);

	this->tid_root_=0;

	// Just define tid_map
	cntr::distributed_timestep_array<double>  gk_all_timesteps_=cntr::distributed_timestep_array<double>(nk_,nt_,ntau_,nrpa_,-1,true);
	this->tid_=gk_all_timesteps_.tid(); 
	this->ntasks_=gk_all_timesteps_.ntasks();
	this->tid_map_=gk_all_timesteps_.data().tid_map();
	this->tid_root_=0;
	// // std::cout << "Init MF 1"<< std::endl;
	// /////////////////////////////////////////////////
	for(int k=0;k<nk_;k++){
		density_k_[k]=kpoint_density<LATTICE>(nt_,ntau_,param.size,beta_,h_,latt_.kpoints_[k],latt_,param.mu,param.den,param.epsilon,param.gamma,param.mix,param.phonontype);
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
double mpi_lattice_step_optical<LATTICE>::step(int tstp,int iter,int kt,double om0,double s,double amp){
  assert(tstp<=nt_);
  int n;
  int n1=(tstp==-1 || tstp>kt ? tstp : 0);
  int n2=(tstp==-1 || tstp>kt ? tstp : kt);
  cdmatrix rtmp;
  // O-th order extrapolation TODO: make it higher order
  if(iter==1){
  	extrapolate_rho(tstp);
  	set_local(tstp);
	set_order(tstp);
  	set_sym(tstp);	
  }
  //TODO: USE just one phonon, remove template and add phonontype as parameter
  this->phonon_=phonon(latt_.omega0_,latt_.g_,this->rho_loc_);
  // gather_gk_timestep(n);
  step_boson(tstp);
  for(int k=0;k<nk_;k++){
    cdmatrix tmp;
    get_Sigma_Hartree(tstp,k,density_k_[k].SHartree_);
    density_k_[k].SHartree_.get_value(tstp,tmp);
    // std::cout <<  "SigmaH "   << n << " " << this->density_k_[k].kk_ << " " << tmp << std::endl;
    get_Sigma_Fock(tstp,k,density_k_[k].SFock_);
    density_k_[k].SFock_.get_value(tstp,tmp);
    // std::cout <<  "SigmaF "   << n << " " << this->density_k_[k].kk_ << " " << tmp << std::endl;
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
  set_order(tstp);
  set_sym(tstp);
  // Plot from [-2,2]
  double domega=eta_/20.0;
  int nomega=int(20.0/domega);

  if(tstp==0 && iter==1 && suscep_){
  	auto start = high_resolution_clock::now();
    CFUNC optics(nomega*2+2,4); //This is not the how it should be used by eigen+ vector are tedious
    CFUNC chi(nomega*2+2,4);
    CFUNC seebeck(nomega*2+2,1);
    CFUNC dos(nomega*2+2,1);
    CFUNC omega(nomega*2+2,1);
    cdmatrix tmpO(4,4),tmpC(4,4),tmpS(1,1),om(1,1),tmpA(1,1);
    chi.set_zero();
    optics.set_zero();
    get_seebeck(seebeck,omega,domega,nomega);
    get_dos(dos,domega,nomega);

    // get_optical0(optics,domega,nomega,om0,s,amp);
 //    for(int i=0;i<2*nomega+2;i++){
 //    	tmpC.setZero();
 //    	for(int k=0;k<4;k++){
	// 		for(int l=0;l<4;l++){
	// 			tmpC(l,k)=get_chi0((i-nomega)*domega,k,l,om0,s,amp);
	// 		}
	// 	}
	// 	chi.set_value(i,tmpC);
	// }
    	
    // get_chi0(chi,domega,nomega,om0,s,amp);


    hid_t file_id = open_hdf5_file("optical0.h5");
    hid_t group_id = create_group(file_id, "opt");
    //optics.write_to_hdf5(group_id,"opt");
    chi.write_to_hdf5(group_id,"chi");
    seebeck.write_to_hdf5(group_id,"seebeck");
    dos.write_to_hdf5(group_id,"dos");
    omega.write_to_hdf5(group_id,"om"); 
    close_hdf5_file(file_id);
    //abort();
    //std::cout << "optical0: " << tmp[0] << std::endl; 
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start); 
	std::cout << "time for susceptibility: " << duration.count()/1000000.0 << " s"   << std::endl;
  }
  return err;		
}


template <class LATTICE>
void mpi_lattice_step_optical<LATTICE>::extrapolate_rho(int tstp){
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
void mpi_lattice_step_optical<LATTICE>::step_boson(int tstp){
	// Set phonons
	cdmatrix Xtmp(1,1),Ptmp(1,1),rtmp(nrpa_,nrpa_);
	// rho_loc_.get_value(tstp,rtmp);
	if(tstp==-1){
		// std::cout << "Bos -1 a" << latt_.omega0_[tstp+1] << " " <<  latt_.g_[tstp+1] << " " << rtmp <<std::endl;
		// PHONON ho(latt_.omega0_[tstp+1],latt_.g_[tstp+1],rtmp);
		this->phonon_.eq(X_,Pi_);
		// Xtmp(0,0)=-4.0*latt_.g_[tstp+1]*rtmp(0,1).real()/latt_.omega0_[tstp+1];
		// X_.set_value(tstp,Xtmp);
		// std::cout << "X distor: " <<tstp << " " <<  Xtmp  <<  std::endl;
		// Ptmp(0,0)=0.0;
		// Pi_.set_value(tstp,Ptmp);
	}else if(tstp==0){
		X_.get_value(-1,Xtmp);
		Pi_.get_value(-1,Ptmp);
		X_.set_value(0,Xtmp);
		Pi_.set_value(0,Ptmp);
	}else{
		// Phonons
		this->phonon_.step(tstp,X_,Pi_,kt_,h_);
		// X_.get_value(tstp,Xtmp);
		// std::cout << "X distor: " <<tstp << " " <<  Xtmp  <<  std::endl;
	}
}

template <class LATTICE>
void mpi_lattice_step_optical<LATTICE>::get_Sigma_Hartree(int tstp,int kk,CFUNC &S){
	cdmatrix stmp(nrpa_,nrpa_),v0;
	if(update_ || tstp<0){
	// kk-independent, nut anyway
	stmp.setZero();
	//// the vertex V_{a1,a2} = vertex(q=0)_{a1,a2}
	vertex_[latt_.G_].get_value(tstp,v0);
	for(int k=0;k<latt_.nk_;k++){
		double wk=latt_.kweight_bz_[k]; 
		int i1,i2,c12,at,ct;
		cdmatrix rtmp;
		// gk_all_timesteps_.G_[latt_.representative_kk(k)].density_matrix(tstp,rtmp);
		density_k_[k].rho_.get_value(tstp,rtmp);
		// std::cout << "hartree inside " << " " << rtmp << std::endl;
		// Electronic
		// for(int i = 0; i < 2; i++) {
  //   		MPI_Barrier(MPI_COMM_WORLD);
  //   		if (i == this->tid_){
  //        		std::cout <<" Inside S hartree " << this->tid_ << " " << kk << " " << latt_.representative_kk(k)  << " " << v0 << " " << std::endl;
		// 		std::cout << rtmp << " " << wk << std::endl;
  //   		}
		// }

		for(int i1=0;i1<nrpa_;i1++){
			for(int itmp=0;itmp<nrpa_;itmp++){
				stmp(i1,i1)+=v0(i1,itmp)*rtmp(itmp,itmp)*wk;
				// std::cout << "Inside hartree " << rtmp << " " << wk <<  " " << v0 << std::endl;
			}
		}
		// stmp*=wk;
	}

	// Electron-phonon interaction
	cdmatrix tmp,tmpP,Vtmp(nrpa_,nrpa_),Xtmp(1,1);
	X_.get_value(tstp,Xtmp);
	// std::cout << "Inside hartree2 " << tstp << " " << Vtmp << " " << Xtmp<< std::endl;
	this->phonon_.hartree(Vtmp,Xtmp);
	// std::cout << "Inside hartree2 " << tstp << " " << Vtmp << " " << Xtmp<< std::endl;
	stmp+=Vtmp;

	S.set_value(tstp,stmp);
	}else{
		S.get_value(-1,stmp);
	}
	S.set_value(tstp,stmp);
	// std::cout << "hartree inside end " << " " << stmp << std::endl;
	// std::cout << "Hartree 3 " << std::endl;
}

template <class LATTICE>
void mpi_lattice_step_optical<LATTICE>::set_local(int tstp){
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
void mpi_lattice_step_optical<LATTICE>::set_order(int tstp){
	cdmatrix loc(nrpa_,nrpa_),tmp(nrpa_,nrpa_),rtmp(nrpa_,nrpa_);
	std::complex<double> I(0.0,1.0);
	// Local density matrix
	// std::cout << "set order 3 " << std::endl;
	loc.setZero();
	for(int k=0;k<latt_.nk_;k++){
		double wk=latt_.kweight_bz_[k];
		density_k_[k].rho_.get_value(tstp,rtmp);
		loc+=rtmp*wk*(1.0+exp(I*density_k_[k].kk_));
	}
	order_.set_value(tstp,loc);
}



template <class LATTICE>
void mpi_lattice_step_optical<LATTICE>::set_sym(int tstp){
	cdmatrix loc(nrpa_,nrpa_),tmp(nrpa_,nrpa_),rtmp(nrpa_,nrpa_);
	// Local density matrix
	loc.setZero();
	double kk;
	for(int k=0;k<latt_.nk_;k++){
		double wk=latt_.kweight_bz_[k];
		density_k_[k].rho_.get_value(tstp,rtmp);
		// kk=density_k_[k].kk_;
		loc+=rtmp*wk*cos(density_k_[k].kk_);
	}
	rho_sym_.set_value(tstp,loc);
}



template <class LATTICE>
void mpi_lattice_step_optical<LATTICE>::get_Sigma_Fock(int tstp,int kk,CFUNC &S){
	cdmatrix stmp(nrpa_,nrpa_),inter_tmp(nrpa_,nrpa_);
	if(update_ || tstp<0){
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
		density_k_[kq].rho_.get_value(tstp,rtmp);
		// kk_functions_[latt_.representative_kk(kq)].rho_.get_value(tstp,rtmp);		
		vertex_[q].get_value(tstp,vtmp);
		// std::cout <<" Inside S fock " << this->tid_ << " " << kk << " " << latt_.representative_kk(q)  << " " << vtmp << std::endl;
		// std::cout <<  " " << rtmp << " " << wk << std::endl;
		// std::cout << "q " <<  rtmp << " " << vtmp << std::endl;
		for(int i1=0;i1<nrpa_;i1++){
			for(int i2=0;i2<nrpa_;i2++){
				stmp(i1,i2)-= vtmp(i1,i2)*rtmp(i1,i2)*wk;
				// std::cout << "Inside Fock " << vtmp(i1,i2)<<std::endl;
				// std::cout << rtmp(i1,i2) << " " <<wk << std::endl;
			}
		}
		// std::cout << "stmp q " <<  stmp << " " << wk << std::endl;
		// std::cout <<" ------------- " << std::endl;
	}
	inter_tmp=stmp;
	// std::cout << "Fock " << kk << " " << stmp << std::endl;
	// std::cout << "After " << tstp << stmp-inter_tmp  << std::endl;
	// std::cout << "After " << tstp << stmp  << std::endl;
	}else{
		S.get_value(-1,stmp);	
	}
	S.set_value(tstp,stmp);
}

template <class LATTICE>
double mpi_lattice_step_optical<LATTICE>::get_ekin(int tstp){
	cdmatrix rtmp,hktmp;
	double ekintmp=0.0,ekin=0.0;
	if(tid_==tid_root_){
		for(int k=0;k<latt_.nk_;k++){
			double wt=latt_.kweight_bz_[k];
			density_k_[k].rho_.get_value(tstp,rtmp);
			latt_.hk(hktmp,tstp,density_k_[k].kk_,12);
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
void mpi_lattice_step_optical<LATTICE>::get_ekin(int tstp,cdmatrix &ekin){
	cdmatrix rtmp,hktmp,vec,tmp1,tmp2,tmp3,hkeff;
	ekin.setZero();
	if(tid_==tid_root_){
		for(int k=0;k<latt_.nk_;k++){
			double wt=latt_.kweight_bz_[k];
			density_k_[k].rho_.get_value(tstp,rtmp);
			latt_.hk(hktmp,tstp,density_k_[k].kk_,12);
			density_k_[k].eigen_vec_.get_value(tstp,vec);

			ekin+=(vec.adjoint()*hktmp*rtmp*vec)*wt;
			// testing
			// tmp1=vec.adjoint()*rtmp*vec;
			// tmp2=vec.adjoint()*hktmp*vec;

			// density_k_[k].hkeff_.get_value(tstp,hkeff);
			
			// density_k_[k].rho_diag_.get_value(tstp,tmp3);
			// std::cout << "rho " << tmp1  << std::endl;
			// std::cout << "rhoDIAG " << tmp3  << std::endl;
			// std::cout << "ekeff " << hkeff  << std::endl;
			// std::cout << "ek " << vec.adjoint()*hkeff*vec  << std::endl;
			// std::cout << "vec " << vec  << std::endl;
		}
	}
	#if FLEX_USE_MPI==1
	MPI_Bcast(&ekin,ekin.size(),MPI_COMPLEX,tid_root_,MPI_COMM_WORLD);
	#endif
}


template <class LATTICE>
double mpi_lattice_step_optical<LATTICE>::get_curr(int tstp){
	cdmatrix rtmp,hktmp,Aktmp,rtmp1,rtmp2,drtmp;
	double vktmp1=0.0,vktmp2=0.0,vk=0.0;

	if(tid_==tid_root_){
	  if(tstp<=0){
	    vk=0.0;
	  }else{
		for(int k=0;k<latt_.nk_;k++){

			double wt=latt_.kweight_bz_[k]; // factor 2: ksum normalized for RBZ
			double sign=(k<=this->latt_.G_ ? -1.0 : 1.0 ); // Parity sym has -1
			// Intraband v_k{\alpha\alpa} \rho_{k,\alpha\alpga}
			density_k_[k].rho_.get_value(tstp,rtmp);
			latt_.vk(hktmp,tstp,latt_.kpoints_[k]);
			vktmp1+=std::real((hktmp*rtmp).trace())*wt;
			// Interband
			latt_.Ak(Aktmp,tstp,latt_.kpoints_[k]);			
			//std::cout << "Ak " << latt_.kpoints_[k] << " " << density_k_[latt_.representative_kk(kk1)].kk_ << " " << Aktmp << std::endl;
			//std::cout << "----------- " << std::endl;
			// Derivative 
			density_k_[k].rho_.get_value(tstp-1,rtmp1);
			density_k_[k].rho_.get_value(tstp,rtmp2);
			drtmp=(rtmp2-rtmp1)/(h_); //todo implement higher order differenciator for function
			//std::cout << "densities " << rtmp1 << " " << rtmp2 << std::endl;
			vktmp2+=std::real((Aktmp*drtmp).trace())*wt;
			
			//std::cout << "momentum " <<Aktmp << " " << drtmp << " " << vktmp2 << std::endl;
		}
		// std::cout << "velocity " << vktmp1 << " " << vktmp2 << std::endl;
		vk=vktmp1+vktmp2;
	  }
	}
	// add all to task tid
	#if FLEX_USE_MPI==1
	    MPI_Bcast(&vk,1,MPI_DOUBLE,tid_root_,MPI_COMM_WORLD);
	#endif
	return vk;
}


template <class LATTICE>
double mpi_lattice_step_optical<LATTICE>::get_curr_peierls(int tstp){
	cdmatrix rtmp,hktmp,Aktmp,rtmp1,rtmp2,drtmp;
	double vktmp1=0.0,vktmp2=0.0,vk=0.0;

	if(tid_==tid_root_){
		for(int k=0;k<latt_.nk_;k++){
			double wt=latt_.kweight_bz_[k]; // factor 2: ksum normalized for RBZ
			double sign=(k<=this->latt_.G_ ? -1.0 : 1.0 ); // Parity sym has -1
			
			// Intraband v_k{\alpha\alpa} \rho_{k,\alpha\alpga}
			density_k_[k].rho_.get_value(tstp,rtmp);
			latt_.vk(hktmp,tstp,latt_.kpoints_[k]);
			//std::cout << "Velocities " << latt_.kpoints_[k] << " " << hktmp << std::endl;
			//latt_.vk(hktmp,tstp,density_k_[latt_.representative_kk(kk1)].kk_);
			vktmp1+=std::real((hktmp*rtmp).trace())*wt;
		}
		//std::cout << "velocity " <<tstp << " "  << vktmp1 << std::endl;
		vk=vktmp1;
	}
	// add all to task tid
	#if FLEX_USE_MPI==1
	    MPI_Bcast(&vk,1,MPI_DOUBLE,tid_root_,MPI_COMM_WORLD);
	#endif
		//std::cout << "velocity FIN: " <<tstp << " "  << vk << std::endl;
	return vk;
}


template <class LATTICE>
double mpi_lattice_step_optical<LATTICE>::get_curr_dip(int tstp){
	cdmatrix rtmp,hktmp,Aktmp,rtmp1,rtmp2,drtmp;
	double vktmp1=0.0,vktmp2=0.0,vk=0.0;

	if(tid_==tid_root_){
	  if(tstp<=0){
	    vk=0.0;
	  }else{
	    for(int k=0;k<latt_.nk_;k++){

			double wt=latt_.kweight_bz_[k]; // factor 2: ksum normalized for RBZ
			// Interband
			latt_.Ak(Aktmp,tstp,latt_.kpoints_[k]);
			// Derivative 
			density_k_[k].rho_.get_value(tstp-1,rtmp1);
			density_k_[k].rho_.get_value(tstp,rtmp2);
			drtmp=(rtmp2-rtmp1)/(h_); //todo implement higher order differenciator for function
			vktmp2+=std::real((Aktmp*drtmp).trace())*wt;
		}
		// std::cout << "velocity " << vktmp1 << " " << vktmp2 << std::endl;
		vk=vktmp2;
	  }
	}
	// add all to task tid
	#if FLEX_USE_MPI==1
	    MPI_Bcast(&vk,1,MPI_DOUBLE,tid_root_,MPI_COMM_WORLD);
	#endif
	return vk;
}


template <class LATTICE>
double mpi_lattice_step_optical<LATTICE>::get_dip(int tstp){
	cdmatrix rtmp,hktmp,Aktmp;
	double vktmp1=0.0,vktmp2=0.0,vk=0.0;

	if(tid_==tid_root_){
		for(int k=0;k<latt_.nk_;k++){

			double wt=latt_.kweight_bz_[k]; // factor 2: ksum normalized for RBZ
			double sign=(k<=this->latt_.G_ ? -1.0 : 1.0 ); // Parity sym has -1
			// Interband
			latt_.Ak(Aktmp,tstp,density_k_[k].kk_);
			// Derivative 
			density_k_[k].rho_.get_value(tstp,rtmp);
			vktmp2+=std::real((Aktmp*rtmp).trace())*wt;
		}
		// std::cout << "velocity " << vktmp1 << " " << vktmp2 << std::endl;
		vk=vktmp2;
	}
	// add all to task tid
	#if FLEX_USE_MPI==1

	    MPI_Bcast(&vk,1,MPI_DOUBLE,tid_root_,MPI_COMM_WORLD);

	#endif
	return vk;
}


// template <class LATTICE>
// std::complex<double> mpi_lattice_step_optical<LATTICE>::get_chi0(CFUNC &chi,double domega,int nomega,double om0,double s,double amp){
//   cdmatrix rtmp(2,2),tmp(2,2);
//   std::complex<double> optical0=std::complex<double>(0.0,0.0);
//   std::vector<cdmatrix> sigma(4);
//   double om;
//   sigma[0].resize(2,2);sigma[1].resize(2,2);sigma[2].resize(2,2);sigma[3].resize(2,2);
//   sigma[0].setZero();sigma[1].setZero();sigma[2].setZero();sigma[3].setZero();
  
//   sigma[0](0,0)=1.0;
//   sigma[0](1,1)=1.0;
//   sigma[1](0,1)=1.0;
//   sigma[1](1,0)=1.0;
//   sigma[2](0,1)=std::complex<double>(0.0,-1.0);
//   sigma[2](1,0)=std::complex<double>(0.0,1.0);
//   sigma[3](0,0)=1.0;
//   sigma[3](1,1)=-1.0;

//   cdmatrix previous(4,4);

//   dvector ek,bk(3);
//   cdmatrix vec,wp(2,2),wm(2,2),vk,Ak;
//   for(int k=0;k<latt_.nk_;k++){
//       double wt=latt_.kweight_bz_[k]; // factor 2: ksum normalized for RBZ
		  
//       // Eigenvalues and structure factors
//       density_k_[k].hkeff_.get_value(-1,rtmp);
//       // std::cout << "Kpoint:  " << k << " " << latt_.kpoints_[k] << " " << rtmp << std::endl;
//       Eigen::SelfAdjointEigenSolver<cdmatrix> eigensolver(rtmp);
//       ek=eigensolver.eigenvalues();
//       vec=eigensolver.eigenvectors();
            
//       bk(0)=std::real(rtmp(0,1));
//       bk(1)=std::imag(rtmp(0,1));
//       bk(2)=std::real(rtmp(1,1));
//       wp=0.5*(sigma[0]+(bk(0)*sigma[1]+bk(1)*sigma[2]+bk(2)*sigma[3])/bk.norm());
//       wm=0.5*(sigma[0]-(bk(0)*sigma[1]+bk(1)*sigma[2]+bk(2)*sigma[3])/bk.norm());
      
//       for(int io=0;io<2*nomega+2;io++){
//     	om=(io-nomega)*domega;
//     	chi.get_value(io,previous);
//     	double f1=cntr::fermi<double>(beta_,ek(1))+amp*gauss(ek(1),om0,s)-amp*gauss(ek(1),-om0,s);
//       	double f0=cntr::fermi<double>(beta_,ek(0))+amp*gauss(ek(0),om0,s)-amp*gauss(ek(0),-om0,s);

//       	for(int mu=0;mu<4;mu++){
// 			for(int nu=0;nu<4;nu++){ 
// 				//Sum +-
// 				// std::cout << "0 : " << previous(mu,nu) << std::endl; 
//     	  		tmp=sigma[nu]*wp*sigma[mu]*wm;
//       			previous(mu,nu)+=wt*tmp.trace()*(f0-f1)/(om-(ek(1)-ek(0))+std::complex<double>(0,eta_));
//       			// std::cout << "1 " << mu << " "<< nu << " " << wt <<" " << tmp.trace() << " " <<  (f0-f1) << " " << previous(mu,nu)  << std::endl;
//       			//Sum -+
//       			// std::cout << "2 : " << previous(mu,nu) << std::endl;
//       			tmp=sigma[nu]*wm*sigma[mu]*wp;
//       			previous(mu,nu)-=wt*tmp.trace()*(f0-f1)/(om-(ek(0)-ek(1))+std::complex<double>(0,eta_));
//       			// std::cout << "3 " << mu << " "<< nu << " " << wt <<" " << tmp.trace() << " " <<  (f0-f1) << " " << previous(mu,nu)  << std::endl;
//       			// std::cout <<"--------------------------" << std::endl;
//       		} //MU loop
//       	} //NU loop
//       	chi.set_value(io,previous);
//       	// std::cout <<"--------------------------" << std::endl;
//       } // Omega loop 
//   }//k loop
// }

template <class LATTICE>
std::complex<double> mpi_lattice_step_optical<LATTICE>::get_chi0(double omega,int mu,int nu,double om0,double s,double amp){
  cdmatrix rtmp(2,2),tmp(2,2);
  std::complex<double> optical0=std::complex<double>(0.0,0.0);
  std::vector<cdmatrix> sigma(4);
  double eta=0.01;
  sigma[0].resize(2,2);sigma[1].resize(2,2);sigma[2].resize(2,2);sigma[3].resize(2,2);
  sigma[0].setZero();sigma[1].setZero();sigma[2].setZero();sigma[3].setZero();
  
  sigma[0](0,0)=1.0;
  sigma[0](1,1)=1.0;
  sigma[1](0,1)=1.0;
  sigma[1](1,0)=1.0;
  sigma[2](0,1)=std::complex<double>(0.0,-1.0);
  sigma[2](1,0)=std::complex<double>(0.0,1.0);
  sigma[3](0,0)=1.0;
  sigma[3](1,1)=-1.0;
  //std::cout << "Sigma0 " << sigma[0] << std::endl; 
  //std::cout << "Sigmax " << sigma[1] << std::endl; 
  //std::cout << "Sigmay " << sigma[2] << std::endl;
  //std::cout << "Sigmaz " << sigma[3] << std::endl;
  if(tid_==tid_root_){
    dvector ek,bk(3);
    cdmatrix vec,wp(2,2),wm(2,2),vk,Ak;
    for(int k=0;k<latt_.nk_;k++){
	    
      double wt=latt_.kweight_bz_[k]; // factor 2: ksum normalized for RBZ
		  
      // Eigenvalues and structure factors
      density_k_[k].hkeff_.get_value(-1,rtmp);
      //std::cout << "Kpoint:  " << k << " " << latt_.kpoints_[k] << " " << rtmp << std::endl;
      Eigen::SelfAdjointEigenSolver<cdmatrix> eigensolver(rtmp);
      ek=eigensolver.eigenvalues();
      vec=eigensolver.eigenvectors();
            
      bk(0)=std::real(rtmp(0,1));
      bk(1)=std::imag(rtmp(0,1));
      bk(2)=std::real(rtmp(1,1));
      wp=0.5*(sigma[0]+(bk(0)*sigma[1]+bk(1)*sigma[2]+bk(2)*sigma[3])/bk.norm());
      wm=0.5*(sigma[0]-(bk(0)*sigma[1]+bk(1)*sigma[2]+bk(2)*sigma[3])/bk.norm());
      double f1=cntr::fermi<double>(beta_,ek(1))+amp*gauss(ek(1),om0,s)-amp*gauss(ek(1),-om0,s);
      double f0=cntr::fermi<double>(beta_,ek(0))+amp*gauss(ek(0),om0,s)-amp*gauss(ek(0),-om0,s);
      //Sum +-
      tmp=sigma[nu]*wp*sigma[mu]*wm;
      //std::cout << "Values chi  " << " " <<  nu << " " << mu << " " << tmp.trace() << std::endl;
      //std::cout << "11 :  " <<  bk(2)*bk(2)/(bk.norm()*bk.norm()) << std::endl;
      //std::cout << "33 :  " <<  bk(0)*bk(0)/(bk.norm()*bk.norm()) << std::endl;
      
      optical0+=wt*tmp.trace()*(f0-f1)/(omega-(ek(1)-ek(0))+std::complex<double>(0,eta));
      
      
      //Sum -+
      tmp=sigma[nu]*wm*sigma[mu]*wp;
      //std::cout << "Values chi  " << " " <<  nu << " " << mu << " " << tmp.trace() << std::endl;
      //std::cout << "11 " << bk(2)*bk(2)/(bk.norm()*bk.norm())  << std::endl;
      //std::cout << "33 :  " << bk(0)*bk(0)/(bk.norm()*bk.norm()) << std::endl;
      
      optical0-=wt*tmp.trace()*(f0-f1)/(omega-(ek(0)-ek(1))+std::complex<double>(0,eta));
      
      //std::cout << "---------------------" << std::endl;
      //std::cout << "---------------------" << std::endl;
    }
  
  }
  return optical0;
  // add all to task tid
  //#if FLEX_USE_MPI==1
  //	MPI_Bcast(vk.data(),1,MPI_DOUBLE,vk.data(),tid_root_);
  //#endif
  //return vk;
}


template <class LATTICE>
std::complex<double> mpi_lattice_step_optical<LATTICE>::get_optical0(CFUNC &optics,double domega,int nomega,double om0,double s,double amp){
  cdmatrix rtmp(2,2),tmp(2,2);
  std::complex<double> optical0=std::complex<double>(0.0,0.0);
  double om;
  std::vector<cdmatrix> sigma(4);
  sigma[0].resize(2,2);sigma[1].resize(2,2);sigma[2].resize(2,2);sigma[3].resize(2,2);
  sigma[0].setZero();sigma[1].setZero();sigma[2].setZero();sigma[3].setZero();
  
  sigma[0](0,0)=1.0;
  sigma[0](1,1)=1.0;
  sigma[1](0,1)=1.0;
  sigma[1](1,0)=1.0;
  sigma[2](0,1)=std::complex<double>(0.0,-1.0);
  sigma[2](1,0)=std::complex<double>(0.0,1.0);
  sigma[3](0,0)=1.0;
  sigma[3](1,1)=-1.0;

  dvector ek,bk(3);
  cdmatrix vec,wp(2,2),wm(2,2),vk,Ak,vk1;
  cdmatrix previous(4,4);

  for(int k=0;k<latt_.nk_;k++){

    double wt=latt_.kweight_bz_[k];
    // Eigenvalues and structure factors
    density_k_[k].hkeff_.get_value(-1,rtmp);
    Eigen::SelfAdjointEigenSolver<cdmatrix> eigensolver(rtmp);
    ek=eigensolver.eigenvalues();
    vec=eigensolver.eigenvectors();
    bk(0)=std::real(rtmp(0,1));
    bk(1)=std::imag(rtmp(0,1));
    bk(2)=std::real(rtmp(1,1));
    wp=0.5*(sigma[0]+(bk(0)*sigma[1]+bk(1)*sigma[2]+bk(2)*sigma[3])/bk.norm());
    wm=0.5*(sigma[0]-(bk(0)*sigma[1]+bk(1)*sigma[2]+bk(2)*sigma[3])/bk.norm());
	//velocities
    latt_.vkFULL(vk,-1,latt_.kpoints_[k]);
    latt_.Ak(Ak,-1,latt_.kpoints_[k]);
	double vnu,vmu;

	for(int io=0;io<2*nomega+2;io++){
    	om=(io-nomega)*domega;
    	optics.get_value(io,previous);
		for(int mu=0;mu<4;mu++){
			for(int nu=0;nu<4;nu++){ 
				// nu direction 
	      		if(nu==0){
					vnu=0.0;
	      		}else if(nu==1){
					vnu=std::real(Ak(0,1));
	      		}else if(nu==2){
					vnu=std::imag(Ak(0,1));
	      		}else if(nu==3){
					vnu=std::real(vk(1,1));
	      		}
	      		// mu direction
	      		if(mu==0){
					vmu=0.0;
	      		}else if(mu==1){
					vmu=std::real(Ak(0,1));
	      		}else if(mu==2){
					vmu=std::imag(Ak(0,1));
	     	 	}else if(mu==3){
				vmu=std::real(vk(1,1));
	      		}

	      		double f1=cntr::fermi<double>(beta_,ek(1))+amp*gauss(ek(1),om0,s)-amp*gauss(ek(1),-om0,s);
	      		double f0=cntr::fermi<double>(beta_,ek(0))+amp*gauss(ek(0),om0,s)-amp*gauss(ek(0),-om0,s);
	      		//Sum +-
				tmp=vnu*vmu*sigma[nu]*wp*sigma[mu]*wm;
	      		previous(mu,nu)+=wt*tmp.trace()*(f0-f1)/(om-(ek(1)-ek(0))+std::complex<double>(0,eta_));

	      		//Sum -+
      			tmp=vnu*vmu*sigma[nu]*wm*sigma[mu]*wp;
      			previous(mu,nu)-=wt*tmp.trace()*(f0-f1)/(om-(ek(0)-ek(1))+std::complex<double>(0,eta_));
      		} //NU LOOP
	    } //MU LOOP
	    optics.set_value(io,previous);
	} //OMEGA LOOP
  }//K LOOP
}


template <class LATTICE>
std::complex<double> mpi_lattice_step_optical<LATTICE>::get_seebeck(CFUNC &seebeck,CFUNC &omega,double domega,int nomega){
  std::complex<double> phi=0.0;
  double om;
  cdmatrix previous(1,1);
  dvector ek(2);
  cdmatrix vec(2,2),G(2,2),diatmp(2,2),vk(2,2),tmp(2,2),rtmp(2,2);
  for(int k=0;k<latt_.nk_;k++){
    double wt=latt_.kweight_bz_[k];
    // Eigenvalues and structure factors
    density_k_[k].hkeff_.get_value(-1,rtmp);
    Eigen::SelfAdjointEigenSolver<cdmatrix> eigensolver(rtmp);
    ek=eigensolver.eigenvalues();
    vec=eigensolver.eigenvectors();
    // Scan through omega
    for(int io=0;io<2*nomega+2;io++){
      om=(io-nomega)*domega;
      for(int i=0;i<nrpa_;i++) diatmp(i,i)=1.0/(om-ek(i)+std::complex<double>(0,eta_));
 	  G=vec*diatmp*vec.adjoint();
  	  G=(-1.0)*G.imag()/3.14159265358979323846;
  	  latt_.vkFULL(vk,-1,latt_.kpoints_[k]);
  	  tmp=vk*G*vk*G;
  	  seebeck.get_value(io,previous);
  	  previous(0,0)-=tmp.trace();
  	  seebeck.set_value(io,previous);
  	  previous(0,0)=om; //Recycle
  	  omega.set_value(io,previous);
 	}
  }
  return phi;
}

template <class LATTICE>
std::complex<double> mpi_lattice_step_optical<LATTICE>::get_dos(CFUNC &dos,double domega,int nomega){
  std::complex<double> phi=0.0;
  cdmatrix previous(1,1);
  double omega;
  dvector ek(2);
  cdmatrix vec(2,2),G(2,2),diatmp(2,2),vk(2,2),tmp(2,2),rtmp(2,2);
  for(int k=0;k<latt_.nk_;k++){
    double wt=latt_.kweight_bz_[k];
    // Eigenvalues and structure factors
    density_k_[k].hkeff_.get_value(-1,rtmp);
    Eigen::SelfAdjointEigenSolver<cdmatrix> eigensolver(rtmp);
    ek=eigensolver.eigenvalues();
    vec=eigensolver.eigenvectors();
    // Scan through omega
    for(int io=0;io<2*nomega+2;io++){
      omega=(io-nomega)*domega;
      for(int i=0;i<nrpa_;i++) diatmp(i,i)=1.0/(omega-ek(i)+std::complex<double>(0,eta_));
 	  G=vec*diatmp*vec.adjoint();
  	  G=(-1.0)*G.imag()/3.14159265358979323846;
  	  dos.get_value(io,previous);
  	  // std::cout << "------------- " << std::endl;
  	  // std::cout << "previous " << previous << std::endl;
  	  previous(0,0)+=G.trace();
  	  // std::cout << "next " << previous << " " << G.trace() << std::endl;
  	  // std::cout << "------------- " << std::endl;
  	  dos.set_value(io,previous);
  	}
  }
  //return phi;
}

#define EXPMAX 100
template <class LATTICE>
double mpi_lattice_step_optical<LATTICE>::dfermi(double omega,double beta){
	double arg = omega * beta;
    if (fabs(arg) > EXPMAX) {
        return 0;
    } else {
        return -1.0*cntr::fermi<double>(beta_,omega)*cntr::fermi<double>(beta_,omega)*exp(arg)*beta;
    }
}


template <class LATTICE>
double mpi_lattice_step_optical<LATTICE>::get_seebeck_boltzmann(double beta){
	double nom=0,mu;
	double denom=0;
	int nk=latt_.nk_;
	double dk=2.0*PI/nk;
	// std::cout << "dk " <<dk << std::endl;

	cdmatrix rtmp1,rtmp2,drtmp,rtmp;
	for(int k=0;k<nk;k++){
    	double wt=latt_.kweight_bz_[k];
    	// Make a derivative in k of eigenvalue
    	density_k_[k].hkeff_eigen_.get_value(-1,rtmp);
    	density_k_[(k+1)%nk].hkeff_eigen_.get_value(-1,rtmp1);
    	density_k_[(k-1+nk)%nk].hkeff_eigen_.get_value(-1,rtmp2);
    	mu=density_k_[k].mu_;
    	drtmp=(rtmp2-rtmp1)/(2.0*dk);
    	// std::cout << "vk " << density_k_[k].kk_ << " " <<  drtmp(0,0) << " " << rtmp(0,0) << " " << dfermi(std::real(rtmp(0,0)-mu),beta)  << " " <<  drtmp(1,1) << " " << rtmp(1,1) << " " << dfermi(std::real(rtmp(1,1)-mu),beta) <<  std::endl;
    	nom+=wt*std::real(drtmp(0,0)*drtmp(0,0)*dfermi(std::real(rtmp(0,0)-mu),beta)*(rtmp(0,0)-mu) + drtmp(1,1)*drtmp(1,1)*dfermi(std::real(rtmp(1,1)-mu),beta)*(rtmp(1,1)-mu));
    	denom+=wt*std::real(drtmp(0,0)*drtmp(0,0)*dfermi(std::real(rtmp(0,0)-mu),beta) + drtmp(1,1)*drtmp(1,1)*dfermi(std::real(rtmp(1,1)-mu),beta));
    }
    // std::cout << "nom " << nom << std::endl;
    // std::cout << "denom " << denom << std::endl;
	return nom*beta/denom;
}

template <class LATTICE>
double mpi_lattice_step_optical<LATTICE>::get_dAk(int tstp){

	cdmatrix rtmp,hktmp;
	double vktmp=0.0,vk=0.0;
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
double mpi_lattice_step_optical<LATTICE>::gauss(double om,double om0,double s){
	return exp(-(om-om0)*(om-om0)/(2*s*s))/(sqrt(2*PI)*s);
}

template <class LATTICE>
double mpi_lattice_step_optical<LATTICE>::get_eneHF(int tstp){
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
