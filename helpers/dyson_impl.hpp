#pragma once

#include "inclusions.hpp"

template <class LATTICE>
kpoint_density<LATTICE>::kpoint_density(void){
}

template <class LATTICE>
kpoint_density<LATTICE>::kpoint_density(int nt,int ntau,int size,double beta,double h,double kk,LATTICE &latt,double mu,double den,double epsilon,double gamma,double mix,int phonontype):
	beta_(beta), h_(h), nt_(nt), ntau_(ntau),
	nrpa_(size),mu_(mu),den_(den),epsilon_(epsilon),phonontype_(phonontype){
		SHartree_=CFUNC(nt_,nrpa_);
		SFock_=CFUNC(nt_,nrpa_);
		rho_=CFUNC(nt_,nrpa_);
		rho_diag_=CFUNC(nt_,nrpa_);
		rho_eq_=cdmatrix(nrpa_,nrpa_);
		hk_=CFUNC(nt_,nrpa_);
		hkeff_=CFUNC(nt_,nrpa_);
		hkeff_eigen_=CFUNC(nt_,nrpa_);
		vertex_=CFUNC(nt_,nrpa_);
		eigen_vec_=CFUNC(nt_,nrpa_);
		gamma_=gamma;
		mix_= mix;
		use_omp(false);
		kk_=kk;
		for(int tstp=-1;tstp<=nt_;tstp++) set_hk(tstp,1,latt);
		for(int tstp=-1;tstp<=nt_;tstp++) set_vertex(tstp,latt);
	}

template <class LATTICE>
void kpoint_density<LATTICE>::use_omp(bool onoff){
	#if FLEX_CAN_USE_OMP==1
		use_omp_=onoff;
		if(use_omp_){
			std::cout << __PRETTY_FUNCTION__<< ":\n using OMP-parallelizaztion over VIE2 equations ..." << std::endl;
			std::cout << omp_get_max_threads() << " threads" << std::endl;
		}
	#else
		std::cerr << __PRETTY_FUNCTION__<< ":\n OMP not allowed " << std::endl;
		exit(0);
	#endif
}

template <class LATTICE>
double kpoint_density<LATTICE>::gauss(double om,double om0,double s){
	return exp(-(om-om0)*(om-om0)/(2*s*s))/(sqrt(2*PI)*s);
}

template <class LATTICE>
double kpoint_density<LATTICE>::heaviside(double ene){
	if(ene<=0){
		return 1;
	}
	else{
		return 0;	
	}
}


template <class LATTICE>
void kpoint_density<LATTICE>::init_rho_free(LATTICE &latt){
	int tstp=-1;
	cdmatrix hkeff;
	dvector ek;
	cdmatrix vec,fk(nrpa_,nrpa_),rtmp(nrpa_,nrpa_);
	latt.hkfree(hkeff,tstp,kk_);
	Eigen::SelfAdjointEigenSolver<cdmatrix> eigensolver(hkeff);
	ek=eigensolver.eigenvalues();
	vec=eigensolver.eigenvectors();
	//  now  vec.adjoint() * hk * vec = diag(ek)
	//  hk = vec * diag(ek) * vec.adjoint() 
	//  rho = vec * diag(fermi(ek)) * vec.adjoint() 
	fk.setZero();

	for(int i=0;i<nrpa_;i++) fk(i,i)=cntr::fermi<double>(beta_,ek(i));	
	
	
	rtmp=vec*fk*vec.adjoint();	
	rho_.set_value(tstp,rtmp);
}


template <class LATTICE>
kpoint_green<LATTICE>::kpoint_green(void){
}

template <class LATTICE>
kpoint_green<LATTICE>::kpoint_green(int nt,int ntau,int size,double beta,double h,double kk,LATTICE &latt,double mu,double epsilon,CFUNC &omega, CFUNC &g,double mix,int migdal):
	beta_(beta), h_(h), nt_(nt), ntau_(ntau), kk_(kk),
	nrpa_(size),mu_(mu),epsilon_(epsilon){
	Sigma_=GREEN(nt_,ntau_,nrpa_,-1);
	G_=GREEN(nt_,ntau_,nrpa_,-1);
	G0_=GREEN(nt_,ntau_,nrpa_,-1);

	G0Sigma_=GREEN(nt_,ntau_,nrpa_,-1);
	SigmaG0_=GREEN(nt_,ntau_,nrpa_,-1);
	chi_=GREEN(nt_,ntau_,nrpa_,+1);
	P_=GREEN(nt_,ntau_,nrpa_,+1);
	PV_=GREEN(nt_,ntau_,nrpa_,+1);
	VP_=GREEN(nt_,ntau_,nrpa_,+1);
	W_=GREEN(nt_,ntau_,nrpa_,+1);
	D0_=GREEN(nt_,ntau_,1,+1);
	migdal_=migdal;
	
	if(migdal_==1){
		Pph_=GREEN(nt_,ntau_,1,+1);
		D0Pph_=GREEN(nt_,ntau_,1,+1);
		PphD0_=GREEN(nt_,ntau_,1,+1);
	}
	D_=GREEN(nt_,ntau_,1,+1);
	mix_=mix;
	gph_=g;

	// cdmatrix gtmp(nrpa_,nrpa_);
	// for(int tstp=-1;tstp<=nt;tstp++){
	// 		gph_.get_value(tstp,gtmp);
	// 		std::cout << "Coupling " << gtmp << std::endl;
	// }
	
	// TODO: Here is the assumption of the constant omega0 -  but the routine is limited to this !
	cdmatrix tmp(nrpa_,nrpa_);
	tmp.setZero();
	cntr::green_from_H(G0_,mu_,tmp,beta_,h_);
	omega.get_value(0,tmp);
	// std::cout << "Phonon " <<  std::real(tmp(0,0)) << " " << beta << " " << h << std::endl;
	cntr::green_single_pole_XX(D0_,std::real(tmp(0,0)),beta,h);
	if(migdal_==0){
		cdmatrix tmp1;
		for(int tstp=-1;tstp<=nt;tstp++){
			gph_.get_value(tstp,tmp1);
			// std::cout << "Phonon2 " << tstp << " " << tmp1<< std::endl;
			D0_.right_multiply(tstp,gph_,1.0);
			D0_.left_multiply(tstp,gph_,1.0);
		}
	}
}

template <class LATTICE>
void kpoint_green<LATTICE>::use_omp(bool onoff){
	#if FLEX_CAN_USE_OMP==1
		use_omp_=onoff;
		if(use_omp_){
			std::cout << __PRETTY_FUNCTION__<< ":\n using OMP-parallelizaztion over VIE2 equations ..." << std::endl;
			std::cout << omp_get_max_threads() << " threads" << std::endl;
		}
	#else
		std::cerr << __PRETTY_FUNCTION__<< ":\n OMP not allowed " << std::endl;
		exit(0);
	#endif
}

template <class LATTICE>
void kpoint_density<LATTICE>::set_hk(int tstp,int iter,LATTICE &latt){
  assert(-1<=tstp && tstp<=nt_);
  cdmatrix hktmp(nrpa_,nrpa_);
  latt.hk(hktmp,tstp,kk_,iter);
  hk_.set_value(tstp,hktmp);
}

template <class LATTICE>
void kpoint_density<LATTICE>::set_vertex(int tstp,LATTICE &latt){
  assert(-1<=tstp && tstp<=nt_);
  cdmatrix vktmp(nrpa_,nrpa_);
  latt.V(vktmp,tstp,kk_);
  vertex_.set_value(tstp,vktmp);
}



template <class LATTICE>
void kpoint_green<LATTICE>::get_Density_matrix(int tstp, kpoint_density<LATTICE> &density){
  cdmatrix tmp(nrpa_,nrpa_);
  G_.density_matrix(tstp,tmp);
  density.rho_.set_value(tstp,tmp);
}


/* MEAN field propagation

Solve G=(idt+mu-hk-Hartree-Fock)^{-1} for the density matrix
		
(1 PROPAGATION) 
d/dt rho (t) = d/dt (-i*Gles(t,t)) = - H(t) Gles(t,t) + Gles(t,t)H(t) = -ii * [ H(t), -ii*Gles(t,t) ] = -ii * [H,rho]
with H = heff_ = hk+Hartree+Fock,
rho(t+dt) = U(t+dt,t) rho(t) U(t+dt,t), U(t,t') = T exp -ii * int_{t'}^t ds H(s)
for the simplest, we just use a lowest order implicit propagator:
		 
U(t+dt,t) = exp( -ii dt (H(t+dt + H(t))/0.5) )
(2 EQUILIBRIUM)
*/

template <class LATTICE>
void kpoint_density<LATTICE>::step_dyson(int tstp,int iter,LATTICE &latt,double om0,double s,double amp){
	if(tstp==-1){
		cdmatrix hk,sh,sf,hkeff,rhok,hkeff_eigen,mu(nrpa_,nrpa_);
		mu.setZero();
		mu(0,0)=mu_;
		mu(1,1)=mu_;
		dmatrix hktmp;
		dvector ek;
		cdmatrix vec,fk(nrpa_,nrpa_),rtmp(nrpa_,nrpa_),rtmpOLD(nrpa_,nrpa_),Xtmp(1,1),Ptmp(1,1);
		set_hk(tstp,iter,latt);
		hk_.get_value(tstp,hk);
		SHartree_.get_value(tstp,sh);
		SFock_.get_value(tstp,sf);
		rho_.get_value(tstp,rhok);
		
		hkeff=hk+sh+sf;
		// std::cout <<  "rho " <<tstp << " " << iter << " " << kk_ << " " << rhok <<std::endl;
		// std::cout <<  "Hartree " <<tstp << " " << iter << " " << kk_ << " " << sh <<std::endl;
		// std::cout <<  "Fock " <<tstp << " " << iter << " " << kk_ << " " << sf <<std::endl;
		// std::cout <<  "hk " <<tstp << " " << iter << " " << kk_ << " " << hk <<std::endl;

		hkeff_.set_value(tstp,hkeff);
		// std::cout <<  "hkeff " <<tstp << " " << iter << " " << kk_ << " " << hkeff <<std::endl;
		// std::cout <<  "heff " <<tstp << " " << iter << " " << kk_ << " " << hkeff <<std::endl;
		Eigen::SelfAdjointEigenSolver<cdmatrix> eigensolver(hkeff);
		ek=eigensolver.eigenvalues();
		vec=eigensolver.eigenvectors();
		eigen_vec_.set_value(tstp,vec);
		//  now  vec.adjoint() * hk * vec = diag(ek)
		hktmp = ek.asDiagonal();
		hk=hktmp.cast<std::complex<double> >();
		hkeff_eigen_.set_value(tstp,hk);
		// std::cout <<  "hkeig " <<tstp << " " << iter << " " << kk_ << " " << hk <<std::endl;
		//  rho = vec * diag(fermi(ek)) * vec.adjoint() 
		fk.setZero();
		
		rho_.get_value(-1,rtmpOLD);
		// std::cout << "Density " << den_ << " " << mu_ << std::endl;
		if(iter<20){
	  		for(int i=0;i<nrpa_;i++) fk(i,i)=cntr::fermi<double>(beta_,ek(i)-mu_*den_);
		}else{
			for(int i=0;i<nrpa_;i++) fk(i,i)=cntr::fermi<double>(beta_,ek(i)-mu_*den_+amp*gauss(ek(i),om0,s)-amp*gauss(ek(i),-om0,s));
		}

		// std::cout <<  "rhok " <<tstp << " " << iter << " " << kk_ << " " << vec*fk*vec.adjoint() << " " << fk <<std::endl;

		rtmp=mix_*rtmpOLD+vec*fk*vec.adjoint()*(1.0-mix_);
		rho_diag_.set_value(tstp,fk);
		// std::cout << "fk " <<  fk  << std::endl;
		rho_.set_value(tstp,rtmp);
		//std::cout <<  "rho " <<ek(0) << " " <<ek(1) << " " << iter << " " << rtmp  <<std::endl;
		// rho_loc_.get_value(tstp,rtmp);
  	}
	if(tstp==0){
		cdmatrix Xtmp(1,1),Ptmp(1,1);
		cdmatrix rtmp(nrpa_,nrpa_),rtmpout(nrpa_,nrpa_);
		std::complex<double> I(0.0,1.0); 
		rho_.get_value(-1,rtmp);
		// TODO: Add the initial condition change for general matrix
		// At the moment I don't know how to do this in general case - example for excitonic two band semiconductor
		rtmpout(0,0)=(rtmp(0,0)+I*epsilon_*(rtmp(1,0)-rtmp(0,1))+epsilon_*epsilon_*rtmp(1,1))/(1.0+epsilon_*epsilon_);
		rtmpout(1,0)=(I*epsilon_*(rtmp(0,0)-rtmp(1,1))+epsilon_*epsilon_*rtmp(0,1)+rtmp(1,0))/(1.0+epsilon_*epsilon_);
		rtmpout(0,1)=(I*epsilon_*(rtmp(1,1)-rtmp(0,0))+epsilon_*epsilon_*rtmp(1,0)+rtmp(0,1))/(1.0+epsilon_*epsilon_);
		rtmpout(1,1)=(rtmp(1,1)+I*epsilon_*(rtmp(0,1)-rtmp(1,0))+epsilon_*epsilon_*rtmp(0,0))/(1.0+epsilon_*epsilon_);
		rho_.set_value(0,rtmpout);
		// std::cout << "rho in " << rtmp <<  std::endl;

		// std::cout << "rho out " << rtmpout <<  std::endl;


		rho_diag_.get_value(-1,rtmp);
		rho_diag_.set_value(0,rtmp);
		
		hk_.get_value(-1,rtmp);
		hk_.set_value(0,rtmp);

		SHartree_.get_value(-1,rtmp);
		SHartree_.set_value(0,rtmp);

		SFock_.get_value(-1,rtmp);
		SFock_.set_value(0,rtmp);

		hkeff_.get_value(-1,rtmp);
		hkeff_.set_value(0,rtmp);

		hkeff_eigen_.get_value(-1,rtmp);
		hkeff_eigen_.set_value(0,rtmp);

		eigen_vec_.get_value(-1,rtmp);
		eigen_vec_.set_value(0,rtmp);

	}
	if(tstp>0){
		cdmatrix hk0,sh0,sf0,hk,sh,sf,hkeff0,hkeff1,hkeff;
		dmatrix hktmp;
		dvector ek;
		cdmatrix vec,rtmp(nrpa_,nrpa_),uktmp(nrpa_,nrpa_),uk(nrpa_,nrpa_),rho0(nrpa_,nrpa_),rhotmp(nrpa_,nrpa_);
		cdmatrix Xtmp(1,1),Xtmp_prev(1,1),Ptmp(1,1),Ptmp_prev(1,1);
		///
		// set_hk(tstp-1,iter,latt);
		hk_.get_value(tstp,hk0);
		SHartree_.get_value(tstp-1,sh0);
		SFock_.get_value(tstp-1,sf0);
		hkeff0=hk0+sh0+sf0;
		//
		set_hk(tstp,iter,latt);
		hk_.get_value(tstp,hk);
		SHartree_.get_value(tstp,sh);
		SFock_.get_value(tstp,sf);
		hkeff1=hk+sh+sf;
		// std::cout << "Difference hk  " << hk0-hk  << " " << hk0  << std::endl;
		// std::cout << "Difference sh  " << sh0-sh  <<  " " << sh0 << std::endl;
		// std::cout << "Difference sF  " << sf0-sf  <<  " " << sf0 << std::endl;
		//
		hkeff=0.5*(hkeff0+hkeff1);
		hkeff_.set_value(tstp,hkeff);
		//
		Eigen::SelfAdjointEigenSolver<cdmatrix> eigensolver(hkeff);
		ek=eigensolver.eigenvalues();
		vec=eigensolver.eigenvectors();
		eigen_vec_.set_value(tstp,vec);


		hktmp = ek.asDiagonal();
		hk=hktmp.cast<std::complex<double> >();
		hkeff_eigen_.set_value(tstp,hk);
		//  now  vec.adjoint() * hk * vec = diag(ek)
		//  hk = vec * diag(ek) * vec.adjoint() 
		//  uk= exp(-i*dt*hk ) = vec * diag(exp(-ii*dt*ek) ) * vec.adjoint() 
		uktmp.setZero();
		for(int i=0;i<nrpa_;i++) uktmp(i,i)=CPLX(cos(ek(i)*h_),-sin(ek(i)*h_));
		rho_diag_.set_value(tstp,uktmp);
		uk=vec*uktmp*vec.adjoint();
		// std:: cout << "t=" << tstp << " " << iter << std::endl;
		// std:: cout << "hkP=" << hkeff0 << std::endl;
		// std:: cout << "hk=" << hkeff << std::endl;
		// std:: cout << "hkN=" << hkeff1 << std::endl;
		// std:: cout << "Diff hk=" << hkeff-hkeff0 << std::endl;
		// std:: cout << "uk =" << uk << std::endl;
		// std:: cout << "uk diff=" << uk-uk.conjugate() << std::endl;
		// std:: cout << ",.............," << std::endl;
			
		rho_.get_value(tstp-1,rho0);
		
		rtmp=uk*rho0*uk.adjoint()-gamma_*h_*(rho0-rho_eq_);

		rho_.set_value(tstp,rtmp);
		// rhotmp=vec.adjoint()*rtmp*vec;
		
		// TEST
		// std::cout << "Rtmp: " << tstp << " " << gamma_ << " " <<  -gamma_*h_*(rho0-rho_eq_) << std::endl;
		// std::cout << "eigenvalues" << vec.adjoint()*hkeff*vec  << std::endl;
		// std::cout << "hkeff " << hkeff  << std::endl;
		// std::cout << "hkeig " << ek(0) << " " << ek(1)  << std::endl;
		// std:: cout << ",.............," << std::endl;
	}
}


template <class LATTICE>
double kpoint_density<LATTICE>::step_dyson_with_error(int tstp,int iter,LATTICE &latt,double om0,double s,double amp){
	cdmatrix rtmp,rtmp1,dr;
	rho_.get_value(tstp,rtmp);
	step_dyson(tstp,iter,latt,om0,s,amp);
	rho_.get_value(tstp,rtmp1);
	dr=rtmp1-rtmp;
	// std::cout << "Error1a " << rtmp << std::endl;
	// std::cout << "Error2a " << rtmp1 << std::endl;   
	// std::cout << " ---------------------- "<< std::endl;
	return (dr*dr.adjoint()).trace().real();
}


// START OF RPA CLASS routines
template <class LATTICE>
void kpoint_green<LATTICE>::init_G_mat_nointeraction(LATTICE &latt,kpoint_density<LATTICE>& density,int kt){
	// set G = (idt + mu - hk(U=V=0=A=V01) )^{-1}
	int tstp=-1;
	cdmatrix hktmp(nrpa_,nrpa_);
	latt.hkfree(hktmp,tstp,kk_);
	cdmatrix tmp(nrpa_,nrpa_);
	density.hkeff_.set_value(-1,hktmp);
	Sigma_.set_timestep_zero(tstp);
	cntr::dyson_mat(G_,Sigma_,mu_,density.hkeff_,integration::I<double>(kt),beta_,0,true);
	
	G_.density_matrix(-1,tmp);
}

template <class LATTICE>
void kpoint_green<LATTICE>::step_chi(int tstp,int kt,LATTICE &latt, kpoint_density<LATTICE>& density){
	// solve chi = P + P*V*chi, assumn P is set on all relevant timesteps
	int n;
	int n1=(tstp==-1 || tstp>kt ? tstp : 0);
	int n2=(tstp==-1 || tstp>kt ? tstp : kt);		
	//for(n=n1;n<=n2;n++) chi_.set_timen,P_);
	//return;
	// set PV=-P*V etc.
	for(n=n1;n<=n2;n++){
		// kpoint_density::set_vertex(n,latt); Already set at init 
		PV_.set_timestep(n,P_);
		PV_.right_multiply(n,density.vertex_);
		PV_.smul(n,-1.0);
		VP_.set_timestep(n,P_);
		VP_.left_multiply(n,density.vertex_);
		VP_.smul(n,-1.0);
	}
	// solve [1+PV]*chi=P
	if(tstp==-1){
		cntr::vie2_mat(chi_,PV_,VP_,P_,beta_,integration::I<double>(kt),6);
	}else if (tstp==0){
	    cntr::set_t0_from_mat(chi_);
	}else if(tstp<=kt){
		cntr::vie2_start(chi_,PV_,VP_,P_,integration::I<double>(kt),beta_,h_);
	}else{
		#if FLEX_CAN_USE_OMP==1
		if(this->use_omp_)
		{
			cntr::vie2_timestep_omp(omp_get_max_threads(),tstp,chi_,PV_,VP_,P_,integration::I<double>(kt),beta_,h_);
		}
		else 
		#endif
		{
			cntr::vie2_timestep(tstp,chi_,PV_,VP_,P_,integration::I<double>(kt),this->beta_,this->h_);
		}
	}
}


template <class LATTICE>
void kpoint_green<LATTICE>::step_W(int tstp,int kt,LATTICE &latt, kpoint_density<LATTICE>& density){
	// solve W = V + V*Pi*W, assuming P is set on all relevant timesteps
	int n;
	int n1=(tstp==-1 || tstp>kt ? tstp : 0);
	int n2=(tstp==-1 || tstp>kt ? tstp : kt);		
	//for(n=n1;n<=n2;n++) chi_.set_timestep(n,P_);
	//return;
	// set PV=-P*V etc.
	for(n=n1;n<=n2;n++){
		// kpoint_density::set_vertex(n,latt); Already set at init 
		PV_.set_timestep(n,P_);
		PV_.right_multiply(n,density.vertex_);
		PV_.smul(n,-1.0);
		VP_.set_timestep(n,P_);
		VP_.left_multiply(n,density.vertex_);
		VP_.smul(n,-1.0);
	}
	// solve [1-UP]*W=U
	CFUNC tmpFsin(nt_,nrpa_);
	tmpFsin.set_zero();
	GREEN Q(nt_,ntau_,nrpa_,1);

	cntr::vie2_timestep_sin(tstp,W_,density.vertex_,VP_,PV_,tmpFsin,Q,density.vertex_,beta_,h_,kt);
}

template <class LATTICE>
void kpoint_green<LATTICE>::step_W2b(int tstp,int kt,LATTICE &latt, kpoint_density<LATTICE>& density){
	// Second born approximation for the W
	// solve W = V + V*Pi*V, assuming P is set on all relevant timesteps
	int n;
	int n1=(tstp==-1 || tstp>kt ? tstp : 0);
	int n2=(tstp==-1 || tstp>kt ? tstp : kt);		
	//for(n=n1;n<=n2;n++) chi_.set_timestep(n,P_);
	//return;
	// set PV=-P*V etc.
	for(n=n1;n<=n2;n++){ 
		W_.set_timestep(n,P_);
		W_.right_multiply(n,density.vertex_);
		W_.left_multiply(n,density.vertex_);
	}
}


template <class LATTICE>
void kpoint_green<LATTICE>::step_D(int tstp,int kt,LATTICE &latt, kpoint_density<LATTICE>& density){
	// solve D = D0 + D0*Pph*D, assuming Pph is set on all relevant timesteps

	if(migdal_==0){
		D_=D0_;
	}else if(migdal_==1){
		// D_=D0_;
		int n1=(tstp==-1 || tstp>kt ? tstp : 0);
		int n2=(tstp==-1 || tstp>kt ? tstp : kt);

		//Make convolution 
		for(int n=n1;n<=n2;n++){
			cntr::convolution_timestep(n,D0Pph_,D0_,D0_,Pph_,Pph_,integration::I<double>(kt),beta_,h_);
			cntr::convolution_timestep(n,PphD0_,Pph_,Pph_,D0_,D0_,integration::I<double>(kt),beta_,h_);
		}
		for(int n=n1;n<=n2;n++){ 
			D0Pph_.smul(n,-1.0);
			PphD0_.smul(n,-1.0);
		}
		// D0Pph_.get_mat(0,tmp);
		// std::cout << "D0Pol " << tmp<< std::endl;

		// PphD0_.get_mat(0,tmp);
		// std::cout << "PolD0 " << tmp<< std::endl;

		for(int n=n1;n<=n2;n++){
			cdmatrix tmp;	
			cntr::vie2_timestep(n,D_,D0Pph_,PphD0_,D0_,integration::I<double>(kt),beta_,h_);
		}

		for(int n=n1;n<=n2;n++){
			D_.right_multiply(n,gph_,1.0);
			D_.left_multiply(n,gph_,1.0);
		}
	}
}


template <class LATTICE>
void kpoint_green<LATTICE>::step_dyson(int tstp,int iter,int kt,LATTICE &latt, kpoint_density<LATTICE> &density){
	// solve G=(-idt+mu-hk-Hartree-Fock-Sigma)^{-1}, assuming Sigma is set
	int n;
	int n1=(tstp==-1 || tstp>kt ? tstp : 0);
	int n2=(tstp==-1 || tstp>kt ? tstp : kt);		
	// set hkeff=hk+Hartree+Fock+Symmetry breaking
	for(n=n1;n<=n2;n++){
		cdmatrix tmp,tmp1,tmp2,rtmp(nrpa_,nrpa_),hkeff,hk,hktmp;
		density.set_hk(n,iter,latt);
		density.SHartree_.get_value(n,tmp);
		density.SFock_.get_value(n,tmp1);
		// std::cout <<  "Hartree " <<tstp << " " << iter << " " << kk_ << " " << tmp <<std::endl;
		// std::cout <<  "Fock " <<tstp << " " << iter << " " << kk_ << " " << tmp1 <<std::endl;
		tmp+=tmp1;
		density.hk_.get_value(n,tmp2);
		
		tmp+=tmp2;
		hkeff=tmp;
		density.hkeff_.set_value(n,tmp);
		// std::cout << tstp << " " << n  <<  "hkeff " << tmp <<std::endl;
		if(tstp==-1){
			dvector ek;
			cdmatrix vec;
			Eigen::SelfAdjointEigenSolver<cdmatrix> eigensolver(hkeff);
			ek=eigensolver.eigenvalues();
			vec=eigensolver.eigenvectors();
			//  now  vec.adjoint() * hk * vec = diag(ek)
			hktmp = ek.asDiagonal();
			hk=hktmp.cast<std::complex<double> >();
			density.hkeff_eigen_.set_value(tstp,hk);
		}
		// Fix Hartree fock for kt steps
		if(n>=0 && n<kt){
			density.SHartree_.get_value(-1,tmp);
			density.SFock_.get_value(-1,tmp1);
			tmp+=tmp1;
			density.hk_.get_value(-1,tmp2);
			tmp+=tmp2;
			density.hkeff_.set_value(n,tmp);
			// std::cout << tstp << " " << n  <<  "hkeff fix " << tmp <<std::endl;
		}
	}
	// solve Dyson
	if(tstp==-1 && iter==1){
		GREEN_TSTP tmp(-1,ntau_,nrpa_,G_.sig());
		cntr::green_from_H(tmp,mu_,density.hkeff_,beta_,h_,kt,4,false);
		G_.set_timestep(-1,tmp);
	}else if(tstp==-1){
		cntr::dyson_mat(G_,Sigma_,mu_,density.hkeff_,integration::I<double>(kt),beta_,0,true);
	}else if (tstp==0){
		cdmatrix tmp,tmp1,tmp2,rtmp(nrpa_,nrpa_),hkeff,hk,hktmp;
		// TODO: Add epsilon to the initial condition to check the linear response
	    cntr::set_t0_from_mat(G_);
	    
		


	}else if(tstp<=kt){
		cntr::dyson_start(G_,mu_,density.hkeff_,Sigma_,integration::I<double>(kt),beta_,h_);
	}else{
		cdmatrix rtmp(nrpa_,nrpa_);
		cdmatrix Xtmp_prev(1,1),Ptmp_prev(1,1);
		cntr::dyson_timestep(tstp,G_,mu_,density.hkeff_,Sigma_,integration::I<double>(kt),beta_,h_);
	}
	cdmatrix test;
	G_.density_matrix(tstp,test);
	// std::cout <<  "Density "   << n << " " << iter << " " << kk_ << " " << test << std::endl;
	for(n=n1;n<=n2;n++) get_Density_matrix(n,density);
}

template <class LATTICE>
void kpoint_green<LATTICE>::step_dyson_integral(int tstp,int iter,int kt,LATTICE &latt, kpoint_density<LATTICE> &density){
	// Solve Dyson equation in the integral form (1-G_0*Sigma)*G=G_0
	int n;
	int n1=(tstp==-1 || tstp>kt ? tstp : 0);
	int n2=(tstp==-1 || tstp>kt ? tstp : kt);		
	
	

	// set hkeff=hk+Hartree+Fock+Symmetry breaking
	for(n=n1;n<=n2;n++){
		cdmatrix tmp,tmp1,tmp2,rtmp(nrpa_,nrpa_);
		density.set_hk(n,iter,latt);
		density.SHartree_.get_value(n,tmp);
		density.SFock_.get_value(n,tmp1);
		
		tmp+=tmp1;
		density.hk_.get_value(n,tmp2);
		tmp+=tmp2;
		density.hkeff_.set_value(n,tmp);
		// std::cout << tstp  << " " << n << " "   <<  "hkeff " << tmp <<std::endl;
		if(tstp==-1){
			dvector ek;
			cdmatrix vec,hk,hkeff_eigen,hktmp;
			Eigen::SelfAdjointEigenSolver<cdmatrix> eigensolver(tmp);
			ek=eigensolver.eigenvalues();
			vec=eigensolver.eigenvectors();
			hktmp = ek.asDiagonal();
			hk=hktmp.cast<std::complex<double> >();
			density.hkeff_eigen_.set_value(tstp,hk);
		}
		// Fix Hartree fock for kt steps
		if(n>=0 && n<kt){
			density.SHartree_.get_value(-1,tmp);
			density.SFock_.get_value(-1,tmp1);
			tmp+=tmp1;
			density.hk_.get_value(-1,tmp2);
			tmp+=tmp2;
			density.hkeff_.set_value(-1,tmp);
			// std::cout << tstp << " " << n  <<  "hkeff fix " << tmp <<std::endl;
		}
	}
	
	cdmatrix one=cdmatrix::Identity(nrpa_,nrpa_);
	CFUNC unity(nt_,nrpa_);
	unity.set_constant(one);
	for(n=n1;n<=n2;n++){
		cntr::convolution_timestep(n,G0Sigma_,G0_,G0_,Sigma_,Sigma_,integration::I<double>(kt),beta_,h_);
		cntr::convolution_timestep(n,SigmaG0_,Sigma_,Sigma_,G0_,G0_,integration::I<double>(kt),beta_,h_);
		G0Sigma_.left_multiply(n,unity,-1.0);
		SigmaG0_.left_multiply(n,unity,-1.0);
	}
	for(n=n1;n<=n2;n++){
		GREEN_TSTP tmp(n,ntau_,nrpa_),tmpcc(n,ntau_,nrpa_);
		G0_.get_timestep(n,tmp);
		tmpcc=tmp;
		tmp.right_multiply(density.hkeff_,1.0);
		tmpcc.left_multiply(density.hkeff_,1.0);
		G0Sigma_.incr_timestep(n,tmp,-1.0);
		SigmaG0_.incr_timestep(n,tmpcc,-1.0);
	}
	GREEN_TSTP Gkmix(tstp,ntau_,nrpa_);
	Gkmix.clear();
	
	if(tstp==-1){
		Gkmix.incr(G_,mix_);
	}
	
	cntr::vie2_timestep(tstp,G_,G0Sigma_,SigmaG0_,G0_,integration::I<double>(kt),beta_,h_,kt);
	
	if(tstp==-1){
		Gkmix.incr(G_,1.0 - mix_);
		G_.set_timestep(tstp,Gkmix);
	}

	//cdmatrix test;
	//G_.density_matrix(tstp,test);
	//std::cout <<  "Density "   << n << " " << iter << " " << kk_ << " " << test << std::endl;
	for(n=n1;n<=n2;n++){
	  get_Density_matrix(n,density);  
	}
}

template <class LATTICE>
double kpoint_green<LATTICE>::step_dyson_with_error(int tstp,int iter,int kt,LATTICE &latt,kpoint_density<LATTICE> &density){
	GREEN_TSTP gtmp;
	G_.get_timestep(tstp,gtmp);
	// std::ostringstream nameGp;
	// nameGp << "gkprev_" << tstp <<"_" << iter << "_" <<density.kk_ <<  ".h5";
	// gtmp.write_to_hdf5(nameGp.str().c_str(),"G");
	step_dyson(tstp,iter,kt,latt,density);
	// std::ostringstream nameG;
	// nameG << "gk_" << tstp <<"_" << iter << "_" <<density.kk_ <<  ".h5";
	// G_.write_to_hdf5(nameG.str().c_str(),"G");

	return cntr::distance_norm2(tstp,gtmp,G_);
}

template <class LATTICE>
double kpoint_green<LATTICE>::step_dyson_with_error_integral(int tstp,int iter,int kt,LATTICE &latt,kpoint_density<LATTICE> &density){
	GREEN_TSTP gtmp;
	G_.get_timestep(tstp,gtmp);
	// std::ostringstream nameGp;
	// nameGp << "gkprev_" << tstp <<"_" << iter << "_" <<density.kk_ <<  ".h5";
	// gtmp.write_to_hdf5(nameGp.str().c_str(),"G");
	
	step_dyson_integral(tstp,iter,kt,latt,density);
	//step_dyson(tstp,iter,kt,latt,density);
	
	// std::ostringstream nameG;
	// nameG << "gk_" << tstp <<"_" << iter << "_" <<density.kk_ <<  ".h5";
	// G_.write_to_hdf5(nameG.str().c_str(),"G");

	return cntr::distance_norm2(tstp,gtmp,G_);
}

// template class kpoint_density<lattice_1d_nofield>;
// template class kpoint_green<lattice_1d_nofield>;
