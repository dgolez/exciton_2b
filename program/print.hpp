#pragma once

#include "../step/step.hpp"

template <class LATTICE>
void kpoint_green<LATTICE>::read_from_hdf5(int tstp,hid_t group_id,LATTICE &latt,kpoint_density<LATTICE> &density){
	std::ostringstream groupname;
	GREEN_TSTP tmp;
	groupname.str("");
	groupname << "G/t" << tstp;
	tmp.read_from_hdf5(group_id,groupname.str().c_str());
	G_.set_timestep(tstp,tmp);
		
	groupname.str("");
	groupname << "Sigma/t" << tstp;
	tmp.read_from_hdf5(group_id,groupname.str().c_str());
	Sigma_.set_timestep(tstp,tmp);

	groupname.str("");
	groupname << "chi/t" << tstp;
	tmp.read_from_hdf5(group_id,groupname.str().c_str());
	chi_.set_timestep(tstp,tmp);

	groupname.str("");
	groupname << "P/t" << tstp;
	tmp.read_from_hdf5(group_id,groupname.str().c_str());
	P_.set_timestep(tstp,tmp);
		
	density.SHartree_.read_from_hdf5(group_id,"SHartree");
	density.SFock_.read_from_hdf5(group_id,"SFock");
	density.hkeff_.read_from_hdf5(group_id,"hkeff");
	density.hkeff_.read_from_hdf5(group_id,"hkeff_eigen");
	density.hk_.read_from_hdf5(group_id,"hk");
	density.rho_.read_from_hdf5(group_id,"rho");
}

template <class LATTICE>
void kpoint_green<LATTICE>::read_from_hdf5(int nt1,const char *filename,LATTICE &latt,kpoint_density<LATTICE> &density){
	hid_t file_id = read_hdf5_file(filename);
	this->read_from_hdf5(nt1,file_id,latt,density);
	close_hdf5_file(file_id);
}

template <class LATTICE>
void kpoint_green<LATTICE>::write_to_hdf5(hid_t group_id,kpoint_density<LATTICE> &density){
	hid_t sub_group_id;
	// -- Impurty parameters
	sub_group_id = create_group(group_id, "parm");
	store_int_attribute_to_hid(sub_group_id, "nt", nt_); 
	store_int_attribute_to_hid(sub_group_id, "ntau", ntau_);
	store_double_attribute_to_hid(sub_group_id, "beta", beta_);
	store_double_attribute_to_hid(sub_group_id, "h", h_);
	store_double_attribute_to_hid(sub_group_id, "mu", mu_);
	store_double_attribute_to_hid(sub_group_id, "kk", kk_);
	close_group(sub_group_id); // End parameters
	// -- Green's functions
	G_.write_to_hdf5(group_id,"G");
	Sigma_.write_to_hdf5(group_id, "Sigma");
	chi_.write_to_hdf5(group_id, "chi");
	P_.write_to_hdf5(group_id, "P");
	density.SHartree_.write_to_hdf5(group_id, "SHartree");
	density.SFock_.write_to_hdf5(group_id, "SFock");
	density.hkeff_.write_to_hdf5(group_id, "hkeff");
	density.hk_.write_to_hdf5(group_id, "hk");
	density.rho_.write_to_hdf5(group_id, "rho");
	density.vertex_.write_to_hdf5(group_id, "vertex");

}

template <class LATTICE>
void kpoint_green<LATTICE>::write_to_hdf5(const char *filename,kpoint_density<LATTICE> &density){
	hid_t file_id = open_hdf5_file(filename);
	this->write_to_hdf5(file_id,density);
	close_hdf5_file(file_id);
}


template <class LATTICE>
void kpoint_density<LATTICE>::write_to_hdf5(hid_t group_id){
	hid_t sub_group_id;
	// -- Impurty parameters
	sub_group_id = create_group(group_id, "parm");
	store_int_attribute_to_hid(sub_group_id, "nt", nt_); 
	store_int_attribute_to_hid(sub_group_id, "ntau", ntau_);
	store_double_attribute_to_hid(sub_group_id, "beta", beta_);
	store_double_attribute_to_hid(sub_group_id, "h", h_);
	store_double_attribute_to_hid(sub_group_id, "mu", mu_);
	store_double_attribute_to_hid(sub_group_id, "kk", kk_);
	close_group(sub_group_id); // End parameters
	SHartree_.write_to_hdf5(group_id, "SHartree");
	SFock_.write_to_hdf5(group_id, "SFock");
	hkeff_.write_to_hdf5(group_id, "hkeff");
	hkeff_eigen_.write_to_hdf5(group_id, "hkeff_eigen");
	hk_.write_to_hdf5(group_id, "hk");
	rho_.write_to_hdf5(group_id, "rho");
	vertex_.write_to_hdf5(group_id, "vertex");
}

template <class LATTICE>
void kpoint_density<LATTICE>::write_to_hdf5(const char *filename){
	hid_t file_id = open_hdf5_file(filename);
	this->write_to_hdf5(file_id);
	close_hdf5_file(file_id);
}


template <class LATTICE>
void kpoint_green<LATTICE>::write_to_hdf5_slices(hid_t group_id,int dt,int tid){
	hid_t sub_group_id;
	// -- Impurity parameters
	sub_group_id = create_group(group_id, "parm");
	store_int_attribute_to_hid(sub_group_id, "nt", this->nt_); 
	store_int_attribute_to_hid(sub_group_id, "ntau", this->ntau_);
	store_double_attribute_to_hid(sub_group_id, "beta", this->beta_);
	store_double_attribute_to_hid(sub_group_id, "h", this->h_);
	store_double_attribute_to_hid(sub_group_id, "mu", this->mu_);
	store_double_attribute_to_hid(sub_group_id, "kk", this->kk_);
	close_group(sub_group_id); // End parameters
	// -- Green's functions
	G_.write_to_hdf5_slices(group_id,"G",dt);
	G_.write_to_hdf5_tavtrel(group_id,"Gtavrel",dt);
	Sigma_.write_to_hdf5_slices(group_id, "Sigma",dt);
	chi_.write_to_hdf5_slices(group_id, "chi",dt);
	chi_.write_to_hdf5_tavtrel(group_id, "chitavtrel",dt);
	P_.write_to_hdf5_slices(group_id, "P",dt);
	// density.SHartree_.write_to_hdf5(group_id, "SHartree");
	// density.SFock_.write_to_hdf5(group_id, "SFock");
	// density.hkeff_.write_to_hdf5(group_id, "hkeff");
	// density.hk_.write_to_hdf5(group_id, "hk");
	// density.rho_.write_to_hdf5(group_id, "rho");
	// density.vertex_.write_to_hdf5(group_id, "vertex");
}

template <class LATTICE>
void kpoint_density<LATTICE>::write_to_hdf5_density(hid_t group_id,int dt,int tid){
	hid_t sub_group_id;
	// -- Impurity parameters
	sub_group_id = create_group(group_id, "parm");
	store_int_attribute_to_hid(sub_group_id, "nt", this->nt_); 
	store_int_attribute_to_hid(sub_group_id, "ntau", this->ntau_);
	store_double_attribute_to_hid(sub_group_id, "beta", this->beta_);
	store_double_attribute_to_hid(sub_group_id, "h", this->h_);
	store_double_attribute_to_hid(sub_group_id, "mu", this->mu_);
	store_double_attribute_to_hid(sub_group_id, "kk", this->kk_);
	close_group(sub_group_id); // End parameters
	// -- Green's functions
	this->SHartree_.write_to_hdf5(group_id, "SHartree");
	this->SFock_.write_to_hdf5(group_id, "SFock");
	this->hkeff_.write_to_hdf5(group_id, "hkeff");
	this->hkeff_eigen_.write_to_hdf5(group_id, "hkeff_eigen");
	this->hk_.write_to_hdf5(group_id, "hk");
	this->rho_.write_to_hdf5(group_id, "rho");
	this->vertex_.write_to_hdf5(group_id, "vertex");
}

template <class LATTICE>
void kpoint_green<LATTICE>::write_to_hdf5_slices(const char *filename,int dt,int tid){
	hid_t file_id = open_hdf5_file(filename);
	this->write_to_hdf5_slices(file_id,dt,tid);
	close_hdf5_file(file_id);
}


template <class LATTICE>
void mpi_lattice_step<LATTICE>::read_from_file_hdf5(int tstp,const char *filename_prefix,int kt){
	// TODO !!
}

template <class LATTICE>
void mpi_lattice_deb_step<LATTICE>::read_from_file_hdf5(int tstp,const char *filename_prefix,int kt){
	// TODO !!
}

template <class LATTICE>
void mpi_lattice_step_cdw<LATTICE>::read_from_file_hdf5(int tstp,const char *filename_prefix,int kt){
	// TODO !!
}

template <class LATTICE>
void mpi_lattice_step_optical<LATTICE>::read_from_file_hdf5(int tstp,const char *filename_prefix,int kt){
	// TODO !!
}

template <class LATTICE>
void mpi_lattice_step_hubbard<LATTICE>::read_from_file_hdf5(int tstp,const char *filename_prefix,int kt){
	// TODO !!
}

template <class LATTICE>
void mpi_lattice_step_2b_optical<LATTICE>::read_from_file_hdf5(int tstp,const char *filename_prefix,int kt){
	// TODO !!
}
template <class LATTICE>
void mpi_lattice_step_hubbard_2b<LATTICE>::read_from_file_hdf5(int tstp,const char *filename_prefix,int kt){
	// TODO !!
}
template <class LATTICE>
void mpi_lattice_step_hubbard_RPA<LATTICE>::read_from_file_hdf5(int tstp,const char *filename_prefix,int kt){
	// TODO !!
}

template <class LATTICE>
void mpi_lattice_step_GKBA<LATTICE>::read_from_file_hdf5(int tstp,const char *filename_prefix,int kt){
	// TODO !!
}

template <class LATTICE>
void mpi_lattice_step_RPA<LATTICE>::read_from_file_hdf5(int tstp,const char *filename_prefix,int kt){
	// reading nt+1 timesteps from files generate by print_to_file_hdf5,
	// named [filename_prefix]_local.out, [filename_prefix]_k0.out, ...
	// reads only the grens-functions etc. does not initialize nt ntau etc .
	char fnametmp[1024];
	assert(std::strlen(filename_prefix)<900);
	
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			std::sprintf(fnametmp,"%s_k%d.out",filename_prefix,k);
		    std::cout << " ... reading kk-functions at k= " << k << std::endl;
		    hid_t file_id = read_hdf5_file(fnametmp);
			green_k_[k].read_from_hdf5(tstp,file_id,this->latt_,this->density_k_[k]);
		}
	}
}

template <class LATTICE>
void mpi_lattice_step_RPA_cdw<LATTICE>::read_from_file_hdf5(int tstp,const char *filename_prefix,int kt){
	// reading nt+1 timesteps from files generate by print_to_file_hdf5,
	// named [filename_prefix]_local.out, [filename_prefix]_k0.out, ...
	// reads only the grens-functions etc. does not initialize nt ntau etc .
	char fnametmp[1024];
	assert(std::strlen(filename_prefix)<900);
	
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			std::sprintf(fnametmp,"%s_k%d.out",filename_prefix,k);
		    std::cout << " ... reading kk-functions at k= " << k << std::endl;
		    hid_t file_id = read_hdf5_file(fnametmp);
			green_k_[k].read_from_hdf5(tstp,file_id,this->latt_,this->density_k_[k]);
		}
	}
}


template <class LATTICE>
void mpi_lattice_step_2b_cdw<LATTICE>::read_from_file_hdf5(int tstp,const char *filename_prefix,int kt){
	// reading nt+1 timesteps from files generate by print_to_file_hdf5,
	// named [filename_prefix]_local.out, [filename_prefix]_k0.out, ...
	// reads only the grens-functions etc. does not initialize nt ntau etc .
	char fnametmp[1024];
	assert(std::strlen(filename_prefix)<900);
	
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			std::sprintf(fnametmp,"%s_k%d.out",filename_prefix,k);
		    std::cout << " ... reading kk-functions at k= " << k << std::endl;
		    hid_t file_id = read_hdf5_file(fnametmp);
			green_k_[k].read_from_hdf5(tstp,file_id,this->latt_,this->density_k_[k]);
		}
	}
}


template <class LATTICE>
void mpi_lattice_step_2b<LATTICE>::read_from_file_hdf5(int tstp,const char *filename_prefix,int kt){
	// reading nt+1 timesteps from files generate by print_to_file_hdf5,
	// named [filename_prefix]_local.out, [filename_prefix]_k0.out, ...
	// reads only the grens-functions etc. does not initialize nt ntau etc .
	char fnametmp[1024];
	assert(std::strlen(filename_prefix)<900);
	
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			std::sprintf(fnametmp,"%s_k%d.out",filename_prefix,k);
		    std::cout << " ... reading kk-functions at k= " << k << std::endl;
		    hid_t file_id = read_hdf5_file(fnametmp);
			green_k_[k].read_from_hdf5(tstp,file_id,this->latt_,this->density_k_[k]);
		}
	}
}

template <class LATTICE>
void mpi_lattice_step_2b_deb<LATTICE>::read_from_file_hdf5(int tstp,const char *filename_prefix,int kt){
	// reading nt+1 timesteps from files generate by print_to_file_hdf5,
	// named [filename_prefix]_local.out, [filename_prefix]_k0.out, ...
	// reads only the grens-functions etc. does not initialize nt ntau etc .
	char fnametmp[1024];
	assert(std::strlen(filename_prefix)<900);
	
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			std::sprintf(fnametmp,"%s_k%d.out",filename_prefix,k);
		    std::cout << " ... reading kk-functions at k= " << k << std::endl;
		    hid_t file_id = read_hdf5_file(fnametmp);
			green_k_[k].read_from_hdf5(tstp,file_id,this->latt_,this->density_k_[k]);
		}
	}
}



template <class LATTICE>
void mpi_lattice_step<LATTICE>::print_to_file_hdf5(const char *filename_prefix,int print_k){
	CFUNC ekin=CFUNC(this->nt_,1);
	CFUNC curr=CFUNC(this->nt_,1);
	CFUNC eHF=CFUNC(this->nt_,1);
	CFUNC eRPA=CFUNC(this->nt_,1);
	CFUNC dAk=CFUNC(this->nt_,1);
	CFUNC HAk=CFUNC(this->nt_,1);
	for(int tstp=-1;tstp<=this->nt_;tstp++){
	  cdmatrix tmp(1,1);
	  tmp(0,0)=get_ekin(tstp);
	  ekin.set_value(tstp,tmp);
	  tmp(0,0)=get_eneHF(tstp);
	  eHF.set_value(tstp,tmp);
	  tmp(0,0)=get_dAk(tstp);
	  dAk.set_value(tstp,tmp);

	  // tmp(0,0)=get_curr(tstp);
	  // curr.set_value(tstp,tmp);

	  // int kt1=(this->nt_>=this->kt_ ? this->kt_ : this->nt_);
	  // int n1=(tstp<kt1 ? kt1 : tstp);
	  // double dHAk=0.0;
	  // if(tstp>0){
	  //   for(int n=0;n<=n1;n++){
	  //     dHAk += integration::I<double>(kt1).gregory_weights(tstp,n)*get_dAk(n)*this->h_;
	  //   }
	  // }
	  // tmp(0,0)=dHAk;
	  // HAk.set_value(tstp,tmp);
	}
	if(this->tid_==this->tid_root_){
		char fnametmp[1000];
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		std::cout << "writing hdf5 data to " << fnametmp << std::endl;
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		hid_t group_id = create_group(file_id, "parm");
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		close_group(group_id); // End parameters
		// -- Green's functions
		group_id = create_group(file_id, "obs");
		rho_loc_.write_to_hdf5(group_id,"rho_loc");
		ekin.write_to_hdf5(group_id,"Ekin");
		eHF.write_to_hdf5(group_id,"eHF");
		dAk.write_to_hdf5(group_id,"dAk");
		HAk.write_to_hdf5(group_id,"HAK");
		curr.write_to_hdf5(group_id,"curr");
		X_.write_to_hdf5(group_id,"X");
		Pi_.write_to_hdf5(group_id,"Pi");
		close_group(group_id);
		close_hdf5_file(file_id);
			
		for(int k=0;k<this->nk_;k++){
			char fnametmp[1000];
			std::sprintf(fnametmp,"%s_k%d.out",filename_prefix,k);
			hid_t file_id = open_hdf5_file(std::string(fnametmp));
			density_k_[k].write_to_hdf5(file_id);
			close_hdf5_file(file_id);
		}
	}
}

template <class LATTICE>
void mpi_lattice_deb_step<LATTICE>::print_to_file_hdf5(const char *filename_prefix,int print_k){
	CFUNC ekin=CFUNC(this->nt_,1);
	CFUNC curr=CFUNC(this->nt_,1);
	CFUNC eHF=CFUNC(this->nt_,1);
	CFUNC eRPA=CFUNC(this->nt_,1);
	CFUNC dAk=CFUNC(this->nt_,1);
	CFUNC HAk=CFUNC(this->nt_,1);
	for(int tstp=-1;tstp<=this->nt_;tstp++){
	  cdmatrix tmp(1,1);
	  tmp(0,0)=get_ekin(tstp);
	  ekin.set_value(tstp,tmp);
	  tmp(0,0)=get_eneHF(tstp);
	  eHF.set_value(tstp,tmp);
	  tmp(0,0)=get_dAk(tstp);
	  dAk.set_value(tstp,tmp);

	  // tmp(0,0)=get_curr(tstp);
	  // curr.set_value(tstp,tmp);

	  // int kt1=(this->nt_>=this->kt_ ? this->kt_ : this->nt_);
	  // int n1=(tstp<kt1 ? kt1 : tstp);
	  // double dHAk=0.0;
	  // if(tstp>0){
	  //   for(int n=0;n<=n1;n++){
	  //     dHAk += integration::I<double>(kt1).gregory_weights(tstp,n)*get_dAk(n)*this->h_;
	  //   }
	  // }
	  // tmp(0,0)=dHAk;
	  // HAk.set_value(tstp,tmp);
	}
	if(this->tid_==this->tid_root_){
		char fnametmp[1000];
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		std::cout << "writing hdf5 data to " << fnametmp << std::endl;
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		hid_t group_id = create_group(file_id, "parm");
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		close_group(group_id); // End parameters
		// -- Green's functions
		group_id = create_group(file_id, "obs");
		rho_loc_.write_to_hdf5(group_id,"rho_loc");
		ekin.write_to_hdf5(group_id,"Ekin");
		eHF.write_to_hdf5(group_id,"eHF");
		dAk.write_to_hdf5(group_id,"dAk");
		HAk.write_to_hdf5(group_id,"HAK");
		curr.write_to_hdf5(group_id,"curr");
		X_.write_to_hdf5(group_id,"X");
		Pi_.write_to_hdf5(group_id,"Pi");
		close_group(group_id);
		close_hdf5_file(file_id);
			
		for(int k=0;k<this->nk_;k++){
			char fnametmp[1000];
			std::sprintf(fnametmp,"%s_k%d.out",filename_prefix,k);
			hid_t file_id = open_hdf5_file(std::string(fnametmp));
			density_k_[k].write_to_hdf5(file_id);
			close_hdf5_file(file_id);
		}
	}
}



template <class LATTICE>
void mpi_lattice_step_cdw<LATTICE>::print_to_file_hdf5(const char *filename_prefix,int print_k){
	CFUNC ekin=CFUNC(this->nt_,1);
	CFUNC curr=CFUNC(this->nt_,1);
	CFUNC eHF=CFUNC(this->nt_,1);
	CFUNC eRPA=CFUNC(this->nt_,1);
	CFUNC dAk=CFUNC(this->nt_,1);
	CFUNC HAk=CFUNC(this->nt_,1);
	for(int tstp=-1;tstp<=this->nt_;tstp++){
	  cdmatrix tmp(1,1);
	  tmp(0,0)=get_ekin(tstp);
	  ekin.set_value(tstp,tmp);
	  tmp(0,0)=get_eneHF(tstp);
	  eHF.set_value(tstp,tmp);
	  tmp(0,0)=get_dAk(tstp);
	  dAk.set_value(tstp,tmp);

	  tmp(0,0)=get_curr(tstp);
	  curr.set_value(tstp,tmp);

	  int kt1=(this->nt_>=this->kt_ ? this->kt_ : this->nt_);
	  int n1=(tstp<kt1 ? kt1 : tstp);
	  double dHAk=0.0;
	  if(tstp>0){
	    for(int n=0;n<=n1;n++){
	      dHAk += integration::I<double>(kt1).gregory_weights(tstp,n)*get_dAk(n)*this->h_;
	    }
	  }
	  tmp(0,0)=dHAk;
	  HAk.set_value(tstp,tmp);
	}
	if(this->tid_==this->tid_root_){
		char fnametmp[1000];
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		std::cout << "writing hdf5 data to " << fnametmp << std::endl;
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		hid_t group_id = create_group(file_id, "parm");
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		close_group(group_id); // End parameters
		// -- Green's functions
		group_id = create_group(file_id, "obs");
		rho_loc_.write_to_hdf5(group_id,"rho_loc");
		ekin.write_to_hdf5(group_id,"Ekin");
		eHF.write_to_hdf5(group_id,"eHF");
		dAk.write_to_hdf5(group_id,"dAk");
		HAk.write_to_hdf5(group_id,"HAK");
		curr.write_to_hdf5(group_id,"curr");
		X_.write_to_hdf5(group_id,"X");
		Pi_.write_to_hdf5(group_id,"Pi");
		close_group(group_id);
		close_hdf5_file(file_id);
			
		for(int k=0;k<this->nk_;k++){
			char fnametmp[1000];
			std::sprintf(fnametmp,"%s_k%d.out",filename_prefix,k);
			hid_t file_id = open_hdf5_file(std::string(fnametmp));
			density_k_[k].write_to_hdf5(file_id);
			close_hdf5_file(file_id);
		}
	}
}


template <class LATTICE>
void mpi_lattice_step_optical<LATTICE>::print_to_file_hdf5(const char *filename_prefix,int print_k){
	CFUNC ekin=CFUNC(this->nt_,1);
	CFUNC ekinMAT=CFUNC(this->nt_,2);
	CFUNC curr=CFUNC(this->nt_,1);
	CFUNC currP=CFUNC(this->nt_,1);
	CFUNC currD=CFUNC(this->nt_,1);
	CFUNC dip=CFUNC(this->nt_,1);
	CFUNC eHF=CFUNC(this->nt_,1);
	CFUNC eRPA=CFUNC(this->nt_,1);
	CFUNC dAk=CFUNC(this->nt_,1);
	CFUNC HAk=CFUNC(this->nt_,1);
	for(int tstp=-1;tstp<=this->nt_;tstp++){
	  cdmatrix tmp(1,1);
	  tmp(0,0)=get_ekin(tstp);
	  ekin.set_value(tstp,tmp);

	  cdmatrix tmpkin(2,2);
	  get_ekin(tstp,tmpkin);
	  ekinMAT.set_value(tstp,tmpkin);


	//   tmp(0,0)=get_eneHF(tstp);
	//   eHF.set_value(tstp,tmp);
	//   if(tstp==-1){
	//     tmp(0,0)=std::complex<double>(0.0,0.0);
	//   }else{
	//     tmp(0,0)=get_curr(tstp);
	//   }
	//   curr.set_value(tstp,tmp);

	//   if(tstp==-1){
	//     tmp(0,0)=std::complex<double>(0.0,0.0);
	//   }else{
	//     tmp(0,0)=get_curr_peierls(tstp);
	//   }
	//   currP.set_value(tstp,tmp);

	//   if(tstp==-1){
	//     tmp(0,0)=std::complex<double>(0.0,0.0);
	//   }else{
	//     tmp(0,0)=get_curr_dip(tstp);
	//   }
	//   currD.set_value(tstp,tmp);

	//   tmp(0,0)=get_dip(tstp);
	//   dip.set_value(tstp,tmp);

	}
	double seebeck=0.0;
	seebeck=this->get_seebeck_boltzmann(this->beta_);

	if(this->tid_==this->tid_root_){
		char fnametmp[1000];
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		std::cout << "writing hdf5 data to " << fnametmp << std::endl;
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		hid_t group_id = create_group(file_id, "parm");
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		store_double_attribute_to_hid(group_id, "mu", this->latt_.mu_);
		store_double_attribute_to_hid(group_id, "U", this->latt_.U_[0]);
		close_group(group_id); // End parameters
		// -- Green's functions
		group_id = create_group(file_id, "obs");
		rho_loc_.write_to_hdf5(group_id,"rho_loc");
		order_.write_to_hdf5(group_id,"order");
		ekin.write_to_hdf5(group_id,"Ekin");
		ekinMAT.write_to_hdf5(group_id,"EkinMAT");
		// eHF.write_to_hdf5(group_id,"eHF");
		// dAk.write_to_hdf5(group_id,"dAk");
		// HAk.write_to_hdf5(group_id,"HAK");
		// curr.write_to_hdf5(group_id,"curr");
		// currP.write_to_hdf5(group_id,"currP");
		// currD.write_to_hdf5(group_id,"currD");
		// dip.write_to_hdf5(group_id,"dip");
		X_.write_to_hdf5(group_id,"X");
		Pi_.write_to_hdf5(group_id,"Pi");
		store_double_attribute_to_hid(group_id, "seebeck", seebeck);
		close_group(group_id);
		close_hdf5_file(file_id);
		
		if(print_k){
		  for(int k=0;k<this->nk_;k++){
		    char fnametmp[1000];
		    std::sprintf(fnametmp,"%s_k%d.out",filename_prefix,k);
		    hid_t file_id = open_hdf5_file(std::string(fnametmp));
		    density_k_[k].write_to_hdf5(file_id);
		    close_hdf5_file(file_id);
		  }
		}
	}
}

template <class LATTICE>
void mpi_lattice_step_hubbard<LATTICE>::print_to_file_hdf5(const char *filename_prefix,int print_k){
	CFUNC ekin=CFUNC(this->nt_,1);
	CFUNC curr=CFUNC(this->nt_,1);
	CFUNC currP=CFUNC(this->nt_,1);
	CFUNC currD=CFUNC(this->nt_,1);
	CFUNC dip=CFUNC(this->nt_,1);
	CFUNC eHF=CFUNC(this->nt_,1);
	CFUNC eRPA=CFUNC(this->nt_,1);
	CFUNC dAk=CFUNC(this->nt_,1);
	CFUNC HAk=CFUNC(this->nt_,1);
	for(int tstp=-1;tstp<=this->nt_;tstp++){
	  cdmatrix tmp(1,1);
	  tmp(0,0)=get_ekin(tstp);
	  ekin.set_value(tstp,tmp);
	  tmp(0,0)=get_eneHF(tstp);
	  eHF.set_value(tstp,tmp);
	  if(tstp==-1){
	    tmp(0,0)=std::complex<double>(0.0,0.0);
	  }else{
	    tmp(0,0)=get_curr(tstp);
	  }
	  curr.set_value(tstp,tmp);

	  if(tstp==-1){
	    tmp(0,0)=std::complex<double>(0.0,0.0);
	  }else{
	    tmp(0,0)=get_curr_peierls(tstp);
	  }
	  currP.set_value(tstp,tmp);

	  if(tstp==-1){
	    tmp(0,0)=std::complex<double>(0.0,0.0);
	  }else{
	    tmp(0,0)=get_curr_dip(tstp);
	  }
	  currD.set_value(tstp,tmp);

	  tmp(0,0)=get_dip(tstp);
	  dip.set_value(tstp,tmp);

	}
	double seebeck=0.0;
	seebeck=this->get_seebeck_boltzmann(this->beta_);

	if(this->tid_==this->tid_root_){
		char fnametmp[1000];
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		std::cout << "writing hdf5 data to " << fnametmp << std::endl;
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		hid_t group_id = create_group(file_id, "parm");
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		store_double_attribute_to_hid(group_id, "mu", this->latt_.mu_);
		store_double_attribute_to_hid(group_id, "U", this->latt_.U_[0]);
		close_group(group_id); // End parameters
		// -- Green's functions
		group_id = create_group(file_id, "obs");
		rho_loc_.write_to_hdf5(group_id,"rho_loc");
		ekin.write_to_hdf5(group_id,"Ekin");
		eHF.write_to_hdf5(group_id,"eHF");
		dAk.write_to_hdf5(group_id,"dAk");
		HAk.write_to_hdf5(group_id,"HAK");
		curr.write_to_hdf5(group_id,"curr");
		currP.write_to_hdf5(group_id,"currP");
		currD.write_to_hdf5(group_id,"currD");
		dip.write_to_hdf5(group_id,"dip");
		X_.write_to_hdf5(group_id,"X");
		Pi_.write_to_hdf5(group_id,"Pi");
		store_double_attribute_to_hid(group_id, "seebeck", seebeck);
		close_group(group_id);
		close_hdf5_file(file_id);
		
		if(print_k){
		  for(int k=0;k<this->nk_;k++){
		    char fnametmp[1000];
		    std::sprintf(fnametmp,"%s_k%d.out",filename_prefix,k);
		    hid_t file_id = open_hdf5_file(std::string(fnametmp));
		    density_k_[k].write_to_hdf5(file_id);
		    close_hdf5_file(file_id);
		  }
		}
	}
}





template <class LATTICE>
void mpi_lattice_step_GKBA<LATTICE>::print_to_file_hdf5(const char *filename_prefix,int print_k){

}


template <class LATTICE>
void mpi_lattice_step_RPA<LATTICE>::print_to_file_hdf5(const char *filename_prefix,int print_k){
  assert(std::strlen(filename_prefix)<900);
	// NOTE on MPI:
	// not using parallel Hdf5, so should write in serial, to ensure coherent file access
	// maybe a bit faster, depending on the network structure:write files for each k
	CFUNC ekin=CFUNC(this->nt_,1);
	CFUNC eHF=CFUNC(this->nt_,1);
	//CFUNC eRPA=CFUNC(this->nt_,1);
	CFUNC dAk=CFUNC(this->nt_,1);
	CFUNC HAk=CFUNC(this->nt_,1);
	for(int tstp=-1;tstp<=this->nt_;tstp++){
	  cdmatrix tmp(1,1);
	  tmp(0,0)=base_type::get_ekin(tstp);
	  ekin.set_value(tstp,tmp);
	  tmp(0,0)=base_type::get_eneHF(tstp);
	  eHF.set_value(tstp,tmp);
	  //tmp(0,0)=get_eneRPA(tstp,this->kt_);
	  //eRPA.set_value(tstp,tmp);
	  tmp(0,0)=base_type::get_dAk(tstp);
	  dAk.set_value(tstp,tmp);
	  int kt1=(this->nt_>=this->kt_ ? this->kt_ : this->nt_);
	  int n1=(tstp<kt1 ? kt1 : tstp);
	  double dHAk=0.0;
	  if(tstp>0){
	    for(int n=0;n<=n1;n++){
	      dHAk += integration::I<double>(kt1).gregory_weights(tstp,n)*base_type::get_dAk(n)*this->h_;
	    }
	  }
	  tmp(0,0)=dHAk;
	  HAk.set_value(tstp,tmp);
	}

	if(this->tid_==this->tid_root_){
		char fnametmp[1000];
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		std::cout << "writing hdf5 data to " << fnametmp << std::endl;
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		hid_t group_id = create_group(file_id, "parm");
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		close_group(group_id); // End parameters
		// -- Green's functions
		group_id = create_group(file_id, "Gloc");
		Gloc_.write_to_hdf5(group_id);
		close_group(group_id);
		group_id = create_group(file_id, "Wloc");
		Wloc_.write_to_hdf5(group_id);
		close_group(group_id);
		group_id = create_group(file_id, "obs");
		this->rho_loc_.write_to_hdf5(group_id,"rho_loc");
		ekin.write_to_hdf5(group_id,"Ekin");
		eHF.write_to_hdf5(group_id,"eHF");
		//eRPA.write_to_hdf5(group_id,"eRPA");
		dAk.write_to_hdf5(group_id,"dAk");
		HAk.write_to_hdf5(group_id,"HAK");
		X_.write_to_hdf5(group_id,"X");
		Pi_.write_to_hdf5(group_id,"Pi");
		close_group(group_id);
		close_hdf5_file(file_id);
	}
}

template <class LATTICE>
void mpi_lattice_step_RPA_cdw<LATTICE>::print_to_file_hdf5(const char *filename_prefix,int print_k){
  assert(std::strlen(filename_prefix)<900);
	// NOTE on MPI:
	// not using parallel Hdf5, so should write in serial, to ensure coherent file access
	// maybe a bit faster, depending on the network structure:write files for each k
	CFUNC ekin=CFUNC(this->nt_,1);
	CFUNC eHF=CFUNC(this->nt_,1);
	CFUNC eRPA=CFUNC(this->nt_,1);
	CFUNC dAk=CFUNC(this->nt_,1);
	CFUNC HAk=CFUNC(this->nt_,1);
	for(int tstp=-1;tstp<=this->nt_;tstp++){
	  cdmatrix tmp(1,1);
	  tmp(0,0)=base_type::get_ekin(tstp);
	  ekin.set_value(tstp,tmp);
	  tmp(0,0)=base_type::get_eneHF(tstp);
	  eHF.set_value(tstp,tmp);
	  tmp(0,0)=get_eneRPA(tstp,this->kt_);
	  eRPA.set_value(tstp,tmp);
	  tmp(0,0)=base_type::get_dAk(tstp);
	  dAk.set_value(tstp,tmp);
	  int kt1=(this->nt_>=this->kt_ ? this->kt_ : this->nt_);
	  int n1=(tstp<kt1 ? kt1 : tstp);
	  double dHAk=0.0;
	  if(tstp>0){
	    for(int n=0;n<=n1;n++){
	      dHAk += integration::I<double>(kt1).gregory_weights(tstp,n)*base_type::get_dAk(n)*this->h_;
	    }
	  }
	  tmp(0,0)=dHAk;
	  HAk.set_value(tstp,tmp);
	}

	if(this->tid_==this->tid_root_){
		char fnametmp[1000];
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		std::cout << "writing hdf5 data to " << fnametmp << std::endl;
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		hid_t group_id = create_group(file_id, "parm");
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		close_group(group_id); // End parameters
		// -- Green's functions
		group_id = create_group(file_id, "Gloc");
		Gloc_.write_to_hdf5(group_id);
		close_group(group_id);
		group_id = create_group(file_id, "Wloc");
		Wloc_.write_to_hdf5(group_id);
		close_group(group_id);
		group_id = create_group(file_id, "obs");
		this->rho_loc_.write_to_hdf5(group_id,"rho_loc");
		ekin.write_to_hdf5(group_id,"Ekin");
		eHF.write_to_hdf5(group_id,"eHF");
		eRPA.write_to_hdf5(group_id,"eRPA");
		dAk.write_to_hdf5(group_id,"dAk");
		HAk.write_to_hdf5(group_id,"HAK");
		X_.write_to_hdf5(group_id,"X");
		Pi_.write_to_hdf5(group_id,"Pi");
		close_group(group_id);
		close_hdf5_file(file_id);
	}
}


template <class LATTICE>
void mpi_lattice_step_2b_cdw<LATTICE>::print_to_file_hdf5(const char *filename_prefix,int print_k){
  assert(std::strlen(filename_prefix)<900);
	// NOTE on MPI:
	// not using parallel Hdf5, so should write in serial, to ensure coherent file access
	// maybe a bit faster, depending on the network structure:write files for each k
	CFUNC ekin=CFUNC(this->nt_,1);
	CFUNC eHF=CFUNC(this->nt_,1);
	CFUNC eRPA=CFUNC(this->nt_,1);
	CFUNC dAk=CFUNC(this->nt_,1);
	CFUNC HAk=CFUNC(this->nt_,1);
	for(int tstp=-1;tstp<=this->nt_;tstp++){
	  cdmatrix tmp(1,1);
	  tmp(0,0)=base_type::get_ekin(tstp);
	  ekin.set_value(tstp,tmp);
	  tmp(0,0)=base_type::get_eneHF(tstp);
	  eHF.set_value(tstp,tmp);
	  tmp(0,0)=get_eneRPA(tstp,this->kt_);
	  eRPA.set_value(tstp,tmp);
	  tmp(0,0)=base_type::get_dAk(tstp);
	  dAk.set_value(tstp,tmp);
	  int kt1=(this->nt_>=this->kt_ ? this->kt_ : this->nt_);
	  int n1=(tstp<kt1 ? kt1 : tstp);
	  double dHAk=0.0;
	  if(tstp>0){
	    for(int n=0;n<=n1;n++){
	      dHAk += integration::I<double>(kt1).gregory_weights(tstp,n)*base_type::get_dAk(n)*this->h_;
	    }
	  }
	  tmp(0,0)=dHAk;
	  HAk.set_value(tstp,tmp);
	}

	if(this->tid_==this->tid_root_){
		char fnametmp[1000];
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		std::cout << "writing hdf5 data to " << fnametmp << std::endl;
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		hid_t group_id = create_group(file_id, "parm");
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		close_group(group_id); // End parameters
		// -- Green's functions
		group_id = create_group(file_id, "Gloc");
		Gloc_.write_to_hdf5(group_id);
		close_group(group_id);
		group_id = create_group(file_id, "Wloc");
		Wloc_.write_to_hdf5(group_id);
		close_group(group_id);
		group_id = create_group(file_id, "obs");
		this->rho_loc_.write_to_hdf5(group_id,"rho_loc");
		ekin.write_to_hdf5(group_id,"Ekin");
		eHF.write_to_hdf5(group_id,"eHF");
		eRPA.write_to_hdf5(group_id,"eRPA");
		dAk.write_to_hdf5(group_id,"dAk");
		HAk.write_to_hdf5(group_id,"HAK");
		X_.write_to_hdf5(group_id,"X");
		Pi_.write_to_hdf5(group_id,"Pi");
		close_group(group_id);
		close_hdf5_file(file_id);
	}
}

template <class LATTICE>
void mpi_lattice_step_2b<LATTICE>::print_to_file_hdf5(const char *filename_prefix,int print_k){
    assert(std::strlen(filename_prefix)<900);
	// NOTE on MPI:
	// not using parallel Hdf5, so should write in serial, to ensure coherent file access
	// maybe a bit faster, depending on the network structure:write files for each k
	// std::cout << "base print 1 " << this->tid_ << std::endl;
	CFUNC ekin=CFUNC(this->nt_,1);
	CFUNC eHF=CFUNC(this->nt_,1);
	CFUNC eRPA=CFUNC(this->nt_,1);
	CFUNC dAk=CFUNC(this->nt_,1);
	CFUNC HAk=CFUNC(this->nt_,1);
	// std::cout << "base print 2 " << this->tid_ << std::endl;
	for(int tstp=-1;tstp<=this->nt_;tstp++){
	  // std::cout << "base print 3 " << this->tid_ << " " << tstp << std::endl;
	  cdmatrix tmp(1,1);
	  tmp(0,0)=base_type::get_ekin(tstp);
	  ekin.set_value(tstp,tmp);
	  // std::cout << "base print 4 " << this->tid_ << " " << tstp << std::endl;
	  tmp(0,0)=base_type::get_eneHF(tstp);
	  eHF.set_value(tstp,tmp);
	  // std::cout << "base print 5 " << this->tid_ << " " << tstp << std::endl;
	  tmp(0,0)=get_eneRPA(tstp,this->kt_);
	  eRPA.set_value(tstp,tmp);
	  // std::cout << "base print 6 " << this->tid_ << " " << tstp << " " << tmp << std::endl;
	  tmp(0,0)=base_type::get_dAk(tstp);
	  dAk.set_value(tstp,tmp);
	  // std::cout << "base print 7 " << this->tid_ << " " << tstp << std::endl;
	  // int kt1=(this->nt_>=this->kt_ ? this->kt_ : this->nt_);
	  // int n1=(tstp<kt1 ? kt1 : tstp);
	  // double dHAk=0.0;
	  // if(tstp>0){
	  //   for(int n=0;n<=n1;n++){
	  //     dHAk += integration::I<double>(kt1).gregory_weights(tstp,n)*base_type::get_dAk(n)*this->h_;
	  //   }
	  // }
	  // tmp(0,0)=dHAk;
	  // HAk.set_value(tstp,tmp);
	}
	
	if(this->tid_==this->tid_root_){
		// std::cout << "base print 8 " << this->tid_ << std::endl;
		char fnametmp[1000];
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		std::cout << "writing hdf5 data to " << fnametmp << std::endl;
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		// std::cout << "base print 9 " << this->tid_ << std::endl;
		hid_t group_id = create_group(file_id, "parm");
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		close_group(group_id); // End parameters
		// std::cout << "base print 10 " << this->tid_ << std::endl;
		// -- Green's functions
		group_id = create_group(file_id, "Gloc");
		Gloc_.write_to_hdf5(group_id);
		close_group(group_id);
		// std::cout << "base print 11 " << this->tid_ << std::endl;
		group_id = create_group(file_id, "Wloc");
		Wloc_.write_to_hdf5(group_id);
		close_group(group_id);
		// std::cout << "base print 12 " << this->tid_ << std::endl;
		group_id = create_group(file_id, "Dloc");
		Dloc_.write_to_hdf5(group_id);
		close_group(group_id);
		// std::cout << "base print 13 " << this->tid_ << std::endl;
		group_id = create_group(file_id, "obs");
		this->rho_loc_.write_to_hdf5(group_id,"rho_loc");
		this->order_.write_to_hdf5(group_id,"order");
		this->rho_sym_.write_to_hdf5(group_id,"order_cos");
		ekin.write_to_hdf5(group_id,"Ekin");
		eHF.write_to_hdf5(group_id,"eHF");
		eRPA.write_to_hdf5(group_id,"eRPA");
		dAk.write_to_hdf5(group_id,"dAk");
		HAk.write_to_hdf5(group_id,"HAK");
		X_.write_to_hdf5(group_id,"X");
		Pi_.write_to_hdf5(group_id,"Pi");
		// std::cout << "base print 14 " << this->tid_ << std::endl;
		close_group(group_id);
		close_hdf5_file(file_id);
	}
	// std::cout << "base print 15 " << this->tid_ << std::endl;
}



template <class LATTICE>
void mpi_lattice_step_2b_deb<LATTICE>::print_to_file_hdf5(const char *filename_prefix,int print_k){
    assert(std::strlen(filename_prefix)<900);
	// NOTE on MPI:
	// not using parallel Hdf5, so should write in serial, to ensure coherent file access
	// maybe a bit faster, depending on the network structure:write files for each k
	// std::cout << "base print 1 " << this->tid_ << std::endl;
	CFUNC ekin=CFUNC(this->nt_,1);
	CFUNC eHF=CFUNC(this->nt_,1);
	CFUNC eRPA=CFUNC(this->nt_,1);
	CFUNC dAk=CFUNC(this->nt_,1);
	CFUNC HAk=CFUNC(this->nt_,1);
	// std::cout << "base print 2 " << this->tid_ << std::endl;
	for(int tstp=-1;tstp<=this->nt_;tstp++){
	  // std::cout << "base print 3 " << this->tid_ << " " << tstp << std::endl;
	  cdmatrix tmp(1,1);
	  tmp(0,0)=base_type::get_ekin(tstp);
	  ekin.set_value(tstp,tmp);
	  // std::cout << "base print 4 " << this->tid_ << " " << tstp << std::endl;
	  tmp(0,0)=base_type::get_eneHF(tstp);
	  eHF.set_value(tstp,tmp);
	  // std::cout << "base print 5 " << this->tid_ << " " << tstp << std::endl;
	  tmp(0,0)=get_eneRPA(tstp,this->kt_);
	  eRPA.set_value(tstp,tmp);
	  // std::cout << "base print 6 " << this->tid_ << " " << tstp << " " << tmp << std::endl;
	  tmp(0,0)=base_type::get_dAk(tstp);
	  dAk.set_value(tstp,tmp);
	  // std::cout << "base print 7 " << this->tid_ << " " << tstp << std::endl;
	  // int kt1=(this->nt_>=this->kt_ ? this->kt_ : this->nt_);
	  // int n1=(tstp<kt1 ? kt1 : tstp);
	  // double dHAk=0.0;
	  // if(tstp>0){
	  //   for(int n=0;n<=n1;n++){
	  //     dHAk += integration::I<double>(kt1).gregory_weights(tstp,n)*base_type::get_dAk(n)*this->h_;
	  //   }
	  // }
	  // tmp(0,0)=dHAk;
	  // HAk.set_value(tstp,tmp);
	}
	
	if(this->tid_==this->tid_root_){
		// std::cout << "base print 8 " << this->tid_ << std::endl;
		char fnametmp[1000];
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		std::cout << "writing hdf5 data to " << fnametmp << std::endl;
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		// std::cout << "base print 9 " << this->tid_ << std::endl;
		hid_t group_id = create_group(file_id, "parm");
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		close_group(group_id); // End parameters
		// std::cout << "base print 10 " << this->tid_ << std::endl;
		// -- Green's functions
		group_id = create_group(file_id, "Gloc");
		Gloc_.write_to_hdf5(group_id);
		close_group(group_id);
		// std::cout << "base print 11 " << this->tid_ << std::endl;
		group_id = create_group(file_id, "Wloc");
		Wloc_.write_to_hdf5(group_id);
		close_group(group_id);
		// std::cout << "base print 12 " << this->tid_ << std::endl;
		group_id = create_group(file_id, "Dloc");
		Dloc_.write_to_hdf5(group_id);
		close_group(group_id);
		// std::cout << "base print 13 " << this->tid_ << std::endl;
		group_id = create_group(file_id, "obs");
		this->rho_loc_.write_to_hdf5(group_id,"rho_loc");
		ekin.write_to_hdf5(group_id,"Ekin");
		eHF.write_to_hdf5(group_id,"eHF");
		eRPA.write_to_hdf5(group_id,"eRPA");
		dAk.write_to_hdf5(group_id,"dAk");
		HAk.write_to_hdf5(group_id,"HAK");
		X_.write_to_hdf5(group_id,"X");
		Pi_.write_to_hdf5(group_id,"Pi");
		// std::cout << "base print 14 " << this->tid_ << std::endl;
		close_group(group_id);
		close_hdf5_file(file_id);
	}
	// std::cout << "base print 15 " << this->tid_ << std::endl;
}

template <class LATTICE>
void mpi_lattice_step_2b_optical<LATTICE>::print_to_file_hdf5(const char *filename_prefix,int print_k){
  assert(std::strlen(filename_prefix)<900);
	// NOTE on MPI:
	// not using parallel Hdf5, so should write in serial, to ensure coherent file access
	// maybe a bit faster, depending on the network structure:write files for each k
	std::cout << "Print 1 " << this->tid_ << std::endl;
	CFUNC ekin=CFUNC(this->nt_,1);
	CFUNC ekinMAT=CFUNC(this->nt_,2);
	CFUNC eHF=CFUNC(this->nt_,1);
	CFUNC eRPA=CFUNC(this->nt_,1);
	CFUNC dAk=CFUNC(this->nt_,1);
	CFUNC HAk=CFUNC(this->nt_,1);
	CFUNC currP=CFUNC(this->nt_,1);
	CFUNC currD=CFUNC(this->nt_,1);
	// std::cout << "Print 2 " << this->tid_ << std::endl;
	for(int tstp=-1;tstp<=this->nt_;tstp++){
	//   std::cout << "Print 2a " << this->tid_ << std::endl;
	  cdmatrix tmp(1,1);
	  tmp(0,0)=base_type::get_ekin(tstp);
	  ekin.set_value(tstp,tmp);
	  
	//   std::cout << "Print 2b " << this->tid_ << std::endl;
	//   cdmatrix tmpkin(2,2);
	//   base_type::get_ekin(tstp,tmpkin);
	//   ekinMAT.set_value(tstp,tmpkin);

	//   std::cout << "Print 2c " << this->tid_ << std::endl;
	//   tmp(0,0)=base_type::get_eneHF(tstp);
	//   eHF.set_value(tstp,tmp);

	//   std::cout << "Print 2d " << this->tid_ << std::endl;
	//   tmp(0,0)=get_eneRPA(tstp,this->kt_);
	//   eRPA.set_value(tstp,tmp);

	  // tmp(0,0)=base_type::get_dAk(tstp);

	  // dAk.set_value(tstp,tmp);
	  // std::cout << "Print 3 " << this->tid_ << std::endl;

	//   int kt1=(this->nt_>=this->kt_ ? this->kt_ : this->nt_);
	//   int n1=(tstp<kt1 ? kt1 : tstp);
	  // double dHAk=0.0;
	  // if(tstp>0){
	  //   for(int n=0;n<=n1;n++){
	  //     dHAk += integration::I<double>(kt1).gregory_weights(tstp,n)*base_type::get_dAk(n)*this->h_;
	  //   }
	  // }
	  // tmp(0,0)=dHAk;
	  // HAk.set_value(tstp,tmp);
	}
	std::cout << "Print 4 " << this->tid_ << std::endl;
	if(this->tid_==this->tid_root_){
		// std::cout << "Print 5 " << this->tid_ << std::endl;
		char fnametmp[1000];
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		std::cout << "writing hdf5 data to " << fnametmp << std::endl;
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		hid_t group_id = create_group(file_id, "parm");
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		close_group(group_id); // End parameters
		// -- Green's functions
		group_id = create_group(file_id, "Gloc");
		Gloc_.write_to_hdf5(group_id);
		close_group(group_id);
		group_id = create_group(file_id, "Wloc");
		Wloc_.write_to_hdf5(group_id);
		group_id = create_group(file_id, "Dloc");
		Dloc_.write_to_hdf5(group_id);
		group_id = create_group(file_id, "D0loc");
		D0loc_.write_to_hdf5(group_id);
		close_group(group_id);
		group_id = create_group(file_id, "obs");
		this->rho_loc_.write_to_hdf5(group_id,"rho_loc");
		this->order_.write_to_hdf5(group_id,"order");
		ekin.write_to_hdf5(group_id,"Ekin");
		// ekinMAT.write_to_hdf5(group_id,"EkinMAT");
		// eHF.write_to_hdf5(group_id,"eHF");
		// eRPA.write_to_hdf5(group_id,"eRPA");
		// dAk.write_to_hdf5(group_id,"dAk");
		// HAk.write_to_hdf5(group_id,"HAK");
		// curr.write_to_hdf5(group_id,"curr");
		// currP.write_to_hdf5(group_id,"currP");
		// currD.write_to_hdf5(group_id,"currD");
		// dip.write_to_hdf5(group_id,"dip");
		X_.write_to_hdf5(group_id,"X");
		Pi_.write_to_hdf5(group_id,"Pi");
		// Pi_.write_to_hdf5(group_id,"Pi_phonon");
		// std::cout << "Print 6 " << std::endl;
		close_group(group_id);
		close_hdf5_file(file_id);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

template <class LATTICE>
void mpi_lattice_step_hubbard_2b<LATTICE>::print_to_file_hdf5(const char *filename_prefix,int print_k){
  assert(std::strlen(filename_prefix)<900);
	// NOTE on MPI:
	// not using parallel Hdf5, so should write in serial, to ensure coherent file access
	// maybe a bit faster, depending on the network structure:write files for each k
	// std::cout << "Print 1 " << this->tid_ << std::endl;
	CFUNC ekin=CFUNC(this->nt_,1);
	CFUNC eHF=CFUNC(this->nt_,1);
	CFUNC eRPA=CFUNC(this->nt_,1);
	CFUNC dAk=CFUNC(this->nt_,1);
	CFUNC HAk=CFUNC(this->nt_,1);
	CFUNC currP=CFUNC(this->nt_,1);
	CFUNC currD=CFUNC(this->nt_,1);
	// std::cout << "Print 2 " << this->tid_ << std::endl;
	for(int tstp=-1;tstp<=this->nt_;tstp++){

	  cdmatrix tmp(1,1);
	  tmp(0,0)=base_type::get_ekin(tstp);
	  ekin.set_value(tstp,tmp);

	  tmp(0,0)=base_type::get_eneHF(tstp);
	  eHF.set_value(tstp,tmp);

	  tmp(0,0)=get_eneRPA(tstp,this->kt_);
	  eRPA.set_value(tstp,tmp);

	  // tmp(0,0)=base_type::get_dAk(tstp);

	  // dAk.set_value(tstp,tmp);
	  // std::cout << "Print 3 " << this->tid_ << std::endl;
	  int kt1=(this->nt_>=this->kt_ ? this->kt_ : this->nt_);
	  int n1=(tstp<kt1 ? kt1 : tstp);
	  // double dHAk=0.0;
	  // if(tstp>0){
	  //   for(int n=0;n<=n1;n++){
	  //     dHAk += integration::I<double>(kt1).gregory_weights(tstp,n)*base_type::get_dAk(n)*this->h_;
	  //   }
	  // }
	  // tmp(0,0)=dHAk;
	  // HAk.set_value(tstp,tmp);
	}
	// std::cout << "Print 4 " << this->tid_ << std::endl;
	if(this->tid_==this->tid_root_){
		// std::cout << "Print 5 " << this->tid_ << std::endl;
		char fnametmp[1000];
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		std::cout << "writing hdf5 data to " << fnametmp << std::endl;
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		hid_t group_id = create_group(file_id, "parm");
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		close_group(group_id); // End parameters
		// -- Green's functions
		group_id = create_group(file_id, "Gloc");
		Gloc_.write_to_hdf5(group_id);
		close_group(group_id);
		group_id = create_group(file_id, "Wloc");
		Wloc_.write_to_hdf5(group_id);
		group_id = create_group(file_id, "Dloc");
		Dloc_.write_to_hdf5(group_id);
		group_id = create_group(file_id, "D0loc");
		D0loc_.write_to_hdf5(group_id);
		close_group(group_id);
		group_id = create_group(file_id, "obs");
		this->rho_loc_.write_to_hdf5(group_id,"rho_loc");
		ekin.write_to_hdf5(group_id,"Ekin");
		eHF.write_to_hdf5(group_id,"eHF");
		eRPA.write_to_hdf5(group_id,"eRPA");
		// dAk.write_to_hdf5(group_id,"dAk");
		// HAk.write_to_hdf5(group_id,"HAK");
		// curr.write_to_hdf5(group_id,"curr");
		// currP.write_to_hdf5(group_id,"currP");
		// currD.write_to_hdf5(group_id,"currD");
		// dip.write_to_hdf5(group_id,"dip");
		X_.write_to_hdf5(group_id,"X");
		Pi_.write_to_hdf5(group_id,"Pi");
		// Pi_.write_to_hdf5(group_id,"Pi_phonon");
		// std::cout << "Print 6 " << std::endl;
		close_group(group_id);
		close_hdf5_file(file_id);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

template <class LATTICE>
void mpi_lattice_step_hubbard_RPA<LATTICE>::print_to_file_hdf5(const char *filename_prefix,int print_k){
  assert(std::strlen(filename_prefix)<900);
	// NOTE on MPI:
	// not using parallel Hdf5, so should write in serial, to ensure coherent file access
	// maybe a bit faster, depending on the network structure:write files for each k
	// std::cout << "Print 1 " << this->tid_ << std::endl;
	CFUNC ekin=CFUNC(this->nt_,1);
	CFUNC eHF=CFUNC(this->nt_,1);
	CFUNC eRPA=CFUNC(this->nt_,1);
	CFUNC dAk=CFUNC(this->nt_,1);
	CFUNC HAk=CFUNC(this->nt_,1);
	CFUNC currP=CFUNC(this->nt_,1);
	CFUNC currD=CFUNC(this->nt_,1);
	// std::cout << "Print 2 " << this->tid_ << std::endl;
	for(int tstp=-1;tstp<=this->nt_;tstp++){

	  cdmatrix tmp(1,1);
	  tmp(0,0)=base_type::get_ekin(tstp);
	  ekin.set_value(tstp,tmp);

	  tmp(0,0)=base_type::get_eneHF(tstp);
	  eHF.set_value(tstp,tmp);

	  tmp(0,0)=get_eneRPA(tstp,this->kt_);
	  eRPA.set_value(tstp,tmp);

	  // tmp(0,0)=base_type::get_dAk(tstp);

	  // dAk.set_value(tstp,tmp);
	  // std::cout << "Print 3 " << this->tid_ << std::endl;
	  int kt1=(this->nt_>=this->kt_ ? this->kt_ : this->nt_);
	  int n1=(tstp<kt1 ? kt1 : tstp);
	  // double dHAk=0.0;
	  // if(tstp>0){
	  //   for(int n=0;n<=n1;n++){
	  //     dHAk += integration::I<double>(kt1).gregory_weights(tstp,n)*base_type::get_dAk(n)*this->h_;
	  //   }
	  // }
	  // tmp(0,0)=dHAk;
	  // HAk.set_value(tstp,tmp);
	}
	// std::cout << "Print 4 " << this->tid_ << std::endl;
	if(this->tid_==this->tid_root_){
		// std::cout << "Print 5 " << this->tid_ << std::endl;
		char fnametmp[1000];
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		std::cout << "writing hdf5 data to " << fnametmp << std::endl;
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		hid_t group_id = create_group(file_id, "parm");
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		close_group(group_id); // End parameters
		// -- Green's functions
		group_id = create_group(file_id, "Gloc");
		Gloc_.write_to_hdf5(group_id);
		close_group(group_id);
		group_id = create_group(file_id, "Wloc");
		Wloc_.write_to_hdf5(group_id);
		group_id = create_group(file_id, "Dloc");
		Dloc_.write_to_hdf5(group_id);
		group_id = create_group(file_id, "D0loc");
		D0loc_.write_to_hdf5(group_id);
		close_group(group_id);
		group_id = create_group(file_id, "obs");
		this->rho_loc_.write_to_hdf5(group_id,"rho_loc");
		ekin.write_to_hdf5(group_id,"Ekin");
		eHF.write_to_hdf5(group_id,"eHF");
		eRPA.write_to_hdf5(group_id,"eRPA");
		// dAk.write_to_hdf5(group_id,"dAk");
		// HAk.write_to_hdf5(group_id,"HAK");
		// curr.write_to_hdf5(group_id,"curr");
		// currP.write_to_hdf5(group_id,"currP");
		// currD.write_to_hdf5(group_id,"currD");
		// dip.write_to_hdf5(group_id,"dip");
		X_.write_to_hdf5(group_id,"X");
		Pi_.write_to_hdf5(group_id,"Pi");
		// Pi_.write_to_hdf5(group_id,"Pi_phonon");
		// std::cout << "Print 6 " << std::endl;
		close_group(group_id);
		close_hdf5_file(file_id);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}



template <class LATTICE>
void mpi_lattice_step<LATTICE>::print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k){
}

template <class LATTICE>
void mpi_lattice_deb_step<LATTICE>::print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k){
}

template <class LATTICE>
void mpi_lattice_step_cdw<LATTICE>::print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k){
}

template <class LATTICE>
void mpi_lattice_step_optical<LATTICE>::print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k){
}

template <class LATTICE>
void mpi_lattice_step_hubbard<LATTICE>::print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k){
}

// 


template <class LATTICE>
void mpi_lattice_step_GKBA<LATTICE>::print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k){
}

template <class LATTICE>
void mpi_lattice_step_RPA<LATTICE>::print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k){
  assert(std::strlen(filename_prefix)<900);
	// NOTE on MPI:
	// not using parallel hdf5, so should write in serial, to ensure coherent file access
	// maybe a bit faster, depending on the network structure:write files for each k
	// std::cout << "slice 1" << std::endl;
	if(this->tid_==this->tid_root_){
		char fnametmp[1000];
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		hid_t group_id = create_group(file_id, "parm");
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		close_group(group_id); // End parameters
		// -- Green's functions
		group_id = create_group(file_id, "Gloc");
		Gloc_.write_to_hdf5_slices(group_id,dt);
		close_group(group_id);
		group_id = create_group(file_id, "Gloc_full");
		Gloc_.write_to_hdf5(group_id);
		close_group(group_id);
		group_id = create_group(file_id, "Gtavrel");
		Gloc_.write_to_hdf5_tavtrel(group_id,dt);
		close_group(group_id);
		close_hdf5_file(file_id);
	}
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			char fnametmp[1000];
			std::sprintf(fnametmp,"%s_denk%d.out",filename_prefix,k);
			std::cout << "writing hdf5 data to " << fnametmp << std::endl;
			hid_t file_id = open_hdf5_file(std::string(fnametmp));
			//hid_t group_id = create_group(file_id,std::string(suffix));
			this->density_k_[k].write_to_hdf5_density(file_id,dt,this->tid_);
			close_hdf5_file(file_id);
		}
	}

	// std::cout << "slice 2 " << this->tid_ << std::endl;
	if(print_k){
		for(int k=0;k<this->nk_;k++){
			if(this->tid_map_[k]==this->tid_){
				// std::cout << "slice 3" << this->tid_ <<  std::endl;
				char fnametmp[1000];
				std::sprintf(fnametmp,"%s_k%d.out",filename_prefix,k);
				std::cout << "writing hdf5 data to " << fnametmp << std::endl;
				hid_t file_id = open_hdf5_file(std::string(fnametmp));
				//hid_t group_id = create_group(file_id,std::string(suffix));
				green_k_[k].write_to_hdf5_slices(file_id,dt,this->tid_);
				close_hdf5_file(file_id);
			}
		}	
	}
	
	if(this->tid_==this->tid_root_){
		// merge???
	}
}

template <class LATTICE>
void mpi_lattice_step_RPA_cdw<LATTICE>::print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k){
  assert(std::strlen(filename_prefix)<900);
	// NOTE on MPI:
	// not using parallel hdf5, so should write in serial, to ensure coherent file access
	// maybe a bit faster, depending on the network structure:write files for each k
	if(this->tid_==this->tid_root_){
		char fnametmp[1000];
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		hid_t group_id = create_group(file_id, "parm");
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		close_group(group_id); // End parameters
		// -- Green's functions
		group_id = create_group(file_id, "Gloc");
		Gloc_.write_to_hdf5_slices(group_id,dt);
		close_group(group_id);
		group_id = create_group(file_id, "Gloc_full");
		Gloc_.write_to_hdf5(group_id);
		close_group(group_id);
		group_id = create_group(file_id, "Gtavrel");
		Gloc_.write_to_hdf5_tavtrel(group_id,dt);
		close_group(group_id);
		close_hdf5_file(file_id);
	}
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			char fnametmp[1000];
			std::sprintf(fnametmp,"%s_denk%d.out",filename_prefix,k);
			std::cout << "writing hdf5 data to " << fnametmp << std::endl;
			hid_t file_id = open_hdf5_file(std::string(fnametmp));
			//hid_t group_id = create_group(file_id,std::string(suffix));
			this->density_k_[k].write_to_hdf5_density(file_id,dt,this->tid_);
			close_hdf5_file(file_id);
		}
	}

	if(print_k){
		for(int k=0;k<this->nk_;k++){
			if(this->tid_map_[k]==this->tid_){
				char fnametmp[1000];
				std::sprintf(fnametmp,"%s_k%d.out",filename_prefix,k);
				std::cout << "writing hdf5 data to " << fnametmp << std::endl;
				hid_t file_id = open_hdf5_file(std::string(fnametmp));
				//hid_t group_id = create_group(file_id,std::string(suffix));
				green_k_[k].write_to_hdf5_slices(file_id,dt,this->tid_);
				close_hdf5_file(file_id);
			}
		}	
	}
	
	if(this->tid_==this->tid_root_){
		// merge???
	}
}



template <class LATTICE>
void mpi_lattice_step_2b_cdw<LATTICE>::print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k){
  assert(std::strlen(filename_prefix)<900);
	// NOTE on MPI:
	// not using parallel hdf5, so should write in serial, to ensure coherent file access
	// maybe a bit faster, depending on the network structure:write files for each k
	if(this->tid_==this->tid_root_){
		char fnametmp[1000];
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		hid_t group_id = create_group(file_id, "parm");
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		close_group(group_id); // End parameters
		// -- Green's functions
		group_id = create_group(file_id, "Gloc");
		Gloc_.write_to_hdf5_slices(group_id,dt);
		close_group(group_id);
		group_id = create_group(file_id, "Gloc_full");
		Gloc_.write_to_hdf5(group_id);
		close_group(group_id);
		group_id = create_group(file_id, "Gtavrel");
		Gloc_.write_to_hdf5_tavtrel(group_id,dt);
		close_group(group_id);
		close_hdf5_file(file_id);
	}
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			char fnametmp[1000];
			std::sprintf(fnametmp,"%s_denk%d.out",filename_prefix,k);
			std::cout << "writing hdf5 data to " << fnametmp << std::endl;
			hid_t file_id = open_hdf5_file(std::string(fnametmp));
			//hid_t group_id = create_group(file_id,std::string(suffix));
			this->density_k_[k].write_to_hdf5_density(file_id,dt,this->tid_);
			close_hdf5_file(file_id);
		}
	}

	if(print_k){
		for(int k=0;k<this->nk_;k++){
			if(this->tid_map_[k]==this->tid_){
				char fnametmp[1000];
				std::sprintf(fnametmp,"%s_k%d.out",filename_prefix,k);
				std::cout << "writing hdf5 data to " << fnametmp << std::endl;
				hid_t file_id = open_hdf5_file(std::string(fnametmp));
				//hid_t group_id = create_group(file_id,std::string(suffix));
				green_k_[k].write_to_hdf5_slices(file_id,dt,this->tid_);
				close_hdf5_file(file_id);
			}
		}	
	}
	
	if(this->tid_==this->tid_root_){
		// merge???
	}
}


template <class LATTICE>
void mpi_lattice_step_2b<LATTICE>::print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k){
  assert(std::strlen(filename_prefix)<900);
	// NOTE on MPI:
	// not using parallel hdf5, so should write in serial, to ensure coherent file access
	// maybe a bit faster, depending on the network structure:write files for each k
	// std::cout << "slice 1" << std::endl;
	if(this->tid_==this->tid_root_){
		char fnametmp[1000];
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		hid_t group_id = create_group(file_id, "parm");
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		close_group(group_id); // End parameters
		// -- Green's functions
		group_id = create_group(file_id, "Gloc");
		Gloc_.write_to_hdf5_slices(group_id,dt);
		close_group(group_id);
		group_id = create_group(file_id, "Gloc_full");
		Gloc_.write_to_hdf5(group_id);
		close_group(group_id);
		group_id = create_group(file_id, "Gtavrel");
		Gloc_.write_to_hdf5_tavtrel(group_id,dt);
		close_group(group_id);
		close_hdf5_file(file_id);
	}
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			char fnametmp[1000];
			std::sprintf(fnametmp,"%s_denk%d.out",filename_prefix,k);
			std::cout << "writing hdf5 data to " << fnametmp << std::endl;
			hid_t file_id = open_hdf5_file(std::string(fnametmp));
			//hid_t group_id = create_group(file_id,std::string(suffix));
			this->density_k_[k].write_to_hdf5_density(file_id,dt,this->tid_);
			close_hdf5_file(file_id);
		}
	}

	// std::cout << "slice 2 " << this->tid_ << std::endl;
	if(print_k){
	  //Print full green for k=0 - test for compression
	  if(this->tid_map_[0]==this->tid_){
	    hid_t file_id = open_hdf5_file("Gk0_full");
	    hid_t group_id = create_group(file_id, "G");
	    green_k_[0].G_.write_to_hdf5(group_id);
	    group_id = create_group(file_id, "Sigma");
	    green_k_[0].Sigma_.write_to_hdf5(group_id);
	    close_group(group_id);
	    close_hdf5_file(file_id);
	  }
	  for(int k=0;k<this->nk_;k++){
	    if(this->tid_map_[k]==this->tid_){
	      // std::cout << "slice 3" << this->tid_ <<  std::endl;
	      char fnametmp[1000];
	      std::sprintf(fnametmp,"%s_k%d.out",filename_prefix,k);
	      std::cout << "writing hdf5 data to " << fnametmp << std::endl;
	      hid_t file_id = open_hdf5_file(std::string(fnametmp));
	      //hid_t group_id = create_group(file_id,std::string(suffix));
	      green_k_[k].write_to_hdf5_slices(file_id,dt,this->tid_);
	      close_hdf5_file(file_id);
	    }
	  }
	}
	
	if(this->tid_==this->tid_root_){
		// merge???
	}
}

template <class LATTICE>
void mpi_lattice_step_2b_deb<LATTICE>::print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k){
  assert(std::strlen(filename_prefix)<900);
	// NOTE on MPI:
	// not using parallel hdf5, so should write in serial, to ensure coherent file access
	// maybe a bit faster, depending on the network structure:write files for each k
	// std::cout << "slice 1" << std::endl;
	if(this->tid_==this->tid_root_){
		char fnametmp[1000];
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		hid_t group_id = create_group(file_id, "parm");
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		close_group(group_id); // End parameters
		// -- Green's functions
		group_id = create_group(file_id, "Gloc");
		Gloc_.write_to_hdf5_slices(group_id,dt);
		close_group(group_id);
		group_id = create_group(file_id, "Gloc_full");
		Gloc_.write_to_hdf5(group_id);
		close_group(group_id);
		group_id = create_group(file_id, "Gtavrel");
		Gloc_.write_to_hdf5_tavtrel(group_id,dt);
		close_group(group_id);
		close_hdf5_file(file_id);
	}
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			char fnametmp[1000];
			std::sprintf(fnametmp,"%s_denk%d.out",filename_prefix,k);
			std::cout << "writing hdf5 data to " << fnametmp << std::endl;
			hid_t file_id = open_hdf5_file(std::string(fnametmp));
			//hid_t group_id = create_group(file_id,std::string(suffix));
			this->density_k_[k].write_to_hdf5_density(file_id,dt,this->tid_);
			close_hdf5_file(file_id);
		}
	}

	// std::cout << "slice 2 " << this->tid_ << std::endl;
	if(print_k){
	  //Print full green for k=0 - test for compression
	  if(this->tid_map_[0]==this->tid_){
	    hid_t file_id = open_hdf5_file("Gk0_full");
	    hid_t group_id = create_group(file_id, "G");
	    green_k_[0].G_.write_to_hdf5(group_id);
	    group_id = create_group(file_id, "Sigma");
	    green_k_[0].Sigma_.write_to_hdf5(group_id);
	    close_group(group_id);
	    close_hdf5_file(file_id);
	  }
	  for(int k=0;k<this->nk_;k++){
	    if(this->tid_map_[k]==this->tid_){
	      // std::cout << "slice 3" << this->tid_ <<  std::endl;
	      char fnametmp[1000];
	      std::sprintf(fnametmp,"%s_k%d.out",filename_prefix,k);
	      std::cout << "writing hdf5 data to " << fnametmp << std::endl;
	      hid_t file_id = open_hdf5_file(std::string(fnametmp));
	      //hid_t group_id = create_group(file_id,std::string(suffix));
	      green_k_[k].write_to_hdf5_slices(file_id,dt,this->tid_);
	      close_hdf5_file(file_id);
	    }
	  }
	}
	
	if(this->tid_==this->tid_root_){
		// merge???
	}
}

// template <class LATTICE>
// void mpi_lattice_step_2b_optical<LATTICE>::print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k){
// }

template <class LATTICE>
void mpi_lattice_step_2b_optical<LATTICE>::print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k){
  assert(std::strlen(filename_prefix)<900);
	// NOTE on MPI:
	// not using parallel hdf5, so should write in serial, to ensure coherent file access
	// maybe a bit faster, depending on the network structure:write files for each k
	// std::cout << "slice 1: " << this->tid_ << std::endl;
	if(this->tid_==this->tid_root_){
		char fnametmp[1000];
		// std::cout << "slice 11 in: " << this->tid_ << std::endl;
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		// std::cout << "slice 111 in: " << this->tid_ << std::endl;
		hid_t group_id = create_group(file_id, "parm");
		// std::cout << "slice 112 in: " << this->tid_ << std::endl;
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		// std::cout << "slice 113 in: " << this->tid_ << std::endl;
		close_group(group_id); // End parameters
		// std::cout << "slice 12 in: " << this->tid_ << std::endl;
		// -- Green's functions
		group_id = create_group(file_id, "Gloc");
		Gloc_.write_to_hdf5_slices(group_id,dt);
		close_group(group_id);
		// std::cout << "slice 13 in: " << this->tid_ << std::endl;
		group_id = create_group(file_id, "Gloc_full");
		Gloc_.write_to_hdf5(group_id);
		close_group(group_id);
		// std::cout << "slice 14 in: " << this->tid_ << std::endl;
		group_id = create_group(file_id, "Gtavrel");
		Gloc_.write_to_hdf5_tavtrel(group_id,dt);
		// std::cout << "slice 15 in: " << this->tid_ << std::endl;
		close_group(group_id);
		close_hdf5_file(file_id);
	}
	// std::cout << "slice 2: " << this->tid_ << std::endl;
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			// std::cout << "slice 21: " << this->tid_ << " " << this->tid_map_[k]  << std::endl;
			char fnametmp[1000];
			// std::cout << "slice 22: " << this->tid_ << " " << this->tid_map_[k]  << std::endl;
			std::sprintf(fnametmp,"%s_denk%d.out",filename_prefix,k);
			// std::cout << "slice 23: " << this->tid_ << " " << this->tid_map_[k]  << std::endl;
			std::cout << "writing hdf5 data to " << fnametmp << std::endl;
			hid_t file_id = open_hdf5_file(std::string(fnametmp));
			//hid_t group_id = create_group(file_id,std::string(suffix));
			// std::cout << "slice 24: " << this->tid_ << " " << this->tid_map_[k]  << std::endl;
			this->density_k_[k].write_to_hdf5_density(file_id,dt,this->tid_);
			close_hdf5_file(file_id);
		}
	}

	// std::cout << "slice 2 " << this->tid_ << std::endl;
	if(print_k){
	  //Print full green for k=0 - test for compression
	  // if(this->tid_map_[0]==this->tid_){
	  //   hid_t file_id = open_hdf5_file("Gk0_full");
	  //   hid_t group_id = create_group(file_id, "G");
	  //   green_k_[0].G_.write_to_hdf5(group_id);
	  //   group_id = create_group(file_id, "Sigma");
	  //   green_k_[0].Sigma_.write_to_hdf5(group_id);
	  //   close_group(group_id);
	  //   close_hdf5_file(file_id);
	  // }
	  for(int k=0;k<this->nk_;k++){
	    if(this->tid_map_[k]==this->tid_){
	      // std::cout << "slice 3" << this->tid_ <<  std::endl;
	      char fnametmp[1000];
	      std::sprintf(fnametmp,"%s_k%d.out",filename_prefix,k);
	      std::cout << "writing hdf5 data to " << fnametmp << std::endl;
	      hid_t file_id = open_hdf5_file(std::string(fnametmp));
	      //hid_t group_id = create_group(file_id,std::string(suffix));
	      green_k_[k].write_to_hdf5_slices(file_id,dt,this->tid_);
	      close_hdf5_file(file_id);
	    }
	  }
	}
	
	if(this->tid_==this->tid_root_){
		// merge???
	}
}

template <class LATTICE>
void mpi_lattice_step_hubbard_2b<LATTICE>::print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k){
  assert(std::strlen(filename_prefix)<900);
	// NOTE on MPI:
	// not using parallel hdf5, so should write in serial, to ensure coherent file access
	// maybe a bit faster, depending on the network structure:write files for each k
	// std::cout << "slice 1: " << this->tid_ << std::endl;
	if(this->tid_==this->tid_root_){
		char fnametmp[1000];
		// std::cout << "slice 11 in: " << this->tid_ << std::endl;
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		// std::cout << "slice 111 in: " << this->tid_ << std::endl;
		hid_t group_id = create_group(file_id, "parm");
		// std::cout << "slice 112 in: " << this->tid_ << std::endl;
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		// std::cout << "slice 113 in: " << this->tid_ << std::endl;
		close_group(group_id); // End parameters
		// std::cout << "slice 12 in: " << this->tid_ << std::endl;
		// -- Green's functions
		group_id = create_group(file_id, "Gloc");
		Gloc_.write_to_hdf5_slices(group_id,dt);
		close_group(group_id);
		// std::cout << "slice 13 in: " << this->tid_ << std::endl;
		group_id = create_group(file_id, "Gloc_full");
		Gloc_.write_to_hdf5(group_id);
		close_group(group_id);
		// std::cout << "slice 14 in: " << this->tid_ << std::endl;
		group_id = create_group(file_id, "Gtavrel");
		Gloc_.write_to_hdf5_tavtrel(group_id,dt);
		// std::cout << "slice 15 in: " << this->tid_ << std::endl;
		close_group(group_id);
		close_hdf5_file(file_id);
	}
	// std::cout << "slice 2: " << this->tid_ << std::endl;
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			// std::cout << "slice 21: " << this->tid_ << " " << this->tid_map_[k]  << std::endl;
			char fnametmp[1000];
			// std::cout << "slice 22: " << this->tid_ << " " << this->tid_map_[k]  << std::endl;
			std::sprintf(fnametmp,"%s_denk%d.out",filename_prefix,k);
			// std::cout << "slice 23: " << this->tid_ << " " << this->tid_map_[k]  << std::endl;
			std::cout << "writing hdf5 data to " << fnametmp << std::endl;
			hid_t file_id = open_hdf5_file(std::string(fnametmp));
			//hid_t group_id = create_group(file_id,std::string(suffix));
			// std::cout << "slice 24: " << this->tid_ << " " << this->tid_map_[k]  << std::endl;
			this->density_k_[k].write_to_hdf5_density(file_id,dt,this->tid_);
			close_hdf5_file(file_id);
		}
	}

	// std::cout << "slice 2 " << this->tid_ << std::endl;
	if(print_k){
	  //Print full green for k=0 - test for compression
	  // if(this->tid_map_[0]==this->tid_){
	  //   hid_t file_id = open_hdf5_file("Gk0_full");
	  //   hid_t group_id = create_group(file_id, "G");
	  //   green_k_[0].G_.write_to_hdf5(group_id);
	  //   group_id = create_group(file_id, "Sigma");
	  //   green_k_[0].Sigma_.write_to_hdf5(group_id);
	  //   close_group(group_id);
	  //   close_hdf5_file(file_id);
	  // }
	  for(int k=0;k<this->nk_;k++){
	    if(this->tid_map_[k]==this->tid_){
	      // std::cout << "slice 3" << this->tid_ <<  std::endl;
	      char fnametmp[1000];
	      std::sprintf(fnametmp,"%s_k%d.out",filename_prefix,k);
	      std::cout << "writing hdf5 data to " << fnametmp << std::endl;
	      hid_t file_id = open_hdf5_file(std::string(fnametmp));
	      //hid_t group_id = create_group(file_id,std::string(suffix));
	      green_k_[k].write_to_hdf5_slices(file_id,dt,this->tid_);
	      close_hdf5_file(file_id);
	    }
	  }
	}
	
	if(this->tid_==this->tid_root_){
		// merge???
	}
}

template <class LATTICE>
void mpi_lattice_step_hubbard_RPA<LATTICE>::print_to_file_hdf5_slice(const char *filename_prefix,int dt,int print_k){
  assert(std::strlen(filename_prefix)<900);
	// NOTE on MPI:
	// not using parallel hdf5, so should write in serial, to ensure coherent file access
	// maybe a bit faster, depending on the network structure:write files for each k
	// std::cout << "slice 1: " << this->tid_ << std::endl;
	if(this->tid_==this->tid_root_){
		char fnametmp[1000];
		// std::cout << "slice 11 in: " << this->tid_ << std::endl;
		std::sprintf(fnametmp,"%s_local.out",filename_prefix);
		hid_t file_id = open_hdf5_file(std::string(fnametmp));
		// std::cout << "slice 111 in: " << this->tid_ << std::endl;
		hid_t group_id = create_group(file_id, "parm");
		// std::cout << "slice 112 in: " << this->tid_ << std::endl;
		store_int_attribute_to_hid(group_id, "nt", this->nt_); 
		store_int_attribute_to_hid(group_id, "ntau", this->ntau_);
		store_double_attribute_to_hid(group_id, "beta", this->beta_);
		store_double_attribute_to_hid(group_id, "h", this->h_);
		// std::cout << "slice 113 in: " << this->tid_ << std::endl;
		close_group(group_id); // End parameters
		// std::cout << "slice 12 in: " << this->tid_ << std::endl;
		// -- Green's functions
		group_id = create_group(file_id, "Gloc");
		Gloc_.write_to_hdf5_slices(group_id,dt);
		close_group(group_id);
		// std::cout << "slice 13 in: " << this->tid_ << std::endl;
		group_id = create_group(file_id, "Gloc_full");
		Gloc_.write_to_hdf5(group_id);
		close_group(group_id);
		// std::cout << "slice 14 in: " << this->tid_ << std::endl;
		group_id = create_group(file_id, "Gtavrel");
		Gloc_.write_to_hdf5_tavtrel(group_id,dt);
		// std::cout << "slice 15 in: " << this->tid_ << std::endl;
		close_group(group_id);
		close_hdf5_file(file_id);
	}
	// std::cout << "slice 2: " << this->tid_ << std::endl;
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			// std::cout << "slice 21: " << this->tid_ << " " << this->tid_map_[k]  << std::endl;
			char fnametmp[1000];
			// std::cout << "slice 22: " << this->tid_ << " " << this->tid_map_[k]  << std::endl;
			std::sprintf(fnametmp,"%s_denk%d.out",filename_prefix,k);
			// std::cout << "slice 23: " << this->tid_ << " " << this->tid_map_[k]  << std::endl;
			std::cout << "writing hdf5 data to " << fnametmp << std::endl;
			hid_t file_id = open_hdf5_file(std::string(fnametmp));
			//hid_t group_id = create_group(file_id,std::string(suffix));
			// std::cout << "slice 24: " << this->tid_ << " " << this->tid_map_[k]  << std::endl;
			this->density_k_[k].write_to_hdf5_density(file_id,dt,this->tid_);
			close_hdf5_file(file_id);
		}
	}

	// std::cout << "slice 2 " << this->tid_ << std::endl;
	if(print_k){
	  //Print full green for k=0 - test for compression
	  // if(this->tid_map_[0]==this->tid_){
	  //   hid_t file_id = open_hdf5_file("Gk0_full");
	  //   hid_t group_id = create_group(file_id, "G");
	  //   green_k_[0].G_.write_to_hdf5(group_id);
	  //   group_id = create_group(file_id, "Sigma");
	  //   green_k_[0].Sigma_.write_to_hdf5(group_id);
	  //   close_group(group_id);
	  //   close_hdf5_file(file_id);
	  // }
	  for(int k=0;k<this->nk_;k++){
	    if(this->tid_map_[k]==this->tid_){
	      // std::cout << "slice 3" << this->tid_ <<  std::endl;
	      char fnametmp[1000];
	      std::sprintf(fnametmp,"%s_k%d.out",filename_prefix,k);
	      std::cout << "writing hdf5 data to " << fnametmp << std::endl;
	      hid_t file_id = open_hdf5_file(std::string(fnametmp));
	      //hid_t group_id = create_group(file_id,std::string(suffix));
	      green_k_[k].write_to_hdf5_slices(file_id,dt,this->tid_);
	      close_hdf5_file(file_id);
	    }
	  }
	}
	
	if(this->tid_==this->tid_root_){
		// merge???
	}
}

