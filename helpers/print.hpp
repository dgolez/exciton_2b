#pragma once

#include "inclusions.hpp"


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
	// -- Impurity parameters
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
void kpoint_green<LATTICE>::write_to_hdf5_slices(hid_t group_id,int outputfrequency,int tid){
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
	G_.write_to_hdf5_slices(group_id,"G",outputfrequency);
	G_.write_to_hdf5_tavtrel(group_id,"Gtavrel",outputfrequency);
	Sigma_.write_to_hdf5_slices(group_id, "Sigma",outputfrequency);
	chi_.write_to_hdf5_slices(group_id, "chi",outputfrequency);
	chi_.write_to_hdf5_tavtrel(group_id, "chitavtrel",outputfrequency);
	P_.write_to_hdf5_slices(group_id, "P",outputfrequency);
}

template <class LATTICE>
void kpoint_density<LATTICE>::write_to_hdf5_density(hid_t group_id,int outputfrequency,int tid){
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
void kpoint_green<LATTICE>::write_to_hdf5_slices(const char *filename,int outputfrequency,int tid){
	hid_t file_id = open_hdf5_file(filename);
	this->write_to_hdf5_slices(file_id,outputfrequency,tid);
	close_hdf5_file(file_id);
}

template <class LATTICE>
void mpi_lattice_step_optical<LATTICE>::read_from_file_hdf5(int tstp,const char *filename_prefix,int kt){
	// TODO !!
}

template <class LATTICE>
void mpi_lattice_step_2b_optical<LATTICE>::read_from_file_hdf5(int tstp,const char *filename_prefix,int kt){
	// TODO !!
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
		X_.write_to_hdf5(group_id,"X");
		Pi_.write_to_hdf5(group_id,"Pi");
		close_group(group_id);
		close_hdf5_file(file_id);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

template <class LATTICE>
void mpi_lattice_step_optical<LATTICE>::print_to_file_hdf5_slice(const char *filename_prefix,int outputfrequency,int print_k){
}

template <class LATTICE>
void mpi_lattice_step_2b_optical<LATTICE>::print_to_file_hdf5_slice(const char *filename_prefix,int outputfrequency,int print_k){
  assert(std::strlen(filename_prefix)<900);
	// NOTE on MPI:
	// not using parallel hdf5, so should write in serial, to ensure coherent file access
	// maybe a bit faster, depending on the network structure:write files for each k
	// std::cout << "slice 1: " << this->tid_ << std::endl;
	if(this->tid_==this->tid_root_){
		//output_local.out          ← file_id
		//└── Gloc                  ← group_id
		//	└── datasets written by Gloc_.write_to_hdf5(...)

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
		Gloc_.write_to_hdf5_slices(group_id,outputfrequency);
		close_group(group_id);
		group_id = create_group(file_id, "Gloc_full");
		Gloc_.write_to_hdf5(group_id);
		close_group(group_id);
		group_id = create_group(file_id, "Gtavrel");
		Gloc_.write_to_hdf5_tavtrel(group_id,outputfrequency);
		close_group(group_id);
		close_hdf5_file(file_id);
	}
	for(int k=0;k<this->nk_;k++){
		if(this->tid_map_[k]==this->tid_){
			char fnametmp[1000];
			std::sprintf(fnametmp,"%s_denk%d.out",filename_prefix,k);
			std::cout << "writing hdf5 data to " << fnametmp << std::endl;
			hid_t file_id = open_hdf5_file(std::string(fnametmp));
			this->density_k_[k].write_to_hdf5_density(file_id,outputfrequency,this->tid_);
			close_hdf5_file(file_id);
		}
	}

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
	      green_k_[k].write_to_hdf5_slices(file_id,outputfrequency,this->tid_);
	      close_hdf5_file(file_id);
	    }
	  }
	}
	
	if(this->tid_==this->tid_root_){
		// merge???
	}
}