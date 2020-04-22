#include "../Headers/Catalog.h"

// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
void Catalog::set_params_catalog(Params params)
{
  this->Output_directory=params._Output_directory();
  this->params=params;
  this->box.masskernel=params._iMAS_Y();
  this->box.Lbox=params._Lbox();
  this->box.Nft=params._Nft();
  this->box.NGRID=static_cast<ULONG>(params._Nft()*params._Nft()*params._Nft());
  this->box.d1=params._d1();
  this->box.d2=params._d2();
  this->box.d3=params._d3();
  this->box.min1=params._xllc();
  this->box.min2=params._yllc();
  this->box.min3=params._zllc();

}

// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //



// ***********************************************************************************************************  //
#ifdef _USE_VELOCITIES_
void Catalog::read_catalog_bin(ULONG n, string filex, string filey, string filez,string filevx, string filevy, string filevz)
#else
void Catalog::read_catalog_bin(ULONG n, string filex, string filey, string filez)
#endif
{
  this->xgal.clear();
  this->ygal.clear();
  this->zgal.clear();
  this->Pxgal.clear();
  this->Pygal.clear();
  this->Pzgal.clear();
  this->Mass.clear();
  this->property.clear();
  



#ifdef _GET_INTERPOLATED_FIELDS_FROM_BIN_FILES_
  cout<<""<<endl;
  So.message_screen("Reading binary files for pos and vel of DM, and interpolating into a grid.");

 #ifdef _GET_VEL_FIELD_
  vector<real_prec>vx_field(this->box.NGRID,0);
  vector<real_prec>vy_field(this->box.NGRID,0);
  vector<real_prec>vz_field(this->box.NGRID,0);
#endif
#ifdef _GET_CIC_DENS_FIELD_
  vector<real_prec>dens_field_cic(this->box.NGRID,0);
#endif
#ifdef _GET_NGP_DENS_FIELD_
  vector<real_prec>dens_field_ngp(this->box.NGRID,0);
#endif
#ifdef _GET_TSC_DENS_FIELD_
  vector<real_prec>dens_field_tsc(this->box.NGRID,0);
#endif

  string bin_file;
  int N_bin_files=32; // this number comes from the number of hdf5 files converted to binary
  vector<ULONG>N_parts(N_bin_files,0);
  string path_hdfc="/net/deimos/scratch1/balaguera/data/Numerics/Minerva_DM/z1/001/";
  ifstream np; np.open(path_hdfc+"number_of_particles.txt");
  for(int ik=0; ik<N_bin_files;++ik)
    np>>N_parts[ik];
  np.close();
      
  for(int ik=0; ik<N_bin_files;++ik)
    {
            
      bin_file=path_hdfc+"minerva_dm_x_file"+to_string(ik)+".dat";
      ULONG Nn=N_parts[ik]; 

      So.message_screen("N particles =", Nn);
      this->xgal.clear();
      this->xgal.shrink_to_fit();
      this->xgal.resize(Nn,0);
      this->File.read_array(bin_file, this->xgal);

      bin_file=path_hdfc+"minerva_dm_y_file"+to_string(ik)+".dat";
      this->ygal.clear();
      this->ygal.shrink_to_fit();
      this->ygal.resize(Nn,0);
      this->File.read_array(bin_file, this->ygal);

      bin_file=path_hdfc+"minerva_dm_z_file"+to_string(ik)+".dat";
      this->zgal.clear();
      this->zgal.shrink_to_fit();
      this->zgal.resize(Nn,0);
      this->File.read_array(bin_file, this->zgal);

#ifdef _GET_VEL_FIELD_
      bin_file=path_hdfc+"minerva_dm_velx_file"+to_string(ik)+".dat";
      this->Pxgal.clear();
      this->Pxgal.shrink_to_fit();
      this->Pxgal.resize(Nn,0);
      this->File.read_array(bin_file, this->Pxgal);
      
      bin_file=path_hdfc+"minerva_dm_vely_file"+to_string(ik)+".dat";
      this->Pygal.clear();
      this->Pygal.shrink_to_fit();
      this->Pygal.resize(Nn,0);
      this->File.read_array(bin_file, this->Pygal);

      bin_file=path_hdfc+"minerva_dm_velz_file"+to_string(ik)+".dat";
      this->Pzgal.clear();
      this->Pzgal.shrink_to_fit();
      this->Pzgal.resize(Nn,0);
      this->File.read_array(bin_file, this->Pzgal);
#endif


#ifdef _GET_NGP_DENS_FIELD_
      So.message_screen("Interpolating density to grid ngp");
      getDensity_NGP(this->box.Nft, this->box.Nft, this->box.Nft,this->box.Lbox, this->box.Lbox,this->box.Lbox, this->box.d1, this->box.d2, this->box.d3,0,0,0,this->xgal,this->ygal,this->zgal,this->zgal, dens_field_ngp,false);
      So.DONE();
#endif

#ifdef _GET_CIC_DENS_FIELD_
      So.message_screen("Interpolating density to grid CIC");
      getDensity_CIC(this->box.Nft, this->box.Nft, this->box.Nft,this->box.Lbox, this->box.Lbox,this->box.Lbox, this->box.d1, this->box.d2, this->box.d3,0,0,0,this->xgal,this->ygal,this->zgal,this->zgal, dens_field_cic,false);
      So.DONE();
#endif

#ifdef _GET_TSC_DENS_FIELD_
      So.message_screen("Interpolating density to grid TSC");
      getDensity_TSC(this->box.Nft, this->box.Nft, this->box.Nft,this->box.Lbox, this->box.Lbox,this->box.Lbox, this->box.d1, this->box.d2, this->box.d3,0,0,0,this->xgal,this->ygal,this->zgal,this->zgal, dens_field_tsc,false);
      So.DONE();
#endif
#ifdef _GET_VEL_FIELD_
      So.message_screen("Interpolating Vx to grid");
      getDensity_CIC(this->box.Nft, this->box.Nft, this->box.Nft,this->box.Lbox, this->box.Lbox,this->box.Lbox, this->box.d1, this->box.d2, this->box.d3,0,0,0,this->xgal,this->ygal,this->zgal,this->Pxgal,vx_field,true);
      So.DONE();
      So.message_screen("Interpolating Vy to grid");
      getDensity_CIC(this->box.Nft, this->box.Nft, this->box.Nft,this->box.Lbox, this->box.Lbox,this->box.Lbox, this->box.d1, this->box.d2, this->box.d3,0,0,0,this->xgal,this->ygal,this->zgal,this->Pygal,vy_field,true);
      So.DONE();
      So.message_screen("Interpolating Vz to grid");
      getDensity_CIC(this->box.Nft, this->box.Nft, this->box.Nft,this->box.Lbox, this->box.Lbox,this->box.Lbox, this->box.d1, this->box.d2, this->box.d3,0,0,0,this->xgal,this->ygal,this->zgal,this->Pzgal,vz_field,true);
      So.DONE();
#endif
  }

#ifdef _GET_VEL_FIELD_
  ULONG ec=0;
#pragma omp parallel for reduction(+:ec)
  for(ULONG i=0;i<this->box.NGRID;++i)
    {
      if(dens_field[i]!=0)
	{
          vx_field[i]/=static_cast<real_prec>(dens_field_cic[i]);
          vy_field[i]/=static_cast<real_prec>(dens_field_cic[i]);
          vz_field[i]/=static_cast<real_prec>(dens_field_cic[i]);
	}
      else
      {
              // this is not extrictly speaking right. We would have to assingn the velocity of the closest dm particle if the cell is empty.
              // However, for practical reasons, we will be using the velocity field in conjunction with the density (e.g., kinetic energy)
	ec++;
	vx_field[i]=0;
	vy_field[i]=0;
	vz_field[i]=0;
       }
    }	  
  So.message_screen("Number of empty cell=", ec);
#endif



  string field=this->Output_directory;
#ifdef _GET_NGP_DENS_FIELD_
  File.write_array(field+"Minerva_DM_density_NGP_Real1_z1_N"+to_string(this->box.Nft),dens_field_ngp);
  dens_field_ngp.clear();
  dens_field_ngp.shrink_to_fit();
#endif

#ifdef _GET_CIC_DENS_FIELD_
  File.write_array(field+"Minerva_DM_density_CICp_Real1_z1_N"+to_string(this->box.Nft),dens_field_cic);
  dens_field_cic.clear();
  dens_field_cic.shrink_to_fit();
#endif

#ifdef _GET_TSC_DENS_FIELD_
  File.write_array(field+"Minerva_DM_density_TSC_Real1_z1_N"+to_string(this->box.Nft),dens_field_tsc);
  dens_field_tsc.clear();
  dens_field_tsc.shrink_to_fit();
#endif

#ifdef _GET_VEL_FIELD_
  File.write_array(field+"Minerva_DM_vx_CIC_Real1_z1_N"+to_string(this->box.Nft),vx_field);
  File.write_array(field+"Minerva_DM_vy_CIC_Real1_z1_N"+to_string(this->box.Nft),vy_field);
  File.write_array(field+"Minerva_DM_vz_CIC_Real1_z1_N"+to_string(this->box.Nft),vz_field);
  vx_field.clear();
  vy_field.clear();
  vz_field.clear();
#endif


#else  
  vector<real_prec>dummy(n,0);
  
  this->File.read_array(filex, dummy);
  for(ULONG i=0;i<n;++i)
    this->xgal.push_back(static_cast<double>(dummy[i]));
  
  this->File.read_array(filey, dummy);
  for(ULONG i=0;i<n;++i)
    this->ygal.push_back(static_cast<double>(dummy[i]));
  
  this->File.read_array(filez, dummy);
  for(ULONG i=0;i<n;++i)
    this->zgal.push_back(static_cast<double>(dummy[i]));

#ifdef _USE_VELOCITIES_
  this->File.read_array(filevx, dummy);
  for(ULONG i=0;i<n;++i)
    this->Pxgal.push_back(static_cast<double>(dummy[i]));
  
  this->File.read_array(filevy, dummy);
  for(ULONG i=0;i<n;++i)
    this->Pygal.push_back(static_cast<double>(dummy[i]));
  
  this->File.read_array(filevz, dummy);
  for(ULONG i=0;i<n;++i)
    this->Pzgal.push_back(static_cast<double>(dummy[i]));

#endif
  
  dummy.clear();
  dummy.shrink_to_fit();


#endif
  
  this->NOBJS=n; 
  this->property.resize(this->NOBJS, 0); 
}

// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //

#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
void Catalog::read_catalog(string input_file, real_prec prop_min)
#elif defined (_USE_MASS_BINS_PK_)
void Catalog::read_catalog(string input_file, real_prec prop_min, real_prec prop_max)
#endif
{



// The input parameter min_mass cann be also VMAX min, according to preproc definitions


  int NTHREADS = omp_get_max_threads();
  omp_set_num_threads(NTHREADS);


#ifdef _USE_MULTISCALE_MASS_ASSIGNMENT_
  this->Mass_threshold=pow(MASS_THRESHOLD, exponent_mass_tracer); // to be deprecated

#ifdef _USE_MULTISCALE_LEVEL_1_
  this->Mass_threshold_multi_scale_1=pow(MASS_THRESHOLD_MULTI_SCALE_1, exponent_mass_tracer);
#endif

#ifdef _USE_MULTISCALE_LEVEL_2_
  this->Mass_threshold_multi_scale_2=pow(MASS_THRESHOLD_MULTI_SCALE_2, exponent_mass_tracer);
#endif
#ifdef _USE_MULTISCALE_LEVEL_3_
  this->Mass_threshold_multi_scale_3=pow(MASS_THRESHOLD_MULTI_SCALE_3, exponent_mass_tracer);
#endif
#ifdef _USE_MULTISCALE_LEVEL_4_
  this->Mass_threshold_multi_scale_4=pow(MASS_THRESHOLD_MULTI_SCALE_4, exponent_mass_tracer);
#endif
#endif



  int i_x, i_y, i_z, i_vx, i_vy, i_vz, i_mass, i_weight, i_mean_density, i_sf, i_vmax;
  if(this->type_of_object=="DM")
  {
    i_x= this->params._i_coord1_dm();
    i_y= this->params._i_coord2_dm();
    i_z= this->params._i_coord3_dm();
    i_vx= this->params._i_v1_dm();
    i_vy= this->params._i_v2_dm();
    i_vz= this->params._i_v3_dm();
    i_mass = this->params._i_mass_dm();
  }
 else if(this->type_of_object=="TRACER" || this->type_of_object=="TRACER_REF" || this->type_of_object=="TRACER_MOCK" || this->type_of_object=="TRACER_MOCK_ONLY_COORDS")
  {
     i_x= this->params._i_coord1_g();
     i_y= this->params._i_coord2_g();
     i_z= this->params._i_coord3_g();
     i_vx= this->params._i_v1_g();
     i_vy= this->params._i_v2_g();
     i_vz= this->params._i_v3_g();
     i_mass = this->params._i_mass_g();
     i_vmax = this->params._i_vmax_g();
     i_weight = this->params._i_weight1_g();
     i_mean_density= this->params._i_mean_density_g();
     i_sf = this->params._i_sf_g();
  }

  else if(this->type_of_object=="RANDOM")
   {
      i_x= this->params._i_coord1_r();
      i_y= this->params._i_coord2_r();
      i_z= this->params._i_coord3_r();
      i_weight = this->params._i_weight1_r();
      i_mean_density= this->params._i_mean_density_r();
   }


  int i_obsevable=0;
  real_prec units_observable=1;
#ifdef _USE_MASS_AS_OBSERVABLE_
  i_obsevable=i_mass;
  units_observable=this->params._MASS_units();
#elif defined _USE_VMAX_AS_OBSERVABLE_
  i_obsevable=i_vmax;
#endif


  int i_setting=0;
  real_prec units_settings=1;
#ifdef _SET_CAT_WITH_MASS_CUT_
  i_setting=i_mass;
  units_settings=this->params._MASS_units();
#elif defined _SET_CAT_WITH_VMAX_CUT_
  i_setting=i_vmax;
#endif

#ifdef _USE_VMAX_AS_OBSERVABLE_
  if(this->type_of_object!="TRACER_MOCK_ONLY_COORDS")
   if(i_vmax<0)
      {
      So.message_warning("No information for Vmax. Check .ini parameter");
      So.message_warning("Code ends here. Sorry");
      exit(0);
      }
#endif


  vector<real_prec> prop;
  ULONG NLINES = this->File.read_file(input_file, prop,NTHREADS);
  
  this->NCOLS=(static_cast<ULONG>(prop.size()/NLINES));

/*
  ofstream tea;
  So.message_screen("Writting new catalog to file ",this->Output_directory+"newcat.txt");
  tea.open(this->Output_directory+"newcat.txt");
  for(ULONG i=0;i<NLINES;++i)
    if(prop[i_mass+i*NCOLS]*params._MASS_units()>= prop_min)
      tea<<prop[i_x+i*NCOLS]<<"  "<<prop[i_y+i*NCOLS]<<"  "<<prop[i_z+i*NCOLS]<<"  "<<prop[i_mass+i*NCOLS]<<"  "<<prop[i_vmax+i*NCOLS]<<"  "<<prop[i_sf+i*NCOLS]<<endl;

  tea.close();
  So.DONE();
  exit(0);
*/


  //************************** JUST COUNTING ABOVE LIMITS OR INTERVALS*********************************
  // The minimum mass (M or VMAX) defines the number of used tracers!!!!!!!!!!!!!

  ULONG count_new_nobj=0;

  if(this->type_of_object!="TRACER_MOCK_ONLY_COORDS")  //This applies for mocks with mass
   {
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:count_new_nobj)
#endif
  for(ULONG i=0;i<NLINES;++i)
    {
       real_prec obser=prop[i_setting+i*this->NCOLS]*units_settings;
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
       if(obser >= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
       if(obser >= prop_min && obser < prop_max)
#endif         
          count_new_nobj++;
    }


#ifdef _SET_CAT_WITH_MASS_CUT_
  So.message_screen("Minimum mass requested =", prop_min, "Ms/h");
#ifdef _USE_MASS_BINS_PK_
  So.message_screen("Maximum mass requested =", prop_max, "Ms/h");
  So.message_screen("Number of in the mass bin =", count_new_nobj);
#else
  So.message_screen("Number of tracers above mass limit =", count_new_nobj);
#endif

#elif defined (_SET_CAT_WITH_VMAX_CUT_)
  So.message_screen("Minimum VMAX requested =", prop_min, "Km/s");
#ifdef _USE_MASS_BINS_PK_
  So.message_screen("Maximum VMAX requested =", prop_max, "Km/s");
  So.message_screen("Number of in the mass bin =", count_new_nobj);
#else
  So.message_screen("Number of tracers above VMAX limit =", count_new_nobj);
#endif

#endif
  }
 else  // Otherwilse, if type is mock with only coords, the total number of objects is the number of lines in the mock
{
      count_new_nobj=NLINES;
      prop_min=-1.0;
  }


  this->NOBJS=count_new_nobj;
  this->Halo.resize(this->NOBJS);

  ULONG iN=0;

  cout<<endl;
  So.message_screen("Getting grid-ID from tracer coordinates");
  for(ULONG i=0;i<NLINES;++i)
    {
      real_prec mass=0;
      if(this->type_of_object!="TRACER_MOCK_ONLY_COORDS")
        mass=prop[i_setting+i*this->NCOLS]*units_settings;

#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
         if(mass>= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
          if(mass>= prop_min && mass< prop_max)
#endif
          {
           real_prec x=prop[i_x+i*NCOLS];
           real_prec y=prop[i_y+i*NCOLS];
           real_prec z=prop[i_z+i*NCOLS];
           this->Halo[iN].coord1=x;
           this->Halo[iN].coord2=y;
           this->Halo[iN].coord3=z;
           this->Halo[iN].GridID = grid_ID(&box, x,y,z);
           this->Halo[iN].observed=true;
#ifdef _USE_NUMBER_OF_SATELLITES_
           this->Halo[iN].n_satellites=prop[i_sf+i*NCOLS];
#endif
#ifdef _USE_VMAX_TRACERS_
         if(this->type_of_object!="TRACER_MOCK_ONLY_COORDS")
            this->Halo[iN].vmax=prop[i_vmax+i*NCOLS];
#endif
           iN++;
	}
    }
So.DONE();


#ifdef _USE_MULTISCALE_MASS_ASSIGNMENT_

 if(type_of_object=="TRACER_REF")
 {

   // These numbers areonly to be determined from the abundance of the reference
   this->N_masses_0   =0;
   this->N_masses_1   =0;
   this->N_masses_2   =0;
   this->N_masses_3   =0;
   this->N_masses_4   =0;

  So.message_screen("Computing numbers for multi-scale mass assgnment");
  for(ULONG i=0;i<NLINES;++i)
    {
       real_prec mass= pow(prop[i_obsevable+i*this->NCOLS]*units_observable,exponent_mass_tracer);

#ifdef _USE_MULTISCALE_LEVEL_4_
       if(mass >= this->Mass_threshold_multi_scale_4)
         this->N_masses_4 ++;
#endif

#ifdef _USE_MULTISCALE_LEVEL_3_  // ensures level 1 is always used
#ifdef _USE_MULTISCALE_LEVEL_4_
      if(mass >= this->Mass_threshold_multi_scale_3 && mass < this->Mass_threshold_multi_scale_4)
#elif !defined (_USE_MULTISCALE_LEVEL_4_)
      if(mass >= this->Mass_threshold_multi_scale_3)
#endif
        this->N_masses_3++;
#endif


#ifdef _USE_MULTISCALE_LEVEL_2_  // ensures level 1 is always used
#ifdef _USE_MULTISCALE_LEVEL_3_
     if(mass >= this->Mass_threshold_multi_scale_2 && mass < this->Mass_threshold_multi_scale_3)
#elif !defined (_USE_MULTISCALE_LEVEL_3_) && defined (_USE_MULTISCALE_LEVEL_4_)
     if(mass >= this->Mass_threshold_multi_scale_2 && mass < this->Mass_threshold_multi_scale_4)
#elif !defined (_USE_MULTISCALE_LEVEL_3_) && !defined (_USE_MULTISCALE_LEVEL_4_)
     if(mass >= this->Mass_threshold_multi_scale_2)
#endif
       this->N_masses_2++;
#endif




#ifdef _USE_MULTISCALE_LEVEL_1_
#if defined (_USE_MULTISCALE_LEVEL_2_)
       if(mass >= this->Mass_threshold_multi_scale_1 && mass < this->Mass_threshold_multi_scale_2)
#elif !defined (_USE_MULTISCALE_LEVEL_2_) && defined (_USE_MULTISCALE_LEVEL_3_)
       if(mass >= this->Mass_threshold_multi_scale_1 && mass < this->Mass_threshold_multi_scale_3)
#elif !defined (_USE_MULTISCALE_LEVEL_2_) && !defined (_USE_MULTISCALE_LEVEL_3_) && defined (_USE_MULTISCALE_LEVEL_4_)
       if(mass >= this->Mass_threshold_multi_scale_1 && mass < this->Mass_threshold_multi_scale_4)
#elif !defined (_USE_MULTISCALE_LEVEL_2_) && !defined (_USE_MULTISCALE_LEVEL_3_) && !defined (_USE_MULTISCALE_LEVEL_4_)
       if(mass >= this->Mass_threshold_multi_scale_1)
#endif
           this->N_masses_1++;
#endif

    }
    So.DONE();
/*
    So.message_screen("N_M1= ", N_masses_1);
    So.message_screen("N_M2= ", N_masses_2);
    So.message_screen("N_M3= ", N_masses_3);
    So.message_screen("N_M4= ", N_masses_4);
*/

 }

#endif





  // If information of weights is available, then
  if(i_weight>0)
  {
      So.message_screen("Retrieving weighs from tracers");
      iN=0;
      for(ULONG i=0;i<NLINES;++i)
        {
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
              if(prop[i_setting+i*this->NCOLS]*units_settings>= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
          if(prop[i_setting+i*this->NCOLS]*units_settings >= prop_min && prop[i_setting+i*this->NCOLS]*units_settings< prop_max)
#endif
            {
              this->Halo[iN].weight1=prop[i_weight+i*NCOLS];
              iN++;
            }
        }
    So.DONE();
  }


  if(i_mean_density>0)
  {
      So.message_screen("Retrieving number density from tracers");
      iN=0;
      for(ULONG i=0;i<NLINES;++i)
        {
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
         if(prop[i_setting+i*this->NCOLS]*units_settings>= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
           if(prop[i_setting+i*this->NCOLS]*units_settings >= prop_min && prop[i_setting+i*this->NCOLS]*units_settings< prop_max)
#endif
           {
              this->Halo[iN].mean_density=prop[i_mean_density+i*NCOLS];
              iN++;
            }
        }
      So.DONE();
  }




// Here we have to implement an if statement in order to select objects within a "property" bin or cut
// before assigning the cat to the structure this->Halo. The quantity this->NOBS must be then
// the remaining number of tracers after the cuts.

#ifdef _USE_VELOCITIES_TRACERS_
  if(i_vx>0) //Best to evaluate only one if than NOBS if's
    {
      iN=0;
      for(ULONG i=0;i<NLINES;++i)
        if(prop[i_mass+i*NCOLS]*params._MASS_units()>= mass_min)
	  {
	    this->Halo[iN].vel1=prop[i_vx+i*NCOLS];
	    iN++;
	  }
    }
  
  if(i_vy>0)
    {
      iN=0;
      for(ULONG i=0;i<NLINES;++i)
	if(prop[i_mass+i*NCOLS]*params._MASS_units()>= mass_min)
	  {
	    this->Halo[iN].vel2=prop[i_vy+i*NCOLS];
	    iN++;
	  }
    }


  if(i_vz>0)
    {
      iN=0;
      for(ULONG i=0;i<NLINES;++i)
	if(prop[i_mass+i*NCOLS]*params._MASS_units()>= mass_min)
	  {
	    this->Halo[iN].vel3=prop[i_vy+i*NCOLS];
	    iN++;
	  }
    }
#endif


  // ********************************************************************************************************************************
 if(this->type_of_object!="TRACER_MOCK_ONLY_COORDS")
#ifdef _USE_MASS_TRACERS_
   if(i_mass>0) //means there is info on the mass
     {
      real_prec mean_m=0;
      iN=0;
      vector<real_prec> aux;
      for(ULONG i=0;i<NLINES;++i)

#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
          if(prop[i_setting+i*this->NCOLS]*units_settings >= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
          if(prop[i_setting+i*this->NCOLS]*units_settings >= prop_min && prop[i_setting+i*this->NCOLS]*units_settings< prop_max)
#endif
              {

                  real_prec proper=prop[i_mass+i*NCOLS];
#ifdef _USE_LOG_MASS_
            mean_m+=log10(proper);
            this->Halo[iN].mass=log10(proper);
#elif defined (_USE_X_ROOT_MASS_)
            this->Halo[iN].mass=pow(proper, exponent_mass_tracer);
            mean_m+=pow(proper, exponent_mass_tracer);
            aux.push_back(pow(proper, exponent_mass_tracer));
#else
	    this->Halo[iN].mass=prop[i_mass+i*NCOLS];
            mean_m+=proper;
            aux.push_back(proper);
#endif
            iN++;
	  }

      mean_m/=static_cast<real_prec>(this->NOBJS);

      So.message_screen("Minimum mass in tracer catalog =", get_min<real_prec>(aux), "Ms/h");
      So.message_screen("Maximum mass in tracer catalog =", get_max<real_prec>(aux), "Ms/h");
      aux.clear(); aux.shrink_to_fit();

#ifdef _USE_LOG_MASS_
      this->So.message_screen("Mean mass of tracer catalogue =", pow(10,mean_m)*this->params._MASS_units(), "Ms/h");
#else
      this->So.message_screen("Mean mass of tracer catalogue =", mean_m*this->params._MASS_units(), "Ms/h");
#endif

      }
  else
   {
     cout<<endl;
     So.message_screen("No Information on the tracer mass available in catalog, according to input parameters");
   }
#endif
   // ********************************************************************************************************************************

 if(this->type_of_object!="TRACER_MOCK_ONLY_COORDS")
#ifdef _USE_VMAX_TRACERS_
   if(i_vmax>0) //means there is info on the mass
     {
      real_prec mean_m=0;
      iN=0;
      vector<real_prec> aux;
      for(ULONG i=0;i<NLINES;++i)

#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
          if(prop[i_setting+i*this->NCOLS]*units_settings >= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
          if(prop[i_setting+i*this->NCOLS]*units_settings >= prop_min && prop[i_setting+i*this->NCOLS]*units_settings< prop_max)
#endif
              {
                real_prec proper=prop[i_vmax+i*NCOLS]*units_observable;
                this->Halo[iN].vmax=proper;
                aux.push_back(proper);
                mean_m+=proper;
               iN++;
          }
      mean_m/=static_cast<real_prec>(this->NOBJS);
      So.message_screen("Minimum VMAX in tracer catalog =", get_min<real_prec>(aux), "km/s");
      So.message_screen("Maximum VMAX in tracer catalog =", get_max<real_prec>(aux), "km/s");
      aux.clear(); aux.shrink_to_fit();

#ifdef _USE_LOG_MASS_
      this->So.message_screen("Mean mass of tracer catalogue =", pow(10,mean_m)*this->params._MASS_units(), "km/s");
#else
      this->So.message_screen("Mean VMAX of tracer catalogue =", mean_m*this->params._MASS_units(), "km/s");
#endif

      }
  else
   {
     cout<<endl;
     So.message_screen("No Information on the tracer mass available in catalog, according to input parameters");
   }
#endif
   // ********************************************************************************************************************************




  
#ifdef _USE_SAT_FRACTION_

  if(i_sf>0)
    {
      iN=0;
      int iN_nosat=0;
      for(ULONG i=0;i<NLINES;++i)
      if(prop[i_mass+i*NCOLS]*params._MASS_units()>= prop_min)
       {
          this->Halo[iN].number_sub_structures=static_cast<int>(prop[i_sf+i*NCOLS]);
	  iN++;
          if(prop[i_sf+i*NCOLS] < 1)
              iN_nosat++;
       }
  double mean_m=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:mean_m)
#endif
  for(ULONG i=0;i<this->NOBJS;++i)
    mean_m+=static_cast<double>(this->Halo[i].number_sub_structures);

  mean_m/=static_cast<double>(this->NOBJS);
  this->So.message_screen("Mean number of satellites per host =", mean_m);
  this->So.message_screen("Fraction of centrals without sub-structure =", 100.0*static_cast<double>(iN_nosat)/static_cast<double>(this->NOBJS),"%");

      }
#endif




  this->mean_number_density=static_cast<real_prec>(this->NOBJS)/pow(this->box.Lbox,3);

  cout<<endl;
  
}
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //

void Catalog::get_density_field_grid(string prop, string output_file)
{
  
  So.message_screen("Interpolating on a grid using MAS", this->params._masskernel());
  vector<real_prec> deltaTR_counts(this->box.NGRID,0);
  MAS_NEW(&box,this->Halo, _COUNTS_, deltaTR_counts);

  vector<real_prec> deltaTR;

  if(_COUNTS_ == prop)
    this->File.write_array(output_file, deltaTR_counts);

  else
    {
      
      deltaTR.resize(this->box.NGRID,0);
      MAS_NEW(&box,this->Halo, prop, deltaTR);
      
      if(_MASS_ == prop)
        {
         So.message_screen("Interpolating tracer on a grid weighting by mass");
/*           for(ULONG i=0 ; i< this->box.NGRID; ++i)
             if(deltaTR_counts[i]!=0)
                deltaTR[i]=static_cast<double>(deltaTR[i])/static_cast<double>(deltaTR_counts[i]); // Mean mass in cells
             else
               deltaTR[i]=0;*/
        }
      else if (_SAT_FRACTION_==prop)
         So.message_screen("Interpolating Number of satellites in a grid");

      this->File.write_array(output_file, deltaTR);
    }
}

// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //

void Catalog::get_density_field_grid(string prop, vector<real_prec>&out)
{
    this->params.masskernel=0;
    if(prop==_COUNTS_)
    So.message_screen("Computing number counts in cells");
    else if (prop==_MASS_)
    So.message_screen("Computing total tracer-mass in cells");
    MAS_NEW(&box,this->Halo,prop, out);
}


// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //

void Catalog::define_property_bins()
{

  //this->NMBINS=this->params._NMASSbins_mf();
#ifdef MBINS
  So.message_screen("Defining bins for abundance:");
#elif defined MCUTS
  So.message_screen("Defining cuts for abundnace:");
#endif

  
  if(this->params.i_mass_g >0)
    {

      this->NMBINS = this->params._NMASSbins_mf() == 1 ? 1 : this->params._NMASSbins_mf();
      So.message_screen("Generating ",this->NMBINS," bin(s) (zeroth -or  first bin- is the full sample) to tracer-mass");
      
      // to account for the zero bin which includes the full sample, except if one asks for only one bin, in which case the full is assumed
      
      real_prec MMin;
      real_prec MMax;

#ifdef _MASS_LOG_
      MMin = this->params._LOGMASSmin();
      MMax = this->params._LOGMASSmax();
      this->So.message_screen("Minimum Mass", pow(10,MMin), "Ms/h");
      this->So.message_screen("Maximum Mass", pow(10,MMax), "Ms/h");
#else
      MMin= pow(10,this->params._LOGMASSmin());
      MMin= this->params._LOGMASSmin();
      this->So.message_screen("Minimum Mass", MMin, "Ms/h");
      this->So.message_screen("Maximum Mass", MMax, "Ms/h");
#endif

      this->logdeltaM=(MMax-MMin)/(static_cast<real_prec>(this->params._NMASSbins_mf()));
      this->MBmin.clear();this->MBmin.shrink_to_fit();

#ifdef _MASS_LOG_
      for(int i=0;i<this->NMBINS;++i)
	this->MBmin.push_back(pow(10,MMin+i*logdeltaM));
#else
      for(int i=0;i<this->NMBINS;++i)
	this->MBmin.push_back(MMin+i*logdeltaM);
#endif
      
      
#ifdef _MASS_LOG_
      this->MBmax.clear();this->MBmax.shrink_to_fit();
      this->MBin.clear();this->MBin.shrink_to_fit();

#ifdef MBINS
      for(int i=0;i<this->NMBINS;++i)
	this->MBmax.push_back(pow(10,MMin+(i+1)*logdeltaM));

      for(int i=0;i<this->NMBINS;++i)
        this->MBin.push_back(this->MBmin[i]+(i+0.5)*(this->MBmax[i]-this->MBmin[i])/static_cast<double>(this->params._NMASSbins_mf()));

#elif defined MCUTS
      for(int i=0;i<this->NMBINS;++i)
	this->MBmax.push_back(pow(10, MMax));
#endif

#else

#ifdef MBINS
      for(int i=0;i<NMBINS;++i)
	this->MBmax.push_back(MMin+(i)*logdeltaM);
#elif defined MCUTS
      for(int i=0;i<NMBINS;++i)
	this->MBmax.push_back(MMax);
#endif
      
#endif
      
      So.DONE();
    }
  else
    {
      So.message_screen("No halo mass info available for type ", this->type_of_object);
      this->NMBINS=1;
      this->logdeltaM=1;
      this->MBmin.resize(NMBINS,-100.0); //  if no mass is available, in the loop over the masses we set mass = 1,
      this->MBmax.resize(NMBINS,100.0);  //
    }




#ifdef _USE_VMAX_AS_OBSERVABLE_


  if(this->params.i_vmax_g >0)
    {
      this->NMBINS = this->params._NMASSbins_mf() == 1 ? 1 : this->params._NMASSbins_mf();
      So.message_screen("Generating ",this->NMBINS," bin(s) (zeroth -or  first bin- is the full sample) to tracer-vmax");

      real_prec MMin;
      real_prec MMax;

#ifdef _MASS_LOG_
      MMin= log10(this->params._VMAXmin());
      MMax= log10(this->params._VMAXmax());
      this->So.message_screen("Minimum VMAX", pow(10,MMin), "km/s");
      this->So.message_screen("Maximum VMAX", pow(10,MMax), "km/s");
#else
      MMin= pow(10,this->params._VMAXmin());
      MMin= this->params._VMAXMASSmin();
      this->So.message_screen("Minimum Mass", MMin, "Ms/h");
      this->So.message_screen("Maximum Mass", MMax, "Ms/h");
#endif
      
      this->logdeltaVMAX=(MMax-MMin)/(static_cast<real_prec>(this->params._NMASSbins_mf()));
      
      this->VMAXBmin.clear();this->VMAXBmin.shrink_to_fit();
#ifdef _MASS_LOG_
      for(int i=0;i<this->NMBINS;++i)
        this->VMAXBmin.push_back(pow(10,MMin+i*logdeltaVMAX));
#else
      for(int i=0;i<this->NMBINS;++i)
        this->VMAXBmin.push_back(MMin+i*logdeltaVMAX);
#endif


#ifdef _MASS_LOG_
      this->VMAXBmax.clear();this->VMAXBmax.shrink_to_fit();
      this->VMAXBin.clear();this->VMAXBin.shrink_to_fit();

#ifdef MBINS
      for(int i=0;i<this->NMBINS;++i)
        this->VMAXBmax.push_back(pow(10,MMin+(i+1)*logdeltaVMAX));

      for(int i=0;i<this->NMBINS;++i)
        this->VMAXBin.push_back(this->VMAXBmin[i]+(i+0.5)*(this->VMAXBmax[i]-this->VMAXBmin[i])/static_cast<double>(this->params._NMASSbins_mf()));

#elif defined MCUTS
      for(int i=0;i<this->NMBINS;++i)
        this->VMAXBmax.push_back(pow(10, MMax));
#endif

#else

#ifdef MBINS
      for(int i=0;i<NMBINS;++i)
        this->VMAXBmax.push_back(MMin+(i)*logdeltaVMAX);
#elif defined MCUTS
      for(int i=0;i<NMBINS;++i)
        this->VMAXBmax.push_back(MMax);
#endif

#endif

      So.DONE();
    }
  else
    {
      So.message_screen("No VMAX info available for type ", this->type_of_object);
      this->NMBINS=1;
      this->logdeltaVMAX=1;
      this->VMAXBmin.resize(NMBINS,-100.0); //  if no mass is available, in the loop over the masses we set mass = 1,
      this->VMAXBmax.resize(NMBINS,100.0);  //
    }




#endif














  
  
}

// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
 void Catalog::get_property_function(string file)
 {

#ifdef accum
   bool accumulate=true;
 #else
   bool accumulate=false;
#endif

     //Her we get by default the mass function and if vmax is used, the v,ax_function

   this->define_property_bins();
   
   if(this->params.i_mass_g>0)
     {
       So.message_screen("Measuring mass function");
       real_prec lm_min=this->params._LOGMASSmin();
       real_prec lm_max=this->params._LOGMASSmax();
       real_prec units_observable_m=this->params._MASS_units();
       this->mass_function.resize(this->params._NMASSbins_mf());
       vector<int >mf_counts(this->params._NMASSbins_mf(),0);
       
       
       int im;
       for(ULONG i=0;i<this->NOBJS;++i)
	 {
	   real_prec property=log10(this->Halo[i].mass)+log10(units_observable_m);
#ifndef accum
	   if(property >=lm_min && property<=lm_max)
	     {
#endif
	       im=get_bin(property,lm_min,this->params._NMASSbins_mf(),this->logdeltaM,accumulate);
	       mf_counts[im]++;
#ifndef accum
	     }
#endif
	 }
       ULONG naux=0;
#pragma omp parallel for reduction(+:naux)
       for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
	 naux+=mf_counts[i];
       So.message_screen("Number of objects used in mass function n(M) = ", naux);
       So.message_screen("Number of objects EXCLUDED in n(M) estimation = ", this->NOBJS-naux);
       
       
       string file_m=file+"_mass";
       ofstream mout (file_m.c_str());
       So.message_screen("Writing mass function written in file ", file_m);
       
       for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
	 {
	   real_prec mint=(this->MBmax[i]-this->MBmin[i]);
	   this->mass_function[i]=static_cast<real_prec>(mf_counts[i])/(pow(this->params._Lbox(),3)*mint);
	   mout<<this->MBmin[i]+(i+0.5)*(this->MBmax[i]-this->MBmin[i])/static_cast<double>(this->params._NMASSbins_mf())<<"\t"<<this->mass_function[i]<<"   "<<mf_counts[i]<<endl;
	 }
       mf_counts.clear(); mf_counts.shrink_to_fit();
       mout.close();
       So.DONE();
       
       vector<real_prec>mcount(this->params._NMASSbins_mf(),0);
       for(int i=0; i< this->params._NMASSbins_mf(); ++i)
	 for(int j=i; j< this->params._NMASSbins_mf(); ++j)
	   mcount[i]+=this->mass_function[j]*(this->MBmax[j]-this->MBmin[j]);
       
       file_m+="_cumulative";
       mout.open(file_m.c_str());
       
       for(int i=0; i< this->params._NMASSbins_mf(); ++i)
	 mout<<this->MBin[i]<<"\t"<<mcount[i]<<endl;
       mout.close();
       So.message_screen("Cumulative Mass function written in file ", file_m);
       So.DONE();
       mcount.clear();mcount.shrink_to_fit();
     }
   
   
   if(this->params.i_vmax_g>0)
     {
        int im=0;
#ifdef _USE_VMAX_AS_OBSERVABLE_
       So.message_screen("Measuring VMAX function");
       real_prec v_min=log10(this->params._VMAXmin());
       real_prec v_max=log10(this->params._VMAXmax());
       real_prec units_observable_v=1;
       this->vmax_function.resize(this->params._NMASSbins_mf());
       vector<int >v_counts(this->params._NMASSbins_mf(),0);
       for(ULONG i=0;i<this->NOBJS;++i)
	 {
	   real_prec property=log10(this->Halo[i].vmax*units_observable_v);
#ifndef accum
	   if(property >=v_min && property<=v_max)
	     {
#endif
	       im=get_bin(property,v_min,this->params._NMASSbins_mf(),this->logdeltaVMAX,accumulate);
	       v_counts[im]++;
#ifndef accum
	     }
#endif
	 }
       
       int naux=0;
#pragma omp parallel for reduction(+:naux)
       for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
	 naux+=v_counts[i];
       So.message_screen("Number of objects used in mass function n(V) = ", naux);
       So.message_screen("Number of objects EXCLUDED in n(V) estimation = ", this->NOBJS-naux);
       string filev=file+"_vmax";
       
       So.message_screen("Writting VMAX-function in file ", filev);
       ofstream mout;
       mout.open(filev.c_str());
       for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
	 {
	   real_prec mint=(this->VMAXBmax[i]-this->VMAXBmin[i]);
	   this->vmax_function[i]=static_cast<real_prec>(v_counts[i])/(pow(this->params._Lbox(),3)*mint);
	   mout<<this->VMAXBmin[i]+(i+0.5)*(this->VMAXBmax[i]-this->VMAXBmin[i])/static_cast<double>(this->params._NMASSbins_mf())<<"\t"<<this->vmax_function[i]<<"   "<<v_counts[i]<<endl;
	 }
       mout.close();
       v_counts.clear();v_counts.shrink_to_fit();
       
       vector<real_prec>Vcount(this->params._NMASSbins_mf(),0);
       for(int i=0; i< this->params._NMASSbins_mf(); ++i)
	 for(int j=i; j< this->params._NMASSbins_mf(); ++j)
	   Vcount[i]+=this->vmax_function[j]*(this->VMAXBmax[j]-this->VMAXBmin[j]);
       
       filev+="_cumulative";
       mout.open(filev.c_str());
       for(int i=0; i< this->params._NMASSbins_mf(); ++i)
	 mout<<this->VMAXBin[i]<<"\t"<<Vcount[i]<<endl;
       mout.close();
       Vcount.clear();Vcount.shrink_to_fit();
       
       So.message_screen("Cumulative Mass function written in file ", filev);
#endif
     }

   
   So.DONE();
      
 }

// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
 // ***********************************************************************************************************  //
  // ***********************************************************************************************************  //
  // ***********************************************************************************************************  //
  // ***********************************************************************************************************  //
  // ***********************************************************************************************************  //
void Catalog:: get_distribution_reduced_mass_in_cell(){
    int NTHREADS = omp_get_max_threads();
    omp_set_num_threads(NTHREADS);

    So.message_screen("Measuring distribution of reduced mass in cells", MAXIMUM_DISTANCE_EXCLUSION, "Mpc/h");

  // In order to gt ths miminum distance for each particle,
  // we nend to allocate for eachcell the coordintas of the particles belonging to it,
  // as we did when collapsing randoms:

  vector<s_cell_info> cell_info_tr(this->box.NGRID);

#pragma omp parallel for
  for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      ULONG ID=this->Halo[i].GridID;
      cell_info_tr[ID].posx_p.push_back(this->Halo[i].coord1);
      cell_info_tr[ID].posy_p.push_back(this->Halo[i].coord2);
      cell_info_tr[ID].posz_p.push_back(this->Halo[i].coord3);
      cell_info_tr[ID].mass.push_back(this->Halo[i].mass);
    }

//  vector<int>dist_reduced_mass(N_BINS_REDUCED_MASS, 0);

  real_prec dmin=MINIMUM_DISTANCE_EXCLUSION;
  real_prec dmax=MAXIMUM_DISTANCE_EXCLUSION;


  vector<vector<int>> dist_reduced_mass(N_BINS_REDUCED_MASS,vector<int>(N_BINS_DIST_EX, 0));

  for(ULONG id=0; id<this->box.NGRID ;++id) // loop over cells
    if(cell_info_tr[id].posx_p.size()>=2)//nos aseguramos de pasar por celdas que tengan al menos dos tracers para los cuales determinar la distancia. Si sólo hay uno, o ninguno, la separación es cero.
      for(int i=0;i <cell_info_tr[id].posx_p.size();++i) // loop overtracers in that cell
        for(int j=i+1;j < cell_info_tr[id].posx_p.size() ; ++j) // loop over ttracer in that very same cell
         {
           real_prec dist=  sqrt(pow(cell_info_tr[id].posx_p[i]-cell_info_tr[id].posx_p[j],2)+pow(cell_info_tr[id].posy_p[i]-cell_info_tr[id].posy_p[j],2)+pow(cell_info_tr[id].posz_p[i]-cell_info_tr[id].posz_p[j],2));
           int index_dist =get_bin(dist,dmin, N_BINS_DIST_EX,(dmax-dmin)/static_cast<double>(N_BINS_DIST_EX),true);
           real_prec mu= cell_info_tr[id].mass[i]*cell_info_tr[id].mass[j]/(cell_info_tr[id].mass[i]+cell_info_tr[id].mass[j]);
           int index_mu =get_bin(mu,MIN_REDUCED_MASS,N_BINS_REDUCED_MASS,DELTA_REDUCED_MASS,true);
           dist_reduced_mass[index_mu][index_dist]++;
         }

  for(int id=0;id<N_BINS_DIST_EX;++id)
      {
        string file=this->Output_directory+"reduced_mass_dist_"+this->type_of_object+"_dbin"+to_string(id)+".txt";
          So.message_screen("Writting to file ", file);
      ofstream sal;
      sal.open(file.c_str());
      sal<<"#Bin info:     Dmin = "<<MINIMUM_DISTANCE_EXCLUSION+(id)*(dmax-dmin)/static_cast<double>(N_BINS_DIST_EX)<<"  Dmax = "<<MINIMUM_DISTANCE_EXCLUSION+(id+1.0)*(dmax-dmin)/static_cast<double>(N_BINS_DIST_EX)<<endl;
      for(int i=0;i<dist_reduced_mass.size();++i)
       sal<<MIN_REDUCED_MASS+(i+0.5)*DELTA_REDUCED_MASS<<" "<<dist_reduced_mass[i][id]<<endl;
    sal.close();
   So.DONE();
      }
  
  
  
}
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //

//#define _test_ms_


// This function computes the minimum separation between pairs in one cell
// Allocates the result in a class member container min_separation_in_cell[ID]=
void Catalog::get_min_separation_in_cell()
{
  
  int NTHREADS = omp_get_max_threads();
  omp_set_num_threads(NTHREADS);
  
  So.message_screen("Measuring minimum separation in cells (Mpc/h) ...");
  So.message_screen("Current type is ", this->type_of_object);
  cout<<endl;
  
  // In order to gt ths miminum distance for each particle,
  // we nend to allocate for eachcell the coordintas of the particles belonging to it,
  // as we did when collapsing randoms:
  
  vector<s_cell_info> cell_info_tr(this->box.NGRID);
  

  So.message_screen("Getting positions inside cells:");
  if("TRACER_MOCK_ONLY_COORDS" != this->type_of_object)
    {
      for(ULONG i=0; i<this->NOBJS; ++i) //loop over the "observed objects", i.e, with cuts already set
	{
          ULONG ID=this->Halo[i].GridID;
	  cell_info_tr[ID].posx_p.push_back(this->Halo[i].coord1);
	  cell_info_tr[ID].posy_p.push_back(this->Halo[i].coord2);
	  cell_info_tr[ID].posz_p.push_back(this->Halo[i].coord3);
#ifdef _test_ms_
#ifdef _USE_VMAX_AS_OBSERVABLE_
	  cell_info_tr[ID].property.push_back(this->Halo[i].vmax);
#elif defined _USE_MASS_AS_OBSERVABLE_
	  cell_info_tr[ID].property.push_back(this->Halo[i].mass);
#endif // endif use_vmax_as _obs
#endif // endif _test_ms_
	}
    }
  else
    {
      for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed objects", i.e, with cuts already set
	{
          ULONG ID=this->Halo[i].GridID;
	  cell_info_tr[ID].posx_p.push_back(this->Halo[i].coord1);
	  cell_info_tr[ID].posy_p.push_back(this->Halo[i].coord2);
	  cell_info_tr[ID].posz_p.push_back(this->Halo[i].coord3);
	}
    }
  So.DONE();
  
  
  
  
#ifdef _test_ms_
  int n_calc=1; // number of measurements of min sep_dist, 1 of no mass bins or cuts through multilevel is desired
  vector<real_prec> min_prop;
  min_prop.push_back(0);
  vector<int> level;
  level.push_back(0);
  
#ifdef _USE_MULTISCALE_LEVEL_1_
  if(this->type_of_object!="TRACER_MOCK_ONLY_COORDS")
    {
      min_prop.push_back(this->Mass_threshold_multi_scale_1);
      level.push_back(1);
      n_calc++;
    }
#endif

#ifdef _USE_MULTISCALE_LEVEL_2_
  if(this->type_of_object!="TRACER_MOCK_ONLY_COORDS")
    {
      min_prop.push_back(this->Mass_threshold_multi_scale_2);
      level.push_back(2);
      n_calc++;
    }
#endif

#ifdef _USE_MULTISCALE_LEVEL_3_
  if(this->type_of_object!="TRACER_MOCK_ONLY_COORDS")
    {
      min_prop.push_back(this->Mass_threshold_multi_scale_3);
      level.push_back(3);
      n_calc++;
    }
#endif

#ifdef _USE_MULTISCALE_LEVEL_4_
  if(this->type_of_object!="TRACER_MOCK_ONLY_COORDS")
    {
      min_prop.push_back(this->Mass_threshold_multi_scale_4);
      level.push_back(4);
      n_calc++;
    }
#endif

#endif  // end of _test_ms_
  
   
  vector<int>dist_min_sep;
  
#ifdef _test_ms_
  for(int ik=n_calc-1 ;ik>=0; ik--) // hacia atrás, de modo que el último que haga sea el general, que es el que se usa
    {
#endif
      
      this->min_separation_in_cell.clear();
      this->min_separation_in_cell.shrink_to_fit();
      this->min_separation_in_cell.resize(this->box.NGRID,0);
      
      // If the tracer has poperties, then ask for them when computing separations. If not, else below
      if("TRACER_MOCK_ONLY_COORDS" != this->type_of_object)
	{
	  
	  
#ifdef _test_ms_
#if defined _USE_MULTISCALE_LEVEL_1_ || defined _USE_MULTISCALE_LEVEL_2_ || defined _USE_MULTISCALE_LEVEL_3_ || defined _USE_MULTISCALE_LEVEL_4_
#ifdef _USE_VMAX_AS_OBSERVABLE_
	  So.message_screen("Getting Min separations for tracers with Vmax > ", min_prop[ik], "km/s: ");
#elif _USE_VMAX_AS_OBSERVABLE_
	  So.message_screen("Getting Min separations for tracers with M > ", min_prop[ik]), "Ms/h : ";
#endif
#endif
#endif
	  
	  for(ULONG id=0; id<this->box.NGRID ;++id) // loop over cells
	    {
	      real_prec aux_min_distance_a;
	      if(cell_info_tr[id].posx_p.size()>=2)//nos aseguramos de pasar por celdas que tengan al menos dos tracers para los cuales determinar la distancia. Si sólo hay uno, o ninguno, la separación es cero.
		{
		  aux_min_distance_a=10.0;
		  for(int i=0;i <cell_info_tr[id].posx_p.size();++i) // loop overtracers in that cell
		    {
#ifdef _test_ms_
#if defined _USE_MULTISCALE_LEVEL_1_ || defined _USE_MULTISCALE_LEVEL_2_ || defined _USE_MULTISCALE_LEVEL_3_ || defined _USE_MULTISCALE_LEVEL_4_
		      if(cell_info_tr[id].property[i]>min_prop[ik])
			{
#endif
#endif
			  real_prec aux_min_distance_b;
			   for(int j=i+1 ; j < cell_info_tr[id].posx_p.size() ; ++j) // loop over ttracer in that very same cell
			     {
#ifdef _test_ms_
#if defined _USE_MULTISCALE_LEVEL_1_ || defined _USE_MULTISCALE_LEVEL_2_ || defined _USE_MULTISCALE_LEVEL_3_ || defined _USE_MULTISCALE_LEVEL_4_
			       if(cell_info_tr[id].property[j]>min_prop[ik])
				 {
#endif
#endif
				   real_prec dist=  sqrt(pow(cell_info_tr[id].posx_p[i]-cell_info_tr[id].posx_p[j],2)+pow(cell_info_tr[id].posy_p[i]-cell_info_tr[id].posy_p[j],2)+pow(cell_info_tr[id].posz_p[i]-cell_info_tr[id].posz_p[j],2));
				   aux_min_distance_b=min(dist, aux_min_distance_a);
				   aux_min_distance_a=aux_min_distance_b;
#ifdef _test_ms_
#if defined _USE_MULTISCALE_LEVEL_1_ || defined _USE_MULTISCALE_LEVEL_2_ || defined _USE_MULTISCALE_LEVEL_3_ || defined _USE_MULTISCALE_LEVEL_4_
				 }
#endif
#endif
			     }
#ifdef _test_ms_
#if defined _USE_MULTISCALE_LEVEL_1_ || defined _USE_MULTISCALE_LEVEL_2_ || defined _USE_MULTISCALE_LEVEL_3_ || defined _USE_MULTISCALE_LEVEL_4_
			 }
#endif
#endif
		     }
#ifdef _USE_LOG_MSIC_
		   this->min_separation_in_cell[id]=aux_min_distance_a>0? log10(aux_min_distance_a):MIN_SEP_IN_CELLS;
#else
		   this->min_separation_in_cell[id]=aux_min_distance_a;
#endif
		 }
	       else
		 this->min_separation_in_cell[id]=MIN_SEP_IN_CELLS;
	     }
	 }
       else // This is menat for tracers with only coordinates
	 {
	   for(ULONG id=0; id<this->box.NGRID ;++id) // loop over cells
	     {
	       real_prec aux_min_distance_a;
	       if(cell_info_tr[id].posx_p.size()>=2)//nos aseguramos de pasar por celdas que tengan al menos dos tracers para los cuales determinar la distancia. Si sólo hay uno, o ninguno, la separación es cero.
		 {
		   aux_min_distance_a=100.0;
		   for(int i=0;i <cell_info_tr[id].posx_p.size();++i) // loop overtracers in that cell
		     {
		       real_prec aux_min_distance_b;
		       for(int j=i+1 ; j < cell_info_tr[id].posx_p.size() ; ++j) // loop over ttracer in that very same cell
			 {
			   real_prec dist=  sqrt(pow(cell_info_tr[id].posx_p[i]-cell_info_tr[id].posx_p[j],2)+pow(cell_info_tr[id].posy_p[i]-cell_info_tr[id].posy_p[j],2)+pow(cell_info_tr[id].posz_p[i]-cell_info_tr[id].posz_p[j],2));
			   aux_min_distance_b=min(dist, aux_min_distance_a);
			   aux_min_distance_a=aux_min_distance_b;
			 }
		     }
#ifdef _USE_LOG_MSIC_
		   this->min_separation_in_cell[id]=aux_min_distance_a>0? log10(aux_min_distance_a):MIN_SEP_IN_CELLS;
#else
		   this->min_separation_in_cell[id]=aux_min_distance_a;
#endif
		 }
	       else
		 this->min_separation_in_cell[id]=MIN_SEP_IN_CELLS;
	     }
	 }
       
       cell_info_tr.clear();cell_info_tr.shrink_to_fit();
       dist_min_sep.clear();
       dist_min_sep.shrink_to_fit();
       dist_min_sep.resize(N_BINS_MIN_SEP_IN_CELLS,0);
       
       for(ULONG id=0; id<this->box.NGRID ;++id) // loop over cells
	 {
	   if(this->min_separation_in_cell[id]>MIN_SEP_IN_CELLS)
	     {
	       int index =get_bin(this->min_separation_in_cell[id],0,N_BINS_MIN_SEP_IN_CELLS,DELTA_MIN_SEP,true);
	       dist_min_sep[index]++;
	     }
	 }
       So.DONE();
       
       this->min_halo_separation=get_min<real_prec>(this->min_separation_in_cell, true);
       
#ifdef _USE_LOG_MSIC_
       So.message_screen("Min value of log min_sep_in_cell= ", this->min_halo_separation);
       So.message_screen("Max value of log min_sep_in_cell= ", get_max<real_prec>(this->min_separation_in_cell));
#else
       So.message_screen("Min value of min_sep_in_cell= ", this->min_halo_separation, "Mpc /h");
       So.message_screen("Max value of min_sep_in_cell= ", get_max<real_prec>(this->min_separation_in_cell), "Mpc /h");
#endif
       
       
       ofstream sal;
       string fileo;
#ifdef _test_ms_
       if(this->type_of_object!="TRACER_MOCK_ONLY_COORDS")
         fileo=this->Output_directory+"min_sep_in_cells_dist_"+this->type_of_object+"_level"+to_string(level[ik])+".txt";
        else
          fileo=this->Output_directory+"min_sep_in_cells_dist_"+this->type_of_object+"_level"+to_string(0)+".txt";
#else
       fileo=this->Output_directory+"min_sep_in_cells_dist_"+this->type_of_object+".txt";
#endif
       So.message_screen("Writting distribution of minimum separations within cells");
       So.message_screen("in file ", fileo);
       sal.open(fileo.c_str());
       for(int i=0;i<dist_min_sep.size();++i)
	 if(dist_min_sep[i]>0.)
	   sal<<(i+0.5)*DELTA_MIN_SEP<<" "<<dist_min_sep[i]<<endl;
       sal.close();
       So.DONE();
       
       //   File.write_array(this->Output_directory+"min_sep", this->min_separation_in_cell);
       
#ifdef _test_ms_
     }
#endif
   
 }



// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// This function identifies the set of tracers around each tracer within a radius of  SCALE_MAX_N_NEIGHBOURS Mpc/h
// The output is allocated in a class member container Number_of_tracers with dimension Nobjects
void Catalog::get_neighbour_tracers(vector<s_nearest_cells>&nearest_cells_to_cell)
{

   int NTHREADS = omp_get_max_threads();
   omp_set_num_threads(NTHREADS);


  So.message_screen("Identifying tracers in neighbouring cells");


   real_prec Rscale = SCALE_MAX_N_NEIGHBOURS;
// we need to cunstrunc the tree.For each grid point, we can build a vector containing the 8 or 29 nearest cells
   
   int max_neigh_per_dim = 1+2*N_CELLS_BACK_FORTH; // Number of neighbouring cells pr dimension, including the same itself;
   int N_Neigh_cells=pow(max_neigh_per_dim,3); // total number of cells to explore around a cell
   
   So.message_screen("Maximum separation probed =",  max_neigh_per_dim*sqrt(3.0)*this->box.Lbox/static_cast<real_prec>(this->box.Nft),"Mpc /h");
   
   vector<s_cell_info> cell_info_tr(this->box.NGRID);
   
#pragma omp parallel for
   for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
     {
       ULONG ID=this->Halo[i].GridID;
       cell_info_tr[ID].posx_p.push_back(this->Halo[i].coord1);
       cell_info_tr[ID].posy_p.push_back(this->Halo[i].coord2);
       cell_info_tr[ID].posz_p.push_back(this->Halo[i].coord3);
       cell_info_tr[ID].gal_index.push_back(i);//allocate the index of each tracer in that cell
     }
   So.DONE();
   
#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
  real_prec aux_min_distance_a;
  real_prec aux_min_distance_b;
#endif


  // This arrawy will allcopate the number of neigh or the bin in distance in which the minn distance to pair falls
  this->Number_of_neighbours.clear();
  this->Number_of_neighbours.shrink_to_fit();
  this->Number_of_neighbours.resize(this->NOBJS,0);
  
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
  So.message_screen("Computing number of neighbours for each tracer");
  So.message_screen("in spheres of radius ", Rscale, "Mpc/h");
#endif


#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
  So.message_screen("Computing minimum separation to other tracer for each tracer");
#endif
  for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      real_prec x_coord=this->Halo[i].coord1;
      real_prec y_coord=this->Halo[i].coord2;
      real_prec z_coord=this->Halo[i].coord3;
      ULONG ID=this->Halo[i].GridID;

#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
      vector<real_prec> min_distances_v(N_Neigh_cells,0);
#endif


      for(int j=0;j<N_Neigh_cells; j++) // loop sobre las celdas cercanas incluida la celda donde esta la particula i
        {
          ULONG ID_NEIGH = nearest_cells_to_cell[ID].close_cell[j];
          real_prec factor_bc_x=nearest_cells_to_cell[ID].bc_x[j]*this->box.Lbox;
          real_prec factor_bc_y=nearest_cells_to_cell[ID].bc_y[j]*this->box.Lbox;
          real_prec factor_bc_z=nearest_cells_to_cell[ID].bc_z[j]*this->box.Lbox;

#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
          aux_min_distance_b=BIG_NUMBER;
          real_prec min_distance_b;
#endif

          for(int k=0;k< cell_info_tr[ID_NEIGH].posx_p.size(); ++k) //loop sobre las particulas en las celdas cercanas
            {
              ULONG index_gal=cell_info_tr[ID_NEIGH].gal_index[k];
              if(index_gal!=i) //this simply avoids to get a 0 in the distance
                {

                  real_prec distx = x_coord - (cell_info_tr[ID_NEIGH].posx_p[k]+factor_bc_x);
                  real_prec disty = y_coord - (cell_info_tr[ID_NEIGH].posy_p[k]+factor_bc_y);
                  real_prec distz = z_coord - (cell_info_tr[ID_NEIGH].posz_p[k]+factor_bc_z);
                  real_prec dist  = sqrt(pow(distx,2)+pow(disty,2)+pow(distz,2));
#ifdef _USE_LOG_DIST_
                  dist=log10(dist);
#endif

#if defined (_USE_NUMBER_OF_NEIGHBOURS_) || defined (_USE_LOCAL_CLUSTERING_)
                  if(dist<Rscale)
                      this->Number_of_neighbours[i]++;
#endif
#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
                  min_distance_b=min(aux_min_distance_b,dist);
                  aux_min_distance_b=min_distance_b;
#endif
                  }
            }
#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
              min_distances_v[j]=aux_min_distance_b;
#endif
        }

#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
        aux_min_distance_a=BIG_NUMBER;
        real_prec min_distance_a;
        for(int j=0;j<N_Neigh_cells; j++)
          if(min_distances_v[j]>0)
             {
               min_distance_a=min(min_distances_v[j],aux_min_distance_a);
               aux_min_distance_a=min_distance_a;
              }
#endif


#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
      int da_bin = get_bin(aux_min_distance_a,MIN_OF_MIN_SEPARATION,N_BINS_MIN_DIST_TO_NEI,MAX_OF_MIN_SEPARATION/static_cast<real_prec>(N_BINS_MIN_DIST_TO_NEI),true);
      this->Number_of_neighbours[i]=da_bin;
#endif
   }
  So.DONE();
  
  
#ifdef _USE_LOCAL_CLUSTERING_
  So.message_screen("Computiong local clustering:");

  real_prec expected_number_of_tracers = 4.0*M_PI*pow(SCALE_MAX_N_NEIGHBOURS,3)*this->mean_number_density/3.0;
  So.message_screen("Expected number of particles in the lcoal volume:", expected_number_of_tracers);

  //  ofstream tea; tea.open("local_cl.txt");
  
#pragma omp parallel for  
  for(ULONG i=0;i<this->NOBJS ;++i)
    {
      real_prec local_xi=log10(static_cast<real_prec>(this->Number_of_neighbours[i])/expected_number_of_tracers+1);
      int index_xi=get_bin(local_xi,MIN_LOCAL_CLUSTERING,N_BINS_MIN_DIST_TO_NEI, (MAX_LOCAL_CLUSTERING-MIN_LOCAL_CLUSTERING)/static_cast<real_prec>(N_BINS_MIN_DIST_TO_NEI),true);
      this->Number_of_neighbours[i]=index_xi;
    }
	
  //	tea.close();
  So.DONE();
#endif

   
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
   int max_nh=get_max<int>(this->Number_of_neighbours);
   So.message_screen("Maximum number of neighbours around a tracer = ", max_nh);
   So.message_screen("Using ", N_NEIGHBOURS_MAX);
#endif


   So.DONE();

}



// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //

void Catalog::get_distribution_min_separations(string data, vector<s_nearest_cells>&nearest_cells_to_cell)
{


// we need to cunstrunc the tree.For each grid point, we can build a vector containing the 8 or 29 nearest cells
   int NTHREADS = omp_get_max_threads();
   omp_set_num_threads(NTHREADS);

   int max_neigh_per_dim = 1+2*N_CELLS_BACK_FORTH; // Number of neighbouring cells pr dimension, including the same itself;

   int N_Neigh_cells=pow(max_neigh_per_dim,3); // total number of cells to explore around a cell

   So.message_screen("Measuring distribution of pairs for ", this->type_of_object);
   cout<<endl;
   So.message_screen("   maximum separation probed =",  max_neigh_per_dim*sqrt(3.0)*this->box.Lbox/static_cast<real_prec>(this->box.Nft),"Mpc /h");
   cout<<endl;
   
   struct s_cell_info{
     vector<real_prec> posx_p;
     vector<real_prec> posy_p;
     vector<real_prec> posz_p;
     vector<ULONG> gal_index;
   };
   
    //DEFINE VECTOR TO ALLOCATE THE NEIGHBOR CELLS OF EACH CELL
   
   // In order to gt ths miminum distance for each particle,
   // we nend to allocate for eachcell the coordintas of the particles belonging to it,
   // as we did when collapsing randoms:
   
   So.message_screen("Identifying tracers in cells");
   
   vector<s_cell_info> cell_info_tr(this->box.NGRID);
  
#pragma omp parallel for
  for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      ULONG ID=this->Halo[i].GridID;
      cell_info_tr[ID].posx_p.push_back(this->Halo[i].coord1);
      cell_info_tr[ID].posy_p.push_back(this->Halo[i].coord2);
      cell_info_tr[ID].posz_p.push_back(this->Halo[i].coord3);
      cell_info_tr[ID].gal_index.push_back(i);//allocate the index of each tracer in that cell
    }
  So.DONE();
  
  // vector<real_prec>min_separations(this->NOBJS,0);
  //  vector<bool>marked(this->NOBS,false);
  
  
  int N_bins_dist=25;
  int N_mass_bins=200;
  real_prec max_dist=25.0;
  vector<vector<real_prec>> min_distance_dist(N_bins_dist,vector<real_prec>(N_mass_bins,0));
  
//  vector<vector<real_prec>> distance_dist(N_bins_dist,vector<real_prec>(N_mass_bins*N_mass_bins,0));


  vector<real_prec> distance_dist(N_bins_dist*N_mass_bins*N_mass_bins,0);

  real_prec lm_min=this->params._LOGMASSmin();
  real_prec lm_max=this->params._LOGMASSmax();
    
  this->logdeltaM=(lm_max-lm_min)/static_cast<real_prec>(N_mass_bins);
  
  real_prec aux_min_distance_a;
  real_prec aux_min_distance_b;
  

  ULONG pair_counting=0;

  So.message_screen("Computing minimum separations");
  // no paralelizar si hay un min o max() dentro del loop
  for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      
      real_prec x_coord=this->Halo[i].coord1;
      real_prec y_coord=this->Halo[i].coord2;
      real_prec z_coord=this->Halo[i].coord3;
      
#ifdef _USE_VMAX_AS_OBSERVABLE_
      real_prec lmass=log10(this->Halo[i].vmax);
#elif defined _USE_MASS_AS_OBSERVABLE_
      real_prec lmass=log10(this->Halo[i].mass)+log10(this->params._MASS_units());
#endif

      int im=get_bin(lmass,lm_min,N_mass_bins,this->logdeltaM,true);

      ULONG ID=this->Halo[i].GridID;
      
      // Allocate the min_distances to i particle form neighbouring cells
      vector<real_prec> min_distances_v(N_Neigh_cells,0);
      
      for(int j=0;j<N_Neigh_cells; j++) // loop sobre las celdas cercanas incluida la celda donde esta la particula i
	{
	  ULONG ID_NEIGH = nearest_cells_to_cell[ID].close_cell[j];
	  
	  real_prec factor_bc_x=nearest_cells_to_cell[ID].bc_x[j]*this->box.Lbox;
	  real_prec factor_bc_y=nearest_cells_to_cell[ID].bc_y[j]*this->box.Lbox;
	  real_prec factor_bc_z=nearest_cells_to_cell[ID].bc_z[j]*this->box.Lbox;
	  
	  aux_min_distance_b=1e5;
	  real_prec min_distance_b;
	  
	  for(int k=0;k< cell_info_tr[ID_NEIGH].posx_p.size(); ++k) //loop sobre las particulas en las celdas cercanas
	    {
	      
              pair_counting++;
	      ULONG index_gal=cell_info_tr[ID_NEIGH].gal_index[k];

              real_prec lmass_k=log10(this->Halo[index_gal].mass)+log10(this->params._MASS_units());
              int im_k=get_bin(lmass_k,lm_min,N_mass_bins,this->logdeltaM,true);

              int index_mass=index_2d(im,im_k,N_mass_bins);

	      if(index_gal!=i) //this simply avoids to get a 0 in the distance
		{
		  
		  real_prec distx = x_coord - (cell_info_tr[ID_NEIGH].posx_p[k]+factor_bc_x);
		  real_prec disty = y_coord - (cell_info_tr[ID_NEIGH].posy_p[k]+factor_bc_y);
		  real_prec distz = z_coord - (cell_info_tr[ID_NEIGH].posz_p[k]+factor_bc_z);
		  real_prec dist  = sqrt(pow(distx,2)+pow(disty,2)+pow(distz,2));
		  
		  min_distance_b=min(aux_min_distance_b,dist);
		  aux_min_distance_b=min_distance_b;
		  
		  int da_bin = get_bin(dist,0.,N_bins_dist,max_dist/static_cast<real_prec>(N_bins_dist),false);
                  ULONG sindex = index_2d(da_bin,index_mass,N_mass_bins);
                  if(dist<max_dist)
                    distance_dist[sindex]++;
		}
	    }// at this point min_distance_b exits as the minimum within *one* neighbour cell
	  
	  min_distances_v[j]=aux_min_distance_b;	  
	  
	}// at this point min_distance_a exits as the minimum among all neighbour cells
      
      aux_min_distance_a=1e5;
      real_prec min_distance_a;
      for(int j=0;j<N_Neigh_cells; j++) 
	{
	  if(min_distances_v[j]>0)
	    {
	      min_distance_a=min(min_distances_v[j],aux_min_distance_a);
	      aux_min_distance_a=min_distance_a;
	    }

	}
      
      int d_bin = get_bin(min_distance_a,0.,N_bins_dist,max_dist/static_cast<real_prec>(N_bins_dist),false);
      if(min_distance_a<max_dist)
        min_distance_dist[d_bin][im]++;

      //if(i==1)cout<<min_distance_a<<endl;
       
     }
   So.DONE();
   

   So.message_screen("Total number of pairs", pair_counting);
   
   ofstream sal;

   sal.open("min_dist_"+data+".txt");
   for(int i=0;i<min_distance_dist.size();++i)
     {
       sal<<(i+0.5)*max_dist/static_cast<real_prec>(N_bins_dist)<<" ";
       for(int j=0;j< N_mass_bins; ++j)
	 sal<<min_distance_dist[i][j]<<" ";
       sal<<endl;
     }
   sal.close();


   So.message_screen("Writting P(r|M,M') ");
   for(int i=0;i<min_distance_dist.size();++i)
     {
       string ff="dist_"+data+"_separation_bin"+to_string(i)+".txt";
       So.message_screen("in file ", ff);
       sal.open(ff.c_str());
       for(int j=0;j< N_mass_bins; ++j)
         for(int k=0;k< N_mass_bins; ++k)
          {
            int index=index_2d(j,k,N_mass_bins);
            ULONG sindex = index_2d(i,index,N_mass_bins);
            sal<<lm_min+(j+0.5)*logdeltaM<<" "<<lm_min+(k+0.5)*logdeltaM<<"  "<<static_cast<real_prec>(distance_dist[sindex])/static_cast<real_prec>(pair_counting)<<endl;
           }
       sal.close();
       }
   So.DONE();
 }



// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// ***********************************************************************************************************  //
// for each bin in prop, we will allocate the list of masses (M1, M2) in pairs with separatIons d1-2 between d_min and 1.5Mpc
// WE do this from the reference as well as from the final mock with assigned masses
// the index im in this->masses_in_cells_min_sep_ref[ih].M1[im] denotes the mass1 of the im-th pair in the theta bin ih separated a distance [dmin, Dmax]
// frin the other partice having mass this->masses_in_cells_min_sep_ref[ih].M2[im]
//

void Catalog::get_masses_of_pairs_in_min_separation_bin_in_theta_bin(real_prec min_sep, vector<s_mass_members> & dm_properties_bins)
{
    int NTHREADS = omp_get_max_threads();
    omp_set_num_threads(NTHREADS);

    So.message_screen("**Getting list of pair of masses in bins of theta with separations");
    So.message_screen("**between the the minimum of the reference and", EXCLUSION_SCALE," Mpc/h:");
  //  So.message_screen("Using ", NTHREADS, "threads");

//#pragma omp parallel for
    for(int ih=0;ih< dm_properties_bins.size();++ih)
       {
        int nobjs=  dm_properties_bins[ih].masses_bin_properties.size();
        if(nobjs>=2)
          {
            for(int j=0;j< nobjs;++j)
              {
//                ULONG Id1= dm_properties_bins[ih].GridID_bin_properties[j];
                for(int k=j+1;k< nobjs;++k)
                  {
  //                  ULONG Id2= dm_properties_bins[ih].GridID_bin_properties[k];
                    real_prec eta_dist=  pow(dm_properties_bins[ih].x_coord_bin_properties[j]-dm_properties_bins[ih].x_coord_bin_properties[k],2);
                    eta_dist+=pow(dm_properties_bins[ih].y_coord_bin_properties[j]-dm_properties_bins[ih].y_coord_bin_properties[k],2);
                    eta_dist+=pow(dm_properties_bins[ih].z_coord_bin_properties[j]-dm_properties_bins[ih].z_coord_bin_properties[k],2);
                    eta_dist=sqrt(eta_dist);
                    if(eta_dist>= min_sep && eta_dist<=EXCLUSION_SCALE)
                     {
                       this->masses_in_cells_min_sep[ih].M1.push_back(dm_properties_bins[ih].masses_bin_properties[j]);
                       this->masses_in_cells_min_sep[ih].M2.push_back(dm_properties_bins[ih].masses_bin_properties[k]);
                     }
                  }
              }
          }
      }

    So.DONE();
}

