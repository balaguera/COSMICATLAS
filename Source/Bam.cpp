//##################################################################################
//##################################################################################
/** @file Bam.cpp
 *
 *  @brief Generation of mock catalogs of DM tracers
 *  based on the BAM method.
 *  @author: Andrés Balaguera-Antolínez, Francisco-Shu Kitaura, IAC, 2017-2019
 */
# include "../Headers/Bam.h"

#ifdef _USE_PYTHON_
void py_plot(vector<real_prec>& field)
{
  return;
}
#endif

//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
void Bam::warnings(){

  this->So.enter(__PRETTY_FUNCTION__);
#if !defined (_USE_DM_IN_BAM_)
      So.message_warning("Warning: _USE_DM_IN_BAM_ is an undefined pre-proc directive");
      So.message_warning("If this code is meant to run with DM (default option), please define it in def.h.");
      So.message_warning("Otherwise, type c to continue");
      string cont;
      cin>>cont;
      if (cont!="c" || cont!="C")
        So.message_warning("COSMICATLAS stops here.");
      exit(0);
#endif
#if !defined (_DO_BAM_CALIBRATION_) && defined (_MODIFY_LIMITS_)
    So.message_warning("Warning: _MODIFY_LIMITS_ is defined under directive _GET_BAM_REALIZATIONS_.");
    So.message_warning("This is a potential source of  bug/collapse of the code. Please check it.");
    exit(0);
#endif
#if defined (_USE_NABLA2DELTA_) && defined (_USE_INVARANT_SHEAR_VFIELD_I_)
    So.message_warning("Warning: _USE_NABLA2DELTA_ and _USE_INVARANT_SHEAR_VFIELD_I_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_S2DELTA_) && defined (_USE_INVARANT_SHEAR_VFIELD_II_)
    So.message_warning("Warning: _USE_S2DELTA_ and _USE_INVARANT_SHEAR_VFIELD_II_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_S3_) && defined (_USE_INVARANT_SHEAR_VFIELD_III_)
    So.message_warning("Warning: _USE_S3_ and _USE_INVARANT_SHEAR_VFIELD_III_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_S2_) && defined (_USE_TIDAL_ANISOTROPY_)
    So.message_warning("Warning: _USE_S2_ and _USE_TIDAL_ANISOTROPY_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_S2_) && defined (_USE_INVARIANT_TIDAL_FIELD_IV_)
    So.message_warning("Warning: _USE_S2_ and _USE_INVARIANT_TIDAL_FIELD_IV_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_TIDAL_ANISOTROPY_) && defined (_USE_INVARIANT_TIDAL_FIELD_IV_)
    So.message_warning("Warning: _USE_TIDAL_ANISOTROPY_ and _USE_INVARIANT_TIDAL_FIELD_IV_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_TIDAL_ANISOTROPY_) && defined (_USE_ELLIPTICITY_)
    So.message_warning("Warning: _USE_TIDAL_ANISOTROPY_ and _USE_ELLIPTICITY_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_TIDAL_ANISOTROPY_) && defined (_USE_PROLATNESS_)
    So.message_warning("Warning: _USE_TIDAL_ANISOTROPY_ and _USE_PROLATNESS_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_PROLATNESS_) && defined (_USE_ELLIPTICITY_)
    So.message_warning("Warning: _USE_PROLATNESS_ and _USE_ELLIPTICITY_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_INVARIANT_TIDAL_FIELD_IV_) && defined (_USE_ELLIPTICITY_)
    So.message_warning("Warning: _USE_TIDAL_ANISOTROPY_ and _USE_ELLIPTICITY_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif

#if defined (_USE_INVARIANT_TIDAL_FIELD_IV_) && defined (_USE_PROLATNESS_)
    So.message_warning("Warning: _USE_TIDAL_ANISOTROPY_ and _USE_PROLATNESS_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_S2_) && defined (_USE_ELLIPTICITY_)
    So.message_warning("Warning: _USE_TIDAL_ANISOTROPY_ and _USE_ELLIPTICITY_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_S2_) && defined (_USE_PROLATNESS_)
    So.message_warning("Warning: _USE_TIDAL_ANISOTROPY_ and _USE_PROLATNESS_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_DELTA3_) && defined (_USE_INVARIANT_TIDAL_FIELD_III_)
    So.message_warning("Warning: _USE_DELTA3_ and _USE_INVARIANT_TIDAL_FIELD_III_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_DELTA2_) && defined (_USE_INVARIANT_TIDAL_FIELD_II_)
    So.message_warning("Warning: _USE_DELTA2_ and _USE_INVARIANT_TIDAL_FIELD_II_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
}
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
void Bam::set_params_bam()
{

  this->So.enter(__PRETTY_FUNCTION__);
#ifdef _FULL_VERBOSE_
  So.message_screen("Loading parameters for BAM");
#endif
  this->NX=this->params._NX();  //!< Parameter NX  //!<
  this->NY=this->params._NY();
  this->NY_MASS=this->params._NY_MASS();
  this->NY_SAT_FRAC=this->params._NY_SAT_FRAC();
  this->ndel_data = this->params._ndel_data();
  this->Output_directory=this->params._Output_directory();
  this->Input_Directory_Y=this->params._Input_Directory_Y();
  this->Name_Catalog_Y=this->params._Name_Catalog_Y();
  this->Name_Catalog_Y_HR=this->params._Name_Catalog_Y_HR();
  //  this->Name_Catalog_Y_MWEIGHTED=this->params._Name_Catalog_Y_MWEIGHTED();
  this->Input_Directory_X=this->params._Input_Directory_X();
  this->Input_Directory_X_REF=this->params._Input_Directory_X_REF();
  this->Input_Directory_X_NEW=this->params._Input_Directory_X_NEW();
  this->XNAME=this->params._XNAME();
  this->Name_Catalog_X=this->params._Name_Catalog_X();
  this->Name_VelFieldx_X=this->params._Name_VelFieldx_X();
  this->Name_VelFieldy_X=this->params._Name_VelFieldy_X();
  this->Name_VelFieldz_X=this->params._Name_VelFieldz_X();
  this->Name_Catalog_X_REF_PDF=this->params._Name_Catalog_X_REF_PDF();
  this->Name_Catalog_X_NEW=this->params._Name_Catalog_X_NEW();
  this->Name_Property_X=this->params._Name_Property_X();
  this->Name_redshift_mask=this->params._Name_redshift_mask();
  this->Name_binary_mask=this->params._Name_binary_mask();
  this->YNAME=this->params._YNAME();
  this->Name_Property_Y=this->params._Name_Property_Y();
  this->iMAS_X=this->params._iMAS_X();
  this->iMAS_X_REF_PDF=this->params._iMAS_X_REF_PDF();
  this->iMAS_X_NEW=this->params._iMAS_X_NEW();
  this->iMAS_Y=this->params._iMAS_Y();
  this->delta_Y_max=this->params._delta_Y_max();
  this->delta_Y_min=this->params._delta_Y_min();
  this->delta_X_max=this->params._delta_X_max();
  this->delta_X_min=this->params._delta_X_min();
  this->ldelta_Y_max=this->params._ldelta_Y_max();
  this->ldelta_Y_min=this->params._ldelta_Y_min();
  this->ldelta_X_max=this->params._ldelta_X_max();
  this->ldelta_X_min=this->params._ldelta_X_min();
  this->Quantity=this->params._Quantity();
  this->NMASSbins=this->params._NMASSbins();
  this->NMASSbins_mf=this->params._NMASSbins_mf();
  this->redshift=this->params._redshift();
  this->smscale=this->params._smscale();
  this->realization=this->params._realization();
  this->Apply_Rankordering=this->params._Apply_Rankordering(); // to be deprecated
  this->Apply_Rankordering_ab_initio=this->params._Apply_Rankordering_ab_initio();
  this->Nft=this->params._Nft();
  this->Nft_low=this->params._Nft_low();
  this->Nft_HR=this->params._Nft_HR();
  this->write_files_for_histograms=this->params._write_files_for_histograms();
  this->Redefine_limits=this->params._Redefine_limits();
  this->Convert_Density_to_Delta_X=this->params._Convert_Density_to_Delta_X();
  this->Convert_Density_to_Delta_Y=this->params._Convert_Density_to_Delta_Y();
  this->lambdath=this->params._lambdath();
  this->Write_PDF_number_counts=this->params._Write_PDF_number_counts();
  this->Scale_X=this->params._Scale_X();
  this->Scale_Y=this->params._Scale_Y();
  this->n_sknot_massbin=this->params._n_sknot_massbin();
  this->n_vknot_massbin=this->params._n_vknot_massbin();
  this->Lbox=this->params._Lbox();
  this->Lbox_low=this->params._Lbox_low();
  this->N_iterations_Kernel=this->params._N_iterations_Kernel();
  this->iteration_ini=this->params._iteration_ini();
  this->N_dm_initial = this->params._N_dm_initial();
  this->N_dm_realizations=this->params._N_dm_realizations();
  this->N_iterations_dm=this->params._N_iterations_dm();
  this->output_at_iteration=this->params._output_at_iteration();
  this->n_cwt=this->params._n_cwt();
  this->n_cwv=this->params._n_cwv();
  //this->aux_kernel = false; // aparentemente no se usa
  this->used_once = true;
  // This is a value used when no itaration is computed, but direact
#ifdef _GET_BAM_REALIZATIONS_
  ifstream nxo;
  nxo.open(file_one_cell);
#ifdef _FULL_VERBOSE_
  if(nxo.is_open())
    So.message_screen("Reading maximum number of reference tracers in cells from file ",file_one_cell);
  else
  {
      So.message_screen("File with max number of tracers in one cell does not exist",file_one_cell);
      So.message_screen("BAM stops here");
      exit(0);
  }
#endif
  int iaux; int iaux2;
  nxo>>iaux>>iaux2;
  this->nmax_y_onecell=iaux;
  this->mean_number_y_onecell=iaux2;
  nxo.close();

#ifdef _FULL_VERBOSE_
  if(iaux==0)
    So.message_warning("Maximum number of cells = 0. End");
  else
    So.message_screen("maximum number of reference tracers in cells =",iaux);
  So.message_screen("mean number of reference tracers in cells =",iaux2);
#endif

#endif
  this->seed=this->params._seed();
  this->dm_already_done=false;
  this->Distance_fraction = this->params._Distance_fraction();
  this->Nft_random_collapse = this->params._Nft_random_collapse();
  // Feed the structure for cosmological parameters
  // This is also done in patchy, so verify that these lines are also in its init par member
  this->s_cosmo_pars.cosmological_redshift=this->params._redshift();
  this->s_cosmo_pars.Hubble=this->params._Hubble();
  this->s_cosmo_pars.hubble=this->params._hubble();
  this->s_cosmo_pars.Om_matter =this->params._om_matter();
  this->s_cosmo_pars.Om_cdm =this->params._om_cdm();
  this->s_cosmo_pars.Om_baryons =this->params._om_baryons();
  this->s_cosmo_pars.Om_radiation =this->params._om_radiation();
  this->s_cosmo_pars.Om_vac =this->params._om_vac();
  this->s_cosmo_pars.Om_k =this->params._om_k();
  this->s_cosmo_pars.n_s =this->params._spectral_index();
  this->s_cosmo_pars.w_eos =this->params._w_eos();
  this->s_cosmo_pars.N_eff =this->params._N_eff();
  this->s_cosmo_pars.sigma8 =this->params._sigma8();
  this->s_cosmo_pars.f_baryon =this->params._om_baryons()/this->params._om_matter();
  this->s_cosmo_pars.use_wiggles =this->params._use_wiggles();
  this->s_cosmo_pars.RR =this->params._RR();
#ifdef _BIN_ACCUMULATE_
  this->bin_accumulate_borders = true;
#else
  this->bin_accumulate_borders = false;
#endif
  So.DONE();
  // Determine the number of properties used to characterize the halo bias
  this->N_bias_properties=1;  // This accounts for the dark amtter density.
  vector<string> bias_properties;
  bias_properties.push_back("Local density");
#ifdef _USE_MASS_KNOTS_
  this->N_bias_properties++;
  bias_properties.push_back("Knot-Mass");
#endif
#ifdef _USE_CWC_
  this->N_bias_properties++;
  bias_properties.push_back("T-WEB");
#endif

#ifdef _USE_VEL_KNOTS_V_
  this->N_bias_properties++;
  bias_properties.push_back("V-dispersion in knots");
#endif
#ifdef _USE_CWC_V_
  this->N_bias_properties++;
  bias_properties.push_back("V-WEB");
#endif
#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_DELTA2_)
  this->N_bias_properties++;
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
  bias_properties.push_back("I-2");
#else
  bias_properties.push_back("ð²");
#endif
#endif
#if defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_DELTA3_)
  this->N_bias_properties++;
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
  bias_properties.push_back("I-3");
#else
  bias_properties.push_back("ð³");
#endif
#endif
#if defined (_USE_INVARIANT_TIDAL_FIELD_IV_) || defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_ELLIPTICITY_)|| defined (_USE_PROLATNESS_)|| defined (_USE_S2_)
  this->N_bias_properties++;
#ifdef _USE_INVARIANT_TIDAL_FIELD_IV_
  bias_properties.push_back("I-4");
#elif defined (_USE_TIDAL_ANISOTROPY_)
  bias_properties.push_back("Tidal anisotropy");
#elif defined (_USE_ELLIPTICITY_)
  bias_properties.push_back("Ellipticity");
#elif defined (_PROLATNESS_)
    bias_properties.push_back("Prolatness");
#else
  bias_properties.push_back("S²");
#endif
#endif
#if defined (_USE_INVARIANT_SHEAR_VFIELD_I_) || defined(_USE_NABLA2DELTA_)
  this->N_bias_properties++;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
  bias_properties.push_back("IV-1");
#else
  bias_properties.push_back("Nabla² ð");
#endif
#endif
#if defined (_USE_INVARIANT_SHEAR_VFIELD_II_) || defined (_USE_S2DELTA_)
  this->N_bias_properties++;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
  bias_properties.push_back("IV-2");
#else
  bias_properties.push_back("s²ð");
#endif
#endif
#if defined (_USE_INVARIANT_SHEAR_VFIELD_III_) || defined (_USE_S3_)
  this->N_bias_properties++;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
  bias_properties.push_back("IV-3");
#else
  bias_properties.push_back("s²");
#endif
#endif
#ifdef _FULL_VERBOSE_
#ifndef _ONLY_PATCHY_
#ifdef _DO_BAM_CALIBRATION_
  So.message_screen("****************************************************************************");
  So.message_screen("BAM is using", this->N_bias_properties,"properties to characterise the bias:");
  for(int i=0;i<bias_properties.size();++i)cout<<YELLOW<<bias_properties[i]<<RESET<<endl;
  So.message_screen("****************************************************************************");
  cout<<endl;
#endif
#endif
#endif
  this->step_mass_assignemt=0; //initialize it
#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_
  this->Prop_threshold=PROP_THRESHOLD; // to be deprecated
#ifdef _USE_MULTISCALE_LEVEL_1_
  this->Prop_threshold_multi_scale_1=(PROP_THRESHOLD_MULTI_SCALE_1);
#endif
#ifdef _USE_MULTISCALE_LEVEL_2_
  this->Prop_threshold_multi_scale_2=(PROP_THRESHOLD_MULTI_SCALE_2);
#endif
#ifdef _USE_MULTISCALE_LEVEL_3_
  this->Prop_threshold_multi_scale_3=(PROP_THRESHOLD_MULTI_SCALE_3);
#endif
#ifdef _USE_MULTISCALE_LEVEL_4_
  this->Prop_threshold_multi_scale_4=(PROP_THRESHOLD_MULTI_SCALE_4);
#endif
#endif
}

//##################################################################################
//##################################################################################
void Bam::get_cosmo()
{
  this->So.enter(__PRETTY_FUNCTION__);

#ifdef _FULL_VERBOSE_
  So.message_screen("Getting cosmological derived-parameters");
#endif

  // Here we fill the structure s_cosmo_info with different cosmological quantities evaluated at the input redshift
  this->s_cosmo_info.scale_factor=1./(this->redshift+1.);
  this->s_cosmo_info.critical_density=this->Cosmo.critical_density(redshift, (void *)&this->s_cosmo_pars);
  this->s_cosmo_info.density_contrast_top_hat=Cosmo.density_contrast_top_hat(redshift, (void *)&this->s_cosmo_pars);
  this->s_cosmo_info.Hubble_parameter=this->Cosmo.Hubble_function(redshift, (void *)&this->s_cosmo_pars);
  this->s_cosmo_info.comoving_distance=this->Cosmo.comoving_distance(redshift, (void *)&this->s_cosmo_pars);
  this->s_cosmo_info.comoving_angular_diameter_distance=Cosmo.comoving_angular_diameter_distance(redshift, (void *)&this->s_cosmo_pars);
  this->s_cosmo_info.mean_matter_density=Cosmo.mean_matter_density(redshift, (void *)&this->s_cosmo_pars);
  this->s_cosmo_info.age_universe=Cosmo.age_universe(redshift, (void *)&this->s_cosmo_pars);
  this->s_cosmo_info.comoving_sound_horizon=Cosmo.comoving_sound_horizon(redshift, (void *)&this->s_cosmo_pars);
  this->s_cosmo_info.Delta_Vir = Cosmo.density_contrast_top_hat(redshift, (void *)&this->s_cosmo_pars);

#ifdef _USE_PATCHY_  // This is repeated here as it is in patchy
  this->s_cosmo_info.growth_factor=this->Cosmo.growth_factor(redshift, (void *)&this->s_cosmo_pars)/this->Cosmo.growth_factor(0.0,(void *)&this->s_cosmo_pars);
  this->s_cosmo_info.growth_index=this->Cosmo.growth_index(redshift, (void *)&this->s_cosmo_pars);
  this->s_cosmo_info.growth_index2=this->Cosmo.growth_index2(redshift, (void *)&this->s_cosmo_pars);
  this->s_cosmo_info.halo_dynamical_time=Cosmo.halo_dynamical_time(redshift, (void *)&this->s_cosmo_pars);
  this->s_cosmo_info.omega_matter=Cosmo.omega_matter(redshift, (void *)&this->s_cosmo_pars);
  this->s_cosmo_info.Distance_Modulus=Cosmo.Distance_Modulus(redshift, (void *)&this->s_cosmo_pars);
  this->s_cosmo_pars.growth_factor=this->s_cosmo_info.growth_factor;

  real_prec fD2=static_cast<real_prec>(pow(this->s_cosmo_info.omega_matter,-1./143.));
  this->s_cosmo_info.D2=static_cast<real_prec>(-(3./7.)*pow(this->s_cosmo_info.growth_factor,2)*fD2);



  PowerSpectrum Pow;
  this->s_cosmo_pars.growth_factor=this->s_cosmo_info.growth_factor;
  this->s_cosmo_pars.pk_normalization=Pow.normalization((void *)&this->s_cosmo_pars);
#endif

#ifdef _USE_PATCHY_  // This is repeated here as it is in patchy
#ifdef _FULL_VERBOSE_
  this->So.write_cosmo_parameters((void *)&this->s_cosmo_pars, (void *)&this->s_cosmo_info);
#endif
#endif

  So.DONE();


  this->cwclass.s_cosmo_info=this->s_cosmo_info;

}


//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
void Bam::show_params()
{
  this->So.enter(__PRETTY_FUNCTION__);

  cout<<BOLDYELLOW;
  cout<<"***************************************************************************"<<endl;
  cout<<"***************************************************************************"<<endl;
  cout<<"BAM                                                                       *"<<endl;
  cout<<"Bias Assignment Method for galaxy/halo mock catalogs                      *"<<endl;
  cout<<"Input values of parameters in parameter file                              *"<<endl;
  cout<<"***************************************************************************"<<endl;
  cout<<"***************************************************************************"<<RESET<<endl;

  cout<<"Redshift                   = "<<this->params._redshift()<<endl;
  cout<<"Lbox                       = "<<this->params._Lbox()<<endl;
  cout<<"NX                         = "<<this->NX<<endl;
  cout<<"NY                         = "<<this->NY<<endl;
  cout<<"Output_directory           = "<<this->params._Output_directory()<<endl;
  cout<<"Input_Directory_X          = "<<this->params._Input_Directory_X()<<endl;
  cout<<"Input_Directory_X_REF      = "<<this->params._Input_Directory_X_REF()<<endl;
  cout<<"XNAME                      = "<<this->params._XNAME()<<endl;
  cout<<"Name_Catalog_X             = "<<this->params._Name_Catalog_X()<<endl;
  cout<<"Name_redshift_mask         = "<<this->params._Name_redshift_mask()<<endl;
  cout<<"Name_Catalog_X_REF_PDF     = "<<this->params._Name_Catalog_X_REF_PDF()<<endl;
  cout<<"Name_Catalog_X_NEW         = "<<this->params._Name_Catalog_X_NEW()<<endl;
  cout<<"Name_Property_X            = "<<this->params._Name_Property_X()<<endl;

  cout<<"iMAS_X                     = "<<this->params._iMAS_X()<<endl;
  cout<<"iMAS_X_REF_PDF             = "<<this->params._iMAS_X_REF_PDF()<<endl;
  cout<<"iMAS_X_NEW                 = "<<this->params._iMAS_X_NEW()<<endl;
  cout<<"Input_Directory_Y          = "<<this->params._Input_Directory_Y()<<endl;
  cout<<"YNAME                      = "<<this->params._YNAME()<<endl;
  cout<<"Name_Catalog_Y             = "<<this->params._Name_Catalog_Y()<<endl;
  cout<<"Name_Catalog_Y_MWEIGHTED   = "<<this->params._Name_Catalog_Y_MWEIGHTED()<<endl;
  cout<<"Name_Catalog_Y_HR          = "<<this->params._Name_Catalog_Y_HR()<<endl;
  cout<<"Name_Property_Y            = "<<this->params._Name_Property_Y()<<endl;
  cout<<"iMAS_Y                     = "<<this->params._iMAS_Y()<<endl;
  cout<<"delta_Y_max                = "<<this->params._delta_Y_max()<<endl;
  cout<<"delta_Y_min                = "<<this->params._delta_Y_min()<<endl;
  cout<<"delta_X_max                = "<<this->params._delta_X_max()<<endl;
  cout<<"delta_X_min                = "<<this->params._delta_X_min()<<endl;
  cout<<"ldelta_Y_max               = "<<this->params._ldelta_Y_max()<<endl;
  cout<<"ldelta_Y_min               = "<<this->params._ldelta_Y_min()<<endl;
  cout<<"ldelta_X_max               = "<<this->params._ldelta_X_max()<<endl;
  cout<<"ldelta_X_min               = "<<this->params._ldelta_X_min()<<endl;
  cout<<"Quantity                   = "<<this->params._Quantity()<<endl;
  cout<<"NMASSbins                  = "<<this->params._NMASSbins()<<endl;
  cout<<"smscale                    = "<<this->params._smscale()<<endl;
  cout<<"realization                = "<<this->params._realization()<<endl;
  cout<<"Apply_Rankordering         = "<<this->params._Apply_Rankordering()<<endl;
  cout<<"Nft                        = "<<this->params._Nft()<<endl;
  cout<<"ndel_data                  = "<<this->params._ndel_data()<<endl;
  cout<<"write_files_for_histograms = "<<this->params._write_files_for_histograms()<<endl;
  cout<<"Redefine_limits            = "<<this->params._Redefine_limits()<<endl;
  cout<<"Convert_Density_to_Delta_X = "<<this->params._Convert_Density_to_Delta_X()<<endl;
  cout<<"Convert_Density_to_Delta_Y = "<<this->params._Convert_Density_to_Delta_Y()<<endl;
  cout<<"lambdath                   = "<<this->params._lambdath()<<endl;
  cout<<"Write_PDF_number_counts    = "<<this->params._Write_PDF_number_counts()<<endl;
  cout<<"Scale_X                    = "<<this->params._Scale_X()<<endl;
  cout<<"Scale_Y                    = "<<this->params._Scale_Y()<<endl;
  cout<<"n_sknot_massbin            = "<<this->params._n_sknot_massbin()<<endl;
  cout<<"n_vknot_massbin            = "<<this->params._n_vknot_massbin()<<endl;
  cout<<"N_iterations_Kernel        = "<<this->params._N_iterations_Kernel()<<endl;
  cout<<"N_iterations_dm            = "<<this->params._N_iterations_dm()<<endl;
  cout<<"N_dm_realizations          = "<<this->params._N_dm_realizations()<<endl;
  cout<<"n_cwt                      = "<<this->params._n_cwt()<<endl;
  cout<<"n_cwv                      = "<<this->params._n_cwv()<<endl;

}

//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
void Bam::show_preproc()
{
  cout<<COLOR_DEFINED<<"DEFINED"<<RESET<<endl;
  cout<<COLOR_UNDEFINED<<"UNDEFINED"<<RESET<<endl;

  this->So.message_screen("PRECISION");
#ifdef SINGLE_PREC
  cout<<COLOR_DEFINED<<"SINGLE_PREC"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"SINGLE_PREC"<<RESET<<RESET<<endl;
#endif
#ifdef DOUBLE_PREC
  cout<<COLOR_DEFINED<<"DOUBLE_PREC"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"DOUBLE_PREC"<<RESET<<RESET<<endl;
#endif

  cout<<endl;
  this->So.message_screen("OMP");
#ifdef _USE_OMP_
  cout<<COLOR_DEFINED<<"_USE_OMP_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_OMP_"<<RESET<<RESET<<endl;
#endif

  cout<<endl;
  this->So.message_screen("VERBOSE");
#ifdef _USE_COLORS_
  cout<<COLOR_DEFINED<<"_USE_COLORS_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_COLORS_"<<RESET<<RESET<<endl;
#endif
#ifdef _FULL_VERBOSE_
  cout<<COLOR_DEFINED<<"_FULL_VERBOSE_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_FULL_VERBOSE_"<<RESET<<RESET<<endl;
#endif

  cout<<endl;
  this->So.message_screen("IO FORMATS");
#ifdef _WRITE_BINARY_BAM_FORMAT_
  cout<<COLOR_DEFINED<<"_WRITE_BINARY_BAM_FORMAT_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_WRITE_BINARY_BAM_FORMAT_"<<RESET<<RESET<<endl;
#endif
#ifdef _READ_BINARY_BAM_FORMAT_
  cout<<COLOR_DEFINED<<"_READ_BINARY_BAM_FORMAT_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_READ_BINARY_BAM_FORMAT_"<<RESET<<RESET<<endl;
#endif
#ifdef _OUTPUT_WITH_HEADERS_
  cout<<COLOR_DEFINED<<"_OUTPUT_WITH_HEADERS_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_OUTPUT_WITH_HEADERS_"<<RESET<<RESET<<endl;
#endif



  cout<<endl;
  this->So.message_screen("COSMOLOGY");
#ifdef _USE_COSMO_PARS_
  cout<<COLOR_DEFINED<<"_USE_COSMO_PARS_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_COSMO_PARS_"<<RESET<<RESET<<endl;
#endif
#ifdef _USE_UNITSIM_COSMOLOGY_
  cout<<COLOR_DEFINED<<"_USE_UNITSIM_COSMOLOGY_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_UNITSIM_COSMOLOGY_"<<RESET<<RESET<<endl;
#endif
#ifdef _USE_PLANCK_COSMOLOGY_
  cout<<COLOR_DEFINED<<"_USE_PLANCK_COSMOLOGY_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_PLANCK_COSMOLOGY_"<<RESET<<RESET<<endl;
#endif


  cout<<endl;
  this->So.message_screen("BAM");
#ifdef BIAS_MODE
  cout<<COLOR_DEFINED<<"BIAS_MODE"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"BIAS_MODE"<<RESET<<endl;
#endif
#ifdef MOCK_MODE
  cout<<COLOR_DEFINED<<"MOCK_MODE"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"MOCK_MODE"<<RESET<<endl;
#endif

#ifdef _GET_BAM_REALIZATIONS_
  cout<<COLOR_DEFINED<<"_GET_BAM_REALIZATIONS_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_GET_BAM_REALIZATIONS_"<<RESET<<endl;
#endif

#ifdef _DO_BAM_CALIBRATION_
  cout<<COLOR_DEFINED<<"_DO_BAM_CALIBRATION_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_DO_BAM_CALIBRATION_"<<RESET<<endl;
#endif


#ifdef _USE_VELOCITIES_
  cout<<COLOR_DEFINED<<"_USE_VELOCITIES_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_VELOCITIES_"<<RESET<<endl;
#endif




cout<<endl;
this->So.message_screen("PATCHY");
#ifdef _USE_PATCHY_
  cout<<COLOR_DEFINED<<"_USE_PATCHY_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_PATCHY_"<<RESET<<endl;
#endif
#ifdef _ONLY_PATCHY_
  cout<<COLOR_DEFINED<<"_ONLY_PATCHY_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_ONLY_PATCHY_"<<RESET<<endl;
#endif
#ifdef _POWER_BOOST_ALPT_
  cout<<COLOR_DEFINED<<"_POWER_BOSST_ALPT_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_POWER_BOOST_ALPT_"<<RESET<<endl;
#endif
#ifdef _GET_VELOCITIES_PATCHY_
  cout<<COLOR_DEFINED<<"_GET_VELOCITIES_PATCHY_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_GET_VELOCITIES_PATCHY_"<<RESET<<endl;
#endif
#ifdef _GET_DENSITY_FIELDS_
  cout<<COLOR_DEFINED<<"_GET_DENSITY_FIELDS_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_GET_DENSITY_FIELDS_"<<RESET<<endl;
#endif


cout<<endl;
this->So.message_screen("CATALOG");



#ifdef _GET_BAM_CAT_
  cout<<COLOR_DEFINED<<"_GET_BAM_CAT_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_GET_BAM_CAT_"<<RESET<<endl;
#endif
#ifdef _READ_REF_CATALOG_
  cout<<COLOR_DEFINED<<"_READ_REF_CATALOG_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_READ_REF_CATALOG_"<<RESET<<endl;
#endif

#ifdef _USE_MASS_TRACERS_
  cout<<COLOR_DEFINED<<"_USE_MASS_TRACERS_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_MASS_TRACERS_"<<RESET<<endl;
#endif

#ifdef _USE_VELOCITIES_TRACERS_
  cout<<COLOR_DEFINED<<"_USE_VELOCITIES_TRACERS_"<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_VELOCITIES_TRACERS_"<<RESET<<endl;
#endif

#ifdef _USE_VMAX_AS_OBSERVABLE_
  cout<<COLOR_DEFINED<<"_USE_VMAX_AS_OBSERVABLE_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_VMAX_AS_OBSERVABLE_"<<RESET<<endl;
#endif


#ifdef _USE_MASS_AS_OBSERVABLE_
  cout<<COLOR_DEFINED<<"_USE_MASS_AS_OBSERVABLE_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_MASS_AS_OBSERVABLE_"<<RESET<<endl;
#endif



#ifdef _SET_CAT_WITH_MASS_CUT_
  cout<<COLOR_DEFINED<<"_SET_CAT_WITH_MASS_CUT_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_SET_CAT_WITH_MASS_CUT_"<<RESET<<endl;
#endif


#ifdef _WRITE_BAM_CATALOGS_
  cout<<COLOR_DEFINED<<"_WRITE_BAM_CATALOGS_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_WRITE_BAM_CATALOGS_"<<RESET<<endl;
#endif


#ifdef _WRITE_BINARY_BAM_FORMAT_
  cout<<COLOR_DEFINED<<"_WRITE_BAM_CATALOGS_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_WRITE_BINARY_BAM_FORMAT_"<<RESET<<endl;
#endif




cout<<endl;
this->So.message_screen("CALIBRATION");

#ifdef _DO_BAM_CALIBRATION_
  cout<<COLOR_DEFINED<<"_DO_BAM_CALIBRATION_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_DO_BAM_CALIBRATION_"<<RESET<<endl;
#endif
#ifdef _GET_BAM_REALIZATIONS_
  cout<<COLOR_DEFINED<<"_GET_BAM_REALIZATIONS_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_GET_BAM_REALIZATIONS_"<<RESET<<endl;
#endif
#ifdef _SMOOTHED_KERNEL_
  cout<<COLOR_DEFINED<<"_SMOOTHED_KERNEL_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_SMOOTHED_KERNEL_"<<RESET<<endl;
#endif

#ifdef _MODIFY_LIMITS_
  cout<<COLOR_DEFINED<<"_MODIFY_LIMITS_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_MODIFY_LIMITS_"<<RESET<<endl;
#endif


#ifdef _RANK_ORDERING_AB_INITIO_
  cout<<COLOR_DEFINED<<"_RANK_ORDERING_AB_INITIO_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_RANK_ORDERING_AB_INITIO_"<<RESET<<endl;
#endif

#ifdef _RANK_ORDERING_AB_INITIO_
#ifdef _RANK_ORDERING_AB_INITIO_
  cout<<COLOR_DEFINED<<"_RO_WITH_DELTA_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_RO_WITH_DELTA_"<<RESET<<endl;
#endif
#endif



#ifdef _USE_DM_IN_BAM_
  cout<<COLOR_DEFINED<<"_USE_DM_IN_BAM_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_DM_IN_BAM_"<<RESET<<endl;
#endif

cout<<endl;
this->So.message_screen("THETA-PROPERTIES");

#ifdef _USE_CWC_
  cout<<COLOR_DEFINED<<"_USE_CWC_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_CWC_"<<RESET<<endl;
#endif
#ifdef _USE_CWC_INSIDE_LOOP_
  cout<<COLOR_DEFINED<<"_USE_CWC_INSIDE_LOOP_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_CWC_INSIDE_LOOP_"<<RESET<<endl;
#endif
#ifdef _USE_MASS_KNOTS_
  cout<<COLOR_DEFINED<<"_USE_MASS_KNOTS_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_MASS_KNOTS_"<<RESET<<endl;
#endif
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
  cout<<COLOR_DEFINED<<"_USE_INVARIANT_TIDAL_FIELD_II_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_INVARIANT_TIDAL_FIELD_II_"<<RESET<<endl;
#endif
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
  cout<<COLOR_DEFINED<<"_USE_INVARIANT_TIDAL_FIELD_III_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_INVARIANT_TIDAL_FIELD_III_"<<RESET<<endl;
#endif
#ifdef _USE_TIDAL_ANISOTROPY_
  cout<<COLOR_DEFINED<<"_USE_TIDAL_ANISOTROPY_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_TIDAL_ANISOTROPY_"<<RESET<<endl;
#endif
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
  cout<<COLOR_DEFINED<<"_USE_INVARIANT_SHEAR_VFIELD_I_"<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_INVARIANT_SHEAR_VFIELD_I_"<<RESET<<endl;
#endif
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
  cout<<COLOR_DEFINED<<"_USE_INVARIANT_SHEAR_VFIELD_II_"<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_INVARIANT_SHEAR_VFIELD_II_"<<RESET<<endl;
#endif
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
  cout<<COLOR_DEFINED<<"_USE_INVARIANT_SHEAR_VFIELD_III_"<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_INVARIANT_SHEAR_VFIELD_III_"<<RESET<<endl;
#endif
#ifdef _USE_S2_
  cout<<COLOR_DEFINED<<"_USE_S2_"<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_S2_"<<RESET<<endl;
#endif
#ifdef _USE_S3_
  cout<<COLOR_DEFINED<<"_USE_S3_"<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_S3_"<<RESET<<endl;
#endif
#ifdef _USE_S2DELTA_
  cout<<COLOR_DEFINED<<"_USE_S2DELTA_"<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_S2DELTA_"<<RESET<<endl;
#endif
#ifdef _USE_DELTA2_
  cout<<COLOR_DEFINED<<"_USE_DELTA2_"<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_DELTA2_"<<RESET<<endl;
#endif
#ifdef _USE_DELTA3_
  cout<<COLOR_DEFINED<<"_USE_DELTA3_"<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_DELTA3_"<<RESET<<endl;
#endif
#ifdef _USE_NABLA2DELTA_
  cout<<COLOR_DEFINED<<"_USE_NABLA2DELTA_"<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_NABLA2DELTA_"<<RESET<<endl;
#endif


 cout<<endl;
this->So.message_screen("ASSIGNMENT");

#ifdef _ASSIGN_PROPERTY_
  cout<<COLOR_DEFINED<<"_ASSIGN_PROPERTY_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_ASSIGN_PROPERTY_"<<RESET<<endl;
#endif


#ifdef _ASSIGN_MASS_POST_
  cout<<COLOR_DEFINED<<"_ASSIGN_MASS_POST_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_ASSIGN_MASS_POST_"<<RESET<<endl;
#endif



#ifdef _ONLY_POST_PROC_
  cout<<COLOR_DEFINED<<"_ONLY_POST_PROC_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_ONLY_POST_PROC_"<<RESET<<endl;
#endif


#ifdef _ONLY_POST_PROC_
#ifdef _DO_NOT_CONVOLVE_
  cout<<COLOR_DEFINED<<"_DO_NOT_CONVOLVE_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_DO_NOT_CONVOLVE_"<<RESET<<endl;
#endif
#endif


#ifdef _ASSIGN_TO_CALIBRATION_
  cout<<COLOR_DEFINED<<"_ASSIGN_TO_CALIBRATION_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_ASSIGN_TO_CALIBRATION_"<<RESET<<endl;
#endif

#ifdef _ASSIGN_TO_REFERENCE_
  cout<<COLOR_DEFINED<<"_ASSIGN_TO_REFERENCE_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_ASSIGN_TO_REFERENCE_"<<RESET<<endl;
#endif


#ifdef _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_
  cout<<COLOR_DEFINED<<"_USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_"<<RESET<<endl;
#endif

#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_
  cout<<COLOR_DEFINED<<"_USE_MULTISCALE_PROPERTY_ASSIGNMENT_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_MULTISCALE_PROPERTY_ASSIGNMENT_"<<RESET<<endl;
#endif


#ifdef _USE_TRACERS_IN_CELLS_
  cout<<COLOR_DEFINED<<"_USE_TRACERS_IN_CELLS_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_TRACERS_IN_CELLS_"<<RESET<<endl;
#endif


#ifdef _USE_NUMBER_OF_NEIGHBOURS_
  cout<<COLOR_DEFINED<<"_USE_NUMBER_OF_NEIGHBOURS_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_NUMBER_OF_NEIGHBOURS_"<<RESET<<endl;
#endif

#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
  cout<<COLOR_DEFINED<<"_USE_MIN_DISTANCE_TO_NEIGHBOURS_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_MIN_DISTANCE_TO_NEIGHBOURS_"<<RESET<<endl;
#endif

#ifdef _USE_LOCAL_CLUSTERING_
  cout<<COLOR_DEFINED<<"_USE_LOCAL_CLUSTERING_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_LOCAL_CLUSTERING_"<<RESET<<endl;
#endif



#ifndef _ASSIGN_TO_REFERENCE_
#ifdef _ASSIGN_TO_CALIBRATION_
  cout<<COLOR_DEFINED<<"_ASSIGN_TO_CALIBRATION_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_ASSIGN_TO_CALIBRATION_"<<RESET<<endl;
#endif
#endif

#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_

#ifdef _USE_MULTISCALE_LEVEL_1_
  cout<<COLOR_DEFINED<<"_USE_MULTISCALE_LEVEL_1_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_EXCLUSION_CELL_LEVEL1_"<<RESET<<endl;
#endif

#ifdef _USE_MULTISCALE_LEVEL_2_
  cout<<COLOR_DEFINED<<"_USE_MULTISCALE_LEVEL_2_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_MULTISCALE_LEVEL_2_"<<RESET<<endl;
#endif


#ifdef _USE_MULTISCALE_LEVEL_3_
  cout<<COLOR_DEFINED<<"_USE_MULTISCALE_LEVEL_3_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_MULTISCALE_LEVEL_3_"<<RESET<<endl;
#endif

#ifdef _USE_MULTISCALE_LEVEL_4_
  cout<<COLOR_DEFINED<<"_USE_MULTISCALE_LEVEL_4_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_MULTISCALE_LEVEL_4_"<<RESET<<endl;
#endif
#ifdef _USE_MULTISCALE_LEVEL_4_

#ifdef _USE_EXCLUSION_CELL_LEVEL4_
  cout<<COLOR_DEFINED<<"_USE_EXCLUSION_CELL_LEVEL4_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_EXCLUSION_CELL_LEVEL4_"<<RESET<<endl;
#endif

#ifdef _USE_EXCLUSION_CELL_LEVEL4_HIGH_RES_
  cout<<COLOR_DEFINED<<"_USE_EXCLUSION_CELL_LEVEL4_HIGH_RES_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_EXCLUSION_CELL_LEVEL4_HIGH_RES_"<<RESET<<endl;
#endif

#endif

#endif


#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
  cout<<COLOR_DEFINED<<"_USE_MIN_SEPARATIONS_IN_CELLS_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_MIN_SEPARATIONS_IN_CELLS_"<<RESET<<endl;
#endif

#ifdef _USE_TRACERS_IN_CELLS_
  cout<<COLOR_DEFINED<<"_USE_TRACERS_IN_CELLS_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_TRACERS_IN_CELLS_"<<RESET<<endl;
#endif





#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
  cout<<COLOR_DEFINED<<"_SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_ "<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_"<<RESET<<endl;
#endif
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#ifdef _TOP_RANDOM_
  cout<<COLOR_DEFINED<<"_TOP_RANDOM_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_TOP_RANDOM_"<<RESET<<endl;
#endif
#ifdef _BOTTOM_RANDOM_
  cout<<COLOR_DEFINED<<"_BOTTOM_RANDOM_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_BOTTOM_RANDOM_"<<RESET<<endl;
#endif


#endif


#ifdef _COLLAPSE_RANDOMS_
  cout<<COLOR_DEFINED<<"_COLLAPSE_RANDOMS_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_COLLAPSE_RANDOMS_"<<RESET<<endl;
#endif
#ifdef _COLLAPSE_RANDOMS_AUX_
  cout<<COLOR_DEFINED<<"_COLLAPSE_RANDOMS_AUX_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_COLLAPSE_RANDOMS_AUX_"<<RESET<<endl;
#endif

#ifdef _COLLAPSE_RANDOMS_USING_EXCLUSION_
  cout<<COLOR_DEFINED<<"_COLLAPSE_RANDOMS_USING_EXCLUSION_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_COLLAPSE_RANDOMS_USING_EXCLUSION_"<<RESET<<endl;
#endif

#ifdef _APPLY_GLOBAL_EXCLUSION_
  cout<<COLOR_DEFINED<<"_APPLY_GLOBAL_EXCLUSION_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_APPLY_GLOBAL_EXCLUSION_"<<RESET<<endl;
#endif



#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
  cout<<COLOR_DEFINED<<"_USE_RS_AS_DERIVED_OBSERVABLE_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_RS_AS_DERIVED_OBSERVABLE_"<<RESET<<endl;
#endif

#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
  cout<<COLOR_DEFINED<<"_USE_SPIN_AS_DERIVED_OBSERVABLE_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_SPIN_AS_DERIVED_OBSERVABLE_"<<RESET<<endl;
#endif



  cout<<endl;
  this->So.message_screen("POWER");
#ifdef _POWER_
  cout<<COLOR_DEFINED<<"_POWER_"<<RESET<<endl;
  cout<<COLOR_DEFINED<<"        This directive is devoted for -p execution (i.e, measure powr spectgrum)_"<<RESET<<endl;
  cout<<COLOR_DEFINED<<"        For other options (e.g., -b), please undefine it."<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_POWER_"<<RESET<<RESET<<endl;
#endif

#ifdef _USE_ALL_PK_
  cout<<COLOR_DEFINED<<"_USE_ALL_PK_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_ALL_PK_"<<RESET<<endl;
#endif

#ifdef _USE_MASS_CUTS_PK_
  cout<<COLOR_DEFINED<<"_USE_MASS_CUTS_PK_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_MASS_CUTS_PK_"<<RESET<<endl;
#endif

#ifdef _USE_MASS_BINS_PK_
  cout<<COLOR_DEFINED<<"_USE_MASS_BINS_PK_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_USE_MASS_BINS_PK_"<<RESET<<endl;
#endif

#ifdef _REDSHIFT_SPACE_
  cout<<COLOR_DEFINED<<"_REDSHIFT_SPACE_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_REDSHIFT_SPACE_"<<RESET<<endl;
#endif
#ifdef _WRITE_MULTIPOLES_
  cout<<COLOR_DEFINED<<"_WRITE_MULTIPOLES_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_WRITE_MULTIPOLES_"<<RESET<<endl;
#endif



  cout<<endl;
  this->So.message_screen("OTHER ACTIONS");

#ifdef _BIN_ACCUMULATE_
  cout<<COLOR_DEFINED<<"_BIN_ACCUMULATE_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_BIN_ACCUMULATE_"<<RESET<<endl;
#endif


#ifdef _NGP2CIC_Y_
  cout<<COLOR_DEFINED<<"_NGP2CICL_Y_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_NGP2CIC_Y_"<<RESET<<endl;
#endif



#ifdef _DYNAMICAL_SAMPLING_
  cout<<COLOR_DEFINED<<"_DYNAMICAL_SAMPLING_"<<RESET<<endl;
#else
  cout<<COLOR_UNDEFINED<<"_DYNAMICAL_SAMPLING_"<<RESET<<endl;
#endif



}

//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
void Bam::set_Fourier_vectors()
{
  this->So.enter(__PRETTY_FUNCTION__);
#ifdef _FULL_VERBOSE_
  So.message_screen("Inizializing some Fourier vectors");
#endif
  ULONG NTT=(this->Nft)*(this->Nft)*(this->Nft/2+1);
  this->Kernel.resize(NTT, 1.0); //initialize to unity
  this->Power_REF.resize(this->Nft/2/this->ndel_data, 0);
  this->Power_NEW.resize(this->Nft/2/this->ndel_data, 0);
  this->power_ratio_unsmoothed.resize(this->Nft/2/this->ndel_data, 10.0);
  this->kvec.resize(this->Nft/2/this->ndel_data, 0);
#ifdef _USE_MASS_KNOTS_
  this->SKNOT_M_info.resize(this->NGRID, 0);
#endif
  this->Comp_conditional_PDF=this->params._Comp_conditional_PDF();
  this->Comp_joint_PDF=this->params._Comp_joint_PDF();
  So.DONE();
}
//##################################################################################
//#################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
void Bam::get_power_spectrum(string type)
{

   this->So.enter(__PRETTY_FUNCTION__);
#ifdef _FULL_VERBOSE_
  So.message_screen("Measuring power spectrum of",type);
#endif
  this->params.mass_assignment_scheme="NGP";
  this->params.MAS_correction=false;

  this->params.dir_output =this->Output_directory;  //This line is important
  // in the case in which the output dir is changing. Otherwise, the O(k) will be written in the dir read initially by the class Params

  this->params.Name_survey=type;

  if(this->step >0 && this->step <=this->N_iterations_Kernel)
    this->params.Name_survey=type+"_iteration"+to_string(this->step);


  if(type=="DM_REF" || type=="DM_REF_NEW")
    {
      this->params.input_type="density_grid";
      this->params.SN_correction=false;
      this->params.MAS_correction=true;

      if(0==this->iMAS_X)
        {
          this->params.mass_assignment_scheme="NGP";
          this->params.MAS_correction=false;
        }
      if(1==this->iMAS_X)
        this->params.mass_assignment_scheme="CIC";
      if(11==this->iMAS_X)
	{
          this->params.mass_assignment_scheme="CIC";
          this->params.MAS_correction=false;
	}
      else if (2==this->iMAS_X)
        this->params.mass_assignment_scheme="TSC";
      else if (3==this->iMAS_X)
        this->params.mass_assignment_scheme="PSC";


      PowerSpectrumF cPSF(this->params);
      cPSF.compute_power_spectrum_grid(this->delta_X_ini);
      So.DONE();
      cPSF.write_power_and_modes();
      this->Power_DM_REF.resize(this->Nft/2/this->ndel_data, 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<this->Power_DM_REF.size(); ++i)
        this->Power_DM_REF[i]=cPSF._pk0(i);
    }

  else if(type=="CROSS_TR_DM")
    {
      this->params.input_type="density_grid";
      this->params.SN_correction=true;
      this->params.mass_assignment_scheme="NGP";
      this->params.MAS_correction=false;
      this->params.measure_cross=true;
      PowerSpectrumF cPSF(this->params);
      cPSF.compute_cross_power_spectrum_grid(true,this->delta_X, this->delta_Y);
      So.DONE();
      cPSF.write_power_and_modes();
    }


  else if(type=="DM_iteration" || type=="DM_real" || type=="DM_KONV"|| type=="DM_RO")
    {
      this->params.input_type="delta_grid";
      this->params.SN_correction=false;
      this->params.mass_assignment_scheme="CIC";
      this->params.MAS_correction=true;
      PowerSpectrumF cPSF(this->params);
      if(type=="DM_RO")
    	{
	      cPSF.compute_power_spectrum_grid(this->delta_X_ini);
	     this->Power_DM_REF.resize(this->Nft/2/this->ndel_data, 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
	  for(int i=0;i<this->Power_DM_REF.size(); ++i)
	    this->Power_DM_REF[i]=cPSF._pk0(i);
	}
      else
        cPSF.compute_power_spectrum_grid(this->delta_X);
      So.DONE();
      cPSF.write_power_and_modes();
    }
  else if(type=="TR_MOCK" || type== "TR_MOCK_REALIZATION")//TR_MOCK refers to the itarations: TR_MOCK_REAL refers to the actuial mocks built from the learning process
    {
      if(this->Name_Property_Y=="COUNTS")
      	{
	        this->params.SN_correction=true;
    	    this->params.MAS_correction=false;
	        this->params.mass_assignment_scheme="NGP";
    	  }
      else
	     {
	       this->params.SN_correction=false;
         this->params.mass_assignment_scheme="CIC";
	       this->params.MAS_correction=false;
	     }    

      if(type== "TR_MOCK_REALIZATION")
        this->params.Name_survey="TR_MOCK";
     
#ifdef _CALIBRATION_WITHOUT_SN_
    	if(_COUNTS_==this->Name_Property_Y )
      	if(type=="TR_MOCK")
     	    if(this->step==0)
	           this->params.SN_correction=false;
#endif
	

      this->params.input_type="density_grid";

      PowerSpectrumF cPSF(this->params);

      cPSF.compute_power_spectrum_grid(this->delta_Y_new);
      this->Power_NEW.resize(this->Nft/2/this->ndel_data, 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<this->Power_NEW.size(); ++i)
        this->Power_NEW[i]=cPSF._pk0(i);  ///     +cPSF.shot_noise; Revisar esto: colocar arriba sn_cor =true, acá se los ponemos, pero lo escribimos con SN desde write_power_modes()

      So.DONE();
      cPSF.write_power_and_modes();
#ifndef _GET_BAM_REALIZATIONS_
      cPSF.write_power_and_modes(this->Output_directory+"power_mock.txt");
#endif
      this->kvec.resize(this->Nft/2/this->ndel_data, 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<this->Power_NEW.size() ;++i)
        this->kvec[i]=cPSF._kvector_data(i);
    }

  else if(type=="TR_MOCK_CATb")
    {
      this->params.SN_correction=true;
      this->params.input_type="density_grid";
      this->params.mass_assignment_scheme="CIC";
      this->params.MAS_correction=true;
      PowerSpectrumF cPSF(this->params);
      cPSF.compute_power_spectrum_grid(this->delta_Y_new);
      cPSF.write_power_and_modes();
    }

  else if(type=="TR_MOCK_CAT")
    {

      this->params.SN_correction=false;
      if(this->Name_Property_Y=="COUNTS")
        this->params.SN_correction=true;

      this->params.input_type="catalog";
      this->params.mass_assignment_scheme="TSC";
      this->params.MAS_correction=true;
      this->params.file_catalogue=this->patchy.fnameTRACERCAT;

      // Change this for the ordering in the param file might not be that of the one used by Patchy to write the catalog
      this->params.i_coord1_g=0;
      this->params.i_coord2_g=1;
      this->params.i_coord3_g=2;

      /*      PowerSpectrumF cPSF(this->params);
              cPSF.tracer_cat=this->tracer;
              cPSF.compute_power_spectrum(false, false);
              So.DONE();
              //the argument false above forces to write explicitely here the write_power and modes:
              cPSF.write_power_and_modes();
      */
    }

  else if(type=="TR_REF")
    {
      this->params.SN_correction=false;
      if(_COUNTS_==this->Name_Property_Y)
      	{
          this->params.SN_correction=true;
          this->params.mass_assignment_scheme="NGP";
          this->params.MAS_correction=false;
        }
       else
         {
           this->params.mass_assignment_scheme="CIC";
           this->params.MAS_correction=false;
         }
      
      
      this->params.input_type="density_grid";
#ifdef _USE_TRACER_HR_
      this->params.Nft=this->Nft_HR;
#endif
      PowerSpectrumF cPSF(this->params);
#ifdef _USE_TRACER_HR_
      cPSF.compute_power_spectrum_grid(this->delta_Y_HR);
#else
      cPSF.compute_power_spectrum_grid(this->delta_Y);
#endif
      // Return to the original value of Nft
      // and write ref power up the the NF of the Nft
#ifdef _USE_TRACER_HR_
      this->params.Nft=this->Nft;
#endif
      So.DONE();
      cPSF.write_power_and_modes();
      this->Power_REF.resize(this->Nft/2/this->ndel_data, 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<this->Nft/2/this->ndel_data;++i)
        this->Power_REF[i]=cPSF._pk0(i);

      this->kvec.resize(this->Nft/2/this->ndel_data, 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<this->Nft/2/this->ndel_data ;++i)
        this->kvec[i]=cPSF._kvector_data(i);
    }
  // For the DM we do not take the kvectors. We take them from the calc of TR power
}


//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//#define _use_more_chis_  //use this option to use the info from adjacen k-bins

void Bam::GetKernel(bool rejection, real_prec exponent)
{
  this->So.enter(__PRETTY_FUNCTION__);


#ifdef _USE_OMP_
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

#ifdef _use_more_chis_
  real_prec weight_back=0.15;
  real_prec weight_forw=0.15;
  real_prec weight_central=0.7;
#endif

#ifdef _use_more_chis_
#ifdef _FULL_VERBOSE_
  this->So.message_screen("Applying Metropolis-Hasting algorithm to generate BAM-Kernel (three-modes based)");
#endif
#else
#ifdef _FULL_VERBOSE_
this->So.message_screen("Applying Metropolis-Hasting algorithm to generate BAM-Kernel (single-mode based)");
#endif
#endif



#ifdef _USE_GNUPLOT_
    // PLOT OF THE POWER SPECTRUM
     vector<pair<real_prec, real_prec> > xy_pts_ref;
     vector<pair<real_prec, real_prec> > xy_pts_new;
     vector<pair<real_prec, real_prec> > xy_pts_dm;
     for(int i=0; i<kvec.size(); ++i) 
       xy_pts_ref.push_back(std::make_pair(this->kvec[i], log10(this->Power_REF[i])));
     for(int i=0; i<kvec.size(); ++i) 
       xy_pts_new.push_back(std::make_pair(this->kvec[i], log10(this->Power_NEW[i])));
     for(int i=0; i<kvec.size(); ++i) 
       xy_pts_dm.push_back(std::make_pair(this->kvec[i], log10(this->Power_DM_REF[i])));
     this->gp_kernel<<"set size 0.5,0.5\n";
     this->gp_kernel<<"set origin 0.5,0.0\n";
     this->gp_kernel<<"set log x \n";
     this->gp_kernel<<"set xlabel 'k [h / Mpc]' font 'Times-Roman,12'\n";
     this->gp_kernel<<"set ylabel 'log P(k) [(Mpc / h)³]' font 'Times-Roman,12'\n";
     this->gp_kernel<<"plot"<<this->gp_kernel.file1d(xy_pts_ref) << "w l lw 2 lt 8 title 'Reference',"<<this->gp_kernel.file1d(xy_pts_new) << " w l lw 2 lt 6 title 'Mock',"<<this->gp_kernel.file1d(xy_pts_dm)<< "w l lw 2 lt 3 title 'DM'"<<endl;
     xy_pts_ref.clear();
     xy_pts_ref.shrink_to_fit();
     xy_pts_new.clear();
     xy_pts_new.shrink_to_fit();
     xy_pts_dm.clear();
     xy_pts_dm.shrink_to_fit();
#endif




  if(false==this->use_iteration_ini) // initialized as false in bamrunner if iteration_ini >0
    {
      int nmodes_f=this->Nft/this->ndel_data/2;
      const gsl_rng_type *  T;
      gsl_rng * rng ;

      vector<real_prec> weight(nmodes_f,1.0);

      vector<real_prec> aux_power(nmodes_f,1.0);
      //#pragma omp parallel private(T, r)
      {
        gsl_rng_env_setup();
        ULONG seed=55457;//time(NULL);
        gsl_rng_default_seed=seed;
        T = gsl_rng_mt19937;//.gsl_rng_ranlux;
        rng = gsl_rng_alloc (T);


        // Select kernel according to the previous step. Update the previous kernel with the current step
        real_prec power_ratio=1.0;
        real_prec partial_ratio=1.00;
        real_prec deltak=2.*M_PI/this->Lbox;
        real_prec vol= pow(this->Lbox,3);
        real_prec sfac=1./(2.*M_PI*deltak*vol); // this factor is 2 / (4pì delta_k V)
        int counter_mh=0;
        int counter_residuals=0;
//#pragma omp parallel for
        real_prec residuals=0;
        real_prec residuals_unsigned=0;
        for(ULONG i=0;i<nmodes_f;++i)
          {
            real_prec kmode=this->kvec[i];
#ifdef _DO_BAM_CALIBRATION_
            real_prec Power_ref=this->Power_REF[i];
            real_prec Power_new=this->Power_NEW[i];
#ifdef _use_more_chis_
            real_prec Power_ref_back= i > 0 ? this->Power_REF[i-1]: 0 ;
            real_prec Power_new_back= i > 0 ? this->Power_NEW[i-1]: 0 ;
            real_prec Power_ref_forw= i< nmodes_f ? this->Power_REF[i+1]:0;
            real_prec Power_new_forw= i< nmodes_f ? this->Power_NEW[i+1]:0;
#endif
#else
            real_prec Power_ref=this->Power_REF_MW[i];
            real_prec Power_new=this->Power_NEW_MW[i];
#endif
            // I make this distiction, for there are approx methods with negative dm power, which lead to "nan" under the square root
            // in the case of the kernel for the DM.
            if(exponent >=0 && exponent<1.0)
              power_ratio=(Power_new == 0 || Power_ref == 0) ? 1.0 : pow(fabs(Power_ref/Power_new), exponent);
            else
              {
                power_ratio = Power_new == 0. ? 1.0 : Power_ref/Power_new;   //idelly for exponent = 1.0
                if(i< N_MODES)
                  partial_ratio+=power_ratio;
              }
            if(false==rejection)
              {
#ifndef _use_random_kernel_
                weight[i]=power_ratio;
#endif
                this->power_ratio_unsmoothed[i]=power_ratio;
              }
            else if(true==rejection)
              {
                real_prec kmode_squared=kmode*kmode;  // k²
                real_prec deltaVK_squared=kmode_squared+(1./12.)*deltak*deltak; // This is (1/Delta_K)²
                real_prec inv_deltaVK_squared=sfac/deltaVK_squared;
                // ***********************************
                // These applies if we use instead the likelihhod and take ratios of it.
                // The approach here resembles more a MCMC with  L=exp(-chi**2) for each kbin.
                // The variance at each k-bin is assumed to be Gaussian, and we neglect here shot-noise
//                real_prec sigma_squared =  Power_ref == 0? 1.0 :  inv_deltaVK_squared*(Power_new*Power_new); // this is sigma²
                real_prec sigma_squared =  Power_ref == 0? 1.0 :  inv_deltaVK_squared*(Power_ref+this->shot_noise_ref)*(Power_ref+this->shot_noise_ref); // this is sigma²

#ifdef _use_more_chis_
                real_prec sigma_back_squared = Power_ref_back == 0? 1.0 : inv_deltaVK_squared*(Power_ref_back+this->shot_noise_ref)*(Power_ref_back+this->shot_noise_ref);
                real_prec sigma_forw_squared = Power_ref_forw == 0? 1.0 : inv_deltaVK_squared*(Power_ref_forw+this->shot_noise_ref)*(Power_ref_forw+this->shot_noise_ref);
                real_prec new_H =0.5*(pow(Power_ref- Power_new, 2)/sigma_squared)*weight_central   ;
                new_H+= 0.5*(pow(Power_ref_back- Power_new_back,2)/sigma_back_squared)*weight_back;
                new_H+= 0.5*(pow(Power_ref_forw- Power_new_forw,2)/sigma_forw_squared)*weight_forw;
#else
                // proposed likelihood in the present step
                real_prec sigma_squared_new =  Power_ref == 0? 1.0 :  inv_deltaVK_squared*(Power_new*Power_new); // this is sigma²
                real_prec new_H = 0.5*(Power_ref- Power_new)*(Power_ref- Power_new)/sigma_squared_new;
#endif
                // power of the previous step obtained from the ratio ref/new saved in the previous iteration
                real_prec Power_old= static_cast<double>(Power_ref)/static_cast<double>(this->power_ratio_unsmoothed[i]);

                real_prec sigma_squared_old =  Power_ref == 0? 1.0 :  inv_deltaVK_squared*(Power_old*Power_old); // this is sigma²


#ifdef _use_more_chis_
                real_prec Power_old_back= i>0 ? static_cast<double>(Power_ref_back)/static_cast<double>(this->power_ratio_unsmoothed[i-1]): 0 ;
                real_prec Power_old_forw= i<nmodes_f? static_cast<double>(Power_ref_forw)/static_cast<double>(this->power_ratio_unsmoothed[i+1]):0 ;
                real_prec old_H = 0.5*(pow(Power_ref-Power_old,2)/sigma_squared)*weight_central;   // likelihood of the previous step
                old_H+= 0.5*(pow(Power_ref_back-Power_old_back,2)/sigma_back_squared)*weight_back;
                old_H+= 0.5*(pow(Power_ref_forw-Power_old_forw,2)/sigma_forw_squared)*weight_forw;
#else
                real_prec old_H = 0.5*(Power_ref-Power_old)*(Power_ref-Power_old)/sigma_squared_old;   // likelihood of the previous step
#endif
                // likelihood of the previous step
  //                double ratio_diff= exp(-static_cast<double>(new_H/10.0))/exp(-static_cast<double>(old_H/10.0)); // ratios between lilekihoods
                double ratio_diff= static_cast<double>(old_H)/static_cast<double>(new_H);  //ratios between chi²

                // ***********************************
                double xran= gsl_rng_uniform (rng);
                if(xran  < min(1.0, ratio_diff))
                  {
                    counter_mh++;
                    weight[i]=power_ratio; 
                    aux_power[i]=power_ratio;    //recorded just to print out
                  }
                else
                  {
#ifdef _use_random_kernel_
                    weight[i]=gsl_rng_uniform (rng);  //in the
//#else
//                    weight[i]=1.0; // this is not needed: if not accepted, the weight remains as inizialized above, i.e, =1. 
#endif
                    aux_power[i]=this->power_ratio_unsmoothed[i]; //recorded just to print out
                  }

                this->power_ratio_unsmoothed[i]=power_ratio;

                if(i>INITIAL_MODE_RESIDUALS)
                 {
                   residuals+=fabs(static_cast<real_prec>(power_ratio)-1.0);
                   residuals_unsigned+=static_cast<real_prec>(power_ratio)-1.0;
                   counter_residuals++;
                 }
            }
        }
        So.DONE();
#ifdef _FULL_VERBOSE_
        So.message_screen("Number of modes upgraded for Kernel = ", counter_mh);
        So.message_screen("Residuals at this iteration (%) = ",100.0*residuals/static_cast<real_prec>(counter_residuals));
        So.message_screen("Average of ratio P(k)_ref / P(k)_new = " , partial_ratio/static_cast<real_prec>(N_MODES));
        So.message_screen("Computed from fundamental mode up to k = ", this->kvec[N_MODES]);
        cout<<endl;
#else
        cout<<BLUE<<"Iteration "<<this->step<<RESET<<endl;
        cout<<BLUE<<"Residuals = "<<100.0*residuals/static_cast<real_prec>(counter_residuals)<<" %  \r"<<RESET<<endl;
        cout<<BLUE<<"Residuals(unsingned) = "<<100.0*residuals_unsigned/static_cast<real_prec>(counter_residuals)<<" %  \r"<<RESET<<endl;
#endif
        residuals_power.push_back(100.0*residuals/static_cast<real_prec>(counter_residuals));
        residuals_power_unsigned.push_back(100.0*residuals_unsigned/static_cast<real_prec>(counter_residuals));
        it_power.push_back(this->step);

    }


#ifdef _DO_BAM_CALIBRATION_
  string kernel_file_or=this->Output_directory+"RatioT_iteration"+to_string(this->step)+".txt";
  this->File.write_to_file(kernel_file_or, this->kvec,aux_power,power_ratio_unsmoothed);
  this->File.write_to_file(this->Output_directory+"RatioT.txt", this->kvec,aux_power,power_ratio_unsmoothed);


#ifdef _USE_GNUPLOT_
  // PLOT OF THE RESIDUALS
  std::vector<std::pair<double, double> > xy_pts_r;
  std::vector<std::pair<double, double> > xy_pts_ru;
  for(int i=0; i<it_power.size(); ++i) 
    xy_pts_r.push_back(std::make_pair(this->it_power[i], this->residuals_power[i]));
  for(int i=0; i<it_power.size(); ++i) 
    xy_pts_ru.push_back(std::make_pair(this->it_power[i], this->residuals_power_unsigned[i]));
  this->gp_kernel<<"set size 0.5,0.5\n";
  this->gp_kernel<<"set origin 0.,0.5\n";
  this->gp_kernel<<"set grid\n";
  this->gp_kernel<<"unset log\n";
  this->gp_kernel<<"set title 'Reference: "<<this->params._seed()<<"'\n";
  this->gp_kernel<<"set xlabel 'Iteration' font 'Times-Roman,12'\n";
  this->gp_kernel<<"set ylabel 'Residuals %' font 'Times-Roman,12'\n";
  this->gp_kernel<<"plot " << gp_kernel.file1d(xy_pts_r) << "w l lw 2 lt 2 title 'Absolute',"<<gp_kernel.file1d(xy_pts_ru)<< "w l lw 2 lt 5 title 'Relative'"<<endl;
  xy_pts_r.clear();
  xy_pts_r.shrink_to_fit();
  xy_pts_ru.clear();
  xy_pts_ru.shrink_to_fit();


  // PLOT OF THE RATIO
  for(int i=0; i<kvec.size(); ++i) 
    xy_pts_r.push_back(std::make_pair(this->kvec[i], aux_power[i]));
  this->gp_kernel<<"set size 0.5,0.5\n";
  this->gp_kernel<<"set origin 0.0,0.0\n";
  this->gp_kernel<<"set grid\n";
  this->gp_kernel<<"unset log y\n";
  this->gp_kernel<<"set log x\n";
  this->gp_kernel<<"set xlabel 'k [h / Mpc]' font 'Times-Roman,12'\n";
  this->gp_kernel<<"set ylabel 'Ratio T(k)' font 'Times-Roman,12'\n";
  this->gp_kernel<<"plot " << gp_kernel.file1d(xy_pts_r) << "w l lw 2 lt 7 notitle"<<endl;
  xy_pts_r.clear();
  xy_pts_r.shrink_to_fit();
#endif



  aux_power.clear(); //to release memory before going out of scope
  aux_power.shrink_to_fit();
#else
  string kernel_file_or=this->Output_directory+"RatioT_mass_assignment_iteration"+to_string(this->step_mass_assignemt)+".txt";
  this->File.write_to_file(kernel_file_or, this->kvec,aux_power,power_ratio_unsmoothed);
  aux_power.clear(); //to release memory before going out of scope
  aux_power.shrink_to_fit();
#endif

#ifdef _DO_BAM_CALIBRATION_
#ifdef _FULL_VERBOSE_
  So.message_screen("Updating BAM-Kernel");
#endif
#else
#ifdef _FULL_VERBOSE_
  So.message_screen("Updating Halo-Mass assignment Kernel");
#endif
#endif

  vector<real_prec> coords(this->Nft,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<coords.size() ;++i)
    coords[i]= (i<=Nft/2? static_cast<real_prec>(i): -static_cast<real_prec>(Nft-i));

  vector<int>nmodes(this->Power_NEW.size(), 0);
  vector<real_prec>kernel_updated(this->Power_NEW.size(), 0);// to print out the shell-averaged kernel
#ifndef _use_random_kernel_
  if(true==rejection)
    {
#endif
#pragma omp parallel for collapse(3)
      for(ULONG i=0;i< this->Nft; ++i)
        for(ULONG j=0;j< this->Nft; ++j)
          for(ULONG k=0;k< this->Nft/2+1; ++k)
            {
              ULONG ind=index_3d(i,j,k, this->Nft, this->Nft/2+1);
              real_prec kxc=coords[i];
              real_prec kyc=coords[j];
              real_prec kzc=coords[k];
              real_prec kv=_get_modulo(kxc,kyc,kzc);// sqrt(kxc*kxc + kyc*kyc + kzc*kzc);
              int kmod=static_cast<int>(floor(kv));
              if(kmod<this->kvec.size())
               {
#ifdef _use_random_kernel_
                 this->Kernel[ind]=weight[i]; // Random kernel
#else

#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                 this->Kernel[ind]*=weight[kmod]; //weight is 1 if no improvement, so the kernel is not updated;
#endif
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                  kernel_updated[kmod]+=this->Kernel[ind]; //weight is 1 of no improvement; the new kernel if it gets closer
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                  nmodes[kmod]++;
                }
            }
           
          
#ifndef _use_random_kernel_
    }
    else
      {
#pragma omp parallel for collapse(3)
        for(ULONG i=0;i< this->Nft; ++i)
          for(ULONG j=0;j< this->Nft; ++j)
            for(ULONG k=0;k< this->Nft/2+1; ++k)
              {
                ULONG ind=index_3d(i,j,k, this->Nft, this->Nft/2+1);
                real_prec kxc=coords[i];
                real_prec kyc=coords[j];
                real_prec kzc=coords[k];
                real_prec kv=_get_modulo(kxc,kyc,kzc);// sqrt(kxc*kxc + kyc*kyc + kzc*kzc);
                int kmod=static_cast<int>(floor(kv));
                if(kmod<this->Nft/2/this->ndel_data)
                  {
                    this->Kernel[ind] = weight[kmod];
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                    nmodes[kmod]++;
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                    kernel_updated[kmod]+=this->Kernel[ind]; //weight is 1 of no improvement; the new kernel if it gets closer
                  }
              }
        }
#endif
// ****************************************************************************************************************************
// Get the shell-averaged version of the Kernel in FOurier space
// ****************************************************************************************************************************
  for(ULONG i=0;i< kernel_updated.size(); ++i)
    kernel_updated[i]/=static_cast<real_prec>(nmodes[i]);
  nmodes.clear(); nmodes.shrink_to_fit();

//   ****************************************************************************************************************************
//   Fix the value of the kernel at the fundamental mode by assigning to that mode an average of the first three values of kernel
// if(this->step == this->N_iterations_Kernel)
   kernel_updated[0]=(kernel_updated[1]+kernel_updated[2]+kernel_updated[3]+kernel_updated[4])/4.0;
//   ****************************************************************************************************************************

#ifdef _DO_BAM_CALIBRATION_
  string file_kernel=this->Output_directory+"Kernel_Fourier_iteration"+to_string(this->step)+".txt";
  this->File.write_to_file(file_kernel, this->kvec, kernel_updated);
  this->File.write_to_file(this->Output_directory+"Kernel_Fourier.txt", this->kvec, kernel_updated);

#ifdef _USE_GNUPLOT_
  // PLOT OF THE KERNEL
  for(int i=0; i<kvec.size(); ++i) 
    xy_pts_r.push_back(std::make_pair(this->kvec[i], kernel_updated[i]));
  this->gp_kernel<<"set size 0.5,0.5\n";
  this->gp_kernel<<"set origin 0.5,0.5\n";
//  this->gp_kernel<<"set title 'Kernel'\n";
  this->gp_kernel<<"unset log y\n";
  this->gp_kernel<<"set log x\n";
  this->gp_kernel<<"set grid\n";
  this->gp_kernel<<"set xlabel 'k [h / Mpc]' font 'Times-Roman,12' \n";
  this->gp_kernel<<"set ylabel 'Kernel K(k)' font 'Times-Roman,12'\n";
  this->gp_kernel<<"plot" << gp_kernel.file1d(xy_pts_r) << "w l lw 2 lt 7 notitle"<<endl;
  xy_pts_r.clear();
  xy_pts_r.shrink_to_fit();
#endif

#else
  string file_kernel=this->Output_directory+"Kernel_Fourier_mass_assignment_iteration"+to_string(this->step_mass_assignemt)+".txt";
  this->File.write_to_file(file_kernel, this->kvec, kernel_updated, weight);
#endif

#ifdef _SMOOTHED_KERNEL_
#ifdef _SMOOTHED_KERNEL_LAST_ITERATION_
  if(this->step==this->N_iterations_Kernel)
    {
#endif
// Smooth twice shell-averaged kernel
      lin_smooth(this->kvec, kernel_updated);
//          lin_smooth(this->kvec, kernel_updated);
//          lin_smooth(this->kvec, kernel_updated);
//          lin_smooth(this->kvec, kernel_updated);
        //This line sets the kernek in the first mode to an average of the first three modes
      kernel_updated[0]=(kernel_updated[1]+kernel_updated[2]+kernel_updated[3])/3.0;
      string kernel_file=this->Output_directory+"Kernel_Fourier_smoothed_iteration"+to_string(this->step)+".txt";
      sal.open(kernel_file.c_str());
#ifdef _FULL_VERBOSE_
      So.message_screen("Writting smoothed version of kernel in file", kernel_file);
#endif
      sal.close();

#ifdef _SMOOTHED_KERNEL_LAST_ITERATION_
    }
#endif
// ****************************************************************************************************************************
// Assign smoothed shell averaged Kernel to 3D kernel:
// ****************************************************************************************************************************
#ifdef _USE_OMP_
#pragma omp parallel for collapse(3)
#endif
  for(ULONG i=0;i< this->Nft; ++i)
    for(ULONG j=0;j< this->Nft; ++j)
      for(ULONG k=0;k< this->Nft/2+1; ++k)
        {
          ULONG ind=index_3d(i,j,k, this->Nft, this->Nft/2+1);
          real_prec kxc=coords[i];
          real_prec kyc=coords[j];
          real_prec kzc=coords[k];
          real_prec kv=_get_modulo(kxc,kyc,kzc);
          int kmod=static_cast<int>(floor(kv));
          if(kmod<this->Nft/2/this->ndel_data)
            this->Kernel[ind]=kernel_updated[kmod];
        }
#endif // end of smooth

// ****************************************************************************************************************************
// All other vectors defined here and not released
// are destroyed her when going out of scope
// ****************************************************************************************************************************

#ifndef _TEST_THRESHOLDS_RESIDUALS_
      if(this->step == this->N_iterations_Kernel)
        this->File.write_array(this->Output_directory+"Bam_Kernel", this->Kernel);
#endif

  if(this->step == this->N_iterations_Kernel)
    {
#ifdef DOUBLE_PREC
     complex_prec *data_out= (complex_prec *)fftw_malloc(2*this->NTT*sizeof(real_prec));
#else
     complex_prec *data_out= (complex_prec *)fftwf_malloc(2*this->NTT*sizeof(real_prec));
#endif
// ****************************************************************************************************************************
// The original BAM code had only the real part assigned, with IMAG set to zero.
// ****************************************************************************************************************************
#pragma omp parallel for
    for(ULONG ind=0;ind< this->NTT ;++ind)
      {
        data_out[ind][REAL]=this->Kernel[ind];
        data_out[ind][IMAG]=0;
      }
    vector<real_prec>aux(this->NGRID,0);
    do_fftw_c2r(this->Nft, data_out,aux);
    this->File.write_array(this->Output_directory+"Bam_Kernel_config_space", aux);
    aux.clear(); aux.shrink_to_fit();

#ifdef DOUBLE_PREC
    fftw_free(data_out);
#else
    fftwf_free(data_out);
#endif
    }
  this->use_iteration_ini=false;
  }
else
    {
      auto out_it = std::find(std::begin(this->output_at_iteration), std::end(this->output_at_iteration), this->step);
      if((this->step==this->N_iterations_Kernel) || (out_it != std::end(this->output_at_iteration)))
      {
#ifdef _FULL_VERBOSE_
        So.message_screen("Reading Kernel from iteration", this->iteration_ini);
#endif
        this->Kernel.clear();
        this->Kernel.resize(this->NTT, 0.0);
        this->File.read_array(this->params._Input_Directory_BIAS_KERNEL()+"Bam_Kernel.dat", this->Kernel);
        this->use_iteration_ini=true; // we set it true for we will need it for the bias
      }
   }
}

//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
void Bam::Konvolve(vector<real_prec> &in, vector<real_prec>&out)
{


  this->So.enter(__PRETTY_FUNCTION__);

#ifdef _USE_OMP_
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

#ifdef _GET_BAM_REALIZATIONS_
#ifdef _FULL_VERBOSE_
  if(this->step_mass_assignemt==0)
    So.message_screen("Generating new DM density field by convolution of input DM with input Kernel");
  else
    So.message_screen("Generating new DM density field by convolution of DM with updated mass-weighterd kernel ");
#endif
#else
#ifdef _FULL_VERBOSE_
  So.message_screen("Generating new DM density field by convolution of input DM with updated Kernel");
#endif
#endif


#ifdef DOUBLE_PREC
  complex_prec *data_out= (complex_prec *)fftw_malloc(2*this->NTT*sizeof(real_prec));
#else
  complex_prec *data_out= (complex_prec *)fftwf_malloc(2*this->NTT*sizeof(real_prec));
#endif

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->NTT;++i)
   {
     data_out[i][REAL]=0;
     data_out[i][IMAG]=0;
   }
  do_fftw_r2c(this->Nft,in, data_out);

/*
double paux=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:paux)
#endif
  for(ULONG i=0;i<this->NGRID;i++)
   paux+=static_cast<double>(in[i]);
   if(true==isinf(paux)){
       So.message_warning("Not defined value found in container at function ",__PRETTY_FUNCTION__);
       So.message_warning("Line" ,__LINE__);
       So.message_warning("Code exits here");
       exit(0);
}
*/
#ifdef _EXTRAPOLATE_VOLUME_
  real_prec correction_factor = 1.0;
#else
  real_prec correction_factor = 1.0;
#endif
  double we=0;

#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:we)
#endif
  for(ULONG ind=0; ind< this->NTT ;++ind)
    {
      real_prec cor = this->Kernel[ind]*correction_factor;
      real_prec reald=data_out[ind][REAL]*cor;
      real_prec imagd=data_out[ind][IMAG]*cor;
      data_out[ind][REAL]=reald;
      data_out[ind][IMAG]=imagd;
      we+=cor;
    }
  do_fftw_c2r(this->Nft, data_out, out);
  // Here I have to correct for the normalization of the kernel
  // for the function  do_fftw_c2r returns the transform normalized by the NGRID, so I divide by NGRID and by multiply by 2 we
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->NGRID;i++)
    {
      real_prec out_s=out[i];
      out[i]=out_s/(static_cast<double>(this->NGRID)/static_cast<double>(2.0*we));
    }
  So.DONE();

#ifdef DOUBLE_PREC
  fftw_free(data_out);
#else
  fftwf_free(data_out);
#endif

  if(this->iteration_ini>0 && this->step!=this->iteration_ini)
    {
      auto out_it = std::find(std::begin(this->output_at_iteration), std::end(this->output_at_iteration), this->step);
      if (out_it != std::end(this->output_at_iteration))
        {
          // Write the kernel interpolated to 3D in Fourier space.
#ifdef DOUBLE_PREC
          complex_prec *kern= (complex_prec *)fftw_malloc(2*this->NTT*sizeof(real_prec));
#else
          complex_prec *kern= (complex_prec *)fftwf_malloc(2*this->NTT*sizeof(real_prec));
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(ULONG i=0;i<this->NTT;++i)
            {
              kern[i][REAL]=this->Kernel[i];
              kern[i][IMAG]=0;
            }

          vector<real_prec>aux(this->NGRID,0);
          do_fftw_c2r(this->Nft,kern,aux);
          // para no escribir lo que ha leído
          string file_kernel=this->Output_directory+"3DKernel_iteration"+to_string(this->step);
          this->File.write_array(file_kernel, aux);
          aux.clear();aux.shrink_to_fit();
#ifdef DOUBLE_PREC
          fftw_free(kern);
#else
          fftwf_free(kern);
#endif
          }
    }



}
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
void Bam::get_new_min_max_properties()
{

    //this->So.enter(__PRETTY_FUNCTION__);

  // These are the limits in the tracer info, which should remain constant
  // Only computed once in the first step (see below, in case we do not have -COUNTS:)

  // Limits on the reference
  this->s_mins.prop0=this->Ymin;
  this->s_maxs.prop0=this->Ymax;
  this->s_deltas.prop0=this->DELTAY;

  // This is only done once, as it is meant for the reference. THIS AHS TO BE REVISED, SINCE, WHEN WE USE CONTINIOUS VARIABLES, MION AND MAX ARE IN ANY CASE REFERRED TO BINS IN TH HISTOGRAMS
  // IN REALITY, this->Name_Property_Y MUST BE ALWAYS COUNTS, EVEN IF IT IS A CONTINIOUS QUANTITY, WE WILL COBNVERT IT TO A HISTOGRAM
  if(this->step==0)
    {
#ifdef _USE_MASS_FIELD_
      this->s_mins.prop0_mass=get_min(this->delta_Y_MASS);
      this->s_maxs.prop0_mass=get_max(this->delta_Y_MASS);
      this->s_deltas.prop0_mass=(this->s_maxs.prop0_mass-this->s_mins.prop0_mass)/static_cast<real_prec>(this->NY_MASS);
      So.message_screen("Maximim mass Tracer = ", this->s_maxs.prop0_mass);
      So.message_screen("Minimum mass  Tracer = ", this->s_mins.prop0_mass);
#endif

    }//close else of if _COUNTS_!=this->Name_Property_Y



#ifdef _MODIFY_LIMITS_

#ifdef _USE_DM_IN_BAM_
  // Gert new extremes for delta_X
  this->s_mins.prop1=get_min(this->delta_X);
  this->s_maxs.prop1=get_max(this->delta_X);
  this->s_deltas.prop1=(this->s_maxs.prop1-this->s_mins.prop1)/static_cast<real_prec>(this->new_nbins_x);
#ifdef _FULL_VERBOSE_
  So.message_screen("Minimum log(2+ð) DM = ", this->s_mins.prop1);
  So.message_screen("Maximum log(2+ð) DM = ", this->s_maxs.prop1);
  //  So.message_screen("N_bins  ð DM = ", this->new_nbins_x);
#endif
#endif


#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
  this->s_mins.prop4=get_min(this->cwclass.Invariant_TF_II);
  this->s_maxs.prop4=get_max(this->cwclass.Invariant_TF_II);
  this->s_deltas.prop4=(this->s_maxs.prop4-this->s_mins.prop4)/static_cast<real_prec>(N_C_BIN1);
#ifdef _FULL_VERBOSE_
  So.message_screen("Maximim InvTF2 = ", this->s_maxs.prop4);
  So.message_screen("Minimum InvTF2 = ", this->s_mins.prop4);
  So.message_screen("Delta InvTF2 = ", this->s_deltas.prop4);
#endif

#elif defined _USE_DELTA2_
  this->s_mins.prop4=get_min(this->cwclass.DELTA2);
  this->s_maxs.prop4=get_max(this->cwclass.DELTA2);
  this->s_deltas.prop4=(this->s_maxs.prop4-this->s_mins.prop4)/static_cast<real_prec>(N_C_BIN1);
  So.message_screen("Maximim ð² = ", this->s_maxs.prop4);
  So.message_screen("Minimum ð² = ", this->s_mins.prop4);
#endif

#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
  this->s_mins.prop5=get_min(this->cwclass.Invariant_TF_III);
  this->s_maxs.prop5=get_max(this->cwclass.Invariant_TF_III);
  this->s_deltas.prop5=(this->s_maxs.prop5-this->s_mins.prop5)/static_cast<real_prec>(N_C_BIN2);
#ifdef _FULL_VERBOSE_
  So.message_screen("Maximim InvTF3 = ", this->s_maxs.prop5);
  So.message_screen("Minimum InvTF3 = ", this->s_mins.prop5);
  So.message_screen("Delta InvTF3 = ", this->s_deltas.prop5);
#endif

#elif defined _USE_DELTA3_
  this->s_mins.prop5=get_min(this->cwclass.DELTA3);
  this->s_maxs.prop5=get_max(this->cwclass.DELTA3);
  this->s_deltas.prop5=(this->s_maxs.prop5-this->s_mins.prop5)/static_cast<real_prec>(N_C_BIN2);
  So.message_screen("Maximim ð³ = ", this->s_maxs.prop5);
  So.message_screen("Minimum ð³ = ", this->s_mins.prop5);
#endif

#ifdef _USE_INVARIANT_TIDAL_FIELD_IV_
  this->s_mins.prop6=get_min(this->cwclass.Invariant_TF_IV);
  this->s_maxs.prop6=get_max(this->cwclass.Invariant_TF_IV);
  this->s_deltas.prop6=(this->s_maxs.prop6-this->s_mins.prop6)/static_cast<real_prec>(N_C_BIN3);
  So.message_screen("Maximim InvTFIII = ", this->s_maxs.prop6);
  So.message_screen("Minimum InvTFII = ", this->s_mins.prop6);
  So.message_screen("Delta InvTFII = ", this->s_deltas.prop6);
#elif _USE_TIDAL_ANISOTROPY_
  this->s_mins.prop6=get_min(this->cwclass.Tidal_Anisotropy);
  this->s_maxs.prop6=get_max(this->cwclass.Tidal_Anisotropy);
  this->s_deltas.prop6=(this->s_maxs.prop6-this->s_mins.prop6)/static_cast<real_prec>(N_C_BIN3);
  So.message_screen("Maximim Tidal Anisotropy = ", this->s_maxs.prop6);
  So.message_screen("Minimum Tidal Anisotropy = ", this->s_mins.prop6);

#elif defined _USE_S2_
  this->s_mins.prop6=get_min(this->cwclass.S2);
  this->s_maxs.prop6=get_max(this->cwclass.S2);
  this->s_deltas.prop6=(this->s_maxs.prop6-this->s_mins.prop6)/static_cast<real_prec>(N_C_BIN3);
  So.message_screen("Maximim s² = ", this->s_maxs.prop6);
  So.message_screen("Minimum s² = ", this->s_mins.prop6);
#endif

#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
  this->s_mins.prop7=get_min(this->cwclass.Invariant_VS_I);
  this->s_maxs.prop7=get_max(this->cwclass.Invariant_VS_I);
  this->s_deltas.prop7=(this->s_maxs.prop7-this->s_mins.prop7)/static_cast<real_prec>(N_CV_BIN1);
  So.message_screen("Maximim InvVS1 = ", this->s_maxs.prop7);
  So.message_screen("Minimum InvVS1 = ", this->s_mins.prop7);
#elif defined _USE_NABLA2DELTA_
  this->s_mins.prop7=get_min(this->cwclass.N2D);
  this->s_maxs.prop7=get_max(this->cwclass.N2D);
  this->s_deltas.prop7=(this->s_maxs.prop7-this->s_mins.prop7)/static_cast<real_prec>(N_CV_BIN1);
  So.message_screen("Maximim Nabla²ð = ", this->s_maxs.prop7);
  So.message_screen("Minimum Nabla²ð = ", this->s_mins.prop7);
#endif

#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
  this->s_mins.prop8=get_min(this->cwclass.Invariant_VS_II);
  this->s_maxs.prop8=get_max(this->cwclass.Invariant_VS_II);
  this->s_deltas.prop8=(this->s_maxs.prop8-this->s_mins.prop8)/static_cast<real_prec>(N_CV_BIN2);
  So.message_screen("Maximim InvVSII = ", this->s_maxs.prop8);
  So.message_screen("Minimum InvVSII = ", this->s_mins.prop8);
#elif defined _USE_S2DELTA_
  this->s_mins.prop8=get_min(this->cwclass.S2DELTA);
  this->s_maxs.prop8=get_max(this->cwclass.S2DELTA);
  this->s_deltas.prop8=(this->s_maxs.prop8-this->s_mins.prop8)/static_cast<real_prec>(N_CV_BIN2);
  So.message_screen("Maximim s²ð = ", this->s_maxs.prop8);
  So.message_screen("Minimum s²ð = ", this->s_mins.prop8);
#endif

#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
  this->s_mins.prop9=get_min(this->cwclass.Invariant_VS_III);  //not yet assigned to a container
  this->s_maxs.prop9=get_max(this->cwclass.Invariant_VS_III);
  this->s_deltas.prop9=(this->s_maxs.prop9-this->s_mins.prop9)/static_cast<real_prec>(N_CV_BIN3);
  So.message_screen("Maximim InvVSIII = ", this->s_maxs.prop9);
  So.message_screen("Minimum InvVSIII = ", this->s_mins.prop9);

#elif defined _USE_S3_
  this->s_mins.prop9=get_min(this->cwclass.S3);
  this->s_maxs.prop9=get_max(this->cwclass.S3);
  this->s_deltas.prop9=(this->s_maxs.prop9-this->s_mins.prop9)/static_cast<real_prec>(N_CV_BIN3);
  So.message_screen("Maximim s³ = ", this->s_maxs.prop9);
  So.message_screen("Minimum s³ = ", this->s_mins.prop9);
#endif

#else    // if _MODIFY_LIMITS_is not defined,


  real_prec xmin_temp=get_min(this->delta_X);
  real_prec xmax_temp=get_max(this->delta_X);
#ifdef _FULL_VERBOSE_
  if(this->Xmax< xmax_temp)
  {
    So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Xmax in delta below the nominal value");
    cout<<"this->Xmax="<<this->Xmax<<"  Current="<<xmax_temp<<endl;
  }
  if(this->Xmin> xmin_temp)
    So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Xmin in delta_X above the nominal value");
#endif

  this->s_mins.prop1=this->Xmin;
  this->s_mins.prop4=C1_MIN;
  this->s_mins.prop5=C2_MIN;
  this->s_mins.prop6=C3_MIN;
  this->s_mins.prop7=CV1_MIN;
  this->s_mins.prop8=CV2_MIN;
  this->s_mins.prop9=CV3_MIN;

  this->s_maxs.prop1=this->Xmax;
  this->s_maxs.prop4=C1_MAX;
  this->s_maxs.prop5=C2_MAX;
  this->s_maxs.prop6=C3_MAX;
  this->s_maxs.prop7=CV1_MAX;
  this->s_maxs.prop8=CV2_MAX;
  this->s_maxs.prop9=CV3_MAX;

  this->s_deltas.prop1=this->DELTAX;
  this->s_deltas.prop4=DELTA_C1;
  this->s_deltas.prop5=DELTA_C2;
  this->s_deltas.prop6=DELTA_C3;
  this->s_deltas.prop7=DELTA_CV1;
  this->s_deltas.prop8=DELTA_CV2;
  this->s_deltas.prop9=DELTA_CV3;
#endif   // end of ifdef _MODIFY_LIMITS_
}

// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################

void  Bam::get_min_max_X_Y()
{
  // This function redefines the values of delta_X_min etc in case. If this is not called, the input values are used
  // The filds delta X and delta Y ahave been already, if requiested, converted to overdensities
    this->So.enter(__PRETTY_FUNCTION__);

  if(true==this->Redefine_limits)
    {
      real_prec num_in_log_x = true==this->Convert_Density_to_Delta_X ? NUM_IN_LOG: 0.;
      real_prec num_in_log_y = true==this->Convert_Density_to_Delta_Y ? NUM_IN_LOG: 0.;

      vector<real_prec>AUX(this->delta_X.size());
      if(this->Scale_X=="linear")
        {
          this->Xmin=get_min<real_prec>(this->delta_X);
          this->Xmax=get_max<real_prec>(this->delta_X);
        }
      else
        if(this->Scale_X=="log"){
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(ULONG i=0;i<this->delta_X.size();++i)AUX[i]=log10(num_in_log_x+this->delta_X[i]);
          this->Xmin=get_min<real_prec>(AUX);
          this->Xmax=get_max<real_prec>(AUX);
        }

      if(this->Scale_Y=="linear")
        {
          this->Ymin=get_min<real_prec>(this->delta_Y);
          this->Ymax=get_max<real_prec>(this->delta_Y);
        }
      else
        if(this->Scale_Y=="log"){

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(ULONG i=0;i<this->delta_Y.size();++i)
            AUX[i]=log10(num_in_log_y+this->delta_Y[i]);
          this->Ymin=get_min<real_prec>(AUX);
          this->Ymax=get_max<real_prec>(AUX);
        }
      AUX.clear();  // to release memory before going out of scope
      AUX.shrink_to_fit();
    }
  else
    {
      if(this->Scale_X=="log")
        {
          this->Xmin=this->ldelta_X_min;
          this->Xmax=this->ldelta_X_max;
        }
      else
        {
          this->Xmin=this->delta_X_min;
          this->Xmax=this->delta_X_max;
        }
      if(this->Scale_Y=="log")
        {
          this->Ymin=this->ldelta_Y_min;
          this->Ymax=this->ldelta_Y_max;
        }
      else
        {
          this->Ymin=this->delta_Y_min;
          this->Ymax=this->delta_Y_max;
        }
    }
}
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################

#ifdef _USE_VELOCITIES_
void Bam::read_bam_files(string file_X, string file_Y, string file_Y_HR, string file_Y_mass, string file_Y_sat_frac, string file_X_ref, string file_Vx, string file_Vy, string file_Vz)
#else
  void Bam::read_bam_files(string file_X, string file_Y, string file_Y_HR, string file_Y_mass, string file_Y_sat_frac, string file_X_ref)
#endif
{

   this->So.enter(__PRETTY_FUNCTION__);


   int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);

#ifdef _FULL_VERBOSE_
  cout<<endl;
  this->So.message_screen("Reading input *reference* files");
  cout<<endl;
#endif

  ULONG NGRID_NEW;
  ULONG NGRID_NEW_HR;
#ifdef _EXTRAPOLATE_VOLUME_
  NGRID_NEW=static_cast<ULONG>(this->Nft_low*this->Nft_low*this->Nft_low);
#else
#ifdef _USE_TRACER_HR_
  NGRID_NEW=this->NGRID;
  NGRID_NEW_HR=static_cast<ULONG>(this->Nft_HR*this->Nft_HR*this->Nft_HR);
#else
  NGRID_NEW=this->NGRID;
  NGRID_NEW_HR=this->NGRID;
#endif
#endif

#ifndef _EXTRAPOLATE_VOLUME_
  this->delta_X.resize(NGRID_NEW,0);
  this->delta_X_ini.resize(NGRID_NEW,0);
  this->File.read_array(file_X,this->delta_X);

    //andres
  cout<<RED<<this->delta_X[1255]<<"   "<<this->delta_X.size()<<"   "<<NGRID_NEW<<RESET<<endl;

#ifdef _CONVERT_CIC_TO_NGP_
  convert_cic_to_ngp(delta_X,delta_X);
#endif
#endif

#ifdef _USE_TRACER_HR_
  this->delta_Y_HR.resize(NGRID_NEW_HR,0);
  this->File.read_array_t<PrecType_Y>(file_Y_HR, this->delta_Y_HR);
#endif

  //  this->delta_Y_new.resize(this->NGRID,0); //Isn't this resized in get_mock?
  this->delta_Y.resize(this->NGRID,0);
  this->File.read_array_t<PrecType_Y>(file_Y, this->delta_Y);


#ifdef _SHOW_EMPTY_CELLS_
  ULONG empty_cells=0;
  for(ULONG i=0;i < this->NGRID; ++i)
    if(this->delta_Y[i]==0)
      empty_cells++;
  So.message_screen("Number of empty cells from number density field =", empty_cells);
#endif

#ifdef _USE_MASS_TRACERS_
#ifdef _USE_MASS_FIELD_
  this->delta_Y_MASS.resize(this->NGRID,0);
  this->File.read_array_t<PrecType_Y>(file_Y_mass, this->delta_Y_MASS);
#endif

#ifdef _USE_LOG_MASS_
  this->Mean_density_Y_MASS=1;
#else
  /*
    this->Mean_density_Y_MASS=get_mean(this->delta_Y_MASS);
    get_overdens(this->delta_Y_MASS,this->Mean_density_Y_MASS, this->delta_Y_MASS);
  */

#ifdef _USE_MASS_FIELD_
#ifdef _SHOW_EMPTY_CELLS_
  empty_cells=0;
  for(ULONG i=0;i < this->NGRID; ++i)
    {
      if(this->delta_Y_MASS[i]<1.0)
        empty_cells++;
      this->delta_Y_MASS[i]=log10(NUM_IN_LOG+ this->delta_Y_MASS[i]/MASS_SCALE);
    }
  So.message_screen("Number of empty cells from mass density field =", empty_cells);
#endif
  So.message_screen("Using log(2+total central mass in cells). Check line ",__LINE__);
#endif

#endif
#endif



#ifdef _NGP2CIC_Y_
  convert_ngp_to_cic(this->delta_Y, this->delta_Y);
#endif
  // this->So.message_screen("Minimum Number of Y", get_min(delta_Y));
  // this->So.message_screen("Maximum Number of Y", get_max(delta_Y));
  // VELOCITIES: TAKEN FROM PATCHY OR READ
#ifdef _USE_VELOCITIES_
  this->Velx_X.resize(this->NGRID,0);
  this->Vely_X.resize(this->NGRID,0);
  this->Velz_X.resize(this->NGRID,0);
  this->File.read_array(file_Vx, this->Velx_X);
  this->File.read_array(file_Vy, this->Vely_X);
  this->File.read_array(file_Vz, this->Velz_X);
#endif

#ifdef _EXCHANGE_X_Y_DENSITY_
  So.message("Exchanging axis X-Y in density field");
  exchange_xy(this->Nft,this->delta_X,this->delta_X);
#endif

#ifdef _USE_VELOCITIES_
#ifdef _EXCHANGE_X_Y_VEL_
  So.message_screen("Exchanging axis X-Y in velocity field");
  exchange_xy(this->Nft,this->Velx_X,this->Velx_X);
  exchange_xy(this->Nft,this->Vely_X,this->Vely_X);
  exchange_xy(this->Nft,this->Velz_X,this->Velz_X);
  So.DONE();
#endif
#endif

#ifdef _KONV_
  vector<real_prec>kernel(this->NTT,0);
  this->File.read_array(this->params._Input_Directory_BIAS_KERNEL()+"Bam_Kernel.dat",kernel);
  convolvek(this->Nft,this->delta_X, kernel,this->delta_X);
#endif
  // Initialize the DM density field with the input (REF) density field
  this->delta_X_ini=this->delta_X;



}
//  #################################################################################################################
//  #################################################################################################################
//  #################################################################################################################
//  #################################################################################################################
//  #################################################################################################################
//  #################################################################################################################
//  #################################################################################################################
//  #################################################################################################################
//  #################################################################################################################


void Bam::analyze_input()
{
  So.enter(__PRETTY_FUNCTION__);

  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);

#ifdef _FULL_VERBOSE_
    So.message_screen("Getting some statistics from input files: ");
    cout<<endl;
#endif

// *********************************************************************************************
// Analyze input file for X variables
// *********************************************************************************************

  real_prec nmean_X=0.;

#ifndef _EXTRAPOLATE_VOLUME_
  nmean_X=get_nobjects(this->delta_X);
  this->N_objects_X=nmean_X;
#ifdef _FULL_VERBOSE_
  if(this->Name_Property_X==_COUNTS_)
    So.message_screen("Total number of X objects =", nmean_X);
#endif

  nmean_X=static_cast<real_prec>(nmean_X)/static_cast<real_prec>(this->NGRID);
#ifdef _FULL_VERBOSE_
  if(_COUNTS_==this->Name_Property_X)
    So.message_screen("Mean number of X objects in cells =", nmean_X);
  else
    So.message_screen("Mean property X in cells =", nmean_X);
#endif

#ifdef _FULL_VERBOSE_
  if(_COUNTS_==this->Name_Property_X)
    So.message_screen("Mean number X density =", nmean_X*static_cast<real_prec>(this->NGRID)/pow(this->Lbox,3), "(Mpc / h )⁻³");
  else
    So.message_screen("Mean property X density =", nmean_X*static_cast<real_prec>(this->NGRID)/pow(this->Lbox,3));
#endif


  if(this->iMAS_X==0) // if DM commes in NGP, we can make a pdf
    {
      int NPART=static_cast<int>(get_max<real_prec>(this->delta_X)); // Number of particles in cells
      this->PDF_NC_X.resize(NPART, 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i< this->delta_X.size();++i)
        if(static_cast<int>(this->delta_X[i])<NPART)
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
          this->PDF_NC_X[static_cast<int>(this->delta_X[i])]++;

      if(true==this->Write_PDF_number_counts)
        {
          string fileX=this->Output_directory+"PDF_NC"+"_X_"+this->XNAME+"_"+this->Name_Property_X+"_MASX"+to_string(this->iMAS_X)+"_Nft"+to_string(this->Nft)+"_SmoothingScale"+to_string(this->smscale)+"_z"+to_string(this->redshift)+"_LambdaTH"+to_string(this->lambdath)+"_CW"+to_string(this->tstruct)+"_CWclass.txt";

          this->File.write_to_file_i(fileX,this->PDF_NC_X);
        }
      this->nmax_x_onecell=NPART;
#ifdef _FULL_VERBOSE_
      So.message_screen("Maximum number of DM particles in one cell =", this->nmax_x_onecell);
      So.message_screen("Estimated Poisson signal-to-noise ratio from DM=", sqrt(this->nmax_x_onecell));
#endif
    }
#ifdef MOCK_MODE
#ifdef _GET_POWER_REFS_
  this->get_power_spectrum("DM_REF");  //gets power from delta_X
#endif
#endif
#endif  // END EXTRAPOLATE

#ifdef _USE_X_COMPLEMENT_I_
#endif

// *********************************************************************************************
// *********************************************************************************************
// Analyze input file for Y variables
// *********************************************************************************************
  real_prec nmean_Y=0;
  nmean_Y=get_nobjects(this->delta_Y);

  if(_COUNTS_==this->Name_Property_Y)
    {
      this->N_objects_Y=nmean_Y;
#ifdef _EXTRAPOLATE_VOLUME_
      So.message_screen("Note that these are properties of the reference simulation (i.e, smaller volume)");
#endif

#ifdef _FULL_VERBOSE_
      cout<<endl;
      So.message_screen("Total number of Y objects =", nmean_Y);
#endif
    }

  real_prec LLBOX;
  ULONG NNGRID;

#ifdef _EXTRAPOLATE_VOLUME_
  LLBOX=this->Lbox_low;
  NNGRID= static_cast<ULONG>(this->Nft_low*this->Nft_low*this->Nft_low);
#else
  LLBOX=this->Lbox;
  NNGRID= this->NGRID;
#endif

  nmean_Y=static_cast<real_prec>(nmean_Y)/static_cast<real_prec>(NNGRID);
#ifdef _FULL_VERBOSE_

  if(_COUNTS_==this->Name_Property_Y)
    So.message_screen("Mean number of Y objects in cells =", nmean_Y);
  else
    So.message_screen("Mean property Y in cells =", nmean_Y);

  this->Mean_density_Y=nmean_Y;

  if(_COUNTS_==this->Name_Property_Y)
    So.message_screen("Mean Y number density =", nmean_Y*static_cast<real_prec>(NNGRID)/pow(LLBOX,3), "(Mpc / h )⁻³");
  else
    {
      So.message_screen("Mean Y property-density =", nmean_Y*static_cast<real_prec>(NNGRID)/pow(LLBOX,3));
    }

#endif
  this->new_Name_Property_Y=this->Name_Property_Y;

  if(_COUNTS_==this->Name_Property_Y && ZERO==this->iMAS_Y)
    {
      this->nmax_y_onecell=static_cast<int>(get_max<real_prec>(this->delta_Y)); // Maximum number of particles in cells
      // Here we get the pdf of the counts in order to measure the mean occupation number
      // thee value is to be allocated in the variable this->mean_number_x_onecell;
      this->mean_number_y_onecell=static_cast<int>(get_mean_occupation_number<real_prec>(this->nmax_y_onecell,  this->delta_Y));
#ifdef _FULL_VERBOSE_
      this->So.message_screen("Maximum number of Y particles in one cell =", this->nmax_y_onecell);
      this->So.message_screen("Average number of Y particles in one cell =", this->mean_number_y_onecell );
      this->So.message_screen("Estimated Maximum Poisson signal-to-noise ratio =", sqrt(this->nmax_y_onecell));
#endif
      ofstream nxo;
        nxo.open(file_one_cell);
        nxo<<this->nmax_y_onecell<<"\t"<<this->mean_number_y_onecell<<endl;
        nxo.close();
#ifdef _USE_MASS_FIELD_
        this->So.message_screen("Maximum tracer mass in one cell =", get_max<real_prec>(this->delta_Y_MASS));
#endif
        if(true==this->Write_PDF_number_counts)
          {

            this->PDF_NC_Y.resize(this->nmax_y_onecell, 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
            for(ULONG i=0;i<NNGRID;++i)
              {
                int nob=static_cast<int>(this->delta_Y[i]);

#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                this->PDF_NC_Y[nob]++;
              }
            string fileY=this->Output_directory+"PDF_NC_Y_REF_"+this->new_Name_Property_Y+"_MASY"+to_string(this->iMAS_Y)+"_Nft"+to_string(this->Nft)+"_z"+to_string(this->redshift)+".txt";
            this->File.write_to_file_i(fileY,this->PDF_NC_Y);
            this->PDF_NC_Y.clear();this->PDF_NC_Y.shrink_to_fit();
	         }

    }

#ifdef _FULL_VERBOSE_
  cout<<endl;
#endif

#ifdef MOCK_MODE
#ifndef _EXTRAPOLATE_VOLUME_
#ifdef _GET_POWER_REFS_
  this->get_power_spectrum("TR_REF");
#endif
#ifdef _USE_TRACER_HR_
  this->delta_Y_HR.clear();
  this->delta_Y_HR.shrink_to_fit();
#endif
#endif
#endif

#ifdef _GET_POWER_REFS_
//  this->get_power_spectrum("CROSS_TR_DM"); // NOT WORKING, MEMCHINK
#endif
#ifdef MOCK_MODE
#ifndef _GET_BAM_REALIZATIONS_
  real_prec lss_bias=0;
  int ncount=0;

  if(this->Power_DM_REF.size()==0 || this->Power_REF.size()==0)
    {
      So.message_warning("Container Power_DM_REF is not resized. Check preprocessor directive _GET_POWER_REFS_. BAM ends here");
      exit(0);
    }

#pragma omp parallel for reduction(+:lss_bias, ncount)
  for(int i=0;i<N_MODES;++i)
    if(this->Power_DM_REF[i]!=0)
      {
      	lss_bias += this->Power_REF[i]/this->Power_DM_REF[i];
	       ncount++;
      }
  lss_bias/=static_cast<real_prec>(ncount);
#ifdef _FULL_VERBOSE_
  So.message_screen("Large-Scale Bias: P(k)_tracer / P(k)_dm  = ", sqrt(lss_bias));
  cout<<endl;
#endif

#ifdef _UNDER_BIASED_
  if(lss_bias<1)
    {
      this->used_once=false;
      So.message_screen("Treating under-biased tracers");
      cout<<endl;
    }
#endif
#endif
#endif
#ifdef _RUN_TEST_
  dark_matter_to_halos_analytical(this->delta_X, this->delta_Y);
  exit(0);
#endif

// *********************************************************************************************
// Preparation of the DM density field
// *********************************************************************************************
  if(this->N_iterations_dm>0)
    {
      So.message_screen("Initial itearation for the DM density field");

      // Leamos el campo de densidad DM de la referencia, medimos el espectro y
      // definimos un power_REF
      this->params.Name_survey="DM_REF";
      this->params.MAS_correction= false;
      PowerSpectrumF dPSF(this->params);

      dPSF.compute_power_spectrum_grid(this->delta_X_REF_PDF); //ojo, el argumento debe ser DENSIDAD
      dPSF.write_power_and_modes();
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<this->Power_REF.size();++i)
        this->Power_REF[i]=dPSF._pk0(i);

      this->params.Name_survey="DM";
      PowerSpectrumF dPSFo(this->params);
      dPSFo.compute_power_spectrum_grid(this->delta_X); //ojo, el argumento debe ser DENSIDAD
      dPSFo.write_power_and_modes();
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<this->Power_REF.size();++i)
        this->Power_NEW[i]=dPSFo._pk0(i);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<this->Power_REF.size() ;++i)
        this->kvec[i]=dPSFo._kvector_data(i);

      for(int iteration_dm = 0; iteration_dm < this->N_iterations_dm ; ++iteration_dm) //THIS LOOP MUST BE DEPRECATED
        {
          So.message_screen("Iteration on the DM field", iteration_dm);
          // i) Get power of DM, getting ready for the kernel. NO rejection in this case, use exponent 0.5
          GetKernel(false, 0.5);
          vector<real_prec>auxv(this->NGRID, 0);
          auxv=this->delta_X;
         // v) Convolve the DM field field with the kernel obtaiende from the reference
          this->Konvolve(auxv,this->delta_X);
          auxv.clear();
          auxv.shrink_to_fit();
          this->get_power_spectrum("DM_iteration");
          this->File.write_array(this->Output_directory+"density_DM_converted_dm_iteration_"+to_string(iteration_dm), this->delta_X);
        }
    }

  // *********************************************************************************************
  // Assign the density field delta_X_ini to *density*  delta_X
  // Below, if delta_X is to be converted to overdensity, we reassign delta_X_ini to overdensity
  // *********************************************************************************************
  this->delta_X_ini=this->delta_X;  // BAM starts with the outcome of this loop. This line is important.
  // *********************************************************************************************
  // *********************************************************************************************
#ifdef _USE_BINARY_MASK_
  vector<real_prec> binary_mask(this->NGRID,0);
  this->File.read_array(this->Name_binary_mask, binary_mask);
#endif
// *********************************************************************************************
  // Convert density to delta for both the reference and the approximated field
  // *********************************************************************************************
  if(true==this->Convert_Density_to_Delta_X)
    {
      this->new_Name_Property_X="DELTA";
#ifdef _ADD_NOISE_
      // With this we add Poisson noise to the DM density field.
      // The noise is added by increasing the mean number density from that of the original filed
      // to a new valie MEAN_NEW, from wihch we construct a new density field
      // MEAN_NEW*1+delta) and use this value as mean to draw Poisson distributed values.
      // Note that, when measuring power spectrum, if we subtract the Poisson shot-noise from this,
      // we obtain the same input field.
      // Hence if we do that operation (subtract Poiss sn) we better use another distribution to add this
      // variance.
      gsl_rng_env_setup();
      gsl_rng_default_seed=75;
      const gsl_rng_type *  T= gsl_rng_ranlux;
      gsl_rng * r = gsl_rng_alloc (T);
      get_overdens(this->delta_X,  nmean_X, this->delta_X);  // get overdensity
      So.message_screen("Converting DENSITY -> Poisson(DENSITY):");
      gsl_real proba=0.8;
      real_prec expo=0.98;
        real_prec mean_new=0;
      for(ULONG i=0;i<delta_X.size();++i)
            mean_new+=pow(1.+delta_X[i],expo);
      mean_new/=static_cast<real_prec>(delta_X.size());

      for(ULONG i=0;i<delta_X.size();++i)
//         delta_X[i]=gsl_ran_poisson(r,MEAN_NEW*(1.+delta_X[i]));   //recover density
         delta_X[i]=gsl_ran_negative_binomial(r,proba, (MEAN_NEW/mean_new)*pow(1.+delta_X[i],expo));   //recover density
    //    delta_X[i]=(nmean_X/mean_new)*pow(1.+delta_X[i],1.8);   //recover density
      So.DONE();
       {
        this->params.input_type="density_grid";
        this->params.Name_survey="DM_Poisson";
        this->params.SN_correction=true;
        this->params.mass_assignment_scheme="CIC";
        this->params.MAS_correction=true;
        PowerSpectrumF dPSFa(this->params);
        dPSFa.compute_power_spectrum_grid(this->delta_X); //ojo, el argumento debe ser DENSIDAD
        dPSFa.write_power_and_modes();
      }
#endif

#ifdef _FULL_VERBOSE_
      So.message_screen("Converting DENSITY into DELTA for X:");
#endif

#ifndef _USE_BINARY_MASK_
#ifdef _ADD_NOISE_
      get_overdens(this->delta_X,  this->delta_X);
#else
      get_overdens(this->delta_X, nmean_X, this->delta_X);
#ifdef _GIVE_POWER_
      So.message_screen("Adding power on small scales to overdensity");
      give_power(this->Lbox, this->Nft , this->delta_X,this->delta_X);
      So.DONE();
      this->params.input_type="delta_grid";
      this->params.Name_survey="DM_extra";
      this->params.SN_correction=false;
      this->params.mass_assignment_scheme="CIC";
      this->params.MAS_correction=true;
      PowerSpectrumF dPSFa(this->params);
      dPSFa.compute_power_spectrum_grid(this->delta_X); //ojo, el argumento debe ser DENSIDAD
      dPSFa.write_power_and_modes();
#endif
#endif
#else
      get_overdens(this->delta_X, binary_mask, this->delta_X,true);
#endif
      // Reassign the delta_X_ini to *overdensities*  delta_X
      this->delta_X_ini=this->delta_X;
#ifdef _USE_BINARY_MASK_
      So.message_screen("Weighting Delta X with selection");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<this->delta_X.size();++i)
        delta_X[i]*=binary_mask[i];
      this->So.DONE();
#endif
    }
  else
    {
      this->new_Name_Property_X=this->Name_Property_X;
#ifdef _USE_BINARY_MASK_
      if(this->Name_Property_X=="DELTA")
        {
          So.message_screen("Weighting Delta X with selection");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(ULONG i=0;i<this->delta_X.size();++i)
            delta_X[i]*=binary_mask[i];
          this->So.DONE();
        }
#endif
    }

// *********************************************************************************************
// Convert the TRACER to delta
// *********************************************************************************************
    if(true==this->Convert_Density_to_Delta_Y)
      {
#ifdef _FULL_VERBOSE_
        So.message_screen("Converting DENSITY into DELTA for Y:");
#endif
        this->new_Name_Property_Y="DELTA";

#ifdef _USE_BINARY_MASK_
        get_overdens(this->delta_Y, binary_mask,this->delta_Y, false);
        binary_mask.clear();
        binary_mask.shrink_to_fit();
#else
        get_overdens(this->delta_Y, this->delta_Y);
#endif
      }
  else // IF NET, THE NEW NAME REMAINS THE SAME
    this->new_Name_Property_Y=this->Name_Property_Y;
}
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
// *********************************************************************************************
// This function converts the DM DELTA fields to log 2+delta and compute the CWC if requested
// The BAM kernel is initialized
  // *********************************************************************************************
void Bam::get_BAM_DM()
{

  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
  string rest_file="_Nft"+to_string(this->Nft)+"_z"+to_string(this->redshift)+".txt";
  real_prec num_in_log_x= NUM_IN_LOG;
  // **********************************************************************
  this->get_min_max_X_Y();// THIS CAN BE REPLACED by finding min and max at converting time
  // **********************************************************************
  vector<real_prec>xbaux(this->NX, 0);
  vector<real_prec>pdf_in(this->NX, 0); //contains the pdf of DM in each iteration, before updating it with teh convolution
#ifdef _WRITE_PDF_
  this->pdf_ref.resize(this->NX, 0);
  this->pdf_ini.resize(this->NX, 0);

  for(int i=0;i<this->NX; ++i)
    xbaux[i]=this->Xmin+static_cast<real_prec>(i+0.5)*(this->Xmax-this->Xmin)/(static_cast<real_prec>(this->NX));

  if(this->step==0)
    {
#ifdef _FULL_VERBOSE_
          So.message_screen("Computing PDF from log(1+delta)  dark matter");
#endif
      calc_pdf("log", this->NGRID,this->NX, this->Xmax, this->Xmin, this->delta_X, this->pdf_ref);
      string filex=this->Output_directory+"pdf_X_"+this->XNAME+"_MASX"+to_string(this->iMAS_X)+rest_file;
      this->File.write_to_file(filex, xbaux,pdf_ref);
    }
#endif
  // *********************************************************************************************
  // At this point, we had already said delta_X_ini =delta_x WITHOUT transforming to logs
  // (we convert delta_X to log10(num_in_log+delta_X) in the next function get_BIAS)
  // *********************************************************************************************
#ifdef _UNDER_BIASED_
  if(false==this->used_once)
    {
      So.message_screen("Using TR as DM for b<1 tests. Trick to treat under-biased populations");
      So.message_screen("Test in line", __LINE__, ", function ",__PRETTY_FUNCTION__);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i = 0;i < this->NGRID ;++i)
        this->delta_X_ini[i] = (this->delta_X[i]<-1 ?  0 :  log10(num_in_log_x+ static_cast<real_prec>(this->delta_Y[i])));
      this->used_once=true;
    }
#endif

  // *********************************************************************************************
  // Do the CWC class for the first iteration, or in the last , when te target field is loaded
  // *********************************************************************************************
  if(this->step==0)
    {
#ifdef _RANK_ORDERING_AB_INITIO_
#ifdef _FULL_VERBOSE_
          So.message_screen("Ready to perform rank ordering of DM to target distribution");
#endif
          this->delta_X_REF_PDF.resize(this->NGRID,0);
          string file_X_ref_pdf=this->Input_Directory_X_REF+this->Name_Catalog_X_REF_PDF;
          this->File.read_array_t<PrecType_X>(file_X_ref_pdf, this->delta_X_REF_PDF); // Read target DM field for RO
#ifdef _RO_WITH_DELTA_
	   	  get_overdens(this->delta_X_REF_PDF, this->delta_X_REF_PDF);
#endif
          pdf_in=this->pdf_ref;  // This was already computed above, so we do not replicate it here
#ifdef _FULL_VERBOSE_
          So.message_screen("Measuring pdf of log(1+delta) DM target reference");
#endif
          calc_pdf("log", this->NGRID,this->NX, this->Xmax, this->Xmin, this->delta_X_REF_PDF, this->pdf_ref);
          So.DONE();
#ifdef _WRITE_PDF_
          string filex=this->Output_directory+"pdf_X_REF_"+this->XNAME+"_MASX"+to_string(this->iMAS_X_REF_PDF)+rest_file;
          this->File.write_to_file(filex, xbaux,this->pdf_ref);
#endif
          this->delta_X_REF_PDF.clear();
          this->delta_X_REF_PDF.shrink_to_fit();
#ifdef _FULL_VERBOSE_
          So.message_screen("Executing rank ordering from DM to DM-target");
#endif
          rankorder(ZERO, xbaux, this->NGRID, this->NX,  this->Xmax, this->Xmin, this->delta_X_ini, pdf_in, this->pdf_ref);
          So.DONE();
          this->get_power_spectrum("DM_RO"); // from here, Power_DM_REF is loaded
          // Get pdf of the rank-ordered and write it to file
#ifdef _FULL_VERBOSE_
          So.message_screen("Measuring pdf of delta DM rank-ordered");
#endif
          calc_pdf("log", this->NGRID,this->NX, this->Xmax, this->Xmin, this->delta_X_ini, pdf_in);
          So.DONE();
#ifdef _WRITE_PDF_
          filex=this->Output_directory+"pdf_rank_ordered_iteration"+to_string(this->step)+"_X_"+this->XNAME+"_MASX"+to_string(this->iMAS_X)+rest_file;
          this->File.write_to_file(filex, xbaux,pdf_in);
#endif
          real_prec lss_bias=0;
          int ncount=0;
          for(int i=0;i<N_MODES;++i)
           if(this->Power_DM_REF[i]>0)
             {
               lss_bias += this->Power_REF[i]/this->Power_DM_REF[i];
               ncount++;
             }
            lss_bias/=static_cast<real_prec>(ncount);
#ifdef _FULL_VERBOSE_
          So.message_screen("New large-Scale Bias: P(k)_tracer / P(k)_dm_RO  = ", sqrt(lss_bias));
          cout<<endl;
#endif
#endif // END RANK ORDERING AB INITIO


#ifdef _USE_CWC_
      this->cwclass.do_CWC(this->delta_X_ini);
#ifdef _USE_MASS_KNOTS_
      this->cwclass.get_Mk_collapsing_regions(this->delta_X_ini,  static_cast<real_prec>(this->N_objects_X)/static_cast<real_prec>(this->NGRID));
#endif

#else
#if defined _USE_MASS_KNOTS_
      this->cwclass.do_CWC(this->delta_X_ini);
      this->cwclass.get_Mk_collapsing_regions(this->delta_X_ini,  static_cast<real_prec>(this->N_objects_X)/static_cast<real_prec>(this->NGRID));
#endif

#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_INVARIANT_TIDAL_FIELD_IV_)|| defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_ELLIPTICITY_) || defined (_USE_PROLATNESS_) || defined (_USE_PWEB_)
      this->cwclass.do_CWC(this->delta_X_ini);  // @at the raw sampling, uses the initial (i.e, the original density field to get the CWC and eigenvalues)
#endif
#endif   //endif _USE_CWC_

#if (defined (_USE_NABLA2DELTA_) || defined (_USE_S2_) || defined (_USE_DELTA2_)  || defined (_USE_S2DELTA_) || defined (_USE_DELTA3_) || defined (_USE_S3_)) && (!defined (_USE_CWC_)) && (!defined (_USE_PWEB_))
      this->cwclass.get_bias_terms(this->delta_X_ini);
#endif
    } // end if step==0
  else    // steps >= 1
    {
// *********************************************************************************************
// Compute the kernel
// *********************************************************************************************
      this->GetKernel(true, KERNEL_INDEX);
// #ifdef _HYDROTEST_
//       // ********hydro test, *******************************
//       this->cwclass.Kernel.clear();
//       this->cwclass.Kernel.shrink_to_fit();
//       this->cwclass.Kernel.resize(this->NTT, 0.0);
//       this->cwclass.Kernel=this->Kernel;
//       // ********************************************
// #endif
// *********************************************************************************************
// Convolve the Dm with the kernel and output delta_X
// with this we always convolve the original overdensity field
// *********************************************************************************************
//       if (this->step==this->N_iterations_Kernel)
//            this->File.write_array(this->Output_directory+"DM_DELTA_NOconvolved_iteration"+to_string(this->step), this->delta_X_ini);

//andres
  cout<<RED<<this->Kernel[1255]<<"  "<<this->delta_X_ini[1255]<<__LINE__<<RESET<<endl;
      this->Konvolve(this->delta_X_ini, this->delta_X);
 //andres
  cout<<RED<<this->Kernel[1255]<<"  "<<this->delta_X[1255]<<__LINE__<<RESET<<endl;


      // Get the PDF of the convolved field


#ifdef _WRITE_PDF_
      this->pdf_ite.resize(this->NX, 0);
#ifdef _FULL_VERBOSE_
      So.message_screen("Measuring pdf of log 1+ delta DM convolved");
#endif
      calc_pdf("log",this->NGRID,this->NX, this->Xmax, this->Xmin, this->delta_X, this->pdf_ite);

      int index= (this->step <=this->N_iterations_Kernel)  ?  this->step  : this->step - (this->N_iterations_Kernel)+1;
      string label_aux = this->step <= this->N_iterations_Kernel ? "_iteration": "_realization";

      string afilex=this->Output_directory+"pdf_convolved"+label_aux+to_string(index)+"_X_"+this->XNAME+"_MASX"+to_string(this->iMAS_X)+rest_file;
      this->File.write_to_file(afilex, xbaux,this->pdf_ite);
#endif



#ifdef _USE_CWC_INSIDE_LOOP_

#ifdef _USE_CWC_
      this->cwclass.do_CWC(this->delta_X);
#ifdef _USE_MASS_KNOTS_
      this->cwclass.get_Mk_collapsing_regions(this->delta_X,  static_cast<real_prec>(this->N_objects_X)/static_cast<real_prec>(this->NGRID));
#endif

#elif !defined _USE_CWC_
#ifdef _USE_MASS_KNOTS_
      this->cwclass.do_CWC(this->delta_X);
      this->cwclass.get_Mk_collapsing_regions(this->delta_X,  static_cast<real_prec>(this->N_objects_X)/static_cast<real_prec>(this->NGRID));
#endif

#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_INVARIANT_TIDAL_FIELD_IV_) || defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_ELLIPTICITY_) || defined (_USE_PROLATNESS_) || defined (_USE_PWEB_)
      this->cwclass.do_CWC(this->delta_X);  // @at iterations,  use the convolved density field to get the CWC
#endif

#endif

#endif

#if (defined (_USE_NABLA2DELTA_) || defined (_USE_S2_)  || defined (_USE_DELTA2_) || defined (_USE_S2DELTA_) || defined (_USE_DELTA3_) || defined (_USE_S3_) ) && (!defined (_USE_CWC_))
      this->cwclass.get_bias_terms(this->delta_X);
#endif

#ifdef _WRITE_DM_DENSITY_FIELD_
      auto out_it = std::find(std::begin(this->output_at_iteration), std::end(this->output_at_iteration), this->step);
      if (this->step==this->N_iterations_Kernel || out_it != std::end(this->output_at_iteration))
        this->File.write_array(this->Output_directory+"DM_DELTA_convolved_iteration"+to_string(this->step), this->delta_X);
      this->File.write_array(this->Output_directory+"DM_DELTA_convolved", this->delta_X);
#endif


#ifdef _GET_POWER_REFS_
    this->get_power_spectrum("DM_KONV");
#endif
    }//end if step>0
  // At the end of this function, dela_X denotes an overdensity

}

//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################

#ifdef MOCK_MODE
void Bam::get_X_function()
{

  this->So.enter(__PRETTY_FUNCTION__);

#ifdef _USE_OMP_
 int NTHREADS=_NTHREADS_;
 omp_set_num_threads(NTHREADS);
#endif

  //HEWRE WE HAVE TO USE tracer_ref to read the reference catalog
#ifdef _FULL_VERBOSE_
  So.message_screen("*************************************************************************");
  So.message_screen("*************************************************************************");
#ifdef _USE_VMAX_AS_OBSERVABLE_
  So.message_screen("***Measuring conditional Vmax Function from reference tracer catalog*****");
#elif defined _USE_MASS_AS_OBSERVABLE_
  So.message_screen("***Measuring conditional Mass Function from reference tracer catalog*****");
#endif
  So.message_screen("****as a function of DM properties of the reference DM density field*****");

  So.message_screen("*************************************************************************");
  cout<<endl;
#endif

  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************

  int N_a=N_C_BIN1;
#ifdef _USE_TOTAL_MASS_IN_CELL_
  N_a = N_BINS_TOTAL_MASS_IN_CELL;
#endif

  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************
  int N_b=N_C_BIN2;
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
  N_b = N_NEIGHBOURS_MAX;
#elif defined (_USE_MIN_DISTANCE_TO_NEIGHBOURS_) || defined (_USE_LOCAL_CLUSTERING_)
  N_b=N_BINS_MIN_DIST_TO_NEI;
#endif

  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************
  int N_c = N_C_BIN3;
#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
  N_c = N_BINS_MIN_SEP_IN_CELLS;
#endif
  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************

  int N_v= N_CV_BIN1;
#ifdef _USE_TRACERS_IN_CELLS_
  N_v = N_BINS_TRACERS_IN_CELLS;
#endif

 // ************************RESIZE VECTORS FOR HISTOGRAMS*************************

  ULONG LENGHT_AB=N_v*N_a* N_b* N_c* N_CV_BIN2*N_CV_BIN3;

  LENGHT_AB*=this->n_sknot_massbin * this->n_cwt * this->n_vknot_massbin * this->n_cwv*this->NX;


  // ********************************************************************************
#ifndef _ASSIGN_TO_REFERENCE_
  int N_x = this->params._NMASSbins();
  ULONG LENGHT_AB_ONE= LENGHT_AB*N_x;
#endif
#ifdef _USE_MASS_AS_OBSERVABLE_
   real_prec lm_min=this->params._LOGMASSmin();
#endif

// *********************************************************************************************
// Get min and max of the differnet properties invovled
// *********************************************************************************************
  this->get_new_min_max_properties();

// *********************************************************************************************
// Open the ref catalog and read it from this->tracer[].Halo.mass
// *********************************************************************************************
#ifdef _FULL_VERBOSE_
  cout<<endl;
  So.message_screen("******** IMPORTANT: the reference is currently being read from input cat file, see line", __LINE__);
#endif

/*
  this->params.i_coord1_g=0;
  this->params.i_coord2_g=1;
  this->params.i_coord3_g=2;  // esto es improtante, pues estamos leyendo el file con mcut ya hecho
  this->params.i_mass_g=7;
#ifdef _USE_VMAX_TRACERS_
  this->params.i_vmax_g=9;
#endif
*/

  this->tracer_ref.type_of_object="TRACER_REF";
  this->tracer_ref.set_params_catalog(this->params);
  string newcat=this->params._file_catalogue();

  // *********************************************************************************************
  //read catalog passing as argument the file and the mininum mass requested
  // *********************************************************************************************

#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
#ifdef _SET_CAT_WITH_MASS_CUT_
     this->tracer_ref.read_catalog(newcat,pow(10,params._LOGMASSmin())*params._MASS_units());
#elif defined (_SET_CAT_WITH_VMAX_CUT_)
  this->tracer.read_catalog(newcat,params._VMAXmin());
#else
  this->tracer_ref.read_catalog(newcat,0);
#endif
#else
#ifdef _SET_CAT_WITH_MASS_CUT_
  this->tracer_ref.read_catalog(newcat,pow(10,params._LOGMASSmin())*params._MASS_units(),static_cast<real_prec>(BIG_NUMBER));
#elif defined (_SET_CAT_WITH_VMAX_CUT_)
  this->tracer_ref.read_catalog(newcat,params._VMAXmin(),static_cast<real_prec>(BIG_NUMBER));
#endif
#endif


#ifdef _ASSIGN_TO_REFERENCE_
     this->tracer.Halo.resize(this->tracer.NOBJS);
     for(ULONG i=0; i<this->tracer.NOBJS; ++i)
       {
         this->tracer.Halo[i].coord1=this->tracer_ref.Halo[i].coord1;
         this->tracer.Halo[i].coord2=this->tracer_ref.Halo[i].coord2;
         this->tracer.Halo[i].coord3=this->tracer_ref.Halo[i].coord3;
         this->tracer.Halo[i].GridID=this->tracer_ref.Halo[i].GridID;
       }
#endif



  this->tracer_ref.get_property_function(this->Output_directory+"tracer_ref_abundance.txt");

#if defined _USE_NEIGHBOURS_ || defined _GET_DIST_MIN_SEP_REF_ || defined _GET_DIST_MIN_SEP_MOCK_
  So.message_screen("Identifying Neighbouring Cells (this is done once)");
  this->ncells_info.resize(this->NGRID);
  get_neighbour_cells(this->Nft, N_CELLS_BACK_FORTH, this->ncells_info);
  So.DONE();
  cout<<endl;
#endif



  // ********************************************************************************
  // ********************************************************************************
  // ********************************************************************************
  // *****************************Deal with reference catalogs and fields ***************

  // Here we open the reference DM catalog: we have to
  // i ) get delta
  // ii) convolve with kernel from BAM
  // iII) get mins and max
  // if inside iterative mass assignment, convolve with kernel computed from mass power spectra
  // iv) do Iweb or Tweb classification
  // v) convert to log(num_in_log+delta)
#ifdef _FULL_VERBOSE_
  cout<<endl;
#endif
  // ********************************************************************************
  this->cwclass_ref.set_params_cwclass(this->params);
  this->cwclass_ref.s_cosmo_pars=this->s_cosmo_pars;

  // Ideally Read the reference only in the first iter  ation, all containers with *aux* are not meant to be kept in memmory, they are like dummy
#ifdef _FULL_VERBOSE_
  So.message_screen("*************************************");
  So.message_screen("********Reading Reference DM*********");
  So.message_screen("*************************************");
#endif
  this->delta_dm_aux_mem.clear();
  this->delta_dm_aux_mem.shrink_to_fit();
  this->delta_dm_aux_mem.resize(this->NGRID,0);  //Keep this untouched, will be used along the iterations
  File.read_array(this->Input_Directory_X+this->Name_Catalog_X,this->delta_dm_aux_mem);
  this->mean_aux=get_mean(this->delta_dm_aux_mem);
  get_overdens(this->delta_dm_aux_mem,this->mean_aux, this->delta_dm_aux_mem);


#ifdef _RANK_ORDERING_AB_INITIO_
#ifdef _FULL_VERBOSE_
       So.message_screen("Executing rank ordering from DM to DM-target");
#endif
   this->delta_X_REF_PDF.resize(this->NGRID,0);
   string file_X_ref_pdf=this->Input_Directory_X_REF+this->Name_Catalog_X_REF_PDF; // TBDep
   this->File.read_array_t<PrecType_X>(file_X_ref_pdf, this->delta_X_REF_PDF); // Read target DM field for RO
#ifdef _RO_WITH_DELTA_
   get_overdens(this->delta_X_REF_PDF, this->delta_X_REF_PDF); //a better bispect is found if the rank ordering is done to the CIC of number counts, not the delta thereof
#endif
   vector<real_prec>xbaux(this->NX, 0);
   vector<real_prec>pdf_in(this->NX, 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int i=0;i<this->NX; ++i)
     xbaux[i]=this->Xmin+static_cast<real_prec>(i+0.5)*(this->Xmax-this->Xmin)/(static_cast<real_prec>(this->NX));
   this->pdf_ref.resize(this->NX, 0);
   calc_pdf("log", this->NGRID,this->NX, this->Xmax, this->Xmin,delta_dm_aux_mem, pdf_in);
   calc_pdf("log", this->NGRID,this->NX, this->Xmax, this->Xmin, this->delta_X_REF_PDF, this->pdf_ref);
   rankorder(ZERO, xbaux, this->NGRID, this->NX,  this->Xmax, this->Xmin, delta_dm_aux_mem, pdf_in, this->pdf_ref);
   this->delta_X_REF_PDF.clear();
   this->delta_X_REF_PDF.shrink_to_fit();

#endif   // end for _RANK_ORDERING_AB_INITIO_


#ifndef _DO_NOT_CONVOLVE_
  this->Konvolve(this->delta_dm_aux_mem,this->delta_dm_aux_mem);  // COnvolve with BAM Kernel
//  File.write_array(this->Output_directory+"DM_field",dm_ref);
  this->Kernel.clear();
  this->Kernel.shrink_to_fit();
#endif


#ifdef _USE_CWC_
      this->cwclass_ref.do_CWC(this->delta_dm_aux_mem);   //get the CW info
#ifdef _USE_MASS_KNOTS_
      this->cwclass_ref.get_Mk_collapsing_regions(this->delta_dm_aux_mem,this->mean_aux);  //get the MK info
#endif //use_mass_knots

#elif !defined _USE_CWC_ 
#if defined (_USE_MASS_KNOTS_)  ||  defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_INVARIANT_TIDAL_FIELD_IV_)|| defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_ELLIPTICITY_) || defined (_USE_PROLATNESS_) || defined (_USE_PWEB_)
   this->cwclass_ref.do_CWC(this->delta_dm_aux_mem);   //get the CW info
#endif
#endif

#ifdef _FULL_VERBOSE_
  cout<<endl;
  So.message_screen("Transforming delta_ref -> log10(2+delta_ref). Line ", __LINE__);
#endif

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i = 0;i < this->NGRID ;++i)  //TRANSFORM DELTA TO LOG10(NUM_IN_LOG + DELTA)
    this->delta_dm_aux_mem[i] = (this->delta_dm_aux_mem[i]<-1 ? 0  :  log10(NUM_IN_LOG+ this->delta_dm_aux_mem[i]));
  So.DONE();

  // ********************************************************************************
#ifdef _USE_TRACERS_IN_CELLS_
  vector<real_prec> REF_DEN_FIELD(this->NGRID,0);
  this->File.read_array(this->params._Input_Directory_Y()+this->params._Name_Catalog_Y(),REF_DEN_FIELD);
//  this->tracer_ref.get_density_field_grid(_COUNTS_,REF_DEN_FIELD);
  int nmax=get_max<real_prec>(REF_DEN_FIELD);
  So.message_screen("Maximum number of tracer in cell", nmax);
#endif

  // ********************************************************************************

#ifdef _USE_TOTAL_MASS_IN_CELL_
  vector<real_prec> REF_MASS_FIELD(this->NGRID,0);
  this->tracer_ref.get_density_field_grid(_MASS_, REF_MASS_FIELD);
  real_prec mmax=get_max<real_prec>(REF_MASS_FIELD);
  So.message_screen("Maximum mass of tracer in cell", mmax);
#endif
  // ********************************************************************************

#ifdef _USE_NUMBER_OF_NEIGHBOURS_
  this->tracer_ref.get_neighbour_tracers(this->ncells_info);
#endif
  // ********************************************************************************
#ifdef _GET_DIST_MIN_SEP_REF_
    this->tracer_ref.get_distribution_min_separations(this->tracer_ref.type_of_object, this->ncells_info);
#endif
// ********************************************************************************

#ifdef _GET_DIST_REDUCED_MASS_
  this->tracer_ref.get_distribution_reduced_mass_in_cell();
#endif
  // ********************************************************************************
#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
  this->tracer_ref.get_min_separation_in_cell();
#ifdef _ASSIGN_TO_REFERENCE_
  this->min_halo_separation=this->tracer_ref.min_halo_separation;
#else
  this->min_halo_separation=MIN_SEP_IN_CELLS;
#endif
  real_prec delta_min_sep = (MAX_SEP_IN_CELLS-this->min_halo_separation)/static_cast<real_prec>(N_BINS_MIN_SEP_IN_CELLS);
#endif
  // ********************************************************************************
  // ********************************************************************************
  // ********************************************************************************

#ifdef _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_

#ifdef _USE_VMAX_AS_OBSERVABLE_
#ifdef _FULL_VERBOSE_
  So.message_screen("Allocating", LENGHT_AB* (sizeof(real_prec)+sizeof(ULONG))/(1e9), "Gb for list of Vmax (reference) in theta-bins");
#endif
#elif defined _USE_MASS_AS_OBSERVABLE_
#ifdef _FULL_VERBOSE_
  So.message_screen("Allocating", LENGHT_AB* (sizeof(real_prec)+sizeof(ULONG))/(1e9), "Gb for list of masses (reference) in theta-bins");
#endif
#endif
  this->dm_properties_bins.clear();            // Container of structures meant to keep the properties of the tracers found in a bin of Theta
  this->dm_properties_bins.shrink_to_fit();
  this->dm_properties_bins.resize(LENGHT_AB);
  So.DONE();


#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
  this->dm_properties_for_randoms_bins.clear();            // Container of structures meant to keep the properties of the tracers found in a bin of Theta
  this->dm_properties_for_randoms_bins.shrink_to_fit();
  this->dm_properties_for_randoms_bins.resize(LENGHT_AB);
  So.DONE();
#endif

#endif

#ifndef _ASSIGN_TO_REFERENCE_
#ifdef _FULL_VERBOSE_
  So.message_screen("Allocating", LENGHT_AB_ONE* (sizeof(ULONG))/(1e9), "Gb for probability distribution");
#endif
  this->ABUNDANCE.clear();
  this->ABUNDANCE.shrink_to_fit();
  this->ABUNDANCE.resize(LENGHT_AB_ONE, 0);
  So.DONE();
#endif

  // ********************************************************************************
#ifdef _GET_HALO_MASS_PROPERTIES_CORRELATIONS_
  real_prec mean_hmass=0;
  real_prec mean_ldm=0;
  real_prec m_dm=0;
  real_prec mean_ntr=0;
  real_prec m_ntr=0;
#endif


#ifndef _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_
  So.message_screen("**Measuring n(M|theta) from the reference  using ", this->tracer_ref.NOBJS," objects");
#else
#ifdef _FULL_VERBOSE_
#ifdef _USE_MASS_AS_OBSERVABLE_
    So.message_screen("**Identifying reference masses in bins of {theta}_ref using ", this->tracer_ref.NOBJS," objects");
#elif defined _USE_VMAX_AS_OBSERVABLE_
  So.message_screen("**Identifying reference values of Vmax in bins of {theta}_ref with ", this->tracer_ref.NOBJS," objects from the reference catalog");
#endif
#endif

#endif


#ifdef  _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
    ULONG counter_vmax=0;
    ULONG counter_vmax_r=0;
#endif

  // *********************************************************************************************
  // *********************************************************************************************
    vector<s_cell_info> cell_info_tr(this->NGRID);

 // *********************************************************************************************
 // This can be improved by doing two loops. In the first loop, Abundance is computed in a prallel way,
  // while the structore allocating the thera_bin ofin which the id where the galaxy ig licves is allocated.
  // This can then be used in the second loop over grid, filling the dm_properties
    // *********************************************************************************************
 // Open loop over the number of tracer objects
// Do not parallelize. There are push_backs inside this loop.


//  cout<<RED<<this->delta_dm_aux_mem[1255]<<"   "<<this->delta_dm_aux_mem.size()<<"   "<<this->NGRID<<RESET<<endl;
  cout<<RED<<this->delta_dm_aux_mem[12255]<<"   "<<cwclass_ref.Invariant_TF_II[12255]<<"  "<<cwclass_ref.Invariant_TF_III[12255]<<endl;

  for(ULONG ig = 0; ig< this->tracer_ref.NOBJS ; ++ig)
    {
     real_prec halo_prop=0;
     int I_Y=0;
#ifdef _USE_MASS_AS_OBSERVABLE_
      halo_prop=log10(this->tracer_ref.Halo[ig].mass/this->params._MASS_units());  // Using the logarithm of Mass
#ifndef _ASSIGN_TO_REFERENCE_
      I_Y= get_bin(halo_prop,lm_min,this->params._NMASSbins_mf(),this->tracer_ref.logdeltaM,this->bin_accumulate_borders);
#endif

#elif defined _USE_VMAX_AS_OBSERVABLE_
      halo_prop=log10(this->tracer_ref.Halo[ig].vmax);   //Using the logarithn of Vmax

#ifndef _ASSIGN_TO_REFERENCE_
      I_Y =get_bin(halo_prop, log10(this->params._VMAXmin()),N_x,log10(this->params._VMAXmax()/this->params._VMAXmin())/static_cast<double>(N_x), true);
#endif

#endif
    // *********************************************************************************************
    // Get ID of this reference dark matter tracers
    // *********************************************************************************************
      ULONG ID=this->tracer_ref.Halo[ig].GridID;
    // *********************************************************************************************
    // Get bin of DM ovedensity
    // *********************************************************************************************
      real_prec xdm  = static_cast<real_prec>(this->delta_dm_aux_mem[ID]);
      int I_X  = get_bin(xdm,this->s_mins.prop1,this->NX,this->s_deltas.prop1,this->bin_accumulate_borders);

      int I_CWT=0;
#ifdef _USE_CWC_
      I_CWT=this->cwclass_ref.get_Tclassification(ID);
#endif

      int I_MK=0;
#ifdef _USE_MASS_KNOTS_
      I_MK= (this->cwclass_ref.cwt_used[I_CWT]== I_KNOT ? this->cwclass_ref.SKNOT_M_info[ID]: 0);
#endif

      int I_CWV=0;
#ifdef _USE_CWC_V_
      I_CWV=cwclass_ref.get_Vclassification(ID);
#endif
      int I_VK=0;
#ifdef _USE_VEL_KNOTS_V_
      I_VK= (this->cwclass_ref.cwv_used[I_CWV]== I_KNOT ? this->cwclass_ref.VDISP_KNOT_info[ID]: 0);
#endif

      int I_C1=0;
#ifdef _USE_TOTAL_MASS_IN_CELL_
      I_CV1=get_bin( log10(REF_MASS_FIELD[ID]), this->params._LOGMASSmin(),N_BINS_TOTAL_MASS_IN_CELL,(log10(MAX_TOTAL_MASS_IN_CELL)-this->params._LOGMASSmin())/static_cast<double>(N_BINS_TOTAL_MASS_IN_CELL), true);
#elif defined _USE_INVARIANT_TIDAL_FIELD_II_
      real_prec C1 = this->cwclass_ref.Invariant_TF_II[ID];
      I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#elif defined _USE_DELTA2_
      real_prec C1 = this->cwclass_ref.DELTA2[ID];
      I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#endif

      int I_C2=0;
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
      I_C2=this->tracer_ref.Number_of_neighbours[ig];
      if(I_C2==N_NEIGHBOURS_MAX)
        I_C2=N_NEIGHBOURS_MAX-1;
#elif defined _USE_INVARIANT_TIDAL_FIELD_III_
      real_prec C2 = this->cwclass_ref.Invariant_TF_III[ID];
      I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#elif defined _USE_DELTA3_
      real_prec C2 = this->cwclass_ref.DELTA3[ID];
      I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#endif


      int I_C3=0;
#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
      real_prec min_sep = this->tracer_ref.min_separation_in_cell[ID];
      I_C3 = get_bin(min_sep, 0, N_c, delta_min_sep , this->bin_accumulate_borders);
#elif defined _USE_INVARIANT_TIDAL_FIELD_IV_
      real_prec C3 = this->cwclass_ref.Invariant_TF_IV[ID];
      I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined _USE_TIDAL_ANISOTROPY_
      real_prec C3 = this->cwclass_ref.Tidal_Anisotropy[ID];
      I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined _USE_S2_
      real_prec C3 = this->cwclass_ref.S2[ID];             // s²
      I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#endif

      int I_CV1=0;
#ifdef _USE_TRACERS_IN_CELLS_
      I_CV1=get_bin(static_cast<int>(REF_DEN_FIELD[ID]),0,N_BINS_TRACERS_IN_CELLS,DELTA_TRACERS_IN_CELLS,true);
//      I_CV1=static_cast<int>(REF_DEN_FIELD[ID]);
#elif defined (_USE_INVARIANT_SHEAR_VFIELD_I_)
      real_prec CV1 = invariant_field_I(cwclass_ref.lambda1_vs[ID],cwc_ref.lambda2_vs[ID],cwc_ref.lambda3_vs[ID]); // Not yet assigned to container
      I_CV1= get_bin(CV1, CV1_MIN, N_CV_BIN1,DELTA_CV1, this->bin_accumulate_borders);
#elif defined _USE_NABLA2DELTA_
      real_prec CV1 = this->cwclass_ref.N2D[ID];      // Nabla² ð
      I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1, s_deltas.prop7,this->bin_accumulate_borders);
#endif

      int I_CV2=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
      real_prec CV2 = invariant_field_II(cwclass_ref.lambda1_vs[ID],cwclass_ref.lambda2_vs[ID],cwclass_ref.lambda3_vs[ID]);// Not yet assigned to container
      I_CV2= get_bin(CV2, CV2_MIN, N_CV_BIN2,DELTA_CV2, this->bin_accumulate_borders);
#elif defined _USE_S2DELTA_
      real_prec CV2 = this->cwclass_ref.S2DELTA[ID];         // s²ð
      I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2, s_deltas.prop8,this->bin_accumulate_borders);
#endif

      int I_CV3=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
      real_prec CV3 = invariant_field_III(cwclass_ref.lambda1_vs[ID],cwclass_ref.lambda2_vs[ID],cwclass_ref.lambda3_vs[ID]);// Not yet assigned to container
      I_CV3= get_bin(CV3, CV3_MIN, N_CV_BIN3,DELTA_CV3, this->bin_accumulate_borders);
#elif defined _USE_S3_
      real_prec CV3 = this->cwclass_ref.S3[ID];                                   // s³
      I_CV3= get_bin(CV3, this->s_mins.prop9, N_CV_BIN3, s_deltas.prop9,this->bin_accumulate_borders);
#endif


#ifdef _GET_HALO_MASS_PROPERTIES_CORRELATIONS_
      mean_ldm+=xdm;
      mean_hmass+=halo_prop;
      m_dm+=halo_prop*xdm;
#ifdef _USE_TRACERS_IN_CELLS_
      mean_ntr+=REF_DEN_FIELD[ID];
      m_ntr+=halo_prop*REF_DEN_FIELD[ID];
#endif
#endif

#ifndef _BIN_ACCUMULATE_
      if(xdm <=this->s_maxs.prop1 && xdm >=this->s_mins.prop1)
#endif
#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
        if(C1<= this->s_maxs.prop4 && C1>= this->s_mins.prop4)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
          if(C2<= this->s_maxs.prop5 && C2>= this->s_mins.prop5)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_TIDAL_FIELD_IV_) || defined(_USE_TIDAL_ANISOTROPY_)
            if(C3<= this->s_maxs.prop6 && C3>= this->s_mins.prop6)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
              if(CV1<= this->s_maxs.prop7 && CV1>= this->s_mins.prop7)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
                if(CV2<= this->s_maxs.prop8 && CV2>= this->s_mins.prop8)
#endif
#endif

#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
                  if(CV3<= this->s_maxs.prop9 && CV3>= this->s_mins.prop9)
#endif
                      {
#endif
                      ULONG index_ant=index_11d(I_X, I_CWT, I_MK, I_CWV,I_VK, I_C1, I_C2, I_C3, I_CV1, I_CV2, I_CV3, this->n_cwt, this->n_sknot_massbin, this->n_cwv, this->n_vknot_massbin,N_a, N_b ,N_c, N_v, N_CV_BIN2,N_CV_BIN3);

                      // The container ABUNDANCE will be used used to assign property to mock tracers below the threshold Xthreshold
#ifndef  _ASSIGN_TO_REFERENCE_
                      ULONG index_prop_ab= index_2d(I_Y,index_ant,LENGHT_AB);
                      this->ABUNDANCE[index_prop_ab]++;
#endif

                      // This will be used for mocks and assign_to_ref for X>=Xthres level 1 only
#ifdef _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_
#ifdef _USE_MASS_AS_OBSERVABLE_
                      this->dm_properties_bins[index_ant].masses_bin_properties.push_back(this->tracer_ref.Halo[ig].mass);
                      this->dm_properties_bins[index_ant].used_mass.push_back(false);
#elif defined _USE_VMAX_AS_OBSERVABLE_

#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_  //this is only done for vmax so far, since we shall  use that property as main property
                     // If the ordered index vmax_index (ordered from bottom-to-top in vmax) for each reference tracer is higher than Nran, means that the vmax for this tracer wll be assigned to DMparticles at assignment campaign
                      // while the lowest this->tracer.Ntracers_ran vmax available will be assigned to tracers generated with random coordinates.
#ifdef _BOTTOM_RANDOM_
                      if(this->tracer_ref.Halo[ig].vmax_index >= this->tracer.Ntracers_ran)// This "if" guarranties that the number of vmax contained here are equal to the number of dm- assigned tracers.
#elif defined (_TOP_RANDOM_)                                                                                     // In the Bottom-random case, the dm particles get the highst Vmax values
                      if(this->tracer_ref.Halo[ig].vmax_index < this->tracer.Ntracers_dm)  // In the Top_Random case , dm particles get the lowest vmax values
#endif
                       {
                           counter_vmax++;
#endif
                           this->dm_properties_bins[index_ant].masses_bin_properties.push_back(this->tracer_ref.Halo[ig].vmax);
                           this->dm_properties_bins[index_ant].used_mass.push_back(false);
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#if defined _TOP_RANDOM_ || defined (_BOTTOM_RANDOM_)
                        }
                       else  // if this->tracer_ref.Halo[ig].vmax_index < this->tracer.Ntracers_ran, we fill withthe value sof Vmax the vectors of structures for the randoms
                        {
                          this->dm_properties_for_randoms_bins[index_ant].masses_bin_properties.push_back(this->tracer_ref.Halo[ig].vmax);
                          this->dm_properties_for_randoms_bins[index_ant].used_mass.push_back(false);
                           counter_vmax_r++;
                        }
#endif // end  #if defined _TOP_RANDOM_ || defined (_BOTTOM_RANDOM_)
#endif // end ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#endif //  end defined _USE_VMAX_AS_OBSERVABLE_
#endif // end _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_

#ifndef _BIN_ACCUMULATE_
                      }
#endif

    }

  this->So.DONE();
//  So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"tracer_ref.Halo.clear() is disabled in order to compute correlation between masses");
 // this->tracer_ref.Halo.clear();
 // this->tracer_ref.Halo.shrink_to_fit();


#ifdef _GET_HALO_MASS_PROPERTIES_CORRELATIONS_
  mean_ldm/=static_cast<real_prec>(this->tracer_ref.NOBJS);
  mean_hmass/=static_cast<real_prec>(this->tracer_ref.NOBJS);
  m_dm/=static_cast<real_prec>(this->tracer_ref.NOBJS);
#ifdef _USE_MASS_AS_OBSERVABLE_
  So.message_screen("Correlation log(M) - DM density =", sqrt(fabs(m_dm-mean_hmass*mean_ldm)));
#elif defined _USE_VMAX_AS_OBSERVABLE_
#ifdef _FULL_VERBOSE_
  So.message_screen("Correlation log(Vmax) - DM density =", sqrt(fabs(m_dm-mean_hmass*mean_ldm)));
#endif
#endif

#ifdef _USE_TRACERS_IN_CELLS_
  mean_ntr/=static_cast<real_prec>(this->tracer_ref.NOBJS);
  m_ntr/=static_cast<real_prec>(this->tracer_ref.NOBJS);
#ifdef _USE_MASS_AS_OBSERVABLE_

#ifdef _FULL_VERBOSE_
  So.message_screen("Halo-mass - N_tracers in cells  =", sqrt(fabs(m_ntr-mean_hmass*mean_ntr)));
#endif
#elif defined _USE_VMAX_AS_OBSERVABLE_
#ifdef _FULL_VERBOSE_
  So.message_screen("Halo-Vmax - N_tracers in cells =", sqrt(fabs(m_ntr-mean_hmass*mean_ldm)));
  cout<<endl;
#endif
#endif
#endif
#endif




#ifdef _USE_TRACERS_IN_CELLS_
  REF_DEN_FIELD.clear();
  REF_DEN_FIELD.shrink_to_fit();
#endif




#ifdef _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_
  real_prec aux_n=-100;
  real_prec aux_m;
  for(int ih=0;ih<LENGHT_AB ;++ih)
     {
       aux_m=max(static_cast<real_prec>(this->dm_properties_bins[ih].masses_bin_properties.size()), aux_n);
       aux_n=aux_m;
      }
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#ifdef _FULL_VERBOSE_
#ifdef _BOTTOM_RANDOM_
  So.message_screen("Maximum number of tracers with Vmax > Vmax_threshold in a bin of DM properties {Theta} = ", static_cast<real_prec>(aux_m));
#else
  So.message_screen("Maximum number of tracers with Vmax < Vmax_threshold in a bin of DM properties {Theta} = ", static_cast<real_prec>(aux_m));
#endif  
#endif  


  aux_n=-100;
  aux_m=0;
  for(int ih=0;ih<LENGHT_AB ;++ih)
     {
       aux_m=max(static_cast<real_prec>(this->dm_properties_for_randoms_bins[ih].masses_bin_properties.size()), aux_n);
       aux_n=aux_m;
      }
#ifdef _FULL_VERBOSE_
#ifdef _BOTTOM_RANDOM_
  So.message_screen("Maximum number of tracers with Vmax <= Vmax_threshold in a bin of DM properties {Theta} = ", static_cast<real_prec>(aux_m));
#else
  So.message_screen("Maximum number of tracers with Vmax >= Vmax_threshold in a bin of DM properties {Theta} = ", static_cast<real_prec>(aux_m));
#endif  
#endif  

#else
#ifdef _FULL_VERBOSE_
  So.message_screen("Maximum number of tracers in a bin of DM properties {Theta} = ", static_cast<real_prec>(aux_m));
  cout<<endl;
#endif
#endif


  ULONG aux_a=0;
#ifdef _FULL_VERBOSE_
  So.message_screen("Checking number of tracers in theta-containers:");
#endif
#pragma omp parallel for reduction(+:aux_a)
  for(ULONG ih=0;ih<LENGHT_AB ;++ih)
    aux_a+=this->dm_properties_bins[ih].used_mass.size();// using masses_bin_properties leads to the same result.


#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
  So.message_screen("Number of tracers_dm in theta-containers =", aux_a);
  So.message_screen("check =", counter_vmax);
  So.message_screen("Expected =", this->tracer.Ntracers_dm);
#else
#ifdef _FULL_VERBOSE_
  So.message_screen("Number of tracers in theta-containers in =", aux_a);
#endif
#endif

#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
  //Here we sum the tracers with Vmax below the threshold identified below which random tracers will take their vmax
  ULONG aux_b=0;
#pragma omp parallel for reduction(+:aux_b)
  for(ULONG ih=0;ih<LENGHT_AB ;++ih)
    aux_b+=this->dm_properties_for_randoms_bins[ih].used_mass.size();

  So.message_screen("Number of tracers_random in theta-containers =", aux_b);
  So.message_screen("check =", counter_vmax_r);
  So.message_screen("Expected =", this->tracer.Ntracers_ran);

  aux_a+=aux_b;
#endif


  if(aux_a<this->tracer_ref.NOBJS)
    {
      So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"The number of tracers counted in vec<dm_props> dm_properties is smaller than the input. Perhaps the range for some of the DMF properties *must* be increased to include all cells. CosmiCatlas stops here");
      cout<<aux_a<<"  "<<this->tracer_ref.NOBJS<<endl;
      exit(0);
    }
  So.DONE();

#ifndef _ASSIGN_TO_REFERENCE_
   aux_a=0;
#ifdef _FULL_VERBOSE_
  So.message_screen("Checking number of tracers in abundance:");
#endif
#pragma omp parallel for reduction(+:aux_a)
  for(ULONG ih=0;ih<this->ABUNDANCE.size() ;++ih)
     aux_a+=ABUNDANCE[ih];
  if(aux_a<this->tracer_ref.NOBJS)
    {
      So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"The number of tracers counted in ABUNDANCE(Y,X) is smaller than the input. Perhaps the range for some of the DMF properties *must* be increased to include all cells. CosmiCatlas stops here");
      cout<<aux_a<<"  "<<this->tracer_ref.NOBJS<<endl;
      exit(0);
    }
  So.DONE();
#endif


#endif


#ifndef _ASSIGN_TO_REFERENCE_
#ifdef _FULL_VERBOSE_
  So.message_screen("Normalizing joint probability distribution");
#endif
  this->ABUNDANCE_normalized.clear();
  this->ABUNDANCE_normalized.shrink_to_fit();
  this->ABUNDANCE_normalized.resize(LENGHT_AB_ONE,0);


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG idm = 0; idm < LENGHT_AB; ++idm)
  {
      long aux_a=-10000;
      long aux_b;

      for(int tr_j = 0; tr_j < N_x; ++tr_j)
        {
          ULONG indexa=index_2d(tr_j,idm,LENGHT_AB);
          long AUX=this->ABUNDANCE[indexa];
          aux_b=max(aux_a, AUX);
          aux_a=aux_b;
        }
      for(int tr_j = 0; tr_j < N_x; ++tr_j)
        {
          ULONG indexa=index_2d(tr_j,idm,LENGHT_AB);
          long AUX=this->ABUNDANCE[indexa];
          this->ABUNDANCE_normalized[indexa]=(aux_b==0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
        }
    }

/*  aux_a=0;
  So.message_screen("Checking sum in abundance_normalized:");
#pragma omp parallel for reduction(+:aux_a)
  for(ULONG ih=0;ih<this->ABUNDANCE.size() ;++ih)
     aux_a+=ABUNDANCE_normalized[ih];

*/


#ifdef _FULL_VERBOSE_
  So.DONE();
  So.message_screen("Freeing memory");
#endif
  this->ABUNDANCE.clear();
  this->ABUNDANCE.shrink_to_fit();
  So.DONE();

#endif




#ifndef _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_
  //* Final check to verify that all tracers have been counted in
  ULONG ncells_used = static_cast<ULONG>(get_nobjects(this->ABUNDANCE));

  if(ncells_used<this->tracer_ref.NOBJS)
    {
      So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"The number of tracers counted in ABUNDANCE(Y,X) is smaller than the input. Perhaps the range for some of the DMF properties *must* be increased to include all cells. CosmiCatlas stops here");
      cout<<"Number of cells used = "<<ncells_used<<"  Number expected = "<<this->tracer_ref.NOBJS<<endl;
    }


   So.message_screen("Assigning memmory space for conditional probability function ABUNDANCE_normalized");
   this->ABUNDANCE_normalized.shrink_to_fit();
   this->ABUNDANCE_normalized.clear();
   this->ABUNDANCE_normalized.resize(LENGHT_AB_ONE, 0);

    ULONG LENGHT_AC  = N_CV_BIN1*N_CV_BIN2*N_CV_BIN3*N_a* N_b* N_c * this->n_sknot_massbin * this->n_cwt*this->n_vknot_massbin*this->n_cwv*this->NX;
    this->NCELLSperDMBIN.clear();
    this->NCELLSperDMBIN.shrink_to_fit();
    this->NCELLSperDMBIN.resize(LENGHT_AC, 0);
    So.DONE();

    So.message_screen("**Normalizing number counts and marginalizing with respect to Y-bins:");


#ifdef _USE_OMP_
#pragma omp parallel for collapse(11)
#endif
    for(int i=0;i < this->NX;++i)// loop sobre bines de dark matter
      for(int sua = 0; sua < this->n_cwt; ++sua)
        for(int k = 0; k < this->n_sknot_massbin ;  ++k)
          for(int vua = 0; vua < this->n_cwv; ++vua)
            for(int vk = 0; vk < this->n_vknot_massbin ;  ++vk)
              for(int l1 = 0; l1 < N_a; ++l1)
                for(int l2 = 0; l2 < N_b ; ++l2)
                  for(int l3 = 0; l3 < N_c; ++l3)
                    for(int lv1 = 0; lv1 < N_CV_BIN1; ++lv1)
                      for(int lv2 = 0; lv2 < N_CV_BIN2; ++lv2)
                        for(int lv3 = 0; lv3 < N_CV_BIN3; ++lv3)
                           {
#ifndef _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_
                                ULONG index_l=index_11d(i,sua,k,vua,vk,l1,l2,l3,lv1,lv2,lv3,this->n_cwt,this->n_sknot_massbin,this->n_cwv,this->n_vknot_massbin,N_a,N_b,N_c, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
#endif

                                long aux_a=-10000;
                                long aux_b;
                                for(int tr_j = 0; tr_j < this->params._NMASSbins_mf(); ++tr_j)
                                  {
                                    ULONG indexa=index_12d(tr_j,i,sua,k,vua,vk,l1,l2,l3,lv1, lv2, lv3,this->NX,this->n_cwt, this->n_sknot_massbin,this->n_cwv, this->n_vknot_massbin, N_a,N_b,N_c, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
                                    long AUX=this->ABUNDANCE[indexa];
                                    aux_b=max(aux_a, AUX);
                                    aux_a=aux_b;

#ifndef _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_
                                    ULONG index_h=index_12d(tr_j,i,sua,k,vua,vk,l1,l2,l3,lv1,lv2,lv3,this->NX,this->n_cwt,this->n_sknot_massbin,this->n_cwv,this->n_vknot_massbin,N_a, N_b,N_c, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
#pragma omp atomic update
                                    this->NCELLSperDMBIN[index_l]+=this->ABUNDANCE[index_h];
#endif
                                  }
                                for(int tr_j = 0; tr_j < this->params._NMASSbins_mf(); ++tr_j)
                                  {
                                    ULONG indexa=index_12d(tr_j,i,sua,k,vua,vk,l1,l2,l3,lv1, lv2, lv3,this->NX,this->n_cwt, this->n_sknot_massbin,this->n_cwv, this->n_vknot_massbin, N_a,N_b,N_c, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
                                    long AUX=this->ABUNDANCE[indexa];
                                    this->ABUNDANCE_normalized[indexa]=(aux_b==0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
                                  }
                           }


    this->So.DONE();

#ifdef _FULL_VERBOSE_
    So.message_screen("Check on number of bins used...");
#endif
    real_prec aux_h=0;

#pragma omp parallel for reduction(+:aux_h)
    for(int ih=0;ih<this->ABUNDANCE_normalized.size();++ih)
      aux_h+=this->ABUNDANCE_normalized[ih];
    if(aux_h<=0)
      {
        So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Normalized Joint probability for number counts ill-defined. CosmicAtlas stops here.");
        exit(0);
      }
    else
      this->So.DONE();


    So.message_screen("Freeing memmory");
//    this->ABUNDANCE.clear();
//    this->ABUNDANCE.shrink_to_fit();
    So.DONE();

    // We do not need here to write these files as they will be inmediatly used in the assign_tracer_mass
    //  this->File.write_array(this->Output_directory+"Bam_Abundance", this->ABUNDANCE);
    // this->File.write_array(this->Output_directory+"Bam_Abundance_Normalized", this->ABUNDANCE_normalized);

#endif



}

//end class member function get_X_function()


//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################


// This function is explicitely meant to be applid to assign masses after having assigne vmax

void Bam::get_X_function_complement(string h_property)
{

    So.enter(__PRETTY_FUNCTION__);


#ifdef _FULL_VERBOSE_
   if(h_property==_MASS_)
     {
    So.message_screen("***Measuring conditional MASS Function from reference tracer catalog*****");
    So.message_screen("****as a function of DM properties of the reference DM density field*****");
     }
    else if(h_property==_RS_)
   {
      So.message_screen("***Measuring conditional Rs Function from reference tracer catalog*****");
      So.message_screen("****as a function of M and Vmax *****");
    }
   else if(h_property==_SPIN_)
  {
     So.message_screen("***Measuring conditional Spin Function from reference tracer catalog*****");
     So.message_screen("****as a function of M and Vmax *****");
   }
#endif


  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);

#ifdef _FULL_VERBOSE_
#if defined (test_vmax_mass) && defined (_add_dm_density_)
#ifdef _add_Xweb_
  So.message_warning("Virial masses assigned from P(M|Vmax,theta) scaling relation");
#else
  if(h_property==_MASS_)
  So.message_warning("Virial masses assigned from P(Mvir|Vmax,delta) scaling relation");
  else if (h_property==_RS_)
      So.message_warning("RS assigned from P(RS|Vmax,Mvir) scaling relation");
  else if (h_property==_SPIN_)
      So.message_warning("Spin assigned from P(RS|Vmax,Mvir) scaling relation");

#endif
#else
  So.message_warning("Virial masses assigned from P(M|Vmax) scaling relation");
#endif
  cout<<endl;
#endif





//#ifndef test_vmax_mass
  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************

  ULONG N_a=1;//N_C_BIN1;  // This function will not be using Theta properties, meant for abudnance P(X|y, delta)

  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************
  ULONG N_b=1;//N_C_BIN2;

  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************
  ULONG N_c = 1; //N_C_BIN3; // We shall not use min_sep for post-mass assignment
  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************

  ULONG N_v=1;// N_CV_BIN1;
  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************
//#endif




 int N_x= N_CV_BIN2;
#ifdef _ASSIGN_MASS_POST_
  N_x = N_VMAX_BINS;
#endif



  ULONG Ntot=N_x;

  if(h_property==_MASS_)
    {    // we will only use dm for mass assignment. For rs, we use M and Vmax for the time being
#ifdef _add_dm_density_
     Ntot*=this->NX;
#ifdef _add_Xweb_
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
     Ntot*=N_C_BIN1;  // We are assuming that Xweb
#endif
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
      Ntot*=N_C_BIN2;  // We are assuming that Xweb
#endif
#endif
#endif
   }
   else if(h_property==_RS_ || h_property==_SPIN_) // If we want to get Rs, we use P(Rs|M,Vmax) so we add a dimention fo M
     Ntot*=N_MASS_BINS;


#ifndef test_vmax_mass
  ULONG LENGHT_AB=N_v*N_x*N_CV_BIN3* N_a* N_b* N_c;
  LENGHT_AB*=this->n_sknot_massbin * this->n_cwt * this->n_vknot_massbin * this->n_cwv*this->NX;
#else
  ULONG LENGHT_AB=Ntot;
#endif


  // ********************************************************************************
  ULONG LENGHT_AC=LENGHT_AB*this->params._NMASSbins(); // Here we explicitly say that the variable to be reconstructed (Rs or Mass) has bin sizes from para file his->params._NMASSbins
  this->ABUNDANCE.clear();
  this->ABUNDANCE.shrink_to_fit();
  this->ABUNDANCE.resize(LENGHT_AC,0);

  // ********************************************************************************
  this->get_new_min_max_properties();

  // ********************************************************************************

#ifdef _USE_TOTAL_MASS_IN_CELL_
  vector<real_prec> REF_MASS_FIELD(this->NGRID,0);
  this->tracer_ref.get_density_field_grid(_MASS_, REF_MASS_FIELD);
#endif
  // ********************************************************************************

#ifndef test_vmax_mass
#ifdef _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_
  So.message_screen("Allocating", LENGHT_AB* (sizeof(real_prec)+sizeof(ULONG))/(1e9), "Gb for list of masses in theta-bins");
  this->dm_properties_bins.clear();
  this->dm_properties_bins.shrink_to_fit();   //container of structures, meant to keep properties of halos living in a given bin of Theta
  this->dm_properties_bins.resize(LENGHT_AB);   //vector de estructiras para guarar las masas que caen en un bin de {Theta}
  So.DONE();
#endif
#endif

  // ********************************************************************************

#ifndef test_vmax_mass
#ifndef _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_
  So.message_screen("**Measuring n(M|theta) from the reference  using ", this->tracer_ref.NOBJS," objects");
#else
  So.message_screen("**Identifying reference Masses in bins of {theta}_ref using ", this->tracer_ref.NOBJS," objects");
#endif
#endif

#ifdef _USE_OMP_
#pragma omp parallel
  {

      vector<real_prec>abundance_parallel(this->ABUNDANCE.size(),0);
#ifdef _USE_OMP_
#pragma omp for nowait
#endif
  for(ULONG ig = 0; ig< this->tracer_ref.NOBJS ; ++ig)
    {

      // Get ID of this tracer
#if !defined test_vmax_mass || defined (_add_dm_density_)
      ULONG ID=this->tracer_ref.Halo[ig].GridID;
#endif

      int I_Y=0;//This index takes the place of the derived variable to assign (e.g. Mvir, Rs)
      if(h_property==_MASS_)
        I_Y  = get_bin(log10(this->tracer_ref.Halo[ig].mass),this->params._LOGMASSmin(),this->params._NMASSbins(),this->tracer_ref.logdeltaM,this->bin_accumulate_borders);
      else if(h_property==_RS_)
        I_Y  = get_bin(log10(this->tracer_ref.Halo[ig].rs),log10(this->params._RSmin()),this->params._NMASSbins(),this->tracer_ref.logdeltaRS,this->bin_accumulate_borders);
      else if(h_property==_SPIN_)
       I_Y  = get_bin(log10(this->tracer_ref.Halo[ig].spin),log10(this->params._SPINmin()),this->params._NMASSbins(),this->tracer_ref.logdeltaSPIN,this->bin_accumulate_borders);
          // Get bin of DM ovedensity
      int I_X=0; // We use this bin as for I_X ->delta for P(M|delta,Vmax) , or I_X -> Mass  for P(Rs|M,Vmax). Vmax will take the index I_CV2. If add_dm_density undef, I_X=0 for P(M|Vmax)

#ifdef _add_dm_density_
      if(h_property==_MASS_) // Using DELTA only allowed so far when assigning mass. FOr Rs we use so fat P(Rs|M,Vmax)
       {
         real_prec xdm  = this->delta_dm_aux_mem[ID]; // this and its I-or T-web  was computed in member function  get_X_function()
         I_X  = get_bin(xdm,this->s_mins.prop1,this->NX,this->s_deltas.prop1,this->bin_accumulate_borders);
       }
      else
#endif
    // If assignment is for RS, IX takes the place of halo mass
      if(h_property==_RS_ || h_property==_SPIN_ ) // Using DELTA only allowed so far when assigning mass. FOr Rs and Spin we use so far P(Rs|M,Vmax) or P(S|M,Vmax)
        I_X =  get_bin(log10(this->tracer_ref.Halo[ig].mass),this->params._LOGMASSmin(),N_MASS_BINS, (this->params._LOGMASSmax()-this->params._LOGMASSmin())/static_cast<double>(N_MASS_BINS)  ,this->bin_accumulate_borders);

       ULONG I_CV2=0;
#ifdef test_vmax_mass
#ifdef _ASSIGN_MASS_POST_
       I_CV2=get_bin(log10(this->tracer_ref.Halo[ig].vmax), log10(this->params._VMAXmin()),N_VMAX_BINS,log10(this->params._VMAXmax()/this->params._VMAXmin())/static_cast<double>(N_VMAX_BINS), true);
#endif
#endif

#ifdef _add_dm_density_
#ifndef _BIN_ACCUMULATE_
       if(xdm <=this->s_maxs.prop1 && xdm >=this->s_mins.prop1)
          {
#endif
#endif
#ifndef test_vmax_mass
                      ULONG index_prop=index_11d(I_X, I_CWT, I_MK, I_CWV,I_VK, I_C1, I_C2, I_C3, I_CV1, I_CV2, I_CV3,this->n_cwt, this->n_sknot_massbin, this->n_cwv, this->n_vknot_massbin,N_a,N_b,N_c, N_v, N_x,N_CV_BIN3);
                      this->dm_properties_bins[index_prop].masses_bin_properties.push_back(this->tracer_ref.Halo[ig].mass);
                      this->dm_properties_bins[index_prop].used_mass.push_back(false);
                      ULONG index_prop_b=index_2d(I_Y,index_prop,LENGHT_AB);
#endif
#ifdef test_vmax_mass
                    ULONG index_dm=0;
                    if(h_property==_MASS_)
                      {
#ifdef _add_dm_density_
                        index_dm=index_2d(I_CV2,I_X,this->NX);
#else
                        index_dm=I_CV2;
#endif
#endif
                     }
                   else if(h_property==_RS_ || h_property==_SPIN_)
                      index_dm=index_2d(I_CV2,I_X,N_MASS_BINS);

                   ULONG index_prop_b=index_2d(I_Y,index_dm,Ntot); // I_Y is the bin of mass
#ifdef _USE_OMP_
                   abundance_parallel[index_prop_b]++;
#else
                   this->ABUNDANCE[index_prop_b]++;
#endif

#ifndef _BIN_ACCUMULATE_
                      }
#endif
    }


#pragma omp critical
  {
    for(ULONG i=0;i <this->ABUNDANCE.size();++i)
       this->ABUNDANCE[i]+=abundance_parallel[i];
  }

  }// close parallel region
#endif

  this->So.DONE();



  ULONG aux_a=0;
#ifdef _FULL_VERBOSE_
  So.message_screen("Checking number of tracers in abundance:");
#endif
#pragma omp parallel for reduction(+:aux_a)
  for(ULONG ih=0;ih<this->ABUNDANCE.size() ;++ih)
     aux_a+=ABUNDANCE[ih];
  if(aux_a<this->tracer_ref.NOBJS)
    {
      So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"The number of tracers counted in ABUNDANCE(Y,X) is smaller than the input. Perhaps the range for some of the DMF properties *must* be increased to include all cells. CosmiCatlas stops here");
      cout<<aux_a<<"  "<<this->tracer_ref.NOBJS<<endl;
      exit(0);
    }
  So.DONE();



 this->delta_dm_aux_mem.clear();delta_dm_aux_mem.shrink_to_fit();

  #if defined (_USE_RS_AS_DERIVED_OBSERVABLE_) && !defined (_USE_SPIN_AS_DERIVED_OBSERVABLE_)
  if(h_property == _RS_) // The tracer info canot be cleared until we have assigned Rs or the last variable
   {
    So.message_screen("**Freeing memmory from tracer_ref and auxiliary containers, line", __LINE__);
    this->tracer_ref.Halo.clear();
    this->tracer_ref.Halo.shrink_to_fit();
   this->So.DONE();
}
#endif

#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
  if(h_property == _SPIN_) // The tracer info canot be cleared until we have assigned Rs or the last variable
   {
#ifdef _FULL_VERBOSE_
    So.message_screen("**Freeing memmory from tracer_ref and auxiliary containers, line", __LINE__);
#endif
    this->tracer_ref.Halo.clear();
    this->tracer_ref.Halo.shrink_to_fit();
   this->So.DONE();
}
#endif


 this->ABUNDANCE_normalized.clear();
 this->ABUNDANCE_normalized.shrink_to_fit();
 this->ABUNDANCE_normalized.resize(LENGHT_AC,0);
#ifdef _FULL_VERBOSE_
  So.message_screen("Normalizing joint probability distribution");
#endif

#ifdef _USE_OMP_
#pragma omp parallel num_threads(NTHREADS)
 {
#pragma omp parallel for
#endif
    for(int ULONG idm = 0; idm < LENGHT_AB; ++idm)
 {
     long aux_a=-10000;
     long aux_b;

     for(int tr_j = 0; tr_j < this->params._NMASSbins(); ++tr_j)
       {
         ULONG indexa=index_2d(tr_j,idm,LENGHT_AB);
         long AUX=this->ABUNDANCE[indexa];
         aux_b=max(aux_a, AUX);
         aux_a=aux_b;
       }
     for(int tr_j = 0; tr_j < this->params._NMASSbins(); ++tr_j)
       {
         ULONG indexa=index_2d(tr_j,idm,LENGHT_AB);
         long AUX=this->ABUNDANCE[indexa];
         this->ABUNDANCE_normalized[indexa]=(aux_b==0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
       }
   }
#ifdef _USE_OMP_
}
#endif

 So.DONE();

 this->ABUNDANCE.clear();
 this->ABUNDANCE.shrink_to_fit();

#ifdef _FULL_VERBOSE_
 cout<<endl;
#endif


#ifndef test_vmax_mass
#ifdef _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_
  real_prec aux_n=-100;
  real_prec aux_m;
  for(int ih=0;ih<LENGHT_AB ;++ih)
     {
       aux_m=max(static_cast<real_prec>(this->dm_properties_bins[ih].masses_bin_properties.size()), aux_n);
       aux_n=aux_m;
      }
  So.message_screen("Maximum number of tracers in a bin of DM properties {Theta} = ", static_cast<real_prec>(aux_m));
  cout<<endl;

  real_prec aux_a;
  So.message_screen("Checking number of tracers in abundance:");
#pragma omp parallel for reduction(+:aux_a)
  for(int ih=0;ih<LENGHT_AB ;++ih)
    aux_a+=this->dm_properties_bins[ih].masses_bin_properties.size();

  if(aux_a<this->tracer_ref.NOBJS)
    {
      So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"The number of tracers counted in ABUNDANCE(Y,X) is smaller than the input. Perhaps the range for some of the DMF properties *must* be increased to include all cells. CosmiCatlas stops here.");
      cout<<aux_a<<"  "<<this->tracer_ref.NOBJS<<endl;
      exit(0);
    }
  So.DONE();
#endif



#endif

}



#endif

//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################

#ifdef MOCK_MODE
 void Bam::get_BIAS(string property)
 {


     So.enter(__PRETTY_FUNCTION__);

    int NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);


     //IOn tbis function the TR and

   real_prec num_in_log_x = true==this->Convert_Density_to_Delta_X ? NUM_IN_LOG: 0.;
   real_prec num_in_log_y = true==this->Convert_Density_to_Delta_Y ? NUM_IN_LOG: 0.;

#ifdef _DO_BAM_CALIBRATION_

#ifdef _FULL_VERBOSE_
   cout<<endl;
   if(this->iteration_ini==0)
     {
       So.message_screen("Measuring Bias: ");
       cout<<endl;
    }
#endif

   if(this->Scale_Y=="log")
     {

       if(this->step==0) // Do this only in the first step, sicne in the iteration process does not change the tracer density field.
        {

#ifdef _FULL_VERBOSE_
       So.message_screen("Transforming delta_Y->log10(num_in_log+delta_Y_ref). Line ", __LINE__);
#endif

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
        for(ULONG i = 0; i < this->NGRID; ++i)
          this->delta_Y[i] = this->delta_Y[i]<-1 ? 0: log10(num_in_log_y+static_cast<real_prec>(this->delta_Y[i]));
          }
      }

   if(this->Scale_X=="log")
     {
#ifdef _FULL_VERBOSE_
       So.message_screen("Transforming delta_X->log10(2+delta_ref). Line ", __LINE__);
#endif


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(ULONG i = 0;i < this->NGRID ;++i)  //TRANSFORM DELTA TO LOG10(NUM_IN_LOG + DELTA). tHIS IS DONE SINCE THE KONVOLUITION ALWAYS GOVES DELTAS, WE HAVE TO TRANSFORMTO LOG(1+DELTA)
	 this->delta_X[i] = this->delta_X[i]<-1 ?  0 : log10(num_in_log_x + static_cast<real_prec>(this->delta_X[i]));
       So.DONE();
     }

#endif  //end if DO_BAM_CALIBRATION

   this->new_nbins_x = this->NX;


   // ******************************************************************************
   // For the mocks we convert always X to delta, so no question mark here. xmin, xmax are not updated
   //  this->new_nbins_x = this->NX;

   if(_COUNTS_==property)
     {
       this->Ymin=NMIN_Y_ONECELL;
       this->Ymax=this->nmax_y_onecell;
       this->new_nbins_y = this->Ymax+1; // One bin per particle, plus 0
     }
   else
     {
       if(this->Scale_Y=="log")
         {
           this->Ymin=this->ldelta_Y_min;
           this->Ymax=this->ldelta_Y_max;

         }
       else
         {
           this->Ymin=this->delta_Y_min;
           this->Ymax=this->delta_Y_max;
         }
       this->new_nbins_y = this->NY;
     }


   this->DELTAX=(this->Xmax-this->Xmin)/static_cast<real_prec>(this->NX);
   this->DELTAY=(this->Ymax-this->Ymin)/static_cast<real_prec>(this->new_nbins_y);

   // ******************************************************************************
   // ************************RESIZE VECTORS FOR HISTOGRAMS*************************

   int NX_NEW=1;
#ifdef _USE_DM_IN_BAM_
   NX_NEW=this->NX;
#endif

   ULONG LENGHT_BIAS_NCOUNT=N_CV_BIN1 * N_CV_BIN2* N_CV_BIN3* N_C_BIN1 * N_C_BIN2* N_C_BIN3* this->n_sknot_massbin * this->n_cwt * this->n_vknot_massbin * this->n_cwv* NX_NEW * this->new_nbins_y;
   ULONG LENGHT_BIAS_NCOUNT_aux=N_CV_BIN1 * N_CV_BIN2* N_CV_BIN3* N_C_BIN1 * N_C_BIN2* N_C_BIN3* this->n_sknot_massbin * this->n_cwt * this->n_vknot_massbin * this->n_cwv* NX_NEW;



   int count_arrays=1;
#ifdef _USE_MASS_TRACERS_
   count_arrays++;
#endif



#if !defined (_ASSIGN_TO_REFERENCE_)

#ifdef _GET_BAM_REALIZATIONS_

#ifdef _FULL_VERBOSE_
   So.message_screen("Allocating", LENGHT_BIAS_NCOUNT*count_arrays*(sizeof(real_prec)+sizeof(ULONG))/(1e9), "Gb for bias");
#endif

#endif


// This section will be done in parallel only if we have more than two threads

  this->BIAS_NCOUNTS.clear();
  this->BIAS_NCOUNTS.shrink_to_fit();
  this->BIAS_NCOUNTS.resize(LENGHT_BIAS_NCOUNT, 0);
  this->BIAS_NCOUNTS_normalized.clear();
  this->BIAS_NCOUNTS_normalized.shrink_to_fit();
  this->BIAS_NCOUNTS_normalized.resize(LENGHT_BIAS_NCOUNT, 0);

#endif


   // ******************************************************************************

#ifdef _GET_BAM_REALIZATIONS_

#ifndef _ONLY_POST_PROC_ // when post_proc is defiened, we do not need the bias for the numbr counts is already generated

#ifdef _FULL_VERBOSE_
   So.message_screen("Reading BIAS from calibration");
#endif


#ifndef __ASSIGN_TO_REFERENCE_
    this->File.read_array(this->params._Input_Directory_BIAS_KERNEL()+"Bam_Bias.dat", this->BIAS_NCOUNTS);
    this->File.read_array(this->params._Input_Directory_BIAS_KERNEL()+"Bam_Bias_Normalized.dat", this->BIAS_NCOUNTS_normalized);
    ULONG ncounts_aux = static_cast<ULONG>(get_nobjects(this->BIAS_NCOUNTS));
#ifdef _FULL_VERBOSE_
    So.message_screen("Number of cells in BIAS = ", ncounts_aux);
#endif
    if(ncounts_aux!=this->NGRID)
       So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Missing cells in get_mock_grid(). Perhaps density range is not wide enough to contain all cells");
   else
    So.DONE();
   real_prec Kd = static_cast<real_prec>(get_nobjects(this->BIAS_NCOUNTS_normalized));
   if(Kd<=0)
       So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Conditional probablity ill-defined.");


#endif  // end of ifndef _MASS_ASSIGNMENT_TO_REFERENCE_

#endif  // endif for  ifndef _ONLY_POST_PROC_


#else  // else for ifdef _GET_BAM_REALIZATIONS_. The following lines are then allowed for CALIBRATION

   if(this->iteration_ini>0 && true==this->use_iteration_ini) //NOT WORKING YET
     {
       So.message_screen("Reading BIAS from iteration", this->iteration_ini);
       this->File.read_array(this->Output_directory+"Bam_Bias.dat", this->BIAS_NCOUNTS);
       this->File.read_array(this->Output_directory+"Bam_Bias_Normalized.dat", this->BIAS_NCOUNTS_normalized);
       this->use_iteration_ini=false; //set it to false such that we will not read this again,
     }
   else
     {
       // ********************************************************************************
       // ********************************************************************************
       this->get_new_min_max_properties();
       // The outputs of this method (private variables) are also used in get_mock
       // ********************************************************************************
       // ********************************************************************************


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(ULONG ig = 0; ig< this->NGRID ; ++ig)
         {
           // Get number counts

#ifdef _DISPLACEMENTS_
           real_prec halo_prop = this->delta_Y[ig]-this->patchy.Displacement[i]; //This line is = Reference Displacement - The displacement computed from patchy
#else
           real_prec halo_prop = this->delta_Y[ig];
#endif

           int I_Y=0;
           if(_COUNTS_== property)
             I_Y=static_cast<int>(halo_prop);
           else
             I_Y=get_bin(halo_prop,this->s_mins.prop0, this->new_nbins_y, this->s_deltas.prop0,this->bin_accumulate_borders);

#ifdef _USE_MASS_FIELD_   //get bin of tracer mass
           real_prec halo_mass_prop=this->delta_Y_MASS[ig];
           int I_Y_MASS=get_bin(halo_mass_prop,this->s_mins.prop0_mass, this->NY_MASS, this->s_deltas.prop0_mass,this->bin_accumulate_borders);
#endif


           int I_X=0;
#ifdef _USE_DM_IN_BAM_
           real_prec xdm    = this->delta_X[ig];
           I_X  = ((this->iMAS_X == 0  && false==this->Convert_Density_to_Delta_X) ? static_cast<int>(this->delta_X[ig]) : get_bin(xdm,this->s_mins.prop1,this->NX,this->s_deltas.prop1,this->bin_accumulate_borders));
#endif

           int I_CWT=0;
#ifdef _USE_CWC_
           I_CWT=this->cwclass.get_Tclassification(ig);
#endif

           int I_MK=0;
#ifdef _USE_MASS_KNOTS_
	   //           I_MK=this->cwclass.SKNOT_M_info[ig];
           I_MK= (this->cwclass.cwt_used[this->cwclass.get_Tclassification(ig)]== I_KNOT ? this->cwclass.SKNOT_M_info[ig]: 0);
#endif


           int I_CWV=0;
#ifdef _USE_CWC_V_
           I_CWV=this->cwclass.get_Vclassification(ig);
#endif

           int I_VK=0;
#ifdef _USE_VEL_KNOTS_V_
           I_VK= (this->cwclass.cwv_used[I_CWV]== I_KNOT ? this->cwclass.VDISP_KNOT_info[ig]: 0);
#endif


           int I_C1=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
           real_prec C1 = this->cwclass.Invariant_TF_II[ig];
           I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#elif defined _USE_DELTA2_
           real_prec C1 = this->cwclass.DELTA2[ig];
           I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#endif


           int I_C2=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
           real_prec C2 = this->cwclass.Invariant_TF_III[ig];
           I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, this->s_deltas.prop5,this->bin_accumulate_borders);
#elif defined _USE_DELTA3_
           real_prec C2 = this->cwclass.DELTA3[ig];
           I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#endif

           int I_C3=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_IV_
           real_prec C3 = this->cwclass.Invariant_TF_IV[ig];
           I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, this->s_deltas.prop6,this->bin_accumulate_borders);
#elif defined (USE_TIDAL_ANISOTROPY_)
           real_prec C3 = this->cwclass.Tidal_Anisotropy[ig];
           I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined _USE_S2_
           real_prec C3 = this->cwclass.S2[ig];             // s²
           I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#endif

           int I_CV1=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
           real_prec CV1 = this->cwclass.Invariant_VS_I[ig];
           I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1, this->s_deltas.prop7, this->bin_accumulate_borders);
#elif defined _USE_NABLA2DELTA_
           real_prec CV1 = this->cwclass.N2D[ig];      // Nabla² ð
           I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1, s_deltas.prop7,this->bin_accumulate_borders);
#endif

           int I_CV2=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
           real_prec CV2 = this->cwclass.Invariant_VS_II[ig];
           I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2,this->s_deltas.prop8, this->bin_accumulate_borders);
#elif defined _USE_S2DELTA_
           real_prec CV2 = this->cwclass.S2DELTA[ig];         // s²ð
           I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2, s_deltas.prop8,this->bin_accumulate_borders);
#endif

           int I_CV3=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
           real_prec CV3 = this->cwclass.Invariant_VS_III[ig];
           I_CV3= get_bin(CV3, this->s_mins.prop9, N_CV_BIN3,this->s_deltas.prop9, this->bin_accumulate_borders);
#elif defined _USE_S3_
           real_prec CV3 = this->cwclass.S3[ig];                                   // s³
           I_CV3= get_bin(CV3, this->s_mins.prop9, N_CV_BIN3, s_deltas.prop9,this->bin_accumulate_borders);
#endif



#ifndef _BIN_ACCUMULATE_
           if(halo_prop <=this->s_maxs.prop0 && halo_prop >=this->s_mins.prop0)
#endif


#ifndef _BIN_ACCUMULATE_
#ifdef _USE_DM_IN_BAM_
             if(xdm <=this->s_maxs.prop1 && xdm >=this->s_mins.prop1)
#endif
#endif

#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
               if(C1<= this->s_maxs.prop4 && C1>= this->s_mins.prop4)
#endif
#endif

#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
                 if(C2<= this->s_maxs.prop5 && C2>= this->s_mins.prop5)
#endif
#endif

#ifndef _BIN_ACCUMULATE_
#if defined (_USE_TIDAL_ANISOTROPY_)  || defined (_USE_INVARIANT_TIDAL_FIELD_III_)
                   if(C3<= this->s_maxs.prop6 && C3>= this->s_mins.prop6)
#endif
#endif

#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
                     if(CV1<= this->s_maxs.prop7 && CV1>= this->s_mins.prop7)
#endif
#endif

#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
                       if(CV2<= this->s_maxs.prop8 && CV2>= this->s_mins.prop8)
#endif
#endif

#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
                         if(CV3<= this->s_maxs.prop9 && CV3>= this->s_mins.prop9)
#endif
#endif

                           {
                             ULONG index_dm=index_11d(I_X, I_CWT, I_MK, I_CWV,I_VK, I_C1, I_C2, I_C3, I_CV1, I_CV2, I_CV3, this->n_cwt, this->n_sknot_massbin, this->n_cwv, this->n_vknot_massbin,N_C_BIN1,N_C_BIN2,N_C_BIN3, N_CV_BIN1, N_CV_BIN2,N_CV_BIN3);
                             ULONG Index=index_2d(I_Y, index_dm, LENGHT_BIAS_NCOUNT_aux);

#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                             this->BIAS_NCOUNTS[Index]++; // Add the number of cells that fall in this bin of Theta Index

			                     }
         }
       this->So.DONE();


       //* Final check done only in case we use counts-in-cells
       // This is only done when we gemnerate mocks with the same volume of the reference
#ifdef _FULL_VERBOSE_
       So.message_screen("Checking the number of grid-cells accounted for in BIAS_NCOUNTS:");
#endif

#ifndef _EXTRAPOLATE_VOLUME_
       ULONG ncells_used = static_cast<ULONG>(get_nobjects(this->BIAS_NCOUNTS));
       if(ncells_used<this->NGRID){
	 So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"The number of cells counted in BIAS(Y,X) is smaller than NGRID. Perhaps the range for some of the DMF properties *must* be increased to include all cells. CosmiCatlas stops here");
       }
       else
	 this->So.DONE();
#endif


#ifdef _FULL_VERBOSE_
       So.message_screen("Normalizing using",NTHREADS,"threads");
#endif

       //     int nbins_y_temp = property==_COUNTS_ ?  this->new_nbins_y: this->new_nbins_y_MW;
       int nbins_y_temp = this->new_nbins_y;


#if !defined (_ASSIGN_TO_REFERENCE_) || defined (_DO_BAM_CALIBRATION_)
       this->BIAS_NCOUNTS_normalized.shrink_to_fit();
       this->BIAS_NCOUNTS_normalized.clear();
       this->BIAS_NCOUNTS_normalized.resize(LENGHT_BIAS_NCOUNT, 0);
#endif



long aux_a_full=-10000;
long aux_b_full;
ULONG aux_b_empty=0;

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
 for(ULONG i=0;i < LENGHT_BIAS_NCOUNT_aux;++i)// loop sobre bines de dark matter
    {
       long aux_a=-10000;
       long aux_b;
       for(int tr_j = 0; tr_j < nbins_y_temp; ++tr_j)
         {
           ULONG inde=index_2d(tr_j, i, LENGHT_BIAS_NCOUNT_aux);
           long AUX=this->BIAS_NCOUNTS[inde];
           if(AUX==0)
            aux_b_empty++;
           aux_b=max(aux_a, AUX);
           aux_a=aux_b;
           aux_b_full=max(aux_a_full, AUX);
           aux_a_full=aux_b_full;
          }
       for(int tr_j = 0; tr_j < nbins_y_temp; ++tr_j)
         {
           ULONG inde=index_2d(tr_j, i, LENGHT_BIAS_NCOUNT_aux);
           long AUX=this->BIAS_NCOUNTS[inde];
           this->BIAS_NCOUNTS_normalized[inde] = (aux_a<=0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
         }
	  }

    double effecive_nbins=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:effecive_nbins)
#endif
   for(ULONG i=0;i < this->BIAS_NCOUNTS.size();++i)// loop sobre bines de dark matter
    effecive_nbins+=this->BIAS_NCOUNTS[i];
  effecive_nbins/=static_cast<double>(aux_b_full);
#ifdef _FULL_VERBOSE_
   So.message_screen("Effective number of bins =",effecive_nbins);
   So.message_screen("Number of empty bins =",aux_b_empty);
   So.message_screen("Number of used bins =",this->BIAS_NCOUNTS.size()-aux_b_empty);
#endif

   this->So.DONE();

   real_prec aux_h=static_cast<real_prec>(get_nobjects(this->BIAS_NCOUNTS_normalized));
   if(aux_h<=0)
     {
       So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Normalized Joint probability for number counts is ill-defined. CosmicAtlas stops here");
       exit(0);
     }


   aux_h=0;


#ifdef _UNDER_BIASED_
  if(0==this->step)
  {

      this->File.write_array(this->Output_directory+"Bam_Raw_Bias", this->BIAS_NCOUNTS);
      this->File.write_array(this->Output_directory+"Bam_Raw_Bias_Normalized", this->BIAS_NCOUNTS_normalized);
  }
#endif


#ifndef _TEST_THRESHOLDS_RESIDUALS_
   auto out_it = std::find(std::begin(this->output_at_iteration), std::end(this->output_at_iteration), this->step);
   if( (this->step==this->N_iterations_Kernel) || (out_it != std::end(this->output_at_iteration)))
     {
       this->File.write_array(this->Output_directory+"Bam_Bias", this->BIAS_NCOUNTS);
       this->File.write_array(this->Output_directory+"Bam_Bias_Normalized", this->BIAS_NCOUNTS_normalized);

       vector<real_prec>BIAS_AUX(LENGHT_BIAS_NCOUNT_aux,0);


#pragma omp parallel for collapse(2)
   for(ULONG i=0;i<LENGHT_BIAS_NCOUNT_aux;++i)
     for(int tr_j = 0; tr_j < nbins_y_temp; ++tr_j)
      {
       ULONG inde_l=index_2d(tr_j,i,LENGHT_BIAS_NCOUNT_aux);
#pragma omp atomic update
       BIAS_AUX[i]+=this->BIAS_NCOUNTS[inde_l];
     }



     this->File.write_array(this->Output_directory+"Bam_Bias_AUX", BIAS_AUX);
     BIAS_AUX.clear(); BIAS_AUX.shrink_to_fit();
    }  // closes  if( (this->step==this->N_iterations_Kernel) || (out_it != std::end(this->output_at_iteration)))
#endif

     }  // closes if(this->iteration_ini>0 && true==this->use_iteration_ini) //NOT WORKING YET, mut be deprecated
#endif // endif for _GET_BAM_REALIZARTIONS_

 }
//end class member function
#endif


//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################

 // This function is used in mode _GET_BAM_REALIZATIONS_ and gets the approximated DMF
 // Action: get new alpt dm filed. Get bam kernel and convolve
 // This fuctionn doies not need to know if we use cwc insde loops, for those loope are meant for the calibration procedure.
 void Bam::get_new_DM_field()
 {
    So.enter(__PRETTY_FUNCTION__);

#ifdef _USE_OMP_TEST_
   int NTHREADS=_NTHREADS_;   // omp_get_max_threads();
   omp_set_num_threads(NTHREADS);
#else
   int NTHREADS=1;   // omp_get_max_threads();
#endif   


#ifndef _ONLY_POST_PROC_
#ifdef _FULL_VERBOSE_
   cout<<endl;
   this->So.message_screen("*************************************************************************");
   this->So.message_screen("****Getting new DM density field in line", __LINE__, "of file",__FILE__);
   this->So.message_screen("*************************************************************************");
   cout<<endl;
#endif
#endif

   real_prec num_in_log_x = (true==this->Convert_Density_to_Delta_X ? NUM_IN_LOG: 0.);

   string file_X;
   string file_Y;
   string file_Vx;
   string file_Vy;
   string file_Vz;


   this->delta_X.clear();
   this->delta_X.shrink_to_fit();
   this->delta_X.resize(this->NGRID,0);
   this->delta_X_ini.clear();
   this->delta_X_ini.shrink_to_fit();
   this->delta_X_ini.resize(this->NGRID,0);

   int n_realization=this->step - this->N_iterations_Kernel + this->N_dm_initial-1;
   real_prec nmean=0;

   vector<real_prec>xbaux(this->NX, 0);
   vector<real_prec>pdf_in(this->NX, 0);

#ifndef _USE_PATCHY_  //If we cannot run patchy, we need to read DM files from somewhere, so we read them from paths hard coded here
   this->File.read_array(this->Input_Directory_X_NEW+this->Name_Catalog_X_NEW, this->delta_X_ini);



#ifdef _USE_VELOCITIES_
   file_Vx="../";
   file_Vy="../";
   file_Vz="../";
   this->Velx_X.resize(this->NGRID,0);
   this->Vely_X.resize(this->NGRID,0);
   this->Velz_X.resize(this->NGRID,0);
   this->File.read_array(file_Vx, this->Velx_X);
   this->File.read_array(file_Vy, this->Vely_X);
   this->File.read_array(file_Vz, this->Velz_X);
#endif

   nmean=get_nobjects(this->delta_X_ini);
#ifdef _FULL_VERBOSE_
   So.message_screen("Total number of X objects (new) =", nmean);
#endif

   nmean/=static_cast<real_prec>(this->NGRID);

#ifdef _FULL_VERBOSE_
   So.message_screen("Average number of X objects (new) =", nmean);
   cout<<endl;
#endif

#ifndef _ASSIGN_TO_REFERENCE_  // just to save tmie while doing the test of mass assingment to tracers
#ifdef _GET_POWER_REFS_
   this->get_power_spectrum("DM_REF_NEW");
#endif
#endif

    //andres
  cout<<RED<<this->delta_X_ini[1255]<<__LINE__<<RESET<<endl;

   get_overdens(this->delta_X_ini, nmean, this->delta_X_ini);
    //andres
  cout<<RED<<this->delta_X_ini[1255]<<__LINE__<<RESET<<endl;


#ifdef _RANK_ORDERING_AB_INITIO_

#ifdef _FULL_VERBOSE_
       So.message_screen("Rank ordering to target field");
#endif

#ifdef _FULL_VERBOSE_
       So.message_screen("Measuring pdf of log(1+delta) of input DM field");
#endif
       calc_pdf("log", this->NGRID,this->NX, this->Xmax, this->Xmin, this->delta_X_ini, pdf_in);
#ifdef _FULL_VERBOSE_
       So.DONE();
#endif
       this->delta_X_REF_PDF.resize(this->NGRID,0);
       string file_X_ref_pdf=this->Input_Directory_X_REF+this->Name_Catalog_X_REF_PDF;
       this->File.read_array_t<PrecType_X>(file_X_ref_pdf, this->delta_X_REF_PDF); // Read target DM field for RO
#ifdef _RO_WITH_DELTA_
       get_overdens(this->delta_X_REF_PDF, this->delta_X_REF_PDF); //a better bispect is found if the rank ordering is done to the CIC of number counts, not the delta thereof
#endif

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(int i=0;i<this->NX; ++i)
         xbaux[i]=this->Xmin+static_cast<real_prec>(i+0.5)*(this->Xmax-this->Xmin)/(static_cast<real_prec>(this->NX));

#ifdef _FULL_VERBOSE_
       So.message_screen("Measuring pdf of log(1+delta) of target reference");
#endif
       this->pdf_ref.resize(this->NX, 0);
       calc_pdf("log", this->NGRID,this->NX, this->Xmax, this->Xmin, this->delta_X_REF_PDF, this->pdf_ref);
#ifdef _FULL_VERBOSE_
       So.DONE();
#endif
       this->delta_X_REF_PDF.clear();
       this->delta_X_REF_PDF.shrink_to_fit();
#ifdef _FULL_VERBOSE_
       So.message_screen("Executing rank ordering from DM to DM-target");
#endif
       rankorder(ZERO, xbaux, this->NGRID, this->NX,  this->Xmax, this->Xmin, this->delta_X_ini, pdf_in, this->pdf_ref);
       So.DONE();

#endif //end if apply rank ordering ab initio

       // we do not need else for the output of rank ordering is again delta_X_ini
#else   //if we have to use patchy
   {

     const gsl_rng_type *rng_t;
#ifdef _USE_OMP_
     gsl_rng **gBaseRand;
#else
     gsl_rng *gBaseRand;
#endif
     gsl_rng_env_setup();
     rng_t = gsl_rng_mt19937;// gsl_rng_default;
#ifdef _USE_OMP_
     int nt=omp_get_max_threads();
     gBaseRand = (gsl_rng **) malloc(nt * sizeof(gsl_rng *));
#endif

#ifdef _USE_OMP_
#pragma omp parallel for num_threads(nt)
#endif
     for(int i=0;i<nt;i++)
       {
         gBaseRand[i] = gsl_rng_alloc(rng_t);
         gsl_rng_set(gBaseRand[i],this->seed);
       }

     time_t start_patchy;
     time(&start_patchy);

     if(false==dm_already_done) //COmpute DM if it has not been done before through a call of this function
       {
#ifdef _USE_OMP_
         this->patchy.get_dm_field(gBaseRand);
#else
         this->patchy.get_dm_field(gBaseRand);
#endif
         this->dm_already_done=true;
         this->So.message_screen("Patchy has created DMDF using ALPT");
         cout<<endl;

       }
     // Patchy has written in files, so now we read them
     file_X=this->patchy.fnameDM+".dat";

#ifdef _USE_TRACER_HR_
     file_Y=this->Input_Directory_Y+this->Name_Catalog_Y_HR;
#else
     file_Y=this->Input_Directory_Y+this->Name_Catalog_Y;
#endif

     this->delta_X_ini.resize(this->NGRID,0);

#ifdef _USE_VELOCITIES_
     //here we have to read the vels of the newly created density field
     file_Vx=this->patchy.fnameVX+".dat";
     file_Vy=this->patchy.fnameVY+".dat";
     file_Vz=this->patchy.fnameVZ+".dat";
     this->Velx_X.resize(this->NGRID,0);
     this->Vely_X.resize(this->NGRID,0);
     this->Velz_X.resize(this->NGRID,0);
     this->File.read_array(file_Vx, this->Velx_X);
     this->File.read_array(file_Vy, this->Vely_X);
     this->File.read_array(file_Vz, this->Velz_X);
#endif // end ifdef _USE_VELOCITIES_

     // read the just created density field
     this->File.read_array(file_X, this->delta_X_ini);


#ifdef _GET_POWER_REFS_
     this->get_power_spectrum("DM_REF_NEW");
#endif
     nmean=get_nobjects(this->delta_X_ini);
     So.message_screen("Total number of X objects (new) =", nmean);
     nmean/=static_cast<real_prec>(this->NGRID);
     So.message_screen("Average number of X objects (new) =", nmean);

     // Tansform input DF to Overdensity field
     get_overdens(this->delta_X_ini, this->delta_X_ini);

   } // end else uf use patchy
#endif


   this->pdf_ini.resize(this->NX, 0);
   string rest_file="_Nft"+to_string(this->Nft)+"_z"+to_string(this->redshift)+".txt";
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int i=0;i<this->NX; ++i)
     xbaux[i]=this->Xmin+static_cast<real_prec>(i+0.5)*(this->Xmax-this->Xmin)/(static_cast<real_prec>(this->NX));

#ifdef _FULL_VERBOSE_
   So.message_screen("Measuring pdf of log(1+delta) DM: ");
#endif
   calc_pdf("log",this->NGRID,this->NX, this->Xmax, this->Xmin, this->delta_X_ini, pdf_in);
   this->pdf_ini=pdf_in;
#ifdef _WRITE_PDF_
   string filex=this->Output_directory+"pdf_X_NEW_"+this->XNAME+"_MASX"+to_string(this->iMAS_X)+rest_file;
   this->File.write_to_file(filex, xbaux,pdf_in);
#endif

   // Now that the new field is loadad, convolve with kernel and then get the T or I WEB

#ifdef _USE_MASS_KNOTS_
   this->cwclass.SKNOT_M_info.resize(this->NGRID, 0);
#endif

   // The order is important: if we do CWC inside each loop,
   // then we convolve the target DF with the updated kernel and then do the CWC
   // Otherwise one does first the CWC once and then convolves in each iteration


   this->delta_X.resize(this->delta_X_ini.size(),0);


#ifndef _DO_NOT_CONVOLVE_
#ifdef _FULL_VERBOSE_
    So.message_screen("Reading Kernel ");
#endif
    this->Kernel.clear();
    this->Kernel.shrink_to_fit();
    this->Kernel.resize(this->NTT, 0.0);
    this->File.read_array(this->params._Input_Directory_BIAS_KERNEL()+"Bam_Kernel.dat", this->Kernel);
    //andres
  cout<<RED<<this->Kernel[1255]<<"  "<<this->delta_X_ini[1255]<<__LINE__<<RESET<<endl;

    this->Konvolve(this->delta_X_ini, this->delta_X);
    //andres
  cout<<RED<<this->Kernel[1255]<<"  "<<this->delta_X[1255]<<__LINE__<<RESET<<endl;



#ifdef _WRITE_DM_DENSITY_FIELD_
        this->File.write_array(this->Output_directory+"DM_DELTA_convolved_realization", this->delta_X);
#endif
#else
        delta_X=delta_X_ini;
#endif



#ifdef _USE_CWC_
   if(this->n_cwt > 1)
     { 
#endif
// Open the kernel


#ifdef _USE_CWC_
       this->cwclass.do_CWC(this->delta_X);   //get the CW info
#ifdef _USE_MASS_KNOTS_
       this->cwclass.get_Mk_collapsing_regions(this->delta_X,nmean);  //get the MK info
#endif //use_mass_knots


#elif !defined _USE_CWC_
#if defined (_USE_MASS_KNOTS_)  
       this->cwclass.do_CWC(this->delta_X);   //get the CW info
       this->cwclass.get_Mk_collapsing_regions(this->delta_X,nmean);  //get the MK info
#endif // use_mass_knots || other stff

#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_INVARIANT_TIDAL_FIELD_IV_)|| defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_ELLIPTICITY_) || defined (_USE_PROLATNESS_) || defined (_USE_PWEB_)
       this->cwclass.do_CWC(this->delta_X);   //get the CW info
#endif // use_mass_knots || other stff
#endif    // !use_cwc




          // If the terms in the bias expansion are to be used, then :
#if (defined (_USE_NABLA2DELTA_) || defined (_USE_S2_) || defined (_USE_S3_) || defined (_USE_S2DELTA_) ) && (!defined (_USE_CWC_))
       this->cwclass.get_bias_terms(this->delta_X);
#endif

#ifdef _USE_CWC_
   }
#endif // end ifdef use_cwc






#ifdef _GET_POWER_REFS_
     this->get_power_spectrum("DM_KONV"); // This has to be done *BEFORE* transfoming to LOg
#endif




#ifdef _FULL_VERBOSE_
          So.message_screen("Measuring pdf of log(1+delta) DM target reference");
#endif


     calc_pdf("log", this->NGRID,this->NX, this->Xmax, this->Xmin, this->delta_X, pdf_in);
     So.DONE();
#ifdef _WRITE_PDF_
     rest_file="_Nft"+to_string(this->Nft)+"_z"+to_string(this->redshift)+".txt";
     filex=this->Output_directory+"pdf_X_realization_"+this->XNAME+"_MASX"+to_string(this->iMAS_X_REF_PDF)+rest_file;
     this->File.write_to_file(filex, xbaux,pdf_in);
#endif
      pdf_in.clear(); pdf_in.shrink_to_fit();


   if(this->Scale_X=="log")
     {
#ifdef _FULL_VERBOSE_
       cout<<endl;
       // If the mass_iterative_nes is nit to be used, here we directly convert to log(numinlog+delta)
       So.message_screen("Transforming new delta->log10(2+delta)");
#endif


//andres
//  cout<<RED<<this->delta_X[1255]<<"  "<<this->cwclass.Invariant_TF_II[1255]<<"  "<<this->cwclass.Invariant_TF_III[1255]<<"  "<<__LINE__<<RESET<<endl;


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(ULONG i = 0;i < this->NGRID ;++i)  //TRANSFORM DELTA TO LOG10(NUM_IN_LOG + DELTA)
         this->delta_X[i] = this->delta_X[i]<-1 ?  0 : log10(num_in_log_x + static_cast<real_prec>(this->delta_X[i]));
       So.DONE();
     }
   // Once the delta field is obtained, get the limits. Note that these limits mut be those of the last iteration
   // of the iterative procedure,
   // but it is not guaranteed that this is the case. So, if in _GET_BAM_REALIZATIONS_mde, Leave the limits fixed for the time being.



   this->get_new_min_max_properties();

 }


 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
#ifdef MOCK_MODE


#ifdef _GET_BAM_REALIZATIONS_
//#define _new_app_
#define _achtung_inside_  // This check inside the mock loop whether the cond probability is >0, such that we do not fall in infinite while loops
#endif

 //#define _use_filled_cells_  Use this when cañlibrating from something else such as sat fraction or mass

 void Bam::get_mock_grid(string property)
 {
   So.enter(__PRETTY_FUNCTION__);

   bool silent=true;
#ifdef _USE_OMP_
   int NTHREADS = _NTHREADS_;  // omp_get_max_threads();
   omp_set_num_threads(NTHREADS);
#endif

#ifdef _FULL_VERBOSE_
#ifdef _USE_OMP_
   So.message_screen("Function get_mock_grid is using",omp_get_max_threads(),"threads");
#endif
   cout<<endl;
#endif
   // ********************************************************************************************************************

    int nx=1;
#ifdef _USE_DM_IN_BAM_
   // Define seed vector for each thread in the paralellized run
   nx=this->NX;  // Number of bins id DM
#endif



#ifdef _use_filled_cells_
   vector<bool>filled_cells(this->NGRID, true);
   vector<real_prec> ncounts(this->NGRID,0); // this is only used in case:_prop=2, whcih was halo mass in cells, deprecated
#endif
   string pname,fname;
   int case_prop=0; int ny=0;
   ULONG size_bias_array=0;

   if(_COUNTS_==property)
     {
       case_prop=1;
       size_bias_array=this->BIAS_NCOUNTS.size();
       ny = this->new_nbins_y;

       if(this->step <=this->N_iterations_Kernel)
         pname ="_iteration"+to_string(this->step);
       else
         pname= "_realization"+to_string(this->params._realization());
       fname=this->Output_directory+"MOCK_TR"+pname+"_"+"MASY"+to_string(this->iMAS_Y)+"_Nft"+to_string(this->Nft)+"_z"+to_string(this->redshift);
     }

   if(_DENSITY_==property)
     {
       case_prop=2;
       size_bias_array=this->BIAS_NCOUNTS.size();
       ny = this->new_nbins_y;
       fname=this->Output_directory+"MOCK_TR"+pname+"_"+"MASY"+to_string(this->iMAS_Y)+"_Nft"+to_string(this->Nft)+"_z"+to_string(this->redshift);
     }
   this->fnameMOCK=fname;


   // Initialize counter for those cells for which no available positions were found
   ULONG counter_orphan=0;
   ULONG counter_out_of_theta=0;
   ULONG counter_orphan_b=0;

   time_t start_mock;
   time(&start_mock);

   // Allocate memmory for the mock density field
   this->delta_Y_new.resize(this->NGRID,0);

#ifdef _DYNAMICAL_SAMPLING_
   // Vector to allocate the Joint distribution updated after assigning Nhalos to a cell.
   vector<ULONG> X_Y_hist_dynamical(size_bias_array);
   vector<real_prec>  X_Y_hist_dynamical_normalized(size_bias_array);
#endif


   // Initialize these vectors with the original distribution in each density bin.

#if defined _USE_MASS_FIELD_
   switch(case_prop)
     {
     case(1):
#endif


#ifdef _DYNAMICAL_SAMPLING_
     X_Y_hist_dynamical=this->BIAS_NCOUNTS;
     X_Y_hist_dynamical_normalized=this->BIAS_NCOUNTS_normalized;


#ifdef _new_counter_
   vector<ULONG> new_X_Y_hist_dynamical(size_bias_array);
   vector<real_prec>  new_X_Y_hist_dynamical_normalized(size_bias_array);
   new_X_Y_hist_dynamical=this->BIAS_NCOUNTS;
   new_X_Y_hist_dynamical_normalized=this->BIAS_NCOUNTS_normalized;
#endif


#if defined _USE_SAT_FRAC_ || defined _USE_MASS_FIELD_
       break;
#endif

#ifdef _USE_MASS_FIELD_
     case(3):
       X_Y_hist_dynamical=this->BIAS_SAT_FRACTION;
       X_Y_hist_dynamical_normalized=this->BIAS_SAT_FRACTION_normalized;
       break;
#endif
#if defined _USE_MASS_FIELD_
   }
#endif

   // Vector allocating the Joint distribution normalized within each Den bin, after having assigned a value of Nhalos to a cell.
   // Initialize these vectors with the original Joint and the normalized distribution in each density bin.

#endif  // enf if def _DYNAMICAL_SAMPLING_

   ULONG lenght_bias=N_CV_BIN1 * N_CV_BIN2* N_CV_BIN3 * N_C_BIN1 * N_C_BIN2* N_C_BIN3 * this->n_sknot_massbin * this->n_cwv * this->n_vknot_massbin * this->n_cwt * nx;


   // Vector containing the total number if cells in a given density bin and CWT and KNOT mass
   this->NCELLSperDMBIN.clear();
   this->NCELLSperDMBIN.shrink_to_fit();
   this->NCELLSperDMBIN.resize(lenght_bias, 0);

   // These two loops replace the 11 loops below. tHIS MODIFICATIONB MUST BE DONE CCONSISTENTLYT EVERYWHERE WHERE BIAS IS LOADED


#ifdef _GET_BAM_REALIZATIONS_

   this->File.read_array(this->params._Input_Directory_BIAS_KERNEL()+"Bam_Bias_AUX.dat", this->NCELLSperDMBIN);

#else

#ifdef _USE_OMP_
#pragma omp parallel for collapse(2)
#endif
   for(ULONG i=0; i< lenght_bias ;++i)
       for(ULONG j=0; j< ny ; ++j)
         {
           ULONG index_h=index_2d(j,i,lenght_bias);
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
           this->NCELLSperDMBIN[i]+= this->BIAS_NCOUNTS[index_h];
          }

#endif // enf ifdef get bam realizations
   {
#ifndef _EXTRAPOLATE_VOLUME_
     if(property==_COUNTS_) // We have only allowed the number counts bias to have all cells.
       {
#ifdef _FULL_VERBOSE_
         So.message_screen("Double-checking the number of grid-cells accounted for in BIAS_NCOUNTS:");
#endif
         ULONG KK = static_cast<ULONG>(get_nobjects(this->BIAS_NCOUNTS));
         if(KK!=this->NGRID)
           {
             So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Missing cells in get_mock_grid(). Perhaps density range is not wide enough to contain all cells", KK);
             exit(0);
            }
         else
             So.DONE();
#endif
         real_prec Kd = static_cast<real_prec>(get_nobjects(this->BIAS_NCOUNTS_normalized));
         if(Kd<=0){
           So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Conditional probablity ill-defined. BAM stops here.", Kd);
           exit(0);
           }
         }
   }

   // Vector to allocate the number of cells in a given density bin during the mapping */
   vector<ULONG>Ncells_density_bin_new(size_bias_array , 0);

   // Vector containing the number of cells in a given density bin, updated everytime a cell has been assigend a value of Nhalos
   vector<ULONG> NCELLSperDMBIN_now(lenght_bias, 0);




#ifdef _GET_BAM_REALIZATIONS_
#ifdef _FULL_VERBOSE_
   So.message_screen("**Generating new Mock (tracer) number counts.");
   cout<<endl;
#endif
#else
#ifdef _FULL_VERBOSE_
   So.message_screen("**Generating Mock (tracer) number counts @ iteration ", this->step);
   cout<<endl;
#endif
#endif

#ifdef _new_app_
    string file_Pdf=this->Output_directory+"PDF_NC_Y_REF_COUNTS_MASY"+to_string(this->iMAS_Y)+"_Nft"+to_string(this->Nft)+"_z"+to_string(this->redshift)+".txt";
    So.message_screen("Reading PDF from reference halo catalog in file", file_Pdf);
    ifstream pdf_file; pdf_file.open(file_Pdf);
    vector<ULONG> pdf_ref(ny+1,0);
    vector<ULONG> pdf_temp(ny+1,0);
    int iaux;
    for(int i =0; i < ny +1;++i)
      pdf_file>>iaux>>pdf_ref[i];
    pdf_file.close();
#ifdef _FULL_VERBOSE_
    So.DONE();
#endif
#endif



#ifdef _GET_BAM_REALIZATIONS_
#ifdef _USE_OMP_
    NTHREADS=1;
    omp_set_num_threads(NTHREADS);
#endif
#endif

#ifdef _FULL_VERBOSE_
#ifdef _USE_OMP_
        So.message_screen("Using",omp_get_max_threads(),"threads");
#endif
#endif
    const gsl_rng_type * rng_t;
    gsl_rng * gBaseRand;
    int jthread=0;
    gsl_rng_env_setup();

#ifdef _USE_OMP_
/*
        //initialize gsl random number generator
  const gsl_rng_type *rng_t;
  gsl_rng **gBaseRand;
  gsl_rng_env_setup();
  rng_t = gsl_rng_mt19937;// gsl_rng_default;

  int nt=NTHREADS;
  int jthread;
  gBaseRand = (gsl_rng **) malloc(nt * sizeof(gsl_rng *));

#pragma omp parallel for num_threads(nt)
  for(int i=0;i<NTHREADS;i++)
    {
      gBaseRand[i] = gsl_rng_alloc(rng_t);
      gsl_rng_set(gBaseRand[i],seed);
    }
*/


   vector<ULONG>vseeds(NTHREADS,0);
   for(int i=0;i<vseeds.size();++i)
     vseeds[i]=35+static_cast<ULONG>(i)*565;
#endif

#ifdef _new_app_

   struct s_scalar_cell_info{
     ULONG theta_index;
     bool assigned;
   };

    vector<s_scalar_cell_info> cell_info(this->NGRID);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->NGRID;++i)
    {
      cell_info[i].assigned=false;
      cell_info[i].theta_index=0;
   }
#endif



#ifdef _USE_OMP_
#pragma omp parallel private (jthread, gBaseRand, rng_t)
   {
#endif

#ifdef _USE_OMP_
     jthread=omp_get_thread_num();
     gsl_rng_default_seed=vseeds[jthread];
#else
     gsl_rng_default_seed=35;
#endif

     rng_t = gsl_rng_mt19937;//_default;
     gBaseRand = gsl_rng_alloc (rng_t);

#ifdef _USE_OMP_
#pragma omp for reduction(+:counter_orphan, counter_out_of_theta)
#endif
     //Start loop over the cells
     for(ULONG i=0;i<this->NGRID;++i)        // First block, meant to construct halo number counts
       {

         int I_X=0;
#ifdef _USE_DM_IN_BAM_
         real_prec dm = this->delta_X[i];
         // Get the bin in the delta dark matter (or log 1+delta) in each cell
         I_X  = ((this->iMAS_X == 0  && false==this->Convert_Density_to_Delta_X) ? static_cast<int>(this->delta_X[i]) : get_bin(dm,this->s_mins.prop1,this->NX,this->s_deltas.prop1,this->bin_accumulate_borders));
#endif

         // **********CWT
         int I_CWT=0;
#ifdef _USE_CWC_
         I_CWT=this->cwclass.get_Tclassification(i);
#endif

         // **********MK
         // Get the corresponding bin on Knot-mass. The current cell is in a cluster with a given mass falling in that particular bin
         int I_MK=0;
#ifdef _USE_MASS_KNOTS_
         I_MK= (this->cwclass.cwt_used[this->cwclass.get_Tclassification(i)]== I_KNOT ? this->cwclass.SKNOT_M_info[i]: 0);
#endif


         // **********CW-V
         int I_CWV=0;
#ifdef _USE_CWC_V_
         I_CWV=this->cwclass.get_Vclassification(i);
#endif

         // **********MV
         int I_VK=0;
#ifdef _USE_VEL_KNOTS_V_
         I_VK= (this->cwclass.cwv_used[I_CWV]== I_KNOT ? this->cwclass.VDISP_KNOT_info[i]: 0);
#endif


         // **********C1
         // Get the corresponding bin in the two invariants of the shear of the tidal field
         int I_C1=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
         real_prec C1 = this->cwclass.Invariant_TF_II[i];
         I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#elif defined _USE_DELTA2_
         real_prec C1 = this->cwclass.DELTA2[i];
         I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#endif

         // **********C2
         int I_C2=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
         real_prec C2 = this->cwclass.Invariant_TF_III[i];
         I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#elif defined _USE_DELTA3_
         real_prec C2 = this->cwclass.DELTA3[i];
         I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#endif

       int I_C3=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_IV_
         real_prec C3 = this->cwclass.Invariant_TF_IV[i];
         I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined (_USE_TIDAL_ANISOTROPY_)
         real_prec C3 = this->cwclass.Tidal_Anisotropy[i];
         I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined _USE_S2_
         real_prec C3 = this->cwclass.S2[i];             // s²
         I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#endif

         int I_CV1=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
         real_prec CV1 = this->cwclass.Invariant_VS_I[i];
         I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1,this->s_deltas.prop7, this->bin_accumulate_borders);
#elif defined _USE_NABLA2DELTA_
         real_prec CV1 = this->cwclass.N2D[i];      // Nabla² ð
         I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1, s_deltas.prop7,this->bin_accumulate_borders);
#endif

         int I_CV2=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
         real_prec CV2 = this->cwclass.Invariant_VS_II[i];
         I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2,this->s_deltas.prop8, this->bin_accumulate_borders);
#elif defined _USE_S2DELTA_
         real_prec CV2 = this->cwclass.S2DELTA[i];         // s²ð
         I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2, s_deltas.prop8,this->bin_accumulate_borders);
#endif


         int I_CV3=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
         real_prec CV3 = this->cwclass.Invariant_VS_III[i];
         I_CV3= get_bin(CV3, this->s_mins.prop9, N_CV_BIN3,this->s_deltas.prop9, this->bin_accumulate_borders);
#elif defined _USE_S3_
         real_prec CV3 = this->cwclass.S3[i];                                   // s³
         I_CV3= get_bin(CV3, this->s_mins.prop9, N_CV_BIN3, s_deltas.prop9,this->bin_accumulate_borders);
#endif

#ifndef _BIN_ACCUMULATE_
         if(dm >=this->s_mins.prop1 && dm <=this->s_maxs.prop1)
#endif

#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_DELTA2_)
           if(C1<= this->s_maxs.prop4 && C1>= this->s_mins.prop4)
#endif
#endif


#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_DELTA3_)
             if(C2<= this->s_maxs.prop5 && C2>= this->s_mins.prop5)
#endif
#endif


#ifndef _BIN_ACCUMULATE_
#if defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_S2_) || defined (_USE_INVARIANT_TIDAL_FIELD_IV_)
               if(C3<= this->s_maxs.prop6 && C3>= this->s_mins.prop6)
#endif
#endif

#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_SHEAR_VFIELD_I_) || defined (_USE_NABLA2DELTA_)
                 if(CV1<= this->s_maxs.prop7 && CV1>= this->s_mins.prop7)
#endif
#endif


#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_SHEAR_VFIELD_II_) || defined (_USE_S2DELTA_)
                   if(CV2<= this->s_maxs.prop8 && CV2>= this->s_mins.prop8)
#endif
#endif


#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_SHEAR_VFIELD_III_) || defined (_USE_S3_)
                     if(CV3<= this->s_maxs.prop9 && CV3>= this->s_mins.prop9)
#endif
#endif

#ifdef _use_filled_cells_
                       if(true==filled_cells[i])
                         {
#endif
                           // Get the bin in DM



                           ULONG index_el=index_11d(I_X,I_CWT,I_MK,I_CWV,I_VK,I_C1,I_C2,I_C3,I_CV1,I_CV2,I_CV3,this->n_cwt,this->n_sknot_massbin,this->n_cwv, this->n_vknot_massbin,N_C_BIN1,N_C_BIN2,N_C_BIN3,N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);

#ifdef _new_app_
                            cell_info[i].theta_index=index_el;
#endif

#ifdef _achtung_inside_
                           // check if prob is zero for all possible valoes of Nh
                           // if aux_h=0, given that the quantity CWT_hist is always > 0, means that all values of prob in this bin are zero
                           // Hence no cells with this DM were found to have halos (for all values of ny).
                           // This is likely to happen when we apply the Bias and the Kernel to a different DM field in order to create a mock.
			                     // During the calibration procedure, this must not happen (by construction)
                           double aux_h=0;
                           for(int ih=0;ih<ny;++ih)
                               aux_h +=  static_cast<double>(this->BIAS_NCOUNTS[index_2d(ih,index_el, lenght_bias)]);//this applyes for cases 1 and 2

                           if(aux_h>0) //This applies in the case of generating mocks. This is likely the reason why parallelization in the mock production is not working, for it works for all cells......
                             {
#endif // end _achtung_inside_


                               // Number of cells in the density bin and CWT.
                               // This is always greater than zero, by construction. A fixed value inside the loop
                               ULONG N_available_positions_in_denbin_cwt_Mk = this->NCELLSperDMBIN[index_el];

                               // Number of used cells for the current bin of DM and CWT
                               ULONG N_used_positions_in_denbin_cwt_Mk = NCELLSperDMBIN_now[index_el];


                               // If the density bin to which this cell belongs to is already filled
                               // then go below to get a value of Nhalos selected according to the original distribution in this density bin
                               // setting a flag=false
                               bool flag = true;

                               if (N_available_positions_in_denbin_cwt_Mk == N_used_positions_in_denbin_cwt_Mk)
                                 flag=false;



                               // If the density bin still has available positions to be assigned, then proceed
                               if(true==flag)
                                 {
                                   bool cell_accepted=false;

                                   int halo_prop;

                                   // This density bin as available positions. Proceed until the cell is assigned a value of Nhalo
                                   while (false==cell_accepted)
                                     {
                                       // Throw one value of Nhalos in the range [0, nmax_y_onecell], where  nmax_y_onecell represents the
                                       // maximum number of tracers in one cell, read from the cell number counts.

                                       // Accept this number according to the normalized distribution of Nhalos in the corresponding density bin.
                                       // Everytime one cell is assigned this value of Nhalo, the normalized distribution is updated (below)
                                       // and computed extracting that position already assigned. This will ensure that more probability
                                       // is given to the remaining available positions.
                                       real_prec prob_ncounts=-10.0;
                                       real_prec ran = 10.0;
                                       ULONG index=static_cast<ULONG>(LARGE_NUMBER);

                                       // Probability for number counts
                                       while (prob_ncounts<ran)
                                         {
#if defined _USE_MASS_FIELD_ 
                                           if(1==case_prop || 2==case_prop)
#endif
                                           halo_prop= gsl_rng_uniform_int(gBaseRand,ny); //draw number counts in cells or bins in the Y distirbution (id _DENSITY_)

#ifdef _USE_MASS_FIELD_
                                           else if(4==case_prop)// this is not expected
                                             halo_prop= static_cast<int>(floor(static_cast<real_prec>(ny)*gsl_rng_uniform(r)));                         // draw bins n the mass distribution
#endif

                                           index  = index_2d(halo_prop,index_el,lenght_bias );
                                           prob_ncounts = static_cast<real_prec>(X_Y_hist_dynamical_normalized[index]);
                                           ran   = gsl_rng_uniform(gBaseRand);

                                       }// closes while(prob_ncounts<ran)

                                       ULONG N_available_positions=this->BIAS_NCOUNTS[index];
                                       ULONG N_used_positions = Ncells_density_bin_new[index];
                                       // Now proceed to assign the value of Nhalos to the current cell i, only in case
                                       // that we have available positions in the current Den-Nh bin.
                                       // If there are no available positions, then cell_accepted will be false and we cannot
                                       // get out pf the while loop. A new Nhalo will be proposed.

                                       //For the other properties such as Mass in a cell or the fraction of satellits, we
                                       // we let the available number of cells be the driver criteria

                                       if(N_available_positions > N_used_positions)
                                         {
                                           this->delta_Y_new[i]=static_cast<real_prec>(halo_prop);

                                           // Claim the current cell as already assigned a value of Nhalos. This breaks the while and continues to the next cell.
                                           cell_accepted=true;
#ifdef _new_app_
                                           cell_info[i].assigned=true;
#endif


#ifdef _new_app_
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                                           pdf_temp[halo_prop]++;
#endif


                                           // Add one (i.e, the accepted) to count the number of positions used in the current den-N bin
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                                           Ncells_density_bin_new[index]++;
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                                           NCELLSperDMBIN_now[index_el]++;

                                           // -----------------------------------------------------------------
                                           // In the current Den-N bin, subtract one (i.e, the accepted) in order
                                           // to upgrade the distribution and give more weight to the remaining available positions
                                           // Attention, do not ask whether this is <0, for it is an unsigned long
                                           if(X_Y_hist_dynamical[index] >=1)
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                                             X_Y_hist_dynamical[index]--;

                                           else
                                             X_Y_hist_dynamical[index]=0;

                                           // Normalize the dynamical histogram,  in each Den bin, to the maximum.
                                           // The resulting distribution, in this den bin, will be used for the next cell
                                           // in order to accept or reject the new value of Nhalos,
                                           // such that more probability is given to the remaining available positions.

                                           long aux_a=-10000;
                                           long lkk;
                                           for(int j=0; j< ny ;++j)// loop over the number of tracers in a mesh cell
                                           {
                                             long AUX=X_Y_hist_dynamical[index_2d(j,index_el, lenght_bias)];
                                             lkk=max(aux_a, AUX);
                                             aux_a=lkk;
                                           }


                                           for(int j=0;j< ny ;++j)
                                           {
                                               ULONG iin=index_2d(j,index_el, lenght_bias);
                                               ULONG AUX=X_Y_hist_dynamical[iin];
                                               X_Y_hist_dynamical_normalized[iin]= (lkk == 0  ? 0.0 :  static_cast<double>(AUX)/static_cast<double>(lkk));
                                            }
                                         }   // end   if(N_available_positions > N_used_positions) .
                                     }// end  while (false==cell_accepted)
                                 } // end of if(true==flag)
#ifndef _new_app_
                               else
                                 // If the density bin has been already filled, we call the full joint distribution in that bin to assign randomly the value Nhalo according to the other properties
                                 {

                                    real_prec prob=-1.0;
                                    real_prec ran=0.0;
                                    int Nhalos_orfan=0;

                                    while(prob<ran)
                                       {
                                          Nhalos_orfan= gsl_rng_uniform_int(gBaseRand,ny); //draw number counts in cells or bins of a density-like quantity
                                          prob = this->BIAS_NCOUNTS_normalized[index_2d(Nhalos_orfan,index_el, lenght_bias)];
                                          ran = gsl_rng_uniform(gBaseRand);
                                       }

                                   this->delta_Y_new[i]= static_cast<real_prec>(Nhalos_orfan);
                                   counter_orphan++;

                                }// end else
#endif
#ifdef _achtung_inside_
                             }// end of if(aux_h>0), available only under  _achtung_inside_

                           else // if in the theta_bin the prob is **ALWAYS zero**, assign Poisson distributed number (or zero) with mean that of the reference;
                             {
                              // if(this->step==this->N_iterations_Kernel)
                                 this->delta_Y_new[i]=0;  // This is radical. Might not allow for cosmic variance
                            //    this->delta_Y_new[i]=gsl_ran_poisson(gBaseRand,this->mean_number_y_onecell); //This gives a Poisson realization with mean=mean_part_in_cell. More realistic perhaps
                                counter_out_of_theta++;

#ifdef _new_app_
                                cell_info[i].assigned=true;
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                                pdf_temp[0]++;
#endif //ende new_app
                               }

#endif // end achtung inside

#ifdef _use_filled_cells_
                   }// end if(dm is in the defined range)
#endif

 }// end loop over cells


  this->Nobjects=static_cast<int>(get_nobjects(this->delta_Y_new));


#ifdef _new_app_

 // Given that the biased mocks are assigned more objects (power is below) we can avoid that by reducing the nmax in each cell
 // usefd
 // get the number of tracers so far assigned
// get the ratio btween the number of objects assigned and the total number of tracers form the reference
  this->delta_Y.resize(this->NGRID,0);
  this->File.read_array_t<PrecType_Y>(this->Input_Directory_Y+this->Name_Catalog_Y, this->delta_Y);
  ULONG Nobjects_ref=static_cast<int>(get_nobjects(this->delta_Y));
/*
  ny= static_cast<int>(floor(ny*static_cast<double>(this->Nobjects)/static_cast<double>(Nobjects_ref)));
  this->delta_Y.clear(); this->delta_Y.shrink_to_fit();
  So.message_screen("Reducing maximum occuipation number to", ny);
  So.message_screen("First loop done");
*/

/*
  So.message_screen("Start second loop");
#ifdef _USE_OMP_
#pragma omp for reduction(+:counter_orphan_b)
#endif
     //Start loop over the cells
     for(ULONG i=0;i<this->NGRID;++i)
         {
          if(false==cell_info[i].assigned)
            {
               ULONG index=cell_info[i].theta_index;
                int Nhalos_orfan;
                real_prec prob=-1.0;
                real_prec ran=0.0;
                while(prob<ran) // This while aims at chosing the right Nh according to the full propability distribution
                 {
                    Nhalos_orfan= gsl_rng_uniform_int(gBaseRand,ny); //draw number counts in cells
                    prob = this->BIAS_NCOUNTS_normalized[index_2d(Nhalos_orfan,index, lenght_bias)];
                    ran = gsl_rng_uniform(gBaseRand);
                 }
                 real_prec sigma_pdf= sqrt(static_cast<real_prec>(pdf_ref[Nhalos_orfan]));
                 real_prec dif_pdf = static_cast<real_prec>(pdf_ref[Nhalos_orfan])-static_cast<real_prec>(pdf_temp[Nhalos_orfan]+1);

                if(dif_pdf >0 )
                   {
                     this->delta_Y_new[i]= static_cast<real_prec>(Nhalos_orfan);
                     cell_info[i].assigned=true;
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                    pdf_temp[Nhalos_orfan]++;
                    counter_orphan_b++;
                  }
            }
         }
     So.message_screen("Second loop done");
*/
     So.message_screen("Start third loop"); // this is the same second part of the mian loop, where the main prob dist is used.

      // One option here is to collect the cells that are not yet assigned and randomize them
     // Then assign number of tracers in cells until the number of tracers reaches a value close to the reference * Some Poisson dispersion.


     vector<ULONG>i_nac;
     for(ULONG i=0;i<this->NGRID;++i)
       if(false==cell_info[i].assigned)
         i_nac.push_back(i);
     gsl_ran_shuffle(gBaseRand,&i_nac[0], i_nac.size(),sizeof(ULONG));


     //Start loop over the cells that are not yeat assigned
     ULONG Ntracer_new=this->Nobjects;
     ULONG Ntracer_ref_poiss = static_cast<ULONG>(gsl_ran_poisson(gBaseRand,Nobjects_ref));
     So.message_screen("Setting threshold at ",Ntracer_ref_poiss );
     So.message_screen("NUmber of reference tracers",Nobjects_ref);

     for(ULONG i=0;i<i_nac.size();++i) // loop ove rrandomized not assigned cells
       {
        if(Ntracer_new<Ntracer_ref_poiss)
         {
          ULONG i_cell = i_nac[i]; // original ID of the not assigned cell
          ULONG index= cell_info[i_cell].theta_index;
          int Nhalos_orfan=0;
          real_prec prob=-1.0;
          real_prec ran=0.0;
          while(prob<ran) // This while aims at chosing the right Nh according to the full propability distribution
            {
              Nhalos_orfan= gsl_rng_uniform_int(gBaseRand,ny); //draw number counts in cells
              prob = this->BIAS_NCOUNTS_normalized[index_2d(Nhalos_orfan,index, lenght_bias)];
              ran = gsl_rng_uniform(gBaseRand);
           }
           counter_orphan++;
           Ntracer_new+=Nhalos_orfan;
           this->delta_Y_new[i_cell]= static_cast<real_prec>(Nhalos_orfan);
        }
      } 
      i_nac.clear();i_nac.shrink_to_fit();

/*
#ifdef _USE_OMP_ 
#pragma omp for reduction(+:counter_orphan)
#endif
     //Start loop over the cells
     for(ULONG i=0;i<this->NGRID;++i)
       {
        if(false==cell_info[i].assigned)
          {
            ULONG index=cell_info[i].theta_index;
            int Nhalos_orfan;
            real_prec prob=-1.0;
            real_prec ran=0.0;
            while(prob<ran) // This while aims at chosing the right Nh according to the full propability distribution
              {
                 Nhalos_orfan= gsl_rng_uniform_int(gBaseRand,ny); //draw number counts in cells
                 prob = this->BIAS_NCOUNTS_normalized[index_2d(Nhalos_orfan,index, lenght_bias)];
                 ran = gsl_rng_uniform(gBaseRand);
              }
              this->delta_Y_new[i]= static_cast<real_prec>(Nhalos_orfan);
              counter_orphan++;
           }
      }
*/

     So.message_screen("Third loop done");

    cell_info.clear();cell_info.shrink_to_fit();
#endif // end_new_app_



#ifdef _USE_OMP_
     gsl_rng_free (gBaseRand);
     }// end parallelized region
#endif


    X_Y_hist_dynamical.clear();X_Y_hist_dynamical.shrink_to_fit();
    X_Y_hist_dynamical_normalized.clear();X_Y_hist_dynamical_normalized.shrink_to_fit();

#ifdef _new_counter_
    new_X_Y_hist_dynamical.clear();new_X_Y_hist_dynamical.shrink_to_fit();
    new_X_Y_hist_dynamical_normalized.clear();new_X_Y_hist_dynamical_normalized.shrink_to_fit();
#endif





    this->NCELLSperDMBIN.clear(); this->NCELLSperDMBIN.shrink_to_fit();
    NCELLSperDMBIN_now.clear();    NCELLSperDMBIN_now.shrink_to_fit();
    Ncells_density_bin_new.clear(); Ncells_density_bin_new.shrink_to_fit();


//     testa.close();
#ifdef _FULL_VERBOSE_
   So.message_time_mock(start_mock);
#endif
   // ******************************************** end assigning ncounts or other prop*******************************

#ifdef _FULL_VERBOSE_
 if(counter_orphan>0)
   So.message_screen("Fraction of cells assigned with global prob. dist =", 100.0*static_cast<double>(counter_orphan)/static_cast<double>(this->NGRID), "%");
 if(counter_orphan_b>0)
   So.message_screen("Fraction of cells assigned with global prob. dist conditional to pdf =", 100.0*static_cast<double>(counter_orphan_b)/static_cast<double>(this->NGRID), "%");
 if(counter_out_of_theta>0)
     So.message_screen("Fraction of cells without repesentation in the ref bias=", 100.0*static_cast<real_prec>(counter_out_of_theta)/static_cast<real_prec>(this->NGRID), "%");

#endif
   if(property==_COUNTS_)
     {

       this->Nobjects=static_cast<int>(get_nobjects(this->delta_Y_new));
#ifdef _FULL_VERBOSE_
       So.message_screen("Number of objects in new mock= ", this->Nobjects);
       So.message_screen("Mean number density = ", static_cast<real_prec>(this->Nobjects)/pow(Lbox, 3), "(Mpc / h )⁻³");
#endif


#ifdef _DO_BAM_CALIBRATION_
#ifndef _EXTRAPOLATE_VOLUME_


#ifdef _FULL_VERBOSE_

       if(this->Nobjects < this->N_objects_Y)
         So.message_screen("->Less objects assigned to mock. Fraction = ", fabs(100.0-100.0*(static_cast<real_prec>(this->Nobjects)/this->N_objects_Y)), "%");
       else if(this->Nobjects  > this->N_objects_Y)
         So.message_screen("->More objects assigned to mock. Fraction = ", fabs(100.0-100.0*(static_cast<real_prec>(this->Nobjects)/this->N_objects_Y)), "%");
#endif

#endif
#endif


       // ----------------------------------------------------------------------------------------------
       if(false==silent)
         {
           real_prec anmean=static_cast<real_prec>(this->Nobjects)/static_cast<real_prec>(this->NGRID);
#ifdef _FULL_VERBOSE_
           So.message_screen("<N> objects in mock =", anmean);
#endif

           int nyy=30;
           real_prec dmax=300.0;

           vector<real_prec>DELTA_HIST(nyy,0);

           real_prec delta_delta= (dmax+1.0)/static_cast<real_prec>(nyy);  //deltamax - deltamin
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
           for(ULONG i=0;i< this->NGRID;++i)
             {
               real_prec delta_mock= (this->delta_Y_new[i]/static_cast<real_prec>(anmean) - 1.0);
               if(delta_mock>=-1.0 && delta_mock<=dmax)
                 {
                   int ID=static_cast<int>(floor(delta_mock +1.0)/delta_delta);
                   if(ID==nyy)ID--;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                   DELTA_HIST[ID]++;
                 }
             }

           real_prec lkk=get_max<real_prec>(DELTA_HIST);
           for(int j=0;j< nyy ;++j)
             DELTA_HIST[j]= lkk == 0  ? 0.0 :  static_cast<real_prec>(DELTA_HIST[j])/static_cast<real_prec>(lkk);

           ofstream aja;
           aja.open("delta_mock_dist.txt");
           for (int i=0;i<DELTA_HIST.size();++i)aja<<i<<"\t"<<DELTA_HIST[i]<<endl;
           aja.close();
           DELTA_HIST.clear();
           DELTA_HIST.shrink_to_fit();

           vector<real_prec>AUX(this->NGRID,0);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
           for(ULONG i=0; i< this->NGRID ;++i)
             AUX[i]=this->delta_Y_new[i]/static_cast<real_prec>(anmean) - 1.0;

           lkk=get_max<real_prec>(AUX);
           So.message_screen("Maximum delta in mock =", lkk);

           lkk=get_min<real_prec>(AUX);
           So.message_screen("Minimum delta in mock =", lkk);
           AUX.clear();
           AUX.shrink_to_fit();

        }

           // ----------------------------------------------------------------------------------
       if(true==this->Write_PDF_number_counts) 	  // Write the PDF of the outcome */
         {
               string fileY=this->Output_directory+"PDF_NC"+"_Y_MOCK"+this->new_Name_Property_Y+"_MASY"+to_string(this->iMAS_Y)+"_Nft"+to_string(this->Nft)+"_z"+to_string(this->redshift)+"_"+this->params._file_power()+".txt";

               So.message_screen("Writting PDF_NC for MOCKS");
               int NPART=ny; // Number of particles in cells

               this->PDF_NC_Y.clear();
               this->PDF_NC_Y.resize(NPART, 0);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
               for(ULONG i=0;i< this->NGRID;++i)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 this->PDF_NC_Y[static_cast<int>(this->delta_Y_new[i])]++;  // Get pdf of tracer number counts


               this->File.write_to_file_i(fileY,this->PDF_NC_Y);// Write pdf to file

               vector<int> cells_with_one;
               for(int i=0;i<this->PDF_NC_Y.size();++i)
                 if(this->PDF_NC_Y[i]>0)
                   cells_with_one.push_back(i);

               So.message_screen("Maximum number of mock-objects in one cell =", get_max<int>(cells_with_one) );

               cells_with_one.clear();
               this->PDF_NC_X.clear();
               cells_with_one.shrink_to_fit();
               this->PDF_NC_X.shrink_to_fit();
	        }
       }

   else if (_DENSITY_==property)
     {
       // if we have Y as density, then we have to sample the bins in Y to rweturn values of a Y field
       gsl_rng_env_setup();
       gsl_rng_default_seed=236655;//time(NULL);
       const gsl_rng_type *  Tn= gsl_rng_ranlux;
       gsl_rng * rn = gsl_rng_alloc (Tn);

#ifdef _USE_OMP_
#pragma omp for
#endif
       for(ULONG i=0;i< this->NGRID;++i) //tbc
	 {  /// cambiar ldelta_Y por le mínimo de delta_Y_new
	   real_prec xr=gsl_rng_uniform (rn);
           real_prec aux_den=this->s_mins.prop0+(static_cast<int>(delta_Y_new[i])+xr)*this->s_deltas.prop0;   // Get value of log10(2+delta_y)
	   this->delta_Y_new[i]=this->Mean_density_Y*(1.0+pow(10,aux_den)-NUM_IN_LOG);                                        // Get nbar*(1+delta_y)
	 }

     }



#ifdef  _EXTRAPOLATE_VOLUME_
#ifdef _WRITE_TR_DENSITY_FIELD_
   this->File.write_array(fname, this->delta_Y_new);
#endif
#endif


   // this is commented for only was applicable to DENSITY case. For the case in which we calibrate ncounts
   // and other properties, those extra fields are written directly as they come in the last iteration
   if(case_prop==4)
     {
       gsl_rng_env_setup();
       gsl_rng_default_seed=236655;//time(NULL);
       const gsl_rng_type *  Tn= gsl_rng_ranlux;
       gsl_rng * rn = gsl_rng_alloc (Tn);

#ifndef _USE_LOG_MASS_


       real_prec mass_min=pow(10,params._LOGMASSmin())*params._MASS_units();
        So.message_warning_ini(__LINE__,__PRETTY_FUNCTION__,__FILE__,"Calibrating from sat fraction. Please check definition of containers vector<int>ncounts and vector<bool>filled_cells");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
        for(ULONG i=0;i< this->NGRID;++i)
         {
           real_prec xr=gsl_rng_uniform (rn);
           real_prec aux=this->s_mins.prop0_mass+(static_cast<int>(delta_Y_new[i])+xr)*this->s_deltas.prop0_mass;   // Get value of log10(2+delta_y)
           real_prec mass_cell = (pow(10,aux)-NUM_IN_LOG)*MASS_SCALE;
           real_prec filled=1;
           int ncounts_m=1;
#ifdef _use_filled_cells_
           filled=static_cast<real_prec>(filled_cells[i]);
           ncounts_m=ncounts[i];
#endif
           this->delta_Y_new[i]=mass_cell*filled;

           int Neff = static_cast<int>(floor(this->delta_Y_new[i]/mass_min));

           // This weights help to upweight the mass at a cell such that, give the number counts it has, it can at least provide mass for all particles
           real_prec weight = (Neff >= ncounts_m ? 1.0 : static_cast<real_prec>(ncounts_m)/static_cast<real_prec>(Neff));
           this->delta_Y_new[i]*=weight;

         }
       gsl_rng_free (rn);
       So.message_screen("Minimum of new mass field = ", get_min(this->delta_Y_new));
       So.message_screen("Maximim of new mass field = ", get_max(this->delta_Y_new));
#endif
     }
#ifdef _use_filled_cells_
   ncounts.clear();
   ncounts.shrink_to_fit();
#endif

   // Writ ethe density fields of properties of mock tracers:.


#ifndef _GET_BAM_REALIZATIONS_
#ifndef _TEST_THRESHOLDS_RESIDUALS_
   auto out_it = std::find(std::begin(this->output_at_iteration), std::end(this->output_at_iteration), this->step);
//   if((this->step==this->N_iterations_Kernel) || (out_it == std::end(this->output_at_iteration)))
//    if(this->step==this->N_iterations_Kernel)
     this->File.write_array(fname, this->delta_Y_new);
#endif
#endif


#ifdef _GET_BAM_REALIZATIONS_
   this->File.write_array(fname, this->delta_Y_new);
   if(case_prop==1 || case_prop==2)
     this->patchy.fname_MOCK_NCOUNTS=fname;
   if(case_prop==4)
     this->patchy.fname_MOCK_MASS=fname;
   if(case_prop==3)
     this->patchy.fname_MOCK_NCOUNTS_SAT=fname;

#endif

   // Get power spectrum only number counts

#ifdef _DISPLACEMENTS_
    this->File.read_array_t<PrecType_Y>(this->Input_Directory_Y+this->Name_Catalog_Y, this->delta_Y);  // read the reference
    this->delta_Y_new=this->delta_Y_new+this->delta_Y;  // Construct psi new = delta Psi + psi_ref
#endif


   if(case_prop==1 || case_prop==2)
   {
#ifdef _GET_BAM_REALIZATIONS_
       this->get_power_spectrum("TR_MOCK_REALIZATION");
    this->delta_Y.resize(this->NGRID,0);
    this->File.read_array_t<PrecType_Y>(this->Input_Directory_Y+this->Name_Catalog_Y, this->delta_Y);
    this->get_power_spectrum("TR_REF");
    this->delta_Y.clear();this->delta_Y.shrink_to_fit();
     vector<pair<real_prec, real_prec> > xy_pts_ref;
     vector<pair<real_prec, real_prec> > xy_pts_new;
     for(int i=0; i<kvec.size(); ++i) 
       xy_pts_ref.push_back(std::make_pair(this->kvec[i], log10(this->Power_REF[i])));
     for(int i=0; i<kvec.size(); ++i) 
       xy_pts_new.push_back(std::make_pair(this->kvec[i], log10(this->Power_NEW[i])));
     this->gp_power<<"set log x \n";
     this->gp_power<<"set border linewidth 2.2\n";
     this->gp_power<<"set title 'Realization:"<<this->params._realization()<<". Reference: "<<this->params._seed()<<"'\n";
     this->gp_power<<"set size square \n";
     this->gp_power<<"set xlabel 'k [h / Mpc]' font 'Times-Roman,14'\n";
     this->gp_power<<"set ylabel 'log P(k) [(Mpc / h)³]'  font 'Times-Roman,12'\n";
#ifdef _USE_GNUPLOT_POWER_PLOT_    // PLOT OF THE POWER SPECTRUM
     this->gp_power<<"plot"<<this->gp_power.file1d(xy_pts_ref) << "w l lw 2.3 lt 8 title 'Reference',"<<this->gp_power.file1d(xy_pts_new) << " w l lw 2.3 lt 6 title 'Mock',"<<endl;
     this->gp_power<<"plot"<<this->gp_power.file1d(xy_pts_ref) << "w l lw 2.3 lt 8 title 'Reference',"<<this->gp_power.file1d(xy_pts_new) << " w l lw 2.3 lt 6 title 'Mock',"<<endl;
#endif
     this->gp_power<<"set terminal pngcairo\n"; 
     this->gp_power<<"set output '"<< this->params._Output_directory()<<"power_comparison_ref"<<this->params._seed()<<"_realization"<<this->params._realization()<<".png'"<<endl;
     this->gp_power<<"replot\n";
     this->gp_power<<"clear\n";
     xy_pts_ref.clear();
     xy_pts_ref.shrink_to_fit();
     xy_pts_new.clear();
     xy_pts_new.shrink_to_fit();



#else

       /* One more test unsucesfullyy tried to solve the SLICS issue
       real_prec aux_mean_Y;
       gsl_rng_env_setup();
       gsl_rng_default_seed=time(NULL);
       const gsl_rng_type *  Tp= gsl_rng_ranlux;
       gsl_rng * rp = gsl_rng_alloc (Tp);

       get_overdens(this->delta_Y_new,  aux_mean_Y, this->delta_Y_new);  // get overdensity
       So.message_screen("Converting DENSITY -> Poisson(DENSITY):");
       for(ULONG i=0;i<delta_Y_new.size();++i)
          delta_Y_new[i]=gsl_ran_poisson(rp,aux_mean_Y*(1.+delta_Y_new[i]));   //recover density
       So.DONE();
*/
       this->get_power_spectrum("TR_MOCK");
#endif
   }



   this->shot_noise_new = pow(this->Lbox,3)/static_cast<real_prec>(this->Nobjects);
   this->shot_noise_ref=this->shot_noise_new;  // This should be the case during the calibration


   this->delta_Y_new.clear();
   this->delta_Y_new.shrink_to_fit();

 }

#endif



 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################
 //  ####################################################################################################################################################################



#ifdef _GET_BAM_CAT_

#ifdef _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_


 // This property assignment explicitely uses the values of the properties taken from the reference and assign them to the mocks provided the underlying properties of the DM field.
 // For Initial assignmetn = true, the initial assignment is done over a property Vmax (which seems to be that correlating the most with DM)
 // Initial_assignment = false, the mass is to be assigned based on the assignment of Vmax done in a previous call of this function
 // This function is used for Vmax assignemnet
 // This function must be in accordance with what Bam::get_X_function() has done


 void Bam::assign_tracer_property_new_new(bool initial_assignment, string h_property)
{

  So.enter(__PRETTY_FUNCTION__);


#ifdef _FULL_VERBOSE_
  So.message_screen("**************************************************");
  if(initial_assignment==true)
#ifdef _USE_VMAX_AS_OBSERVABLE_
      So.message_screen("*Assigning Vmax from reference values to mock");
#elif defined _USE_MASS_AS_OBSERVABLE_
      So.message_screen("*Assigning Halo Mass from reference*");
#endif
  So.message_screen("Number of tracers in reference = ", this->tracer_ref.NOBJS);
  So.message_screen("Number of tracers in mock      = ", this->tracer.NOBJS);
  So.message_screen("**************************************************");
#endif

  ULONG cumulative_counter=0; // this counts all

  // ********************************************************************************************************************
  // *********************************************OMP stuff***********************************************************************
  int NTHREADS=_NTHREADS_;
#ifdef _FULL_VERBOSE_
  So.message_screen("Using ",NTHREADS," threads");
#endif
  omp_set_num_threads(NTHREADS);


  // These are to be used out of parallel regions
  gsl_rng_env_setup();
  const gsl_rng_type *Tn;
  gsl_rng_default_seed=10016515;
  Tn =  gsl_rng_mt19937 ;//gsl_rng_default;
  gsl_rng *rn = gsl_rng_alloc (Tn);

  // ********************************************************************************************************************
  ULONG LENGHTdm=this->dm_properties_bins.size();
  //Vector of structures aimed at keeping the values of the tracer properties (vmax, or halos) which fall in a bin of {Theta}
  //This information will be modified below, so we make a copy ni order not to measure that again above in get_X_function()


#ifndef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_
  vector<ULONG>number_in_theta_ref(LENGHTdm,0);
  // container to track the number of available masses in each theta-bin
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<LENGHTdm; ++i)
    number_in_theta_ref[i]=this->dm_properties_bins[i].masses_bin_properties.size();


#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_

  ULONG LENGHTdm_for_randoms=this->dm_properties_for_randoms_bins.size(); //In principe this number is equal to this->dm_properties_bins.size();
  vector<ULONG>number_in_theta_ref_for_randoms(LENGHTdm_for_randoms,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<LENGHTdm_for_randoms; ++i)
    number_in_theta_ref_for_randoms[i]=this->dm_properties_for_randoms_bins[i].masses_bin_properties.size();
#endif


#endif  // end of ifndef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_

  // ********************************************************************************************************************
  // ********************************************************************************************************************
  // *******************************Mass, Vmax function stufff************************************************************************************
  int nbins_mf=0;
  real_prec lm_min, lm_max;
  gsl_interp_accel *acc;
  gsl_spline *spline ;

#ifdef _USE_VMAX_AS_OBSERVABLE_
  if(true==initial_assignment)
    nbins_mf=this->tracer_ref.vmax_function.size();
  else if(h_property==_MASS_)
     nbins_mf=this->tracer_ref.mass_function.size();
  else if (h_property==_RS_)
     nbins_mf=this->tracer_ref.rs_function.size();
  else if (h_property==_SPIN_)
     nbins_mf=this->tracer_ref.s_function.size();

#elif defined _USE_MASS_AS_OBSERVABLE_
      nbins_mf=this->tracer_ref.mass_function.size();
#endif

   if(nbins_mf==0)
     {
       So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Reference mass function has not be allocated and perhaps not measured. Code exits here.");
       exit(0);
     }
#ifdef _FULL_VERBOSE_
   So.message_screen("Preparing arrays for global X-function from reference: ");
#endif
   this->tracer.define_property_bins();
   this->mfunc.resize(nbins_mf,0);

#ifdef _USE_VMAX_AS_OBSERVABLE_
     if(true==initial_assignment)
      {
        this->prop_min=this->tracer_ref.VMAXBin[0];
        this->prop_max=this->tracer_ref.VMAXBin[nbins_mf-1];
     }
     else
       {
        if(h_property==_MASS_)
         {
           this->prop_min=this->tracer_ref.MBin[0];
           this->prop_max=this->tracer_ref.MBin[nbins_mf-1];
         }
        if(h_property==_RS_)
           {
             this->prop_min=this->tracer_ref.RSBin[0];
             this->prop_max=this->tracer_ref.RSBin[nbins_mf-1];
           }
         else if (h_property==_SPIN_)
            {
             this->prop_min=this->tracer_ref.SPINBin[0];
             this->prop_max=this->tracer_ref.SPINBin[nbins_mf-1];
             }

      }
#elif defined _USE_MASS_AS_OBSERVABLE_
      this->prop_min=this->tracer_ref.MBin[0];
      this->prop_max=this->tracer_ref.MBin[nbins_mf-1];
#endif

      // I cannot use here the this->tracer MBmin and max for it might happen that the reference mass function has been measured with a different number of bins
      // different to the current this>tracer_ref_NMBINS
      vector<real_prec>delta_prop_aux(nbins_mf,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0; i< nbins_mf; ++i)
        {
          real_prec lmmin=log10(this->prop_min)+i*log10(this->prop_max/this->prop_min)/static_cast<real_prec>(nbins_mf);
          real_prec lmmax=log10(this->prop_min)+(i+1)*log10(this->prop_max/this->prop_min)/static_cast<real_prec>(nbins_mf);
          delta_prop_aux[i]=pow(10,lmmax)-pow(10,lmmin);
        }

#ifdef _USE_VMAX_AS_OBSERVABLE_
      if(true==initial_assignment)
        {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(int i=0; i< nbins_mf; ++i)
            this->mfunc[i]=this->tracer_ref.vmax_function[i]*delta_prop_aux[i];
        }
        else
          if(h_property==_MASS_)
            {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
             for(int i=0; i< nbins_mf; ++i)
              this->mfunc[i]=this->tracer_ref.mass_function[i]*delta_prop_aux[i];
        }
        else if(h_property==_RS_)
          {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
            for(int i=0; i< nbins_mf; ++i)
              this->mfunc[i]=this->tracer_ref.rs_function[i]*delta_prop_aux[i];
          }
        else if(h_property==_SPIN_)
          {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(int i=0; i< nbins_mf; ++i)
             this->mfunc[i]=this->tracer_ref.s_function[i]*delta_prop_aux[i];
        } // closes  if(true==initial_assignment)

#elif defined _USE_MASS_AS_OBSERVABLE_

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0; i< nbins_mf; ++i)
        this->mfunc[i]=this->tracer_ref.mass_function[i]*delta_prop_aux[i];
#endif

      delta_prop_aux.clear(); delta_prop_aux.shrink_to_fit();

      real_prec max_mf=static_cast<real_prec>(get_max(mfunc));
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0; i< nbins_mf; ++i)
        {
          real_prec mf=this->mfunc[i];
          this->mfunc[i]=mf/max_mf;
        }

      acc = gsl_interp_accel_alloc ();

#ifdef _USE_VMAX_AS_OBSERVABLE_
    if(true==initial_assignment) // for true, we initially assign vmax
      {
        spline    = gsl_spline_alloc (gsl_interp_linear, this->tracer_ref.vmax_function.size());
        gsl_spline_init (spline, &(this->tracer_ref.VMAXBin[0]), &(this->mfunc[0]), this->tracer_ref.vmax_function.size());
        lm_min=log10(this->params._VMAXmin());
        lm_max=log10(this->params._VMAXmax());
      }
      else
       {
        if(h_property==_MASS_)
         {
           spline    = gsl_spline_alloc (gsl_interp_linear, this->tracer_ref.mass_function.size());
           gsl_spline_init (spline, &(this->tracer_ref.MBin[0]), &(this->mfunc[0]), this->tracer_ref.mass_function.size());
           lm_min=this->params._LOGMASSmin();
           lm_max=this->params._LOGMASSmax();
         }
         else if(h_property==_RS_)
          {
           spline    = gsl_spline_alloc (gsl_interp_linear, this->tracer_ref.rs_function.size());
           gsl_spline_init (spline, &(this->tracer_ref.RSBin[0]), &(this->mfunc[0]), this->tracer_ref.rs_function.size());
           lm_min=log10(this->params._RSmin());
           lm_max=log10(this->params._RSmax());
          }
         else if(h_property==_SPIN_)
          {
           spline    = gsl_spline_alloc (gsl_interp_linear, this->tracer_ref.s_function.size());
           gsl_spline_init (spline, &(this->tracer_ref.SPINBin[0]), &(this->mfunc[0]), this->tracer_ref.s_function.size());
           lm_min=log10(this->params._SPINmin());
           lm_max=log10(this->params._SPINmax());
          }

        }
#elif defined _USE_MASS_AS_OBSERVABLE_
      spline    = gsl_spline_alloc (gsl_interp_linear, this->tracer_ref.mass_function.size());
      gsl_spline_init (spline, &(this->tracer_ref.MBin[0]), &(this->mfunc[0]), this->tracer_ref.mass_function.size());
      lm_min=this->params._LOGMASSmin();
      lm_max=this->params._LOGMASSmax();
#endif

  // ********************************************************************************************************************
  // ************************************************RENAME SOME BIN PROPPERTIES*****************************************
  // ********************************************************************************************************************

  int N_a=N_C_BIN1;


  if(true==initial_assignment) // This applies for the assignment of Vmax with the different techniques
    {

#ifdef _USE_TOTAL_MASS_IN_CELL_  //check this, things need to be fixed inside, related to paths

      vector<real_prec> MOCK_MASS_FIELD;
      N_a = N_BINS_TOTAL_MASS_IN_CELL;
      // Can happen that the used mock properties fall in a theta bin in which the reference has nothing, hence  dm_properties_bins[index_bins].masses_bin_properties.size()=0
      // This has been visible when using a total mass density field sampled from the calibration (as it must be) instead of using the reference
      // This is likely to happen if we do not use the sampled total mass field that correspond to the final bumber density field that one uses to assign positions
      // Hence, if we are assigning as test_mode to the reference, use the mass density field from the reference
      MOCK_MASS_FIELD.resize(this->NGRID,0);
#ifdef _ASSIGN_TO_REFERENCE_ 
      So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"I am passing the mass density field from the reference, in order toto do tests. Ideally this should be the one sampled in the last step of the calibration");
      this->tracer_ref.get_density_field_grid(_MASS_, MOCK_MASS_FIELD); // ojo que estoy cargando la ref, no el mock:esto debería ser el mock con la masa total dada por la calibración
#else
      File.read_array(this->Output_directory+"MOCK_TR_MASS_iteration200_MASY0_Nft256_z1.124.dat",MOCK_MASS_FIELD);
#endif

      real_prec mmax=get_max<real_prec>(MOCK_MASS_FIELD);
#ifdef _FULL_VERBOSE_
      So.message_screen("Maximum mass of tracer in mock-cell", mmax);
#endif
      mmax=get_min<real_prec>(MOCK_MASS_FIELD);
#ifdef _FULL_VERBOSE_
      So.message_screen("Minimum mass of tracer in mock-cell", mmax);
#endif
#endif //end totalmass in cell

    }

  // ************************************************RENAME SOME BIN PROPPERTIES********************************************************************

  int N_b=N_C_BIN2;
  if(true==initial_assignment)
    {
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
      this->tracer.get_neighbour_tracers(this->ncells_info);
      N_b = N_NEIGHBOURS_MAX;
#elif defined (_USE_MIN_DISTANCE_TO_NEIGHBOURS_) || defined (_USE_LOCAL_CLUSTERING_)
      N_b=N_BINS_MIN_DIST_TO_NEI;
#endif
    }
  // ************************************************RENAME SOME BIN PROPPERTIES********************************************************************
  int N_c = N_C_BIN3;
#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
  real_prec delta_min_sep=0;
#endif

  if(true==initial_assignment)  // Ww shall not use his info for post-mass assignment
    {
#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
      this->tracer.get_min_separation_in_cell();
      N_c = N_BINS_MIN_SEP_IN_CELLS;
      delta_min_sep = (MAX_SEP_IN_CELLS-this->min_halo_separation)/static_cast<real_prec>(N_BINS_MIN_SEP_IN_CELLS);
#endif
    }

  // ********************************************************************************
  int N_v= N_CV_BIN1;

#ifdef _USE_TRACERS_IN_CELLS_
  vector<real_prec> MOCK_DEN_FIELD;
#endif

  if(true==initial_assignment)  // Ww shall not use his info for post-mass assignment
    {
#ifdef _USE_TRACERS_IN_CELLS_
      N_v = N_BINS_TRACERS_IN_CELLS;
      MOCK_DEN_FIELD.resize(this->NGRID,0);
      this->File.read_array(this->fnameMOCK,MOCK_DEN_FIELD); // fname mock has been defined in bamrunner before calling makecat

//     vector<real_prec> MOCK_DEN_FIELD;
//      MOCK_DEN_FIELD.resize(this->NGRID,0);
//#ifdef _FULL_VERBOSE_
//      So.message_screen("Using number of tracers in cell");
//#endif

//     this->tracer.get_density_field_grid(_COUNTS_, MOCK_DEN_FIELD);
      int nmax=get_max<real_prec>(MOCK_DEN_FIELD);
#ifdef _FULL_VERBOSE_
      So.message_screen("Maximum number of tracer in cell", nmax);
#endif

#endif // end USE_TRACERS_IN_CELLS_
    }
  // ********************************************************************************
  // ********************************************************************************************************************
  // ********************************************************************************************************************
  // ********************************************************************************************************************

#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_
  // Number of tracers with properties in each of the multi-scale levels defined
  ULONG N_props_0=0;
  ULONG N_props_1=0;
  ULONG N_props_2=0;
  ULONG N_props_3=0;
  ULONG N_props_4=0;
#endif

  int N_x= N_CV_BIN2;

#ifdef _USE_VMAX_AS_OBSERVABLE_
#ifdef _ASSIGN_MASS_POST_
  if(false==initial_assignment) // if initial_assignem t is false, we proceed to assign mass once the vmax is already assigned. This will be false upon a second call of this function.
    N_x = N_VMAX_BINS;
#endif
#endif

  ULONG counter_fmf=0;

  // ********************************************************************************************************************
  // ********************************************************************************************************************
#if defined _USE_VMAX_AS_OBSERVABLE_ || defined _USE_MASS_AS_OBSERVABLE_ 
  if(true==initial_assignment) // This applies for the assignment of vmax or mass, reading from the reference catalogs, whether we are using or not multi-scale approach
    {
#endif

#if defined (_USE_MULTISCALE_PROPERTY_ASSIGNMENT_) || defined (_SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_)
// This container of structures can in principle hold the same info as the vector tracer[].galThetaID.
// but since it is useful when loops are performed over the grid (where we tracer[].galThetaID is not accesible)
// we use this:
    vector<s_cell_info> cell_info_tr(this->NGRID);
#endif

// ************************************************************************************
// ************************************************************************************
// ************************************************************************************
//  The following loop aims at collecting the values of Vmax in theta bins.
// ************************************************************************************
// ************************************************************************************



  counter_fmf=0;
  ULONG count_dm=0;

#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_  //this is only done for vmax so far, since we shall  use that property as main property
  ULONG count_dm_full=0;
  ULONG count_ran_full=0;
#endif

  // Note:
  // If the number of neighbours is to be used, we do a loop over the particles.
  // Also, when multiscale *is not* used, we do a loop over the particles.
  // Else (e.g. using multiscale or using neighbour info), we do a loop over the grid cells, since the some of
  // the used quantities are computed at the cells.
  // In both cases, parallelization can be done, taking care of defining the random objects inside parallel region


      // Start parallel region. Problems with this paralellization
#ifdef _USE_OMP_TEST_

  int jthread=0;
  const gsl_rng_type *Trn;
  gsl_rng *randa;
  vector<ULONG>vseeds(NTHREADS,0);
  for(int i=0;i<vseeds.size();++i)
    vseeds[i]=35+static_cast<ULONG>(i+14)*56045;


#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_
#ifdef _FULL_VERBOSE_
  this->So.message_screen("Collecting info for multuscaling");
#endif
#elif defined _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#ifdef _FULL_VERBOSE_
  this->So.message_screen("Assinging with threhold in Vmax");
#endif
#endif

#pragma omp parallel private (jthread, Trn, randa)// this is causing problem
  {
    jthread=omp_get_thread_num();
    gsl_rng_default_seed=vseeds[jthread];
    Trn = gsl_rng_mt19937;//_default;
    randa= gsl_rng_alloc (Trn);
#endif

#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_


#ifdef _USE_OMP_TEST_
#pragma omp for reduction(+:count_dm)
#endif
#if defined _USE_NUMBER_OF_NEIGHBOURS_ || defined _ASSIGN_MASS_POST_
    for(ULONG ig=0; ig < this->tracer.NOBJS;++ig)
#else
    for(ULONG id=0; id < this->NGRID;++id)
#endif

#else   //else for _USE_MULTISCALE_PROPERTY_ASSIGNMENT_




     
      ULONG  count_test=0;
      for(ULONG i=0;i<this->dm_properties_bins.size();++i)
          for(ULONG j=0;j<this->dm_properties_bins[i].used_mass.size();++j)
            if(false==this->dm_properties_bins[i].used_mass[j])
              count_test++;
      cout<<"Check Previous to assignment. Remaining number of properties in theta_dm = "<<count_test<<endl;
      cout<<this->tracer.NOBJS<<endl;
      count_test=0;

/*
      count_test=0;
      for(ULONG i=0;i<this->dm_properties_for_randoms_bins.size();++i)
          for(ULONG j=0;j<this->dm_properties_for_randoms_bins[i].used_mass.size();++j)
            if(false==this->dm_properties_for_randoms_bins[i].used_mass[j])
              count_test++;
      cout<<"remaining theta_dm_for_randoms "<<count_test<<endl;
*/



// andres
//  cout<<RED<<this->delta_X[1255]<<"   "<<this->delta_X.size()<<"   "<<this->NGRID<<RESET<<endl;
  cout<<RED<<this->delta_X[12255]<<"   "<<cwclass.Invariant_TF_II[12255]<<"  "<<cwclass_ref.Invariant_TF_III[12255]<<endl;


#ifdef _USE_OMP_TEST_
#pragma omp for reduction(+:counter_fmf, count_dm)
#endif
    for(ULONG ig=0; ig < this->tracer.NOBJS; ++ig)
#endif    // endif for _USE_MULTISCALE_PROPERTY_ASSIGNMENT_
    	{

        

       // Get the cell ID where the mock tracer is located
#if defined _USE_MULTISCALE_PROPERTY_ASSIGNMENT_ || defined (_SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_)
#if defined _USE_NUMBER_OF_NEIGHBOURS_ || defined _ASSIGN_MASS_POST_
	      ULONG id=this->tracer.Halo[ig].GridID;
#endif
#else

        // ANDres
        ULONG id=this->tracer.Halo[ig].GridID;
#endif

        int I_X=0;
       // Get the bin in the delta dark matter (or log 1+delta) in each cell
        real_prec xdm = static_cast<real_prec>(this->delta_X[id]);
        I_X  = get_bin(xdm,this->s_mins.prop1,this->NX,this->s_deltas.prop1,this->bin_accumulate_borders);

        int I_CWT=0;
#ifdef _USE_CWC_
    	  I_CWT=this->cwclass.get_Tclassification(id);
#endif

	  // Get the corresponding bin on Knot-mass. The current cell is in a cluster with a given mass falling in that particular bin
       int I_MK=0;
#ifdef _USE_MASS_KNOTS_
       I_MK= (this->cwclass.cwt_used[I_CWT]== I_KNOT ? this->cwclass.SKNOT_M_info[id]: 0);
#endif

       int I_CWV=0;
#ifdef _USE_CWC_V_
       I_CWV=this->cwclass.get_Vclassification(id);
#endif

       // Get the corresponding bin on Knot-mass. The current cell is in a cluster with a given mass falling in that particular bin
       int I_VK=0;
#ifdef _USE_VEL_KNOTS_V_
       I_VK= (this->cwclass.cwt_used[I_CWV]== I_KNOT ? this->cwclass.VDISP_KNOT_info[id]: 0);
#endif

       // Get the corresponding bin in the two invariants of the shear of the tidal field
       int I_C1=0;
#ifdef _USE_TOTAL_MASS_IN_CELL_
       if(true==initial_assignment)
         I_C1=get_bin(log10(MOCK_MASS_FIELD[id]),this->params._LOGMASSmin(),N_BINS_TOTAL_MAS  S_IN_CELL,(log10(MAX_TOTAL_MASS_IN_CELL)-this->params._LOGMASSmin())/static_cast<double>(N_BINS_TOTAL_MASS_IN_CELL), true);

#elif defined _USE_INVARIANT_TIDAL_FIELD_II_
       real_prec C1 = this->cwclass.Invariant_TF_II[id];
       I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#elif defined _USE_DELTA2_
       real_prec C1 = this->cwclass.DELTA2[id];
       I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#endif

       int I_C2=0;
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
       if(true==initial_assignment)
         {
           I_C2=this->tracer.Number_of_neighbours[ig];
           if(I_C2>=N_NEIGHBOURS_MAX)
              I_C2=N_NEIGHBOURS_MAX-1;
          }
#elif defined _USE_INVARIANT_TIDAL_FIELD_III_
          real_prec C2 = this->cwclass.Invariant_TF_III[id];
          I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#elif defined _USE_DELTA3_
          real_prec C2 = this->cwclass.DELTA3[id];
          I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#endif

          int I_C3=0;
#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
          if(true==initial_assignment)
           {
             real_prec min_sep = this->tracer.min_separation_in_cell[id];
             I_C3 = get_bin(min_sep, 0, N_c, delta_min_sep , this->bin_accumulate_borders);
           }
#elif defined (_USE_INVARIANT_TIDAL_FIELD_IV_)
          real_prec C3 = this->cwclass.Invariant_TF_IV[id];
          I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined _USE_TIDAL_ANISOTROPY_
          real_prec C3 = this->cwclass.Tidal_Anisotropy[id];
          I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined _USE_S2_
          real_prec C3 = this->cwclass.S2[id];             // s²
          I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#endif
          int I_CV1=0;
#ifdef _USE_TRACERS_IN_CELLS_
          if(true==initial_assignment)
              I_CV1=get_bin(static_cast<int>(MOCK_DEN_FIELD[id]),0,N_BINS_TRACERS_IN_CELLS,DELTA_TRACERS_IN_CELLS,true);

#elif defined (_USE_INVARIANT_SHEAR_VFIELD_I_)
          real_prec CV1 = invariant_field_I(this->cwclass.lambda1_vs[id],this->cwclass.lambda2_vs[id],this->cwclass.lambda3_vs[id]); // Not yet assigned to container
          I_CV1= get_bin(CV1, CV1_MIN, N_CV_BIN1,DELTA_CV1, this->bin_accumulate_borders);
#elif defined _USE_NABLA2DELTA_
          real_prec CV1 = this->cwclass.N2D[id];      // Nabla² ð
          I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1, s_deltas.prop7,this->bin_accumulate_borders);
#endif

         int I_CV2=0;
#ifdef _ASSIGN_MASS_POST_
         if(false==initial_assignment)
           I_CV2=get_bin(log10(this->tracer.Halo[ig].vmax), log10(this->params._VMAXmin()),N_x,log10(this->params._VMAXmax()/this->params._VMAXmin())/static_cast<double>(N_x), true);
#elif defined (_USE_INVARIANT_SHEAR_VFIELD_II_)
         real_prec CV2 = invariant_field_II(this->cwclass.lambda1_vs[id],this->cwclass.lambda2_vs[id],this->cwclass.lambda3_vs[id]);// Not yet assigned to container
         I_CV2= get_bin(CV2, CV2_MIN, N_CV_BIN2,DELTA_CV2, this->bin_accumulate_borders);
#elif defined _USE_S2DELTA_
         real_prec CV2 = this->cwclass.S2DELTA[id];         // s²ð
         I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2, s_deltas.prop8,this->bin_accumulate_borders);
#endif

        int I_CV3=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
        real_prec CV3 = invariant_field_III(this->cwclass.lambda1_vs[id],this->cwclass.lambda2_vs[id],this->cwclass.lambda3_vs[id]);// Not yet assigned to container
        I_CV3= get_bin(CV3, CV3_MIN, N_CV_BIN3,DELTA_CV3, this->bin_accumulate_borders);
#elif defined _USE_S3_
        real_prec CV3 = this->cwclass.S3[id];                                   // s³
        I_CV3= get_bin(CV3, this->s_mins.prop9, N_CV_BIN3, s_deltas.prop9,this->bin_accumulate_borders);
#endif




#ifndef _BIN_ACCUMULATE_
        if(xdm >=this->s_mins.prop1 && xdm <=this->s_maxs.prop1)
#endif

#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_DELTA2_)
          if(C1<= this->s_maxs.prop4 && C1>= this->s_mins.prop4)
#endif
#endif

#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_DELTA3_)
            if(C2<= this->s_maxs.prop5 && C2>= this->s_mins.prop5)
#endif 
#endif

#ifndef _BIN_ACCUMULATE_
#if defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_S2_) || (_USE_INVARIANT_TIDAL_FIELD_I_)
              if(C3<= this->s_maxs.prop6 && C3>= this->s_mins.prop6)
#endif
#endif

#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_SHEAR_VFIELD_I_) || defined (_USE_NABLA2DELTA_)
                if(CV1<= this->s_maxs.prop7 && CV1>= this->s_mins.prop7)
#endif
#endif

#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_SHEAR_VFIELD_II_) || defined (_USE_S2DELTA_)
                  if(CV2<= this->s_maxs.prop8 && CV2>= this->s_mins.prop8)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_SHEAR_VFIELD_III_) || defined (_USE_S3_)
                   if(CV3<= this->s_maxs.prop9 && CV3>= this->s_mins.prop9)
#endif
                     {
#endif


                       ULONG index_bins = index_11d(I_X, I_CWT, I_MK,I_CWV,I_VK,I_C1,I_C2, I_C3, I_CV1,I_CV2, I_CV3, this->n_cwt, this->n_sknot_massbin,this->n_cwv, this->n_vknot_massbin, N_a, N_b, N_c,N_v,N_x,N_CV_BIN3);


#if defined (_USE_MULTISCALE_PROPERTY_ASSIGNMENT_) || defined (_SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_)
                       cell_info_tr[id].Theta_bin=index_bins;     // New_Theta-bin in which the cell ID has been identified
#endif



#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_  //this is only done for vmax so far, since we shall  use that property as main property
                       if(this->tracer.Halo[ig].identity >0 )  // Select tracer associated to a DM particle
			                   {
                           count_dm_full++;
#endif

#ifndef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_
#ifdef  _USE_GLOBAL_MASS_FUNCTION_
#ifdef _MASS_ASSIGNMENT_TO_REFERENCE_
			                     int N_props_left=0;
#endif
#else
                  			   ULONG N_tracers_in_theta_bin = this->dm_properties_bins[index_bins].masses_bin_properties.size();
			                     ULONG N_tracers_left =  number_in_theta_ref[index_bins]; // this container is updated below
#endif
			   // If we have available properties in the theta-bin to assign, proceed 
                			     if(N_tracers_left>0) //when we use the ref as a mock, this condition will be always satisfied by construction
			                       {
                			         bool flag=false;
			                         while(false == flag) // este while obliga a elegir una masa -no elegida antes- de entre las disponibles
                      		 		  {
#ifdef _USE_OMP_TEST_
                                  int i_halo_index= gsl_rng_uniform_int(randa,N_tracers_in_theta_bin);
#else
                                  int i_halo_index= gsl_rng_uniform_int(rn,N_tracers_in_theta_bin);
#endif
                                  bool used_tracer = this->dm_properties_bins[index_bins].used_mass[i_halo_index]; //true or false if the mass was already chosen or not

                                  ULONG itest=0;
                                  if(false == used_tracer)// if the property of the tracer (e.g., vmax, mass) has not been used before, then assign that property to the current particle i
                      				     {
#ifdef _USE_VMAX_AS_OBSERVABLE_
                                      this->tracer.Halo[ig].vmax= this->dm_properties_bins[index_bins].masses_bin_properties[i_halo_index]; //prop_to_assign;
#elif defined _USE_MASS_AS_OBSERVABLE_
                                      this->tracer.Halo[ig].mass= this->dm_properties_bins[index_bins].masses_bin_properties[i_halo_index];  //prop_to_assign;
#endif
                      		   		      this->tracer.Halo[ig].observed=true; // Mark this tarcer as already assigned a property. This will be used below.
                                      this->dm_properties_bins[index_bins].used_mass[i_halo_index] = true; //mark this tracer as already assigned a property. This is used withhin this loop.
                     				          flag = true;
				                              counter_fmf++;

#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_  //this is only done for vmax so far, since we shall  use that property as main property
                                      count_dm++;
#endif


                                      if(number_in_theta_ref[index_bins]>1)
#ifdef _USE_OMP_TEST_
#pragma omp atomic update
#endif
                                        number_in_theta_ref[index_bins]--;
                                      else
                                        number_in_theta_ref[index_bins]=0;

                        				    }

                        				 }
			                       }
                                    // andres
//                         cout<<ig<<"   "<<N_tracers_left<<"  "<<i_halo_index<<"    "<<counter_fmf<<"  "<<N_tracers_in_theta_bin<<"  "<<N_tracers_left<<"  "<<this->tracer.Halo[ig].mass<<"  "<<itest<<endl;
//                      if(index_bins==322232)   cout<<ig<<"   "<<id<<"   "<<N_tracers_left<<"  "<<index_bins<<"  "<<counter_fmf<<"  "<<N_tracers_in_theta_bin<<"  "<<this->tracer.Halo[ig].mass<<endl;



#endif  // end of ifndef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_ line 7068

#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_  //this is only done for vmax so far, since we shall  use that property as main property
                          }
#endif


#ifndef _BIN_ACCUMULATE_
                     }  // end of ifs for theta-bins
#endif

#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#if defined(_FULL_VERBOSE_) && !defined (_USE_OMP_TEST_)
                   So.message_screen_flush("Number of dm-like tracers with assigned property = ",static_cast<int>(counter_fmf));
#endif
#else
#ifndef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_
#if defined _FULL_VERBOSE_ && !defined (_USE_OMP_TEST_)
//                  So.message_screen_flush("Number of tracers with assigned property = ",static_cast<int>(counter_fmf));
#endif
#endif
#endif
              }// end loop over tracers or grid
#ifdef _FULL_VERBOSE_
            cout<<endl;
#endif

#ifdef _USE_OMP_TEST_
        }  // end of parallel region
#endif

#if !defined (_SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_) && !defined (_USE_MULTISCALE_PROPERTY_ASSIGNMENT_)
#ifdef _FULL_VERBOSE_
    this->So.DONE();
    So.message_screen("Number of tracers with assigned property = ",static_cast<int>(counter_fmf));
    cout<<endl;
#endif
#endif
    // ***********************


#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
    vector<real_prec>left_overs_properties;
    vector<bool>left_overs_properties_used;



    for(ULONG i = 0; i < this->dm_properties_bins.size(); ++i)
      for(ULONG j = 0; j < this->dm_properties_bins[i].used_mass.size(); ++j)
        if(false==this->dm_properties_bins[i].used_mass[j])
          {
            left_overs_properties.push_back(this->dm_properties_bins[i].masses_bin_properties[j]);
            left_overs_properties_used.push_back(false);
          }
#ifdef _FULL_VERBOSE_
       So.message_screen("Assigning randomly to remainig set -dm:");
       cout<<endl;
#endif
    for(ULONG ig=0; ig < this->tracer.NOBJS; ++ig)
      if( this->tracer.Halo[ig].identity>0  && false==this->tracer.Halo[ig].observed)// proceed if this object is random
        {
          bool flag=false;
          while(false==flag){
            int jj = gsl_rng_uniform_int(rn,left_overs_properties.size());
            if(false==left_overs_properties_used[jj])
              {
#ifdef _USE_MASS_AS_OBSERVABLE_
                this->tracer.Halo[ig].mass=left_overs_properties[jj] ;
#elif defined _USE_VMAX_AS_OBSERVABLE_
                this->tracer.Halo[ig].vmax= left_overs_properties[jj] ;
#endif
                this->tracer.Halo[ig].observed=true;
                left_overs_properties_used[jj]=true;
                count_dm++;
                flag=true;
              }
           }
#ifdef _FULL_VERBOSE_
       So.message_screen_flush("Number of ran-like tracers with assigned property (remaining) = ",static_cast<int>(count_dm));
#endif
         }
    cout<<endl;
    left_overs_properties.clear();
    left_overs_properties.shrink_to_fit();
    left_overs_properties_used.clear();
    left_overs_properties_used.shrink_to_fit();


    // This part is meant for the random distributed tracers



#ifdef _FULL_VERBOSE_
    So.message_screen("Assigning property to random tracers. Expected ", this->tracer.Ntracers_ran);
#endif
    ULONG count_ran=0;

    // This loop will not work perfectly: The container this->dm_properties_for_randoms_bins[theta_bin].masses_bin_properties[]
    // has been filled with the reference tracers below a threshold Vmax, but this doest not guarranty that
    // if a object here is a random living in a theta bin with Nrans inside, we can authomatically find "Nrans"
    // in a theta bin of the this->dm_properties_for_randoms_bins[theta_bin].masses_bin_properties[]
    // There are some tracers marked as "not-observed" and "random" living in a theta_bin where the reference has none of these "randoms".
    // Thes objeects are theferefre not assigned a value of vmax. But the container this->dm_properties_for_randoms_bins has more,
    // so we will just assign vmax to the remainnhg orphan objects randomly from the left-overs.
    for(ULONG ig=0; ig < this->tracer.NOBJS; ++ig)
      {
        if( this->tracer.Halo[ig].identity<0)// proceed if this object is random
          {
            count_ran_full++;
            ULONG id=this->tracer.Halo[ig].GridID;
            ULONG index_bins = cell_info_tr[id].Theta_bin;
            ULONG N_tracers_in_theta_bin_r = this->dm_properties_for_randoms_bins[index_bins].masses_bin_properties.size();
            ULONG N_tracers_left_r =  number_in_theta_ref_for_randoms[index_bins]; // this container is updated below
           // If we have available properties in the theta-bin
        //    cout<<ig<<" "<<this->tracer.Halo[ig].identity<<"   "<<this->tracer.Halo[ig].observed<<"   "<<N_tracers_left_r<<endl;
            if(N_tracers_left_r>0) //when we use the ref as a mock, this condition will be always satisfied by construction
              {
                bool flag=false;
                while(flag == false) // este while obliga a elegir una masa -no elegida antes- de entre las disponibles
                 {
                   int i_halo_index_r= gsl_rng_uniform_int(rn,N_tracers_in_theta_bin_r);
                   bool used_tracer = this->dm_properties_for_randoms_bins[index_bins].used_mass[i_halo_index_r]; //true or false if the mass was already chosen or not
                   if(false == used_tracer)// if the property of the tracer (e.g., vmax, mass) has not been used before, then assign that property to the current particle i
                    {
                      this->tracer.Halo[ig].vmax= this->dm_properties_for_randoms_bins[index_bins].masses_bin_properties[i_halo_index_r]; //prop_to_assign;
                      this->tracer.Halo[ig].observed = true; // Mark this tarcer as already assigned a property. This will be used below.
                      this->dm_properties_for_randoms_bins[index_bins].used_mass[i_halo_index_r] = true; //mark this tracer as already assigned a property. This is used within this loop.
                      flag = true;
                      counter_fmf++;
                      count_ran++;
                      number_in_theta_ref_for_randoms[index_bins]--;
                    }
                 }
           }
        }
#ifdef _FULL_VERBOSE_
       So.message_screen_flush("Number of ran-like tracers with assigned property = ",static_cast<int>(count_ran));
#endif
      }
    cout<<endl;
        // Left overs:

    for(ULONG i = 0; i < this->dm_properties_for_randoms_bins.size(); ++i)
      for(ULONG j = 0; j < this->dm_properties_for_randoms_bins[i].used_mass.size(); ++j)
        if(false==this->dm_properties_for_randoms_bins[i].used_mass[j])
          {
            left_overs_properties.push_back(this->dm_properties_for_randoms_bins[i].masses_bin_properties[j]);
            left_overs_properties_used.push_back(false);
          }

    this->dm_properties_for_randoms_bins.clear();
    this->dm_properties_for_randoms_bins.shrink_to_fit();
#ifdef _FULL_VERBOSE_
       So.message_screen("Assigning randomly remainig set (random):");
       cout<<endl;
#endif
    for(ULONG ig=0; ig < this->tracer.NOBJS; ++ig)
      if( this->tracer.Halo[ig].identity<0  && false==this->tracer.Halo[ig].observed)// proceed if this object is random
        {
          bool flag=false;
          while(false==flag){
            int jj = gsl_rng_uniform_int(rn,left_overs_properties.size());
            if(false==left_overs_properties_used[jj])
              {
#ifdef _USE_MASS_AS_OBSERVABLE_
                this->tracer.Halo[ig].mass=left_overs_properties[jj] ;
#elif defined _USE_VMAX_AS_OBSERVABLE_
                this->tracer.Halo[ig].vmax= left_overs_properties[jj] ;
#endif
                this->tracer.Halo[ig].observed=true;
                left_overs_properties_used[jj]=true;
                count_ran++;
                flag=true;
              }
           }
#ifdef _FULL_VERBOSE_
       So.message_screen_flush("Number of ran-like tracers with assigned property (remaining) = ",static_cast<int>(count_ran));
#endif
         }
    left_overs_properties.clear();
    left_overs_properties.shrink_to_fit();
    left_overs_properties_used.clear();
    left_overs_properties_used.shrink_to_fit();

    cout<<endl;
    So.message_screen("Number of tracers_dm with assigned property = ",static_cast<int>(count_dm));
    So.message_screen("Number of tracers_dm with assigned property_full = ",static_cast<int>(count_dm_full));//to check, ok
    So.message_screen("Number of tracers_ran with assigned property = ",static_cast<int>(count_ran));
    So.message_screen("Number of tracers_ran with assigned property_full = ",static_cast<int>(count_ran_full));//to check, ok



#endif


#ifndef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_
      number_in_theta_ref.clear();
      number_in_theta_ref.shrink_to_fit();
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
      number_in_theta_ref_for_randoms.clear();
      number_in_theta_ref_for_randoms.shrink_to_fit();
#endif


#endif

#ifdef _USE_TOTAL_MASS_IN_CELL_
      MOCK_MASS_FIELD.clear();
      MOCK_MASS_FIELD.shrink_to_fit();
#endif



#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_
      // *********************************** EXPLANATION NOTE*************************************************************************
      // In the following loop we identify cells with obtjects above the X threshold Mth, imposing already positions for the masses
      // But we can avoid that threshold here, and then we would be counting all cells.
      // Note that counter_high_mass and N_props_multiscale are merely consequence of the X (vmax or mass) function
      // These numbers are computed from the reference.
      // With the number of tracers in each interval (level), we then can pass to assign the largest X values onto the low resolution mesh.
      // The fact that we assig the larger X-values lies in the fact that we use the X values *ordered* within the bins of theta and we incrase the index
      // every time we assign a mass, such that at the next assignment, the code will use the next less "massive" object in the corresponding theta-bin.
      // Note that this applies for vmax and Mass, as these display a positive correlation. If a halo property displays an anticorrelation (but still strong)
      // with the underlying dm field and hence with the mass, the assignment would have to be done from bottom-to-top.
      // Regimes of multi-scale mass assignment:
   // i) M>M_multiscale_1, assign with NFFT_MS_1
   // ii) M_multiscale_2<M< M_multiscale_1 assign with NFFT_MS_2
   // iii) M_multiscale_3<M< M_multiscale_2 assign with NFFT_MS_3
   // iv) M< M_multiscale_3 done at the particle level, assigning the remaining masses to tracers according to the theta bin of their respective cell
   // Steps:
   // i) Count number of masses available in the different regimes O
   N_props_0=0;
   N_props_1=0;
   N_props_2=0;
   N_props_3=0;
   N_props_4=0;

#ifdef _FULL_VERBOSE_
   So.message_screen("Number of available masses in the different multi-scale levels obtained from the reference: ");
#endif
// *******************************************************************************************************************************
#ifdef _USE_MULTISCALE_LEVEL_4_
   N_props_4 = this->tracer_ref.N_props_4; // number of masses in level 4   M > M_th4
#ifndef _ASSIGN_TO_CALIBRATION_
   N_props_4 = static_cast<ULONG>(floor(TOLERANCE_FACTOR_L4*N_props_4));
#endif
#ifdef _FULL_VERBOSE_
   So.message_screen("Number of masses in level 4 =", N_props_4);
   So.message_screen("Suggested grid size for this population = ", floor(pow(N_props_4,1./3.)));
   So.message_screen("Using", NFT_LOW_4);
#endif
   real_prec d1_l4=this->Lbox/static_cast<real_prec>(NFT_LOW_4);		/* grid spacing x-direction */
   vector<int>number_in_cells_aux_l4(NGRID_MS_4,0); // this would be for level 4
   if(NFT_LOW_4>this->Nft)
     So.message_error("NFT for level 4 smaller than original value");
#endif
// *******************************************************************************************************************************

#ifdef _USE_MULTISCALE_LEVEL_3_
   N_props_3 = this->tracer_ref.N_props_3; // number of masses in level 3   Mth3 < M < M_th4
#ifndef _ASSIGN_TO_CALIBRATION_
   N_props_3 = static_cast<ULONG>(floor(TOLERANCE_FACTOR_L3*N_props_3));
#endif
#ifdef _FULL_VERBOSE_
   So.message_screen("Number of masses in level 3 =", N_props_3);
   So.message_screen("Suggested grid size for this population = ", floor(pow(N_props_3,1./3.)));
   So.message_screen("Using", NFT_LOW_3);
#endif
   real_prec d1_l3=this->Lbox/static_cast<real_prec>(NFT_LOW_3);		/* grid spacing x-direction */
   vector<int>number_in_cells_aux_l3(NGRID_MS_3,0);
   if(NFT_LOW_3>this->Nft)
     So.message_error("NFT for level 3 smaller than original value");
#endif
   // *******************************************************************************************************************************

#ifdef _USE_MULTISCALE_LEVEL_2_
   N_props_2 = this->tracer_ref.N_props_2; // number of masses in level 2   Mth2 < M < M_th3
#ifndef _ASSIGN_TO_CALIBRATION_
   N_props_2 = static_cast<ULONG>(floor(TOLERANCE_FACTOR_L2*N_props_2));
#endif
#ifdef _FULL_VERBOSE_
   So.message_screen("Number of masses in level 2 =", N_props_2);
   So.message_screen("Suggested grid size for this population = ", floor(pow(N_props_2,1./3.)));
   So.message_screen("Using", NFT_LOW_2);
#endif
   real_prec d1_l2=this->Lbox/static_cast<real_prec>(NFT_LOW_2);		/* grid spacing x-direction */
   vector<int>number_in_cells_aux_l2(NGRID_MS_2,0); // this would be for level 2
   if(NFT_LOW_2>this->Nft)
     So.message_error("NFT for level 2 smaller than original value");
#endif
   // *******************************************************************************************************************************
#ifdef _USE_MULTISCALE_LEVEL_1_
   N_props_1 = this->tracer_ref.N_props_1; // number of masses in level 1   Mth1 < M < Mth2
#ifndef _ASSIGN_TO_CALIBRATION_
   N_props_1 = static_cast<ULONG>(floor(TOLERANCE_FACTOR_L1*N_props_1));
#endif
#ifdef _FULL_VERBOSE_
   So.message_screen("Number of masses in level 1 =", N_props_1);
   So.message_screen("Suggested grid size for this population = ", floor(pow(N_props_1,1./3.)));
   So.message_screen("Using the orignal,", this->Nft);
   cout<<endl;
#endif
#endif
   // *******************************************************************************************************************************
// PARTICLE LEVEL
   vector<int>number_in_cells_aux(this->NGRID,0); // this is needed regardless the level defined
   ULONG ep_ncounts=this->tracer_ref.NOBJS-this->tracer.NOBJS;
   if(this->tracer_ref.NOBJS-this->tracer.NOBJS<0)
       ep_ncounts*=-1.0;

   // Using ep_ncounts, N0 has then the info of the remainig set of properties to assign in the mock, not in the reference.
   N_props_0=this->tracer_ref.NOBJS-ep_ncounts-(N_props_1+N_props_2+N_props_3+N_props_4);

#ifdef _FULL_VERBOSE_
   So.message_screen("Total to be assigned from grid-approach = ", N_props_1+N_props_2+N_props_3+N_props_4);
   cout<<endl;
#endif

   if(this->tracer_ref.NOBJS<N_props_1+N_props_2+N_props_3+N_props_4)
    {
       So.message_warning("Error in line", __LINE__);
       throw std::invalid_argument( "ULONG defined variable was assigned a negative value" );
       So.message_warning("This might be caused by pre-proc declarations defining 'mass_cuts' or 'mass_bins'");
       So.message_screen("N_props=", N_props_0);
       So.message_screen("Assigned =",N_props_1+N_props_2+N_props_3+N_props_4);
       So.message_warning("BAM stops here to avoid propagation of wrong numbers");
       exit(0);
    }
#ifdef _FULL_VERBOSE_
   So.message_screen("Number of reference masses to be assigned at particle level (level 0) = ",N_props_0);
   cout<<endl;
   So.message_screen("Identifying galaxy id to get X-values in the different levels");
#endif

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG ig=0;ig < this->tracer.NOBJS; ++ig)
       this->tracer.Halo[ig].observed=false;  // this will be used below, so let keep it inzialized like that

   // This loop is valid as long as X-values have been assigned before, even if they will be reassigned.
   // Since this has a push_bvack inside, I leave it without parallelization
   for(ULONG ig=0;ig < this->tracer.NOBJS; ++ig)
     {
       ULONG id=this->tracer.Halo[ig].GridID;
       cell_info_tr[id].gal_index.push_back(ig);  // allocate the galaxy ID "ig" in the corresponding grid cell "id"
    }

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG ig=0;ig < this->tracer.NOBJS; ++ig)
     {
       ULONG id=this->tracer.Halo[ig].GridID;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
       number_in_cells_aux[id]++;    // this is the number counts in the ref grid and will be updated in every level

#if defined _USE_MULTISCALE_LEVEL_4_ || defined _USE_MULTISCALE_LEVEL_3_ || defined _USE_MULTISCALE_LEVEL_2_
       ULONG il,jl,kl,id_low;
       real_prec x = this->tracer.Halo[ig].coord1;
       real_prec y = this->tracer.Halo[ig].coord2;
       real_prec z = this->tracer.Halo[ig].coord3;
#endif

       // Tracers in level 4 (grid-based)
#ifdef _USE_MULTISCALE_LEVEL_4_
       il=static_cast<ULONG>(floor(x/d1_l4));
       jl=static_cast<ULONG>(floor(y/d1_l4));
       kl=static_cast<ULONG>(floor(z/d1_l4));
       if(il==NFT_LOW_4)
      	 il--;
       if(jl==NFT_LOW_4)
	       jl--;
       if(kl==NFT_LOW_4)
	       kl--;
       id_low = index_3d(il,jl,kl,NFT_LOW_4,NFT_LOW_4);
#ifdef _USE_OMP_
#pragma omp atomic
#endif
       number_in_cells_aux_l4[id_low]++;     //Counting the number of tracers in the level-4 mesh
       this->tracer.Halo[ig].GridID_l4=id_low;  //Allocate the level4-mesh ID in which each tracer resides
#endif


       // Tracers in level 3 (grid-based)
#ifdef _USE_MULTISCALE_LEVEL_3_  // ensures level 1 is always used
       il=static_cast<ULONG>(floor(x/d1_l3));
       jl=static_cast<ULONG>(floor(y/d1_l3));
       kl=static_cast<ULONG>(floor(z/d1_l3));
       if(il==NFT_LOW_3)
      	 il--;
       if(jl==NFT_LOW_3)
	       jl--;
       if(kl==NFT_LOW_3)
      	 kl--;
       id_low = index_3d(il,jl,kl,NFT_LOW_3,NFT_LOW_3);
#ifdef _USE_OMP_
#pragma omp atomic
#endif
       number_in_cells_aux_l3[id_low]++;          //Counting the number of tracers in the level-3 mesh
       this->tracer.Halo[ig].GridID_l3=id_low;  //Allocate the level3-mesh ID in which each tracer resides
#endif

       // Number of objects in level 2 (grid-based)
#ifdef _USE_MULTISCALE_LEVEL_2_  // ensures level 1 is always used
       il=static_cast<ULONG>(floor(x/d1_l2));
       jl=static_cast<ULONG>(floor(y/d1_l2));
       kl=static_cast<ULONG>(floor(z/d1_l2));
       if(il==NFT_LOW_2)
      	 il--;
       if(jl==NFT_LOW_2)
      	 jl--;
       if(kl==NFT_LOW_2)
      	 kl--;
       id_low = index_3d(il,jl,kl,NFT_LOW_2,NFT_LOW_2);
#ifdef _USE_OMP_
#pragma omp atomic
#endif
       number_in_cells_aux_l2[id_low]++;
       this->tracer.Halo[ig].GridID_l2=id_low;
#endif
       // Number of objects in level 1: (grid- based)
     }
   So.DONE();

   // This follows the HADRON-approach, in which
   // we order masses in each theta*-bin from the max to the minimum and assign
   //them in that order to tracers located in randomly selected cells (within the same theta-bin)

#ifdef _FULL_VERBOSE_
#ifdef _USE_MASS_AS_OBSERVABLE_
   So.message_screen("Sorting masses in bins of theta");
#elif defined _USE_VMAX_AS_OBSERVABLE_
   So.message_screen("Sorting Vmax in bins of theta");
#endif
   So.message_screen("used in the multiscale approach: ");
#endif

#ifdef _USE_OMP_
//#pragma omp parallel for
#endif
   for(ULONG it=0; it < LENGHTdm; ++it)
     {
       int ll=this->dm_properties_bins[it].masses_bin_properties.size();
       if(ll>0)
	      {
           gsl_vector *aux_properties;
           aux_properties=gsl_vector_alloc(ll);
	   for(int im=0; im< ll; ++im )
             gsl_vector_set(aux_properties,im,this->dm_properties_bins[it].masses_bin_properties[im]);
           gsl_sort_vector(aux_properties);
	   for(int im=0; im< ll; ++im )
             this->dm_properties_bins[it].masses_bin_properties[im]=gsl_vector_get(aux_properties,ll-im-1);
           gsl_vector_free(aux_properties);
	 }
     }

   So.DONE();

   // This is independent of the level; will be used for all levels in order to pinpoint the position of the ordered X-properties

   // ****************************************************************************************
   // ****************************************************************************************
   // ****************************************************************************************
   // ****************************************************************************************
   // ****************************************************************************************
   // ****************************MULTISCALE Assignment **************************************

   vector<ULONG>count_aux_cell(LENGHTdm,0);  // This must be outside of the do-while loop
   vector<ULONG>cells_id_still_to_assign;
   ULONG counter_multiscale=0;
   bool do_flag=false;
   int loop_counter=0;
   // ****************************************************************************************
   // ****************************************************************************************
   // ***************************************************Level 4 *****************************

#ifdef _USE_MULTISCALE_LEVEL_4_

   //Container to allocate the ID of the cells with particles still to get assigned properties
   //This has to be cleared, resized and updated for every level
   vector<s_cell_info>cell_info_cell(NGRID_MS_4);

   vector<int>ncount_low_cell(NGRID_MS_4);


   // This function returns the list of grid id in each ID_low (all cells included, regarless with or without tracers)
   get_low_res_id_from_high_res_id(this->Nft, NFT_LOW_4, cell_info_cell);

#ifdef _FULL_VERBOSE_
   So.message_screen("Identifying level-4 empty cells");
#else
   So.message_screen("Level-4");
#endif

   // THis is meant to be shuffled in order to run randomly the cells
   for(ULONG id_low=0; id_low< NGRID_MS_4; ++id_low)
     if(number_in_cells_aux_l4[id_low]>0)// Pinpoint the low_res cells with tracers. These tracers do not have an assigned mass
       cells_id_still_to_assign.push_back(id_low); //
   So.DONE();

   number_in_cells_aux_l4.clear();
   number_in_cells_aux_l4.shrink_to_fit();


#ifdef _USE_EXCLUSION_CELL_LEVEL4_
   So.message_screen("Identifying neighbouring cells in level 4 mesh");
   vector<s_nearest_cells> ID_nearest_cells_lowres(NGRID_MS_4);
   get_neighbour_cells(NFT_LOW_4,NUMBER_OF_NEIGHBOUR_CELL_FOR_EXCLUSION_LEVEL4,ID_nearest_cells_lowres);
   So.DONE();
   vector<bool>rejected_cell_lowres;

#endif


#ifdef _USE_EXCLUSION_CELL_LEVEL4_HIGHRES_
   vector<bool>visited_highres(this->NGRID,false);
#endif

   ULONG N_cells=cells_id_still_to_assign.size();

#ifdef _FULL_VERBOSE_
   So.message_screen("Shuffling order of level-4 cells containing tracers without assigned masses");
#endif
   gsl_ran_shuffle(rn,&cells_id_still_to_assign[0], cells_id_still_to_assign.size() ,sizeof(ULONG));

#ifdef _FULL_VERBOSE_
   So.DONE();
   So.message_screen("*****Going through LEVEL 4 mesh to assign ", N_props_4, "properties");
#endif

   counter_multiscale=0;

   loop_counter=0;  // this couts the number of for loops wr ned to perform in order to get the expected number of tracers at this level
   do_flag=false;

   do{
     // Assignment campaing, going through level4-cells and assigning a value of X-property *to one tracer per cell*.
     // For that we select low_res cells with tracers without assigned X-value,
     // and do a loop over these cells. For each cell (previously randomly sorted), we randomly select one of the 
     // original resolution cells, and randomly pick up one tracer therefrom.
     // We assign Xvalues top-to-bottom (vmax-mass correlate positively. For negative correlation, bottom-to-top must be applied)
     // until we reach the number N_props_multiscales, which corresponds to the X-threshold_multilscale
     // The do loop guarranties that the number of assigned X-values equals the number of tracers in this multl-scale level.
     // When this condition is achieved, the do loop ends.


      ULONG loop_assinged=0;
#ifdef _USE_EXCLUSION_CELL_LEVEL4_
      rejected_cell_lowres.resize(NGRID_MS_4,false);
#endif

#ifdef _USE_OMP_BARR_L4_
#pragma omp parallel for  // watch out, there is a random generator iside these loops!
#endif
     for(ULONG id_ini=0; id_ini< N_cells; ++id_ini) // loop over low-res cells containing tracers without assigned X-vales
       {
#ifdef _USE_OMP_BARR_L4_
         if(do_flag) continue;
#endif

    	   ULONG id_low = cells_id_still_to_assign[id_ini];  // Low Grid ID. These come randomized
         // THis must be raandomized only once, such that we guarantte that we pass through every vell once (in a random way)

  	     ULONG N_cells_in_cells= cell_info_cell[id_low].gal_index.size();  // Number of high-res cells in a low-res cell
         ULONG index_cell_in_low_res_cell= gsl_rng_uniform_int(rn,N_cells_in_cells); //select randomy one of the high-res cells contained in the low-res cell
         ULONG id=cell_info_cell[id_low].gal_index[index_cell_in_low_res_cell]; // Retrieve the id. This chosen id will be used below to assign mass

	      if(number_in_cells_aux[id]>0)// if have available tracers in the high-res cell:
	       {

	         ULONG index_bins = cell_info_tr[id].Theta_bin; // Get the theta-bin in which the cell with "id" has been classified. This has been calculated in the step M<M*
#ifdef _USE_EXCLUSION_CELL_LEVEL4_HIGHRES_
           if(false==visited_highres[id]) // if the high-res cell hsa not been visited, proceed
            {
#endif

#ifdef _USE_EXCLUSION_CELL_LEVEL4_
       	     bool used_previous_cell=false;  // ask whether the previos cell is one of the neighbours of the current one
             bool rejected_cell=false;  // ask whether the previos cell is one of the neighbours of the current one
	           ULONG N_neighbour_cells=ID_nearest_cells_lowres[id_low].close_cell.size(); // Number of neighbours computed for the low_res
             if(id_ini > 0) // do not ask the fiorst visited cell
	  	         {
               for(int in = 0 ; in < N_neighbour_cells; ++in)// Loop over the neighbour cells of the current id_ini cell
                 {
		               if(ID_nearest_cells_lowres[id_low].close_cell[in]==cells_id_still_to_assign[id_ini-1]) //If the current cell has as neighbour the previous visited cell, break
		               {
		                 used_previous_cell=true; // if any of the neighbour cells was visited before, 
                     rejected_cell_lowres[id_ini]=true;
		                 break;                   // break the loop over the neighbout cells and report it.
     		           }
                 }
               }
       // Proceed to assign if the previous cell is not one of the closest neighbours of the current cell
	     // or if the previus  cell id-1 was rejected because being neighbour of the id-2
   
          // if the current cell is not neighbour of the previous visited, or if the previous cell was rejected, I can use this one even if the current cell is neighbour of the previous one
            if( false==used_previous_cell ||  (true==rejected_cell_lowres[id_ini-1] && true==used_previous_cell) ) 
             { 
        
#endif  // end of use_exclusion in cells

               ULONG N_tracers_without_prop_in_cell=cell_info_tr[id].gal_index.size();
               ULONG jk= gsl_rng_uniform_int(rn,N_tracers_without_prop_in_cell); //choose a random integer in  [0,N_props_in_cell_without_mass). N_props_withopur mass must be the one fixed
           		 ULONG ig=cell_info_tr[id].gal_index[jk];                    // get the Galaxy ID for the cell id in the position jk of tracers without mass
               ULONG i_prop_halo_label=count_aux_cell[index_bins];  // get the (top-bottom sorted) id of the property in this theta-bin

//               cout<<id<<"  "<<index_bins<<"   "<<ig<<endl;
#ifndef _ASSIGN_TO_CALIBRATION_
               if(this->dm_properties_bins[index_bins].masses_bin_properties.size()>0) // This will be always the case when the new density filed is the same usd for the calibration
                  {
#endif

                    real_prec assigned_property=this->dm_properties_bins[index_bins].masses_bin_properties[i_prop_halo_label]; // Read the X-value of the ig tracer in the theta-bin index_bins

                   // if this tracer has not been used and its property is above the threshold for this level,AND if we have trracers available in this theta bin
                    if(false==this->tracer.Halo[ig].observed && assigned_property>=this->Prop_threshold_multi_scale_4 )
                     {


#ifdef _USE_MASS_AS_OBSERVABLE_
                       this->tracer.Halo[ig].mass=assigned_property ;
#elif defined _USE_VMAX_AS_OBSERVABLE_
                       this->tracer.Halo[ig].vmax=assigned_property ;
#endif
              		     this->tracer.Halo[ig].observed=true;
                       this->dm_properties_bins[index_bins].used_mass[i_prop_halo_label]=true;

#ifdef _USE_OMP_BARR_L4_
#pragma omp atomic update
#endif
                      loop_assinged++;  // keep track of assigned in this do loop.

#ifdef _USE_OMP_BARR_L4_
#pragma omp atomic update
#endif
        		     counter_multiscale++; // keepe track of all assigned
#ifdef _USE_OMP_BARR_L4_
#pragma omp atomic update
#endif
                     number_in_cells_aux[id]--; // Remove one tracer from the id cell
#ifdef _USE_OMP_BARR_L4_
#pragma omp atomic update
#endif
                     count_aux_cell[index_bins]++;  // This keeps track of the X-properties used. This is IMPORTANT for it increases 1 as long as a prop is assigned, such that the following assignment will poin to the next less massive objectand will be used below for campaign at proceeding levels
                     ncount_low_cell[id_ini]++;   // this counts the number of assigned tracers in each low_res cell
#ifdef _USE_EXCLUSION_CELL_LEVEL4_HIGHRES_
                     visited_highres[id]=true;  // Keep track of the high res cell not to step in the same twice or thrice
#endif

		   }
#ifndef _ASSIGN_TO_CALIBRATION_
                     }
#endif


#ifdef _USE_EXCLUSION_CELL_LEVEL4_
                     	       }//closes if(false==used_previous_cell)
	     else
	       rejected_cell_lowres[id_low]=true;
#endif

#ifdef _USE_EXCLUSION_CELL_LEVEL4_HIGHRES_
         } // closes if(visited high res is false)
#endif
	     } //end if number_in cells >0


         if(counter_multiscale == N_props_4)
#ifdef _USE_OMP_BARR_L4_
               do_flag=true;
#else
             break;
#endif
      }// end loop over cells

#ifndef _USE_OMP_BARR_
#ifdef _FULL_VERBOSE_
       So.message_screen_flush("Number of masses assigned   = ",static_cast<int>(counter_multiscale));
#endif
#endif



#ifdef _USE_EXCLUSION_CELL_LEVEL4_
//    if(counter_multiscale>= N_props_4-N_MASSES_ABOVE_THRESHOLD_LEVEL4)  //if ig gets exhausted, randomize id of cells
 //     {
//	         So.message_screen("Reshufling cells");
//	         gsl_ran_shuffle(rn,&cells_id_still_to_assign[0],cells_id_still_to_assign.size(),sizeof(ULONG));
	     //    rejected_cell_lowres.resize(this->NGRID,false);
	   //      reshuffle_once=false;
	 //     }
#endif

     loop_counter++;

  //   cout<<"Loop "<<loop_counter<<"   Assigned "<<loop_assinged<<endl;


#ifdef _USE_OMP_BARR_L4_
   }while(flase=do_flag);
   //   N_props_4=counter_multiscale;// update in case the do loop is broken by tolerant factor
#else
       }while(counter_multiscale< N_props_4-N_MASSES_ABOVE_THRESHOLD_LEVEL4);

#endif


   ULONG check_n=0;
   ULONG check_one=0;
   ULONG check_two=0;
   ULONG check_three=0;
   ULONG check_cero=0;
   for(int i=0;i<ncount_low_cell.size();++i)
    {
     check_n+=ncount_low_cell[i];
    if(ncount_low_cell[i]==1)
      check_one++;
    if(ncount_low_cell[i]==2)
      check_two++;
    if(ncount_low_cell[i]==3)
      check_three++;
    if(ncount_low_cell[i]==0)
      check_cero++;
     }
  cout<<check_n<<endl;
  cout<<check_cero<< " cells with 0 tracers"<<endl;

  cout<<check_one<< " cells with 1 tracers"<<endl;
  cout<<check_two<< " cells with 2 tracers"<<endl;
  cout<<check_three<< " cells with 3 tracers"<<endl;


   ncount_low_cell.clear();
   ncount_low_cell.shrink_to_fit();


#ifdef _FULL_VERBOSE_
   cout<<endl;
#ifdef _USE_OMP_BARR_
   So.message_screen("Number of masses assigned   = ",static_cast<int>(counter_multiscale));
#endif
   So.message_screen("Used ", loop_counter, "loops over L4 MESH to assign property");
   cout<<endl;
#endif

   cumulative_counter+=counter_multiscale;
#ifdef _FULL_VERBOSE_
   So.message_screen("Partial number of properties assigned = ",cumulative_counter);
#endif


   So.DONE();
#ifdef _FULL_VERBOSE_
   So.message_screen("Freeing memmory in line", __LINE__);
#endif
   //     used_cell.clear();used_cell.shrink_to_fit();
   cell_info_cell.clear();
   cell_info_cell.shrink_to_fit();
   cells_id_still_to_assign.clear();
   cells_id_still_to_assign.shrink_to_fit();


#ifdef _USE_EXCLUSION_CELL_LEVEL4_
   rejected_cell_lowres.clear();
   rejected_cell_lowres.shrink_to_fit();
   ID_nearest_cells_lowres.clear();
   ID_nearest_cells_lowres.shrink_to_fit();
#endif

   So.DONE();


   // test:
   vector<ULONG>ncounts_aux(NGRID_MS_4,0);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0; i<this->tracer.NOBJS;++i)
     if(true==this->tracer.Halo[i].observed)
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
       ncounts_aux[this->tracer.Halo[i].GridID_l4]++;

#ifdef _FULL_VERBOSE_
   So.message_screen("Maximum number of objects in cells from level-4 mesh = ", get_max<ULONG>(ncounts_aux));
   So.message_screen("Check: number of tracers = ", get_nobjects(ncounts_aux));
#endif

#endif



   // ******************************************************************************************************************************************
   // ******************************************************************************************************************************************
   // ******************************Level 3****************************************************************************************************
#ifdef _USE_MULTISCALE_LEVEL_3_

#ifdef _USE_MULTISCALE_LEVEL_4_
   cell_info_cell.resize(NGRID_MS_3);
#else
   cells_id_still_to_assign.clear();
   cells_id_still_to_assign.shrink_to_fit();
   vector<s_cell_info>cell_info_cell(NGRID_MS_3);
#endif

   // This function returns the list of id in each ID_low (all cells included, regarless with or without mass)
   get_low_res_id_from_high_res_id(this->Nft, NFT_LOW_3,cell_info_cell);

#ifdef _FULL_VERBOSE_
   cout<<endl;
   So.message_screen("Identifying level 3 empty cells");
#else
   So.message_screen("Level-3");
#endif
   for(ULONG id_low=0; id_low< NGRID_MS_3; ++id_low)
     if(number_in_cells_aux_l3[id_low]>0)// Pinpoint the low_res cells with tracers. These tracers do not have an assigned mass
       cells_id_still_to_assign.push_back(id_low); //

   number_in_cells_aux_l3.clear();
   number_in_cells_aux_l3.shrink_to_fit();


#ifdef _USE_EXCLUSION_CELL_LEVEL3_
   So.message_screen("Identifying neighbouring cells in level 3 mesh");

#ifdef _USE_EXCLUSION_CELL_LEVEL4_
   ID_nearest_cells_lowres.resize(NGRID_MS_3);
#else
   vector<s_nearest_cells> ID_nearest_cells_lowres(NGRID_MS_3);
#endif


   get_neighbour_cells(NFT_LOW_3,NUMBER_OF_NEIGHBOUR_CELL_FOR_EXCLUSION_LEVEL3,ID_nearest_cells_lowres);
   So.DONE();
#endif
#ifdef _FULL_VERBOSE_
   So.message_screen("Shuffling order of level 3 cells containing tracers without assigned masses");
#endif

   gsl_ran_shuffle(rn,&cells_id_still_to_assign[0], cells_id_still_to_assign.size() ,sizeof(ULONG));
   So.DONE();

#ifdef _USE_EXCLUSION_CELL_LEVEL3_
#ifdef _USE_EXCLUSION_CELL_LEVEL4_
   rejected_cell_lowres.resize(NGRID_MS_3,false);
#else
   vector<bool>rejected_cell_lowres(NGRID_MS_3,false);
#endif
#endif

#ifdef _FULL_VERBOSE_
   So.message_screen("*****Going through LEVEL 3 mesh to assign ", N_props_3, "properties");
#endif

   counter_multiscale=0;
   do_flag=false;
   do{

      ULONG N_cells=cells_id_still_to_assign.size();

#ifdef _USE_OMP_BARR_L3_
#pragma omp parallel for
#endif
     for(ULONG id_ini=0; id_ini< N_cells; ++id_ini) // loop over low-res cells containing tracers without assigned masses
       {
#ifdef _USE_OMP_BARR_L3_
         if(do_flag) continue;
#endif

      	 ULONG id_low = cells_id_still_to_assign[id_ini];  // Low Grid ID. These come randomized
      	 // ojo que acá estan las high-res celdas que incluso no tienen tracers
#ifdef _USE_EXCLUSION_CELL_LEVEL3_
      	 ULONG N_neighbour_cells=ID_nearest_cells_lowres[id_low].close_cell.size(); // Number of neighbpurs computed for the low_res
#endif
	 int N_cells_in_cells= cell_info_cell[id_low].gal_index.size();  // Number of high-res cells in a low-res cell
	 int index_cell_in_low_res_cell= gsl_rng_uniform_int(rn,N_cells_in_cells); //select randomy one of the high-id's in the low res cell
	 ULONG id=cell_info_cell[id_low].gal_index[index_cell_in_low_res_cell]; // This chosen id will be used below to assign mass
	 ULONG index_bins = cell_info_tr[id].Theta_bin; // Get the theta-bin in which this cell has been classified. This has been calculated in the step M<M*
	 if(number_in_cells_aux[id]>0)
	   {
#ifdef _USE_EXCLUSION_CELL_LEVEL3_
	     bool used_previous_cell=false;  // ask whether the previos cell is one of the neighbours of the current one
	     bool rejected_previous_cell=false;   // to ask whether the previously id-1 cell was rejected because of being neighbour of the cell id-2
	     if(id_ini > 0)
	       {
		 for(int in = 0 ; in < N_neighbour_cells; ++in)
		   if(ID_nearest_cells_lowres[id_low].close_cell[in] == cells_id_still_to_assign[id_ini-1]) /// this is satisfied only once (at most) inside the loop!
		     {
		       used_previous_cell=true;
		       break;
		     }
	       }
	     // Proceed to assign if the previous cell is not one of the closest neighbours of the current cell
	     // or if the previus  cell id-1 was rejected because being neighbour of the id-2
	     if(false==used_previous_cell || (true==rejected_previous_cell && true==used_previous_cell)  ) // if the previous cell was rejected, I can use this one wven if the current cell is neighbour of the previous one
	       {
#endif  // end of use_exclusion in cells
                 ULONG N_props_without_mass_in_cell=cell_info_tr[id].gal_index.size(); //if masses are to be assigned with this, this line can be MOCK_DEN_FIELD[id];
                 ULONG jk= gsl_rng_uniform_int(rn,N_props_without_mass_in_cell); //choose a random integer in  [0,N_props_in_cell_without_mass). N_props_withopur mass must be the one fixed
		 ULONG ig=cell_info_tr[id].gal_index[jk];                       // get the Galaxy ID for the cell id in the position jk of tracers without mass
		 if(false==this->tracer.Halo[ig].observed)
		   {
                     ULONG i_mass_halo_label=count_aux_cell[index_bins];
#ifndef _ASSIGN_TO_CALIBRATION_
                     if(this->dm_properties_bins[index_bins].masses_bin_properties.size()>0) // This will be always the case when the new density filed is the same usd for the calibration
                       {
#endif
                       real_prec assigned_property=this->dm_properties_bins[index_bins].masses_bin_properties[i_mass_halo_label];
#ifdef _USE_MULTISCALE_LEVEL_4_
                       if(assigned_property >= this->Prop_threshold_multi_scale_3 && assigned_property < this->Prop_threshold_multi_scale_4)
#else
                       if(assigned_property >= this->Prop_threshold_multi_scale_3)
#endif
                         {

#ifdef _USE_MASS_AS_OBSERVABLE_
                           this->tracer.Halo[ig].mass=assigned_property ;
#elif defined _USE_VMAX_AS_OBSERVABLE_
                           this->tracer.Halo[ig].vmax=assigned_property ;
#endif
                           this->tracer.Halo[ig].observed=true;
                           this->dm_properties_bins[index_bins].used_mass[i_mass_halo_label]=true;
#ifdef _USE_OMP_BARR_L3_
#pragma omp atomic update
#endif
                           counter_multiscale++;
#ifdef _USE_OMP_BARR_L3_
#pragma omp atomic update
#endif
                           number_in_cells_aux[id]--; // Remove one tracer from the id cell
#ifdef _USE_OMP_BARR_L3_
#pragma omp atomic update
#endif
                           count_aux_cell[index_bins]++;  // This keeps track of the masses used, and will be used below for next-scale campaign
                         }
#ifndef _ASSIGN_TO_CALIBRATION_
                     } // closes  if(this->dm_properties_bins[index_bins].masses_bin_properties.size()>0)
#endif
                     }

#ifdef _USE_EXCLUSION_CELL_LEVEL3_
	       }
	     else
	       rejected_cell_lowres[id_low]=true;
#endif

	   } //end if number_in cells >0
#ifdef _USE_OMP_BARR_L3_
         if(counter_multiscale == N_props_3)
               do_flag=true;
#else
           if(counter_multiscale==N_props_3) break;
#endif

       }// end loop over cells


#ifndef _USE_OMP_BARR_
#ifdef _FULL_VERBOSE_
       So.message_screen_flush("Number of properties assigned   = ",static_cast<int>(counter_multiscale));
#endif
#endif
#ifdef _USE_OMP_BARR_L3_
   }while(false==do_flag);
   N_props_3=counter_multiscale;// update in case the do loop is broken by tolerance factor
#else
       }while(counter_multiscale< N_props_3);
#endif


   cumulative_counter+=counter_multiscale;
#ifdef _FULL_VERBOSE_
   cout<<endl;
#ifdef _USE_OMP_BARR_
       So.message_screen("Number of masses assigned   = ",static_cast<int>(counter_multiscale));
#endif

   So.message_screen("Partial number of properties assigned = ",cumulative_counter);
#endif

   So.DONE();
#ifdef _FULL_VERBOSE_
   So.message_screen("Freeing memmory");
 #endif
   cell_info_cell.clear();
   cell_info_cell.shrink_to_fit();
   cells_id_still_to_assign.clear();
   cells_id_still_to_assign.shrink_to_fit();
   So.DONE();

#ifdef _USE_EXCLUSION_CELL_LEVEL3_
   rejected_cell_lowres.clear();
   rejected_cell_lowres.shrink_to_fit();
   ID_nearest_cells_lowres.clear();
   ID_nearest_cells_lowres.shrink_to_fit();
#endif

#endif

   // ******************************************************************************************************************************************
   // ******************************************************************************************************************************************
   // ******************************Level 2 ****************************************************************************************************
#ifdef _USE_MULTISCALE_LEVEL_2_


#if defined _USE_MULTISCALE_LEVEL_3_ || defined _USE_MULTISCALE_LEVEL_4_
   cell_info_cell.resize(NGRID_MS_2);
#else
   cells_id_still_to_assign.clear();
   cells_id_still_to_assign.shrink_to_fit();
#ifdef _HIGHEST_RES_LEVEL2_
   vector<ULONG>id_info(NGRID_MS_2, 0);
#else
   vector<s_cell_info>cell_info_cell(NGRID_MS_2);
#endif


#endif

   // This function returns the list of id in each ID_low (all cells included, regarless with or without mass)

#ifdef _HIGHEST_RES_LEVEL2_
   get_high_res_id_from_low_res_id(NFT_LOW_2, this->Nft, id_info);// in this case, nFT_LOW_2 > this->NFT
#else
   get_low_res_id_from_high_res_id(this->Nft, NFT_LOW_2,cell_info_cell);
#endif

#ifdef _FULL_VERBOSE_
   So.message_screen("Identifying level 3 empty cells");
#endif

   for(ULONG id_low=0; id_low< NGRID_MS_2; ++id_low)
     if(number_in_cells_aux_l2[id_low]>0)// Pinpoint the low_res cells with tracers. These tracers do not have an assigned mass
       cells_id_still_to_assign.push_back(id_low); //

#ifndef _HIGHEST_RES_LEVEL2_
   number_in_cells_aux_l2.clear();
   number_in_cells_aux_l2.shrink_to_fit();
#endif

#ifdef _USE_EXCLUSION_CELL_LEVEL2_
   So.message_screen("Identifying neighbouring cells in level 2 mesh");

#if defined _USE_EXCLUSION_CELL_LEVEL3_ || defined (_USE_EXCLUSION_CELL_LEVEL4_)
   ID_nearest_cells_lowres.resize(NGRID_MS_3);
#else
   vector<s_nearest_cells> ID_nearest_cells_lowres(NGRID_MS_2);
#endif

   get_neighbour_cells(NFT_LOW_2,NUMBER_OF_NEIGHBOUR_CELL_FOR_EXCLUSION_LEVEL2,ID_nearest_cells_lowres);
   So.DONE();
#endif

#ifdef _FULL_VERBOSE_
   So.message_screen("Shuffling order of level 2 cells containing tracers without assigned masses");
#endif
   gsl_ran_shuffle(rn,&cells_id_still_to_assign[0], cells_id_still_to_assign.size() ,sizeof(ULONG));
   So.DONE();


#ifdef _USE_EXCLUSION_CELL_LEVEL2_
#if defined  _USE_EXCLUSION_CELL_LEVEL3_ || defined (_USE_EXCLUSION_CELL_LEVEL4_)
   rejected_cell_lowres.resize(NGRID_MS_2,false);
#else
   vector<bool>rejected_cell_lowres(NGRID_MS_2,false);
#endif
#endif

#ifdef _FULL_VERBOSE_
   So.message_screen("*****Going through LEVEL 2 mesh to assign ", N_props_2, "properties");
#else
   So.message_screen("Level-2");
#endif

   counter_multiscale=0; // this is reset to zero at every leel's beginning

   loop_counter=0;  // this couts the number of for loops wr ned to perform in order to get the expected number of tracers at this level
   do_flag=false;

   do{

     ULONG N_cells=cells_id_still_to_assign.size();

#ifdef _USE_OMP_BARR_L2_
#pragma omp parallel for
#endif
     for(ULONG id_ini=0; id_ini< N_cells; ++id_ini) // loop over low-res cells containing tracers without assigned masses
       {
#ifdef _USE_OMP_BARR_L2_
         if(do_flag) continue;
#endif
         ULONG id_low = cells_id_still_to_assign[id_ini];  // Low Grid ID. These come randomized
	 // ojo que acá estan las high-res celdas que incluso no tienen tracers
#ifdef _USE_EXCLUSION_CELL_LEVEL2_
	 ULONG N_neighbour_cells=ID_nearest_cells_lowres[id_low].close_cell.size(); // Number of neighbpurs computed for the low_res
#endif
#ifdef _HIGHEST_RES_LEVEL2_
	 ULONG id=id_info[id_low];
#else
	 int N_cells_in_cells= cell_info_cell[id_low].gal_index.size();  // Number of high-res cells in a low-res cell
	 int index_cell_in_low_res_cell= gsl_rng_uniform_int(rn,N_cells_in_cells); //select randomy one of the high-id's in the low res cell
	 ULONG id=cell_info_cell[id_low].gal_index[index_cell_in_low_res_cell]; // This chosen id will be used below to assign mass
#endif

	 ULONG index_bins = cell_info_tr[id].Theta_bin; // Get the theta-bin in which this cell has been classified. This has been calculated in the step M<M*

#ifdef _HIGHEST_RES_LEVEL2_
	 if(number_in_cells_aux_l2[id_low]>0)
#else
	   if(number_in_cells_aux[id]>0)
#endif
	     {
#ifdef _USE_EXCLUSION_CELL_LEVEL2_
	       bool used_previous_cell=false;  // ask whether the previos cell is one of the neighbours of the current one
	       if(id_ini > 0)
		 {
		   for(int in = 0 ; in < N_neighbour_cells; ++in)
		     if(ID_nearest_cells_lowres[id_low].close_cell[in] == cells_id_still_to_assign[id_ini-1]) /// this is satisfied only once (at most) inside the loop!
		       {
			 used_previous_cell=true;
			 break;
		       }
		 }
	       // Proceed to assign if the previous cell is not one of the closest neighbours of the current cell
	       // or if the previus  cell id-1 was rejected because being neighbour of the id-2
	       if(false==used_previous_cell) // if the previous cell was rejected, I can use this one wven if the current cell is neighbour of the previous one
		 {
#endif  // end of use_exclusion in cells
                   ULONG N_props_without_mass_in_cell=cell_info_tr[id].gal_index.size(); //if masses are to be assigned with this, this line can be MOCK_DEN_FIELD[id];
                   ULONG jk= gsl_rng_uniform_int(rn,N_props_without_mass_in_cell); //choose a random integer in  [0,N_props_in_cell_without_mass). N_props_withopur mass must be the one fixed
		   ULONG ig=cell_info_tr[id].gal_index[jk];                       // get the Galaxy ID for the cell id in the position jk of tracers without mass
		   if(false==this->tracer.Halo[ig].observed)
		     {
		       ULONG i_mass_halo_label=count_aux_cell[index_bins];
#ifndef _ASSIGN_TO_CALIBRATION_
                       if(this->dm_properties_bins[index_bins].masses_bin_properties.size()>0) // This will be always the case when the new density filed is the same usd for the calibration
                         {
#endif
                       real_prec assigned_property=this->dm_properties_bins[index_bins].masses_bin_properties[i_mass_halo_label];
#ifdef _USE_MULTISCALE_LEVEL_3_
                       if(assigned_property >= this->Prop_threshold_multi_scale_2 && assigned_property < this->Prop_threshold_multi_scale_3)
#elif !defined (_USE_MULTISCALE_LEVEL_3_) && defined (_USE_MULTISCALE_LEVEL_4_)
                         if(assigned_property >= this->Prop_threshold_multi_scale_2 && assigned_property < this->Prop_threshold_multi_scale_4)
#elif !defined (_USE_MULTISCALE_LEVEL_3_) && !defined (_USE_MULTISCALE_LEVEL_4_)
                           if(assigned_property >= this->Prop_threshold_multi_scale_2)
#endif
			     {

#ifdef _USE_MASS_AS_OBSERVABLE_
                               this->tracer.Halo[ig].mass=assigned_property ;
#elif defined _USE_VMAX_AS_OBSERVABLE_
                               this->tracer.Halo[ig].vmax=assigned_property ;
#endif
			       this->tracer.Halo[ig].observed=true;
                               this->dm_properties_bins[index_bins].used_mass[i_mass_halo_label]=true;
#ifdef _USE_OMP_BARR_L2_
#pragma omp atomic update
#endif
			       counter_multiscale++;
#ifdef _USE_OMP_BARR_L2_
#pragma omp atomic update
#endif
			       number_in_cells_aux[id]--; // Remove one tracer from the id cell
#ifdef _USE_OMP_BARR_L2_
#pragma omp atomic update
#endif
			       count_aux_cell[index_bins]++;  // This keeps track of the masses used, and will be used below for next-scale campaign
			     }
#ifndef _ASSIGN_TO_CALIBRATION_
                        }
#endif
                     }

#ifdef _USE_EXCLUSION_CELL_LEVEL2_
		 }
	   else
		    rejected_cell_lowres[id_low]=true;
#endif

	     } //end if number_in cells >0

           if(counter_multiscale == N_props_2)
#ifdef _USE_OMP_BARR_L2_
               do_flag=true;
#else
               break;
#endif
         }// end loop over cells


#ifndef _USE_OMP_BARR_L2_
#ifdef _FULL_VERBOSE_
       So.message_screen_flush("Number of properties assigned   = ",static_cast<int>(counter_multiscale));
#endif
#endif
#ifdef _USE_OMP_BARR_L2_
   }while(false==do_flag);
   N_props_2=counter_multiscale; // update in case the do loop is broken by tolerance factor
#else
       }while(counter_multiscale < N_props_2);
#endif

   cumulative_counter+=counter_multiscale;
#ifdef _FULL_VERBOSE_
   cout<<endl;
#ifdef _USE_OMP_BARR_
       So.message_screen("Number of properties assigned   = ",static_cast<int>(counter_multiscale));
#endif
   So.message_screen("Partial number of properties assigned = ",cumulative_counter);
#endif


#ifdef _FULL_VERBOSE_
   So.DONE();
   So.message_screen("Freeing memmory in line", __LINE__);
#endif


#ifdef _HIGHEST_RES_LEVEL2_
   id_info.clear();
   id_info.shrink_to_fit();
#else
   cell_info_cell.clear();
   cell_info_cell.shrink_to_fit();
   number_in_cells_aux_l2.clear();
   number_in_cells_aux_l2.shrink_to_fit();
#endif

   cells_id_still_to_assign.clear();
   cells_id_still_to_assign.shrink_to_fit();
   So.DONE();

#ifdef _USE_EXCLUSION_CELL_LEVEL2_
   rejected_cell_lowres.clear();
   rejected_cell_lowres.shrink_to_fit();
   ID_nearest_cells_lowres.clear();
   ID_nearest_cells_lowres.shrink_to_fit();
#endif

#endif

   // ******************************************************************************************************************************************
   // ******************************************************************************************************************************************
   // ******************************************************************************************************************************************
   // ******************************************************************************************************************************************
   // ******************************Level 1 ****************************************************************************************************
#ifdef _USE_MULTISCALE_LEVEL_1_

   cells_id_still_to_assign.clear();
   cells_id_still_to_assign.shrink_to_fit();

   // *****************Assignment campaign with original reslution******************************************************************************
   //  Level 1 is always defined, unless _USE_MULTISCALE_PROPERTY_ASSIGNMENT_ is undef
   // UPDATE info of cells with tracers after having done the low-res assignment.
   // Soms of the high res cells might be empty now

#ifdef _FULL_VERBOSE_
   So.message_screen("Updating cells with particles without assigned property");
#endif
   for(ULONG id=0; id<this->NGRID; ++id)
     if(number_in_cells_aux[id]>0)// Pinpoint the cells with tracers which do not have an assigned mass
       cells_id_still_to_assign.push_back(id);
   So.DONE();

   number_in_cells_aux.clear();
   number_in_cells_aux.shrink_to_fit();


   // We randomize here the entries of cells_id_still_to assign. This is important
#ifdef _FULL_VERBOSE_
   So.message_screen("Shuffling order of cells containing tracers without assigned property");
#endif
   gsl_ran_shuffle(rn,&cells_id_still_to_assign[0],cells_id_still_to_assign.size(),sizeof(ULONG));
   So.DONE();


#ifdef _USE_EXCLUSION_CELL_LEVEL1_
   So.message_screen("Identifying neighbouring cells for level1");
   vector<s_nearest_cells> ID_nearest_cells(this->NGRID);
   get_neighbour_cells(this->Nft,NUMBER_OF_NEIGHBOUR_CELL_FOR_EXCLUSION_LEVEL1,ID_nearest_cells);
   So.DONE();
   reshuffle_once=true;
#endif

#ifdef _FULL_VERBOSE_
   So.message_screen("*****Going through LEVEL 1 mesh to assign ", N_props_1, "properties");
#else
   So.message_screen("Level-1");
#endif

#ifdef  _USE_EXCLUSION_WITH_MASS_BINS_
   // Container allocating the number of objects in a id cell in a mass bin
   vector<int>mass_bins_exc(this->NGRID*N_MASS_BINS_EXC,0);
#endif


   counter_multiscale=0;
   do_flag=false;
   do{

#ifdef _USE_OMP_BARR_L1_
#pragma omp parallel for
#endif
     for(ULONG id_ini=0; id_ini< cells_id_still_to_assign.size(); ++id_ini) // loop over cells containing tracers without assigned masses
       {

#ifdef _USE_OMP_BARR_L1_
         if(do_flag) continue;
#endif

         ULONG id = cells_id_still_to_assign[id_ini];

#ifdef _USE_EXCLUSION_CELL_LEVEL1_
	 ULONG N_neighbour_cells=ID_nearest_cells[id].close_cell.size();
	 bool used_previous_cell=false;  // ask whether the previos cell is one of the neighbours of the current one
	 if(id_ini>0)                  // This part tarkes into account the memmory from the previus assignment campaign
	   {
	     for(int in = 0 ; in < N_neighbour_cells; ++in)
	       if(ID_nearest_cells[id].close_cell[in] == cells_id_still_to_assign[N_neighbour_cells-1]) /// this is satisfied only once (at most) inside the loop!
		        {
		          used_previous_cell=true;
		          break;
		         }
	   }

	 // Proceed to assign if the previous cell is not one of the closest neighbours of the current cell
	 // or if the previus  cell id-1 was rejected because being neighbour of the id-2
	 if(false==used_previous_cell) // if the previous cell was rejected, I can use this one wven if the current cell is neighbour of the previous one
	   {
#endif  // end of use_exclusion in cells
	     ULONG index_bins = cell_info_tr[id].Theta_bin; // Get the theta-bin in which this cell has been classified. This has been calculated in the step M<M*
	     ULONG N_mocks_without_mass_in_cell=cell_info_tr[id].gal_index.size();
       ULONG jk= gsl_rng_uniform_int(rn,N_mocks_without_mass_in_cell); //choose a random integer in  [0,N_props_in_cell_without_mass). N_props_withopur mass must be the one fixed
	     ULONG ig=cell_info_tr[id].gal_index[jk];        // get the Galaxy ID for the cell id in the position jk of tracers without mass
	     ULONG i_mass_halo_label=count_aux_cell[index_bins];
	     if(false==this->tracer.Halo[ig].observed)
	       {
#ifndef _ASSIGN_TO_CALIBRATION_
            if(this->dm_properties_bins[index_bins].masses_bin_properties.size()>0)// this if must be added a theta bin filled from the mock might be empty in the ref
               {
#endif
                       real_prec assigned_property = this->dm_properties_bins[index_bins].masses_bin_properties[i_mass_halo_label];
#ifdef _USE_EXCLUSION_WITH_MASS_BINS_
                        ULONG index_mass_bin=get_bin(assigned_property,log10(this->prop_min),N_MASS_BINS_EXC,log10(this->prop_max/this->prop_min)/static_cast<double>(N_MASS_BINS_EXC),true);
                        ULONG index_mass_id=index_2d(id, index_mass_bin,N_MASS_BINS_EXC);
                        if(mass_bins_exc[index_mass_id]< MAX_NUMBER_OF_MASSES_IN_CELL_SAME_MBIN) // allow only three masses from the same mass bin in a cell
                          {
#endif
#if defined (_USE_MULTISCALE_LEVEL_2_)
                            if(assigned_property >= this->Prop_threshold_multi_scale_1 && assigned_property < this->Prop_threshold_multi_scale_2)
#elif defined (_USE_MULTISCALE_LEVEL_2_) && defined (_USE_MULTISCALE_LEVEL_3_) && defined (_USE_MULTISCALE_LEVEL_4_)
                            if(assigned_property >= this->Prop_threshold_multi_scale_1 && assigned_property < this->Prop_threshold_multi_scale_2)
#elif !defined (_USE_MULTISCALE_LEVEL_2_) && defined (_USE_MULTISCALE_LEVEL_3_) && defined (_USE_MULTISCALE_LEVEL_4_)
                            if(assigned_property >= this->Prop_threshold_multi_scale_1 && assigned_property < this->Prop_threshold_multi_scale_3)
#elif !defined (_USE_MULTISCALE_LEVEL_2_) && !defined (_USE_MULTISCALE_LEVEL_3_) && defined (_USE_MULTISCALE_LEVEL_4_)
                            if(assigned_property >= this->Prop_threshold_multi_scale_1 && assigned_property < this->Prop_threshold_multi_scale_4)
#elif !defined (_USE_MULTISCALE_LEVEL_2_) && !defined (_USE_MULTISCALE_LEVEL_3_) && !defined (_USE_MULTISCALE_LEVEL_4_)
                            if(assigned_property >= this->Prop_threshold_multi_scale_1)
#endif
                              {
#ifdef _USE_MASS_AS_OBSERVABLE_
                               this->tracer.Halo[ig].mass=assigned_property ;
#elif defined _USE_VMAX_AS_OBSERVABLE_
#ifdef _ASSIGN_MASS_POST_
                               if(true==initial_assignment)
#endif
                                  this->tracer.Halo[ig].vmax=assigned_property ;
#ifdef _ASSIGN_MASS_POST_
                               else
                                  this->tracer.Halo[ig].mass=assigned_property ;
#endif
#endif

#ifdef _USE_OMP_BARR_L1_
#pragma omp atomic update
#endif
                                 counter_multiscale++;
#ifdef _USE_OMP_BARR_L1_
#pragma omp atomic update
#endif
                                 count_aux_cell[index_bins]++;
                                 this->tracer.Halo[ig].observed=true;
                                 this->dm_properties_bins[index_bins].used_mass[i_mass_halo_label] = true;  //mark this mass as already assigned
#ifdef _USE_EXCLUSION_WITH_MASS_BINS_
                          mass_bins_exc[index_mass_id]++;
#endif
                                }

#ifdef _USE_EXCLUSION_WITH_MASS_BINS_
                       }
#endif
#ifndef _ASSIGN_TO_CALIBRATION_
                     }// if no available particles in this theta bin, use the global mass function
#endif

#if !defined _ASSIGN_TO_CALIBRATION_

/* This line was origianlly meant to be
#if !defined _MASS_ASSIGNMENT_TO_REFERENCE_ || !defined _ASSIGN_TO_CALIBRATION_
Nevertheless, it is not celar why the || is not working ghere
*/

                 else // if no particles available, use the global mass function. THis should not happen if assignmet to reference or to same Ncount field from calibration
                      {
                        real_prec prob=-10.0;
                        real_prec ran=10.0;
                        real_prec mass_tracer;
                        while(prob<ran)
                          {
                            real_prec xr = static_cast<real_prec>(gsl_rng_uniform(rn));
                            real_prec fraction_mass = xr*log10(this->prop_max/this->prop_min);
                            mass_tracer=pow(10,log10(this->prop_min)+fraction_mass);
                            real_prec aux_mass=0;
#ifdef _USE_MASS_AS_OBSERVABLE_
                            aux_mass= (mass_tracer <= this->tracer_ref.MBmin[0]? this->tracer_ref.MBmin[0]: mass_tracer);
                            aux_mass= (mass_tracer >= this->tracer_ref.MBmin[nbins_mf-1]? this->tracer.MBin[nbins_mf-1]: mass_tracer);
#elif defined _USE_VMAX_AS_OBSERVABLE_
                              aux_mass= (mass_tracer <= this->tracer_ref.VMAXBmin[0]? this->tracer_ref.VMAXBmin[0]: mass_tracer);
                              aux_mass= (mass_tracer >= this->tracer_ref.VMAXBmin[nbins_mf-1]? this->tracer.VMAXBin[nbins_mf-1]: mass_tracer);
#endif
                              if(aux_mass<this->Prop_threshold_multi_scale_1)
                                prob=0.0;
                              else
                                {
                                  ran  = static_cast<real_prec>(gsl_rng_uniform(rn));
                                  prob = static_cast<real_prec>(gsl_spline_eval (spline, aux_mass, acc));
                                }
                              mass_tracer=aux_mass;
                            }
#ifdef _USE_MASS_AS_OBSERVABLE_
                            this->tracer.Halo[ig].mass=mass_tracer;
#elif defined _USE_VMAX_AS_OBSERVABLE_
                            this->tracer.Halo[ig].vmax=mass_tracer;
#endif
                            this->tracer.Halo[ig].observed=true;
#ifdef _USE_OMP_BARR_L1_
#pragma omp atomic update
#endif
                            counter_multiscale++;
                       }
#endif   // END IFNDEF_MASS_ASSIGNMENT_TO_REFERENCE_

                 }// end if false==observed
#ifdef _USE_EXCLUSION_CELL_LEVEL1_
	         }
#endif

          if(counter_multiscale == N_props_1)
#ifdef _USE_OMP_BARR_L1_
              do_flag=true;
#else
              break;
#endif
      }// end loop over cells


#ifndef _USE_OMP_BARR_L1_
#ifdef _FULL_VERBOSE_
       So.message_screen_flush("Number of properties assigned   = ",static_cast<int>(counter_multiscale));
#endif
#endif

#ifdef _USE_EXCLUSION_CELL_LEVEL1_
     if(reshuffle_once==true)
       if(counter_multiscale>= N_props_1-N_props_ABOVE_THRESHOLD)
	      {
	         gsl_ran_shuffle(rn,&cells_id_still_to_assign[0],cells_id_still_to_assign.size(),sizeof(ULONG));
	   //	     rejected_cell.resize(this->NGRID,false);
	         reshuffle_once=false;
	       }
#endif

#ifdef _USE_EXCLUSION_WITH_MASS_BINS_
     if(counter_multiscale>=N_props_1-N_props_ABOVE_M_EX)
       break;
#endif
#ifdef _USE_OMP_BARR_L1_
   }while(false==do_flag);   // do untill all masses have acquisted mass
   N_props_1=counter_multiscale; // reassigne in case the do loop was stoped because of non being able to get all requested numbers
#else
    }while(counter_multiscale < N_props_1);   // do untill all masses have acquisted mass
#endif


   cumulative_counter+=counter_multiscale;
#ifdef _USE_EXCLUSION_WITH_MASS_BINS_
   N_props_0=this->tracer_ref.NOBJS-((N_props_1-N_props_ABOVE_M_EX)+N_props_2+N_props_3+N_props_4);
#endif

#ifdef _FULL_VERBOSE_
   cout<<endl;
#ifdef _USE_OMP_BARR_
       So.message_screen("Number of properties assigned   = ",static_cast<int>(counter_multiscale));
#endif
   So.message_screen("Partial number of properties assigned = ",cumulative_counter);
   cout<<endl;
   So.message_screen("Freeing memory in line", __LINE__);
#endif   // end of level 1


   cells_id_still_to_assign.clear();
   cells_id_still_to_assign.shrink_to_fit();
   So.DONE();

#endif   // end of _USE_MULTISCALE_LEVEL_1_

/*
#ifndef _MASS_ASSIGNMENT_TO_REFERENCE_
   So.message_screen("Freeing memory");
   this->dm_properties_bins.clear(); this->dm_properties_bins.shrink_to_fit();
   So.DONE();
#endif
*/
// ***************************************************************************************************************
// ***************************************************************************************************************
// ***************************************************************************************************************
// **********************************************Level 0 *********************************************************
// ***************************************Assign at particle level for X<X threshold *****************************

#ifdef _USE_VMAX_AS_OBSERVABLE_
#ifdef _FULL_VERBOSE_
  So.message_screen("*****Assigning Vmax at level 0 (particle by particle based on {Theta} properties)");
#else
   So.message_screen("Level-0");
 #endif

#ifdef _MASS_ASSIGNMENT_TO_REFERENCE_
   So.message_screen("Available = ", N_props_0);
#endif

#else
   So.message_screen("Assigning Mhalo at particle level. Available = ", N_props_0);
#endif

   ULONG counter_masses_b=0; // Number of tracers without assigned property
   ULONG counter_masses_0=0; // NUmber of tracers with property assigned: This will be updated inside the loop
   ULONG counter_from_ref=0; // Number of tracers read from the reference

#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:counter_masses_b)
#endif
   for(ULONG ig=0; ig < this->tracer.NOBJS;++ig)
     if(false==this->tracer.Halo[ig].observed)
       counter_masses_b++;
#ifdef _FULL_VERBOSE_
   So.message_screen("Requested(from non-observed tracer)  = ", counter_masses_b);
   So.message_screen("Requested(check Nobs-cumulative) = ", this->tracer.NOBJS-cumulative_counter);
#endif

   ULONG counter_masses_b_tolerance=static_cast<ULONG>(counter_masses_b*TOLERANCE_FACTOR_L0);
#ifdef _FULL_VERBOSE_
   So.message_screen("Requested with tolerance factor = ", counter_masses_b_tolerance);
#endif


   int bin_threshold=get_bin(log10(this->Prop_threshold_multi_scale_1), log10(this->params._VMAXmin()),this->params._NMASSbins(),log10(this->params._VMAXmax()/this->params._VMAXmin())/static_cast<double>(this->params._NMASSbins()), false);

   // count if there are still available masses from the reference:
  // If we still have, then we select them from dm_properties_bins. If not, we use the global mass function.
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:counter_from_ref)
#endif
   for(ULONG i = 0; i < this->dm_properties_bins.size(); ++i)
       for(ULONG j = 0; j < this->dm_properties_bins[i].used_mass.size(); ++j)
         if(false==this->dm_properties_bins[i].used_mass[j])
            counter_from_ref++;
#ifdef _FULL_VERBOSE_
   So.message_screen("Number of tracers available form reference = ", counter_from_ref);
#endif

#ifdef _MASS_ASSIGNMENT_TO_REFERENCE_
   if(N_props_0>0)
#else
     if(counter_masses_b>0) // If we need to assign, proceed
#endif
       {



#ifdef _USE_OMP_BARR_L0_
          int jthread=0;
          vector<ULONG>vseeds(NTHREADS,0);
          for(int i=0;i<vseeds.size();++i)
             vseeds[i]=35+static_cast<ULONG>(i+14)*56045;

           const gsl_rng_type *Trn;
           gsl_rng *rna;

#pragma omp parallel private (jthread, Trn, rna)// this is causing problem
        {
#endif


#if !defined _ASSIGN_TO_CALIBRATION_
         volatile bool do_flag=false;
         ULONG aux_counter=0;
         ULONG aux_counter_prev=0;
         int aux_repetition=0;

         do{ // DO this part untill all available tracers in ref are used. When sampling from the DM of the calibration, this do-loop is not needed.
#endif




#ifdef _USE_OMP_BARR_L0_
#pragma omp for
#endif
           for(ULONG ig=0; ig < this->tracer.NOBJS;++ig)
           {
#ifndef _ASSIGN_TO_CALIBRATION_
              if(do_flag) continue;
#endif
              if(false==this->tracer.Halo[ig].observed)
                {
                 ULONG id=this->tracer.Halo[ig].GridID;
                 ULONG index_bins = cell_info_tr[id].Theta_bin;
                 if(counter_from_ref>counter_masses_0) // counter_masses_0 is updated inside this loop. This "if" statment ensures that we try to use all refernece properties
                   {
                     ULONG N_props_in_bin = this->dm_properties_bins[index_bins].masses_bin_properties.size();
                     bool flag=false;
#ifndef _ASSIGN_TO_CALIBRATION_
                     if(N_props_in_bin>0) // This will be always the case when the new density filed is the same usd for the calibration
                       {
#endif
                         while(flag == false) // este while obliga a elegir una masa -no elegida antes- de entre las disponibles
                          {
                            int i_mass_halo_label= gsl_rng_uniform_int(rn,N_props_in_bin);
                            real_prec assigned_property=dm_properties_bins[index_bins].masses_bin_properties[i_mass_halo_label] ;
                            bool used_prop = this->dm_properties_bins[index_bins].used_mass[i_mass_halo_label]; //true or false if the mass was already chosen or not

#ifndef _ASSIGN_TO_CALIBRATION_
                            // When assigning properties to the tracer field generated with the DM field used in the calibration, we do not have conflicts here.
                            // However, when sampling a new DM field, it can happen that we reach this point with used_prop=true (which is ead form the reference)
                            // and this->tracer.Halo[ig].observed=false (which is read from the newly sampled tracer number count)
                            // thus inducing anb infinite while loop (because flag will be always flag=false).
                            // Hence, here we ask if used_prop is true (i.e, this property linked to a particle in this theta bins has been already used),
                            // in which case we set the flag=true and leave the while loop without assigning anything
                          if(true==used_prop)
                             flag=true;
                          else//         If used_prop=false (the property has not been used before), then assign that mass to the current particle i
                            {
#else
                          if(false==used_prop) // if the prop has not been used before, then assign that mass to the current particle i
                            {
#endif

#ifdef _USE_MASS_AS_OBSERVABLE_
                             this->tracer.Halo[ig].mass=assigned_property ;
#elif defined _USE_VMAX_AS_OBSERVABLE_
                             this->tracer.Halo[ig].vmax=assigned_property ;
#endif
                             this->tracer.Halo[ig].observed=true;
                             this->dm_properties_bins[index_bins].used_mass[i_mass_halo_label] = true; //mark this mass as already assigned
                             flag = true;
#ifdef _USE_OMP_BARR_L0_
#pragma omp atomic update
#endif
                             counter_masses_0++;
                            }
                          }  // close while
#ifndef _ASSIGN_TO_CALIBRATION_
                         } // closes  if(N_props_in_bin>0)
#endif
                     } //closes counter_from_ref>counter_mass_0
#ifndef _ASSIGN_TO_CALIBRATION_
                 else
                   {
                    real_prec aux_h=0;//
                    for(int iy=0;iy< bin_threshold+1;++iy) // we only need to go up to the threshold mass
                       aux_h+=this->ABUNDANCE_normalized[index_bins+iy*LENGHTdm];

                   if(aux_h>0)
                     {
                       bool flag=false;
                       while(false==flag)
                         {
                           real_prec prob=-10.0;
                           real_prec ran=10.0;
                           int i_mass_halo_index;
                           while(prob<ran)
                             {
                               i_mass_halo_index= gsl_rng_uniform_int(rn,bin_threshold+1); // we only need to go up to the threshold mass
                               ULONG index_or=index_bins+i_mass_halo_index*LENGHTdm;//     index_2d(i_mass_halo_index,index_bins,LENGHTdm);
                               prob = this->ABUNDANCE_normalized[index_or];
                               ran = static_cast<real_prec>(gsl_rng_uniform(rn));
                             }
                           real_prec xr = static_cast<real_prec>(gsl_rng_uniform(rn));
#ifdef _USE_MASS_AS_OBSERVABLE_
                           real_prec fraction_mass= xr*log10(this->tracer_ref.MBmax[i_mass_halo_index]/this->tracer_ref.MBmin[i_mass_halo_index]);
                           real_prec lmass_halo = log10(this->tracer_ref.MBmin[i_mass_halo_index])+fraction_mass ;
                           real_prec mass_assigned=pow(10,lmass_halo);
                           if(mass_assigned<this->Prop_threshold_multi_scale_1) // just to be sure
                             {
                               this->tracer.Halo[ig].mass = mass_assigned;
                               this->tracer.Halo[ig].observed=true;
                               counter_masses_0++;
                               flag=true;
                             }
#elif  defined _USE_VMAX_AS_OBSERVABLE_
                           real_prec fraction_mass= xr*log10(this->tracer_ref.VMAXBmax[i_mass_halo_index]/this->tracer_ref.VMAXBmin[i_mass_halo_index]);
                           real_prec lmass_halo = log10(this->tracer_ref.VMAXBmin[i_mass_halo_index])+fraction_mass ;
                           real_prec prop_assigned=pow(10,lmass_halo);
                           if(prop_assigned<this->Prop_threshold_multi_scale_1) // just to be sure that we are below the threshold of multilevel 1
                             {
                               this->tracer.Halo[ig].vmax = prop_assigned;
                               this->tracer.Halo[ig].observed=true;
                               flag=true;
#ifdef _USE_OMP_BARR_L0_
#pragma omp atomic update
#endif
                               counter_masses_0++;
                             }
#endif
                         }

                     }
               } //closes else
#endif

          }// closes if(false==this->tracer.Halo[ig].observed)



#ifdef _FULL_VERBOSE_
#ifdef _USE_MASS_AS_OBSERVABLE_
          So.message_screen_flush("Number of masses assigned   = ",static_cast<int>(counter_masses_0));
#elif defined _USE_VMAX_AS_OBSERVABLE_
#ifndef _USE_OMP_BARR_L0_
              So.message_screen_flush("Number of properties assigned   = ", counter_masses_0);
#endif
#endif
#endif


#if !defined _ASSIGN_TO_CALIBRATION_
          if(counter_masses_0 == counter_masses_b_tolerance)
            do_flag=false;
#endif

        }// end loop over tracers

#if !defined _ASSIGN_TO_CALIBRATION_
         }while(counter_masses_0 < counter_masses_b_tolerance); // this will be never exact, for the inner loop cannot comunicate outside the current nu ber of accepted `
#endif

#ifdef _USE_OMP_BARR_L0_
             }  // closes parallel region
#endif

         } // end if N_props_0 > 0
     So.DONE();

#ifdef _FULL_VERBOSE_
#ifdef _USE_OMP_BARR_
#ifdef _USE_VMAX_AS_OBSERVABLE_
     So.message_screen("Number of Vmax assigned   = ", counter_masses_0);
#elif defined _USE_MASS_AS_OBSERVABLE_
     So.message_screen("Number of Masses assigned   = ", counter_masses_0);
#endif

#endif
#endif

     cumulative_counter+=counter_masses_0;
#ifdef _FULL_VERBOSE_
   cout<<endl;
   So.message_screen("Partial number of properties assigned = ",cumulative_counter);
   So.message_screen_flush("Number of properties to assign = ",this->tracer.NOBJS-cumulative_counter);
   cout<<endl;
#endif


   if(this->tracer.NOBJS-cumulative_counter<0)
       {

   counter_from_ref=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:counter_from_ref)
#endif
    for(ULONG i = 0; i < this->dm_properties_bins.size(); ++i)
       for(ULONG j = 0; j < this->dm_properties_bins[i].used_mass.size(); ++j)
         if(false==this->dm_properties_bins[i].used_mass[j])
            counter_from_ref++;
      // When sampling the same DM field used for the calibration, at this point there should be no available tracers in the reference from which properties can be read.
#ifdef _FULL_VERBOSE_
    So.message_screen("CHECK: available number of proparties to be assigned, read from the reference = ", counter_from_ref);
    cout<<endl;
#endif


// If we are assigning to the DM used for calibration, things should end up here. However, at this point we still have tracers without propertes and lefts overs from the ref catalog in cae the DM is another
// so we can simply assign the missing ones as randomly from the lefts overs of the refereece
// In oder ot to that, we collect in a single arraay the lefst overs
#ifdef _FULL_VERBOSE_
    So.message_screen("Assigning remaining properties randomly selecting from the remainig set of the reference:");
#endif
    vector<real_prec>left_overs_properties;
    vector<bool>left_overs_properties_used;

    for(ULONG i = 0; i < this->dm_properties_bins.size(); ++i)
       for(ULONG j = 0; j < this->dm_properties_bins[i].used_mass.size(); ++j)
          if(false==this->dm_properties_bins[i].used_mass[j])
            {
              left_overs_properties.push_back(this->dm_properties_bins[i].masses_bin_properties[j]  );
              left_overs_properties_used.push_back(false);
             }
#ifdef _FULL_VERBOSE_
    So.message_screen("CHECK: size of container for these objects = ", left_overs_properties.size());
    cout<<endl;
#endif
    ULONG nc_test=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:nc_test)
#endif
    for(ULONG ig=0; ig < this->tracer.NOBJS;++ig)
       if(false==this->tracer.Halo[ig].observed)
          nc_test++;
#ifdef _FULL_VERBOSE_
    So.message_screen("Number of properties to assign, taken from from 'non-observed' tracers = ", nc_test);
    So.message_screen("Number of tracers to assign with tolerance = ", static_cast<ULONG>(nc_test*TOLERANCE_FACTOR_L0b));
#endif

    if(nc_test*TOLERANCE_FACTOR_L0b >left_overs_properties.size())// In this case, we have more requested than available
     {
#ifdef _FULL_VERBOSE_
        So.message_screen("Changing the number to be assigned to the number of available refs");
#endif
       nc_test=left_overs_properties.size();
     }

    ULONG counter_masses_left=0;
    do{
       for(ULONG ig=0; ig < this->tracer.NOBJS;++ig)
         {
           if(false==this->tracer.Halo[ig].observed)
             {
              int jj = gsl_rng_uniform_int(rn,left_overs_properties.size());
              if(false==left_overs_properties_used[jj])
                 {
#ifdef _USE_MASS_AS_OBSERVABLE_
                   this->tracer.Halo[ig].mass=left_overs_properties[jj] ;
#elif defined _USE_VMAX_AS_OBSERVABLE_
                   this->tracer.Halo[ig].vmax= left_overs_properties[jj] ;
#endif
                   this->tracer.Halo[ig].observed=true;
                   this->tracer.Halo[ig].multi_scale_level=8;
                   left_overs_properties_used[jj]=true;
                   counter_masses_left++;
                   counter_masses_0++;
                 }
              }
           }
#ifdef _FULL_VERBOSE_
           So.message_screen_flush("Assigned from left-overs= ", counter_masses_left);
#endif
           if(counter_masses_left==static_cast<ULONG>(TOLERANCE_FACTOR_L0b*nc_test))break;
        }while(counter_masses_left < static_cast<ULONG>(TOLERANCE_FACTOR_L0b*nc_test));
#ifdef _FULL_VERBOSE_
      So.DONE();
      So.message_screen("Total Assigned at Level 0 = ", counter_masses_0);
#endif
      cumulative_counter+=counter_masses_left;
#ifdef _FULL_VERBOSE_
     So.message_screen("Partial number of properties assigned = ",cumulative_counter);
#endif
#ifdef _FULL_VERBOSE_
     cout<<endl;
      So.message_screen("Freeing memory in line", __LINE__);
#endif
      left_overs_properties.clear();
      left_overs_properties.shrink_to_fit();
      left_overs_properties_used.clear();
      left_overs_properties_used.shrink_to_fit();
      So.DONE();

      // *************************** end of the left-opver assignment, if requested

   }//  close if(this->tracer.NOBJS-cumulative_counter<0)


     N_props_0=counter_masses_0;  // update this number. Asked below

// **********************************************END Level 0 *********************************************************
#endif // end for _USE_MULTISCALE_PROPERTY_ASSIGNMENT_


//#ifndef _MASS_ASSIGNMENT_TO_REFERENCE_  // this was moved from a slot before the loop, because I changed a ndef by a def inside the looip in orer to use dm _prop
#ifdef _FULL_VERBOSE_
   So.message_screen("Freeing memory in line", __LINE__);
#endif
   this->dm_properties_bins.clear();
   this->dm_properties_bins.shrink_to_fit();


   So.DONE();
//#endif

#if defined _USE_VMAX_AS_OBSERVABLE_ || defined _USE_MASS_AS_OBSERVABLE_ 
    }  // close if(true==initial_assignment): 
#endif



#ifdef _USE_VMAX_AS_OBSERVABLE_
  else
#endif
   if (false==initial_assignment)//to assign MASSES based on the information of  v max already assigned
    {
      ULONG N_props_0=this->tracer.NOBJS;
#ifdef _FULL_VERBOSE_
      if(h_property==_MASS_)
          So.message_screen("Assigning Mvir.  Expected = ", this->tracer.NOBJS);
      else if(h_property==_RS_)
          So.message_screen("Assigning RS.  Expected = ", this->tracer.NOBJS);
      else if(h_property==_SPIN_)
          So.message_screen("Assigning SPIN.  Expected = ", this->tracer.NOBJS);
#endif

#ifndef test_vmax_mass
      vector<ULONG>number_in_theta_ref(LENGHTdm,0);
      // container to track the number of available masses in each theta-bin
#pragma omp parallel for
      for(ULONG i=0;i<LENGHTdm; ++i)
	number_in_theta_ref[i]=this->dm_properties_bins[i].masses_bin_properties.size();
#endif

      ULONG counter_masses_0=0;
      ULONG counter_masses_or=0;
      ULONG counter_masses_global=0;

      ULONG Ntot=0;
#ifdef _USE_VMAX_AS_OBSERVABLE_
      Ntot=N_VMAX_BINS;
#endif


    if(h_property==_MASS_)
        {
#ifdef _add_dm_density_
        Ntot*=this->NX;
#ifdef _add_Xweb_
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
        Ntot*=N_C_BIN1;  // We are assuming that Xweb
#endif

#ifdef _USE_INVARIANT_TIDAL_FIELD_III
         Ntot*=N_C_BIN2  // We are assuming that Xweb
#endif
#endif
#endif
   }
 else if(h_property==_RS_)
   Ntot*=N_MASS_BINS;
 else if(h_property==_SPIN_)
  Ntot*=N_MASS_BINS;


// The following loops are meant to inizialize the tracer containers
  if(h_property==_MASS_)
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG ig=0; ig < this->tracer.NOBJS;++ig)
        {
          this->tracer.Halo[ig].observed=false;
          this->tracer.Halo[ig].mass=0;
        }

  if(h_property==_RS_)
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG ig=0; ig < this->tracer.NOBJS;++ig)
        {
          this->tracer.Halo[ig].observed=false;
          this->tracer.Halo[ig].rs=0;
        }

  if(h_property==_SPIN_)
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG ig=0; ig < this->tracer.NOBJS;++ig)
        {
          this->tracer.Halo[ig].observed=false;
          this->tracer.Halo[ig].spin=0;
        }


 // ************************************************************************************
 // ************************************************************************************
 //  The following loop aims at measuriung P(M|Vmax, theta) .
 // ************************************************************************************
 // ************************************************************************************


#ifdef _USE_OMP_TEST2_
  vector<ULONG>vseeds(NTHREADS,0);
     for(int i=0;i<vseeds.size();++i)
         vseeds[i]=554565+static_cast<ULONG>(i)*565;

     int jthread=0;
     const gsl_rng_type *Trn;
     gsl_rng *randa;

#pragma omp parallel private(jthread, randa, Trn)
    {
         jthread=omp_get_thread_num();
         gsl_rng_default_seed=vseeds[jthread];
         Trn = gsl_rng_default;
         randa = gsl_rng_alloc (Trn);

#pragma omp parallel for reduction(+:counter_masses_0, counter_masses_or, counter_masses_global)
#endif
   for(ULONG ig=0; ig < this->tracer.NOBJS;++ig)
       {
         ULONG id=this->tracer.Halo[ig].GridID;

          int I_X=0;
#ifdef _add_dm_density_
          if(h_property==_MASS_)
           {
             real_prec xdm = static_cast<real_prec>(this->delta_X[id]);
             I_X  = get_bin(xdm,this->s_mins.prop1,this->NX,this->s_deltas.prop1,this->bin_accumulate_borders);
          }
          else
#endif
          if(h_property==_RS_ || h_property==_SPIN_)
             I_X= get_bin(log10(this->tracer.Halo[ig].mass),this->params._LOGMASSmin(),N_MASS_BINS, (this->params._LOGMASSmax()-this->params._LOGMASSmin())/static_cast<double>(N_MASS_BINS)  ,this->bin_accumulate_borders);

#ifdef _add_Xweb_
          int I_C1=0;
          int I_C2=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
          real_prec C1 = this->cwclass.Invariant_TF_II[id];
          I_C1=get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#endif
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
          real_prec C2 = this->cwclass.Invariant_TF_III[id];
          I_C2=get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#endif

#endif


#ifndef test_vmax_mass
          ULONG index_bins = cell_info_tr[id].Theta_bin;
          ULONG N_props_left =  number_in_theta_ref[index_bins];
          ULONG N_props_in_bin = this->dm_properties_bins[index_bins].masses_bin_properties.size();
          if(N_props_left>0)
            {
              bool flag=false;
              while(flag == false) // este while obliga a elegir una masa -no elegida antes- de entre las disponibles
                {
                  int i_mass_halo_label= gsl_rng_uniform_int(rn,N_props_in_bin);
                  bool used_mass = this->dm_properties_bins[index_bins].used_mass[i_mass_halo_label]; //true or false if the mass was already chosen or not
                  if(false==used_mass)// if the mass has not been used before, then assign that mass to the current particle i
                    {
                      this->tracer.Halo[ig].mass=dm_properties_bins[index_bins].masses_bin_properties[i_mass_halo_label] ;
                      this->tracer.Halo[ig].observed=true;
                      this->dm_properties_bins[index_bins].used_mass[i_mass_halo_label] = true; //mark this mass as already assigned
                      flag = true;
#ifdef _USE_OMP_TEST2_
#pragma omp atomic
                      number_in_theta_ref[index_bins]--;
#endif
                      counter_masses_0++;
                    }
                }
            }
          else
            {
#endif

              ULONG I_CV2=0;
#ifdef test_vmax_mass
              I_CV2=get_bin(log10(this->tracer.Halo[ig].vmax), log10(this->params._VMAXmin()),N_VMAX_BINS,log10(this->params._VMAXmax()/this->params._VMAXmin())/static_cast<double>(N_VMAX_BINS), true);
#endif

              ULONG index_dm=0;

              if(h_property==_MASS_)
               {
#ifdef _add_Xweb_
                 index_dm=index_4d(index_bins,I_X,I_C1, I_C2, this->NX, N_C_BIN1,N_C_BIN2); //I_CV2 takes here the place of Vmax
#endif
#if defined (_add_dm_density_) && !defined (_add_Xweb_)
                 index_dm=index_2d(I_CV2,I_X,this->NX);
#endif

#ifndef _add_dm_density_
                 index_dm=I_CV2;
#endif
               }
              else if(h_property==_RS_ || h_property==_SPIN_)
                index_dm=index_2d(I_CV2,I_X,N_MASS_BINS);

	      // Control loop over bins of mass
              double aux_h=0;
              for(int iy=0;iy<this->params._NMASSbins();++iy)
#ifdef test_vmax_mass
                aux_h+=static_cast<double>(this->ABUNDANCE_normalized[index_2d(iy,index_dm,Ntot)]);
#else
                aux_h+=this->ABUNDANCE_normalized[index_2d(iy,index_bins,LENGHTdm)];
#endif

              // here we assign mass according to the P(M|Vmax,ð) relation obtained from the reference
              if(aux_h>0)
                {
                  real_prec prob=-10.0;
                  real_prec ran=10.0;
                  int i_halo_index;
                  while(prob<ran)
                    {
#ifdef _USE_OMP_TEST2_
                      i_halo_index= gsl_rng_uniform_int(randa,this->params._NMASSbins()); // we have chosen to use NMASSBIns also for RS, Vmax etc when these are to be assigned
#else 
                      i_halo_index= gsl_rng_uniform_int(rn,this->params._NMASSbins()); // we have chosen to use NMASSBIns also for RS, Vmax etc when these are to be assigned
#endif

#ifdef test_vmax_mass
                      ULONG index_or=index_2d(i_halo_index,index_dm,Ntot);
#else
                      ULONG index_or=index_2d(i_mass_halo_index,index_bins,LENGHTdm);
#endif
              		    prob = this->ABUNDANCE_normalized[index_or];
#ifdef _USE_OMP_TEST2_
                      ran = gsl_rng_uniform(randa);
#else
                      ran = gsl_rng_uniform(rn);
#endif
                    }
#ifdef _USE_OMP_TEST2_
                    real_prec xr = static_cast<real_prec>(gsl_rng_uniform(randa));
#else
                    real_prec xr = static_cast<real_prec>(gsl_rng_uniform(rn));
#endif
                    if(_MASS_ == h_property)
                      {
                        real_prec fraction_mass= xr*log10(this->tracer_ref.MBmax[i_halo_index]/this->tracer_ref.MBmin[i_halo_index]);
                        real_prec lmass_halo = log10(this->tracer_ref.MBmin[i_halo_index])+fraction_mass ;
                        this->tracer.Halo[ig].mass = pow(10,lmass_halo);
                      }
                    else if(_RS_== h_property)
                      {
                        real_prec fraction_mass= xr*log10(this->tracer_ref.RSBmax[i_halo_index]/this->tracer_ref.RSBmin[i_halo_index]);
                        real_prec lmass_halo = log10(this->tracer_ref.RSBmin[i_halo_index])+fraction_mass ;
                        this->tracer.Halo[ig].rs = pow(10,lmass_halo);
                      }
                    else if(_SPIN_== h_property)
                      {
                        real_prec fraction_mass= xr*log10(this->tracer_ref.SPINBmax[i_halo_index]/this->tracer_ref.SPINBmin[i_halo_index]);
                        real_prec lmass_halo = log10(this->tracer_ref.SPINBmin[i_halo_index])+fraction_mass ;
                        this->tracer.Halo[ig].spin = pow(10,lmass_halo);
                      }

                    this->tracer.Halo[ig].observed=true;
              		  counter_masses_0++;
		                counter_masses_or++;
		              }
#ifndef _MASS_ASSIGNMENT_TO_REFERENCE_
                  else{
		                real_prec prob=-10.0;
		                real_prec ran=10.0;
		                real_prec mass_tracer;
                    // here we assign mass according to the measured abundance
                    while(prob<ran)
                		  {
#ifdef _USE_OMP_TEST2_
                        real_prec xr = static_cast<real_prec>(gsl_rng_uniform(randa));
#else
                        real_prec xr = static_cast<real_prec>(gsl_rng_uniform(rn));
#endif
                        real_prec fraction_mass = xr*log10(this->prop_max/this->prop_min);
                        mass_tracer=pow(10,log10(this->prop_min)+fraction_mass);
                        if(h_property==_MASS_)
                          {
                            real_prec aux_mass= (mass_tracer <= this->tracer_ref.MBin[0]? this->tracer_ref.MBin[0]: mass_tracer);
                            aux_mass= (mass_tracer >= this->tracer_ref.MBin[nbins_mf-1]? this->tracer_ref.MBin[nbins_mf-1]: mass_tracer);
                            if(aux_mass<=this->tracer_ref.MBin[0] || aux_mass >= this->tracer_ref.MBin[nbins_mf-1])
                              prob=0.0;
                            else
                              {
#ifdef _USE_OMP_TEST2_
                                ran  = static_cast<real_prec>(gsl_rng_uniform(randa));
#else
                                ran  = static_cast<real_prec>(gsl_rng_uniform(rn));
#endif
                                prob = static_cast<real_prec>(gsl_spline_eval (spline, aux_mass, acc));
                               }
                            mass_tracer=aux_mass;
                            this->tracer.Halo[ig].mass=mass_tracer;
                            this->tracer.Halo[ig].observed=true;
                          }
                         else if(h_property==_RS_)
                          {
                            real_prec aux_mass= (mass_tracer <= this->tracer_ref.RSBin[0]? this->tracer_ref.RSBin[0]: mass_tracer);
                            aux_mass= (mass_tracer >= this->tracer_ref.RSBin[nbins_mf-1]? this->tracer_ref.RSBin[nbins_mf-1]: mass_tracer);
                            if(aux_mass<=this->tracer_ref.RSBin[0] || aux_mass >= this->tracer_ref.RSBin[nbins_mf-1])
                              prob=0.0;
                            else
                             {
#ifdef _USE_OMP_TEST2_
                               ran  = static_cast<real_prec>(gsl_rng_uniform(randa));
#else
                               ran  = static_cast<real_prec>(gsl_rng_uniform(rn));
#endif
                               prob = static_cast<real_prec>(gsl_spline_eval (spline, aux_mass, acc));
                             }
                            mass_tracer=aux_mass;
                            this->tracer.Halo[ig].rs=mass_tracer;
                            this->tracer.Halo[ig].observed=true;
                          }
                         else if(h_property==_SPIN_)
                          {
                            real_prec aux_mass= (mass_tracer <= this->tracer_ref.SPINBin[0]? this->tracer_ref.SPINBin[0]: mass_tracer);
                            aux_mass= (mass_tracer >= this->tracer_ref.SPINBin[nbins_mf-1]? this->tracer_ref.SPINBin[nbins_mf-1]: mass_tracer);
                            if(aux_mass<=this->tracer_ref.SPINBin[0] || aux_mass >= this->tracer_ref.SPINBin[nbins_mf-1])
                              prob=0.0;
                            else
                             {
#ifdef _USE_OMP_TEST2_
                               ran  = static_cast<real_prec>(gsl_rng_uniform(randa));
#else 
                               ran  = static_cast<real_prec>(gsl_rng_uniform(rn));
#endif 
                               prob = static_cast<real_prec>(gsl_spline_eval (spline, aux_mass, acc));
                            }
                            mass_tracer=aux_mass;
                            this->tracer.Halo[ig].spin=mass_tracer;
                            this->tracer.Halo[ig].observed=true;
                          }
                       }
                      counter_masses_global++;
                      counter_masses_0++;
	                }
#endif   //end of ifndef _MASS_ASSIGNMENT_TO_REFERENCE_

#ifndef test_vmax_mass
            }
#endif

#ifdef _FULL_VERBOSE_
#ifndef _USE_OMP_TEST2_
          So.message_screen_flush("Number of properties assigned   = ",static_cast<int>(counter_masses_0));
#endif
#endif
      } // end loop over tracers


#ifdef _USE_OMP_TEST2_
     gsl_rng_free(randa);
    } // closes parallel region
 #endif



#ifdef _FULL_VERBOSE_
    cout<<endl;

#ifdef _USE_OMP_
#ifdef _USE_MASS_AS_OBSERVABLE_
          So.message_screen("Number of masses assigned   = ",static_cast<int>(counter_masses_0));
#elif defined _USE_VMAX_AS_OBSERVABLE_
#ifdef _FULL_VERBOSE_
        So.message_screen("Number of property values assigned = ",static_cast<int>(counter_masses_0));
#endif
#endif
#endif

    if(h_property==_MASS_)
      So.message_screen("Number of Mvir values assigned (using joint prob-dist) = ",static_cast<int>(counter_masses_or));
    else if(h_property==_RS_)
       So.message_screen("Number of Rs values assigned (using joint prob-dist) = ",static_cast<int>(counter_masses_or));
    else if(h_property==_SPIN_)
      So.message_screen("Number of Spin values assigned (using joint prob-dist) = ",static_cast<int>(counter_masses_or));

    So.message_screen("Number of properties assigned with global abundance = ",static_cast<int>(counter_masses_global));

    So.DONE();
    cout<<endl;
#endif

#ifndef test_vmax_mass
    number_in_theta_ref.clear();
    number_in_theta_ref.shrink_to_fit();
#endif

    N_props_0=counter_masses_0;
  } // closes  else if (false==initial_assignment)


#ifdef _FULL_VERBOSE_
  So.message_screen("Freeing memmory in line",__LINE__);
#endif
  this->ABUNDANCE_normalized.clear(); this->ABUNDANCE_normalized.shrink_to_fit();
#ifndef test_vmax_mass
  this->dm_properties_bins.clear();
  this->dm_properties_bins.shrink_to_fit();
#endif
#ifdef _USE_TRACERS_IN_CELLS_
  MOCK_DEN_FIELD.clear();
  MOCK_DEN_FIELD.shrink_to_fit();
#endif
  So.DONE();



#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_
  ULONG assigned_propertyes_mf=N_props_0+N_props_1+N_props_2+N_props_3+N_props_4;
#else
  ULONG assigned_propertyes_mf=counter_fmf;
#endif


  ULONG Ncount_ab=0;


  if(true==initial_assignment)  //this applies for vmax assignment when the global vmax function  is needed
    {
      if(assigned_propertyes_mf  < this->tracer.NOBJS)
      	{
#ifdef _FULL_VERBOSE_
          So.message_screen("Number of tracers with *no property* assigned so far = ", this->tracer.NOBJS-assigned_propertyes_mf);
#endif

          ULONG nc_test=0;
#pragma parallel for reduction(+:nc_test)
          for(ULONG ig=0; ig < this->tracer.NOBJS; ++ig)
            if(false==this->tracer.Halo[ig].observed)
              nc_test++;



#ifdef _FULL_VERBOSE_
          So.message_screen("Number of tracers with *no property* assigned so far, from 'non-observed tracer' = ", nc_test);
          cout<<endl;
          So.message_screen("Assigning property with global *property* function");
#endif
          NTHREADS=_NTHREADS_;
          omp_set_num_threads(NTHREADS);

#ifdef _USE_OMP_TEST_
          int jthread=0;
          vector<ULONG>vseeds(NTHREADS,0);
          for(int i=0;i<vseeds.size();++i)
             vseeds[i]=35+static_cast<ULONG>(i+14)*56045;
#endif
          const gsl_rng_type *Trn;
          gsl_rng *rna;

#ifdef _USE_OMP_TEST_
#pragma omp parallel private (jthread, Trn, rna)// this is causing problem
        {
              jthread=omp_get_thread_num();
              gsl_rng_default_seed=vseeds[jthread];
#else
             gsl_rng_default_seed=1555;
#endif
            Trn = gsl_rng_mt19937;//_default;
            rna= gsl_rng_alloc (Trn);


#ifdef _USE_OMP_TEST_
#pragma omp for reduction(+:Ncount_ab) // this parallelization was causing a problem
#endif
            for(ULONG ig=0; ig < this->tracer.NOBJS; ++ig)
                {
                  if(false==this->tracer.Halo[ig].observed)
                      {
                         real_prec prob=-10.0;
                         real_prec ran=10.0;
                         real_prec property_tracer;
                         while(prob<ran)
                            {
                              real_prec xr = static_cast<real_prec>(gsl_rng_uniform(rna));
                              real_prec fraction_prop = xr*log10(this->prop_max/this->prop_min);
                              property_tracer=pow(10,log10(this->prop_min)+fraction_prop);
                              real_prec aux_prop= (property_tracer < this->prop_min? this->prop_min : property_tracer);
                              aux_prop= (property_tracer > this->prop_max? this->prop_max: property_tracer);

#ifdef _USE_VMAX_AS_OBSERVABLE_
                              if(aux_prop<this->tracer_ref.VMAXBin[0] || aux_prop>this->tracer_ref.VMAXBin[nbins_mf-1])
#elif defined _USE_MASS_AS_OBSERVABLE_
                              if(aux_prop<=this->tracer_ref.MBmin[0] || aux_prop>=this->tracer_ref.MBmin[nbins_mf-1])
#endif
                                 prob=0.0;
                              else
                                {
                                  ran  = static_cast<real_prec>(gsl_rng_uniform(rna));
//                                  cout<<aux_prop<<endl;
                                  prob = static_cast<real_prec>(gsl_spline_eval (spline, aux_prop, acc));
                                }
                              property_tracer=aux_prop;
                              }
#ifdef _USE_MASS_AS_OBSERVABLE_
                	 	        this->tracer.Halo[ig].mass=property_tracer;
#elif defined _USE_VMAX_AS_OBSERVABLE_
                            this->tracer.Halo[ig].vmax=property_tracer;
#endif
                            this->tracer.Halo[ig].observed=true;
                            this->tracer.Halo[ig].multi_scale_level=-1;
                            Ncount_ab++;
                      }//closes if
               }// closes loop over tracers
          gsl_rng_free(rna);
#ifdef _USE_OMP_TEST_
             }  // closes parallel region
#endif
          So.DONE();
#ifdef _FULL_VERBOSE_
      So.message_screen("Number of props assigned using global abundance function =",Ncount_ab);
      cout<<endl;
#endif
            }
    }
#ifndef test_vmax_mass
  else  // if false==initial_assignment, so this section will assign the mass from the mass function
    { // global mass assignment once vmax was done
      if(assigned_propertyes_mf  < this->tracer.NOBJS)
      	{
          So.message_screen("Number of tracers with no property assigned = ", this->tracer.NOBJS-assigned_propertyes_mf);
          cout<<endl;
          So.message_screen("Assigning mass with global *property* function");
          int counter_m=0;
      	  for(ULONG ig=0; ig < this->tracer.NOBJS; ++ig)
	         {
	           if(false==this->tracer.Halo[ig].observed)
          		{
		            real_prec prob=-10.0;
          		  real_prec ran=10.0;
		            real_prec mass_tracer;
          		  while(prob<ran)
		              {
		                real_prec xr = static_cast<real_prec>(gsl_rng_uniform(rn));
          		      real_prec fraction_mass = xr*log10(this->prop_max/this->prop_min);
		                mass_tracer=pow(10,log10(this->prop_min)+fraction_mass);
		                real_prec aux_mass= (mass_tracer <= this->prop_min ? this->prop_min : mass_tracer);
		                aux_mass= (mass_tracer >= this->prop_max ? this->prop_max: mass_tracer);
          		      if(aux_mass<  this->prop_min  || aux_mass>=this->prop_max)
			                 prob=0.0;
		                else
                			{
                  		  ran  = static_cast<real_prec>(gsl_rng_uniform(rn));
               	  		  prob = static_cast<real_prec>(gsl_spline_eval (spline, aux_mass, acc));
			                }
	                  mass_tracer=aux_mass;
		              }
           		  this->tracer.Halo[ig].mass=mass_tracer;
          		  counter_m++;
		          }
#ifdef _FULL_VERBOSE_
            So.message_screen_flush("Number of masses assigned   = ",static_cast<int>(counter_m));
#endif

  	      }
#ifdef _FULL_VERBOSE_
	  So.DONE();
#endif
      	}
    }
#endif

  gsl_rng_free (rn);
  cumulative_counter+=Ncount_ab;

  NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);



  string fname_mass_function_Y = this->Output_directory+"tracer_mock_abundance.txt";

#ifdef _USE_VMAX_AS_OBSERVABLE_

  this->tracer.aux_flag=false; // This avoids the calculation of vmax threshods from the mocks, since it is not necessary


  if(true==initial_assignment)
    {
      this->tracer.params.i_vmax_g=5; //this allows the tracer to get the vmax function
      this->tracer.params.i_mass_g=-5; //no mass information here
      this->tracer.get_property_function(fname_mass_function_Y);

      real_prec residuals=0;
      int count_b=0;
#pragma omp parallel for reduction(+:count_b,residuals)
      for(int i=0;i<this->params._NMASSbins_mf() ;++i)
      	if(this->tracer.vmax_function[i]>0)
	       {
    	     count_b++;
	         residuals+= fabs(this->tracer_ref.vmax_function[i]/this->tracer.vmax_function[i] -1.0 );
     	   }
        residuals/=static_cast<real_prec>(count_b)/100.0;
#ifdef _FULL_VERBOSE_
        So.message_screen("Residuals from  vmax-function = ", residuals, "%");
        cout<<endl;
#endif
    }
  else
    {
     if(h_property==_MASS_)
      {
        this->tracer.params.i_vmax_g=-5; //this allows the abvid the measurement of the vmax function again
        this->tracer.params.i_mass_g=4; //this allow the tracer to get the mass function
        this->tracer.params.i_rs_g=-5; //no rs information here
        this->tracer.params.i_spin_g=-6; //this allow the tracer to get the mass function
        this->tracer.get_property_function(fname_mass_function_Y);
        real_prec residuals=0;
        int count_b=0;
#pragma omp parallel for reduction(+:count_b,residuals)
        for(int i=0;i<this->params._NMASSbins_mf() ;++i)
        	if(this->tracer.mass_function[i]>0)
        	  {
       	      count_b++;
	            residuals+= fabs(this->tracer_ref.mass_function[i]/this->tracer.mass_function[i]-1.0);
        	  }
          residuals/=static_cast<real_prec>(count_b)/100.0;
#ifdef _FULL_VERBOSE_
          So.message_screen("Residuals from mass-function = ", residuals, "%");
          cout<<endl;
#endif
      }

    else if(h_property==_RS_)
      {
        this->tracer.params.i_vmax_g=-5; //this allows the avoid the measurement of the vmax function again
        this->tracer.params.i_mass_g=-4; //this allow aviud the tracer to get the mass function
        this->tracer.params.i_rs_g=5; //this allow the tracer to get the mass function
        this->tracer.params.i_spin_g=-6; //this allow the tracer to get the mass function
        this->tracer.get_property_function(fname_mass_function_Y);
        real_prec residuals=0;
        int count_b=0;
#pragma omp parallel for reduction(+:count_b,residuals)
        for(int i=0;i<this->params._NMASSbins_mf() ;++i)
          if(this->tracer.rs_function[i]>0)
            {
              count_b++;
              residuals+= fabs(this->tracer_ref.rs_function[i]/this->tracer.rs_function[i]-1.0);
            }
        residuals/=static_cast<real_prec>(count_b)/100.0;
#ifdef _FULL_VERBOSE_
        So.message_screen("Residuals from Rs-function = ", residuals, "%");
        cout<<endl;
#endif
      }
      else if(h_property==_SPIN_)
       {
        this->tracer.params.i_vmax_g=-5; //this allows the avoid the measurement of the vmax function again
        this->tracer.params.i_mass_g=-4; //this allow aviud the tracer to get the mass function
        this->tracer.params.i_rs_g=-5; //this allow the tracer to get the mass function
        this->tracer.params.i_spin_g=6; //this allow the tracer to get the mass function
        this->tracer.get_property_function(fname_mass_function_Y);
        real_prec residuals=0;
        int count_b=0;
#pragma omp parallel for reduction(+:count_b,residuals)
        for(int i=0;i<this->params._NMASSbins_mf() ;++i)
          if(this->tracer.s_function[i]>0)
            {
              count_b++;
              residuals+= fabs(this->tracer_ref.s_function[i]/this->tracer.s_function[i]-1.0);
            }
        residuals/=static_cast<real_prec>(count_b)/100.0;
#ifdef _FULL_VERBOSE_
        So.message_screen("Residuals from Spin-function = ", residuals, "%");
        cout<<endl;
#endif
      }
    }
#elif defined _USE_MASS_AS_OBSERVABLE_
  real_prec residuals=0;

      this->tracer.get_property_function(fname_mass_function_Y);

#pragma omp parallel for reduction(+:residuals)
  for(int i=0;i<this->params._NMASSbins_mf() ;++i)
    if(this->tracer.mass_function[i]>0)
      residuals+= fabs(this->tracer_ref.mass_function[i]/this->tracer.mass_function[i]-1.0);
  residuals/=static_cast<real_prec>(this->params._NMASSbins_mf())/100.0;
  So.message_screen("Residuals from mass-function = ", residuals, "%");
  cout<<endl;
#endif


}
#endif
#endif // nd if mock mode

//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################

   //  ####################################################################################################################################################################
   //  ####################################################################################################################################################################
   //  ####################################################################################################################################################################
   //  ####################################################################################################################################################################
   //  ####################################################################################################################################################################
   //  ####################################################################################################################################################################

   void Bam::sample_mock()
   {

       So.enter(__PRETTY_FUNCTION__);
     this->So.message_screen("Sampling halo density field with random positions within cells");
     const gsl_rng_type *  T;
     gsl_rng * r ;

     gsl_rng_env_setup();
     gsl_rng_default_seed=1152;
     T = gsl_rng_ranlux;
     r = gsl_rng_alloc (T);

     string pname;
     if(this->step <=this->N_iterations_Kernel)
       pname ="_iteration"+to_string(this->step);
     else
       pname= "_realization"+to_string(this->params._realization());
     string fname=this->Output_directory+"MOCK_TR"+pname+"_CAT_z"+to_string(this->redshift)+".txt";

     vector<real_prec>x_cart;
     vector<real_prec>y_cart;
     vector<real_prec>z_cart;

     real_prec delta=this->Lbox/static_cast<real_prec>(this->Nft);
     for(ULONG i=0;i<this->delta_Y_new.size();++i)
       {
         int Ngal=this->delta_Y_new[i];
         int counter=0;
         ULONG x_min, y_min, z_min;
         index2coords(i,this->Nft, x_min, y_min, z_min);
         x_min*=delta;
         y_min*=delta;
         z_min*=delta;

         while(counter<=Ngal)
           {
             x_cart.push_back(x_min+delta*gsl_rng_uniform (r));
             y_cart.push_back(y_min+delta*gsl_rng_uniform (r));
             z_cart.push_back(z_min+delta*gsl_rng_uniform (r));
             counter++;
           }
       }
     So.message_screen("Sampled with",x_cart.size(),"objects");
     this->File.write_to_file(fname, x_cart,y_cart,z_cart);
     So.DONE();

   }






   //  ####################################################################################################################################################################
   //  ####################################################################################################################################################################
   //  ####################################################################################################################################################################
   //  ####################################################################################################################################################################
   //  ####################################################################################################################################################################
   //  ####################################################################################################################################################################

#ifdef _SEVERAL_REAL_
   void Bam::get_pdf(int real)
#elif !defined _SEVERAL_REAL_
     void Bam::get_pdf()
#endif
   {

      So.enter(__PRETTY_FUNCTION__);

     this->So.message_screen("BAM in BIAS mode");
     this->So.message_screen("Computing statistics from input density fields");

     // The int sua indicates 0, 1, 2.,  ... n_cwt usados. SI tomamos knots y el resto, sua =0, 1
     string type_X=this->Scale_X;
     string type_Y=this->Scale_Y;

     // real_prec num_in_log_x = true==this->Convert_Density_to_Delta_X ? NUM_IN_LOG: 0.;
     // real_prec num_in_log_y = true==this->Convert_Density_to_Delta_Y ? NUM_IN_LOG: 0.;

     real_prec num_in_log_x=0;
     if(type_X=="log")
       num_in_log_x = NUM_IN_LOG;

     real_prec num_in_log_y=0;
     if(type_Y=="log")
       num_in_log_y = NUM_IN_LOG;



#ifdef _SEVERAL_REAL_
     string prop_X="_X"+this->XNAME+"_ScaleX"+type_X+"_MASX"+to_string(this->iMAS_X)+"_Realization"+to_string(real);
     string prop_Y="Y"+this->YNAME+"_ScaleY"+type_Y+"_MASY"+to_string(this->iMAS_Y);
#else
     string prop_X="_X"+this->XNAME+"_ScaleX"+type_X+"_MASX"+to_string(this->iMAS_X);
     string prop_Y="Y"+this->YNAME+"_ScaleY"+type_Y+"_MASY"+to_string(this->iMAS_Y)+"_"+this->params._extra_info();
#endif

     // ******************************************************************************

#ifdef _USE_CWC_
     this->cwclass.do_CWC(this->delta_X_ini); // CWClass done with the full delta
#ifdef _USE_MASS_KNOTS_
     this->cwclass.get_Mk_collapsing_regions(this->delta_X_ini,  static_cast<real_prec>(this->N_objects_X)/static_cast<real_prec>(this->NGRID));
#endif

#elif !defined _USE_CWC_
#ifdef _USE_MASS_KNOTS_
     this->cwclass.do_CWC(this->delta_X_ini); // CWClass done with the full delta
     this->cwclass.get_Mk_collapsing_regions(this->delta_X_ini,  static_cast<real_prec>(this->N_objects_X)/static_cast<real_prec>(this->NGRID));
#endif

#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_PWEB_)
      this->cwclass.do_CWC(this->delta_X_ini);
#endif

#endif

#if (defined (_USE_NABLA2DELTA_) || defined (_USE_S2_)  || defined (_USE_S3_) || defined (_USE_S2DELTA_)) && (!defined (_USE_CWC_)) && !(defined (_USE_PWEB_))
     this->cwclass.get_bias_terms(this->delta_X_ini);
#endif


     vector<real_prec> delta_Y_ini(this->NGRID,0);
     delta_Y_ini=this->delta_Y;
     vector<real_prec> delta_X_ini(this->NGRID,0);
     delta_X_ini=this->delta_X;
     // ******************************************************************************


     int sua=0;
#ifdef _USE_CWC_
     for(sua=0;sua<this->n_cwt;++sua) // Loop over the different CWC requested in the parameter file
       {
#endif

         this->tstruct=0;
#ifdef _USE_CWC_
         this->tstruct=this->cwclass.cwt_used[sua];
#endif

//         string file=this->Output_directory+"2dbias"+prop_X+"_"+prop_Y+"_Nft"+to_string(this->Nft)+"_z"+to_string(this->redshift)+"_LambdaTH"+to_string(lambdath)+"_CW"+to_string(this->tstruct);
         string file=this->Output_directory+"2dbias"+prop_X+"_"+prop_Y+"Nft"+to_string(this->Nft);

         // ******************************************************************************
         // So.message_screen("Computing variances for cosmic web type ",this->tstruct);
         // real_prec var_XX=gsl_stats_variance(&this->delta_X[0],1, this->NGRID);
         // real_prec var_YY=gsl_stats_variance(&this->delta_Y[0],1, this->NGRID);
         // real_prec corr_XY=gsl_stats_correlation(&this->delta_X[0],1, &this->delta_Y[0], 1, this->NGRID);
         // So.message_screen("Variance DM =",var_XX);
         // So.message_screen("Variance TR =",var_YY);
         // So.message_screen("Correlation X-Y =",corr_XY);

         this->get_min_max_X_Y();// THIS CAN BE REPLACED by finding min and max at converting time

         // ******************************************************************************
         if(this->tstruct!=0)
           {
             // Here we need to redefine the overdensities according to the mean for each CWT
#ifdef _FULL_VERBOSE_
             So.message_screen("Transforming to overdensities for CWT",this->tstruct);
#endif



             ULONG new_nobjects=0;
             real_prec beta=0;

             if(true==this->Convert_Density_to_Delta_Y)
               {
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:new_nobjects)
#endif
                 for(ULONG i = 0; i < this->NGRID; ++i)
#ifdef _USE_CWC_
                   if(true==this->cwclass.get_cell_classified(sua, i))
#endif
                     new_nobjects+=(this->N_objects_Y/static_cast<real_prec>(this->NGRID))*(1.0+delta_Y_ini[i]);

                 beta=static_cast<real_prec>(new_nobjects)/static_cast<real_prec>(this->N_objects_Y);



#ifdef _USE_OMP_
#pragma omp parallel for
#endif
                 for(ULONG i = 0; i < this->NGRID; ++i)
                   this->delta_Y[i]=(1./beta)*delta_Y_ini[i]-(beta-1)/beta;

                 So.message_screen("New nobs",  new_nobjects);
                 So.message_screen("Original nobs",  this->N_objects_Y);

                 So.message_screen("Factor beta for TR", beta);
               }


             if(true==this->Convert_Density_to_Delta_X)
               {
                 new_nobjects=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:new_nobjects)
#endif
                 for(ULONG i = 0; i < this->NGRID; ++i)
#ifdef _USE_CWC_
                   if(true==this->cwclass.get_cell_classified(sua, i))
#endif
                     new_nobjects+=(this->N_objects_X/static_cast<real_prec>(this->NGRID))*(1.0+delta_X_ini[i]);

                 beta=static_cast<real_prec>(new_nobjects)/static_cast<real_prec>(this->N_objects_X);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
                 for(ULONG i = 0; i < this->NGRID; ++i)
                   this->delta_X[i]=(1./beta)*this->delta_X[i]-(beta-1)/beta;
                 So.message_screen("Factor beta for DM", beta);

               }
           }

         // ******************************************************************************
         this->So.message_screen("Minimum of delta X", get_min(delta_X));
         this->So.message_screen("Maximum of delta X", get_max(delta_X));
         this->So.message_screen("Minimum of delta Y", get_min(delta_Y));
         this->So.message_screen("Maximum of delta Y", get_max(delta_Y));

         if(type_Y=="log")
           {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
            for(ULONG i = 0; i < this->NGRID; ++i)
#ifdef _USE_CWC_
               if(true==this->cwclass.get_cell_classified(sua, i))
#endif
                 this->delta_Y[i]=  log10(num_in_log_y+static_cast<real_prec>(delta_Y_ini[i]));
#ifdef _USE_CWC_
               else
                 delta_Y[i]=NOCELL;
#endif
           }
#ifdef _USE_CWC_
         else
           {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
             for(ULONG i = 0;i < this->NGRID ;++i)
               if(false!=this->cwclass.get_cell_classified(sua, i))
                 delta_Y[i]=NOCELL;
           }
#endif




         if(type_X=="log")
           {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
             for(ULONG i = 0;i < this->NGRID ;++i)
               {
#ifdef _USE_CWC_
                 if(true==this->cwclass.get_cell_classified(sua, i))
                   {
#endif
                     delta_X[i]= log10(num_in_log_x+static_cast<real_prec>(delta_X_ini[i]));

#ifdef _USE_CWC_
                   }
                 else
                   delta_X[i]=NOCELL;
#endif
               }
           }
#ifdef _USE_CWC_
         else
           {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
             for(ULONG i = 0;i < this->NGRID ;++i)
               if(false!=this->cwclass.get_cell_classified(sua, i))
                 delta_X[i]=NOCELL;
           }
#endif




         if(type_Y=="log")
           {
             this->So.message_screen("Minimum of log10 1 + delta Y", get_min(delta_Y));
             this->So.message_screen("Maximum of log10 1 + delta Y", get_max(delta_Y));
           }


         if(type_X=="log")
           {
             this->So.message_screen("Minimum of log10 1 + delta X", get_min(delta_X));
             this->So.message_screen("Maximum of log10 1 + delta X", get_max(delta_X));
           }


         // ******************************************************************************
#ifdef _FULL_VERBOSE_
         So.message_screen("Getting X and Y bins used for the bias information");
#endif
         // ******************************************************************************
         // Aca redefinimos los bines en X y Y.
         // Hacemos una distinición importante. Cuando usamos NGP,
         // los bines estarán identificados con el número de particulas
         // Para otros tipos de interpolaciones, hacemos bines propiamente.
         // nmax_x_onecell is the maximum number of objects in one cell.

         if(this-> iMAS_X==0)
           {
             if(false==Convert_Density_to_Delta_X)
               if(type_X=="linear")
                 {
                   this->new_nbins_x = this->nmax_x_onecell+1;
                   this->Xmin=0;
                   this->Xmax=static_cast<int>(nmax_x_onecell);
                 }
               else{
                 this->new_nbins_x = this->NX;
                 this->Xmin=this->ldelta_X_min;
                 this->Xmax=this->ldelta_X_max;
               }
             else
               {
                 if(type_X=="linear")
                   {
                     this->new_nbins_x =this->NX;
                     this->Xmin=this->delta_X_min;
                     this->Xmax=this->delta_X_max;
                   }
                 else{
                   this->new_nbins_x = this->NX;
                   this->Xmin=this->ldelta_X_min;
                   this->Xmax=this->ldelta_X_max;
                 }
               }
           }


         else
           {
             if(true==Convert_Density_to_Delta_X)
               {
                 if(type_X=="linear")
                   {
                     this->new_nbins_x = this->NX;
                     this->Xmin=this->delta_X_min;
                     this->Xmax=this->delta_X_max;
                   }
                 else if(type_X=="log")
                   {
                     this->new_nbins_x = this->NX;
                     this->Xmin=this->ldelta_X_min;
                     this->Xmax=this->ldelta_X_max;
                   }
               }

             else
               {
                 if(type_X=="linear")
                   {
                     this->new_nbins_x = this->NX;
                     this->Xmin=this->delta_X_min;
                     this->Xmax=this->delta_X_max;
                   }
                 else if(type_X=="log")
                   {
                     this->new_nbins_x = this->NX;
                     this->Xmin=this->ldelta_X_min;
                     this->Xmax=this->ldelta_X_max;
                   }
               }
           }





         // ******************************************************************************
         if(this-> iMAS_Y==0)
           {
             if(false==Convert_Density_to_Delta_Y)
               if(type_Y=="linear")
                 {
                   this->new_nbins_y = this->nmax_y_onecell+1;
                   this->Ymin=0;
                   this->Ymax=static_cast<int>(nmax_y_onecell);
                 }
               else{
                 this->new_nbins_y = this->NY;
                 this->Ymin=this->ldelta_Y_min;
                 this->Ymax=this->ldelta_Y_max;
               }
             else
               {
                 if(type_Y=="linear")
                   {
                     this->new_nbins_y =this->NY;
                     this->Ymin=this->delta_Y_min;
                     this->Ymax=this->delta_Y_max;
                   }
                 else{
                   this->new_nbins_y = this->NY;
                   this->Ymin=this->ldelta_Y_min;
                   this->Ymax=this->ldelta_Y_max;
                 }
               }
           }


         else{  // if CIC or anything higher
           if(true==Convert_Density_to_Delta_Y)
             {
               if(type_Y=="linear")
                 {
                   this->Ymin=this->delta_Y_min;
                   this->Ymax=this->delta_Y_max;
                 }
               else if(type_Y=="log")
                 {
                   this->Ymin=this->ldelta_Y_min;
                   this->Ymax=this->ldelta_Y_max;
                 }
             }
           else
             {
               if(type_Y=="linear")
                 this->Ymax=this->delta_Y_max;
               else if(type_Y=="log")
                 this->Ymax=this->ldelta_Y_max;
               if(type_Y=="linear")
                 this->Ymin=this->delta_Y_min;
               else if(type_Y=="log")
                 this->Ymin=this->ldelta_Y_min;
           }
           this->new_nbins_y = this->NY;

         }
         So.DONE();


/*
         this->Ymax=get_max(this->delta_Y);
         this->Ymin=get_min(this->delta_Y);
         this->Xmax=get_max(this->delta_X);
         this->Xmin=get_min(this->delta_X);
*/

         // ******************************************************************************
         this->DELTAX=(this->Xmax-this->Xmin)/static_cast<real_prec>(this->NX);
         this->DELTAY=(this->Ymax-this->Ymin)/static_cast<real_prec>(this->new_nbins_y);
         cout<<DELTAX<<endl;
#ifdef _nomas_
         this->So.message_screen("Minimum  Y", this->Ymin);
         this->So.message_screen("Maximum  Y", this->Ymax);
         this->So.message_screen("Minimum  X", this->Xmin);
         this->So.message_screen("Maximum  X", this->Xmax);
#endif
         // ******************************************************************************
         // ******************************************************************************

         this->X_bins.resize(new_nbins_x,0);


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
         for(ULONG i = 0; i < new_nbins_x; ++i)
           X_bins[i]= (iMAS_X == 0 && false==this->Convert_Density_to_Delta_X && type_X=="linear") ? static_cast<real_prec>(i) :  this->Xmin+(i+0.5)*(this->Xmax-this->Xmin)/static_cast<real_prec>(new_nbins_x);


         ofstream sx; sx.open("xbins.txt");
         for(ULONG i = 0; i < new_nbins_x; ++i)
           sx<<i<<"  "<<X_bins[i]<<endl;
         sx.close();

         this->Y_bins.resize(new_nbins_y,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
//         for(ULONG i = 0; i < new_nbins_y; ++i)
 //          Y_bins[i]= (iMAS_Y == 0  && false==this->Convert_Density_to_Delta_Y  && type_Y=="linear") ? static_cast<real_prec>(i) :  this->Ymin+(i+0.5)*(this->Ymax-this->Ymin)/static_cast<real_prec>(new_nbins_y);

         for(ULONG i = 0; i < new_nbins_y; ++i)
           Y_bins[i]=   this->Ymin+(i+0.5)*(this->Ymax-this->Ymin)/static_cast<real_prec>(new_nbins_y);


         ofstream sy; sy.open("ybins.txt");
         for(ULONG i = 0; i < new_nbins_y; ++i)sy<<i<<"  "<<Y_bins[i]<<endl;
         sy.close();


         // ******************************************************************************
         // ************************RESIZE VECTORS FOR HISTOGRAMS*************************
         // ******************************************************************************
         // Resize the vector to allocate the number of cells in the Den-Nh bins

         // This vectors will contain P(X,Y) only used for analyizing , not to create mocks


         BIAS_NCOUNTS.resize(this->n_sknot_massbin*this->n_cwt*new_nbins_x*new_nbins_y*N_REDSHIFT_BINS, 0);
         BIAS_NCOUNTS_normalized.resize(this->n_sknot_massbin*this->n_cwt*new_nbins_x*new_nbins_y*N_REDSHIFT_BINS, 0);

         // ******************************************************************************
         // ******************************************************************************
#ifdef _FULL_VERBOSE_
         So.message_screen("Computing BIAS(X,Y, CWT, ...) ");
#endif


#ifdef _USE_REDSHIFT_MASK_
         vector<real_prec> redshift_mask(this->NGRID,0);
         this->File.read_array(this->Name_redshift_mask,redshift_mask);

         cout<<"Max redshift in mask = "<<get_max(redshift_mask)<<endl;
         cout<<"Min redshift in mask = "<<get_min(redshift_mask)<<endl;


         vector<int> redshift_mask_bin(this->NGRID,0);
         vector<int> ncells_zbin(N_REDSHIFT_BINS,0);
#endif


#ifdef _USE_BINARY_MASK_
         vector<real_prec> binary_mask(this->NGRID,0);
         this->File.read_array(this->Name_binary_mask, binary_mask);

#endif


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
         for(ULONG ig = 0; ig< this->NGRID ; ++ig)
           {

             int redshift_bin=0;
             real_prec redshift_mask_cell=0;
#ifdef _USE_REDSHIFT_MASK_
#ifdef _USE_BINARY_MASK_
             if(static_cast<real_prec>(binary_mask[ig])>0)
               {
#endif
                 redshift_mask_cell=redshift_mask[ig];

                 if(redshift_mask_cell<REDSHIFT_MAX && redshift_mask_cell>=REDSHIFT_MIN)
                   {
                     redshift_bin = static_cast<int>(floor((redshift_mask_cell-REDSHIFT_MIN)/DELTA_Z));
                     redshift_mask_bin[ig]=redshift_bin;


#pragma omp atomic update
                     ncells_zbin[redshift_bin]++;

#endif

                     int ict=0;
#ifdef _USE_MASS_KNOTS_
                     ict=this->cwclass.SKNOT_M_info[ig];
#endif
#ifdef _USE_CWC_
                     if(true==this->cwclass.get_cell_classified(sua, ig))
                       {
#endif
                         if(delta_X[ig]!= NOCELL && delta_Y[ig]!= NOCELL)   // This means, use cells being classified.
                           if((delta_X[ig]>=this->Xmin && delta_X[ig]<=this->Xmax) && (delta_Y[ig]>=this->Ymin && delta_Y[ig]<=this->Ymax))
                             {
                               int IX= get_bin(delta_X[ig], this->Xmin, this->NX,this->DELTAX, false);
                               int IY= get_bin(delta_Y[ig], this->Ymin, this->NY,this->DELTAY, false);

#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                               this->BIAS_NCOUNTS[index_5d(IX, IY, sua, ict, redshift_bin, new_nbins_y, this->n_cwt, this->n_sknot_massbin, N_REDSHIFT_BINS)]++;


                           }
#ifdef _USE_CWC_
                       }
#endif

#ifdef _USE_REDSHIFT_MASK_
                   }
#ifdef _USE_BINARY_MASK_
               }//close if binary mask
#endif

#endif
           }





#ifdef _USE_REDSHIFT_MASK_
         ULONG NCELLS=0;
         for(int i=0;i<ncells_zbin.size();++i)
           NCELLS+=ncells_zbin[i];

         cout<<"Used "<<NCELLS<<" out of "<<this->NGRID<<"  ("<<100.*abs((this->NGRID-NCELLS)/static_cast<real_prec>(this->NGRID))<<" %)"<<endl;

         redshift_mask.clear();
         redshift_mask.shrink_to_fit();
#endif


         // DO THIS ONLY FOR KNOTS, SINCE THIS IS A LOOP OVER THE MKNOTS
         //Now I normalize the 4d array such that I can also do contours for the knots in different bins of MK


         int iz=0;
#ifdef _USE_REDSHIFT_MASK_
         for(iz=0;iz<N_REDSHIFT_BINS;++iz)
           {
#endif

             int mk=0;
#ifdef _USE_MASS_KNOTS_
             for(mk=0;mk<this->n_sknot_massbin;++mk)
               {
#endif
                 vector<LONG>aux_v(new_nbins_x*new_nbins_y, 0);
                 for(int lx=0;lx<new_nbins_x;++lx)
                   for(int ly=0;ly<new_nbins_y;++ly)
                     aux_v[index_2d(lx,ly,new_nbins_y)]=this->BIAS_NCOUNTS[index_5d(lx,ly,sua,mk,iz,new_nbins_y, this->n_cwt, this->n_sknot_massbin, N_REDSHIFT_BINS)];

                 ULONG lkk=get_max<LONG>(aux_v);

                 vector<real_prec>aux_n(new_nbins_x*new_nbins_y, 0);
                 for(int i=0;i< aux_v.size() ;++i)
                   aux_n[i]= lkk==0 ? 0 : static_cast<real_prec>(aux_v[i])/static_cast<real_prec>(lkk);
                 aux_v.clear(); aux_v.shrink_to_fit();

                 So.message_screen("Writting the joint distribution P(X,Y)");
                 vector<vector<real_prec> > Vaux;
                 Vaux.resize(new_nbins_x);
                 for(int i=0;i<new_nbins_x;++i)
                   for(int j=0;j<new_nbins_y;++j)
                     Vaux[i].push_back(aux_n[index_2d(i,j,new_nbins_y)]);
                 aux_n.clear(); aux_n.shrink_to_fit();

                 string filek=file;
#if defined(_USE_MASS_KNOTS_) && !defined(_USE_REDSHIFT_MASK_)
                 filek+="_MKbin"+to_string(mk);
#endif

#if defined(_USE_REDSHIFT_MASK_) && !defined(_USE_MASS_KNOTS_)
                 filek+="_zbin"+to_string(iz);
#endif

#if  defined(_USE_REDSHIFT_MASK_) && defined(_USE_MASS_KNOTS_)
                 filek+="_zbin"+to_string(iz)+"_MKbin"+to_string(mk);
#endif

                 this->File.write_to_file(filek+".txt",X_bins,Y_bins,Vaux);

                 mcmc.get_contour_levels(filek+"_contour_levels",Vaux);

                 Vaux.clear(); Vaux.shrink_to_fit();
#ifdef _USE_MASS_KNOTS_
               }
#endif


#ifdef _USE_REDSHIFT_MASK_
           }
#endif



         // ******************************************************************************
         // ******************************************************************************


         So.message_screen("Allocating the values of Y found in the bins of X.");

         // if(this->cwt_used[sua]==I_KNOT)
         //   {

#ifdef _nomas_

#ifdef _USE_REDSHIFT_MASK_
         for(int iz=0;iz<N_REDSHIFT_BINS;++iz)
           {
#endif


#ifdef _USE_MASS_KNOTS_
             for(int mk=0;mk<this->n_sknot_massbin;++mk)
               {
#endif
                 vector<vector<real_prec> >DELTAXY;
                 DELTAXY.resize(this->new_nbins_x);

                 // Get the values of Y in each bin of X
                 for(ULONG i=0;i<this->NGRID ;++i)
                   {
#ifdef _USE_BINARY_MASK_
                     if(static_cast<real_prec>(binary_mask[i])>0)
                       {
#endif

#ifdef _USE_MASS_KNOTS_
                             if(this->cwclass.SKNOT_M_info[i]==mk)
                               {
#endif

#ifdef _USE_REDSHIFT_MASK_
                                 if(redshift_mask_bin[i]==iz)
                                   {
#endif

                                     // This last condition is important, specially when log10(0) is involved.
                                     // For Y, we ask to do statistics with cells that have at least ione halo, the > imposes it
                                     {
                                      if((delta_X[i]>=this->Xmin && delta_X[i]<=this->Xmax) && (delta_Y[i]>=this->Ymin && delta_Y[i]<=this->Ymax))
                                      {
                                        int IX=get_bin(delta_X[i],this->Xmin,this->new_nbins_x,DELTAX,this->bin_accumulate_borders);
                                         DELTAXY[IX].push_back(delta_Y[i]);// IX ES EL BIN DE X: allocate all values of Y in the bins of dm.
                                       }
                                      }
#ifdef _USE_REDSHIFT_MASK_
                                   }
#endif

#ifdef _USE_MASS_KNOTS_
                               }
#endif

#ifdef _USE_REDSHIFT_MASK_
                           }
#endif



                   }


                 So.message_screen("Done");
                 this->mean_Y.clear();
                 this->mean_Y.shrink_to_fit();
                 this->mean_Y.resize(new_nbins_x,0);
                 vector<real_prec> var_Y(new_nbins_x,0);
                 vector<real_prec> skew_Y(new_nbins_x,0);
                 vector<real_prec> kurt_Y(new_nbins_x,0);

                 So.message_screen("Getting the PDF  P(Y|X) (mean, rms, higher moments).");

                 for(int i=0;i<new_nbins_x;++i)// loop sobre bines de dark matter
                   {
                     // Ddeefine this vector to perform some statistics
                     vector<real_prec> vEPSILON;
                     for(int j=0;j<DELTAXY[i].size();++j)
                       vEPSILON.push_back(DELTAXY[i][j]);


//                     for(int j=0;j<vEPSILON.size();++j)
 //                       cout<<vEPSILON[j]<<endl;

                     // Mean of Y in the current X bin
                     mean_Y[i]= vEPSILON.size()>0 ?  get_mean(vEPSILON): 0 ;

                     // RMS in the current X bin
                     var_Y[i]=vEPSILON.size()> num_1 ? sqrt(get_var(vEPSILON)) : 0.  ;

                     // Skewness
                     //                  skew_Y[i]=gsl_stats_skew(&vEPSILON[0],num_1, vEPSILON.size());

                     // Kurtossis
                     //                  kurt_Y[i]=gsl_stats_kurtosis(&vEPSILON[0],num_1, vEPSILON.size());


                     // THIS IS TO WRITE AND DO HISTOGRAMS, BUT IF WE KEEP TRACK OF THE MOMENTS, WE COULD GET THE DISTRIBUTION
                     if (true==this->write_files_for_histograms && vEPSILON.size()>0)
                       {
                         string file_X_bins=file+"_dist_bin"+to_string(i)+".txt";;
#if defined(_USE_MASS_KNOTS_) && !defined(_USE_REDSHIFT_MASK_)
                         file_X_bins=file+"_dist_bin"+to_string(i)+"_MKbin"+to_string(mk)+".txt";
#endif


#if defined(_USE_REDSHIFT_MASK_) && !defined(_USE_MASS_KNOTS_)
                         file_X_bins=file+"_dist_bin"+to_string(i)+"_zbin"+to_string(iz)+".txt";
#endif

#if defined(_USE_REDSHIFT_MASK_) && defined(_USE_MASS_KNOTS_)
                         file_X_bins=file+"_dist_bin"+to_string(i)+"_zbin"+to_string(iz)+"_MKbin"+to_string(mk)+".txt";
#endif

                         this->sal.open(file_X_bins.c_str());
                         So.message_screen("Writting values of Y in bins of X in file ",file_X_bins);
                         for(int ie=0;ie<vEPSILON.size();++ie)sal<<vEPSILON[ie]<<endl;
                         sal.close();
                       }
                     vEPSILON.clear();
                     vEPSILON.shrink_to_fit();
                   }


                 string out2=file+"_mean"+XNAME+"bins.txt";
#if defined(_USE_MASS_KNOTS_) && !defined(_USE_REDSHIFT_MASK_)
                 out2=file+"_mean"+XNAME+"bins_MKbin"+to_string(mk)+".txt";
#endif

#if defined(_USE_REDSHIFT_MASK_) && !defined(_USE_MASS_KNOTS_)
                 out2=file+"_mean"+XNAME+"bins_zbin"+to_string(iz)+".txt";
#endif

#if  defined(_USE_REDSHIFT_MASK_) && defined(_USE_MASS_KNOTS_)
                 out2=file+"_mean"+XNAME+"bins_zbin"+to_string(iz)+"bins_MKbin"+to_string(mk)+".txt";
#endif


                 this->sal.open(out2.c_str());
                 int lnc=0;
                 for(int i=0;i<X_bins.size();++i)
                   if(mean_Y[i]!=0)
                     {
                       lnc++;
                       this->sal<<X_bins[i]<<"\t"<<mean_Y[i]<<"\t"<<var_Y[i]<<endl;//"\t"<<skew_Y[i]<<"\t"<<kurt_Y[i]<<endl;
                     }
                 So.message_screen("Wrote output in file", out2, ". Number of lines =", lnc);
                 this->sal.close();
#ifdef _USE_MASS_KNOTS_
               }
#endif


#ifdef _USE_REDSHIFT_MASK_
           }
#endif



#ifdef _USE_CWC_
       }
#endif

#endif   //edif del nomas

#ifdef _USE_CWC_
}
#endif   //edif del nomas


}


   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************

// *******************************************************************
#ifdef _SEARCH_CLOSEST_DM_IN_NEIGH_CELLS_

struct s_particle_neighbours
{
  vector<ULONG>global_id_candidate;
  vector<ULONG>id_dist_candidate;
  ULONG sort_candidates(){
    ULONG rank=0;
    sort_1d_vectors<ULONG>(this->id_dist_candidate,rank);
    return this->global_id_candidate[rank];
  }

};
#endif

// *******************************************************************

#ifdef MOCK_MODE
   void Bam::collapse_randoms()
   {
     So.enter(__PRETTY_FUNCTION__);
#ifdef _FULL_VERBOSE_
     this->So.message_screen("**Collapsing randoms towards the DM particles**");
     cout<<endl;
#endif

     int NTHREADS = _NTHREADS_;  // omp_get_max_threads();
     omp_set_num_threads(NTHREADS);

     s_params_box_mas box_collps;
     box_collps.min1=this->params._xllc();
     box_collps.min2=this->params._yllc();
     box_collps.min3=this->params._zllc();
     box_collps.Lbox=this->params._Lbox();
     box_collps.Nft=this->Nft_random_collapse;
     box_collps.d1= box_collps.Lbox/static_cast<real_prec>(box_collps.Nft);		/* grid spacing x-direction */
     box_collps.d2= box_collps.d1;
     box_collps.d3= box_collps.d1;
     box_collps.NGRID=(box_collps.Nft*box_collps.Nft*box_collps.Nft);


     ULONG Ntracers=this->tracer.NOBJS;


     vector<int>dm_count(box_collps.NGRID,0);
     vector<real_prec> x_dm_pos;
     vector<real_prec> y_dm_pos;
     vector<real_prec> z_dm_pos;
     vector<real_prec> mass_dm;
     vector<ULONG> dm_id;
     vector<ULONG> dm_id_global;

     vector<real_prec> x_random_pos;
     vector<real_prec> y_random_pos;
     vector<real_prec> z_random_pos;
     vector<real_prec> mass_random;
     vector<ULONG> ran_id;
     vector<ULONG> ran_id_global;

#ifdef _FULL_VERBOSE_
     this->So.message_screen("Separating DM and random:");
#endif

     ULONG N_dms=0;
     ULONG N_rand=0;

     for(ULONG i=0; i< Ntracers; ++i)
       {
         ULONG id=grid_ID(&box_collps, this->tracer.Halo[i].coord1,this->tracer.Halo[i].coord2,this->tracer.Halo[i].coord3);
         if(this->tracer.Halo[i].identity>0)
           {
             x_dm_pos.push_back(this->tracer.Halo[i].coord1);
             y_dm_pos.push_back(this->tracer.Halo[i].coord2);
             z_dm_pos.push_back(this->tracer.Halo[i].coord3);
             mass_dm.push_back(this->tracer.Halo[i].mass);
             dm_id.push_back(id);// For each dm, keep the ID of the cell in which it s located
             dm_id_global.push_back(i);
             dm_count[id]++;     // Count the numnber of DM particles in this cell
             N_dms++;            // Count the numnber of DM particles
           }
         else if(this->tracer.Halo[i].identity<0)
           {
             x_random_pos.push_back(this->tracer.Halo[i].coord1);
             y_random_pos.push_back(this->tracer.Halo[i].coord2);
             z_random_pos.push_back(this->tracer.Halo[i].coord3);
             mass_random.push_back(this->tracer.Halo[i].mass);
             ran_id.push_back(id);// For each random, keep the ID of the cell in which it s located
             ran_id_global.push_back(i);
             N_rand++;            // Count the numnber of random particles
           }
       }
     this->So.DONE();

     int max_count=get_max(dm_count); // get the maximum number of dm particles in one cell.
#ifdef _FULL_VERBOSE_
     So.message_screen("Maximum number of DM particles found in one cell (low resolution) = ", max_count);
#endif


     // Initialize this to NO_NUM, a negative number, for not all elements of the container
     // dm_index_cell will be filled and initializing with 0 makes confusion, for 0 is a valid/used entry
#define NO_NUM -999
     ULONG Mdm_index_cell=box_collps.NGRID*max_count;


     vector<ULONG> dm_index_cell(Mdm_index_cell, NO_NUM);
#ifdef _SEARCH_CLOSEST_DM_IN_NEIGH_CELLS_
     vector<ULONG>dm_particle_id(Mdm_index_cell,0);
#endif

#ifdef _FULL_VERBOSE_
     So.message_screen("Getting ids of dm particles in low-res cells");
#endif
     ULONG idm=0;
     dm_count.clear(); // We need to count again in order to assign index in cells
     dm_count.shrink_to_fit();
     dm_count.resize(box_collps.NGRID,0);
     for(ULONG i=0;i<Ntracers;++i) // Loop over the DM particles
     {
       if(this->tracer.Halo[i].identity>0) // Loop over the DM particles
         {
           ULONG idc=dm_id[i];         // Get the ID of the cell where this DM particle is located
           ULONG dmc=dm_count[idc];    // Retrieve the current number of dm particles (local id) in this IDcell, 0, 1, 2, ..., Ndm(cell) in the cell
           dm_index_cell[index_2d(idc, dmc, max_count)]=idm; // Associate the partial id in 0, 1, 2, ...Ndm-1 of the dm particle in the cell 
#ifdef _SEARCH_CLOSEST_DM_IN_NEIGH_CELLS_
           dm_particle_id[index_2d(idc, dmc, max_count)]=i; // Associate the global i in (0,Ntracers) of the dm particle in the cell 
#endif
           dm_count[idc]++;              // Count the number of dm particle sin this particular lowres cell
           idm++;                        // Count dm particles, to get the partial id
         }
      }
     dm_id.clear();
     dm_id.shrink_to_fit();
     this->So.DONE();

     vector<ULONG>dm_index_closer_tot(N_rand,0);

#ifdef _SEARCH_CLOSEST_DM_IN_NEIGH_CELLS_
     vector<s_particle_neighbours> dm_closest_to_ran(N_rand);
#endif

#ifdef _FULL_VERBOSE_
     So.message_screen("Identifying closest DM particles for randoms IN-cell");
#endif


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0; i< N_rand; ++i)//loop over randoms
       {
         ULONG id_tracer=ran_id[i]; // identify the ID (low-res cell) where this random lives. 
         real_prec xra=x_random_pos[i];
         real_prec yra=y_random_pos[i];
         real_prec zra=z_random_pos[i];

         vector<ULONG>i_r_to_dm_dist; // define array to allocate distance in form of index
         vector<ULONG>n_index_tot;    // define array to allocate the label of the dm particle in the cell

         ULONG N_dm_cell=dm_count[id_tracer];// Get the number of dm particles in this cell

         if(N_dm_cell>0)// if there are dm particles in the cell where this random is located, then
           for(ULONG jc = 0; jc < N_dm_cell; ++jc)//loop over the dm particles in the cell where the current random is 
              {
                ULONG jdm=dm_index_cell[index_2d(id_tracer,jc,max_count)]; //indice entre (0,N_dm-1) que tiene cada particula de dm dentro de la celda id_ran
                if(jdm!=NO_NUM)
                   {
                     real_prec xda=xra-x_dm_pos[jdm];// x_random - x_dm 
                     real_prec yda=yra-y_dm_pos[jdm];
                     real_prec zda=zra-z_dm_pos[jdm];
                     ULONG i_dist=_get_modulo_squared(xda,yda,zda);  //distance² between rand and each dm particle, converted to floor
                     i_r_to_dm_dist.push_back(i_dist);  // allocate the index for the distance
                     n_index_tot.push_back(jdm);        // allocat the id of the dm particle
                   }
              }
         else
           So.message_screen("No dm particles found in the cell corresponding to the random tracer ", i, ". You might want to increase the cell-size");
            // Now sort and get the identification to the closest dm particle to the current i-random
            // Sort with respect to i_r_to_dm_dist; the n_index_tot vector is sorted correspndingly such that
            // its first element is the id of the nearest dm particle to the random i. 
            //The function sort_1dvectors is modified to avoid loopws and returns the ero element of the sorted array,

            // Note that here we allocate jdm in n_index_tot and  dm_particle_id[s_ind] to the vector if structures  dm_closest_to_ran[i]
            // In the loop OUTcell below we allocate instead dm_particle_id[s_ind] directly to n_index_tot in order to sort wrt distances

         ULONG rank=0;
         sort_1d_vectors<ULONG>(i_r_to_dm_dist, rank); // sort the container with distance indices
         ULONG s_ind=n_index_tot[rank];                // get the if dm in (0,.Ndm-1) particle associated to the first element of the sorted distance

#ifdef _SEARCH_CLOSEST_DM_IN_NEIGH_CELLS_
         ULONG s_val=i_r_to_dm_dist[rank];             // get the distance index
         dm_closest_to_ran[i].id_dist_candidate.push_back(s_val); // This is the distance_id of the closest dm particle in the cell. Keep this value to compare with those from the outcell search
         dm_closest_to_ran[i].global_id_candidate.push_back(s_ind);              // keep the global id of this dm candidate

#else
         dm_index_closer_tot[i]=s_ind;// get the identification (within this cell) of the closest dm particle                             
#endif
       }
      So.DONE();


// Now look for in neighbour cells for closest dm particles

#ifdef _SEARCH_CLOSEST_DM_IN_NEIGH_CELLS_

#ifdef _FULL_VERBOSE_
     So.message_screen("Identifying neighbouring cells");
#endif
     vector<s_nearest_cells>nearest_cells(box_collps.NGRID);
     get_neighbour_cells(box_collps.Nft,1,nearest_cells);
     So.DONE();

#ifdef _FULL_VERBOSE_
     So.message_screen("Identifying closest DM particles for randoms OUT-cell");
#endif

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0; i< N_rand; ++i)//loop over randoms
       {
         ULONG id_tracer=ran_id[i]; // identify the ID (low-res cell) where this random lives. 
         real_prec xra=x_random_pos[i];
         real_prec yra=y_random_pos[i];
         real_prec zra=z_random_pos[i];
         vector<ULONG>i_r_to_dm_dist; // define array to allocate distance in form of index
         vector<ULONG>n_index_tot;    // define array to allocate the label of the dm particle in the cell

         int N_neighbour_cells=nearest_cells[id_tracer].close_cell.size(); // Number of neighbour cells oif the current cell
         for (int ic=0; ic < N_neighbour_cells; ++ic)  // Loop over neighbour cells
           {
             
             ULONG id_neighbour=nearest_cells[id_tracer].close_cell[ic]; // ID of the current neighbout cell
             ULONG N_dm_cell=dm_count[id_neighbour];                     // Get the number of dm particles in this cell
          
             if(N_dm_cell > 0)// for this "if" I do not write an else, for if there are no dm in neighbourr cells, we do not care
              {
                for(ULONG jc = 0; jc < N_dm_cell; ++jc)//loop over the dm particles in the neighbour cell id_n
                  {
                    ULONG jdm  = dm_index_cell[index_2d(id_neighbour, jc, max_count)]; //index jdm in(0,N_dm-1) for each dm particle in the cell id_n
                    if(jdm!=NO_NUM)
                     {
                       real_prec xda = xra-x_dm_pos[jdm];
                       real_prec yda = yra-y_dm_pos[jdm];
                       real_prec zda = zra-z_dm_pos[jdm];
                       ULONG i_dist= _get_modulo_squared(xda,yda,zda); //100000 * distance² between rand and each dm particle, converted to floor
                       ULONG id_dm = dm_particle_id[index_2d(id_neighbour, jc, max_count)]; // global index in (0,N_tracer-1) for each dm particle in the cell id_n
                       i_r_to_dm_dist.push_back(i_dist);  // allocate the distance index          
                       n_index_tot.push_back(id_dm);// allocate the global id of the dm particle 
                     }
                  }               
               }
            }


          ULONG rank=0;
          //For all distances in each neighbour cell, sort the distances and get the closest dm
          sort_1d_vectors<ULONG>(i_r_to_dm_dist, rank); 
          dm_closest_to_ran[i].id_dist_candidate.push_back(i_r_to_dm_dist[rank]);       // keep this value to compare with those formn the incell search
          dm_closest_to_ran[i].global_id_candidate.push_back(n_index_tot[rank]);     // keep the global id of this dm candidate
       }
      So.DONE();

// Here find the minimum between the two sets in and out cell
#ifdef _FULL_VERBOSE_
      So.message_screen("Identifying id of closest DM particle");
#endif

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0; i< N_rand; ++i)//loop over randoms
       dm_index_closer_tot[i]=dm_closest_to_ran[i].sort_candidates();// get the global identification                
 So.DONE();
 nearest_cells.clear();nearest_cells.shrink_to_fit();
 dm_closest_to_ran.clear();dm_closest_to_ran.shrink_to_fit();

#endif // end of SEARCH_CLOSEST_DM_IN_NEIGH_CELLS_


#ifdef _FULL_VERBOSE_
     So.message_screen("Freeing memmory in line", __LINE__);
#endif
   dm_count.clear();
   dm_count.shrink_to_fit();
   ran_id.clear();ran_id.shrink_to_fit();
   So.DONE();


#ifdef _FULL_VERBOSE_
     this->So.message_screen("Collapsing randoms using fraction of distance to closest dm particle = ", this->Distance_fraction);
#endif

#if defined _COLLAPSE_RANDOMS_USING_EXCLUSION_ || defined (_APPLY_GLOBAL_EXCLUSION_)
        // Estimate mean density of dark matter halos as rho_mean_background times Delta_vir
     real_prec density=this->Cosmo.mean_matter_density(this->params._redshift(), (void *)&this->s_cosmo_pars)*this->s_cosmo_info.Delta_Vir;
#endif


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG counter=0; counter< N_rand; ++counter)
      {
         //total index (i.e, from 0 to N_dms) of the closest dm particle to the ramdom i. Useful to retrieve coordinates
         // When closesd tm in only searched in the cell (incell), this gives the id of the closest dm tracer in the cell
         // Otherwise, this gives the global id of the closest dm particle
         ULONG index_dm_closer_a = dm_index_closer_tot[counter]; 

         // Read coordinates of random particles
         real_prec new_x = x_random_pos[counter];
         real_prec new_y = y_random_pos[counter];
         real_prec new_z = z_random_pos[counter];

#ifdef _SEARCH_CLOSEST_DM_IN_NEIGH_CELLS_
         real_prec xdm = this->tracer.Halo[index_dm_closer_a].coord1;
         real_prec ydm = this->tracer.Halo[index_dm_closer_a].coord2;
         real_prec zdm = this->tracer.Halo[index_dm_closer_a].coord3;
#else
         real_prec xdm = x_dm_pos[index_dm_closer_a]; 
         real_prec ydm = y_dm_pos[index_dm_closer_a];
         real_prec zdm = z_dm_pos[index_dm_closer_a];
#endif

         // redefine ran coords to the ref sistem of its closest dm particle:
         new_x -= xdm;
         new_y -= ydm;
         new_z -= zdm;

#ifdef _COLLAPSE_RANDOMS_USING_EXCLUSION_
         real_prec dm_mass = mass_dm[index_dm_closer_a];              // Mass of dm particle
         real_prec radius_dm = pow(3.0*dm_mass/(4.*M_PI*density), 1./3.);
         real_prec ran_mass = mass_random[counter];                                // Mass of random
         real_prec radius_ran = pow(3.0*ran_mass/(4.*M_PI*density), 1./3.);             // get an estimate of the radius of the
         real_prec RR=radius_dm+radius_ran;              // sum of the radii of the two halos
#endif

         // get the distance bewteen randoms and their closest DM particle:
         real_prec dist_random_to_dm = _get_modulo(new_x,new_y,new_z);
         //get the angular coordinates:
         real_prec theta = acos(new_z/dist_random_to_dm);
         real_prec phi   = atan2(new_y,new_x);
         real_prec new_sep = this->Distance_fraction*dist_random_to_dm;

             // The separation between halos must be greater or equal than the sum of their radii
#ifdef _COLLAPSE_RANDOMS_USING_EXCLUSION_
         real_prec sep_with_exclusion=max(new_sep,RR);
#else
         real_prec sep_with_exclusion=new_sep;
#endif
             // transfor to cartesian given a reduced distance and return ot origin of coords:
         this->tracer.Halo[ran_id_global[counter]].coord1=sep_with_exclusion*sin(theta)*cos(phi)+xdm;
         this->tracer.Halo[ran_id_global[counter]].coord2=sep_with_exclusion*sin(theta)*sin(phi)+ydm;
         this->tracer.Halo[ran_id_global[counter]].coord3=sep_with_exclusion*cos(theta)+zdm;
       }
     this->So.DONE();


#ifdef _APPLY_GLOBAL_EXCLUSION_


#ifdef _FULL_VERBOSE_
     this->So.message_screen("Applying random shift to objects above ", PROP_THRESHOLD_MULTI_SCALE_4);
#endif



#ifdef _USE_OMP_

     vector<ULONG>vseeds(NTHREADS,0);
     for(int i=0;i<vseeds.size();++i)
         vseeds[i]=35+static_cast<ULONG>(i)*565;

      const gsl_rng_type *  rng_t;
      gsl_rng * gBaseRand;
      int jthread=0;
      gsl_rng_env_setup();

#pragma omp parallel private(jthread, rng_t,gBaseRand)
  {

      jthread=omp_get_thread_num();
      gsl_rng_default_seed=vseeds[jthread];

      rng_t = gsl_rng_mt19937;//_default;
      gBaseRand = gsl_rng_alloc (rng_t);

#pragma omp for
#endif
     for(ULONG i=0; i< this->tracer.NOBJS; ++i)
      {
         if(this->tracer.Halo[i].vmax>PROP_THRESHOLD_MULTI_SCALE_4)
         {
          real_prec x = this->tracer.Halo[i].coord1;
          real_prec y = this->tracer.Halo[i].coord2;
          real_prec z = this->tracer.Halo[i].coord3;    
          real_prec mass = this->tracer.Halo[i].mass;              // Mass of dm particle
          real_prec radius = 6.0*pow(3.0*mass/(4.*M_PI*density), 1./3.);
          this->tracer.Halo[i].coord1= x + gsl_ran_gaussian(gBaseRand, radius);
          this->tracer.Halo[i].coord2= y + gsl_ran_gaussian(gBaseRand, radius);
          this->tracer.Halo[i].coord3= z + gsl_ran_gaussian(gBaseRand, radius);
         }
      }

#ifdef _USE_OMP_
}
#endif
     this->So.DONE();


#endif




   }


   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   void Bam::collapse_randoms_isolated()
   {
       So.enter(__PRETTY_FUNCTION__);
     cout<<endl;

#ifdef _FULL_VERBOSE_
     this->So.message_screen("**Collapsing randoms towards the DM particles**");
     cout<<endl;
#endif

     int NTHREADS = _NTHREADS_;  // omp_get_max_threads();
     omp_set_num_threads(NTHREADS);


     vector<real_prec> prop;
     this->patchy.fnameTRACERCAT=this->Output_directory+"CAT_realization"+to_string(1)+"_"+this->patchy.stradd+string(".txt");

     s_params_box_mas box;
     box.min1=this->params._xllc();
     box.min2=this->params._yllc();
     box.min3=this->params._zllc();
     box.Lbox=this->params._Lbox();
     box.Nft=this->Nft_random_collapse;
     box.d1= box.Lbox/static_cast<real_prec>(box.Nft);		/* grid spacing x-direction */
     box.d2= box.d1;
     box.d3= box.d1;
     box.NGRID=(box.Nft*box.Nft*box.Nft);


     ULONG Ntracers=File.read_file(this->patchy.fnameTRACERCAT,prop, omp_get_max_threads());
     int NCOLS=(static_cast<ULONG>(prop.size()/Ntracers));

     vector<real_prec> x_dm_pos;
     vector<real_prec> y_dm_pos;
     vector<real_prec> z_dm_pos;

     vector<real_prec> x_random_pos;
     vector<real_prec> y_random_pos;
     vector<real_prec> z_random_pos;
     vector<ULONG> ran_id;
     vector<ULONG> dm_id;
     vector<int>dm_count(box.NGRID,0);
#ifdef _FULL_VERBOSE_
     this->So.message_screen("Separating DM and random:");
#endif

     ULONG N_dms=0;
     ULONG N_rand=0;

     // struct s_cell_info{
     //   vector<real_prec> id_within_cell;
     // };

     // vector<s_cell_info> dm_inf(box.NGRID);

     int identificator=7;  // column where we put the flag identifying if the particle is dm or random

     for(ULONG i=0; i< Ntracers; ++i)
       {
         ULONG id=grid_ID(&box, prop[0+i*NCOLS],prop[1+i*NCOLS],prop[2+i*NCOLS]);
         if(prop[identificator+i*NCOLS]>0)
           {
             x_dm_pos.push_back(prop[0+i*NCOLS]);
             y_dm_pos.push_back(prop[1+i*NCOLS]);
             z_dm_pos.push_back(prop[2+i*NCOLS]);
             dm_id.push_back(id);
             dm_count[id]++;
             N_dms++;
           }
         else
           {
             x_random_pos.push_back(prop[0+i*NCOLS]);
             y_random_pos.push_back(prop[1+i*NCOLS]);
             z_random_pos.push_back(prop[2+i*NCOLS]);
             ran_id.push_back(id);
             N_rand++;
           }
       }
     this->So.DONE();
#ifdef _FULL_VERBOSE_

     So.message_screen("Number of tracers associated to DM particles =", N_dms);
     So.message_screen("(", 100.0*static_cast<real_prec>(N_dms)/(static_cast<real_prec>(N_dms)+static_cast<real_prec>(N_rand)), "%)");

     So.message_screen("Number of tracers associated to random particles =", N_rand);
     So.message_screen("(", 100.0*static_cast<real_prec>(N_rand)/(static_cast<real_prec>(N_dms)+static_cast<real_prec>(N_rand)), "%)");
#endif
     int max_count=get_max(dm_count); // get the maximum number of dm particles in one cell.
#ifdef _FULL_VERBOSE_
     So.message_screen("Maximum number of DM particles found in one cell (low resolution) = ", max_count);
#endif
     // Initialize this to NO_NUM, a negative number, for not all elements of the container
     // dm_index_cell will be filled and initializing with 0 makes confusion, for 0 is a valid/used entry



#define NO_NUM -999
     vector<vector<int> > dm_index_cell(box.NGRID, vector<int>(max_count,NO_NUM));

     dm_count.clear();
     dm_count.resize(box.NGRID,0);
#ifdef _FULL_VERBOSE_
     So.message_screen("Getting ids in cells");
#endif

#pragma omp parallel for
     for(ULONG idm=0;idm<N_dms;++idm)
       {
         ULONG id=dm_id[idm];
         dm_index_cell[id][dm_count[id]]=idm;
#pragma atomic update
         dm_count[id]++;
       }
     this->So.DONE();
     dm_id.clear();dm_id.shrink_to_fit();

     vector<int>dm_index_closer_tot(N_rand,0);
#ifdef _FULL_VERBOSE_
     So.message_screen("Identifying closest DM particles for randoms");
#endif

#pragma omp parallel for
     for(ULONG i=0; i< N_rand; ++i)//loop over randoms
       {
         int count =0;
         int id_ran=ran_id[i];
         int N_dm_cell=dm_count[id_ran];
         vector<ULONG>i_r_to_dm_dist;
         vector<ULONG>n_index_tot;
         if(N_dm_cell>0)// if there are dm particles in the cell where this random is located, then
           {
             for(int jc = 0; jc < N_dm_cell; ++jc)//loop over the dm in that cell where the random is
               {
                 int jdm=dm_index_cell[id_ran][jc]; //indice entre (0,Ndm-1) que tiene cada particula de dm dentro de la celda id_ran
                 if(jdm!=NO_NUM)
                   {
                     real_prec xdd=x_random_pos[i]-x_dm_pos[jdm];
                     real_prec ydd=y_random_pos[i]-y_dm_pos[jdm];
                     real_prec zdd=z_random_pos[i]-z_dm_pos[jdm];
                     real_prec dist_dm_r=_get_modulo(xdd,ydd,zdd) ; //distance between rand and teach central
                     int i_dist=static_cast<int>(floor(100.0*dist_dm_r)); // done in order to sort vectors below. 100 is arbitrary
                     i_r_to_dm_dist.push_back(i_dist);
                     n_index_tot.push_back(jdm);
                     count++;
                   }
               }
             sort_1dvectors(i_r_to_dm_dist, n_index_tot);
             dm_index_closer_tot[i]=n_index_tot[0]; // el primer elemento (0) es el menor
           }
         else
           So.message_screen("No dm particles found in the cell correspoding to random ", i, ". You might want to increase the cell-size");
       }
     So.DONE();
     dm_count.clear();
     dm_count.shrink_to_fit();
     ran_id.clear();ran_id.shrink_to_fit();

#ifdef _FULL_VERBOSE_
     this->So.message_screen("Collapsing randoms now:");
     this->So.message_screen("Using fraction of distance to closest dm particle = ", this->Distance_fraction);
#endif

     int counter=0;
     for(ULONG i=0; i<Ntracers; ++i)
       {
         if(prop[identificator+i*NCOLS]<0)  //randoms
           {
             int index_dm_closer_a=dm_index_closer_tot[counter]; //total index (i.e, from 0 to N_dms) of the closest dm particle to the ramdom i. Usefil to retrieve coordinates
             real_prec xdm=x_dm_pos[index_dm_closer_a];
             real_prec ydm=y_dm_pos[index_dm_closer_a];
             real_prec zdm=z_dm_pos[index_dm_closer_a];
             // redefine ran coords to the ref sistem of its closest dm particle:
             real_prec new_x=x_random_pos[counter]-xdm;
             real_prec new_y=y_random_pos[counter]-ydm;
             real_prec new_z=z_random_pos[counter]-zdm;
             // get the distance:
             real_prec dist_random_to_dm=_get_modulo(new_x,new_y,new_z);
             //get the angular coordinates:
             real_prec theta=acos(new_z/dist_random_to_dm);
             real_prec phi=atan2(new_y,new_x);
             real_prec new_distance = this->Distance_fraction*dist_random_to_dm;
             // transfor to cartesian given a reduced distance and return ot origin of coords:
             prop[0+i*NCOLS]=new_distance*sin(theta)*cos(phi)+xdm;
             prop[1+i*NCOLS]=new_distance*sin(theta)*sin(phi)+ydm;
             prop[2+i*NCOLS]=new_distance*cos(theta)+zdm;
             counter++;
           }
       }
     this->So.DONE();
     x_random_pos.clear();
     x_random_pos.shrink_to_fit();
     y_random_pos.clear();
     y_random_pos.shrink_to_fit();
     z_random_pos.clear();
     z_random_pos.shrink_to_fit();
     x_dm_pos.clear();
     x_dm_pos.shrink_to_fit();
     y_dm_pos.clear();
     y_dm_pos.shrink_to_fit();
     z_dm_pos.clear();
     z_dm_pos.shrink_to_fit();

#ifdef _FULL_VERBOSE_
     this->So.message_screen("Writing new catalog");
#endif
     ofstream trout;
     trout.precision(8);
     trout.setf(ios::showpoint);
     trout.setf(ios::scientific);
     string newcat=this->patchy.fnameTRACERCAT;
     this->patchy.fnameTRACERCAT=newcat;
     trout.open(newcat);
     for(ULONG i = 0; i<Ntracers; ++i)
       {
         for(int j = 0; j<NCOLS; ++j)
           trout<<prop[j+i*NCOLS]<<"\t";
         trout<<endl;
       }
     trout.close();
     this->So.DONE();
     prop.clear();
     prop.shrink_to_fit();
   }
#endif


   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************


 void Bam::bamrunner()
  {
    cout<<endl;
    time_t time_bam;
    time(&time_bam);

#ifdef _GET_BAM_REALIZATIONS_
    So.message_screen("Bulding BAM catalogs");
    So.message_screen("Realization ", this->params._IC_index());
#endif
    this->So.initial_time=time_bam;
    So.enter(__PRETTY_FUNCTION__);
    int NTHREADS=_NTHREADS_;   // omp_get_max_threads();
    omp_set_num_threads(NTHREADS);

       this->warnings();
#ifdef MOCK_MODE
#ifdef _GET_BAM_REALIZATIONS_
#ifdef _FULL_VERBOSE_
       cout<<endl;
     So.message_screen("GENERATING NEW TRACER DENSITY FIELD FROM BIAS AND KERNEL");
     cout<<endl;
#endif
#else
#ifdef _FULL_VERBOSE_
     cout<<endl;
     So.message_screen("CALIBRATING BIAS AND KERNEL FROM REFERNECE SIMULATION");
     cout<<endl;
#endif
#endif
#else
       So.message_screen("Statistics of halo bias");
#endif
     // **************************************************************************************************
     // Initialize cosmological functions using the input redshift
     this->tstruct=0;
     this->im=0; //mas

     // **************************************************************************************************
     this->get_cosmo();
     // Get the cosmological quantities derived from input cosmological parameters.
     // **************************************************************************************************
     this->tracer_ref.type_of_object="TRACER_REF";
     this->tracer.type_of_object="TRACER_MOCK";

     this->tracer.set_params_catalog(this->params);
     // **************************************************************************************************
     // **************************************************************************************************

     // Define strings for input DF
     string file_X, file_X_ref_pdf, file_Y,file_Y_HR;
     string file_Vx, file_Vy, file_Vz;



     // **************************************************************************************************
     // **************************************************************************************************
     // IF DEFINED, READ INPUT FILE WITH INFO FROM THE TRACER



     // If READ_REF_CAT is undef or GET_BAM_RALIZATIONS is def , the code will read these fields and analyze then
#ifdef _READ_REF_CATALOG_
     string file_density_field_tracer=this->Output_directory+"TR_DENS_FIELD";
     string file_Y_mass=this->Output_directory+"TR_MASS_DENS_FIELD";
     string file_Y_sat_frac=this->Output_directory+"TR_SAT_FRACTION_FIELD";
#else
     string file_density_field_tracer=this->Output_directory+"TR_DENS_FIELD.dat";
     string file_Y_mass=this->Output_directory+"TR_MASS_DENS_FIELD.dat";
     string file_Y_sat_frac=this->Output_directory+"TR_SAT_FRACTION_FIELD.dat";
#endif

#ifdef _DO_BAM_CALIBRATION_

#ifdef _READ_REF_CATALOG_
     // Here we shall read the reference and produce two density fields,
     // one with number counts and other with mass weighted number counts
     // The one with number counts will be written in the file pointed to
     // by the parameter file_Y, such that it can be read later below

     this->So.message_screen("Reading the Reference Catalog of tracers");
#ifdef _SET_CAT_WITH_MASS_CUT_
     this->tracer.read_catalog(this->params._file_catalogue(),pow(10,params._LOGMASSmin())*params._MASS_units());
#elif defined (_SET_CAT_WITH_VMAX_CUT_)
     this->tracer.read_catalog(this->params._file_catalogue(),params._VMAXmin());
#endif

     this->tracer.get_density_field_grid(_COUNTS_, file_density_field_tracer);


#ifdef _USE_MASS_FIELD_
     this->tracer.get_density_field_grid(_MASS_, file_Y_mass);
     string fname_mass_function_Y = this->Output_directory+"tracer_ref_abundance.txt";
     this->tracer.get_property_function(fname_mass_function_Y);
#endif


#endif

#endif

#ifdef _GET_INTERPOLATED_FIELDS_FROM_BIN_FILES_
     this->tracer.read_catalog_bin(this->NGRID,"A","A","A");
     exit(0);
#endif
     // **************************************************************************************************
     // **************************************************************************************************

     // HERE WE CAN USE GET_NW_DM() ALSO


#if defined(_USE_PATCHY_) || defined (_GET_BAM_CAT_)

     const gsl_rng_type *rng_t;
#ifdef _USE_OMP_
     gsl_rng **gBaseRand;
#else
     gsl_rng *gBaseRand;
#endif

     gsl_rng_env_setup();
     rng_t = gsl_rng_mt19937;// gsl_rng_default;
     int nt=omp_get_max_threads();
     gBaseRand = (gsl_rng **) malloc(nt * sizeof(gsl_rng *));

#ifdef _USE_OMP_
#pragma omp parallel for num_threads(nt)
#endif
     for(int i=0;i<nt;i++)
       {
         gBaseRand[i] = gsl_rng_alloc(rng_t);
         gsl_rng_set(gBaseRand[i],this->seed);
       }


     // This is used for the calibration *and* the construction of mocks
     this->patchy.s_cosmo_info=this->s_cosmo_info;
     //share the s_cosmo_info with Patchy

     this->patchy.set_params_patchy(this->params);
     //Patchy reads its parameters from params class

     this->patchy.set_fnames();
     // set the output name files for patchyt

#endif    // end use patchy and get_bam_cat line 8905



     // If we need to do the calibration,
#ifdef _DO_BAM_CALIBRATION_

     //  we might want to use patchy either to creat the DM, or assign position to partices in mock, or both
#if defined(_USE_PATCHY_) || defined (_GET_BAM_CAT_)

     time_t start_patchy;
     time(&start_patchy);
#endif // end use patchy and get_bam_cat line 8942



#ifdef _USE_PATCHY_


#ifdef _USE_OMP_
     this->patchy.get_dm_field(gBaseRand);
#else
     this->patchy.get_dm_field(gBaseRand);  // Run Patchy!
#endif

     this->So.message_screen("Patchy has created DMDF using ALPT");
     So.message_time2(start_patchy);
     cout<<endl;


#endif  // end use patchy line 8950


#if defined (_USE_PATCHY_) & !defined(_READ_VELOCITIES_)
     // Read the file names generated in Patchy
     file_X =this->patchy.fnameDM+".dat";
     file_Vx=this->patchy.fnameVX+".dat";
     file_Vy=this->patchy.fnameVY+".dat";
     file_Vz=this->patchy.fnameVZ+".dat";
#endif


#ifdef _FULL_VERBOSE_
     this->So.message_screen("BAM:");
#endif


#ifdef _READ_VELOCITIES_
     file_Vx=this->Input_Directory_X+this->Name_VelFieldx_X;
     file_Vy=this->Input_Directory_X+this->Name_VelFieldy_X;
     file_Vz=this->Input_Directory_X+this->Name_VelFieldz_X;
#endif



#endif //ifdef do_bam_calibration


#ifdef _USE_PATCHY_
     file_X=this->patchy.fnameDM+".dat";
#else
     file_X=this->Input_Directory_X+this->Name_Catalog_X;
#endif

     // Input file coptaining the reference tracer
#ifndef _READ_REF_CATALOG_
     file_Y_HR=this->Input_Directory_Y+this->Name_Catalog_Y_HR;
     file_Y=this->Input_Directory_Y+this->Name_Catalog_Y;
#else
     file_Y = file_density_field_tracer+".dat";
     file_Y_HR = file_density_field_tracer+".dat"; //in thie case we pass the same file
     file_Y_mass += ".dat";
     file_Y_sat_frac += ".dat";
#endif


     this->step=0;


     // **************************************************************************************************
#ifdef _DO_BAM_CALIBRATION_
#ifdef _USE_VELOCITIES_
     this->read_bam_files(file_X, file_Y, file_Y_HR, file_Y_mass, file_Y_sat_frac, file_X_ref_pdf, file_Vx, file_Vy, file_Vz);
#endif // end ifdef _USE_VELOCITIES_
#endif

#ifdef _ONLY_PATCHY_  // This region is devoted to mesure the power of the DM generated in Patchy when only Patchy is used.
         this->delta_X_ini.resize(this->NGRID,0);
         this->File.read_array(file_X,delta_X_ini);
         this->get_power_spectrum("DM_REF");
         this->delta_X_ini.clear(); delta_X.shrink_to_fit();
         exit(0);
#endif

#ifdef _FULL_VERBOSE_
         So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Be aware of this ifdef, I sometimes *comment it* when I need to *read* a reference");
#endif

#if defined (_DO_BAM_CALIBRATION_) || defined (BIAS_MODE)
#ifndef _test_mass_assign_
         {
#ifndef _USE_VELOCITIES_
           this->read_bam_files(file_X, file_Y, file_Y_HR, file_Y_mass, file_Y_sat_frac, file_X_ref_pdf);
           // Read all input files (density fields in binary )
#endif
         }
#endif
#endif


         // **************************************************************************
         this->set_Fourier_vectors();          // Intialize arrays for power spectrum of input fields

#if defined (_DO_BAM_CALIBRATION_)  || defined (BIAS_MODE)
         this->analyze_input(); // analyze the input references, requested also to create mock for limits
#endif
         // Get numbers from the input density fields and their power spectrum.
         // In this function, if requested from parameter file, an iterative process is performed in order to make the approx method
         // DM field match the reference DM power (if N_iterations_dm >0)



         // *******************************************************************************************************
         // *******************************************************************************************************
         // load parameters for the cwclass
#if defined (_USE_CWC_) || defined (_USE_MASS_KNOTS_) || defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_INVARIANT_TIDAL_FIELD_IV_) || defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_ELLIPTICITY_)|| defined (_USE_PROLATNESS_) || defined (_USE_INVARIANT_SHEAR_VFIELD_I_) || defined (_USE_INVARIANT_SHEAR_VFIELD_II_) || defined (_USE_INVARIANT_SHEAR_VFIELD_III_) || defined (_USE_NABLA2DELTA_) || defined (_USE_S2DELTA_) || defined (_USE_S3_)  || defined (_USE_S2_) ||  defined (_USE_DELTA3_) || defined (_USE_PWEB_)
         this->cwclass.set_params_cwclass(this->params);
         this->cwclass.s_cosmo_pars=this->s_cosmo_pars;
#endif


         // *******************************************************************************************************
         //  Clean and initialize arrays for power spectrum  and kernels
         this->set_Fourier_vectors();

         // *******************************************************************************************************
         // *******************************************************************************************************

#ifdef BIAS_MODE
#ifdef _SEVERAL_REAL_BIAS_
         this->get_pdf(id);
#elif !defined _SEVERAL_REAL_BIAS_
         this->get_pdf();
#endif // end _SEVERAL_REAL_
#endif // end ifdef BIAS


#ifdef _TEST_THRESHOLDS_RESIDUALS_

//         this->file_residuals=this->Output_directory+"Resuduals_thresholds_original_vel.txt";
         this->file_residuals=this->Output_directory+"chis_thresholds_original_vel.txt";
         this->output_res.open(this->file_residuals.c_str());

         for(int ilt=0; ilt<100;ilt++)
           for(int ilv=0; ilv<100;ilv++)
             {
               cout<<endl;
               this->lambdath_v=-1.0+2.0*(static_cast<real_prec>(ilv))/static_cast<real_prec>(100.0);
               this->lambdath  =(1.0)*(static_cast<real_prec>(ilt))/static_cast<real_prec>(100.0);
               this->cwclass.lambdath=this->lambdath;
               this->cwclass.lambdath_v=this->lambdath_v;
#endif


               // DO the V-classification from the velocity field.
               // This is done once, since the Vel field won't change during the iterative process.
               // Otherwise this should have been done within the iterations, in the get_BAM_DM method.
#if defined (_USE_CWC_V_) || defined (_USE_INVARIANT_SHEAR_VFIELD_I_)  || defined (_USE_INVARIANT_SHEAR_VFIELD_II_) || defined (_USE_INVARIANT_SHEAR_VFIELD_III_)
               this->cwclass.do_CWC_V(this->Velx_X, this->Vely_X,this->Velz_X);


#ifdef _USE_VEL_KNOTS_V_
               this->cwclass.get_SigmaVel_collapsing_regions(this->delta_X_ini,  this->Velx_X, this->Vely_X,this->Velz_X, static_cast<real_prec>(this->N_objects_X)/static_cast<real_prec>(this->NGRID));
#endif



#endif


               // *******************************************************************************************************
               // *******************************************************************************************************
               // ****************************ITERATIVE BAM AND MOCK GENERATION******************************************
               // *******************************************************************************************************
#ifdef MOCK_MODE
               // Get the total number of steps adding those of the calibration of the kernal and the number of
               // DM realizations used to get the same number of mock density fields

               int N_steps=this->N_iterations_Kernel+(this->N_dm_realizations-this->N_dm_initial+1);

               int init=this->iteration_ini;

#ifdef _TEST_THRESHOLDS_RESIDUALS_
               init=0;
               this->N_iterations_Kernel=0;
#endif




               // If we read the Kernel, we start directly creating the mocks

               // **************************************************************************************************
               this->use_iteration_ini=false;

               if(this->iteration_ini>0)
                 this->use_iteration_ini=true; //used ask the code to read kernel and bias from an iteration, in case the code breaks down during calibration


#ifdef _GET_BAM_REALIZATIONS_
               init = this->N_iterations_Kernel+1;
               this->get_min_max_X_Y();

               // in get bias we read the bias and the kernel ,only when we want to get the realizations
               this->get_BIAS(this->Name_Property_Y);
#endif

               // i=0, zero order approach
               // from i=1, to i=N_iterations_Kernel, we calibrate the Kernel.
               // From i=N_iterations_Kernel+1 to i=N_steps, we create mocks based on independent realization of approx DM fields

#ifdef _FULL_VERBOSE_
               time_t start_end;
               time_t start_aux;
               time(&start_end);
               start_aux=start_end;
#endif

#ifdef _DO_BAM_CALIBRATION_

               // Loop over the number if iterations demanded to perform the calibration of the Kernel and the Bias


#ifdef _USE_GNUPLOT_
              this->gp_kernel<<"set size 1.0, 1.0\n";
              this->gp_kernel<<"set origin 0.0, 0.0\n";
              this->gp_kernel<<"set style function lines\n";
#endif


#ifdef _DISPLACEMENTS_
              this->patchy.Displacement.resize(this->NGRID,0);
#endif

               for(int i=init; i<=this->N_iterations_Kernel ;++i)
                 {

#ifdef _USE_GNUPLOT_
              this->gp_kernel<<"set multiplot\n";
#endif



#ifdef _FULL_VERBOSE_
                   time_t start;
                   time(&start);
                   cout<<endl;
                   cout<<CYAN<<"**************************************************************************"<<RESET<<endl;
                   cout<<CYAN<<"**************************************************************************"<<RESET<<endl;
#endif
                   this->step=i;
#if defined _USE_CWC_ || defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_INVARIANT_TIDAL_FIELD_III_)
                   this->cwclass.step=i;
#endif
                   // Calculations performed during the iterative process.

#ifdef _FULL_VERBOSE_
                   if(i==0)
                     So.message_screen("BAM raw mapping" , i, start, start_aux);
                   else
                     So.message_screen("Iteration ", i, start,start_aux);
                   cout<<CYAN<<"**************************************************************************"<<RESET<<endl;
#endif
                   // ****************************************************************************

#ifdef _DISPLACEMENTS_
                   this->GetKernel(true,1);// Get the kernel from the ratio of P_disp_new/P_disp_ref
                   this->Konvolve(this->delta_X,this->delta_X);
                   this->patchy.get_displacement(gBaseRand,this->delta_X);// delta_X acts as IC. This function generates the container vector<real_prec>patchy.displacement

#else
                   this->get_BAM_DM();
#endif
                   // step i) Do the CWC if requested
                   // Get the ratio T, update kernel K, convolve K with original DM field.
                   // The kernel is computed with the Power aspectrum or the mass weighted power spectrum, according to the preproc def _USE_MASS_WEIGHTED_KERNEL_

                   // ****************************************************************************
                   this->get_BIAS(this->Name_Property_Y);
                   // Step ii)
                   // Compute the halo bias from Numnber counts reference and DM from step i)
                   // ****************************************************************************
          		    if(this->Name_Property_Y==_COUNTS_)
		               this->get_mock_grid(_COUNTS_);
          		   else
		               this->get_mock_grid(_DENSITY_);

		   // Step iii)
                   // Generate the Halo number density field by sampling the DM from step i) using the information of the bias from step ii)
                   // argument false indicates that the DM used in the one of the reference (either original Nbody or approximated)
                   // Also gets the power spectrum of the mock


                   // Since the info of the mass distribution and the sat fraction is not used for the clustering analysis to calibrate the bias,
                   // we can compute them in the last step of the iteration, with the DM already transformed with the Bam kernel.

                   // ****************************************************************************
#ifdef _FULL_VERBOSE_
                   start_aux=start;
#endif
                   } // closes iteration

#endif

#ifdef _USE_GNUPLOT_
                   this->gp_kernel<<"unset multiplot\n";
#endif



#ifdef _TEST_THRESHOLDS_RESIDUALS_
               real_prec residuals=0;
               int ncounts=0;

               /*
#pragma omp parallel for reduction(+:residuals,ncounts)
               for(int i=0;i<this->Power_REF.size();++i)
                if(this->Power_NEW[i]>0)
                 {
                   ncounts++;
                  residuals+=fabs(this->Power_REF[i]/this->Power_NEW[i]-1.0);
                 }
               residuals/=static_cast<real_prec>(ncounts);
               So.message_screen("Residuals at this iteration (%) =",100.0*residuals);
               this->output_res<<this->lambdath<<"\t"<<this->lambdath_v<<"\t"<<100.0*residuals<<endl;
             }
               So.message_screen("Residuals at this iteration (%) =",100*residuals);
               this->output_res<<this->lambdath<<"\t"<<this->lambdath_v<<"\t"<<100.0*residuals<<endl;

 */





               real_prec deltak=this->kvec[1]-this->kvec[0];
               for(int i=0;i<this->Power_REF.size();++i)
                if(this->Power_NEW[i]>0)
                 {
                   ncounts++;
                   real_prec sigma_squared = (1./(2.*M_PI*deltak*pow(this->kvec[i],2)*pow(this->Lbox,3))) *(this->Power_NEW[i]*this->Power_NEW[i]); // ESTE PARA EL KERNEL DE LA DM
                   residuals+=pow(this->Power_REF[i]-this->Power_NEW[i],2)/(sigma_squared);
                 }
               residuals/=static_cast<real_prec>(ncounts);

               /* not rady yet:

               if(this->cwclass_ref.volume_knots>0)
                 residuals+=pow(this->cwclass.volume_knots-this->cwclass_ref.volume_knots,2)/sqrt(this->cwclass_ref.volume_knots);

               if(this->cwclass_ref.volume_filaments>0)
                 residuals+=pow(this->cwclass.volume_filaments-this->cwclass_ref.volume_filaments,2)/sqrt(this->cwclass_ref.volume_filaments);

               if(this->cwclass_ref.volume_sheets>0)
                  residuals+=pow(this->cwclass.volume_sheets-this->cwclass_ref.volume_sheets,2)/sqrt(this->cwclass_ref.volume_sheets);

               if(this->cwclass_ref.volume_voids>0)
                 residuals+=pow(this->cwclass.volume_voids-this->cwclass_ref.volume_voids,2)/sqrt(this->cwclass_ref.volume_voids);
                */

               So.message_screen("chi² =",residuals);
               this->output_res<<this->lambdath<<"\t"<<this->lambdath_v<<"\t"<<residuals<<endl;

           }



         this->output_res.close();
#endif


#if defined(_SEVERAL_REAL_BIAS_) || defined(_SEVERAL_REAL_CAL_)
       }
#endif






#ifdef _GET_BAM_REALIZATIONS_
#ifdef _FULL_VERBOSE_
     time(&start_end);
     start_aux=start_end;
#endif
     // Sample the bias measured from the reference into other realzations of the Approximated density field.

     // Loop over the new DM fields

     int i=init;

#ifndef _ONLY_POST_PROC_
     for(i=init; i<=N_steps ;++i)
       {
#endif
#ifdef _FULL_VERBOSE_

             time_t start;
         time(&start);
#endif
         this->step=init;

         int i_realization=i-this->N_iterations_Kernel+this->N_dm_initial-1;

         this->step=i;

#ifndef _ONLY_POST_PROC_
#ifdef _FULL_VERBOSE_
         So.message_screen("Creating mock");
#endif
#else
#ifdef _FULL_VERBOSE_
         So.message_screen("Assigning properties to mock");
#endif
#endif

         // Get the new dm Df from ALPT, convolve with the kernel and  and get its properties ({theta}). This will be used for the conditional assignment of properties
         this->get_new_DM_field();



#ifndef _ONLY_POST_PROC_   // These Lines are to be deprecated
   if(this->Name_Property_Y==_COUNTS_)
     this->get_mock_grid(_COUNTS_);
   else
     this->get_mock_grid(_DENSITY_);
#endif

         // I) Using the information of the bias and kernel, apply these two to a new DM field
         // in order to generate a new halo number denisty field. Measure Power spectrum pof the mock.


         // **************************************************************************** CREAT CATALOG *********************


#ifdef _GET_BAM_CAT_

#ifdef _ASSIGN_TO_CALIBRATION_
#ifdef _ASSIGN_TO_REFERENCE_
         this->fnameMOCK=this->params._Input_Directory_X()+this->params._Name_Catalog_Y();
#else
         this->fnameMOCK=this->Output_directory+"MOCK_TR_realization"+to_string(0)+"_MASY"+to_string(this->iMAS_Y)+"_Nft"+to_string(this->Nft)+"_z"+to_string(this->redshift);
#endif
#else
         this->fnameMOCK=this->Output_directory+"MOCK_TR_realization"+to_string(this->params._realization())+"_MASY"+to_string(this->iMAS_Y)+"_Nft"+to_string(this->Nft)+"_z"+to_string(this->redshift);
#endif


         this->makecat(this->patchy.stradd,this->fnameMOCK,gBaseRand,this->params._IC_index());

#endif //end get_bam_cat


#ifdef _FULL_VERBOSE_
         start_aux=start;
#endif

#ifndef _ONLY_POST_PROC_
       }
#endif


#endif //end of GET_BAM_REALIZATIONS

#endif //end of MOCK_MODE

     this->So.message_time(time_bam);
 }


   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************



#ifdef _GET_BAM_CAT_
#ifdef _USE_OMP_
void Bam::makecat(string stradd,string fnameMOCK,gsl_rng ** gBaseRand,int ir)
#else
void Bam::makecat(string stradd,string fnameMOCK,gsl_rng * gBaseRand,int ir)
#endif
{
    So.enter(__PRETTY_FUNCTION__);

    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);

    // In this function we asign coordinates from DM + random particles based on the mock number counts
     // Masses are also assigned and the collapse of randoms towards dm (to correct small scale clusterin) is also performed after that
     // Catalog is then written
#ifdef _FULL_VERBOSE_
     cout<<endl;
     So.message_screen("***********************************************************************");
     So.message_screen("********************GENERATING MOCK CATALOG****************************");
     So.message_screen("***********************************************************************");
     So.message_screen("***********************************************************************");
     cout<<endl;
#endif

#ifndef _ASSIGN_TO_REFERENCE_

#ifdef _FULL_VERBOSE_
     cout<<endl;
     this->So.message_screen("*****Assigning position and velocities to tracers using DM particles****");
     cout<<endl;
#endif

     real_prec redshift=this->s_cosmo_pars.cosmological_redshift;

     s_params_box_mas box;
     box.min1=this->params._xllc();
     box.min2=this->params._yllc();
     box.min3=this->params._zllc();
     box.Lbox=this->params._Lbox();
     box.Nft=this->Nft;
     box.d1= box.Lbox/static_cast<real_prec>(box.Nft);		/* grid spacing x-direction */
     box.d2= box.d1;
     box.d3= box.d1;
     box.NGRID=this->NGRID;


     ULONG N_dms=0;
     ULONG N_rand=0;


#ifdef _COLLAPSE_RANDOMS_AUX_
     s_params_box_mas box_collps;
     box_collps.min1=this->params._xllc();
     box_collps.min2=this->params._yllc();
     box_collps.min3=this->params._zllc();
     box_collps.Lbox=this->params._Lbox();
     box_collps.Nft=this->Nft_random_collapse;
     box_collps.d1= box_collps.Lbox/static_cast<real_prec>(box_collps.Nft);		/* grid spacing x-direction */
     box_collps.d2= box_collps.d1;
     box_collps.d3= box_collps.d1;
     box_collps.NGRID=(box_collps.Nft*box_collps.Nft*box_collps.Nft);
     vector<int>dm_count(box_collps.NGRID,0);
     vector<real_prec> x_dm_pos;
     vector<real_prec> y_dm_pos;
     vector<real_prec> z_dm_pos;
     vector<ULONG> dm_id;
     vector<real_prec> x_random_pos;
     vector<real_prec> y_random_pos;
     vector<real_prec> z_random_pos;
     vector<ULONG> ran_id;
#endif


     real_prec d1=this->params._d1();
     real_prec d2=this->params._d2();
     real_prec d3=this->params._d3();


     // ****************************DM POSITIONS *******************************************
     ULONG N_dm=this->NGRID;
     // Note that here the number of dm particles is that of the Ngrid

     vector<real_prec> posx(N_dm,0),posy(N_dm,0),posz(N_dm,0);
     vector<ULONG>index(N_dm,0);

     // Filenames of the bninary files containing the positions of the dark matter particles
     this->File.read_array(this->patchy.fnamePOSX+".dat",posx);
     this->File.read_array(this->patchy.fnamePOSY+".dat",posy);
     this->File.read_array(this->patchy.fnamePOSZ+".dat",posz);

     vector<real_prec> velx(N_dm,0),vely(N_dm,0),velz(N_dm,0);
     this->File.read_array(this->patchy.fnameVXpart+".dat",velx);
     this->File.read_array(this->patchy.fnameVYpart+".dat",vely);
     this->File.read_array(this->patchy.fnameVZpart+".dat",velz);

     real_prec vel_bias_dm=1.00;
     real_prec vel_bias_ran=1.00;

/*
     for(ULONG i=0;i<velz.size();++i)
         cout<<velx[i]<<"  "<<vely[i]<<"  "<<velz[i]<<endl;
     exit(0);
*/
     struct s_cell_info{
       vector<real_prec> posx_p;
       vector<real_prec> posy_p;
       vector<real_prec> posz_p;
       vector<real_prec> velx_p;
       vector<real_prec> vely_p;
       vector<real_prec> velz_p;
       vector<real_prec> mass_p;
       vector<ULONG> id_p;
     };

#endif   // end of ifndef _ASSIGN_TO_REFERENCE_

// we make a break here in the ifndef assign to reference in order to open the reference number count fiule

     vector<real_prec> MOCK_DEN_FIELD(this->NGRID,0);
     this->File.read_array(fnameMOCK+".dat",MOCK_DEN_FIELD);
     ULONG Nobjects_mock=get_nobjects(MOCK_DEN_FIELD);
     ULONG Ntracers = Nobjects_mock;  //reduntant, but still
     this->tracer.NOBJS=Nobjects_mock;


#ifndef _ASSIGN_TO_REFERENCE_

     vector<int> aux_cont(this->NGRID,0);
#ifdef _USE_OMP_
#pragma omp parallel num_threads(4)
     {
         vector<int> aux_cont_par(this->NGRID,0);
#pragma omp for nowait
#endif
     for(ULONG i=0;i<N_dm;++i)
       {
         real_prec x = static_cast<real_prec>(posx[i]);
         real_prec y = static_cast<real_prec>(posy[i]);
         real_prec z = static_cast<real_prec>(posz[i]);
         ULONG ind=grid_ID(&box, x,y,z);
         index[i]=ind;
         aux_cont_par[ind]++;  //count the number of dm particle in a cell
       }

#ifdef _USE_OMP_
#pragma omp critical
     for(ULONG i=0;i<this->NGRID;++i)
         aux_cont[i]+=aux_cont_par[i];  //count the number of dm particle in a cell
      }
#endif



     vector<bool>dm_used(this->NGRID,false);
     vector<int>dm_cases(this->NGRID,0);

     //ACA DEBEMOS CONTEMPLAR 4 CASOS:

     // I) EN UNA CELDA HAY MAS TRACERS QUE DM
     // II) EN UNA CELDA HAY DM QUE TRACERS
     // III) En una celda hay Tracers pero no hay DM
     // IV) Celdas vacias, deben permanecer vacias



#ifdef _FULL_VERBOSE_
     So.message_screen("Total number of tracers = ", Ntracers);
     So.message_screen("Total number of DM particles = ", posz.size());
     cout<<endl;
#endif


     ULONG empty_cells_original=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:empty_cells_original)
#endif
     for(ULONG id=0;id<this->NGRID;++id)
       if(MOCK_DEN_FIELD[id]==0)
         {
           empty_cells_original++;
           dm_cases[id]=4;
           dm_used[id]=false;
         }

     //   So.message_screen("Original number of empty cells ", empty_cells_original);

#ifdef _FULL_VERBOSE_
    So.message_screen("Assigning coordinates and velocities to tracers from coordiantes of ALPT-DM particles");
#endif
     vector<int> aux_cont1(this->NGRID,0);
     vector<int> aux_cont2(this->NGRID,0);

     // define container of structure to allocate coordinates within each cell
     vector<s_cell_info> cell_inf_dm(this->NGRID);

     // caso I: mas trcers que dm
     for(ULONG i=0;i<N_dm;++i)
       {
         real_prec x = static_cast<real_prec>(posx[i]);
         real_prec y = static_cast<real_prec>(posy[i]);
         real_prec z = static_cast<real_prec>(posz[i]);
         real_prec vx= static_cast<real_prec>(vel_bias_dm*velx[i]);
         real_prec vy= static_cast<real_prec>(vel_bias_dm*vely[i]);
         real_prec vz= static_cast<real_prec>(vel_bias_dm*velz[i]);
         ULONG id   =  index[i];    //identify the cell where this particle lives
#ifdef _COLLAPSE_RANDOMS_AUX_
         ULONG id_collapse = grid_ID(&box_collps, x,y,z);
#endif
         if(MOCK_DEN_FIELD[id]>0 &&  aux_cont[id]>0)  // si hay tracers en esta celda y dm tambien
           {
             if(MOCK_DEN_FIELD[id]> aux_cont[id]) //si  hay mas o igual número tracers que dm (y hay dm), tomar todas las dm que hay en cada celda. Al haber mas tracers que dm, vendrá la necesidad de tener randoms
               {
                 cell_inf_dm[id].posx_p.push_back(x);
                 cell_inf_dm[id].posy_p.push_back(y);
                 cell_inf_dm[id].posz_p.push_back(z);
                 cell_inf_dm[id].velx_p.push_back(vx);
                 cell_inf_dm[id].vely_p.push_back(vy);
                 cell_inf_dm[id].velz_p.push_back(vz);
                 cell_inf_dm[id].id_p.push_back(id);
                 dm_used[id]=true;
                 dm_cases[id]=1;
                 aux_cont1[id]++;
#ifdef _COLLAPSE_RANDOMS_AUX_
                 x_dm_pos.push_back(x);
                 y_dm_pos.push_back(y);
                 z_dm_pos.push_back(z);
                 dm_id.push_back(id_collapse);
                 dm_count[id_collapse]++;
                 N_dms++;
#endif
               }
             else
               { // caso II: se necesitan menor numero de tracers que el numero de dm en la celda. Ntr<=Ndm
                 if(aux_cont2[id]< MOCK_DEN_FIELD[id]) // el "if" es para tomar sólo los que necesitamos
                   {
                     cell_inf_dm[id].posx_p.push_back(x);
                     cell_inf_dm[id].posy_p.push_back(y);
                     cell_inf_dm[id].posz_p.push_back(z);
                     cell_inf_dm[id].velx_p.push_back(vx);
                     cell_inf_dm[id].vely_p.push_back(vy);
                     cell_inf_dm[id].velz_p.push_back(vz);
                     cell_inf_dm[id].id_p.push_back(id);
                     dm_used[id]=true;
                     dm_cases[id]=2;
                     aux_cont2[id]++;
#ifdef _COLLAPSE_RANDOMS_AUX_
                     x_dm_pos.push_back(x);
                     y_dm_pos.push_back(y);
                     z_dm_pos.push_back(z);
                     dm_id.push_back(id_collapse);
                     dm_count[id_collapse]++;
                     N_dms++;
#endif
                   }
               }
           }
       }
     index.clear();
     index.shrink_to_fit();

     // Get the number of randoms from case i and ii
     vector<int>Nrandom_tracers(this->NGRID,0);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG id=0;id<this->NGRID;++id)
       if(MOCK_DEN_FIELD[id]>0)
         if(1==dm_cases[id] || 2==dm_cases[id])
           Nrandom_tracers[id]=MOCK_DEN_FIELD[id]-cell_inf_dm[id].posx_p.size();


     ULONG Ndm_used_I=0;
     ULONG Nrandoms1=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:Ndm_used_I, Nrandoms1)
#endif
     for(ULONG id=0;id<this->NGRID;++id)
       if(dm_cases[id]==1)
         {
           Ndm_used_I+=aux_cont1[id];
           Nrandoms1+=Nrandom_tracers[id];
         }
     // So.message_screen("Number of dm used in case I",Ndm_used_I);
     // So.message_screen("Number of randoms demanded from case I ", Nrandoms1);


     //    ULONG Ndm_used_II=0;
     // #pragma omp parallel for reduction(+:Ndm_used_II)
     //    for(ULONG id=0;id<this->NGRID;++id)
     //      Ndm_used_II+=aux_cont2[id];
     // So.message_screen("Number of dm used in case II",Ndm_used_II);
     // So.message_screen("(case II demands no randoms)");

     posx.clear();posx.shrink_to_fit();
     posy.clear();posy.shrink_to_fit();
     posz.clear();posz.shrink_to_fit();

     // So.message_screen("Total Number of dm used",Ndm_used_I+Ndm_used_II);


     // caso 3, celdas con tracers pero sin dm particles -> Todas random
     //   So.message_screen("Getting randoms from case 3");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG id=0;id<this->NGRID;++id)
       {
         if(MOCK_DEN_FIELD[id]>0 && 0==aux_cont[id])
           if(dm_cases[id]!=1 && dm_cases[id]!=2)
             {
               dm_used[id]=false;
               Nrandom_tracers[id]=MOCK_DEN_FIELD[id]; //todas randoms
               dm_cases[id]=3;
             }
       }

     ULONG Nrandoms3=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:Nrandoms3)
#endif
     for(ULONG id=0;id<this->NGRID;++id)
       if(MOCK_DEN_FIELD[id]>0 && aux_cont[id]==0)
         if(3==dm_cases[id])
           Nrandoms3+=Nrandom_tracers[id];

     // So.message_screen("Total number of randoms requested from case 3", Nrandoms3);
     // So.message_screen("Total number of randoms requested", Nrandoms1+Nrandoms3);
     // So.message_screen("Randoms requested + Ndm", Nrandoms3+Nrandoms1+Ndm_used_I+Ndm_used_II);
     // So.message_screen("Number of original tracers",Nobjects_mock);
     Nrandoms1+=Nrandoms3;

     ULONG ncells_check=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ncells_check)
#endif
     for(ULONG id=0;id<this->NGRID;++id)
       {
         if(dm_cases[id]>=1 && dm_cases[id]<=4)
           ncells_check++;
       }

     if(ncells_check!=this->NGRID)
       {
         So.message_screen("Cells in cases", ncells_check);
         So.message_screen("Total number of cells", this->NGRID);
         exit(0);
       }


     int jthread=1;
     vector<s_cell_info> cell_inf_ran(this->NGRID);
     vector<bool>random_used(this->NGRID,false);

     if(Nrandoms1>0)
       {
#ifdef _FULL_VERBOSE_
         So.message_screen("Generating coordinates and velocities for random tracers");
#endif
         for(ULONG i=0; i<this->Nft; ++i)
           for(ULONG j=0; j<this->Nft; ++j)
             for(ULONG k=0; k<this->Nft; ++k)
               {
                 ULONG id= index_3d(i,j,k,this->Nft,this->Nft);
                 if(Nrandom_tracers[id]>0)
                   {
                     cell_inf_ran[id].posx_p.resize(Nrandom_tracers[id],0);
                     cell_inf_ran[id].posy_p.resize(Nrandom_tracers[id],0);
                     cell_inf_ran[id].posz_p.resize(Nrandom_tracers[id],0);
                     cell_inf_ran[id].velx_p.resize(Nrandom_tracers[id],0);
                     cell_inf_ran[id].vely_p.resize(Nrandom_tracers[id],0);
                     cell_inf_ran[id].velz_p.resize(Nrandom_tracers[id],0);
                     random_used[id]=true;
                     for(int ir=0;ir<Nrandom_tracers[id];++ir)
                       {
#ifdef _USE_OMP_
                         real_prec rx= static_cast<real_prec>(gsl_rng_uniform(gBaseRand[jthread]));
                         real_prec ry= static_cast<real_prec>(gsl_rng_uniform(gBaseRand[jthread]));
                         real_prec rz= static_cast<real_prec>(gsl_rng_uniform(gBaseRand[jthread]));
#else
                         real_prec rx= static_cast<real_prec>(gsl_rng_uniform(gBaseRand));
                         real_prec ry= static_cast<real_prec>(gsl_rng_uniform(gBaseRand));
                         real_prec rz= static_cast<real_prec>(gsl_rng_uniform(gBaseRand));
#endif


                         real_prec x = d1*static_cast<real_prec>(i) + d1*rx;
                         real_prec y = d2*static_cast<real_prec>(j) + d2*ry;
                         real_prec z = d3*static_cast<real_prec>(k) + d3*rz;
                         cell_inf_ran[id].posx_p[ir]=x;
                         cell_inf_ran[id].posy_p[ir]=y;
                         cell_inf_ran[id].posz_p[ir]=z;
                         cell_inf_ran[id].velx_p[ir]=vel_bias_ran*this->patchy.linInterp(x,y,z,velx);
                         cell_inf_ran[id].vely_p[ir]=vel_bias_ran*this->patchy.linInterp(x,y,z,vely);
                         cell_inf_ran[id].velz_p[ir]=vel_bias_ran*this->patchy.linInterp(x,y,z,velz);
                         cell_inf_ran[id].id_p.push_back(id);
#ifdef _COLLAPSE_RANDOMS_AUX_
                         x_random_pos.push_back(x);
                         y_random_pos.push_back(y);
                         z_random_pos.push_back(z);
                         ULONG id_collapse = grid_ID(&box_collps, x,y,z);
                         ran_id.push_back(id_collapse);
                         N_rand++;
#endif

                       }
                   }
               }
         So.DONE();
       }
     velx.clear();velx.shrink_to_fit();
     vely.clear();vely.shrink_to_fit();
     velz.clear();velz.shrink_to_fit();



#ifdef _COLLAPSE_RANDOMS_AUX_
     int max_count=get_max(dm_count); // get the maximum number of dm particles in one cell.
#ifdef _FULL_VERBOSE_
     So.message_screen("Maximum number of DM particles found in one cell (low resolution) = ", max_count);
#endif

     // Initialize this to NO_NUM, a negative number, for not all elements of the container
     // dm_index_cell will be filled and initializing with 0 makes confusion, for 0 is a valid/used entry


#define NO_NUM -999
     vector<vector<int> > dm_index_cell(box_collps.NGRID, vector<int>(max_count,NO_NUM));

     dm_count.clear();
     dm_count.resize(box_collps.NGRID,0);
#ifdef _FULL_VERBOSE_
     So.message_screen("Getting ids in cells");
#endif


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG idm=0;idm<N_dms;++idm)
       {
         ULONG id=dm_id[idm];
         dm_index_cell[id][dm_count[id]]=idm;
#ifdef _USE_OMP_
#pragma atomic update
#endif
         dm_count[id]++;
       }
     this->So.DONE();
     dm_id.clear();dm_id.shrink_to_fit();

     vector<int>dm_index_closer_tot(N_rand,0);
#ifdef _FULL_VERBOSE_
     So.message_screen("Identifying closest DM particles for randoms");
#endif

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0; i< N_rand; ++i)//loop over randoms
       {
         ULONG count =0;
         int id_ran=ran_id[i];
         int N_dm_cell=dm_count[id_ran];
         vector<int>i_r_to_dm_dist;
         vector<int>n_index_tot;
         if(N_dm_cell>0)// if there are dm particles in the low res cell where this random is located, then
           {
             for(ULONG jc = 0; jc < N_dm_cell; ++jc)//loop over the dm in that cell where the random is
               {
                 ULONG jdm=dm_index_cell[id_ran][jc]; //indice entre (0,Ndm-1) que identifica a cada particula DM dentro de la celda id_ran
                 if(jdm!=NO_NUM)
                   {
                     real_prec xdd=x_random_pos[i]-x_dm_pos[jdm];
                     real_prec ydd=y_random_pos[i]-y_dm_pos[jdm];
                     real_prec zdd=z_random_pos[i]-z_dm_pos[jdm];
                     real_prec dist_dm_r=_get_modulo(xdd,ydd,zdd); //distance between rand and each central
                     ULONG i_dist=static_cast<ULONG>(floor(100.0*dist_dm_r)); // done in order to sort vectors below. 100 is arbitrary
                     i_r_to_dm_dist.push_back(i_dist);
                     n_index_tot.push_back(jdm);
                     count++;
                   }
               }
             sort_1dvectors(i_r_to_dm_dist, n_index_tot);
             dm_index_closer_tot[i]=n_index_tot[0]; // el primer elemento (0) es el menor
           }
         else
           So.message_screen("No dm particles found in the cell correspoding to random ", i, ". You might want to increase the cell-size");
       }


     So.DONE();
     dm_count.clear();
     dm_count.shrink_to_fit();
     ran_id.clear();ran_id.shrink_to_fit();

#endif  // end collapse_random_aux


     // Before write to file we assign masses. This replaces the call of the function in bamrunner
     // Here we pass the coordinates to a prop vctor and pass it as vector& to the mass assignment function,
     // to then pass it to collapse randoms

     ULONG Nobjects_mock_rc=0;



     this->tracer.Halo.resize(this->tracer.NOBJS);
#ifdef _FULL_VERBOSE_
     this->So.message_screen("Filling structure with dm positions:");
#endif

// Problems with this parallelization. LEAVE COMMENTED!!!!!!!!!1
     for(ULONG id=0;id<this->NGRID;++id)
       {
         if(MOCK_DEN_FIELD[id]>0)
           if(true==dm_used[id])
             for(int in=0;in < cell_inf_dm[id].posx_p.size(); ++in)
               {
                 this->tracer.Halo[Nobjects_mock_rc].coord1 = cell_inf_dm[id].posx_p[in];
                 this->tracer.Halo[Nobjects_mock_rc].coord2 = cell_inf_dm[id].posy_p[in];
                 this->tracer.Halo[Nobjects_mock_rc].coord3 = cell_inf_dm[id].posz_p[in];
                 this->tracer.Halo[Nobjects_mock_rc].vel1 = cell_inf_dm[id].velx_p[in];
                 this->tracer.Halo[Nobjects_mock_rc].vel2 = cell_inf_dm[id].vely_p[in];
                 this->tracer.Halo[Nobjects_mock_rc].vel3 = cell_inf_dm[id].velz_p[in];
                 this->tracer.Halo[Nobjects_mock_rc].identity = 1 ;
                 this->tracer.Halo[Nobjects_mock_rc].GridID = cell_inf_dm[id].id_p[in];
                 Nobjects_mock_rc++;
               }
       }
     So.DONE();
#ifdef _FULL_VERBOSE_
     So.message_screen("Number of tracers associated to DM particles =", Nobjects_mock_rc);
#endif
     this->tracer.Ntracers_dm=Nobjects_mock_rc;
     this->tracer_ref.Ntracers_dm=Nobjects_mock_rc;
     this->tracer.fraction_tracer_from_dm=static_cast<real_prec>(Nobjects_mock_rc)/(static_cast<real_prec>(Nobjects_mock_rc)+static_cast<real_prec>(Nrandoms1)) ;
     this->tracer_ref.fraction_tracer_from_dm=this->tracer.fraction_tracer_from_dm;

#ifdef _FULL_VERBOSE_
     So.message_screen("(", 100.0*this->tracer.fraction_tracer_from_dm, "%)");
     So.message_screen("Freeing memory");
#endif
     cell_inf_dm.clear(); cell_inf_dm.shrink_to_fit();
     So.DONE();


     ULONG N_dms_b= Nobjects_mock_rc;

     if(Nrandoms1>0)
       {
#ifdef _FULL_VERBOSE_
         this->So.message_screen("Adding information of random particles:");
#endif
//#ifdef _USE_OMP_
//#pragma omp parallel for reduction (+:Nobjects_mock_rc)  // lio
//#endif
         for(ULONG id=0;id<this->NGRID;++id)
           {
             if(MOCK_DEN_FIELD[id]>0)
               if(true==random_used[id])
                 for(ULONG in=0;in<cell_inf_ran[id].posx_p.size();++in)
                   {
                     this->tracer.Halo[Nobjects_mock_rc].coord1=cell_inf_ran[id].posx_p[in];
                     this->tracer.Halo[Nobjects_mock_rc].coord2=cell_inf_ran[id].posy_p[in];
                     this->tracer.Halo[Nobjects_mock_rc].coord3=cell_inf_ran[id].posz_p[in];
                     this->tracer.Halo[Nobjects_mock_rc].vel1=cell_inf_ran[id].velx_p[in];
                     this->tracer.Halo[Nobjects_mock_rc].vel2=cell_inf_ran[id].vely_p[in];
                     this->tracer.Halo[Nobjects_mock_rc].vel3=cell_inf_ran[id].velz_p[in];
                     this->tracer.Halo[Nobjects_mock_rc].identity=-1 ;
                     this->tracer.Halo[Nobjects_mock_rc].GridID=cell_inf_ran[id].id_p[in];
                     Nobjects_mock_rc++; // This counter keeps on counting starting from the number of dm particles above
                   }
           }
#ifdef _FULL_VERBOSE_
         So.DONE();
         So.message_screen("Number of tracers associated to random particles =", Nrandoms1);
         this->tracer.Ntracers_ran=Nrandoms1;
         this->tracer_ref.Ntracers_ran=Nrandoms1;

         this->tracer.fraction_tracer_from_random=static_cast<real_prec>(Nrandoms1)/(static_cast<real_prec>(N_dms_b)+static_cast<real_prec>(Nrandoms1));
         this->tracer_ref.fraction_tracer_from_random=this->tracer.fraction_tracer_from_random;

         So.message_screen("(", 100.0*this->tracer.fraction_tracer_from_random, "%)");

         So.message_screen("Freeing memory");
#endif
         cell_inf_ran.clear(); cell_inf_ran.shrink_to_fit();
         random_used.clear();
         random_used.shrink_to_fit();
#ifdef _FULL_VERBOSE_
         So.DONE();
#endif
       }



#ifdef _COLLAPSE_RANDOMS_AUX_

     // THE GrdID kept in the sturcture Tracer is the *ORIGINAL*, not the one computed after the collapse of randoms

     cout<<endl;
#ifdef _FULL_VERBOSE_
     this->So.message_screen("Using fraction of distance to closest dm particle = ", this->Distance_fraction);
     this->So.message_screen("Computing new position of randoms:");
#endif

     ULONG counter=0;
//#ifdef _USE_OMP_
//#pragma omp parallel for reduction(+:counter)  // problem with paralelization and counter, leave commented
//#endif
     for(ULONG i=0; i<Ntracers; ++i)
       {
         if(this->tracer.Halo[i].identity<0)  //This means, use the randoms
           {

             //total index (i.e, from 0 to N_dms) of the closest dm particle to the ramdom i. Usefil to retrieve coordinates
             ULONG index_dm_closer_a=dm_index_closer_tot[counter];
             //Cartesian coordinates of the closest DM particle
             real_prec xdm=x_dm_pos[index_dm_closer_a];
             real_prec ydm=y_dm_pos[index_dm_closer_a];
             real_prec zdm=z_dm_pos[index_dm_closer_a];
             // redefine ran coords to the ref sistem of its closest dm particle:
             real_prec new_x=x_random_pos[counter]-xdm;
             real_prec new_y=y_random_pos[counter]-ydm;
             real_prec new_z=z_random_pos[counter]-zdm;
             real_prec dist_random_to_dm=_get_module(new_x,new_y,new_z);             // get the distance:
             //get the angular coordinates:
             real_prec theta=acos(new_z/dist_random_to_dm);
             real_prec phi=atan2(new_y,new_x);
             // get the new distance between the random and its closest dm particle
             real_prec new_distance = dist_random_to_dm*this->Distance_fraction;
             // transfor to cartesian given a reduced distance and return ot origin of coords:
             // Assign the new cartesian coordinates
             this->tracer.Halo[i].coord1=new_distance*sin(theta)*cos(phi)+xdm;
             this->tracer.Halo[i].coord2=new_distance*sin(theta)*sin(phi)+ydm;
             this->tracer.Halo[i].coord3=new_distance*cos(theta)+zdm;
             counter++;
             }
       }

#ifdef _FULL_VERBOSE_
     this->So.DONE();
#endif
     x_random_pos.clear();
     x_random_pos.shrink_to_fit();
     y_random_pos.clear();
     y_random_pos.shrink_to_fit();
     z_random_pos.clear();
     z_random_pos.shrink_to_fit();
     x_dm_pos.clear();
     x_dm_pos.shrink_to_fit();
     y_dm_pos.clear();
     y_dm_pos.shrink_to_fit();
     z_dm_pos.clear();
     z_dm_pos.shrink_to_fit();

#endif // end of _COLLAPSE_RANDOMS_AUX_



#else   // else for ifndef assign_to_reference, meant to allocate the vector of struturures tracer_ref.Halo[] with the same elements of the tracer.Halo
     // In order to save time, we are assuming that the mock number count used in the test "assign to reference" is the same obtained from the catalog reference
     // THis allocation cannot be done here, for the reference has not been read. It will be read in the function get_X_function()
//     so we make it there interanally.

#endif  //end of ifndef _ASSIGN_TO_REFERENCE




#ifdef _WRITE_BINARY_BAM_FORMAT_
     string outputFileName=this->Output_directory+"CAT_realization"+to_string(this->params._IC_index())+"_"+this->patchy.stradd+string(".dat");
#else
     string outputFileName=this->Output_directory+"CAT_realization"+to_string(this->params._IC_index())+"_"+this->patchy.stradd+string(".txt");
#endif

/*
         Params params_aux=this->params;
         this->tracer.type_of_object="TRACER_MOCK_ONLY_COORDS"; // THis name is important.
         params_aux.mass_assignment_scheme ="TSC";
         params_aux.MAS_correction = true;
         params_aux.SN_correction = true;
         params_aux.vel_units_g = "kmps";

         params_aux.Nft = 800;
         params_aux.input_type ="catalog";
#if defined (_COLLAPSE_RANDOMS_AUX_) || defined (_COLLAPSE_RANDOMS_)
         params_aux.Name_survey="BAM_CAT_collapsed_f"+to_string(this->params._Distance_fraction());
#else
         params_aux.Name_survey="BAM_CAT";

#endif
         PowerSpectrumF Power(params_aux);
         Power.compute_power_spectrum(false,this->tracer.Halo);
*/



#ifdef _ASSIGN_PROPERTY_

     // Get info from the reference
     this->get_X_function();

            // If Im ahave thte reference with the mass cut alrady done:
         //     this->tracer.type_of_object="TRACER_MOCK_ONLY_COORDS"; // THis name is important.
    // if not, use this:
       this->tracer.type_of_object="TRACER_MOCK"; // THis name is important.


/**
#ifdef _MASS_ASSIGNMENT_TO_REFERENCE_
     cout<<endl;

     So.message_warning("***NOTE: The reference is currently being used as mock in order to re-assign masses anc do tests, see line", __LINE__);
//     this->params.file_catalogue=this->Output_directory+"newcat.txt";
     this->params.file_catalogue=this->params._file_catalogue();
//    this->params.i_coord1_g=0;
//     this->params.i_coord2_g=1;
//     this->params.i_coord3_g=2;
 //    this->params.i_vmax_g=-1;
 //    this->params.i_mass_g=-1;
 //    this->tracer.set_params_catalog(this->params);

#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
#ifdef _SET_CAT_WITH_MASS_CUT_
     this->tracer.read_catalog(this->params._file_catalogue(),pow(10,params._LOGMASSmin())*params._MASS_units());
#elif defined (_SET_CAT_WITH_VMAX_CUT_)
     this->tracer.read_catalog(this->params._file_catalogue(),params._VMAXmin());
#else
     So.message_warning("No property for cut defining sample. Check preproc definitions. Code ends here.");
     exit(0);
#endif

#else

#ifdef _SET_CAT_WITH_MASS_CUT_
     this->tracer.read_catalog(this->params._file_catalogue(),pow(10,params._LOGMASSmin())*params._MASS_units(),static_cast<real_prec>(BIG_NUMBER));
#elif defined (_SET_CAT_WITH_VMAX_CUT_)
     this->tracer.read_catalog(this->params._file_catalogue(),params._VMAXmin(),static_cast<real_prec>(BIG_NUMBER));
#else
     So.message_warning("No property for cut defining sample. Check preproc definitions. Code ends here.");
     exit(0);
#endif

#endif

#endif  // end if _MASS_ASSIGNMENT_TO_REFERENCE_


**/






#ifdef _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_

#ifdef _USE_VMAX_AS_OBSERVABLE_
     this->assign_tracer_property_new_new(true, _VMAX_);
#else
     this->assign_tracer_property_new_new(true, _MASS_);
#endif
#else
     this->assign_tracer_mass_new();
#endif
     this->So.DONE();
     this->tracer.type_of_object="TRACER_MOCK"; // Once masses are assigned, change the identity of this->tracer. This name is important.

#ifdef _ASSIGN_MASS_POST_
     // Here we can now complement the property assignment, using the ingormation of vmax already assigned
#ifdef _FULL_VERBOSE_
     So.message_screen("******************************************************");
     So.message_screen("Assigning halo masses using VMAX information::");
     So.message_screen("******************************************************");
     cout<<endl;
     cout<<endl;
#endif
     this->params.i_vmax_g=8;// this number only needs to be positive
     this->tracer.type_of_object="TRACER_MOCK"; // THis name is important.
     this->get_X_function_complement(_MASS_);   // Learn P(M|V,delta) from reference
     this->assign_tracer_property_new_new(false, _MASS_);  // Drawn M from P(M|V=Vref, delta=delta_ref)
#endif

#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
#ifdef _FULL_VERBOSE_
     So.message_screen("******************************************************");
     So.message_screen("Assigning Rs using VMAX, Mvir information::");
     So.message_screen("******************************************************");
#endif
     this->params.i_mass_g=9;// this number only needs to be positive
     this->tracer.type_of_object="TRACER_MOCK"; // THis name is important.
     this->get_X_function_complement(_RS_); // Learn P(RS|M,V) from reference
     this->assign_tracer_property_new_new(false, _RS_);
#endif

#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
#ifdef _FULL_VERBOSE_
     So.message_screen("******************************************************");
     So.message_screen("Assigning Spin using VMAX, Mvir information::");
     So.message_screen("******************************************************");
#endif
     this->params.i_spin_g=10;// this number only needs to be positive
     this->tracer.type_of_object="TRACER_MOCK"; // THis name is important.
     this->get_X_function_complement(_SPIN_); // Learn P(SPIN|M,V)
     this->assign_tracer_property_new_new(false, _SPIN_);
#endif
#endif // end if ASSIGN_PROPERTY


     // Collapse randoms:  this applies if collapse_randoms_aux is undef
#ifdef _COLLAPSE_RANDOMS_
     this->collapse_randoms();
#endif

     this->patchy.fnameTRACERCAT=outputFileName;


#ifdef _FULL_VERBOSE_
     So.message_warning("Power Spectrum of catalog is commented @ line ", __LINE__);
#endif

      Params params_aux=params;
      params_aux.file_catalogue=outputFileName;
      params_aux.Nft=360;
      params_aux.Name_survey="BAM_CAT";
      params_aux.mass_assignment_scheme="TSC";
      params_aux.weight_with_mass=false;

#ifdef _MASS_WEIGHT_POWER_
      params_aux.weight_with_mass=true;
      params_aux.file_power="MOCK_mass_weight";
      params_aux.SN_correction=false;
#else
      params_aux.file_power="MOCK"+to_string(this->params._IC_index());
      params_aux.SN_correction=true;

#endif
      PowerSpectrumF Power(params_aux);
      Power.compute_power_spectrum(false,this->tracer.Halo);


/*
      // If cuts are done, the residuals are to be computed using the last mcut
      this->Power_NEW_MW.resize(this->Nft/2/this->ndel_data,0);
      #ifdef _USE_OMP_
      #pragma omp parallel for
      #endif
      for(int i=0;i<this->Nft/2/this->ndel_data;++i)
      this->Power_NEW_MW[i]=Power._pk0(i);

      real_prec residuals=0;
      int count_p=0;

      for(int i=0;i<this->Nft/2/this->ndel_data;++i)
      if(this->Power_REF_MW[i]>0 && this->kvec[i]<KMAX_RESIDUALS_low)
      {
      count_p++;
      residuals+= fabs(this->Power_NEW_MW[i]/this->Power_REF_MW[i]-1.0);
      }
      residuals/=static_cast<real_prec>(count_p)/100.0;
      So.message_screen("At KMAX_RES = ", KMAX_RESIDUALS_low,"  h / Mpc");
      So.message_screen("Residuals from power spectrum = ", residuals, "%");


      residuals=0;
      count_p=0;
      for(int i=0;i<this->Nft/2/this->ndel_data;++i)
      if(this->Power_REF_MW[i]>0 && this->kvec[i]<KMAX_RESIDUALS_high)
      {
      count_p++;
      residuals+= fabs(this->Power_NEW_MW[i]/this->Power_REF_MW[i]-1.0);
      }
      residuals/=static_cast<real_prec>(count_p)/100.0;
      So.message_screen("At KMAX_RES = ", KMAX_RESIDUALS_high,"  h / Mpc");
      So.message_screen("Residuals from power spectrum = ", residuals, "%");



      residuals=0;
      count_p=0;111
      for(int i=0;i<this->Nft/2/this->ndel_data;++i)
      if(this->Power_REF_MW[i]>0)
      {
      count_p++;
      residuals+= fabs(this->Power_NEW_MW[i]/this->Power_REF_MW[i]-1.0);
      }
      residuals/=static_cast<real_prec>(count_p)/100.0;
      So.message_screen("At Niquyst frequency" );
      So.message_screen("Residuals from power spectrum  = ", residuals, "%");
     */


#ifdef _GET_DIST_MIN_SEP_MOCK_
    this->tracer_ref.get_distribution_min_separations(this->tracer.type_of_object, this->ncells_info);
#endif



#ifdef _WRITE_BAM_CATALOGS_
#ifdef _WRITE_BINARY_BAM_FORMAT_
    this->tracer.write_catalog_bin(outputFileName.c_str());
#else
    this->tracer.write_catalog_ascii(outputFileName.c_str());
#endif

#ifdef _FULL_VERBOSE_
	 So.message_screen("Freeing memmory from tracer");
#endif
     this->tracer.Halo.clear();
     this->tracer.Halo.shrink_to_fit();
     So.DONE();
#endif

}
#endif



// ********************************************************************************************************************************************************************************
// ********************************************************************************************************************************************************************************
// ********************************************************************************************************************************************************************************
// ********************************************************************************************************************************************************************************
#ifdef MOCK_MODE
#ifndef _DEPRECATED_
void Bam::correct_for_exclusion(ULONG LENGHTdm){
    So.enter(__PRETTY_FUNCTION__);

    //A bin in DM properties have contributions from different sections of the volume (i.e, differnet IDGrids).
    // We study here the distribution of separations within objects in a THETA bin.
    // The idea is to move the masses in that particular theta_bin such that the distribution of the mock follow that of the reference, computed here.

    // An alternative idea (used in Hadron) is to  define a Mass threshold such that masses in a theta-bin above that mass threshold
    //are to be associated to distant cells.

    const gsl_rng_type *Tn;
    gsl_rng *rn ;
    Tn = gsl_rng_default;
    rn = gsl_rng_alloc (Tn);

    this->tracer.masses_in_cells_min_sep.clear();
    this->tracer.masses_in_cells_min_sep.shrink_to_fit();
    this->tracer.masses_in_cells_min_sep.resize(LENGHTdm);

    vector<int>number_in_theta_mock(LENGHTdm,0);
 #pragma omp parallel for
    for(ULONG i=0;i<LENGHTdm; ++i)
      number_in_theta_mock[i]=this->dm_properties_bins_mock[i].masses_bin_properties.size();
    So.DONE();

    // Obtain the list of pairs of masses allocated in a bin of theta that are within a distance min_halo_sep, EXCLUSION_SCALE. output in  this->tracer.masses_in_cells_min_sep
    this->tracer.get_masses_of_pairs_in_min_separation_bin_in_theta_bin(this->min_halo_separation, this->dm_properties_bins_mock);

    for(int ih=0;ih< LENGHTdm;++ih) //Loop over the bins in Theta
      {
        real_prec M1_max_ref; // min and max of reference masses in theta_bin in the min_separation bin
        real_prec M2_max_ref;
        real_prec M1_max_aux=-1e5; // min and max of reference masses in theta_bin in the min_separation bin
        real_prec M2_max_aux=-1e5;

        for(int j=0; j< this->tracer_ref.masses_in_cells_min_sep[ih].M1.size();++j)
          {
            M1_max_ref=max(this->tracer_ref.masses_in_cells_min_sep[ih].M1[j], M1_max_aux);
            M1_max_aux=M1_max_ref;
            M2_max_ref=max(this->tracer_ref.masses_in_cells_min_sep[ih].M2[j], M1_max_aux);
            M2_max_aux=M2_max_ref;
          }
        // search for the masses M1 above the max in the bin. It is enough to look at one particle. These are the masses to be swaped

        for(int j=0; j< this->tracer.masses_in_cells_min_sep[ih].M1.size();++j)
          {
            real_prec M1=this->tracer.masses_in_cells_min_sep[ih].M1[j];
            if(M1>M1_max_ref) //allocate the galID of this mass
              this->tracer.masses_in_cells_min_sep[ih].mass_to_swap.push_back(this->dm_properties_bins_mock[ih].GalID_bin_properties[j]);
          }
     }



    So.message_screen("Swaping masses...");

    for(int ih=0;ih< LENGHTdm;++ih)//Loop over the bins in Theta
      {
        int N_props_in_bin = number_in_theta_mock[ih]; //Available Masses in bin of THETA, regardless of the separation between objects

        for(int j=0; j< this->tracer.masses_in_cells_min_sep[ih].mass_to_swap.size();++j)//loopover the masses to be swaped
          {

           ULONG old_ig=this->tracer.masses_in_cells_min_sep[ih].mass_to_swap[j]; // galaxy index IDg [0,NOBJS) of the mass requested to be swaped
           real_prec old_mass=this->tracer.Halo[old_ig].mass; //mass requested to be swaped

           bool baux=false;
           while(baux==false)
            {
              int i_mass_halo_label= gsl_rng_uniform_int(rn,N_props_in_bin); // pick-up randomly some other mass in the same bin-theta
              ULONG new_ig=this->dm_properties_bins_mock[ih].GalID_bin_properties[i_mass_halo_label];  // galaxy index IDg [0,NOBJS) of an randomly selected object
              real_prec new_mass=this->tracer.Halo[new_ig].mass; //mass of the randomly selected object

              //swap masses:
              if(new_ig!=this->tracer.masses_in_cells_min_sep[ih].mass_to_swap[j])
               {
                 this->tracer.Halo[old_ig].mass=new_mass; // assign mass
                  this->tracer.Halo[new_ig].mass=old_mass;
                  number_in_theta_mock[ih]--;
                  baux=true;
               }
             }
          }
     }

    So.DONE();
    number_in_theta_mock.clear();number_in_theta_mock.shrink_to_fit();
    this->dm_properties_bins_mock.clear();
    this->dm_properties_bins_mock.shrink_to_fit();

}
#endif
#endif
