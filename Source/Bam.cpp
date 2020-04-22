//##################################################################################
//##################################################################################
/** @file Bam.cpp
 *
 *  @brief Generation of mock catalogs of DM tracers
 *  based on the BAM method.
 *  @author: Andrés Balaguera-Antolínez, Francisco-Shu Kitaura, IAC, 2017-2019
 */
# include "../Headers/Bam.h"

//##################################################################################
//##################################################################################
void Bam::set_params_bam()
{

  So.message_screen("Loading parameters for BAM");
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
  this->Apply_Rankordering=this->params._Apply_Rankordering();
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
  this->nmax_y_onecell=50;
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

#ifdef _USE_MASS_KNOTS_
  this->N_bias_properties++;
#endif

#ifdef _USE_CWC_
  this->N_bias_properties++;
#endif

#ifdef _USE_VEL_KNOTS_V_
  this->N_bias_properties++;
#endif

#ifdef _USE_CWC_V_
  this->N_bias_properties++;
#endif


#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_DELTA2_)
  this->N_bias_properties++;
#endif

#if defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_DELTA3_)
  this->N_bias_properties++;
#endif

#if defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_S2_)
  this->N_bias_properties++;
#endif


#if defined (_USE_INVARIANT_SHEAR_VFIELD_I_) || defined(_USE_NABLA2DELTA_)
  this->N_bias_properties++;
#endif

#if defined (_USE_INVARIANT_SHEAR_VFIELD_II_) || defined (_USE_S2DELTA_)
  this->N_bias_properties++;
#endif

#if defined (_USE_INVARIANT_SHEAR_VFIELD_III_) || defined (_USE_S3_)
  this->N_bias_properties++;
#endif


  So.message_screen("BAM is using", this->N_bias_properties,"properties to characterise the bias");

  this->step_mass_assignemt=0; //initialize it

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


}

//##################################################################################
//##################################################################################
void Bam::get_cosmo()
{
  So.message_screen("Getting cosmological derived-parameters");

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


#ifdef _USE_PATCHY_  // This is repeated here as it is in patchy
  this->s_cosmo_info.growth_factor=this->Cosmo.growth_factor(redshift, (void *)&this->s_cosmo_pars)/this->Cosmo.growth_factor(0.0,(void *)&this->s_cosmo_pars);
  this->s_cosmo_info.growth_index=this->Cosmo.growth_index(redshift, (void *)&this->s_cosmo_pars);
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
#ifdef _VERBOSE_
  this->So.write_cosmo_parameters((void *)&this->s_cosmo_pars, (void *)&this->s_cosmo_info);
#endif
#endif

  So.DONE();


  this->cwclass.s_cosmo_info=this->s_cosmo_info;

}



//##################################################################################
//##################################################################################
void Bam::show_params()
{
  cout<<BOLDYELLOW;
  cout<<"***************************************************************************"<<endl;
  cout<<"***************************************************************************"<<endl;
  cout<<"BAM                                                                       *"<<endl;
  cout<<"Bias Assignment Method for galaxy/halo mock catalogs                      *"<<endl;
  cout<<"Input values of parameters in parameter file                              *"<<endl;
  cout<<"***************************************************************************"<<endl;
  cout<<"***************************************************************************"<<RESET<<endl;

  cout<<"redshift                   = "<<this->params._redshift()<<endl;
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

  cout<<BOLDYELLOW<<"***************************************************************************"<<RESET<<endl;
  cout<<"Preprocessor directives"<<endl;
#ifdef _USE_OMP_
  cout<<"_USE_OMP_                         : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_USE_OMP_                         : "<<RED<<"undefined"<<RESET<<endl;
#endif


#ifdef _READ_REF_CATALOG_
  cout<<"_READ_REF_CATALOG_               : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_READ_REF_CATALOG_               : "<<RED<<"undefined"<<RESET<<endl;
#endif


#ifdef _USE_MASS_TRACERS_
  cout<<"_USE_MASS_TRACERS_               : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_USE_MASS_TRACERS_               : "<<RED<<"undefined"<<RESET<<endl;
#endif

#ifdef _USE_VELOCITIES_TRACERS_
  cout<<"_USE_VELOCITIES_TRACERS_         : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_USE_VELOCITIES_TRACERS_         : "<<RED<<"undefined"<<RESET<<endl;
#endif


#ifdef _ASSIGN_PROPERTY_
  cout<<"_ASSIGN_PROPERTY_                  : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_ASSIGN_PROPERTY_                  : "<<RED<<"undefined"<<RESET<<endl;
#endif


#ifdef _COLLAPSE_RANDOMS_
  cout<<"_COLLAPSE_RANDOMS_               : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_COLLAPSE_RANDOMS_               : "<<RED<<"undefined"<<RESET<<endl;
#endif




#ifdef BIAS_MODE
  cout<<"BIAS_MODE                         : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"BIAS_MODE                         : "<<RED<<"undefined"<<RESET<<endl;
#endif
#ifdef MOCK_MODE
  cout<<"MOCK_MODE                         : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"MOCK_MODE                         : "<<RED<<"undefined"<<RESET<<endl;
#endif

#ifdef _USE_PATCHY_
  cout<<"_USE_PATCHY_                      : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_USE_PATCHY_                      : "<<RED<<"undefined"<<RESET<<endl;
#endif


#ifdef _DO_BAM_CALIBRATION_
  cout<<"_DO_BAM_CALIBRATION_              : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_DO_BAM_CALIBRATION_              : "<<RED<<"undefined"<<RESET<<endl;
#endif


#ifdef _GET_BAM_REALIZATIONS_
  cout<<"_GET_BAM_REALIZATIONS_            : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_GET_BAM_REALIZATIONS_            : "<<RED<<"undefined"<<RESET<<endl;
#endif


#ifdef _RUN_TEST_
  cout<<"_RUN_TEST_                          : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_RUN_TEST_                          : "<<RED<<"undefined"<<RESET<<endl;
#endif
#ifdef _DYNAMICAL_SAMPLING_
  cout<<"_DYNAMICAL_SAMPLING_              : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_DYNAMICAL_SAMPLING_              : "<<RED<<"undefined"<<RESET<<endl;
#endif
#ifdef _USE_CWC_
  cout<<"_USE_CWC_                         : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_USE_CWC_                         : "<<RED<<"undefined"<<RESET<<endl;
#endif
#ifdef _USE_CWC_INSIDE_LOOP_
  cout<<"_USE_CWC_INSIDE_LOOP_             : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_USE_CWC_INSIDE_LOOP_             : "<<RED<<"undefined"<<RESET<<endl;
#endif

#ifdef _USE_MASS_KNOTS_
  cout<<"_USE_MASS_KNOTS_                  : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_USE_MASS_KNOTS_                  : "<<RED<<"undefined"<<RESET<<endl;
#endif
#ifdef _USE_VELOCITIES_
  cout<<"_USE_VELOCITIES_                  : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_USE_VELOCITIES_                  : "<<RED<<"undefined"<<RESET<<endl;
#endif


#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
  cout<<"_USE_INVARIANT_TIDAL_FIELD_II_     : "<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<"_USE_INVARIANT_TIDAL_FIELD_II_     : "<<RED<<"undefined"<<RESET<<endl;
#endif
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
  cout<<"_USE_INVARIANT_TIDAL_FIELD_III_    : "<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<"_USE_INVARIANT_TIDAL_FIELD_III_    : "<<RED<<"undefined"<<RESET<<endl;
#endif

#ifdef _USE_TIDAL_ANISOTROPY_
  cout<<"_USE_TIDAL_ANISOTROPY_            : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_USE_TIDAL_ANISOTROPY_            : "<<RED<<"undefined"<<RESET<<endl;
#endif


#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
  cout<<"_USE_INVARIANT_SHEAR_VFIELD_I_    : "<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<"_USE_INVARIANT_SHEAR_VFIELD_I_    : "<<RED<<"undefined"<<RESET<<endl;
#endif

#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
  cout<<"_USE_INVARIANT_SHEAR_VFIELD_II_   : "<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<"_USE_INVARIANT_SHEAR_VFIELD_II_   : "<<RED<<"undefined"<<RESET<<endl;
#endif

#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
  cout<<"_USE_INVARIANT_SHEAR_VFIELD_III_  : "<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<"_USE_INVARIANT_SHEAR_VFIELD_III_  : "<<RED<<"undefined"<<RESET<<endl;
#endif

#ifdef _USE_S2_
  cout<<"_USE_S2_                          : "<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<"_USE_S2_                          : "<<RED<<"undefined"<<RESET<<endl;
#endif

#ifdef _USE_S3_
  cout<<"_USE_S3_                          : "<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<"_USE_S3_                          : "<<RED<<"undefined"<<RESET<<endl;
#endif



#ifdef _USE_S2DELTA_
  cout<<"_USE_S2DELTA_                     : "<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<"_USE_S2DELTA_                     : "<<RED<<"undefined"<<RESET<<endl;
#endif

#ifdef _USE_DELTA2_
  cout<<"_USE_DELTA2_                      : "<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<"_USE_DELTA2_                      : "<<RED<<"undefined"<<RESET<<endl;
#endif

#ifdef _USE_DELTA3_
  cout<<"_USE_DELTA3_                      : "<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<"_USE_DELTA3_                      : "<<RED<<"undefined"<<RESET<<endl;
#endif

#ifdef _USE_NABLA2DELTA_
  cout<<"_USE_NABLA2DELTA_                 : "<<BOLDGREEN<<"defined"<<endl;
#else
  cout<<"_USE_NABLA2DELTA_                 : "<<RED<<"undefined"<<RESET<<endl;
#endif



#ifdef _GET_BAM_CAT_
  cout<<"_GET_BAM_CAT_                     : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_GET_BAM_CAT_                     : "<<RED<<"undefined"<<RESET<<endl;
#endif


#ifdef _NGP2CIC_Y_
  cout<<"_NGP2CICL_Y_                      : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_NGP2CIC_Y_                       : "<<RED<<"undefined"<<RESET<<endl;
#endif


#ifdef _MODIFY_LIMITS_
  cout<<"_MODIFY_LIMITS_                   : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_MODIFY_LIMITS_                   : "<<RED<<"undefined"<<RESET<<endl;
#endif


#ifdef _MASS_ASSIGNMENT_TO_REFERENCE_
  cout<<"_MASS_ASSIGNMENT_TO_REFERENCE_    : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_MASS_ASSIGNMENT_TO_REFERENCE_    : "<<RED<<"undefined"<<RESET<<endl;
#endif

#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
  cout<<"_USE_MIN_SEPARATIONS_IN_CELLS_         : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_USE_MIN_SEPARATIONS_IN_CELLS_    : "<<RED<<"undefined"<<RESET<<endl;
#endif


#ifdef _CORRECT_FOR_EXCLUSION_
  cout<<"_CORRECT_FOR_EXCLUSION_         : "<<BOLDGREEN<<"defined"<<RESET<<endl;
#else
  cout<<"_CORRECT_FOR_EXCLUSION_    : "<<RED<<"undefined"<<RESET<<endl;
#endif



  cout<<BOLDYELLOW<<"***************************************************************************"<<RESET<<endl;
  cout<<BOLDYELLOW<<"***************************************************************************"<<RESET<<endl;
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

  So.message_screen("Inizializing some Fourier vectors");
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
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
void Bam::get_power_spectrum(string type)
{

  So.message_screen("Measuring power spectrum of",type);

  this->params.mass_assignment_scheme="NGP";
  this->params.MAS_correction=false;

  this->params.dir_output =this->Output_directory;  //This line is important
  // in the case in which the output dir is changing. Otherwise, the O(k) will be written in the dir read initially by the class Params


  if(this->step==0)
    this->params.Name_survey=type;

  if(this->step <=this->N_iterations_Kernel)
    this->params.Name_survey=type+"_iteration"+to_string(this->step);
  else
    this->params.Name_survey=type+"_realization"+to_string(this->step - (this->N_iterations_Kernel)+this->N_dm_initial-1);

  if(type=="DM_REF" || type=="DM_REF_NEW")
    {
      this->params.input_type="density";
      this->params.SN_correction=false;
      this->params.MAS_correction=true;

      if(0==this->iMAS_X)
        {
          this->params.mass_assignment_scheme="NGP";
          this->params.MAS_correction=false;
        }
      if(1==this->iMAS_X)
        this->params.mass_assignment_scheme="CIC";
      else if (2==this->iMAS_X)
        this->params.mass_assignment_scheme="TSC";
      else if (3==this->iMAS_X)
        this->params.mass_assignment_scheme="PSC";


      PowerSpectrumF cPSF(this->params);
      cPSF.compute_power_spectrum_grid(this->delta_X_ini);
      So.DONE();
      cPSF.write_power_and_modes();
      this->Power_DM_REF.resize(this->Nft/2/this->ndel_data, 0);
      for(int i=0;i<this->Power_DM_REF.size(); ++i)
        this->Power_DM_REF[i]=cPSF._pk0(i);

    }
  else if(type=="DM_iteration" || type=="DM_real" || type=="DM_KONV")
    {
      this->params.input_type="delta";
      this->params.SN_correction=false;
      this->params.mass_assignment_scheme="CIC";
      this->params.MAS_correction=true;
      PowerSpectrumF cPSF(this->params);
      cPSF.compute_power_spectrum_grid(this->delta_X);
      So.DONE();
      cPSF.write_power_and_modes();
    }
  else if(type=="TR_MOCK")
    {

      this->params.SN_correction=false;
      //      if(this->step>0) // Do not correct for SN at the first one, for it goes below the ref in some cases and the kernel might get >1
      if(_COUNTS_==this->Name_Property_Y)
        this->params.SN_correction=true;

      this->params.input_type="density";

      PowerSpectrumF cPSF(this->params);

      cPSF.compute_power_spectrum_grid(this->delta_Y_new);
      this->Power_NEW.resize(this->Nft/2/this->ndel_data, 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<this->Power_NEW.size(); ++i)
        this->Power_NEW[i]=cPSF._pk0(i);


      So.DONE();
      cPSF.write_power_and_modes();

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
      this->params.input_type="density";
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
      if(this->Name_Property_Y=="COUNTS" || this->Name_Property_Y=="DENSITY")
        this->params.SN_correction=true;

      this->params.mass_assignment_scheme="NGP";
      this->params.MAS_correction=false;

      this->params.input_type="density";

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
void Bam::GetKernel(bool rejection, real_prec exponent)
{


  this->So.message_screen("Metropolis-Hasting (mode-by-mode) for kernel");

  if(false==this->use_iteration_ini) // initialized as false in bamrunner if iteration_ini >0
    {

      int NTHREADS = omp_get_max_threads();
      const gsl_rng_type *  T;
      gsl_rng * r ;

      vector<int>vseeds(this->N_iterations_Kernel,0);
      for(int i=0;i<vseeds.size();++i)vseeds[i]=(i+1);

      int nmodes_f=this->Nft/this->ndel_data/2;

      vector<real_prec> weight(nmodes_f,1.0);
      vector<real_prec> aux_power(nmodes_f,1.0);



      //#pragma omp parallel private(T, r)
      {
        gsl_rng_env_setup();
        gsl_rng_default_seed=1+521563+vseeds[this->step_mass_assignemt];
        //    gsl_rng_default_seed=vseeds[omp_get_thread_num()];

        T = gsl_rng_ranlux;
        r = gsl_rng_alloc (T);

        // Select kernel according to the previous step. Update the previous kernel witht eh current step

        real_prec power_ratio;
        real_prec partial_ratio=0;
        real_prec deltak=this->kvec[1]-this->kvec[0];

        //#pragma omp parallel for


        for(ULONG i=0;i<nmodes_f;++i)
          {

#ifdef _DO_BAM_CALIBRATION_
            real_prec Power_ref=this->Power_REF[i];
            real_prec Power_new=this->Power_NEW[i];
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
                power_ratio=(Power_new == 0 || Power_ref == 0) ? 1.0 : pow(Power_ref/Power_new,exponent);   //idelly for exponent = 1.0
                if(i< N_MODES)
                  partial_ratio+=power_ratio;
              }

            if(false==rejection)
              {
                weight[i]=power_ratio;
                this->power_ratio_unsmoothed[i]=power_ratio;
              }
            else
              {
                // ***********************************
                // These applies if we use instead the likelihhod and take ratios of it.
                // The approach here resembles more a MCMC with  L=exp(-chi**2) for each kbin.
                // The variance at each k-bin is assumed to be Gaussian, and we neglect here shot-noise
                real_prec sigma = Power_ref == 0? 1.0 : sqrt((4.*M_PI)/(deltak*pow(this->kvec[i],2)*pow(this->Lbox,3)))*(Power_ref); // ESTE PARA EL KERNEL DE LA DM
                real_prec new_H = exp(-0.5*pow((Power_ref- Power_new)/sigma, 2)); // proposed likelihood
                real_prec Power_old= static_cast<double>(Power_ref)/static_cast<double>(this->power_ratio_unsmoothed[i]);  // power of the previous step
                real_prec old_H = exp(-0.5*pow((Power_ref-Power_old)/sigma, 2));   // likelihood of the previous step
                real_prec ratio_diff=static_cast<double>(new_H)/static_cast<double>(old_H);
                // ***********************************
                real_prec x= gsl_rng_uniform (r);
                if(x  < min(1.0, static_cast<double>(ratio_diff)))
                  {
                    aux_power[i]=power_ratio;    //recorded just to print out
                    weight[i]= power_ratio; // accepted. Use this ratio (at this spherical shell) in this step to update the kernel below*/
                  }
                else
                  {
                    aux_power[i]=this->power_ratio_unsmoothed[i]; //recorded just to print out
                    weight[i]=1.0; // Rejected. Will use the ratio of the last step when updating the kernel below
                  }
                this->power_ratio_unsmoothed[i]=power_ratio;
              }

          }
        So.DONE();


        real_prec residuals;
        for(int i=0;i<this->power_ratio_unsmoothed.size();++i)
          residuals+=fabs(this->power_ratio_unsmoothed[i]-1.0)/static_cast<real_prec>(this->power_ratio_unsmoothed.size());
        So.message_screen("Residuals at this iteration (%) =",100.0*residuals);

        //        this->output_res<<this->lambdath<<"\t"<<this->lambdath_v<<"\t"<<100.0*residuals<<endl;
        So.message_screen("Average of ratio P(k)_ref / P(k)_new = " , partial_ratio/static_cast<real_prec>(N_MODES));
        cout<<endl;
      }


      vseeds.clear();
      vseeds.shrink_to_fit();

#ifdef _DO_BAM_CALIBRATION_
      string kernel_file_or=this->Output_directory+"RatioT_iteration"+to_string(this->step)+".txt";
      this->File.write_to_file(kernel_file_or, this->kvec,aux_power,power_ratio_unsmoothed);
      aux_power.clear(); //to release memory before going out of scope
      aux_power.shrink_to_fit();
#else
      string kernel_file_or=this->Output_directory+"RatioT_mass_assignment_iteration"+to_string(this->step_mass_assignemt)+".txt";
      this->File.write_to_file(kernel_file_or, this->kvec,aux_power,power_ratio_unsmoothed);
      aux_power.clear(); //to release memory before going out of scope
      aux_power.shrink_to_fit();
#endif

#ifdef _DO_BAM_CALIBRATION_
      So.message_screen("Updating BAM Kernel");
#else
      So.message_screen("Updating Halo-Mass assignment Kernel");
#endif

      vector<real_prec> coords(this->Nft,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<Nft ;++i)
        coords[i]= (i<=Nft/2? static_cast<real_prec>(i): -static_cast<real_prec>(Nft-i));

      if(true==rejection)
        {
#pragma omp parallel for collapse(3)
          for(ULONG i=0;i< this->Nft; ++i)
            for(ULONG j=0;j< this->Nft; ++j)
              for(ULONG k=0;k< this->Nft/2+1; ++k)
                {
                  ULONG ind=index_3d(i,j,k, this->Nft, this->Nft/2+1);
                  real_prec kv=sqrt(pow(coords[i],2)+pow(coords[j],2)+pow(coords[k],2));
                  int kmod=static_cast<int>(floor(kv));
                  if(kmod<this->Nft/2/this->ndel_data)
                    this->Kernel[ind]*=weight[kmod]; //weight is 1 if no improvement;
                }
        }
      else
        {
#pragma omp parallel for collapse(3)
          for(ULONG i=0;i< this->Nft; ++i)
            for(ULONG j=0;j< this->Nft; ++j)
              for(ULONG k=0;k< this->Nft/2+1; ++k)
                {
                  ULONG ind=index_3d(i,j,k, this->Nft, this->Nft/2+1);
                  real_prec kv=sqrt(pow(coords[i],2)+pow(coords[j],2)+pow(coords[k],2));
                  int kmod=static_cast<int>(floor(kv));
                  if(kmod<this->Nft/2/this->ndel_data)
                    this->Kernel[ind] = weight[kmod];
                }
        }

//      weight.clear(); // to release memory before going out of scope
//      weight.shrink_to_fit();

      vector<real_prec>kernel_updated(this->Power_NEW.size(), 0);// to print out the shell-averaged kernel
      vector<int>nmodes(this->Power_NEW.size(), 0);

#pragma omp parallel for collapse(3)
      for(ULONG i=0;i< this->Nft; ++i)
        for(ULONG j=0;j< this->Nft; ++j)
          for(ULONG k=0;k< this->Nft/2+1; ++k)
            {
              ULONG ind=index_3d(i,j,k, this->Nft, this->Nft/2+1);
              real_prec kv=sqrt(pow(coords[i],2)+pow(coords[j],2)+pow(coords[k],2));
              int kmod=static_cast<int>(floor(kv));
              if(kmod< nmodes.size())
                {
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
      So.DONE();


#pragma omp parallel for
      // Get the shell-averaged version of the Kernel in FOurier space
      for(ULONG i=0;i< kernel_updated.size(); ++i)
        kernel_updated[i]/=static_cast<real_prec>(nmodes[i]);

#ifdef _DO_BAM_CALIBRATION_
      string file_kernel=this->Output_directory+"Kernel_Fourier_iteration"+to_string(this->step)+".txt";
      this->File.write_to_file(file_kernel, this->kvec, kernel_updated);
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
          log_smooth(this->kvec, kernel_updated);
          log_smooth(this->kvec, kernel_updated);
          log_smooth(this->kvec, kernel_updated);
          log_smooth(this->kvec, kernel_updated);
          string kernel_file=this->Output_directory+"Kernel_Fourier_smoothed_iteration"+to_string(this->step)+".txt";
          sal.open(kernel_file.c_str());
          for(ULONG i=0;i<Nft/2;++i)sal<<this->kvec[i]<<"  "<<kernel_updated[i]<<endl;
          So.message_screen("Writting smootehd version of kernel in file", kernel_file);
          sal.close();

#ifdef _SMOOTHED_KERNEL_LAST_ITERATION_
        }
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      // Assign smoothed shell averaged Kernel to 3D kernel
      for(ULONG i=0;i< this->Nft; ++i)
        for(ULONG j=0;j< this->Nft; ++j)
          for(ULONG k=0;k< this->Nft/2+1; ++k)
            {
              ULONG ind=index_3d(i,j,k, this->Nft, this->Nft/2+1);
              real_prec kv=sqrt(pow(coords[i],2)+pow(coords[j],2)+pow(coords[k],2));
              int kmod=static_cast<int>(floor(kv));
              if(kmod<this->Nft/2/this->ndel_data)
                this->Kernel[ind]=kernel_updated[kmod];
            }
#endif


      // All other vectors defined here and not released
      // are destroyed her when going out of scope

#ifndef _TEST_THRESHOLDS_RESIDUALS_
      if(this->step == this->N_iterations_Kernel)
        this->File.write_array(this->Output_directory+"Bam_Kernel", this->Kernel);
#endif



  if(this->step == this->N_iterations_Kernel)
    {
     complex_prec *data_out= (complex_prec *)fftw_malloc(2*this->NTT*sizeof(real_prec));
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
     fftw_free(data_out);

    }
    this->use_iteration_ini=false;
  }





  else if(this->iteration_ini>0 && true==this->use_iteration_ini)
    {
      So.message_screen("Reading Kernel from iteration", this->iteration_ini);
      this->Kernel.clear();
      this->Kernel.resize(this->NTT, 0.0);
      this->File.read_array(this->Output_directory+"Bam_Kernel.dat", this->Kernel);
      this->use_iteration_ini=true; // we set it true for we will need it for the bias
    }

}

//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################


void Bam::Konvolve(vector<real_prec> &in, vector<real_prec>&out, string type)
{
#ifdef _GET_BAM_REALIZATIONS_
  if(this->step_mass_assignemt==0)
    So.message_screen("Generating new DM density field by convolution of input DM with input Kernel");
  else
    So.message_screen("Generating new DM density field by convolution of DM with upated mass-weighterd kernel ");
#else
  So.message_screen("Generating new DM density field by convolution of input DM with updated Kernel");
#endif


  complex_prec *data_out= (complex_prec *)fftw_malloc(2*this->NTT*sizeof(real_prec));
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->NTT;++i){
    data_out[i][REAL]=0;
    data_out[i][IMAG]=0;
  }

  do_fftw_r2c(this->Nft,in, data_out);


#ifdef _EXTRAPOLATE_VOLUME_
  real_prec correction_factor = 1.0;
#else
  real_prec correction_factor = 1.00;
#endif

  real_prec we=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:we)
#endif
  for(ULONG ind=0;ind< this->NTT ;++ind)
    {
      real_prec cor = this->Kernel[ind]*correction_factor;

#ifdef _USE_OMP_
#pragma omp atomic update
#endif
      data_out[ind][REAL]*=cor;

#ifdef _USE_OMP_
#pragma omp atomic update
#endif
      data_out[ind][IMAG]*=cor;
      we+=cor;
    }

  vector<real_prec>aux(this->NGRID,0);
  do_fftw_c2r(this->Nft, data_out,aux);

  // Here I have to correct for the normalization of the kernel
  // for the function  do_fftw_c2r returns the transform normalized by the NGRID, so I divide by NGRID and by multiply by 2 we
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->NGRID;i++)
    aux[i]/=(static_cast<real_prec>(this->NGRID)/static_cast<real_prec>(2.0*we));

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->NGRID;i++)
    out[i]=aux[i];

  So.DONE();

  if(this->iteration_ini>0 && this->step!=this->iteration_ini)
    {
      auto out_it = std::find(std::begin(this->output_at_iteration), std::end(this->output_at_iteration), this->step);
      if (out_it != std::end(this->output_at_iteration))
        {
          // Write the kernel interpolated to 3D in Fourier space.
          complex_prec *kern= (complex_prec *)fftw_malloc(2*this->NTT*sizeof(real_prec));
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(ULONG i=0;i<this->NTT;++i)
            {
              kern[i][REAL]=this->Kernel[i];
              kern[i][IMAG]=0;
            }


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(ULONG i=0;i<this->NGRID;i++)
            aux[i]=0;

          do_fftw_c2r(this->Nft,kern,aux);

          // para no escribir lo que ha leído
          string file_kernel=this->Output_directory+"3DKernel_iteration"+to_string(this->step);
          this->File.write_array(file_kernel, aux);

          fftw_free(kern);
        }
    }




  fftw_free(data_out);

}

// //##################################################################################
// //##################################################################################
// ##################################################################################
void Bam::get_new_min_max_properties()
{

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

#ifdef _USE_SAT_FRACTION_
      this->s_mins.prop0_sf=get_min(this->delta_Y_SAT_FRACTION);
      this->s_maxs.prop0_sf=get_max(this->delta_Y_SAT_FRACTION);
      this->s_deltas.prop0_mass=(this->s_maxs.prop0_sf-this->s_mins.prop0_sf)/static_cast<real_prec>(this->NY_SAT_FRAC);
      So.message_screen("Maximim mass sat_fraction = ", this->s_maxs.prop0_sf);
      So.message_screen("Minimum mass sat_fraction = ", this->s_mins.prop0_sf);
#endif
    }//close else of if _COUNTS_!=this->Name_Property_Y


#ifdef _MODIFY_LIMITS_

  // Gert new extremes for delta_X
  this->s_mins.prop1=get_min(this->delta_X);
  this->s_maxs.prop1=get_max(this->delta_X);
  this->s_deltas.prop1=(this->s_maxs.prop1-this->s_mins.prop1)/static_cast<real_prec>(this->new_nbins_x);
  So.message_screen("Minimum log(2+ð) DM = ", this->s_mins.prop1);
  So.message_screen("Maximum log(2+ð) DM = ", this->s_maxs.prop1);
  //  So.message_screen("N_bins  ð DM = ", this->new_nbins_x);


#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
  this->s_mins.prop4=get_min(this->cwclass.Invariant_TF_II);
  this->s_maxs.prop4=get_max(this->cwclass.Invariant_TF_II);
  this->s_deltas.prop4=(this->s_maxs.prop4-this->s_mins.prop4)/static_cast<real_prec>(N_C_BIN1);
  So.message_screen("Maximim InvTF2 = ", this->s_maxs.prop4);
  So.message_screen("Minimum InvTF2 = ", this->s_mins.prop4);

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
  So.message_screen("Maximim InvTFIII = ", this->s_maxs.prop5);
  So.message_screen("Minimum InvTFII = ", this->s_mins.prop5);
  So.message_screen("Delta InvTFII = ", this->s_deltas.prop5);

#elif defined _USE_DELTA3_
  this->s_mins.prop5=get_min(this->cwclass.DELTA3);
  this->s_maxs.prop5=get_max(this->cwclass.DELTA3);
  this->s_deltas.prop5=(this->s_maxs.prop5-this->s_mins.prop5)/static_cast<real_prec>(N_C_BIN2);
  So.message_screen("Maximim ð³ = ", this->s_maxs.prop5);
  So.message_screen("Minimum ð³ = ", this->s_mins.prop5);
#endif


#ifdef _USE_TIDAL_ANISOTROPY_
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

#ifdef _DO_BAM_CALIBRATION_
  if(this->step==0)
    {
#endif

#ifdef _USE_SAT_FRACTION_
      this->s_mins.prop0_sf=get_min(this->delta_Y_SAT_FRAC);
      this->s_maxs.prop0_sf=get_max(this->delta_Y_SAT_FRAC);
      //         this->NY_SAT_FRAC = this->s_maxs.prop0_sf; // THIS HAS BEEN ALREADY DONE IN LINE 1371
      this->s_deltas.prop0_mass=(this->s_maxs.prop0_sf-this->s_mins.prop0_sf)/static_cast<real_prec>(this->NY_SAT_FRAC);
      So.message_screen("Maximim mass sat_fraction = ", this->s_maxs.prop0_sf);
      So.message_screen("Minimum mass sat_fraction = ", this->s_mins.prop0_sf);
#endif
#ifdef _DO_BAM_CALIBRATION_
    }
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

void  Bam::get_min_max_X_Y()
{
  // This function redefines the values of delta_X_min etc in case. If this is not called, the input values are used
  // The filds delta X and delta Y ahave been already, if requiested, converted to overdensities

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



  cout<<endl;
  this->So.message_screen("Reading input *reference* files");
  cout<<endl;

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
  this->File.read_array_t<PrecType_X>(file_X,this->delta_X);
#ifdef _CONVERT_CIC_TO_NGP_
  convert_cic_to_ngp(delta_X,delta_X);
#endif
#endif



#ifdef _USE_TRACER_HR_
  this->delta_Y_HR.resize(NGRID_NEW_HR,0);
  this->File.read_array_t<PrecType_Y>(file_Y_HR, this->delta_Y_HR);
#endif

  //  this->delta_Y_new.resize(this->NGRID,0); //Isn't this resized in get_mock?
  ULONG empty_cells=0;
  this->delta_Y.resize(this->NGRID,0);
  this->File.read_array_t<PrecType_Y>(file_Y, this->delta_Y);
  for(ULONG i=0;i < this->NGRID; ++i)
    if(this->delta_Y[i]==0)
      empty_cells++;
  So.message_screen("Number of empty cells from number density field =", empty_cells);


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
  empty_cells=0;
  for(ULONG i=0;i < this->NGRID; ++i)
    {
      if(this->delta_Y_MASS[i]<1.0)
        empty_cells++;
      this->delta_Y_MASS[i]=log10(NUM_IN_LOG+ this->delta_Y_MASS[i]/MASS_SCALE);
    }
  So.message_screen("Number of empty cells from mass density field =", empty_cells);
  So.message_screen("Using log(2+total central mass in cells). Check line ",__LINE__);
#endif

#endif
#endif


#ifdef _USE_SAT_FRACTION_
  this->delta_Y_SAT_FRAC.resize(this->NGRID,0);
  this->File.read_array_t<PrecType_Y>(file_Y_sat_frac, this->delta_Y_SAT_FRAC);

  this->nmax_y_sat_onecell=get_max(delta_Y_SAT_FRAC);
  this->NY_SAT_FRAC=this->nmax_y_sat_onecell+1;
  //  this->Mean_density_Y_SAT_FRAC=get_mean(this->delta_Y_SAT_FRAC);
  //  get_overdens(this->delta_Y_SAT_FRAC,this->Mean_density_Y_SAT_FRAC, this->delta_Y_SAT_FRAC);
  //for(ULONG i=0;i < this->NGRID; ++i)
  //      this->delta_Y_SAT_FRAC[i]=log10(NUM_IN_LOG+this->delta_Y_SAT_FRAC[i]);

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



  if(true==this->Apply_Rankordering || this->N_iterations_dm>0)
    {
      this->delta_X_REF_PDF.resize(NGRID_NEW,0);
      string file_X_ref_pdf=this->Input_Directory_X_REF+this->Name_Catalog_X_REF_PDF; // TBDep
      this->File.read_array_t<PrecType_X>(file_X_ref_pdf, this->delta_X_REF_PDF);
    }


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
  this->File.read_array(this->Output_directory+"Bam_Kernel.dat",kernel);
  convolvek(this->Nft,this->delta_X, kernel,this->delta_X);
#endif

  // Initialize the DM density field with the input (REF) density field
  this->delta_X_ini=this->delta_X;

}

//  #################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  #################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################


void Bam::analyze_input()
{
  So.message_screen("Getting some statistics from input files: ");
  cout<<endl;

  // X:

  real_prec nmean_X=0.;

#ifndef _EXTRAPOLATE_VOLUME_
  nmean_X=get_nobjects(this->delta_X);
  this->N_objects_X=nmean_X;
  if(this->Name_Property_X==_COUNTS_)
    So.message_screen("Total number of X objects =", nmean_X);

  nmean_X=static_cast<real_prec>(nmean_X)/static_cast<real_prec>(this->NGRID);
  if(this->Name_Property_X=="COUNTS")
    So.message_screen("Mean number of X objects in cells =", nmean_X);
  else
    So.message_screen("Mean property X in cells =", nmean_X);

  if(this->Name_Property_X=="COUNTS")
    So.message_screen("Mean number density =", nmean_X*static_cast<real_prec>(this->NGRID)/pow(this->Lbox,3), "(Mpc / h )⁻³");
  else
    So.message_screen("Mean property density =", nmean_X*static_cast<real_prec>(this->NGRID)/pow(this->Lbox,3));


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

      So.message_screen("Maximum number of X particles in one cell =", this->nmax_x_onecell);
      So.message_screen("Estimated Poisson signal-to-noise ratio =", sqrt(this->nmax_x_onecell));

    }


#ifdef MOCK_MODE
#ifdef _GET_POWER_REFS_
  this->get_power_spectrum("DM_REF");  //gets power from delta_X

#endif
#endif

#endif  // END EXTRAPOLATE

  // ----------------------------------------------------------------------------------------------------------------------
  // Y:


  real_prec nmean_Y=0;
  nmean_Y=get_nobjects(this->delta_Y);

  if(_COUNTS_==this->Name_Property_Y)
    {
      this->N_objects_Y=nmean_Y;
#ifdef _EXTRAPOLATE_VOLUME_
      So.message_screen("Note that these are properties of the reference simulation (i.e, smaller volume)");
#endif

      cout<<endl;
      So.message_screen("Total number of Y objects =", nmean_Y);
    }


  real_prec LLBOX;
  ULONG NNGRID;
#ifdef _VERBOSE_

#ifdef _EXTRAPOLATE_VOLUME_
  LLBOX=this->Lbox_low;
  NNGRID= static_cast<ULONG>(this->Nft_low*this->Nft_low*this->Nft_low);
#else
  LLBOX=this->Lbox;
  NNGRID= this->NGRID;
#endif

  nmean_Y=static_cast<real_prec>(nmean_Y)/static_cast<real_prec>(NNGRID);

  if(this->Name_Property_Y=="COUNTS")
    So.message_screen("Mean number of Y objects in cells =", nmean_Y);
  else
    So.message_screen("Mean property Y in cells =", nmean_Y);

  this->Mean_density_Y=nmean_Y;

  if(this->Name_Property_Y=="COUNTS")
    So.message_screen("Mean number density =", nmean_Y*static_cast<real_prec>(NNGRID)/pow(LLBOX,3), "(Mpc / h )⁻³");
  else
    So.message_screen("Mean property-density =", nmean_Y*static_cast<real_prec>(NNGRID)/pow(LLBOX,3));


#endif
  this->new_Name_Property_Y=this->Name_Property_Y;

  if(_COUNTS_==this->Name_Property_Y)
    if(this->iMAS_Y==0)
      {
        this->nmax_y_onecell=static_cast<int>(get_max<real_prec>(this->delta_Y)); // Maximum number of particles in cells
        this->So.message_screen("Maximum number of Y particles in one cell =", this->nmax_y_onecell);
        this->So.message_screen("Estimated Maximum Poisson signal-to-noise ratio =", sqrt(this->nmax_y_onecell));

#ifdef _USE_MASS_FIELD_
        this->So.message_screen("Maximum tracer mass in one cell =", get_max<real_prec>(this->delta_Y_MASS));
#endif

#ifdef _USE_SAT_FRACTION_
        this->So.message_screen("Maximum number of satellites in one cell =", get_max<real_prec>(this->delta_Y_SAT_FRAC));
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
            string fileY=this->Output_directory+"PDF_NC"+"_Y_"+this->YNAME+"_"+this->new_Name_Property_Y+"_MASY"+to_string(this->iMAS_Y)+"_Nft"+to_string(this->Nft)+"_SmoothingScale"+to_string(this->smscale)+"_z"+to_string(this->redshift)+"_LambdaTH"+to_string(this->lambdath)+"_CW"+to_string(this->tstruct)+"_CWclass.txt";
            this->File.write_to_file_i(fileY,this->PDF_NC_Y);
          }
      }



  cout<<endl;
#ifdef MOCK_MODE
#ifndef _EXTRAPOLATE_VOLUME_
  this->get_power_spectrum("TR_REF");
#ifdef _USE_TRACER_HR_
  this->delta_Y_HR.clear();
  this->delta_Y_HR.shrink_to_fit();
#endif

#endif



#endif


#ifdef MOCK_MODE
#ifndef _GET_BAM_REALIZATIONS_


  real_prec lss_bias=0;
  for(int i=0;i<N_MODES;++i)
    if(this->Power_DM_REF[i]!=0)
      lss_bias += this->Power_REF[i]/this->Power_DM_REF[i];

  lss_bias/=static_cast<real_prec>(N_MODES);
  So.message_screen("Large-Scale Bias: P(k)_tracer / P(k)_dm  = ", sqrt(lss_bias));
  cout<<endl;

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




  // **********************************************************************
  // **********************************************************************

#ifdef _RUN_TEST_
  dark_matter_to_halos_analytical(this->delta_X, this->delta_Y);
  exit(0);
#endif



  // **********************************************************************
  // **********************************************************************
  // **********************************************************************
  // Preparation of the DM density field

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


      for(int iteration_dm = 0; iteration_dm < this->N_iterations_dm ; ++iteration_dm)
        {
          So.message_screen("Iteration on the DM field", iteration_dm);

          // i) Get power of DM, getting ready for the kernel. NO rejection in this case, use exponent 0.5
          GetKernel(false, 0.5);
          vector<real_prec>auxv(this->NGRID, 0);
          auxv=this->delta_X;

          // v) Convolve the DM field field with the kernel obtaiende from the reference
          this->Konvolve(auxv,this->delta_X, "DENSITY");

          auxv.clear();
          auxv.shrink_to_fit();

          this->get_power_spectrum("DM_iteration");
          this->File.write_array(this->Output_directory+"density_DM_converted_dm_iteration_"+to_string(iteration_dm), this->delta_X);

        }

    }

  // **********************************************************************
  // **********************************************************************
  // **********************************************************************
  // **********************************************************************

  // Assign the density field delta_X_ini to *density*  delta_X
  // Bellow, if delta_X is to be converted to overdensity, we reassign delta_X_ini to overdensity
  this->delta_X_ini=this->delta_X;  // BAM starts with the outcome of this loop. This line is important.


  // **********************************************************************
  // **********************************************************************
  // **********************************************************************
  // **********************************************************************
  // **********************************************************************
  // **********************************************************************
#ifdef _USE_BINARY_MASK_
  vector<real_prec> binary_mask(this->NGRID,0);
  this->File.read_array(this->Name_binary_mask, binary_mask);
#endif

  // Convert density to delta for both the reference and the approximated field
  if(true==this->Convert_Density_to_Delta_X)
    {
      this->new_Name_Property_X="DELTA";

#ifdef _ADD_POISSON_TO_CIC_
      gsl_rng_env_setup();
      gsl_rng_default_seed=75;
      const gsl_rng_type *  T= gsl_rng_ranlux;
      gsl_rng * r = gsl_rng_alloc (T);
      get_overdens(this->delta_X,  nmean_X, this->delta_X);  // get overdensity
      So.message_screen("Converting DENSITY -> Poisson(DENSITY):");
      for(ULONG i=0;i<delta_X.size();++i)
         delta_X[i]=gsl_ran_poisson(r,MEAN_NEW*(1.+delta_X[i]));   //recover density
      So.DONE();
       {
        this->params.input_type="density";
        this->params.Name_survey="DM_Poisson";
        this->params.SN_correction=false;
        this->params.mass_assignment_scheme="CIC";
        this->params.MAS_correction=true;
        PowerSpectrumF dPSFa(this->params);
        dPSFa.compute_power_spectrum_grid(this->delta_X); //ojo, el argumento debe ser DENSIDAD
        dPSFa.write_power_and_modes();
      }
#endif


      So.message_screen("Converting DENSITY into DELTA for X:");
#ifndef _USE_BINARY_MASK_
#ifdef _ADD_POISSON_TO_CIC_
      get_overdens(this->delta_X,  this->delta_X);
#else
      get_overdens(this->delta_X, nmean_X, this->delta_X);
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

      if(true==this->Apply_Rankordering) // Only if RO is to be applied on the DM, convert the reference to delta
        get_overdens(this->delta_X_REF_PDF,this->delta_X_REF_PDF);
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


  // Convert the TRACER to delta

  if(true==this->Convert_Density_to_Delta_Y)
    {
      So.message_screen("Converting DENSITY into DELTA for Y:");
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
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################
//  ####################################################################################################################################################################


void Bam::get_BAM_DM()
{

  // This function converts the DM DELTA fields to log 2+delta and compute the CWC if requested
  // The BAM kernel is initialized

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
      So.message_screen("Computing PDF from dark matter");
      calc_pdf("log", this->NGRID,this->NX, this->Xmax, this->Xmin, this->delta_X, this->pdf_ref);
      string filex=this->Output_directory+"pdf_X_"+this->XNAME+"_MASX"+to_string(this->iMAS_X)+rest_file;
      this->File.write_to_file(filex, xbaux,pdf_ref);
    }
#endif
  // **********************************************************************
  // At this point, we had already said delta_X_ini =delta_x WITHOUT transforming to logs
  // (we convert delta_X to log10(num_in_log+delta_X) in the next function get_BIAS)
#ifdef _UNDER_BIASED_
  if(false==this->used_once)
    {
      So.message_screen("Using TR as DM for b<1 tests. Trick to treat under-biased populations");
      So.message_screen("Test in line", __LINE__, ", function get_BAM_DM()");

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i = 0;i < this->NGRID ;++i)
        this->delta_X_ini[i] = (this->delta_X[i]<-1 ?  0 :  log10(num_in_log_x+ static_cast<real_prec>(this->delta_Y[i])));

      this->used_once=true;
    }
#endif

  // **********************************************************************
  // Do the CWC class for the first iteration, or in the last , when te target field is loaded
  if(this->step==0)
    {

#ifdef _USE_CWC_
      this->cwclass.do_CWC(this->delta_X_ini);
#ifdef _USE_MASS_KNOTS_
      this->cwclass.get_Mk_collapsing_regions(this->delta_X_ini,  static_cast<real_prec>(this->N_objects_X)/static_cast<real_prec>(this->NGRID));
#endif

#elif !defined _USE_CWC_
#if defined _USE_MASS_KNOTS_
      this->cwclass.do_CWC(this->delta_X_ini);
      this->cwclass.get_Mk_collapsing_regions(this->delta_X_ini,  static_cast<real_prec>(this->N_objects_X)/static_cast<real_prec>(this->NGRID));
#endif

#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_INVARIANT_TIDAL_FIELD_III_)
      this->cwclass.do_CWC(this->delta_X_ini);
#endif

#endif


#if (defined (_USE_NABLA2DELTA_) || defined (_USE_S2_) || defined (_USE_DELTA2_)  || defined (_USE_S2DELTA_) || defined (_USE_DELTA3_) || defined (_USE_S3_)) && (!defined (_USE_CWC_))
      this->cwclass.get_bias_terms(this->delta_X);
#endif

      if(true==this->Apply_Rankordering)
        {
          So.message_screen("Measuring pdf of log(1+delta) DM");
          calc_pdf("log", this->NGRID,this->NX, this->Xmax, this->Xmin, this->delta_X_ini, pdf_in);
          this->pdf_ini=pdf_in;
#ifdef WRITE_PDF
          string filex=this->Output_directory+"pdf_X_"+this->XNAME+"_MASX"+to_string(this->iMAS_X)+rest_file;
          this->File.write_to_file(filex, xbaux,pdf_in);
#endif

          So.message_screen("Measuring pdf of log(1+delta) DM reference");
          calc_pdf("log", this->NGRID,this->NX, this->Xmax, this->Xmin, this->delta_X_REF_PDF, this->pdf_ref);
#ifdef WRITE_PDF
          filex=this->Output_directory+"pdf_X_REF_"+this->XNAME+"_MASX"+to_string(this->iMAS_X_REF_PDF)+rest_file;
          this->File.write_to_file(filex, xbaux,this->pdf_ref);
#endif

          So.message_screen("Executing rank ordering");
          rankorder(this->step, xbaux, this->NGRID, this->NX,  this->Xmax, this->Xmin, this->delta_X_ini, pdf_in, this->pdf_ref);

          // Get pdf of the rank-ordered and write it to file
#ifdef _VERBOSE_
          So.message_screen("Measuring pdf of delta DM rank-ordered");
#endif
          calc_pdf("log", this->NGRID,this->NX, this->Xmax, this->Xmin, this->delta_X_ini, pdf_in);
#ifdef WRITE_PDF
          filex=this->Output_directory+"pdf_rank_ordered_iteration"+to_string(this->step)+"_X_"+this->XNAME+"_MASX"+to_string(this->iMAS_X)+rest_file;
          this->File.write_to_file(filex, xbaux,pdf_in);
#endif

        }//end if apply rank ordering
    } // end if step==0

  else
    {

      // Compute the kernel
      this->GetKernel(true, 1.0);

      // Convolve the Dm with the kernel and output delta_X
      this->Konvolve(this->delta_X_ini, this->delta_X, "DELTA");   // with this we always convolve the original overdensity field


      // Get the PDF of the convolved field
#ifdef WRITE_PDF
      this->pdf_ite.resize(this->NX, 0);
      So.message_screen("Measuring pdf of log 1+ delta DM convolved");
      calc_pdf("log",this->NGRID,this->NX, this->Xmax, this->Xmin, this->delta_X, this->pdf_ite);

      int index= (this->step <=this->N_iterations_Kernel)  ?  this->step  : this->step - (this->N_iterations_Kernel)+1;
      string label_aux = this->step <= this->N_iterations_Kernel ? "_iteration": "_realization";

      string afilex=this->Output_directory+"pdf_convolved"+label_aux+to_string(index)+"_X_"+this->XNAME+"_MASX"+to_string(this->iMAS_X)+rest_file;
      this->File.write_to_file(afilex, xbaux,this->pdf_ite);
#endif

      if(true==this->Apply_Rankordering)
        {
          /* // Get the pdf of the convolved  */

          So.message_screen("Executing rank ordering");
          rankorder(this->step, xbaux, this->NGRID, this->NX,  this->Xmax, this->Xmin, this->delta_X, pdf_in, this->pdf_ref);
          // Get pdf of the rank-ordered and write it to file

          So.message_screen("Measuring pdf of delta DM rank-ordered");
          calc_pdf("log",this->NGRID,this->NX, this->Xmax, this->Xmin, this->delta_X, pdf_in);

#ifdef WRITE_PDF
          string filex=this->Output_directory+"pdf_rank_ordered_iteration"+to_string(this->step)+"_X_"+this->XNAME+"_MASX"+to_string(this->iMAS_X)+rest_file;
          this->File.write_to_file(filex, xbaux,pdf_in);
#endif

        }


      // Once the new density filed is produced, do again the CWClassification
      //#ifdef _USE_MASS_KNOTS_
      //      this->cwclass.SKNOT_M_info.resize(this->NGRID, 0);
      //#endif

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
#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_INVARIANT_TIDAL_FIELD_III_)
      this->cwclass.do_CWC(this->delta_X_ini);
#endif

#endif

#endif


#if (defined (_USE_NABLA2DELTA_) || defined (_USE_S2_)  || defined (_USE_DELTA2_) || defined (_USE_S2DELTA_) || defined (_USE_DELTA3_) || defined (_USE_S3_) ) && (!defined (_USE_CWC_))
      this->cwclass.get_bias_terms(this->delta_X);
#endif


#ifdef _WRITE_DM_DENSITY_FIELD_
      auto out_it = std::find(std::begin(this->output_at_iteration), std::end(this->output_at_iteration), this->step);
      if (out_it != std::end(this->output_at_iteration))
        this->File.write_array(this->Output_directory+"DM_DELTA_convolved_iteration"+to_string(this->step), this->delta_X);
#endif

      // #ifdef _GET_POWER_REFS_
      //       this->get_power_spectrum("DM_CONV");
      // #endif


    }//end if step>0


  // At the end of this function, delta_X denotes an overdensity

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
  //HEWRE WE HAVE TO USE tracer_ref to read the reference catalog
  So.message_screen("*************************************************************************");
  So.message_screen("*************************************************************************");
#ifdef _USE_VMAX_AS_OBSERVABLE_
  So.message_screen("***Measuring conditional Vmax Function from reference tracer catalog*****");
#else
  So.message_screen("***Measuring conditional Mass Function from reference tracer catalog*****");

#endif
  So.message_screen("****as a function of DM properties of the reference DM density field*****");

  So.message_screen("*************************************************************************");
  cout<<endl;
  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************

  int N_a=N_C_BIN1;
#ifdef _USE_TRACERS_IN_CELLS_
  N_a = N_TRACERS_IN_CELLS_MAX;
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
#ifdef _USE_TOTAL_MASS_IN_CELL_
  N_v = N_BINS_TOTAL_MASS_IN_CELL;
#endif

 // ************************RESIZE VECTORS FOR HISTOGRAMS*************************

  ULONG LENGHT_AB=N_v*N_CV_BIN2*N_CV_BIN3* N_a* N_b* N_c;
  LENGHT_AB*=this->n_sknot_massbin * this->n_cwt * this->n_vknot_massbin * this->n_cwv*this->NX;


  // ********************************************************************************
#ifndef _MASS_ASSIGNMENT_TO_REFERENCE_
  int N_x = this->params._NMASSbins();
  ULONG LENGHT_AB_ONE= LENGHT_AB*N_x;
  So.DONE();
#endif

#ifdef _USE_MASS_AS_OBSERVABLE_
   real_prec lm_min=this->params._LOGMASSmin();
#endif

  // ********************************************************************************
  // ********************************************************************************
  this->get_new_min_max_properties();
  // ********************************************************************************
  // ********************************************************************************

  // Open the ref catalog and read it from this->tracer[].Halo.mass
  cout<<endl;
  So.message_screen("******** IMPORTANT: the reference is currently being read from file newcat.txt, see line", __LINE__);
  this->params.i_coord1_g=0;
  this->params.i_coord2_g=1;
  this->params.i_coord3_g=2;  // esto es improtante, pues estamos leyendo el file con mcut ya hecho
  this->params.i_mass_g=3;
#ifdef _USE_VMAX_TRACERS_
  this->params.i_vmax_g=4;
#endif

  this->tracer_ref.type_of_object="TRACER_REF";
  this->tracer_ref.set_params_catalog(this->params);
  string newcat=this->Output_directory+"newcat.txt";

  //read catalog passing as argument the file and the mininum mass requested
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
#ifdef _SET_CAT_WITH_MASS_CUT_
     this->tracer_ref.read_catalog(newcat,pow(10,params._LOGMASSmin())*params._MASS_units());
#elif defined (_SET_CAT_WITH_VMAX_CUT_)
  this->tracer.read_catalog(newcat,params._VMAXmin());
#endif
#else
#ifdef _SET_CAT_WITH_MASS_CUT_
  this->tracer_ref.read_catalog(newcat,pow(10,params._LOGMASSmin())*params._MASS_units(),static_cast<real_prec>(BIG_NUMBER));
#elif defined (_SET_CAT_WITH_VMAX_CUT_)
  this->tracer_ref.read_catalog(newcat,params._VMAXmin(),static_cast<real_prec>(BIG_NUMBER));
#endif

#endif

//  this->tracer_ref.define_property_bins();
  this->tracer_ref.get_property_function(this->Output_directory+"tracer_ref_abundance.txt");


#if defined _USE_NEIGHBOURS_ || defined _GET_DIST_MIN_SEP_REF_ || defined _GET_DIST_MIN_SEP_MOCK_
  So.message_screen("Identifying Neighbouring Cells (this is done once)");
  this->ncells_info.resize(this->NGRID);
  get_neighbour_cells(this->Nft, N_CELLS_BACK_FORTH, this->ncells_info);
  So.DONE();
  cout<<endl;
#endif

  // ************************ready for power *************************

  Params params_aux=params;
  params_aux.file_catalogue=newcat;
  params_aux.weight_with_mass=false;
  params_aux.i_mass_g=3;
#ifdef _USE_VMAX_TRACERS_
  params_aux.i_vmax_g=4;
#endif
#ifdef _MASS_WEIGHT_POWER_
  params_aux.weight_with_mass=true;
  params_aux.file_power="REF_mass_weight";
  params_aux.SN_correction=false;
#else
  params_aux.file_power="REF";
  params_aux.SN_correction=true;
#endif

  params_aux.file_catalogue=newcat;
  params_aux.Name_survey="UNITSIM_HOSTS_REF";


  PowerSpectrumF Power(params_aux);
  Power.compute_power_spectrum(false,this->tracer_ref.Halo);

  this->Power_REF_MW.resize(this->Nft/2/this->ndel_data,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<this->Nft/2/this->ndel_data;++i)
    this->Power_REF_MW[i]=Power._pk0(i);

  this->kvec.resize(this->Nft/2/this->ndel_data, 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<this->Nft/2/this->ndel_data ;++i)
    this->kvec[i]=Power._kvector_data(i);

  // ********************************************************************************
  // *****************************Deal with reference cats and fields ***************

  // Here we open the reference DM catalog: we have to
  // i ) get delta
  // ii) convolve with kernel from BAM
  // iII) get mins and max
  // if isnide iterative mass assignment, convolve with kernel computed from mass power spectra
  // iv) do CWC classification
  // v) convert to log(num_in_log+delta)

  cout<<endl;

  // ********************************************************************************
  this->cwclass_ref.set_params_cwclass(this->params);
  this->cwclass_ref.s_cosmo_pars=this->s_cosmo_pars;

  // Ideally Read the reference only in the first iteration, all containers with *aux* are not meant to be kept in memmory, they are like dummy

  So.message_screen("*************************************");
  So.message_screen("***********Reading Reference DM******");
  So.message_screen("*************************************");
  this->delta_dm_aux_mem.clear();
  this->delta_dm_aux_mem.shrink_to_fit();
  this->delta_dm_aux_mem.resize(this->NGRID,0);  //Keep this untouched, will be used along the iterations
  File.read_array(this->Output_directory+this->Name_Catalog_X,this->delta_dm_aux_mem);
  this->mean_aux=get_mean(this->delta_dm_aux_mem);
  get_overdens(this->delta_dm_aux_mem,this->mean_aux, this->delta_dm_aux_mem);
  this->Konvolve(this->delta_dm_aux_mem,this->delta_dm_aux_mem,"DELTA");  // COnvolve with BAM Kernel
      // Now delete the enters of th Kernel container, for so far it was the BAm Kernel. It si important to make clear and shrink to fit (resize does not delete enters in the container, it just adds the val specified if the new size is larger than previous)
//  File.write_array(this->Output_directory+"DM_field",dm_ref);

  this->Kernel.clear();
  this->Kernel.shrink_to_fit();

#ifdef _USE_CWC_
      this->cwclass_ref.do_CWC(this->delta_dm_aux_mem);   //get the CW info
#ifdef _USE_MASS_KNOTS_
      this->cwclass_ref.get_Mk_collapsing_regions(this->delta_dm_aux_mem,this->mean_aux);  //get the MK info
#endif //use_mass_knots
#endif

  cout<<endl;
  So.message_screen("Transforming delta_ref -> log10(2+delta_ref). According to preproc, done in line ", __LINE__);
#pragma omp parallel for
  for(ULONG i = 0;i < this->NGRID ;++i)  //TRANSFORM DELTA TO LOG10(NUM_IN_LOG + DELTA)
    this->delta_dm_aux_mem[i] = (this->delta_dm_aux_mem[i]<-1 ? 0  :  log10(NUM_IN_LOG+ this->delta_dm_aux_mem[i]));
  So.DONE();

  // ********************************************************************************
#ifdef _USE_TRACERS_IN_CELLS_
  vector<real_prec> REF_DEN_FIELD(this->NGRID,0);
  this->tracer_ref.get_density_field_grid(_COUNTS_,REF_DEN_FIELD);
  int nmax=get_max<real_prec>(REF_DEN_FIELD);
  So.message_screen("Maximum number of tracer in cell", nmax);
  if(nmax>=N_TRACERS_IN_CELLS_MAX)
    {
      So.message_warning("Please increase the maximum number of tracer in one cell to " ,nmax);
      exit(0);
    }
  cout<<endl;
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
#ifdef _MASS_ASSIGNMENT_TO_REFERENCE_
  this->min_halo_separation=this->tracer_ref.min_halo_separation;
#else
  this->min_halo_separation=MIN_SEP_IN_CELLS;
#endif
  real_prec delta_min_sep = (MAX_SEP_IN_CELLS-this->min_halo_separation)/static_cast<real_prec>(N_BINS_MIN_SEP_IN_CELLS);
#endif


  // ********************************************************************************
  // ********************************************************************************
  // ********************************************************************************

#ifdef _USE_MASS_ASSIGNMENT_READING_REF_MASSES_

#ifdef _USE_VMAX_AS_OBSERVABLE_
  So.message_screen("Allocating", LENGHT_AB* (sizeof(real_prec)+sizeof(ULONG))/(1e9), "Gb for list of Vmax (reference) in theta-bins");
#elif defined _USE_MASS_AS_OBSERVABLE_
  So.message_screen("Allocating", LENGHT_AB* (sizeof(real_prec)+sizeof(ULONG))/(1e9), "Gb for list of masses (reference) in theta-bins");
#endif
  this->dm_properties_bins.clear();
  this->dm_properties_bins.shrink_to_fit();   //vector de estructiras para guarar las masas que caen en un bin de {Theta}
  this->dm_properties_bins.resize(LENGHT_AB);   //vector de estructiras para guarar las masas que caen en un bin de {Theta}
  So.DONE();
#endif


#ifndef _MASS_ASSIGNMENT_TO_REFERENCE_
  So.message_screen("Allocating", LENGHT_AB_ONE* (sizeof(ULONG))/(1e9), "Gb for probability distribution");
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

#ifndef _USE_MASS_ASSIGNMENT_READING_REF_MASSES_
  So.message_screen("**Measuring n(M|theta) from the reference  using ", this->tracer_ref.NOBJS," objects");
#else
#ifdef _USE_MASS_AS_OBSERVABLE_
    So.message_screen("**Identifying reference masses in bins of {theta}_ref using ", this->tracer_ref.NOBJS," objects");
#elif defined _USE_VMAX_AS_OBSERVABLE_
  So.message_screen("**Identifying reference Vmax in bins of {theta}_ref using ", this->tracer_ref.NOBJS," objects");
#endif

#endif

  for(ULONG ig = 0; ig< this->tracer_ref.NOBJS ; ++ig)
    {

#ifdef _USE_MASS_AS_OBSERVABLE_
      real_prec halo_mass=log10(this->tracer_ref.Halo[ig].mass/this->params._MASS_units());
#ifndef _MASS_ASSIGNMENT_TO_REFERENCE_
      int I_Y= get_bin(halo_mass,lm_min,this->params._NMASSbins_mf(),this->tracer_ref.logdeltaM,this->bin_accumulate_borders);
#endif

#elif defined _USE_VMAX_AS_OBSERVABLE_
      real_prec halo_mass=log10(this->tracer_ref.Halo[ig].vmax);

#ifndef _MASS_ASSIGNMENT_TO_REFERENCE_
      int I_Y =get_bin(halo_mass, log10(this->params._VMAXmin()),N_x,log10(this->params._VMAXmax()/this->params._VMAXmin())/static_cast<double>(N_x), true);
#endif

#endif


      // Get ID of this tracer
      ULONG ID=this->tracer_ref.Halo[ig].GridID;

      // Get bin of DM ovedensity
      real_prec xdm  = this->delta_dm_aux_mem[ID];
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
#ifdef _USE_TRACERS_IN_CELLS_
      I_C1=static_cast<int>(REF_DEN_FIELD[ID]);
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
#elif defined _USE_TIDAL_ANISOTROPY_
      real_prec C3 = this->cwclass_ref.Tidal_Anisotropy[ID];
      I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined _USE_S2_
      real_prec C3 = this->cwclass_ref.S2[ID];             // s²
      I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#endif

      int I_CV1=0;
#ifdef _USE_TOTAL_MASS_IN_CELL_
      I_CV1=get_bin( log10(REF_MASS_FIELD[ID]), this->params._LOGMASSmin(),N_BINS_TOTAL_MASS_IN_CELL,(log10(MAX_TOTAL_MASS_IN_CELL)-this->params._LOGMASSmin())/static_cast<double>(N_BINS_TOTAL_MASS_IN_CELL), true);
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
      mean_hmass+=halo_mass;
      m_dm+=halo_mass*xdm;
#ifdef _USE_TRACERS_IN_CELLS_
      mean_ntr+=REF_DEN_FIELD[ID];
      m_ntr+=halo_mass*REF_DEN_FIELD[ID];
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
#ifdef _USE_TIDAL_ANISOTROPY_
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
                        // This will be used for the mocks at X<Xthreshold
#ifndef  _MASS_ASSIGNMENT_TO_REFERENCE_
                      ULONG index_prop_ab= index_ant + LENGHT_AB * I_Y; // index_2d(I_Y,index_ant,LENGHT_AB); I use this form because the argumenrs are ULONG, not int as the index_xd expects
                      this->ABUNDANCE[index_prop_ab]++;
#endif

                      // This will be used for mocks and assign_to_ref for X>=Xthres level 1 only
#ifdef _USE_MASS_ASSIGNMENT_READING_REF_MASSES_
#ifdef _USE_MASS_AS_OBSERVABLE_
                      this->dm_properties_bins[index_ant].masses_bin_properties.push_back(this->tracer_ref.Halo[ig].mass);
#elif defined _USE_VMAX_AS_OBSERVABLE_
                      this->dm_properties_bins[index_ant].masses_bin_properties.push_back(this->tracer_ref.Halo[ig].vmax);
#endif
                      this->dm_properties_bins[index_ant].used_mass.push_back(false);
#endif



#ifndef _BIN_ACCUMULATE_
                      }
#endif

    }


  this->So.DONE();

  // So.message_warning("trafer_ref.Halo.clear() is disabled in order to compute correlation between masses, line", __LINE__);
 // this->tracer_ref.Halo.clear();
 // this->tracer_ref.Halo.shrink_to_fit();


#ifdef _GET_HALO_MASS_PROPERTIES_CORRELATIONS_
  mean_ldm/=static_cast<real_prec>(this->tracer_ref.NOBJS);
  mean_hmass/=static_cast<real_prec>(this->tracer_ref.NOBJS);
  m_dm/=static_cast<real_prec>(this->tracer_ref.NOBJS);
#ifdef _USE_MASS_AS_OBSERVABLE_
  So.message_screen("Correlation log(M) - DM density =", sqrt(fabs(m_dm-mean_hmass*mean_ldm)));
#elif defined _USE_VMAX_AS_OBSERVABLE_
  So.message_screen("Correlation log(Vmax) - DM density =", sqrt(fabs(m_dm-mean_hmass*mean_ldm)));
#endif

#ifdef _USE_TRACERS_IN_CELLS_
  mean_ntr/=static_cast<real_prec>(this->tracer_ref.NOBJS);
  m_ntr/=static_cast<real_prec>(this->tracer_ref.NOBJS);
#ifdef _USE_MASS_AS_OBSERVABLE_

  So.message_screen("Halo-mass - N_tracers in cells  =", sqrt(fabs(m_ntr-mean_hmass*mean_ntr)));
#elif defined _USE_VMAX_AS_OBSERVABLE_
  So.message_screen("Halo-Vmax - N_tracers in cells =", sqrt(fabs(m_ntr-mean_hmass*mean_ldm)));
#endif
#endif
#endif
  cout<<endl;




#ifdef _USE_TRACERS_IN_CELLS_
  REF_DEN_FIELD.clear();
  REF_DEN_FIELD.shrink_to_fit();
#endif



#ifdef _USE_MASS_ASSIGNMENT_READING_REF_MASSES_
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
      So.message_warning("The number of tracers counted in ABUNDANCE(Y,X) is smaller than the input. Perhaps the range for some of the DMF properties *must* be increased to include all cells. CosmiCatlas stops here. Line", __LINE__);
      cout<<aux_a<<"  "<<this->tracer_ref.NOBJS<<endl;
      exit(0);
    }
  So.DONE();
#endif




#ifndef _MASS_ASSIGNMENT_TO_REFERENCE_
  So.message_screen("Normalizing probability distribution");
  this->ABUNDANCE_normalized.clear();
  this->ABUNDANCE_normalized.shrink_to_fit();
  this->ABUNDANCE_normalized.resize(LENGHT_AB_ONE,0);
 #ifdef _USE_OMP_
 #pragma omp parallel for
 #endif
  for(int ULONG idm = 0; idm < LENGHT_AB; ++idm)
  {
      long aux_a=-10000;
      long aux_b;

      for(int tr_j = 0; tr_j < N_x; ++tr_j)
        {
//          ULONG indexa=index_2d(tr_j,idm,LENGHT_AB);
          long AUX=this->ABUNDANCE[idm+tr_j*LENGHT_AB];
          aux_b=max(aux_a, AUX);
          aux_a=aux_b;
        }
      for(int tr_j = 0; tr_j < N_x; ++tr_j)
        {
          ULONG indexa=idm+tr_j*LENGHT_AB;//index_2d(tr_j,idm,LENGHT_AB);
          long AUX=this->ABUNDANCE[indexa];
          this->ABUNDANCE_normalized[indexa]=(aux_b==0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
        }
    }

  So.DONE();
  So.message_screen("Freeing memory");
  this->ABUNDANCE.clear();
  this->ABUNDANCE.shrink_to_fit();
  So.DONE();

#endif




#ifndef _USE_MASS_ASSIGNMENT_READING_REF_MASSES_
  //* Final check to verify that all tracers have been counted in
  ULONG ncells_used = static_cast<ULONG>(get_nobjects(this->ABUNDANCE));

  if(ncells_used<this->tracer_ref.NOBJS)
    {
      So.message_warning("The number of tracers counted in ABUNDANCE(Y,X) is smaller than the input. Perhaps the range for some of the DMF properties *must* be increased to include all cells. CosmiCatlas stops here. Line", __LINE__);
      cout<<ncells_used<<"  "<<this->tracer_ref.NOBJS<<endl;
    }


   So.message_screen("Assigning memmory space for conditional probability function");
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
#ifndef _USE_MASS_ASSIGNMENT_READING_REF_MASSES_
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

#ifndef _USE_MASS_ASSIGNMENT_READING_REF_MASSES_
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


    So.message_screen("Check on number of bins used...");
    real_prec aux_h=0;

#pragma omp parallel for reduction(+:aux_h)
    for(int ih=0;ih<this->ABUNDANCE_normalized.size();++ih)
      aux_h+=this->ABUNDANCE_normalized[ih];
    if(aux_h<=0)
      {
        So.message_warning("Normalized Joint probability for number counts ill-defined. CosmicAtlas stops here. Line ",__LINE__);
        exit(0);
      }
    else
      this->So.DONE();


    So.message_screen("Freeing memmory");
    this->ABUNDANCE.clear();
    this->ABUNDANCE.shrink_to_fit();
    So.DONE();

    // We do not need her to write these files as they will be inmediatly used in the assign_tracer_mass
    //  this->File.write_array(this->Output_directory+"Bam_Abundance", this->ABUNDANCE);
    // this->File.write_array(this->Output_directory+"Bam_Abundance_Normalized", this->ABUNDANCE_normalized);

#endif



}

//end class member function


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
void Bam::get_X_function_complement()
{
  So.message_screen("***Measuring conditional Mass Function from reference tracer catalog*****");
  So.message_screen("****as a function of DM properties of the reference DM density field*****");



#ifdef test_vmax_mass
  So.message_warning("The pre-proc directive test_vmax_mass is defined");
  So.message_warning("This implies that, once Vmax has been assigned (with all DM properties), BAM assign masses using only the P(M|Vmax) scaling relation");
#endif

  cout<<endl;
#ifndef test_vmax_mass
  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************

  int N_a=N_C_BIN1;

  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************
  int N_b=N_C_BIN2;

  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************
  int N_c = N_C_BIN3; // We shall not use min_sep for post-mass assignment
  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************

  int N_v= N_CV_BIN1;
  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************
#endif

  int N_x= N_CV_BIN2;
#ifdef _ASSIGN_MASS_POST_
  N_x = N_VMAX_BINS;
#endif

  ULONG Ntot=N_x;
#ifdef _add_dm_density_
 Ntot*=this->NX;
#endif


#ifndef test_vmax_mass
  ULONG LENGHT_AB=N_v*N_x*N_CV_BIN3* N_a* N_b* N_c;
  LENGHT_AB*=this->n_sknot_massbin * this->n_cwt * this->n_vknot_massbin * this->n_cwv*this->NX;
#else
  ULONG LENGHT_AB=Ntot;
#endif
  // ********************************************************************************
  this->ABUNDANCE.clear();
  this->ABUNDANCE.shrink_to_fit();
  ULONG LENGHT_AC=LENGHT_AB*this->params._NMASSbins();
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
#ifdef _USE_MASS_ASSIGNMENT_READING_REF_MASSES_
  So.message_screen("Allocating", LENGHT_AB* (sizeof(real_prec)+sizeof(ULONG))/(1e9), "Gb for list of masses in theta-bins");
  this->dm_properties_bins.clear();
  this->dm_properties_bins.shrink_to_fit();   //vector de estructiras para guarar las masas que caen en un bin de {Theta}
  this->dm_properties_bins.resize(LENGHT_AB);   //vector de estructiras para guarar las masas que caen en un bin de {Theta}
  So.DONE();
#endif
#endif

  // ********************************************************************************




#ifndef _USE_MASS_ASSIGNMENT_READING_REF_MASSES_
  So.message_screen("**Measuring n(M|theta) from the reference  using ", this->tracer_ref.NOBJS," objects");
#else
  So.message_screen("**Identifying reference Masses in bins of {theta}_ref using ", this->tracer_ref.NOBJS," objects");
#endif


  for(ULONG ig = 0; ig< this->tracer_ref.NOBJS ; ++ig)
    {

      // Get ID of this tracer
#if !defined test_vmax_mass || defined (_add_dm_density_)
      ULONG ID=this->tracer_ref.Halo[ig].GridID;
#endif

      int I_Y  = get_bin(log10(this->tracer_ref.Halo[ig].mass),this->params._LOGMASSmin(),this->params._NMASSbins(),this->tracer_ref.logdeltaM,this->bin_accumulate_borders);

      // Get bin of DM ovedensity
#if  !defined (test_vmax_mass) || defined (_add_dm_density_)
      int I_X=0;
      real_prec xdm  = this->delta_dm_aux_mem[ID]; // this and its CWC was cmputed in get_X_function()
      I_X  = get_bin(xdm,this->s_mins.prop1,this->NX,this->s_deltas.prop1,this->bin_accumulate_borders);
#endif


#ifndef test_vmax_mass
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
      I_CWV=this->cwclass_ref.get_Vclassification(ID);
#endif

      int I_VK=0;
#ifdef _USE_VEL_KNOTS_V_
      I_VK= (this->cwclass_ref.cwv_used[I_CWV]== I_KNOT ? cwc_ref.VDISP_KNOT_info[ID]: 0);
#endif


      int I_C1=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
      real_prec C1 = this->cwclass_ref.Invariant_TF_II[ID];
      I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#elif defined _USE_DELTA2_
      real_prec C1 = this->cwclass_ref.DELTA2[ID];
      I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#endif

      int I_C2=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
      real_prec C2 = cwclass_ref.Invariant_TF_III[ID];
      I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#elif defined _USE_DELTA3_
      real_prec C2 = cwclass_ref.DELTA3[ID];
      I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#endif

      int I_C3=0;
#ifdef _USE_TIDAL_ANISOTROPY_
      real_prec C3 = cwclass_ref.Tidal_Anisotropy[ID];
      I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined _USE_S2_
      real_prec C3 = cwclass_ref.S2[ID];             // s²
      I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#endif



      int I_CV1=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
      real_prec CV1 = invariant_field_I(cwclass_ref.lambda1_vs[ID],cwc_ref.lambda2_vs[ID],cwc_ref.lambda3_vs[ID]); // Not yet assigned to container
      I_CV1= get_bin(CV1, CV1_MIN, N_CV_BIN1,DELTA_CV1, this->bin_accumulate_borders);
#elif defined _USE_NABLA2DELTA_
      real_prec CV1 = this->cwclass_ref.N2D[ID];      // Nabla² ð
      I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1, s_deltas.prop7,this->bin_accumulate_borders);
#endif
#endif

      int I_CV2=0;
#ifdef _ASSIGN_MASS_POST_
      I_CV2=get_bin(log10(this->tracer_ref.Halo[ig].vmax), log10(this->params._VMAXmin()),N_x,log10(this->params._VMAXmax()/this->params._VMAXmin())/static_cast<double>(N_x), true);
#elif defined (_USE_INVARIANT_SHEAR_VFIELD_II_)
      real_prec CV2 = invariant_field_II(cwclass_ref.lambda1_vs[ID],cwc_ref.lambda2_vs[ID],cwc_ref.lambda3_vs[ID]);// Not yet assigned to container
      I_CV2= get_bin(CV2, CV2_MIN, N_CV_BIN2,DELTA_CV2, this->bin_accumulate_borders);
#elif defined _USE_S2DELTA_
      real_prec CV2 = this->cwclass_ref.S2DELTA[ID];         // s²ð
      I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2, s_deltas.prop8,this->bin_accumulate_borders);
#endif


#ifndef test_vmax_mass
      int I_CV3=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
      real_prec CV3 = invariant_field_III(cwclass_ref.lambda1_vs[ID],cwc_ref.lambda2_vs[ID],cwc_ref.lambda3_vs[ID]);// Not yet assigned to container
      I_CV3= get_bin(CV3, CV3_MIN, N_CV_BIN3,DELTA_CV3, this->bin_accumulate_borders);
#elif defined _USE_S3_
      real_prec CV3 = this->cwclass_ref.S3[ID];                                   // s³
      I_CV3= get_bin(CV3, this->s_mins.prop9, N_CV_BIN3, s_deltas.prop9,this->bin_accumulate_borders);
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
#ifdef _USE_TIDAL_ANISOTROPY_
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

#ifndef test_vmax_mass
                      ULONG index_prop=index_11d(I_X, I_CWT, I_MK, I_CWV,I_VK, I_C1, I_C2, I_C3, I_CV1, I_CV2, I_CV3,this->n_cwt, this->n_sknot_massbin, this->n_cwv, this->n_vknot_massbin,N_a,N_b,N_c, N_v, N_x,N_CV_BIN3);
                      this->dm_properties_bins[index_prop].masses_bin_properties.push_back(this->tracer_ref.Halo[ig].mass);
                      this->dm_properties_bins[index_prop].used_mass.push_back(false);
                      ULONG index_prop_b=index_2d(I_Y,index_prop,LENGHT_AB);
#endif

#ifdef test_vmax_mass
                      ULONG index_dm;
#ifdef _add_dm_density_
                      index_dm=index_2d(I_CV2,I_X,this->NX);
#else
                      index_dm=I_CV2;
#endif
                      ULONG index_prop_b=index_2d(I_Y,index_dm,Ntot);
#endif
                      this->ABUNDANCE[index_prop_b]++;
#ifndef _BIN_ACCUMULATE_
                      }
#endif
    }


 this->So.DONE();

 So.message_screen("**Freeing memmory from tracer_ref and auxiliary containers, line", __LINE__);
 this->delta_dm_aux_mem.clear();delta_dm_aux_mem.shrink_to_fit();
 this->tracer_ref.Halo.clear();
 this->tracer_ref.Halo.shrink_to_fit();
 this->So.DONE();


 So.message_screen("Normalizing");
 this->ABUNDANCE_normalized.clear();
 this->ABUNDANCE_normalized.shrink_to_fit();
 this->ABUNDANCE_normalized.resize(LENGHT_AC,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
 for(int ULONG idm = 0; idm < LENGHT_AB; ++idm)
 {
     long aux_a=-10000;
     long aux_b;

     for(int tr_j = 0; tr_j < this->params._NMASSbins_mf(); ++tr_j)
       {
         ULONG indexa=index_2d(tr_j,idm,LENGHT_AB);
         long AUX=this->ABUNDANCE[indexa];
         aux_b=max(aux_a, AUX);
         aux_a=aux_b;
       }
     for(int tr_j = 0; tr_j < this->params._NMASSbins_mf(); ++tr_j)
       {
         ULONG indexa=index_2d(tr_j,idm,LENGHT_AB);
         long AUX=this->ABUNDANCE[indexa];
         this->ABUNDANCE_normalized[indexa]=(aux_b==0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
       }
   }

 So.DONE();

 this->ABUNDANCE.clear();
 this->ABUNDANCE.shrink_to_fit();

 cout<<endl;


#ifndef test_vmax_mass
#ifdef _USE_MASS_ASSIGNMENT_READING_REF_MASSES_
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
      So.message_warning("The number of tracers counted in ABUNDANCE(Y,X) is smaller than the input. Perhaps the range for some of the DMF properties *must* be increased to include all cells. CosmiCatlas stops here. Line", __LINE__);
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


   //IOn tbis function the TR and

   real_prec num_in_log_x = true==this->Convert_Density_to_Delta_X ? NUM_IN_LOG: 0.;
   real_prec num_in_log_y = true==this->Convert_Density_to_Delta_Y ? NUM_IN_LOG: 0.;

#ifdef _DO_BAM_CALIBRATION_
   cout<<endl;

   if(this->iteration_ini==0)
     So.message_screen("Measuring Bias");


   if(this->Scale_Y=="log")

       if(this->step==0) // Do this only in the first step, sicne in the iteration process does not change the tracer density field.
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
        for(ULONG i = 0; i < this->NGRID; ++i)
          this->delta_Y[i] = log10(num_in_log_y+static_cast<real_prec>(this->delta_Y[i]));


   if(this->Scale_X=="log")
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i = 0;i < this->NGRID ;++i)  //TRANSFORM DELTA TO LOG10(NUM_IN_LOG + DELTA)
       this->delta_X[i] = this->delta_X[i]<-1 ?  0 : log10(num_in_log_x + static_cast<real_prec>(this->delta_X[i]));
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
   ULONG LENGHT_BIAS_NCOUNT=N_CV_BIN1 * N_CV_BIN2* N_CV_BIN3* N_C_BIN1 * N_C_BIN2* N_C_BIN3* this->n_sknot_massbin * this->n_cwt * this->n_vknot_massbin * this->n_cwv* NX * new_nbins_y;

#ifdef _USE_SAT_FRACTION_
   ULONG LENGHT_BIAS_SAT=N_CV_BIN1 *1 N_CV_BIN2* N_CV_BIN3* N_C_BIN1 * N_C_BIN2* N_C_BIN3* this->n_sknot_massbin * this->n_cwt * this->n_vknot_massbin * this->n_cwv* NX * this->NY_SAT_FRAC;
#endif


   int count_arrays=1;
#ifdef _USE_MASS_TRACERS_
   count_arrays++;
#endif
#ifdef _USE_SAT_FRACTION_
   count_arrays++;
#endif



#if !defined (_MASS_ASSIGNMENT_TO_REFERENCE_) || defined (_DO_BAM_CALIBRATION_)
   So.message_warning("Recall to change DM_NEW from that of the reference to the new DM field with which the used catalog has been constructed");
   So.message_screen("Allocating", LENGHT_BIAS_NCOUNT*count_arrays*(sizeof(real_prec)+sizeof(ULONG))/(1e9), "Gb for bias");

   this->BIAS_NCOUNTS.clear();
   this->BIAS_NCOUNTS.shrink_to_fit();
   this->BIAS_NCOUNTS.resize(LENGHT_BIAS_NCOUNT, 0);
#endif


#ifdef _USE_SAT_FRACTION_
   this->BIAS_SAT_FRACTION.clear();
   this->BIAS_SAT_FRACTION.shrink_to_fit();
   this->BIAS_SAT_FRACTION.resize(LENGHT_BIAS_SAT, 0);
   this->BIAS_SAT_FRACTION_normalized.shrink_to_fit();
   this->BIAS_SAT_FRACTION_normalized.clear();
   this->BIAS_SAT_FRACTION_normalized.resize(LENGHT_BIAS_SAT, 0);
#endif

   // ******************************************************************************

#ifdef _GET_BAM_REALIZATIONS_

   So.message_screen("Reading BIAS and Kernel from calibration");
   this->Kernel.clear();
   this->Kernel.shrink_to_fit();
   this->Kernel.resize(this->NTT, 0.0);
   this->File.read_array(this->Output_directory+"Bam_Kernel.dat", this->Kernel);

#ifndef _MASS_ASSIGNMENT_TO_REFERENCE_
   this->File.read_array(this->Output_directory+"Bam_Bias.dat", this->BIAS_NCOUNTS);
   ULONG ncounts_aux=0;
#pragma parallel for reduction(+:ncounts_aux)
   for(ULONG i=0;i <this->BIAS_NCOUNTS.size(); ++i)
     ncounts_aux+=this->BIAS_NCOUNTS[i];
   So.message_screen("Number of cells BIAS = ", ncounts_aux);

   this->File.read_array(this->Output_directory+"Bam_Bias_Normalized.dat", this->BIAS_NCOUNTS_normalized);



#ifdef _USE_SAT_FRACTION_
   this->File.read_array(this->Output_directory+"Bam_Bias_SatFraction.dat", this->BIAS_SAT_FRACTION);
   this->File.read_array(this->Output_directory+"Bam_Bias_SatFraction_Normalized.dat", this->BIAS_SAT_FRACTION_normalized);

#endif
#endif

#else

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
           real_prec halo_prop = this->delta_Y[ig];


           int I_Y;
           if(_COUNTS_== property)
             I_Y=static_cast<int>(halo_prop);
           else
             I_Y=get_bin(halo_prop,this->s_mins.prop0, this->new_nbins_y, this->s_deltas.prop0,this->bin_accumulate_borders);

#ifdef _USE_MASS_FIELD_   //get bin of tracer mass
           real_prec halo_mass_prop=this->delta_Y_MASS[ig];
           int I_Y_MASS=get_bin(halo_mass_prop,this->s_mins.prop0_mass, this->NY_MASS, this->s_deltas.prop0_mass,this->bin_accumulate_borders);
#endif


#ifdef _USE_SAT_FRACTION_   //get # of  satellites, that is why we leave it as it comes instead of computing the histogram
           int I_Y_SFRAC=static_cast<int>(this->delta_Y_SAT_FRAC[ig]);

           //get_bin(this->delta_Y_SAT_FRAC[ig],this->s_mins.prop0_sf, this->NY_SAT_FRAC, this->s_deltas.prop0_sf,this->bin_accumulate_borders);
#endif


           real_prec xdm    = this->delta_X[ig];
           int I_X  = ((this->iMAS_X == 0  && false==this->Convert_Density_to_Delta_X) ? static_cast<int>(this->delta_X[ig]) : get_bin(xdm,this->s_mins.prop1,this->NX,this->s_deltas.prop1,this->bin_accumulate_borders));



           int I_CWT=0;
#ifdef _USE_CWC_
           I_CWT=this->cwclass.get_Tclassification(ig);
#endif

           int I_MK=0;
#ifdef _USE_MASS_KNOTS_
           I_MK=this->cwclass.SKNOT_M_info[ig];
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
#ifdef _USE_TIDAL_ANISOTROPY_
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
#ifdef _USE_TIDAL_ANISOTROPY_
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
			     ULONG Index = index_12d(I_Y,I_X, I_CWT, I_MK, I_CWV,I_VK, I_C1, I_C2, I_C3, I_CV1, I_CV2, I_CV3, this->NX, this->n_cwt, this->n_sknot_massbin, this->n_cwv, this->n_vknot_massbin,N_C_BIN1,N_C_BIN2,N_C_BIN3, N_CV_BIN1, N_CV_BIN2,N_CV_BIN3);
#pragma omp atomic update
			     this->BIAS_NCOUNTS[Index]++;
			     
			     
#ifdef _USE_SAT_FRACTION_
#pragma omp atomic update
			     this->BIAS_SAT_FRACTION[index_12d(I_Y_SFRAC,I_X, I_CWT, I_MK,I_CWV,I_VK, I_C1, I_C2, I_C3, I_CV1, I_CV2, I_CV3, this->NX, this->n_cwt, this->n_sknot_massbin, this->n_cwv, this->n_vknot_massbin,N_C_BIN1,N_C_BIN2,N_C_BIN3, N_CV_BIN1, N_CV_BIN2,N_CV_BIN3)]++;
#endif
			   }
         }
       this->So.DONE();
       
       
       //* Final check done only in case we use counts-in-cells
       // This is only done when we gemnerate mocks with the same volume of the reference
       So.message_screen("Check: number of cells in BIAS");

#ifndef _EXTRAPOLATE_VOLUME_
       ULONG ncells_used = static_cast<ULONG>(get_nobjects(this->BIAS_NCOUNTS));
       if(ncells_used<this->NGRID)
           So.message_warning("The number of cells counted in BIAS(Y,X) is smaller than NGRID. Perhaps the range for some of the DMF properties *must* be increased to include all cells. CosmiCatlas stops here. Line", __LINE__);
        else
           this->So.DONE();
#endif
       
       
       
       So.message_screen("Normalizing:");
       //     int nbins_y_temp = property==_COUNTS_ ?  this->new_nbins_y: this->new_nbins_y_MW;
       int nbins_y_temp = this->new_nbins_y;
       

#if !defined (_MASS_ASSIGNMENT_TO_REFERENCE_) || defined (_DO_BAM_CALIBRATION_)
       this->BIAS_NCOUNTS_normalized.shrink_to_fit();
       this->BIAS_NCOUNTS_normalized.clear();
       this->BIAS_NCOUNTS_normalized.resize(LENGHT_BIAS_NCOUNT, 0);
#endif
              
#ifdef _USE_OMP_
#pragma omp parallel for collapse(11)
#endif
   for(int i=0;i < this->NX;++i)// loop sobre bines de dark matter
     for(int sua = 0; sua < this->n_cwt; ++sua)
       for(int vua = 0; vua < this->n_cwv; ++vua)
	 for(int k = 0; k < this->n_sknot_massbin ;  ++k)
           for(int kv = 0; kv < this->n_vknot_massbin ;  ++kv)
	     for(int l1 = 0; l1 < N_C_BIN1; ++l1)
	       for(int l2 = 0; l2 < N_C_BIN2; ++l2)
		 for(int l3 = 0; l3 < N_C_BIN3; ++l3)
		   for(int lv1 = 0; lv1 < N_CV_BIN1; ++lv1)
		     for(int lv2 = 0; lv2 < N_CV_BIN2; ++lv2)
		       for(int lv3 = 0; lv3 < N_CV_BIN3; ++lv3)
			 {
                           long aux_a=-10000;
                           long aux_b;
                           for(int tr_j = 0; tr_j < nbins_y_temp; ++tr_j)
                            {
                              ULONG inde=index_12d(tr_j,i,sua,k,vua,kv,l1,l2,l3,lv1, lv2, lv3,this->NX,this->n_cwt, this->n_sknot_massbin, this->n_cwv, this->n_vknot_massbin,N_C_BIN1,N_C_BIN2,N_C_BIN3, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
                              long AUX=this->BIAS_NCOUNTS[inde];
                              aux_b=max(aux_a, AUX);
                              aux_a=aux_b;
                            }
                           for(int tr_j = 0; tr_j < nbins_y_temp; ++tr_j)
                           {
                             ULONG inde=index_12d(tr_j,i,sua,k,vua,kv,l1,l2,l3,lv1, lv2, lv3,this->NX,this->n_cwt, this->n_sknot_massbin, this->n_cwv, this->n_vknot_massbin,N_C_BIN1,N_C_BIN2,N_C_BIN3, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
                             long AUX=this->BIAS_NCOUNTS[inde];
                             this->BIAS_NCOUNTS_normalized[inde] = (aux_a<=0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
                           }

#ifdef _USE_SAT_FRACTION_
			   AUX.clear();AUX.shrink_to_fit();
			   AUX.resize(this->NY_SAT_FRAC,0);
			   for(int tr_j = 0; tr_j < this->NY_SAT_FRAC; ++tr_j)
			     AUX[tr_j]=this->BIAS_SAT_FRACTION[index_12d(tr_j,i,sua,k,l1,l2,l3,lv1, lv2, lv3,this->NX,this->n_cwt, this->n_sknot_massbin, this->n_cwv, this->n_vknot_massbin,N_C_BIN1,N_C_BIN2,N_C_BIN3, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3)];
			   lkk=get_max<long>(AUX);
			   for(int tr_j = 0; tr_j < this->NY_SAT_FRAC; ++tr_j)
			     this->BIAS_SAT_FRACTION_normalized[index_12d(tr_j,i,sua,k,l1,l2,l3,lv1, lv2, lv3, this->NX,this->n_cwt, this->n_sknot_massbin,this->n_cwv, this->n_vknot_massbin,N_C_BIN1,N_C_BIN2,N_C_BIN3, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3)] = (lkk==0 ? 0. : static_cast<real_prec>(AUX[tr_j])/static_cast<real_prec>(lkk));
#endif
			   
			 }


   
   this->So.DONE();

   real_prec aux_h=0;
#pragma omp parallel for reduction(+:aux_h)
   for(int ih=0;ih<this->BIAS_NCOUNTS_normalized.size();++ih)
     aux_h+=this->BIAS_NCOUNTS_normalized[ih];



   if(aux_h<=0)
     {
       cout<<"Sum of conditional probability = "<<aux_h<<endl;
       So.message_warning("Normalized Joint probability for number counts is ill-defined. CosmicAtlas stops here (Bam.cpp), line ",__LINE__);
       exit(0);
     }
   
   
   
#ifdef _USE_SAT_FRACTION_
   aux_h=0;
   for(int ih=0;ih<this->BIAS_SAT_FRACTION.size();++ih)
     aux_h+=this->BIAS_SAT_FRACTION_normalized[ih];
   if(aux_h<=0)
     {
       So.message_warning("Normalized Joint probability for sat_fraction ill-defined. CosmicAtlas stops here. Line ",__LINE__);
       exit(0);
     }
#endif
   
   
   
   
#ifndef _TEST_THRESHOLDS_RESIDUALS_
   auto out_it = std::find(std::begin(this->output_at_iteration), std::end(this->output_at_iteration), this->step);
   if( (this->step==this->N_iterations_Kernel) || (out_it != std::end(this->output_at_iteration)))
     {
       this->File.write_array(this->Output_directory+"Bam_Bias", this->BIAS_NCOUNTS);
       this->File.write_array(this->Output_directory+"Bam_Bias_Normalized", this->BIAS_NCOUNTS_normalized);
       
#ifdef _USE_SAT_FRACTION_
       this->File.write_array(this->Output_directory+"Bam_Bias_SatFraction", this->BIAS_SAT_FRACTION);
       this->File.write_array(this->Output_directory+"Bam_Bias_SatFraction_Normalized", this->BIAS_SAT_FRACTION_normalized);
#endif
       
     }
#endif
   
     }
#endif // endif _GET_BAM_REALIZARTIONS_
   
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

 // This function is used in mode _GET_BAM_REALIZATIONS_ and gets the approximated DMF , convolves it with the BAM kernel
 void Bam::get_new_DM_field()
 {

   cout<<endl;
   this->So.message_screen("*************************************************************************");
   this->So.message_screen("*************************************************************************");
   this->So.message_screen("****Getting new DM density field in line", __LINE__, "of file",__FILE__);
   this->So.message_screen("*************************************************************************");
   this->So.message_screen("*************************************************************************");
   cout<<endl;

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

   // #ifndef _USE_PATCHY_  //If we cannot run patchy, we need to read DM files from somewhere, so we read them from paths hard coded here
   //       string dir_dm;
   //       this->So.message("Reading DM fields from ALPT+IC Minerva in /scratch/marcos/DMforBAM/");
   //       if(n_realization<10)
   // 	dir_dm= "/scratch/marcos/DMforBAM/DM_dens_ALPT_000"+to_string(n_realization)+"/";
   //       else if(n_realization>=10 && n_realization<100)
   // 	dir_dm= "/scratch/marcos/DMforBAM/DM_dens_ALPT_00"+to_string(n_realization)+"/";
   //       else
   // 	dir_dm= "/scratch/marcos/DMforBAM/DM_dens_ALPT_0"+to_string(n_realization)+"/";

   //       this->File.read_array(dir_dm+this->Name_Catalog_X_NEW, this->delta_X_ini);



#ifndef _USE_PATCHY_  //If we cannot run patchy, we need to read DM files from somewhere, so we read them from paths hard coded here
   this->File.read_array(this->Output_directory+this->Name_Catalog_X_NEW, this->delta_X_ini);



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
   So.message_screen("Total number of X objects (new) =", nmean);
   nmean/=static_cast<real_prec>(this->NGRID);
   So.message_screen("Average number of X objects (new) =", nmean);
   cout<<endl;

#ifndef _MASS_ASSIGNMENT_TO_REFERENCE_   // just to save tmie while doing the test of mass assingment to tracers

#ifndef _MASS_WEIGHT_POWER_
   this->get_power_spectrum("DM_REF_NEW");
#endif
#endif

   get_overdens(this->delta_X_ini, this->delta_X_ini);


#else
   {

     const gsl_rng_type *rng_t;
     gsl_rng **gBaseRand;
     gsl_rng_env_setup();
     rng_t = gsl_rng_mt19937;// gsl_rng_default;
     int nt=omp_get_max_threads();
     gBaseRand = (gsl_rng **) malloc(nt * sizeof(gsl_rng *));

#pragma omp parallel for num_threads(nt)
     for(int i=0;i<nt;i++)
       {
         gBaseRand[i] = gsl_rng_alloc(rng_t);
         gsl_rng_set(gBaseRand[i],this->seed);
       }

     /*
     //Construct strings with names for output fiels in patchy.
     //this->patchy.set_fnames();// we did this before. If needed, update names by accessing the public variables

     string dir_ic; //path to initial conditions (White Noise)
     string dir_icdelta; //path to delta (innitial delta with input power)
     if(n_realization<10)
     {
     dir_ic= this->params._dir()+this->params._ic_WN_dir()+"c_renorm_000"+to_string(n_realization)+"_500_delta";
     dir_icdelta=this->params._dir()+this->params._ic_WN_dir()+"c_renorm_000"+to_string(n_realization)+"_500_deltaICField";
     }
     else if(n_realization>=10 && n_realization<100)
     {
     dir_ic=  this->params._dir()+this->params._ic_WN_dir()+"c_renorm_00"+to_string(n_realization)+"_500_delta";
     dir_icdelta=this->params._dir()+this->params._ic_WN_dir()+"c_renorm_00"+to_string(n_realization)+"_500_deltaICField";
     }
     else
     {
     dir_ic=  this->params._dir()+this->params._ic_WN_dir()+"c_renorm_0"+to_string(n_realization)+"_500_delta";
     dir_icdelta=this->params._dir()+this->params._ic_WN_dir()+"c_renorm_0"+to_string(n_realization)+"_500_deltaICField";
     }


     // Set the new name of the intial conditions
     this->patchy.fnameIC=dir_ic;

     // Set the new name of the DM filed generated by Patchy
     this->patchy.fnameDM += "realization"+to_string(n_realization);

     this->patchy.fnameICDELTA=dir_icdelta;


     */


     time_t start_patchy;
     time(&start_patchy);

     if(false==dm_already_done) //COmpute DM if it has not been done before through a call of this function
       {
#ifdef OMPPARRAN
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


     //#ifdef _GET_POWER_REFS_
     this->get_power_spectrum("DM_REF_NEW");
     //#endif
     nmean=get_nobjects(this->delta_X_ini);
     So.message_screen("Total number of X objects (new) =", nmean);
     nmean/=static_cast<real_prec>(this->NGRID);
     So.message_screen("Average number of X objects (new) =", nmean);

     // Tansform input DF to Overdensity field
     get_overdens(this->delta_X_ini, this->delta_X_ini);

   }


#endif


   vector<real_prec>xbaux(this->NX, 0);
   vector<real_prec>pdf_in(this->NX, 0);
   this->pdf_ref.resize(this->NX, 0);
   this->pdf_ini.resize(this->NX, 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int i=0;i<this->NX; ++i)
     xbaux[i]=this->Xmin+static_cast<real_prec>(i+0.5)*(this->Xmax-this->Xmin)/(static_cast<real_prec>(this->NX));

   //Here we attemtp to apply rank ordering in order to map the pdf of the used DM to that of the DM used
   // to calibrate the bias and the kernel
   if(true==this->Apply_Rankordering)
     {

       string rest_file="_Nft"+to_string(this->Nft)+"_z"+to_string(this->redshift)+".txt";

       ULONG NGRID_NEW;
#ifdef _EXTRAPOLATE_VOLUME_
       NGRID_NEW=static_cast<ULONG>(this->Nft_low*this->Nft_low*this->Nft_low);
#else
       NGRID_NEW=this->NGRID;
#endif

       So.message_screen("Measuring pdf of log(1+delta) DM: ");
       calc_pdf("log",this->NGRID,this->NX, this->Xmax, this->Xmin, this->delta_X_ini, pdf_in);
       this->pdf_ini=pdf_in;
#ifdef WRITE_PDF
       string filex=this->Output_directory+"pdf_X_"+this->XNAME+"_MASX"+to_string(this->iMAS_X)+rest_file;
       this->File.write_to_file(filex, xbaux,pdf_in);
#endif


       string file_X_ref=this->Input_Directory_X_REF+this->Name_Catalog_X_REF_PDF; // TBDep
       this->delta_X_REF_PDF.resize(NGRID_NEW,0);
       this->File.read_array(file_X_ref, this->delta_X_REF_PDF);
       get_overdens(this->delta_X_REF_PDF, this->delta_X_REF_PDF);
       So.message_screen("Measuring pdf of log(1+delta) DM reference : ");
       calc_pdf("log",NGRID_NEW,this->NX, this->Xmax, this->Xmin, this->delta_X_REF_PDF, this->pdf_ref);
       So.DONE();
#ifdef WRITE_PDF
       filex=this->Output_directory+"pdf_X_"+this->XNAME+"_REF_MASX"+to_string(this->iMAS_X)+rest_file;
       this->File.write_to_file(filex, xbaux,this->pdf_ref);
#endif

       So.message_screen("Executing rank ordering: ");
       rankorder(this->step, xbaux, this->NGRID, this->NX,  this->Xmax, this->Xmin, this->delta_X_ini, pdf_in, this->pdf_ref);
       So.DONE();

       // Get pdf of the rank-ordered and write it to file
       So.message_screen("Measuring pdf of delta DM rank-ordered");

       calc_pdf("log",this->NGRID,this->NX, this->Xmax, this->Xmin, this->delta_X_ini, pdf_in);
       So.DONE();

#ifdef WRITE_PDF
       filex=this->Output_directory+"pdf_rank_ordered_X_"+this->XNAME+"_MASX"+to_string(this->iMAS_X)+rest_file;
       this->File.write_to_file(filex, xbaux,pdf_in);
#endif


     } // end if apply rank ordering





#ifdef _USE_MASS_KNOTS_
   this->cwclass.SKNOT_M_info.resize(this->NGRID, 0);
#endif

   // The order is important: if we do CWC inside each loop,
   // then we convolve the target DF with the updated kernel and then do the CWC
   // Otherwise one does forst the CWC once and then convolves in each iteration
#ifdef _USE_CWC_
   if(this->n_cwt > 1)
     {
#endif


#ifdef _USE_CWC_INSIDE_LOOP_
#ifndef _DO_NOT_CONVOLVE_
       this->Konvolve(this->delta_X_ini, this->delta_X, "DELTA");
#endif

#ifdef _USE_CWC_
       this->cwclass.do_CWC(this->delta_X);   //get the CW info
#ifdef _USE_MASS_KNOTS_
       this->cwclass.get_Mk_collapsing_regions(this->delta_X,nmean);  //get the MK info
#endif //use_mass_knots

#elif !defined _USE_CWC_
#if defined (_USE_MASS_KNOTS_)  ||  defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_INVARIANT_TIDAL_FIELD_III_)
       this->cwclass.do_CWC(this->delta_X);   //get the CW info
#endif // use_mass_knots

#if defined (_USE_MASS_KNOTS_)
       this->cwclass.get_Mk_collapsing_regions(this->delta_X,nmean);  //get the MK info
#endif // use_mass_knots

#endif	  // !use_cwc


          // If the terms in the bias expansion are to be used, then :
#if (defined (_USE_NABLA2DELTA_) || defined (_USE_S2_) || defined (_USE_S3_) || defined (_USE_S2DELTA_) ) && (!defined (_USE_CWC_))
       this->cwclass.get_bias_terms(this->delta_X);
#endif

       // If we do CWC but only in the initial iteration, then, for the new density field (target)
       // we do the CWC and the convolve with the kernel.

#elif !defined _USE_CWC_INSIDE_LOOP_

#ifdef _USE_CWC_
       this->cwclass.do_CWC(this->delta_X_ini);   //get the CW info
#ifdef _USE_MASS_KNOTS_
       this->cwclass.get_Mk_collapsing_regions(this->delta_X_ini,nmean);  //get the MK info
#endif
#elif !defined _USE_CWC_

#ifdef _USE_MASS_KNOTS_
       this->cwclass.do_CWC(this->delta_X);   //get the CW info
       this->cwclass.get_Mk_collapsing_regions(this->delta_X,nmean);  //get the MK info
#endif	  // use mass knots

#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_INVARIANT_TIDAL_FIELD_III_)
      this->cwclass.do_CWC(this->delta_X_ini);
#endif

#endif // end use cwc


#if (defined (_USE_NABLA2DELTA_) || defined (_USE_S2_) || defined (_USE_S3_) || defined (_USE_S2DELTA_)) && (!defined (_USE_CWC_))
       this->cwclass.get_bias_terms(this->delta_X_ini);
#endif

       this->Konvolve(this->delta_X_ini, this->delta_X, "DELTA");
#endif	  // endif !defined use cwc inside



#ifdef _USE_CWC_
   }

   else // if n_cwt=1, anyway, do the convolution with the updated Kernel
     this->Konvolve(this->delta_X_ini, this->delta_X, "DELTA");
#endif


   // #ifdef _GET_POWER_REFS_
   //  this->get_power_spectrum("DM_KONV"); // This has to be done *BEFORE* transfoming to LOg
   // #endif


   if(this->Scale_X=="log")
     {
       cout<<endl;
       // If the mass_iterative_nes is nit to be used, here we directly convert to log(numinlog+delta)
       So.message_screen("Transforming delta -> log10(2+delta)");
#pragma omp parallel for
       for(ULONG i = 0;i < this->NGRID ;++i)  //TRANSFORM DELTA TO LOG10(NUM_IN_LOG + DELTA)
         this->delta_X[i] = this->delta_X[i]<-1 ?  0 :  log10(num_in_log_x+ static_cast<real_prec>(this->delta_X[i]));
       So.DONE();
     }
   // Once the delta field is obtained, get the limits. Note that these limits mut be those of the last iteration
   // of the iterative procedure,
   // but it is not guaranteed that this is the case. So, if in _GET_BAM_REALIZATIONS_mde, Leave the limits fixed for the time being.



   /*
     real_prec aux_mn=get_min(this->delta_X);
     real_prec aux_mx=get_max(this->delta_X);
     #ifdef _USE_OMP_
     #pragma omp parallel for
     #endif
     for(ULONG i = 0;i < this->NGRID ;++i)  //TRANSFORM DELTA TO LOG10(NUM_IN_LOG + DELTA)
     this->delta_X[i] =  (this->delta_X[i]-aux_mn)/(aux_mx-aux_mn);

   */
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
#ifdef MOCK_MODE

 void Bam::get_mock_grid(string property)
 {

   bool silent=true;
   int NTHREADS;

#ifdef _USE_OMP_
#ifdef _GET_BAM_REALIZATIONS_  //if bam mocks
   NTHREADS = 1;  //This is because I have a problem with paralelization when we apply BAM to other DM fields

#ifdef _SHOW_ISSUES_
   this->So.message_screen("Issue at line", __LINE__, ": omp is not working when kernel and bias are applied to new DM field");
#endif

#elif !defined _GET_BAM_REALIZATIONS_    //if calibration
   NTHREADS = omp_get_max_threads();
#endif
   omp_set_num_threads(NTHREADS);
#elif !defined _USE_OMP_
   NTHREADS=1;
#endif


//   NTHREADS=1;

   So.message_screen("Using ",NTHREADS," threads");

   // ********************************************************************************************************************
   // ********************************************************************************************************************
   // ********************************************************************************************************************
   this->So.message_screen(" ");



   // Define seed vector for each thread in the paralellized run
   int nx=this->NX;  // Number of bins id DM

   int case_prop, ny;
   ULONG size_bias_array;


   vector<bool>filled_cells(this->NGRID, true);
   string pname,fname;

   if(property==_COUNTS_)
     {
       case_prop=1;
       size_bias_array=this->BIAS_NCOUNTS.size();
       ny = this->new_nbins_y;
       So.message_screen("**Generating Mock (tracer) number counts: ");
       cout<<endl;

       if(this->step <=this->N_iterations_Kernel)
         pname ="_iteration"+to_string(this->step);
       else
         pname= "_realization"+to_string(this->step - (this->N_iterations_Kernel)+this->N_dm_initial-1);
       fname=this->Output_directory+"MOCK_TR"+pname+"_"+"MASY"+to_string(this->iMAS_Y)+"_Nft"+to_string(this->Nft)+"_z"+to_string(this->redshift);

       this->fnameMOCK=fname;
     }


   vector<real_prec> ncounts(this->NGRID,0);

#ifdef _USE_SAT_FRACTION_
   if(property==_SAT_FRACTION_)
     {
       case_prop=3;
       size_bias_array=this->BIAS_SAT_FRACTION.size();
       ny=this->NY_SAT_FRAC;
       So.message_screen("**Generating satellite number counts**");
       if(this->step <=this->N_iterations_Kernel)
         pname ="_iteration"+to_string(this->step);
       else
         pname= "_realization"+to_string(this->step - (this->N_iterations_Kernel)+this->N_dm_initial-1);
       fname=this->Output_directory+"MOCK_TR_SAT_FRACTION"+pname+"_"+"MASY"+to_string(this->iMAS_Y)+"_Nft"+to_string(this->Nft)+"_z"+to_string(this->redshift);

       // Should we here do the same of filled_cells as woth case 2?

     }
#endif


   // Initialize counter for those cells for which no available positions were found
   ULONG counter_orphan=0;

   const gsl_rng_type *  T;
   gsl_rng * r ;

   vector<int>vseeds(NTHREADS,0);
   for(int i=0;i<vseeds.size();++i)vseeds[i]=(i*665 + (this->step+1)*1151);
   int jthread=0;
   gsl_rng_env_setup();


   time_t start_mock;
   time(&start_mock);


   // Allocate memmory for the mock density field
   this->delta_Y_new.resize(this->NGRID,0);



#ifdef _DYNAMICAL_SAMPLING_
   // Vector to allocate the Joint distribution updated after assigning Nhalos to a cell.
   vector<ULONG> X_Y_hist_dynamical(size_bias_array);

   vector<real_prec>  X_Y_hist_dynamical_normalized(size_bias_array);

   // Initialize these vectors with the original distribution in each density bin.
#if defined _USE_SAT_FRAC_ || defined _USE_MASS_FIELD_
   switch(case_prop)
     {
     case(1):
#endif
       X_Y_hist_dynamical=this->BIAS_NCOUNTS;
       X_Y_hist_dynamical_normalized=this->BIAS_NCOUNTS_normalized;
#if defined _USE_SAT_FRAC_ || defined _USE_MASS_FIELD_
       break;
#endif
#ifdef _USE_SAT_FRAC_
     case(2):
       break;
#endif
#ifdef _USE_MASS_FIELD_
     case(3):
       X_Y_hist_dynamical=this->BIAS_SAT_FRACTION;
       X_Y_hist_dynamical_normalized=this->BIAS_SAT_FRACTION_normalized;
       break;
#endif
#if defined _USE_SAT_FRAC_ || defined _USE_MASS_FIELD_
   }
#endif

   // Vector allocating the Joint distribution normalized within each Den bin, after having assigned a value of Nhalos to a cell.
   // Initialize these vectors with the original Joint and the normalized distribution in each density bin.

#endif


   ULONG lenght_bias=N_CV_BIN1 * N_CV_BIN2* N_CV_BIN3 * N_C_BIN1 * N_C_BIN2* N_C_BIN3 * this->n_sknot_massbin * this->n_cwv * this->n_vknot_massbin * this->n_cwt * nx;

   // Vector to allocate the number of cells in a given density bin during the mapping */
   vector<ULONG>Ncells_density_bin_new(size_bias_array , 0);

   // Vector containing the number of cells in a given density bin, updated everytime a cell has been assigend a value of Nhalos
   vector<ULONG> NCELLSperDMBIN_now(lenght_bias, 0);

   // Vector containing the total number if cells in a given density bin and CWT and KNOT mass
   this->NCELLSperDMBIN.clear();
   this->NCELLSperDMBIN.shrink_to_fit();
   this->NCELLSperDMBIN.resize(lenght_bias, 0);


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int i=0; i< this->NX ;++i)
     for(int w=0; w< this->n_cwt ; ++w)
       for(int v=0; v< this->n_cwv ; ++v)
         for(int k=0; k< this->n_sknot_massbin; ++k)
           for(int kv=0; kv< this->n_vknot_massbin; ++kv)
             for(int l1=0; l1< N_C_BIN1 ; ++l1)
               for(int l2=0; l2< N_C_BIN2 ; ++l2)
                 for(int l3=0; l3< N_C_BIN3 ; ++l3)
                   for(int lv1=0; lv1< N_CV_BIN1 ; ++lv1)
                     for(int lv2=0; lv2< N_CV_BIN2 ; ++lv2)
                       for(int lv3=0; lv3< N_CV_BIN3 ; ++lv3)
                         for(int j=0; j< ny ; ++j)
                           {
                             ULONG index_h=index_12d(j,i,w,k,v,kv,l1,l2,l3,lv1,lv2,lv3,this->NX,this->n_cwt,this->n_sknot_massbin,this->n_cwv,this->n_vknot_massbin,N_C_BIN1,N_C_BIN2,N_C_BIN3, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
                             ULONG index_l=index_11d(i,w,k,v,kv,l1,l2,l3,lv1,lv2,lv3, this->n_cwt,this->n_sknot_massbin,this->n_cwv,this->n_vknot_massbin,N_C_BIN1, N_C_BIN2,N_C_BIN3, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                             this->NCELLSperDMBIN[index_l]+=this->BIAS_NCOUNTS[index_h];
                           }



   {
#ifndef _EXTRAPOLATE_VOLUME_

     if(property==_COUNTS_) // We have only allowed the number counts bias to have all cells.
       {
         So.message_screen("Checking the number of cells accounted for in the bias (in function get_mock)");
         ULONG KK = static_cast<ULONG>(get_nobjects(this->BIAS_NCOUNTS));

         if(KK<this->NGRID)
           So.message_warning("Missing cells in get_mock_grid(). Perhaps density range is not wide enough to contain all cells");
         else
           So.DONE();
#endif
         real_prec Kd = static_cast<real_prec>(get_nobjects(this->BIAS_NCOUNTS_normalized));
         if(Kd<=0)
           So.message_warning("Conditional probablity ill-defined");
       }
   }


#ifdef _USE_OMP_
#pragma omp parallel private (jthread, r, T)
   {
#endif

     T = gsl_rng_default;
     r = gsl_rng_alloc (T);

#ifdef _USE_OMP_
     jthread=omp_get_thread_num();
#endif

     gsl_rng_default_seed=vseeds[jthread];



     // First block, meant to construct halo number counts



#ifdef _USE_OMP_
#pragma omp for reduction(+:counter_orphan)
#endif
     //Start loop over the cells
     for(ULONG i=0;i<this->NGRID;++i)
       {
         real_prec dm = static_cast<real_prec>(this->delta_X[i]);

         int I_X  = get_bin(dm,this->s_mins.prop1,this->NX,this->s_deltas.prop1,this->bin_accumulate_borders);
         // Get the bin in the delta dark matter (or log 1+delta) in each cell


         // **********CWT

         int I_CWT=0;
#ifdef _USE_CWC_
         I_CWT=this->cwclass.get_Tclassification(i);
#endif

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

         int I_VK=0;
#ifdef _USE_VEL_KNOTS_V_
         I_VK= (this->cwclass.cwv_used[I_CWV]== I_KNOT ? this->cwclass.VDISP_KNOT_info[i]: 0);
#endif


         // Get the corresponding bin in the two invariants of the shear of the tidal field
         int I_C1=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
         real_prec C1 = this->cwclass.Invariant_TF_II[i];
         I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#elif defined _USE_DELTA2_
         real_prec C1 = this->cwclass.DELTA2[i];
         I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#endif



         int I_C2=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
         real_prec C2 = this->cwclass.Invariant_TF_III[i];

         I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#elif defined _USE_DELTA3_
         real_prec C2 = this->cwclass.DELTA3[i];
         I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#endif

         int I_C3=0;
#ifdef _USE_TIDAL_ANISOTROPY_
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
#if defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_S2_)
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

                       if(true==filled_cells[i])

                         {

                           // Get the bin in DM


                           // check if prob is zero for all possible valoes of Nh
                           // if aux_h=0, given that the quantity CWT_hist is always > 0, means that all values of prob in this bin are zero
                           // It means then that no cells with this DM were found to have at least one halo, with respect to the BIAS_NCOUNTS
                           // This is likely to happen when we apply the Bias and the Kernel to a different DM field in order toe create a mock.
                           real_prec aux_h=0;
                           switch(case_prop)
                             {
                             case(1):
                               for(int ih=0;ih<ny;++ih)
                                 aux_h+=this->BIAS_NCOUNTS_normalized[index_12d(ih,I_X,I_CWT,I_MK,I_CWV,I_VK, I_C1,I_C2,I_C3,I_CV1,I_CV2,I_CV3,this->NX,this->n_cwt, this->n_sknot_massbin,this->n_cwv, this->n_vknot_massbin,N_C_BIN1,N_C_BIN2,N_C_BIN3,N_CV_BIN1,N_CV_BIN2,N_CV_BIN3)];
                               break;

#ifdef _USE_SAT_FRACTION_
                             case(3):
                               for(int ih=0;ih<ny;++ih)
                                 aux_h+=this->BIAS_SAT_FRACTION_normalized[index_12d(ih,I_X,I_CWT,I_MK,I_CWV,I_VK,I_C1,I_C2,I_C3,I_CV1,I_CV2,I_CV3,this->NX,this->n_cwt, this->n_sknot_massbin,this->n_cwv, this->n_vknot_massbin,N_C_BIN1,N_C_BIN2,N_C_BIN3,N_CV_BIN1,N_CV_BIN2,N_CV_BIN3)];
                               break;
#endif
                           }

                           int Nhalos_orfan;
                           if(aux_h<0)
                             exit(0);



                           if(aux_h>0)
                             {

                               int halo_prop;
                               ULONG index_el=index_11d(I_X,I_CWT,I_MK,I_CWV,I_VK,I_C1,I_C2,I_C3,I_CV1,I_CV2,I_CV3,this->n_cwt,this->n_sknot_massbin,this->n_cwv, this->n_vknot_massbin,N_C_BIN1,N_C_BIN2,N_C_BIN3,N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
                               ULONG N_available_positions_in_denbin_cwt_Mk = this->NCELLSperDMBIN[index_el];

                               // Number of cells in the density bin and CWT.
                               // This is always greater than zero, by construction. A fixed value inside the loop
                               ULONG N_used_positions_in_denbin_cwt_Mk = NCELLSperDMBIN_now[index_el];

                               // Number of used cells for the current bin of DM and CWT

                               bool flag = true;
                               // If the density bin to which this cell belongs to is already filled
                               // then go below to get a value of Nhalos selected according to the original distribution in this density bin
                               // setting a flag=false

                               if (N_available_positions_in_denbin_cwt_Mk == N_used_positions_in_denbin_cwt_Mk)
                                 flag=false;


                               // If the density bin still has available positions to be assigned, then proceed
                               if(true==flag)
                                 {
                                   bool cell_accepted=false;

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

                                       // Probability for number counts
                                       while (prob_ncounts<ran)
                                         {
#if defined _USE_MASS_FIELD_ || defined _USE_SAT_FRACTION_
                                           if(1==case_prop)
#endif

//                                               halo_prop= static_cast<int>(floor((1.0+static_cast<real_prec>(this->nmax_y_onecell))*gsl_rng_uniform(r))); //draw number counts in cells
                                              halo_prop= gsl_rng_uniform_int(r,this->nmax_y_onecell+1); //draw number counts in cells
#ifdef _USE_MASS_FIELD_
                                           else if(2==case_prop)
                                             halo_prop= static_cast<int>(floor(static_cast<real_prec>(ny)*gsl_rng_uniform(r)));                         // draw bins n the mass distribution
#endif
#ifdef _USE_SAT_FRACTION_
                                           else if(3==case_prop)
                                             halo_prop= static_cast<int>(floor((1.0+static_cast<real_prec>(this->nmax_y_sat_onecell))*gsl_rng_uniform(r)));  //drwa number counts of satellites in cells
#endif
                                           ran   = gsl_rng_uniform(r);
#ifdef _DYNAMICAL_SAMPLING_
                                           prob_ncounts  = static_cast<real_prec>(X_Y_hist_dynamical_normalized[index_12d(halo_prop,I_X,I_CWT,I_MK,I_CWV,I_VK,I_C1,I_C2,I_C3,I_CV1,I_CV2,I_CV3,this->NX,this->n_cwt,this->n_sknot_massbin,this->n_cwv,this->n_vknot_massbin,N_C_BIN1,N_C_BIN2,N_C_BIN3,N_CV_BIN1,N_CV_BIN2,N_CV_BIN3)]);

#else //this case is just as the one usnig mass, isn't it?
                                           prob_ncounts  = static_cast<real_prec>(this->BIAS_NCOUNTS_normalized[index_12d(halo_prop,I_X,I_CWT,I_MK,I_CWV,I_VK,I_C1,I_C2,I_C3,I_CV1,I_CV2,I_CV3,this->NX,this->n_cwt,this->n_sknot_massbin,this->n_cwv,this->n_vknot_massbin,N_C_BIN1,N_C_BIN2,N_C_BIN3, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3)]);

#endif

                                       }

                                       ULONG index=0;
                                       ULONG index_low=0;
                                       index      =index_12d(halo_prop,I_X,I_CWT, I_MK, I_CWV, I_VK,I_C1,I_C2, I_C3, I_CV1,I_CV2,I_CV3, nx, this->n_cwt, this->n_sknot_massbin,this->n_cwv, this->n_vknot_massbin ,N_C_BIN1,N_C_BIN2,N_C_BIN3,N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
                                       index_low  =index_11d(I_X, I_CWT, I_MK, I_CWV, I_VK,I_C1,I_C2, I_C3, I_CV1,I_CV2, I_CV3, this->n_cwt, this->n_sknot_massbin, this->n_cwv, this->n_vknot_massbin,N_C_BIN1,N_C_BIN2,N_C_BIN3,N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);

                                       ULONG N_available_positions;
#ifdef _USE_SAT_FRACTION_
                                       switch(case_prop)
                                         {
                                         case(1):
#endif
                                           N_available_positions=this->BIAS_NCOUNTS[index];
#ifdef _USE_SAT_FRACTION_
                                           break;
                                         case(3):
                                           N_available_positions=this->BIAS_SAT_FRACTION[index];
                                           break;
                                         }
#endif

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

                                           // Add one (i.e, the accepted) to count the number of positions used in the current den-N bin
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                                           Ncells_density_bin_new[index]++;

#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                                           NCELLSperDMBIN_now[index_low]++;

                                           // -----------------------------------------------------------------
                                           // In the current Den-N bin, subtract one (i.e, the accepted) in order
                                           // to upgrade the distribution and give more weight to the remaining available positions
                                           // Attention, do not ask whether this is <0, for it is an unsigned long

#ifdef _DYNAMICAL_SAMPLING_
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

                                           vector<long>AUX(ny,0);

                                           for(int j=0; j< ny ;++j)

                                             // Here we CANNOT use the definitions of index and index_low, for those depend on Nhalos, while here they depend on j
                                             AUX[j]=X_Y_hist_dynamical[index_12d(j,I_X,I_CWT,I_MK,I_CWV, I_VK,I_C1,I_C2,I_C3,I_CV1,I_CV2,I_CV3, nx,this->n_cwt,this->n_sknot_massbin,this->n_cwv,this->n_vknot_massbin,N_C_BIN1,N_C_BIN2,N_C_BIN3, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3)];


                                           long lkk=get_max<long>(AUX);

                                           for(int j=0;j< ny ;++j)
                                             // Here we CANNOT use the definitions of index and index_low, for those depend on Nhalos, while here they depend on j
                                             X_Y_hist_dynamical_normalized[index_12d(j,I_X,I_CWT,I_MK,I_CWV, I_VK,I_C1,I_C2,I_C3,I_CV1,I_CV2,I_CV3, this->NX,this->n_cwt,this->n_sknot_massbin,this->n_cwv,this->n_vknot_massbin,N_C_BIN1,N_C_BIN2,N_C_BIN3,N_CV_BIN1,N_CV_BIN2,N_CV_BIN3)]= (lkk == 0  ? 0.0 :  static_cast<real_prec>(AUX[j])/static_cast<real_prec>(lkk));

                                           AUX.clear();
                                           AUX.shrink_to_fit();
#endif

                                         }   // end if(N_cells > N_cells_now).

                                     }// end while(cell_accepted==false && true==flag[I_X])
                                 } // end of if(true==flag)
                               else
                                 // If the density bin has been already filled, we call the distribution in that bin to assign randomly the value Nhalo according to the other properties
                                 {
                                   real_prec prob=-1.0;
                                   real_prec ran=0.0;
                                   while(prob<ran)
                                     {
#ifdef _USE_SAT_FRAC_
                                       switch(case_prop)
                                         {
                                         case(1):
#endif
//                                           Nhalos_orfan= static_cast<int>(floor((1.0+static_cast<real_prec>(this->nmax_y_onecell))*gsl_rng_uniform(r)));
                                           Nhalos_orfan= gsl_rng_uniform_int(r,this->nmax_y_onecell+1); //draw number counts in cells

                                           prob = this->BIAS_NCOUNTS_normalized[index_12d(Nhalos_orfan,I_X,I_CWT,I_MK, I_CWV, I_VK, I_C1,I_C2, I_C3, I_CV1,I_CV2, I_CV3,this->NX,this->n_cwt, this->n_sknot_massbin,this->n_cwv, this->n_vknot_massbin,N_C_BIN1,N_C_BIN2, N_C_BIN3, N_CV_BIN1,N_CV_BIN2, N_CV_BIN3)];
#ifdef _USE_SAT_FRAC_
                                           break;
                                         case(3):
                                           Nhalos_orfan= static_cast<int>(floor((1.0+static_cast<real_prec>(this->nmax_y_sat_onecell))*gsl_rng_uniform(r)));
                                           prob = this->BIAS_SAT_FRACTION_normalized[index_12d(Nhalos_orfan,I_X,I_CWT,I_MK, I_CWV, I_VK,I_C1,I_C2, I_C3, I_CV1,I_CV2, I_CV3,this->NX,this->n_cwt, this->n_sknot_massbin,this->n_cwv, this->n_vknot_massbin,N_C_BIN1,N_C_BIN2, N_C_BIN3, N_CV_BIN1,N_CV_BIN2, N_CV_BIN3)];
                                           break;
                                         }
#endif
                                       ran = gsl_rng_uniform(r);
                                       counter_orphan++;
                                     }

                                   this->delta_Y_new[i]= static_cast<real_prec>(Nhalos_orfan);

                                 }// end else
                             }// end of if(aux_h>0)
                           else // if in the theta_bin the prob is zero always, then assign random; this won't solve the issue of less assigned objects
                             {
//                               if(this->step==this->N_iterations_Kernel)
                                 this->delta_Y_new[i]=gsl_rng_uniform_int(r,this->nmax_y_onecell+1);
                             }
                   }// end if(dm is in the edfined range)
       }// end loop over cells
     gsl_rng_free (r);
#ifdef _USE_OMP_
   }// end parallelized region
#endif


   So.message_time_mock(start_mock);

   // ******************************************** end assigning ncounts or other prop*******************************


   if(counter_orphan>0)
     So.message_screen("Fraction of orfan cells =", 100.0*static_cast<real_prec>(counter_orphan)/static_cast<real_prec>(this->NGRID), "%");

   if(property==_COUNTS_)
     {

       this->Nobjects=get_nobjects(this->delta_Y_new);
       So.message_screen("Number of objects = ", this->Nobjects);
       So.message_screen("Mean number denisty = ", this->Nobjects/pow(Lbox, 3), "(Mpc / h )⁻³");

#ifdef _DO_BAM_CALIBRATION_
#ifndef _EXTRAPOLATE_VOLUME_
       if(this->Nobjects < this->N_objects_Y)
         So.message_screen("Less objects assigned to mock. Fraction = ", fabs(100.0-100.0*(static_cast<real_prec>(this->Nobjects)/this->N_objects_Y)), "%");
       if(this->Nobjects  > this->N_objects_Y)
         So.message_screen("More objects assigned to mock. Fraction = ", fabs(100.0-100.0*(static_cast<real_prec>(this->Nobjects)/this->N_objects_Y)), "%");
#endif
#endif


       // ----------------------------------------------------------------------------------------------
       if(false==silent)
         {
           real_prec anmean=static_cast<real_prec>(this->Nobjects)/static_cast<real_prec>(this->NGRID);
           So.message_screen("<N> objects in mock =", anmean);


           int nyy=30;
           real_prec dmax=300.0;

           vector<real_prec>DELTA_HIST(nyy,0);

           real_prec delta_delta= (dmax+1.0)/static_cast<real_prec>(nyy);  //deltamax - deltamin
#ifdef _USE_OMP_
#pragma omp for
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
#pragma omp for
#endif
           for(ULONG i=0; i< this->NGRID ;++i)
             AUX[i]=this->delta_Y_new[i]/static_cast<real_prec>(anmean) - 1.0;

           lkk=get_max<real_prec>(AUX);
           So.message_screen("Maximum delta in mock =", lkk);

           lkk=get_min<real_prec>(AUX);
           So.message_screen("Minimum delta in mock =", lkk);
           AUX.clear();
           AUX.shrink_to_fit();

           // ----------------------------------------------------------------------------------
           if(true==this->Write_PDF_number_counts) 	  // Write the PDF of the outcome */
             {
               string fileY=this->Output_directory+"PDF_NC"+"_Y_MOCK_"+this->new_Name_Property_Y+"_MASY"+to_string(this->iMAS_Y)+"_Nft"+to_string(this->Nft)+"_z"+to_string(this->redshift)+".txt";

               So.message_screen("Writting PDF_NC for Y in file ", fileY);
               int NPART=20; // Number of particles in cells

               this->PDF_NC_Y.clear();
               this->PDF_NC_Y.resize(NPART, 0);

#ifdef _USE_OMP_
#pragma omp for
#endif
               for(ULONG i=0;i< this->NGRID;++i)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 this->PDF_NC_Y[static_cast<int>(this->delta_Y_new[i])]++;


               this->File.write_to_file_i(fileY,this->PDF_NC_Y);

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

     }




#ifdef  _EXTRAPOLATE_VOLUME_
#ifdef _WRITE_TR_DENSITY_FIELD_
   this->File.write_array(fname, this->delta_Y_new);
#endif
#endif


   // this is commented for only was applicable to DENSITY case. For the case in which we calibrate ncounts
   // and other properties, those extra fields are written directly as they come in the last iteration
   if(case_prop==2)
     {
       gsl_rng_env_setup();
       gsl_rng_default_seed=75;
       const gsl_rng_type *  Tn= gsl_rng_ranlux;
       gsl_rng * rn = gsl_rng_alloc (Tn);

#ifndef _USE_LOG_MASS_


       real_prec mass_min=pow(10,params._LOGMASSmin())*params._MASS_units();

#pragma omp parallel for
       for(ULONG i=0;i< this->NGRID;++i)
         {
           real_prec xr=gsl_rng_uniform (rn);
           real_prec aux=this->s_mins.prop0_mass+(static_cast<int>(delta_Y_new[i])+xr)*this->s_deltas.prop0_mass;   // Get value of log10(2+delta_y)
           real_prec mass_cell = (pow(10,aux)-NUM_IN_LOG)*MASS_SCALE;
           this->delta_Y_new[i]=mass_cell*filled_cells[i];

           int Neff = static_cast<int>(floor(this->delta_Y_new[i]/mass_min));

           // This weights help to upweight the mass at a cell such that, give the number counts it has, it can at least provide mass for all particles
           real_prec weight = (Neff >= ncounts[i] ? 1.0 : static_cast<real_prec>(ncounts[i])/static_cast<real_prec>(Neff));
           this->delta_Y_new[i]*=weight;

         }
       gsl_rng_free (rn);
       So.message_screen("Minimum of new mass field = ", get_min(this->delta_Y_new));
       So.message_screen("Maximim of new mass field = ", get_max(this->delta_Y_new));
#endif

     }

   ncounts.clear();
   ncounts.shrink_to_fit();

   //  if(case_prop==3)  // do nothing, we deal with counts of satellites
   //   {
   //#ifdef _USE_OMP_
   //#pragma omp for
   //#endif
   //    for(ULONG i=0;i< this->NGRID;++i) //tbc
   //      {  /// cambiar ldelta_Y por le mínimo de delta_Y_new
   //         this->delta_Y_new[i]=this->s_mins.prop0_sf+(static_cast<int>(delta_Y_new[i])+0.5)*this->s_deltas.prop0_sf;   // Get value of log10(2+delta_y)
   //        this->delta_Y_new[i]=this->Mean_density_Y_SAT_FRAC*(1.0+pow(10,delta_Y_new[i])-NUM_IN_LOG);                                        // Get nbar*(1+delta_y)
   //   }
   // }

   // Writ ethe density fields of properties of mock tracers:.



#ifndef _TEST_THRESHOLDS_RESIDUALS_
   auto out_it = std::find(std::begin(this->output_at_iteration), std::end(this->output_at_iteration), this->step);
   if((this->step==this->N_iterations_Kernel) || (out_it != std::end(this->output_at_iteration)))
     this->File.write_array(fname, this->delta_Y_new);
#endif

#ifdef _GET_BAM_REALIZATIONS_

   this->File.write_array(fname, this->delta_Y_new);


   if(case_prop==1)
     this->patchy.fname_MOCK_NCOUNTS=fname;
   if(case_prop==2)
     this->patchy.fname_MOCK_MASS=fname;
   if(case_prop==3)
     this->patchy.fname_MOCK_NCOUNTS_SAT=fname;

#endif

   // Get power spectrum only number counts


   if(case_prop==1)
     this->get_power_spectrum("TR_MOCK");


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

 void Bam::assign_tracer_mass_new()
 {

   So.message_screen("**Assigning masses");
   cout<<endl;

   int NTHREADS;


#ifdef _USE_OMP_
#ifdef _GET_BAM_REALIZATIONS_  //if bam mocks
 NTHREADS = 1;  //This is because I have a problem with paralelization when we apply BAM to other DM fields
 //NTHREADS = omp_get_max_threads();


#elif !defined _GET_BAM_REALIZATIONS_    //if calibration
    NTHREADS = omp_get_max_threads();
 #endif
    omp_set_num_threads(NTHREADS);
 #elif !defined _USE_OMP_
    NTHREADS=1;
#endif

   const gsl_rng_type *  T;
   gsl_rng * r ;

   // ********************************************************************************************************************
   // ********************************************************************************************************************
   // ************************************************RENAME SOME BIN PROPPERTIES********************************************************************
   int N_a=N_C_BIN1;
#ifdef _USE_TRACERS_IN_CELLS_
 vector<real_prec> MOCK_DEN_FIELD(this->NGRID,0);
 this->tracer.get_density_field_grid(_COUNTS_,MOCK_DEN_FIELD);
 int nmax=get_max<real_prec>(MOCK_DEN_FIELD);
 So.message_screen("Maximum number of tracer in cell", nmax);
 if(nmax>=N_TRACERS_IN_CELLS_MAX){
   So.message_warning("Please increase the maximum number of tracer in one cell to " ,nmax);
   exit(0);
 }
 N_a = N_TRACERS_IN_CELLS_MAX;
#endif

 // ************************************************RENAME SOME BIN PROPPERTIES********************************************************************
 int N_b=N_C_BIN2;
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
 this->tracer.get_neighbour_tracers(this->ncells_info);
 N_b = N_NEIGHBOURS_MAX;
#elif defined (_USE_MIN_DISTANCE_TO_NEIGHBOURS_) || defined (_USE_LOCAL_CLUSTERING_)
 N_b=N_BINS_MIN_DIST_TO_NEI;
#endif

 // ************************************************RENAME SOME BIN PROPPERTIES********************************************************************

 int N_c = N_C_BIN3;
#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
 this->tracer.get_min_separation_in_cell();
 N_c = N_BINS_MIN_SEP_IN_CELLS;
#ifdef _MASS_ASSIGNMENT_TO_REFERENCE_
  real_prec delta_min_sep = (MAX_SEP_IN_CELLS-this->min_halo_separation)/static_cast<real_prec>(N_BINS_MIN_SEP_IN_CELLS);
#else
  real_prec delta_min_sep = DELTA_MIN_SEP;
#endif
#endif

   // ********************************************************************************************************************

   // Define seed vector for each thread in the paralellized run
   int ny_mf = this->params._NMASSbins_mf();

   vector<int>vseeds(NTHREADS,0);
   for(int i=0;i<vseeds.size();++i)vseeds[i]=i*66 + (this->step+1)*11;
   int jthread=0;
   gsl_rng_env_setup();


   time_t start_mock;
   time(&start_mock);

   ULONG lenght_bias_full= N_CV_BIN1 * N_CV_BIN2* N_CV_BIN3 * N_a * N_C_BIN2* N_c * this->n_sknot_massbin * this->n_cwt * this->n_vknot_massbin * this->n_cwv* this->NX*ny_mf;
   ULONG lenght_bias     = N_CV_BIN1 * N_CV_BIN2* N_CV_BIN3 * N_a * N_C_BIN2* N_c * this->n_sknot_massbin * this->n_cwt*this->n_vknot_massbin*this->n_cwv*this->NX;


   // Vector to allocate the Joint distribution updated after assigning Nhalos to a cell.
   // Initialize these vectors with the original distribution in each density bin.
   vector<ULONG> X_Y_hist_dynamical(lenght_bias_full);
   X_Y_hist_dynamical=this->ABUNDANCE;

   vector<real_prec>  X_Y_hist_dynamical_normalized(lenght_bias_full);
   X_Y_hist_dynamical_normalized=this->ABUNDANCE_normalized;

   // Vector to allocate the number of cells in a given density bin during the mapping */
   vector<ULONG>Ncells_density_bin_new(X_Y_hist_dynamical.size() , 0);

   // Vector containing the number of cells in a given density bin, updated everytime a cell has been assigend a value of Nhalos
   vector<ULONG> NCELLSperDMBIN_now(lenght_bias, 0);


   /*
   // Vector containing the total number if cells in a given density bin and CWT and KNOT mass
   vector<ULONG> NCELLSperDMBIN(lenght_bias, 0);

   So.message_screen("Marginalizing wrt mass bins...");
#pragma omp parallel for collapse(11)
   for(int i=0; i< this->NX ;++i)
     for(int w=0; w< this->n_cwt ; ++w)
       for(int k=0; k< this->n_sknot_massbin; ++k)
         for(int v=0; v< this->n_cwv ; ++v)
           for(int kv=0; kv< this->n_vknot_massbin; ++kv)
             for(int l1=0; l1< N_a ; ++l1)
               for(int l2=0; l2< N_C_BIN2 ; ++l2)
                 for(int l3=0; l3< N_c ; ++l3)
                   for(int lv1=0; lv1< N_CV_BIN1 ; ++lv1)
                     for(int lv2=0; lv2< N_CV_BIN2 ; ++lv2)
                       for(int lv3=0; lv3< N_CV_BIN3 ; ++lv3)
                         {
                           ULONG index_l=index_11d(i,w,k,v,kv,l1,l2,l3,lv1,lv2,lv3,this->n_cwt,this->n_sknot_massbin,this->n_cwv,this->n_vknot_massbin,N_a,N_C_BIN2,N_c, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
                           for(int j=0; j< ny_mf ; ++j)
                             {
                               ULONG index_h=index_12d(j,i,w,k,v,kv,l1,l2,l3,lv1,lv2,lv3,this->NX,this->n_cwt,this->n_sknot_massbin,this->n_cwv,this->n_vknot_massbin,N_a,N_C_BIN2,N_c, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
#pragma omp atomic update
NCELLSperDMBIN[index_l]+=X_Y_hist_dynamical[index_h];
                             }
                         }
    So.DONE();

    */
   // ********************************************************************************************************************

   int nbins_mf=this->tracer_ref.mass_function.size();

       this->So.message_screen("Using reference mass function to complement assignment.");
       this->tracer.define_property_bins();

       /*
         this->So.message_screen("Using reference mass function to complement assignment.");
         this->tracer.type_of_object="TRACER";
         this->tracer.set_params_catalog(this->params);
         this->tracer.define_property_bins();
         string massf_file=this->Output_directory+"tracer_ref_mass_function.txt";

         vector<real_prec>massf;
         int nbins_mf=this->File.read_file(massf_file,massf,1);
         this->mass_bins.resize(nbins_mf,0);
         this->mfunc.resize(nbins_mf,0);

         #pragma omp parallel for
         for(int i=0; i< nbins_mf; ++i)
         this->mass_bins[i]=massf[0+i*2];

         this->mass_min=this->mass_bins[0];
         this->mass_max=this->mass_bins[nbins_mf-1];
       */

       this->mfunc.resize(nbins_mf,0);
       this->mass_min=this->tracer_ref.MBin[0];
       this->mass_max=this->tracer_ref.MBin[nbins_mf-1];

       // I cannot use here the this->tracer MBmin and max for it might happen that the reference mass function has been measured with a different number of bins
       // different to the current this>tracer_ref_NMBINS
       vector<real_prec>delta_mass_aux(nbins_mf,0);
       for(int i=0; i< nbins_mf; ++i)
         {
           real_prec lmmin=log10(mass_min)+i*log10(mass_max/mass_min)/static_cast<real_prec>(nbins_mf);
           real_prec lmmax=log10(mass_min)+(i+1)*log10(mass_max/mass_min)/static_cast<real_prec>(nbins_mf);
           delta_mass_aux[i]=pow(10,lmmax)-pow(10,lmmin);
         }

#pragma omp parallel for
       for(int i=0; i< nbins_mf; ++i)
         this->mfunc[i]=this->tracer_ref.mass_function[i]*delta_mass_aux[i];
       delta_mass_aux.clear(); delta_mass_aux.shrink_to_fit();

       real_prec max_mf=static_cast<real_prec>(get_max(mfunc));
#pragma omp parallel for
       for(int i=0; i< nbins_mf; ++i)
         this->mfunc[i]/=max_mf;



   gsl_interp_accel *acc = gsl_interp_accel_alloc ();
   gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, this->tracer_ref.mass_function.size());
   gsl_spline_init (spline, &(this->tracer_ref.MBin[0]), &(this->mfunc[0]), this->tracer_ref.mass_function.size());

   real_prec lm_min=this->params._LOGMASSmin();


   // ********************************************************************************************************************

   ULONG naux=0;
#pragma omp parallel for reduction(+:naux)
   for(ULONG i =0; i< X_Y_hist_dynamical.size(); ++i)
     naux+=X_Y_hist_dynamical[i];
   So.message_screen("Number of tracers in abundance (from reference) =", naux);
   So.message_screen("Number of tracers in catalog =", this->tracer.NOBJS);
   So.DONE();
   // ********************************************************************************************************************


   // This vector is defiend in order to allocate the convolution of the DM field with the mass-weighted kernel

#ifdef _USE_ITERATIVE_MASS_ASSIGNMENT_NEW_
   this->delta_dm_aux.clear();
   this->delta_dm_aux.shrink_to_fit();
   this->delta_dm_aux.resize(this->NGRID,0);
   // Here we take the new DM field, (i.e, the new Dm convlvede with tthe BAM kernel) and convolve it with the kernel
   // computed from the mass weighted power spectrum
   if(this->step_mass_assignemt>0)
     this->Konvolve(this->delta_X,this->delta_X,"DELTA"); //recall that if you comment this line, must also comment line 2389

   So.message_screen("Transforming delta -> log10(2+delta)");
   for(ULONG i = 0;i < this->NGRID ;++i)  //TRANSFORM DELTA TO LOG10(NUM_IN_LOG + DELTA)
     this->delta_dm_aux[i] = this->delta_X[i]<-1 ?  0 :  log10(NUM_IN_LOG+ static_cast<real_prec>(this->delta_X[i]));
   So.DONE();
#endif


   So.message_screen("Going through mock catalog ...");

   // Initialize counter for those cells for which no available positions were found
   int counter_orphan=0;
   int counter_mf=0;
   int counter_massless=0;
   int massless=0;
   int Nhalos_orfan=0;
   int counter_in=0;

#ifdef _USE_OMP_
#pragma omp parallel private (jthread, r, T)
   {
#endif
     T = gsl_rng_default;
     r = gsl_rng_alloc (T);

#ifdef _USE_OMP_
     jthread=omp_get_thread_num();
#endif
     gsl_rng_default_seed=vseeds[jthread];

#ifdef _USE_OMP_
#pragma omp for reduction(+:counter_orphan, counter_mf,counter_massless, counter_in)
#endif
     //Start loop over tracer catalog
     for(ULONG i=0;i<this->tracer.NOBJS;++i)
       {

         ULONG id=this->tracer.Halo[i].GridID;  // Get the cell ID where the tracer is located

         #ifdef _USE_ITERATIVE_MASS_ASSIGNMENT_NEW_
         real_prec dm = static_cast<real_prec>(this->delta_dm_aux[id]);
#else
         real_prec dm = static_cast<real_prec>(this->delta_X[id]);
#endif

         int I_X  = get_bin(dm,this->s_mins.prop1,this->NX,this->s_deltas.prop1,this->bin_accumulate_borders);

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
#ifdef _USE_TRACERS_IN_CELLS_
         I_C1 = static_cast<int>(MOCK_DEN_FIELD[id]);
#elif defined _USE_INVARIANT_TIDAL_FIELD_II_
         real_prec C1 = this->cwclass.Invariant_TF_II[id];
         I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#elif defined _USE_DELTA2_
         real_prec C1 = this->cwclass.DELTA2[id];
         I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#endif


         int I_C2=0;
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
         I_C2=this->tracer.Number_of_neighbours[i];
         if(I_C2==N_NEIGHBOURS_MAX)
           I_C2=N_NEIGHBOURS_MAX-1;
#elif defined _USE_INVARIANT_TIDAL_FIELD_III_
         real_prec C2 = this->cwclass.Invariant_TF_III[id];
         I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#elif defined _USE_DELTA3_
         real_prec C2 = this->cwclass.DELTA3[id];
         I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#endif


         int I_C3=0;
#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
         real_prec min_sep = this->tracer.min_separation_in_cell[id];
         I_C3 = get_bin(min_sep, 0, N_c, delta_min_sep , this->bin_accumulate_borders);
#elif defined _USE_TIDAL_ANISOTROPY_
         real_prec C3 = this->cwclass.Tidal_Anisotropy[id];
         I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined _USE_S2_
         real_prec C3 = this->cwclass.S2[id];             // s²
         I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#endif

         int I_CV1=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
         real_prec CV1 = invariant_field_I(this->cwclass.lambda1_vs[id],this->cwclass.lambda2_vs[id],this->cwclass.lambda3_vs[id]); // Not yet assigned to container
         I_CV1= get_bin(CV1, CV1_MIN, N_CV_BIN1,DELTA_CV1, this->bin_accumulate_borders);
#elif defined _USE_NABLA2DELTA_
         real_prec CV1 = this->cwclass.N2D[id];      // Nabla² ð
         I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1, s_deltas.prop7,this->bin_accumulate_borders);
#endif

         int I_CV2=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
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
         if(dm >=this->s_mins.prop1 && dm <=this->s_maxs.prop1)
#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_DELTA2_)
           if(C1<= this->s_maxs.prop4 && C1>= this->s_mins.prop4)
#endif
#if defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_DELTA3_)
             if(C2<= this->s_maxs.prop5 && C2>= this->s_mins.prop5)
#endif
#if defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_S2_)
               if(C3<= this->s_maxs.prop6 && C3>= this->s_mins.prop6)
#endif
#if defined (_USE_INVARIANT_SHEAR_VFIELD_I_) || defined (_USE_NABLA2DELTA_)
                 if(CV1<= this->s_maxs.prop7 && CV1>= this->s_mins.prop7)
#endif
#if defined (_USE_INVARIANT_SHEAR_VFIELD_II_) || defined (_USE_S2DELTA_)
                   if(CV2<= this->s_maxs.prop8 && CV2>= this->s_mins.prop8)
#endif
#if defined (_USE_INVARIANT_SHEAR_VFIELD_III_) || defined (_USE_S3_)
                     if(CV3<= this->s_maxs.prop9 && CV3>= this->s_mins.prop9)
#endif
                       {
                         ULONG index_low  =  index_11d(I_X,I_CWT,I_MK,I_CWV,I_VK,I_C1,I_C2,I_C3,I_CV1,I_CV2,I_CV3, this->n_cwt,this->n_sknot_massbin,this->n_cwv, this->n_vknot_massbin,N_a,N_b,N_c,N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
                         ULONG N_available_positions_in_denbin_cwt_Mk =this->NCELLSperDMBIN[index_low];
                         ULONG N_used_positions_in_denbin_cwt_Mk = NCELLSperDMBIN_now[index_low];
                         int halo_prop;

                         if(N_available_positions_in_denbin_cwt_Mk>0) //as long as we have available positions
                           {

                             bool flag = true;
                             if(N_available_positions_in_denbin_cwt_Mk == N_used_positions_in_denbin_cwt_Mk)
                               flag=false;

                             // If the density bin still has available positions to be assigned, then proceed
                             if(true==flag)
                               {
                                 bool cell_accepted=false;

                                 // This density bin as available positions. Proceed until the cell is assigned a value of Nhalo
                                 while (false==cell_accepted)
                                   {
                                     real_prec prob_ncounts=-10.0;
                                     real_prec ran = 10.0;

                                     while (prob_ncounts<ran)
                                       {
                                         halo_prop= gsl_rng_uniform_int(r,ny_mf);// draw mass bins used to measure the conditional masas function abundance
                                         ULONG index_ran=index_12d(halo_prop,I_X,I_CWT,I_MK,I_CWV,I_VK,I_C1,I_C2,I_C3,I_CV1,I_CV2,I_CV3,this->NX,this->n_cwt,this->n_sknot_massbin,this->n_cwv,this->n_vknot_massbin,N_a,N_b,N_c,N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
                                         prob_ncounts  = static_cast<real_prec>(X_Y_hist_dynamical_normalized[index_ran]);
                                         ran   = gsl_rng_uniform(r);
                                       }

                                     ULONG index = index_12d(halo_prop,I_X,I_CWT, I_MK, I_CWV,I_VK,I_C1,I_C2, I_C3, I_CV1,I_CV2,I_CV3, this->NX, this->n_cwt, this->n_sknot_massbin,this->n_cwv,this->n_vknot_massbin,N_a,N_b,N_c,N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
                                     ULONG N_available_positions=this->ABUNDANCE[index];
                                     ULONG N_used_positions = Ncells_density_bin_new[index];

                                     if(N_available_positions > N_used_positions)
                                       {
                                         // the abundance is measured in bins of log Mass, and hence the bin index halo_prop is in that space
                                         // Note that MBmin and MBmax comes already in M, (not in log M) from Catalog class.
                                         real_prec xr = static_cast<real_prec>(gsl_rng_uniform(r));
                                         real_prec fraction_mass= xr*log10(this->tracer.MBmax[halo_prop]/this->tracer.MBmin[halo_prop]);
                                         real_prec lmass_halo = log10(this->tracer.MBmin[halo_prop])+fraction_mass;
                                         this->tracer.Halo[i].mass = pow(10,lmass_halo);
                                         cell_accepted=true;
                                         counter_in++;
#pragma omp atomic update
                                         Ncells_density_bin_new[index]++;
#pragma omp atomic update
                                         NCELLSperDMBIN_now[index_low]++;

                                         if(X_Y_hist_dynamical[index] >=1)
#pragma omp atomic update
                                           X_Y_hist_dynamical[index]--;
                                         else
                                           X_Y_hist_dynamical[index]=0;

                                         vector<long>AUX(ny_mf,0);
                                         for(int j=0; j< ny_mf ;++j)
                                           AUX[j]=X_Y_hist_dynamical[index_12d(j,I_X,I_CWT,I_MK,I_CWV,I_VK,I_C1,I_C2,I_C3,I_CV1,I_CV2,I_CV3, this->NX,this->n_cwt,this->n_sknot_massbin,this->n_cwv,this->n_vknot_massbin,N_a,N_b,N_c, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3)];
                                         long lkk=get_max<long>(AUX);

                                         for(int j=0;j< ny_mf ;++j)
                                           X_Y_hist_dynamical_normalized[index_12d(j,I_X,I_CWT,I_MK,I_CWV,I_VK,I_C1,I_C2,I_C3,I_CV1,I_CV2,I_CV3, this->NX,this->n_cwt,this->n_sknot_massbin,this->n_cwv,this->n_vknot_massbin,N_a,N_b,N_c,N_CV_BIN1,N_CV_BIN2,N_CV_BIN3)]= (lkk == 0  ? 0.0 :  static_cast<real_prec>(AUX[j])/static_cast<real_prec>(lkk));

                                         AUX.clear();
                                         AUX.shrink_to_fit();

                                       }   // end if(N_cells > N_cells_now).

                                   }// end while(cell_accepted==false)
                               } // end of if(true==flag)
                             else  // If the density bin has been already filled, we call the distribution in that bin to assign randomly the value Nmass-BIN according to the other properties
                               {
                                 real_prec prob=-10.0;
                                 real_prec ran=10.0;
                                 ULONG index_or;
                                 while(prob<ran)
                                   {
                                     Nhalos_orfan= gsl_rng_uniform_int(r,ny_mf);
                                     index_or=index_12d(Nhalos_orfan,I_X,I_CWT,I_MK,I_CWV,I_VK, I_C1,I_C2, I_C3, I_CV1,I_CV2, I_CV3,this->NX,this->n_cwt, this->n_sknot_massbin, this->n_cwv, this->n_vknot_massbin,N_a,N_b, N_c, N_CV_BIN1,N_CV_BIN2, N_CV_BIN3);
                                     prob = this->ABUNDANCE_normalized[index_or];
                                     ran = gsl_rng_uniform(r);
                                   }
                                 real_prec xr = static_cast<real_prec>(gsl_rng_uniform(r));
                                 real_prec fraction_mass= xr*log10(this->tracer.MBmax[Nhalos_orfan]/this->tracer.MBmin[Nhalos_orfan]);
                                 real_prec lmass_halo = log10(this->tracer.MBmin[halo_prop])+fraction_mass ;
                                 this->tracer.Halo[i].mass = pow(10,lmass_halo);
                                 counter_orphan++;
                               }
                           }
                         else  // if no information in that particular theta-bin was in the reference, then assign using global mass function
                           {
                             real_prec prob=-10.;
                             real_prec newxr=10.;
                             real_prec mass_tracer;
                             while(prob<newxr)
                               {
                                 real_prec xr = static_cast<real_prec>(gsl_rng_uniform(r));
                                 real_prec fraction_mass = xr*log10(this->mass_max/this->mass_min);
                                 mass_tracer=pow(10,log10(this->mass_min)+fraction_mass);
                                 real_prec aux_mass= (mass_tracer <= this->tracer.MBmin[0]? this->tracer.MBmin[0]: mass_tracer);
                                 aux_mass= (mass_tracer >= this->tracer.MBmin[nbins_mf-1]? this->tracer.MBin[nbins_mf-1]: mass_tracer);
                                 if(aux_mass<=1.001*this->mass_min || aux_mass>=this->tracer.MBmin[nbins_mf-1])
                                   prob=0.0;
                                 else
                                   {
                                     newxr = static_cast<real_prec>(gsl_rng_uniform(r));
                                     prob  = static_cast<real_prec>(gsl_spline_eval (spline, aux_mass, acc));
                                   }
                              }
                             this->tracer.Halo[i].mass=mass_tracer;
                             counter_mf++;
                           }
                       }

       }// end loop over tracers
#ifdef _USE_OMP_
     gsl_rng_free (r);
   }// end parallelized region
#endif
   So.DONE();

   So.message_screen("Objects with masses following dynamical distribution =", counter_in);
#ifndef _MASS_ASSIGNMENT_TO_REFERENCE_
   So.message_screen("Objects with masses following full distribution n(M|{theta}) (orphans) =", counter_orphan);
   So.message_screen("Objects with masses following global n(M) mass function =", counter_mf);
#endif

   X_Y_hist_dynamical.clear();X_Y_hist_dynamical.shrink_to_fit();
   X_Y_hist_dynamical_normalized.clear();X_Y_hist_dynamical_normalized.shrink_to_fit();
   ABUNDANCE.clear();ABUNDANCE.shrink_to_fit();
   ABUNDANCE_normalized.clear();ABUNDANCE_normalized.shrink_to_fit();

   So.message_time_mock(start_mock);

   string fname_mass_function_Y = this->Output_directory+"tracer_mock_abundance.txt";
   this->tracer.get_property_function(fname_mass_function_Y);

   real_prec residuals=0;
   int aux_c=0;
#pragma omp parallel for reduction(+:residuals, aux_c)
   for(int i=0;i<this->params._NMASSbins_mf() ;++i)
     if(this->tracer.mass_function[i]!=0)
       {
         aux_c++;
         residuals+= fabs(this->tracer_ref.mass_function[i]/this->tracer.mass_function[i]-1.0);
         //        cout<<this->tracer_ref.mass_function[i]<<"  "<<this->tracer.mass_function[i]<<endl;
       }
   residuals/=static_cast<real_prec>(aux_c)/100.0;
   So.message_screen("Residuals from mass function = ", residuals, "%");
   cout<<endl;

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

 // This mas assignment explicitely uses the masses from the reference and assign them to the mocks provided the underlying properties of the DM field.
#ifdef MOCK_MODE
void Bam::assign_tracer_property_new_new(bool initial_assignment)
{
  // INitial assignmetn= false means that this functionis to be used after another primitive property has been alrady assigned


  So.message_screen("**************************************************");
#ifdef _USE_VMAX_AS_OBSERVABLE_
#ifdef _ASSIGN_MASS_POST_
  if(initial_assignment==true)
      So.message_screen("*Assigning Vmax from reference values to mock");
#endif

#elif defined _USE_MASS_AS_OBSERVABLE_
  So.message_screen("*Assigning Halo Mass from reference*");
#endif

  So.message_screen("Number of tracers in reference = ", this->tracer_ref.NOBJS);
  So.message_screen("Number of tracers in mock = ", this->tracer.NOBJS);
  So.message_screen("**************************************************");




  // ********************************************************************************************************************
  // *********************************************OMP stuff***********************************************************************
  int NTHREADS=omp_get_max_threads();
  //   So.message_screen("Using ",NTHREADS," threads");
  omp_set_num_threads(NTHREADS);
  // Define seed vector for each thread in the paralellized run
  gsl_rng_env_setup();

  const gsl_rng_type *Tn = gsl_rng_default;
  gsl_rng *rn = gsl_rng_alloc (Tn);
  gsl_rng_default_seed=155;

  // ********************************************************************************************************************
  ULONG LENGHTdm=this->dm_properties_bins.size();
  //vector de estructiras para guardar las masas que caen en un bin de {Theta}
  //This information is to be modified below, so we make a copy ni order not to measure that again above in get_X_function()


#ifndef _USE_MULTISCALE_MASS_ASSIGNMENT_
  vector<ULONG>number_in_theta_ref(LENGHTdm,0);
  // container to track the number of available masses in each theta-bin
#pragma omp parallel for
  for(ULONG i=0;i<LENGHTdm; ++i)
    number_in_theta_ref[i]=this->dm_properties_bins[i].masses_bin_properties.size();
  //   So.DONE();
#endif

  // ********************************************************************************************************************
  // ********************************************************************************************************************
  // *******************************Mass function stufff************************************************************************************
  int nbins_mf=0;
  real_prec lm_min, lm_max;
  gsl_interp_accel *acc;
  gsl_spline *spline ;

#ifdef _USE_VMAX_AS_OBSERVABLE_
      if(true==initial_assignment)
        nbins_mf=this->tracer_ref.vmax_function.size();
      else
        nbins_mf=this->tracer_ref.mass_function.size();
#else
      nbins_mf=this->tracer_ref.mass_function.size();
#endif

      if(nbins_mf==0)
        {
          So.message_warning("Reference mass function has not be allocated and perhaps not measured. Stop at line", __LINE__);
          exit(0);
        }

      So.message_screen("Preparing arrays for global mass function from reference: ");
      this->tracer.define_property_bins();
      this->mfunc.resize(nbins_mf,0);

#ifdef _USE_VMAX_AS_OBSERVABLE_
      if(true==initial_assignment)
        {
          this->mass_min=this->tracer_ref.VMAXBin[0];
          this->mass_max=this->tracer_ref.VMAXBin[nbins_mf-1];
        }
      else
        {
          this->mass_min=this->tracer_ref.MBin[0];
          this->mass_max=this->tracer_ref.MBin[nbins_mf-1];
        }
#else
      this->mass_min=this->tracer_ref.MBin[0];
      this->mass_max=this->tracer_ref.MBin[nbins_mf-1];
#endif

      // I cannot use here the this->tracer MBmin and max for it might happen that the reference mass function has been measured with a different number of bins
      // different to the current this>tracer_ref_NMBINS
      vector<real_prec>delta_mass_aux(nbins_mf,0);

      for(int i=0; i< nbins_mf; ++i)
        {
          real_prec lmmin=log10(mass_min)+i*log10(mass_max/mass_min)/static_cast<real_prec>(nbins_mf);
          real_prec lmmax=log10(mass_min)+(i+1)*log10(mass_max/mass_min)/static_cast<real_prec>(nbins_mf);
          delta_mass_aux[i]=pow(10,lmmax)-pow(10,lmmin);
        }

#ifdef _USE_VMAX_AS_OBSERVABLE_
      if(true==initial_assignment)
#pragma omp parallel for
        for(int i=0; i< nbins_mf; ++i)
          this->mfunc[i]=this->tracer_ref.vmax_function[i]*delta_mass_aux[i];
      else
#pragma omp parallel for
        for(int i=0; i< nbins_mf; ++i)
          this->mfunc[i]=this->tracer_ref.mass_function[i]*delta_mass_aux[i];
#else
#pragma omp parallel for
      for(int i=0; i< nbins_mf; ++i)
        this->mfunc[i]=this->tracer_ref.mass_function[i]*delta_mass_aux[i];
#endif


      delta_mass_aux.clear(); delta_mass_aux.shrink_to_fit();

      real_prec max_mf=static_cast<real_prec>(get_max(mfunc));
#pragma omp parallel for
      for(int i=0; i< nbins_mf; ++i)
        this->mfunc[i]/=max_mf;

      acc = gsl_interp_accel_alloc ();


#ifdef _USE_VMAX_AS_OBSERVABLE_
      if(true==initial_assignment)
        {
          spline    = gsl_spline_alloc (gsl_interp_linear, this->tracer_ref.vmax_function.size());
          gsl_spline_init (spline, &(this->tracer_ref.VMAXBin[0]), &(this->mfunc[0]), this->tracer_ref.vmax_function.size());
          lm_min=log10(this->params._VMAXmin());
          lm_max=log10(this->params._VMAXmax());
        }
      else{
        spline    = gsl_spline_alloc (gsl_interp_linear, this->tracer_ref.mass_function.size());
        gsl_spline_init (spline, &(this->tracer_ref.MBin[0]), &(this->mfunc[0]), this->tracer_ref.mass_function.size());
        lm_min=this->params._LOGMASSmin();
        lm_max=this->params._LOGMASSmax();
      }
#else
      spline    = gsl_spline_alloc (gsl_interp_linear, this->tracer_ref.mass_function.size());
      gsl_spline_init (spline, &(this->tracer_ref.MBin[0]), &(this->mfunc[0]), this->tracer_ref.mass_function.size());
      lm_min=this->params._LOGMASSmin();
      lm_max=this->params._LOGMASSmax();

#endif


  // ********************************************************************************************************************
  // ************************************************RENAME SOME BIN PROPPERTIES********************************************************************
  So.message_screen("*Assessing mock properties:");
  cout<<endl;
  // ********************************************************************************************************************

  int N_a=N_C_BIN1;
  vector<real_prec> MOCK_DEN_FIELD;
  if(true==initial_assignment)
    {
#ifdef _USE_TRACERS_IN_CELLS_
      MOCK_DEN_FIELD.resize(this->NGRID,0);
      this->tracer.get_density_field_grid(_COUNTS_, MOCK_DEN_FIELD);
      int nmax=get_max<real_prec>(MOCK_DEN_FIELD);
      So.message_screen("Maximum number of tracer in cell", nmax);
      if(nmax>=N_TRACERS_IN_CELLS_MAX){
        So.message_warning("Please increase the maximum number of tracer in one cell to " ,nmax);
        exit(0);
      }
      N_a = N_TRACERS_IN_CELLS_MAX;
#endif
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
  real_prec delta_min_sep=0;
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
  vector<real_prec> MOCK_MASS_FIELD;
  if(true==initial_assignment)  // Ww shall not use his info for post-mass assignment
    {
#ifdef _USE_TOTAL_MASS_IN_CELL_
      N_v = N_BINS_TOTAL_MASS_IN_CELL;

      // Can happen that the used mock properties fall in a theta bin in which the reference has nothing, hence  dm_properties_bins[index_bins].masses_bin_properties.size()=0
      // This has been visible when using a total mass density field sampled from the calibration (as it must be) instead of using the reference
      // This is likely to happen if we do not use the sampled total mass field that correspond to the final bumber density field that one uses to assign positions
      // Hence, if we are assigning as test_mode to the reference, use the mass density field from the reference

      MOCK_MASS_FIELD.resize(this->NGRID,0);

#ifdef _MASS_ASSIGNMENT_TO_REFERENCE_
      So.message_warning("I am passing the mass density field from the reference, in order toto do tests. Ideally this should be the one sampled in the last step of the calibration, LINE", __LINE__);
      this->tracer_ref.get_density_field_grid(_MASS_, MOCK_MASS_FIELD); // ojo que estoy cargando la ref, no el mock:esto debería ser el mock con la masa total dada por la calibración
#else
      File.read_array(this->Output_directory+"MOCK_TR_MASS_iteration200_MASY0_Nft256_z1.124.dat",MOCK_MASS_FIELD);
#endif

      real_prec mmax=get_max<real_prec>(MOCK_MASS_FIELD);
      So.message_screen("Maximum mass of tracer in mock-cell", mmax);
      mmax=get_min<real_prec>(MOCK_MASS_FIELD);
      So.message_screen("Minimum mass of tracer in mock-cell", mmax);
#endif
    }
  // ********************************************************************************

  int N_x= N_CV_BIN2;


#ifdef _ASSIGN_MASS_POST_
  if(false==initial_assignment) // if initial_assignem t is false, we proceed to assign mass once the vmax is already assigned
    N_x = N_VMAX_BINS;
#endif


  // ********************************************************************************************************************
  // ********************************************************************************************************************
  cout<<endl;
#ifdef _USE_MULTISCALE_MASS_ASSIGNMENT_
  vector<s_cell_info> cell_info_tr(this->NGRID);
#endif

  // If the number of neighbours is used, we do a loop over the particles.
  // If not, we do a loop over the grid cells, for the rest of quantities are already computed for the cells
  // When multiscale is not used, we have to do a loop over the particles here

#ifdef _USE_MULTISCALE_MASS_ASSIGNMENT_

#if defined _USE_NUMBER_OF_NEIGHBOURS_ || defined _ASSIGN_MASS_POST_
   for(ULONG ig=0; ig < this->tracer.NOBJS;++ig)
#else
#pragma omp parallel for
     for(ULONG id=0; id < this->NGRID;++id)
#endif

#else
       ULONG counter_fmf=0;
   for(ULONG ig=0; ig < this->tracer.NOBJS;++ig)
#endif
     {
       // Get the cell ID where the tracer is located
#ifdef _USE_MULTISCALE_MASS_ASSIGNMENT_
#if defined _USE_NUMBER_OF_NEIGHBOURS_ || defined _ASSIGN_MASS_POST_
       ULONG id=this->tracer.Halo[ig].GridID;
#endif
#else
       ULONG id=this->tracer.Halo[ig].GridID;
#endif

       // Get the bin in the delta dark matter (or log 1+delta) in each cell
       real_prec xdm = static_cast<real_prec>(this->delta_X[id]);

       int I_X  = get_bin(xdm,this->s_mins.prop1,this->NX,this->s_deltas.prop1,this->bin_accumulate_borders);

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
#ifdef _USE_TRACERS_IN_CELLS_
       if(true==initial_assignment)
          I_C1 = static_cast<int>(MOCK_DEN_FIELD[id]);
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
#elif defined _USE_TIDAL_ANISOTROPY_
       real_prec C3 = this->cwclass.Tidal_Anisotropy[id];
       I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined _USE_S2_
       real_prec C3 = this->cwclass.S2[id];             // s²
       I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#endif

       int I_CV1=0;
#ifdef _USE_TOTAL_MASS_IN_CELL_
       if(true==initial_assignment)
         I_CV1=get_bin(log10(MOCK_MASS_FIELD[id]),this->params._LOGMASSmin(),N_BINS_TOTAL_MASS_IN_CELL,(log10(MAX_TOTAL_MASS_IN_CELL)-this->params._LOGMASSmin())/static_cast<double>(N_BINS_TOTAL_MASS_IN_CELL), true);
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
#if defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_S2_)
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
#ifdef _USE_MULTISCALE_MASS_ASSIGNMENT_
                       cell_info_tr[id].Theta_bin=index_bins;     // Theta-bin in which the cell ID has been assigned
#endif

#ifndef _USE_MULTISCALE_MASS_ASSIGNMENT_
#ifdef  _USE_GLOBAL_MASS_FUNCTION_
#ifdef _MASS_ASSIGNMENT_TO_REFERENCE_
                       int N_masses_left=0;
#endif
#else
                       int N_masses_in_bin = this->dm_properties_bins[index_bins].masses_bin_properties.size();
                       int N_masses_left =  number_in_theta_ref[index_bins]; // this container is updated below
#endif

                       // If we have available masses in the theta-bin
                       if(N_masses_left>0) //when we use the ref as a mock, this condition will be always satisfied by construction
                         {
                           bool flag=false;
                           while(flag == false) // este while obliga a elegir una masa -no elegida antes- de entre las disponibles
                             {
                               int i_mass_halo_label= gsl_rng_uniform_int(rn,N_masses_in_bin);
                               bool used_mass = this->dm_properties_bins[index_bins].used_mass[i_mass_halo_label]; //true or false if the mass was already chosen or not
                               real_prec prop_to_assign = this->dm_properties_bins[index_bins].masses_bin_properties[i_mass_halo_label];
                               if(false==used_mass)// if the mass has not been used before, then assign that mass to the current particle i
                                 {
#ifdef _USE_VMAX_AS_OBSERVABLE_
                                   this->tracer.Halo[ig].vmax=prop_to_assign;
#elif defined _USE_MASS_AS_OBSERVABLE_
                                   this->tracer.Halo[ig].mass=prop_to_assign;
#endif
                                   this->tracer.Halo[ig].observed=true;
                                   counter_fmf++;
                                   this->dm_properties_bins[index_bins].used_mass[i_mass_halo_label] = true; //mark this mass as already assigned
                                   flag = true;
                                   number_in_theta_ref[index_bins]--;
                                 }
                             }
                         }
#endif
#ifndef _BIN_ACCUMULATE_
                     }  // end of ifs for theta-bins
#endif

#ifndef _USE_MULTISCALE_MASS_ASSIGNMENT_
       So.message_screen_flush("Number of masses assigned   = ",static_cast<int>(counter_fmf));
#endif

     }// end loop over tracers or grid


#ifndef _USE_MULTISCALE_MASS_ASSIGNMENT_
   number_in_theta_ref.clear();
   number_in_theta_ref.shrink_to_fit();
#endif

#ifdef _USE_TOTAL_MASS_IN_CELL_
   MOCK_MASS_FIELD.clear();
   MOCK_MASS_FIELD.shrink_to_fit();
#endif


#ifdef _USE_MULTISCALE_MASS_ASSIGNMENT_
   // *********************************** EXPLANATION NOTE*************************************************************************
   // Note that counter_high_mass and N_masses_multiscale are merely consequence of the mass function
   // These numbers can be computed from the refernece, but as we have assigned these to the mock, we get them from the mock here.
   // In the folliwnig loop we identify cells with obtjects above the mass threshold Mth, imposing olready positions for the masses
   // But we can avoid that threshold here, and then we would be counting all cells.
   // With the number of masses in each interval, we than can pass to assign the largest masses onto the low resolution mesh.
   // The fact that we assig the larger masses is in teh fact that we use the masses ordered in bins of theta and we incrase the index
   // everythime we assign, such that the nex assignment wiill use the next less massive obejec in the theta bin.
   // Regimes of multi-scale mass assignment:
   // i) M>M_multiscale_1, assign with NFFT_MS_1
   // ii) M_multiscale_2<M< M_multiscale_1 assign with NFFT_MS_2
   // iii) M_multiscale_3<M< M_multiscale_2 assign with NFFT_MS_3
   // iv) M< M_multiscale_3 done at the particle level, assigning the remaining masses to tracers according to the theta bin of their respective cell
   // Steps:
   // i) Count number of masses available in the different regimes
   ULONG N_masses_0=0;
   ULONG N_masses_1=0;
   ULONG N_masses_2=0;
   ULONG N_masses_3=0;
   ULONG N_masses_4=0;


   if(true==initial_assignment) // This if implies that mass assignment using vmax info is done at a particle level
     {

       real_prec factor_nobjects = 1; //static_cast<real_prec>(this->tracer.NOBJS)/static_cast<real_prec>(this->tracer_ref.NOBJS);
       So.message_screen("Correcting fraction of objects with a factor ",   factor_nobjects);

       // *******************************************************************************************************************************
       So.message_screen("Number of available masses in the different multi-scale levels obtained from the reference: ");
#ifdef _USE_MULTISCALE_LEVEL_4_
       N_masses_4=this->tracer_ref.N_masses_4; // number of masses in level 4   M > M_th4
       N_masses_4 = static_cast<ULONG>(floor(factor_nobjects*N_masses_4));
       So.message_screen("Number of masses in level 4 = ", N_masses_4);
       So.message_screen("Suggested grid size for this population = ", floor(pow(N_masses_4,1./3.)));
       So.message_screen("Using ", NFT_LOW_4);
       real_prec d1_l4=this->Lbox/static_cast<real_prec>(NFT_LOW_4);		/* grid spacing x-direction */
       vector<int>number_in_cells_aux_l4(NGRID_MS_4,0); // this would be for level 4
       cout<<endl;
       if(NFT_LOW_4>this->Nft)
         So.message_error("NFT for level 4 smaller than original value");
#endif

#ifdef _USE_MULTISCALE_LEVEL_3_
       N_masses_3=this->tracer_ref.N_masses_3; // number of masses in level 3   Mth3 < M < M_th4
       N_masses_3 = static_cast<ULONG>(floor(factor_nobjects*N_masses_3));
       So.message_screen("Number of masses in level 3 = ", N_masses_3);
       So.message_screen("Suggested grid size for this population = ", floor(pow(N_masses_3,1./3.)));
       So.message_screen("Using ", NFT_LOW_3);
       real_prec d1_l3=this->Lbox/static_cast<real_prec>(NFT_LOW_3);		/* grid spacing x-direction */
       vector<int>number_in_cells_aux_l3(NGRID_MS_3,0); // this would be for level 4
       cout<<endl;
       if(NFT_LOW_3>this->Nft)
         So.message_error("NFT for level 3 smaller than original value");
#endif

#ifdef _USE_MULTISCALE_LEVEL_2_
       N_masses_2=this->tracer_ref.N_masses_2; // number of masses in level 2   Mth2 < M < M_th3
       N_masses_2 = static_cast<ULONG>(floor(factor_nobjects*N_masses_2));
       So.message_screen("Number of masses in level 2 = ", N_masses_2);
       So.message_screen("Suggested grid size for this population = ", floor(pow(N_masses_2,1./3.)));
       So.message_screen("Using ", NFT_LOW_2);
       real_prec d1_l2=this->Lbox/static_cast<real_prec>(NFT_LOW_2);		/* grid spacing x-direction */
       vector<int>number_in_cells_aux_l2(NGRID_MS_2,0); // this would be for level 4
       cout<<endl;
       //     if(NFT_LOW_2>this->Nft)
       //     So.message_error("NFT for level 2 smaller than original value");
#endif

#ifdef _USE_MULTISCALE_LEVEL_1_
       N_masses_1=this->tracer_ref.N_masses_1; // number of masses in level 1   Mth1 < M < Mth2
       N_masses_1 = static_cast<ULONG>(floor(factor_nobjects*N_masses_1));
       So.message_screen("Number of masses in level 1 = ", N_masses_1);
       So.message_screen("Suggested grid size for this population = ", floor(pow(N_masses_1,1./3.)));
       So.message_screen("Using the orignal, ", this->Nft);
       cout<<endl;
#endif




       vector<int>number_in_cells_aux(this->NGRID,0); // this is needed regardless the level defined

       N_masses_0=this->tracer_ref.NOBJS-(N_masses_1+N_masses_2+N_masses_3+N_masses_4);

       So.message_screen("Total to be assigned from grid-approach = ", N_masses_1+N_masses_2+N_masses_3+N_masses_4);
       cout<<endl;
       So.message_screen("Number of reference masses to be assigned at particle level (level 0) = ",N_masses_0);
       cout<<endl;


       So.message_screen("Identifying galaxy id to get masses in the different regimes");

       // This loop is valid as long as masses have been assigned before, even if they will be reassigned.
       for(ULONG ig=0;ig < this->tracer.NOBJS; ++ig)
         {
           ULONG id=this->tracer.Halo[ig].GridID;
           this->tracer.Halo[ig].observed=false;  // this will be used below, se let keep it inzialized like that

           cell_info_tr[id].gal_index.push_back(ig);
           number_in_cells_aux[id]++;    // this is the number counts in the ref grid and will be updated in every level


#if defined _USE_MULTISCALE_LEVEL_4_ || defined _USE_MULTISCALE_LEVEL_3_ || defined _USE_MULTISCALE_LEVEL_2_
           ULONG il,jl,kl,id_low;
           real_prec x = this->tracer.Halo[ig].coord1;
           real_prec y = this->tracer.Halo[ig].coord2;
           real_prec z = this->tracer.Halo[ig].coord3;
#endif


           // Number of masses in level 4 (grid-based)
#ifdef _USE_MULTISCALE_LEVEL_4_
           il=static_cast<ULONG>(floor((x)/d1_l4));
           jl=static_cast<ULONG>(floor((y)/d1_l4));
           kl=static_cast<ULONG>(floor((z)/d1_l4));
           if(il==NFT_LOW_4)
             il--;
           if(jl==NFT_LOW_4)
             jl--;
           if(kl==NFT_LOW_4)
             kl--;
           id_low = index_3d(il,jl,kl,NFT_LOW_4,NFT_LOW_4);
           number_in_cells_aux_l4[id_low]++;
           this->tracer.Halo[ig].GridID_l4=id_low;
#endif



           // Number of masses in level 3 (grid-based)
#ifdef _USE_MULTISCALE_LEVEL_3_  // ensures level 1 is always used
           il=static_cast<ULONG>(floor((x)/d1_l3));
           jl=static_cast<ULONG>(floor((y)/d1_l3));
           kl=static_cast<ULONG>(floor((z)/d1_l3));
           if(il==NFT_LOW_3)
             il--;
           if(jl==NFT_LOW_3)
             jl--;
           if(kl==NFT_LOW_3)
             kl--;
           id_low = index_3d(il,jl,kl,NFT_LOW_3,NFT_LOW_3);
           number_in_cells_aux_l3[id_low]++;
           this->tracer.Halo[ig].GridID_l3=id_low;
#endif

           // Number of objects in level 2 (grid-based)
#ifdef _USE_MULTISCALE_LEVEL_2_  // ensures level 1 is always used
           il=static_cast<ULONG>(floor((x)/d1_l2));
           jl=static_cast<ULONG>(floor((y)/d1_l2));
           kl=static_cast<ULONG>(floor((z)/d1_l2));
           if(il==NFT_LOW_2)
             il--;
           if(jl==NFT_LOW_2)
             jl--;
           if(kl==NFT_LOW_2)
             kl--;
           id_low = index_3d(il,jl,kl,NFT_LOW_2,NFT_LOW_2);
           number_in_cells_aux_l2[id_low]++;
           this->tracer.Halo[ig].GridID_l2=id_low;
#endif
           // Number of objects in level 1: (grid- based)
         }
       So.DONE();

       // This follows the HADRON-approach, in which
       // we order masses in each theta*-bin from the max to the minimum and assign
       //them in that order to tracers located in randomly selected cells (within the same theta-bin)



#ifdef _USE_MASS_AS_OBSERVABLE_
       So.message_screen("Sorting masses in bins of theta");
#elif defined _USE_VMAX_AS_OBSERVABLE_
       So.message_screen("Sorting Vmax in bins of theta");
#endif

       So.message_screen("used in the multiscale approach: ");

       for(ULONG it=0; it < LENGHTdm; ++it)
         {
           int ll=this->dm_properties_bins[it].masses_bin_properties.size();
           if(ll>0)
             {
               gsl_vector *aux_masses;
               aux_masses=gsl_vector_alloc(ll);
               for(int im=0; im< ll; ++im )
                 gsl_vector_set(aux_masses,im,this->dm_properties_bins[it].masses_bin_properties[im]);
               gsl_sort_vector(aux_masses);
               for(int im=0; im< ll; ++im )
                 this->dm_properties_bins[it].masses_bin_properties[im]=gsl_vector_get(aux_masses,ll-im-1);
               gsl_vector_free(aux_masses);
             }
         }

       So.DONE();

       // This is indenepdent of the level, will be used for all in order to pinpoint the position of the ordered masses



       // ****************************************************************************************
       // ****************************************************************************************
       // ****************************************************************************************
       // ****************************************************************************************
       // ****************************************************************************************
       // ***************************************************multiscale Assignment *******************************************************

       vector<ULONG>count_aux_cell(LENGHTdm,0);  // This must be outside of the do-while loop
       vector<ULONG>cells_id_still_to_assign;
       ULONG counter_multiscale=0;
       bool reshuffle_once=true;

       // ****************************************************************************************
       // ****************************************************************************************
       // ***************************************************Level 4 *******************************************************


       //container to allocate the ID of the cells with particles still to get masses
       //This has to be clear, resized and updated for every level
#ifdef _USE_MULTISCALE_LEVEL_4_

       vector<s_cell_info>cell_info_cell(NGRID_MS_4);
       // This function returns the list of id in each ID_low (all cells included, regarless with or without mass)
       get_low_res_id_from_high_res_id(this->Nft, NFT_LOW_4,cell_info_cell);

       So.message_screen("Identifying level-4 empty cells");
       // THis is ment to be shiffled in order to run randomly the cells
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
       vector<bool>rejected_cell_lowres(NGRID_MS_4,false);
       reshuffle_once=true;
#endif

       ULONG N_cells=cells_id_still_to_assign.size();

       So.message_screen("Shuffling order of level-4 cells containing tracers without assigned masses");
       gsl_ran_shuffle(rn,&cells_id_still_to_assign[0], cells_id_still_to_assign.size() ,sizeof(ULONG));
       So.DONE();

       So.message_screen("Going through LEVEL 4 mesh to assign ", N_masses_4, "masses");

       counter_multiscale=0;

       int loop_counter=0;
       do{  // Assignment campaing, going through cells and assigning one tracer mass per cell


         // For that we select low_res cells with tracers without mass,
         //  and do a loop over them. For each cell (randomly ordered), we randomly select one of the original resolution cells, and pick up one tracer thereof.
         // We assign masses top-to-bottom untill we reach the number N_masses_multiscales, which corresponds to the mass threshold_multilscale
         // Once that first campaing is done, we start using the original mesh. We can use two campaings with the low-res cells, perhaps.
         // In this first part all links to ordered masses are ommited, since this works by consctruction under that definition.


         for(ULONG id_ini=0; id_ini< N_cells; ++id_ini) // loop over low-res cells containing tracers without assigned masses
           {

             ULONG id_low = cells_id_still_to_assign[id_ini];  // Low Grid ID. These come randomized

             int N_cells_in_cells= cell_info_cell[id_low].gal_index.size();  // Number of high-res cells in a low-res cell
             int index_cell_in_low_res_cell= gsl_rng_uniform_int(rn,N_cells_in_cells); //select randomy one of the high-id's in the low res cell
             ULONG id=cell_info_cell[id_low].gal_index[index_cell_in_low_res_cell]; // This chosen id will be used below to assign mass

             if(number_in_cells_aux[id]>0)
               {
                 ULONG index_bins = cell_info_tr[id].Theta_bin; // Get the theta-bin in which this cell has been classified. This has been calculated in the step M<M*


#ifdef _USE_EXCLUSION_CELL_LEVEL4_
                 bool used_previous_cell=false;  // ask whether the previos cell is one of the neighbours of the current one
                 bool rejected_previous_cell=false;   // to ask whether the previously id-1 cell was rejected because of being neighbour of the cell id-2
                 ULONG N_neighbour_cells=ID_nearest_cells_lowres[id_low].close_cell.size(); // Number of neighbours computed for the low_res

                 if(true==reshuffle_once)
                   if(id_ini > 0)
                     for(int in = 0 ; in < N_neighbour_cells; ++in)
                       if(ID_nearest_cells_lowres[id_low].close_cell[in]==cells_id_still_to_assign[id_ini-1]) /// this is satisfied only once (at most) inside the loop!
                         {
                           used_previous_cell=true;
                           break;
                         }
                       else  // this applies after the situation in which the cells are reshuffled and ow masses are assigned randomly reagardless whether the prevois cell was used or not. Forget the past
                         {
                           rejected_previous_cell=false;
                           used_previous_cell=false;
                         }
                 // Proceed to assign if the previous cell is not one of the closest neighbours of the current cell
                 // or if the previus  cell id-1 was rejected because being neighbour of the id-2

                 //               if(false==used_previous_cell || (true==rejected_previous_cell && true==used_previous_cell)  ) // if the previous cell was rejected, I can use this one wven if the current cell is neighbour of the previous one
                 if(false==used_previous_cell) // if the previous cell was rejected, I can use this one wven if the current cell is neighbour of the previous one
                   {
#endif  // end of use_exclusion in cells

                     ULONG N_masses_without_mass_in_cell=cell_info_tr[id].gal_index.size(); //if masses are to be assigned with this, this line can be MOCK_DEN_FIELD[id];
                     ULONG jk= gsl_rng_uniform_int(rn,N_masses_without_mass_in_cell); //choose a random integer in  [0,N_masses_in_cell_without_mass). N_masses_withopur mass must be the one fixed
                     ULONG ig=cell_info_tr[id].gal_index[jk];                       // get the Galaxy ID for the cell id in the position jk of tracers without mass
                     ULONG i_mass_halo_label=count_aux_cell[index_bins];
                     real_prec assigned_mass=this->dm_properties_bins[index_bins].masses_bin_properties[i_mass_halo_label];

                     if(false==this->tracer.Halo[ig].observed && assigned_mass>=this->Mass_threshold_multi_scale_4)
                       {
#ifdef _USE_MASS_AS_OBSERVABLE_
                         this->tracer.Halo[ig].mass=assigned_mass ;
#elif defined _USE_VMAX_AS_OBSERVABLE_
                         this->tracer.Halo[ig].vmax=assigned_mass ;
#endif
                         this->tracer.Halo[ig].observed=true;
                         this->dm_properties_bins[index_bins].used_mass[i_mass_halo_label]=true;
                         counter_multiscale++;
                         number_in_cells_aux[id]--; // Remove one tracer from the id cell
                         count_aux_cell[index_bins]++;  // This keeps track of the masses used, and will be used below for next-scale campaign
                       }
#ifdef _USE_EXCLUSION_CELL_LEVEL4_
                   }

                 else
                   rejected_cell_lowres[id_low]=true;
#endif

               } //end if number_in cells >0

             if(counter_multiscale==N_masses_4)
               break;

           }// end loop over cells

         So.message_screen_flush("Number of masses assigned = ",counter_multiscale);

#ifdef _USE_EXCLUSION_CELL_LEVEL4_
         if(reshuffle_once==true)
           if(counter_multiscale>= N_masses_4-N_MASSES_ABOVE_THRESHOLD_LEVEL4)  //
             {
               So.message_screen("Reshufling cells");
               gsl_ran_shuffle(rn,&cells_id_still_to_assign[0],cells_id_still_to_assign.size(),sizeof(ULONG));
               rejected_cell_lowres.resize(this->NGRID,false);
               reshuffle_once=false;
             }
#endif

         loop_counter++;

       }while(counter_multiscale< N_masses_4);

       cout<<endl;
       So.message_screen("Used ", loop_counter, "loops over L4 MESH to assign property");
       cout<<endl;

       So.DONE();
       So.message_screen("Freeing memmory");

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
       for(ULONG i=0; i<this->tracer.NOBJS;++i)
         if(true==this->tracer.Halo[i].observed)
           ncounts_aux[this->tracer.Halo[i].GridID_l4]++;
       So.message_screen("Maximum number of objects in cells from level-4 mesh = ", get_max<ULONG>(ncounts_aux));
       So.message_screen("Check: number of tracers = ", get_nobjects(ncounts_aux));
       /*
         vector<ULONG>ncounts_auxv(NGRID_MS_4,0);
         for(ULONG i=0; i<this->tracer_ref.NOBJS;++i)
         #ifdef _USE_MASS_AS_OBSERVABLE_
         if(this->tracer_ref.Halo[i].mass>= this->Mass_threshold_multi_scale_4)
         #elif defined _USE_VMAX_AS_OBSERVABLE_
         if(this->tracer_ref.Halo[i].vmax>= this->Mass_threshold_multi_scale_4)
         ncounts_auxv[this->tracer.Halo[i].GridID_l4]++;
         #endif


         So.message_screen("Maximum number of objects in cells from level-4 mesh from ref= ", get_max<ULONG>(ncounts_aux));
         So.message_screen("NUmber of tracers = ", get_nobjects(ncounts_aux));
         ncounts_aux.clear();ncounts_aux.shrink_to_fit();
         ncounts_auxv.clear();ncounts_auxv.shrink_to_fit();
       */

       cout<<endl;
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

       So.message_screen("Identifying level 3 empty cells");
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

       So.message_screen("Shuffling order of level 3 cells containing tracers without assigned masses");
       gsl_ran_shuffle(rn,&cells_id_still_to_assign[0], cells_id_still_to_assign.size() ,sizeof(ULONG));
       So.DONE();

#ifdef _USE_EXCLUSION_CELL_LEVEL3_
#ifdef _USE_EXCLUSION_CELL_LEVEL4_
       rejected_cell_lowres.resize(NGRID_MS_3,false);
#else
       vector<bool>rejected_cell_lowres(NGRID_MS_3,false);
#endif
#endif

       So.message_screen("Going through LEVEL 3 mesh to assign ", N_masses_3, "masses");

       counter_multiscale=0;

       do{  // Assignment campaing, going through cells and assigning one tracer mass per cell


         // For that we select low_res cells with tracers without mass,
         //  and do a loop over them. For each cell (randomly ordered), we randomly select one of the original resolution cells, and pick up one tracer thereof.
         // We assign masses top-to-bottom untill we reach the number N_masses_multiscales, which corresponds to the mass threshold_multilscale
         // Once that first campaing is done, we start using the original mesh. We can use two campaings with the low-res cells, perhaps.
         // In this first part all links to ordered masses are ommited, since this works by consctruction under that definition.
         ULONG N_cells=cells_id_still_to_assign.size();


         for(ULONG id_ini=0; id_ini< N_cells; ++id_ini) // loop over low-res cells containing tracers without assigned masses
           {
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
                     ULONG N_masses_without_mass_in_cell=cell_info_tr[id].gal_index.size(); //if masses are to be assigned with this, this line can be MOCK_DEN_FIELD[id];
                     ULONG jk= gsl_rng_uniform_int(rn,N_masses_without_mass_in_cell); //choose a random integer in  [0,N_masses_in_cell_without_mass). N_masses_withopur mass must be the one fixed
                     ULONG ig=cell_info_tr[id].gal_index[jk];                       // get the Galaxy ID for the cell id in the position jk of tracers without mass
                     if(false==this->tracer.Halo[ig].observed)
                       {
                         ULONG i_mass_halo_label=count_aux_cell[index_bins];
                         real_prec assigned_mass=this->dm_properties_bins[index_bins].masses_bin_properties[i_mass_halo_label];
#ifdef _USE_MULTISCALE_LEVEL_4_
                         if(assigned_mass >= this->Mass_threshold_multi_scale_3 && assigned_mass < this->Mass_threshold_multi_scale_4)
#elif !defined (_USE_MULTISCALE_LEVEL_4_)
                           if(assigned_mass >= this->Mass_threshold_multi_scale_3)
#endif
                             {

#ifdef _USE_MASS_AS_OBSERVABLE_
                               this->tracer.Halo[ig].mass=assigned_mass ;
#elif defined _USE_VMAX_AS_OBSERVABLE_
                               this->tracer.Halo[ig].vmax=assigned_mass ;
#endif
                               this->tracer.Halo[ig].observed=true;
                               this->dm_properties_bins[index_bins].used_mass[i_mass_halo_label]=true;
                               counter_multiscale++;
                               number_in_cells_aux[id]--; // Remove one tracer from the id cell
                               count_aux_cell[index_bins]++;  // This keeps track of the masses used, and will be used below for next-scale campaign
                             }
                       }

#ifdef _USE_EXCLUSION_CELL_LEVEL3_
                   }
                 else
                   rejected_cell_lowres[id_low]=true;
#endif

               } //end if number_in cells >0

             if(counter_multiscale==N_masses_3)
               break;

           }// end loop over cells

         So.message_screen_flush("Number of masses assigned = ",counter_multiscale);

       }while(counter_multiscale< N_masses_3);

       So.DONE();
       So.message_screen("Freeing memmory");
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

       So.message_screen("Identifying level 3 empty cells");
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

       So.message_screen("Shuffling order of level 2 cells containing tracers without assigned masses");
       gsl_ran_shuffle(rn,&cells_id_still_to_assign[0], cells_id_still_to_assign.size() ,sizeof(ULONG));
       So.DONE();


#ifdef _USE_EXCLUSION_CELL_LEVEL2_
#if defined  _USE_EXCLUSION_CELL_LEVEL3_ || defined (_USE_EXCLUSION_CELL_LEVEL4_)
       rejected_cell_lowres.resize(NGRID_MS_2,false);
#else
       vector<bool>rejected_cell_lowres(NGRID_MS_2,false);
#endif
#endif

       So.message_screen("Going through LEVEL 2 mesh to assign ", N_masses_2, "masses");

       counter_multiscale=0;

       do{  // Assignment campaing, going through cells and assigning one tracer mass per cell


         // For that we select low_res cells with tracers without mass,
         //  and do a loop over them. For each cell (randomly ordered), we randomly select one of the original resolution cells, and pick up one tracer thereof.
         // We assign masses top-to-bottom untill we reach the number N_masses_multiscales, which corresponds to the mass threshold_multilscale
         // Once that first campaing is done, we start using the original mesh. We can use two campaings with the low-res cells, perhaps.
         // In this first part all links to ordered masses are ommited, since this works by consctruction under that definition.
         ULONG N_cells=cells_id_still_to_assign.size();


         for(ULONG id_ini=0; id_ini< N_cells; ++id_ini) // loop over low-res cells containing tracers without assigned masses
           {
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
                       ULONG N_masses_without_mass_in_cell=cell_info_tr[id].gal_index.size(); //if masses are to be assigned with this, this line can be MOCK_DEN_FIELD[id];
                       ULONG jk= gsl_rng_uniform_int(rn,N_masses_without_mass_in_cell); //choose a random integer in  [0,N_masses_in_cell_without_mass). N_masses_withopur mass must be the one fixed
                       ULONG ig=cell_info_tr[id].gal_index[jk];                       // get the Galaxy ID for the cell id in the position jk of tracers without mass
                       if(false==this->tracer.Halo[ig].observed)
                         {
                           ULONG i_mass_halo_label=count_aux_cell[index_bins];
                           real_prec assigned_mass=this->dm_properties_bins[index_bins].masses_bin_properties[i_mass_halo_label];
#ifdef _USE_MULTISCALE_LEVEL_3_
                           if(assigned_mass >= this->Mass_threshold_multi_scale_2 && assigned_mass < this->Mass_threshold_multi_scale_3)
#elif !defined (_USE_MULTISCALE_LEVEL_3_) && defined (_USE_MULTISCALE_LEVEL_4_)
                             if(assigned_mass >= this->Mass_threshold_multi_scale_2 && assigned_mass < this->Mass_threshold_multi_scale_4)
#elif !defined (_USE_MULTISCALE_LEVEL_3_) && !defined (_USE_MULTISCALE_LEVEL_4_)
                               if(assigned_mass >= this->Mass_threshold_multi_scale_2)
#endif
                                 {

#ifdef _USE_MASS_AS_OBSERVABLE_
                                   this->tracer.Halo[ig].mass=assigned_mass ;
#elif defined _USE_VMAX_AS_OBSERVABLE_
                                   this->tracer.Halo[ig].vmax=assigned_mass ;
#endif
                                   this->tracer.Halo[ig].observed=true;
                                   this->dm_properties_bins[index_bins].used_mass[i_mass_halo_label]=true;
                                   counter_multiscale++;
                                   number_in_cells_aux[id]--; // Remove one tracer from the id cell
                                   count_aux_cell[index_bins]++;  // This keeps track of the masses used, and will be used below for next-scale campaign
                                 }
                         }

#ifdef _USE_EXCLUSION_CELL_LEVEL2_
                     }
                   else
                     rejected_cell_lowres[id_low]=true;
#endif

                 } //end if number_in cells >0

             if(counter_multiscale==N_masses_2)
               break;

           }// end loop over cells

         So.message_screen_flush("Number of masses assigned = ",counter_multiscale);

       }while(counter_multiscale< N_masses_2);

       So.DONE();
       So.message_screen("Freeing memmory");

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
       // ******************************************************************************************************************************************
       // ******************************************************************************************************************************************
       // ******************************Level 1 ****************************************************************************************************
#ifdef _USE_MULTISCALE_LEVEL_1_

       cells_id_still_to_assign.clear();
       cells_id_still_to_assign.shrink_to_fit();

       // *****************Assignment campaign with original reslution******************************************************************************
       //  Level 1 is always defined, unless _USE_MULTISCALE_MASS_ASSIGNMENT_ is undef
       // UPDATE info of cells with tracers after having done the low-res assignment.
       // Soms of the high res cells might be empty now


       So.message_screen("Updating cells with particles without assigned property");

       for(ULONG id=0; id<this->NGRID; ++id)
         if(number_in_cells_aux[id]>0)// Pinpoint the cells with tracers which do not have an assigned mass
           cells_id_still_to_assign.push_back(id);
       So.DONE();

       number_in_cells_aux.clear();
       number_in_cells_aux.shrink_to_fit();


       // We randomize here the entries of cells_id_still_to assign. This is important
       So.message_screen("Shuffling order of cells containing tracers without assigned property");
       gsl_ran_shuffle(rn,&cells_id_still_to_assign[0],cells_id_still_to_assign.size(),sizeof(ULONG));
       So.DONE();


#ifdef _USE_EXCLUSION_CELL_LEVEL1_
       So.message_screen("Identifying neighbouring cells for level1");
       vector<s_nearest_cells> ID_nearest_cells(this->NGRID);
       get_neighbour_cells(this->Nft,NUMBER_OF_NEIGHBOUR_CELL_FOR_EXCLUSION_LEVEL1,ID_nearest_cells);
       So.DONE();
       reshuffle_once=true;
#endif

       So.message_screen("Going through LEVEL 1 mesh to assign properties to ", N_masses_1, " tracers");

#ifdef  _USE_EXCLUSION_WITH_MASS_BINS_
       // Container allocating the number of objects in a id cell in a mass bin
       vector<int>mass_bins_exc(this->NGRID*N_MASS_BINS_EXC,0);
#endif


       counter_multiscale=0;


       do{


         for(ULONG id_ini=0; id_ini< cells_id_still_to_assign.size(); ++id_ini) // loop over cells containing tracers without assigned masses
           {

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
		 
                 ULONG jk= gsl_rng_uniform_int(rn,N_mocks_without_mass_in_cell); //choose a random integer in  [0,N_masses_in_cell_without_mass). N_masses_withopur mass must be the one fixed
		 
                 ULONG ig=cell_info_tr[id].gal_index[jk];        // get the Galaxy ID for the cell id in the position jk of tracers without mass
		 
                 ULONG i_mass_halo_label=count_aux_cell[index_bins];
                 if(false==this->tracer.Halo[ig].observed)
                   {

                  if(this->dm_properties_bins[index_bins].masses_bin_properties.size()>0)// this if must be added a theta bin filled from the mock might be empty in the ref
                    {

			 real_prec assigned_mass = this->dm_properties_bins[index_bins].masses_bin_properties[i_mass_halo_label];
#ifdef _USE_EXCLUSION_WITH_MASS_BINS_
			 ULONG index_mass_bin=get_bin(assigned_mass,log10(this->mass_min),N_MASS_BINS_EXC,log10(this->mass_max/this->mass_min)/static_cast<double>(N_MASS_BINS_EXC),true);
			 ULONG index_mass_id=index_2d(id, index_mass_bin,N_MASS_BINS_EXC);
			 if(mass_bins_exc[index_mass_id]< MAX_NUMBER_OF_MASSES_IN_CELL_SAME_MBIN) // allow only three masses from the same mass bin in a cell
			   {
#endif


#if defined (_USE_MULTISCALE_LEVEL_2_) && defined (_USE_MULTISCALE_LEVEL_3_) && defined (_USE_MULTISCALE_LEVEL_4_)
                         if(assigned_mass >= this->Mass_threshold_multi_scale_1 && assigned_mass < this->Mass_threshold_multi_scale_2)
#elif !defined (_USE_MULTISCALE_LEVEL_2_) && defined (_USE_MULTISCALE_LEVEL_3_) && defined (_USE_MULTISCALE_LEVEL_4_)
                         if(assigned_mass >= this->Mass_threshold_multi_scale_1 && assigned_mass < this->Mass_threshold_multi_scale_3)
#elif !defined (_USE_MULTISCALE_LEVEL_2_) && !defined (_USE_MULTISCALE_LEVEL_3_) && defined (_USE_MULTISCALE_LEVEL_4_)
                         if(assigned_mass >= this->Mass_threshold_multi_scale_1 && assigned_mass < this->Mass_threshold_multi_scale_4)
#elif !defined (_USE_MULTISCALE_LEVEL_2_) && !defined (_USE_MULTISCALE_LEVEL_3_) && !defined (_USE_MULTISCALE_LEVEL_4_)
                         if(assigned_mass >= this->Mass_threshold_multi_scale_1)
#endif
                           {
#ifdef _USE_MASS_AS_OBSERVABLE_
                              this->tracer.Halo[ig].mass=assigned_mass ;
#elif defined _USE_VMAX_AS_OBSERVABLE_
#ifdef _ASSIGN_MASS_POST_
                              if(true==initial_assignment)
#endif
                                this->tracer.Halo[ig].vmax=assigned_mass ;
#ifdef _ASSIGN_MASS_POST_
                              else
                                this->tracer.Halo[ig].mass=assigned_mass ;
#endif
#endif
                                   counter_multiscale++;
                                   count_aux_cell[index_bins]++;
                                   this->tracer.Halo[ig].observed=true;
                                   this->dm_properties_bins[index_bins].used_mass[i_mass_halo_label] = true;                                   //mark this mass as already assigned
#ifdef _USE_EXCLUSION_WITH_MASS_BINS_
                                   mass_bins_exc[index_mass_id]++;
#endif
                                 }

#ifdef _USE_EXCLUSION_WITH_MASS_BINS_
                       }
#endif
                 }// if no available particles in this theta bin, use the global mass function
#ifndef _MASS_ASSIGNMENT_TO_REFERENCE_
                else // if no particles available, use the global mass function
                  {
                    real_prec prob=-10.0;
                    real_prec ran=10.0;
                    real_prec mass_tracer;
                    while(prob<ran)
                      {
                         real_prec xr = static_cast<real_prec>(gsl_rng_uniform(rn));
                         real_prec fraction_mass = xr*log10(this->mass_max/this->mass_min);
                         mass_tracer=pow(10,log10(this->mass_min)+fraction_mass);
#ifdef _USE_MASS_AS_OBSERVABLE_
                         real_prec aux_mass= (mass_tracer <= this->tracer_ref.MBmin[0]? this->tracer_ref.MBmin[0]: mass_tracer);
                         aux_mass= (mass_tracer >= this->tracer_ref.MBmin[nbins_mf-1]? this->tracer.MBin[nbins_mf-1]: mass_tracer);
#elif defined _USE_VMAX_AS_OBSERVABLE_
                         real_prec aux_mass= (mass_tracer <= this->tracer_ref.VMAXBmin[0]? this->tracer_ref.VMAXBmin[0]: mass_tracer);
                         aux_mass= (mass_tracer >= this->tracer_ref.VMAXBmin[nbins_mf-1]? this->tracer.VMAXBin[nbins_mf-1]: mass_tracer);
#endif

                         if(aux_mass<this->Mass_threshold_multi_scale_1)
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
                      counter_multiscale++;
                 }
#endif   // END IFNDEF_MASS_ASSIGNMENT_TO_REFERENCE_
          }// end if false==observed



#ifdef _USE_EXCLUSION_CELL_LEVEL1_
               }
#endif

             if(counter_multiscale==N_masses_1)  // comment when the do-while is also commented
               break;
           }// end loop over cells

         So.message_screen_flush("Number of masses assigned   = ",static_cast<int>(counter_multiscale));

#ifdef _USE_EXCLUSION_CELL_LEVEL1_
         if(reshuffle_once==true)
           if(counter_multiscale>= N_masses_1-N_MASSES_ABOVE_THRESHOLD)
             {
               gsl_ran_shuffle(rn,&cells_id_still_to_assign[0],cells_id_still_to_assign.size(),sizeof(ULONG));
               //	     rejected_cell.resize(this->NGRID,false);
               reshuffle_once=false;
             }
#endif

#ifdef _USE_EXCLUSION_WITH_MASS_BINS_
         if(counter_multiscale>=N_masses_1-N_MASSES_ABOVE_M_EX)
           break;
#endif
       }while(counter_multiscale<N_masses_1);   // do untill all masses have acquisted mass

#ifdef _USE_EXCLUSION_WITH_MASS_BINS_
       N_masses_0=this->tracer_ref.NOBJS-((N_masses_1-N_MASSES_ABOVE_M_EX)+N_masses_2+N_masses_3+N_masses_4);
#endif
       cout<<endl;

       So.message_screen("Freeing memory");
       cells_id_still_to_assign.clear();
       cells_id_still_to_assign.shrink_to_fit();
       So.DONE();
#endif   // end of level 1

              

#ifndef _MASS_ASSIGNMENT_TO_REFERENCE_
       So.message_screen("Freeing memory");
       this->dm_properties_bins.clear(); this->dm_properties_bins.shrink_to_fit();
       So.DONE();
#endif

       // **********************************************Level 0 *********************************************************
       // ***************************************Assign at particle level for X<X threshold *****************************
       cout<<endl;
       cout<<endl;
#ifdef _USE_VMAX_AS_OBSERVABLE_
       So.message_screen("Assigning Vmax at particle level");

#ifdef _MASS_ASSIGNMENT_TO_REFERENCE_
       So.message_screen("Available = ", N_masses_0);
#endif

#else
       So.message_screen("Assigning Mhalo at particle level. Available = ", N_masses_0);
#endif

       ULONG counter_masses_0=0;
       ULONG counter_masses_b=0;

       for(ULONG ig=0; ig < this->tracer.NOBJS;++ig)
         if(false==this->tracer.Halo[ig].observed)
           counter_masses_b++;

       So.message_screen("Requested = ", counter_masses_b);
       So.message_screen("Needed = ", this->tracer.NOBJS-N_masses_1);

       int bin_threshold=get_bin(log10(this->Mass_threshold_multi_scale_1), log10(this->params._VMAXmin()),this->params._NMASSbins(),log10(this->params._VMAXmax()/this->params._VMAXmin())/static_cast<double>(this->params._NMASSbins()), false);



       counter_masses_0=0;
       
#ifdef _MASS_ASSIGNMENT_TO_REFERENCE_
       if(N_masses_0>0)
#else
       if(counter_masses_b>0)
#endif
           {
           for(ULONG ig=0; ig < this->tracer.NOBJS;++ig)
             {
	       
               if(false==this->tracer.Halo[ig].observed)
                 {
                   ULONG id=this->tracer.Halo[ig].GridID;
                   ULONG index_bins = cell_info_tr[id].Theta_bin;
		   
		   
#ifdef _MASS_ASSIGNMENT_TO_REFERENCE_
                   ULONG N_masses_in_bin = this->dm_properties_bins[index_bins].masses_bin_properties.size();
                   bool flag=false;
                   while(flag == false) // este while obliga a elegir una masa -no elegida antes- de entre las disponibles
		     {
		       int i_mass_halo_label= gsl_rng_uniform_int(rn,N_masses_in_bin);
		       real_prec assigned_mass=dm_properties_bins[index_bins].masses_bin_properties[i_mass_halo_label] ;
		       bool used_mass = this->dm_properties_bins[index_bins].used_mass[i_mass_halo_label]; //true or false if the mass was already chosen or not
		       if(false==used_mass)// if the mass has not been used before, then assign that mass to the current particle i
			 {
#ifdef _USE_MASS_AS_OBSERVABLE_
			   this->tracer.Halo[ig].mass=assigned_mass ;
#elif defined _USE_VMAX_AS_OBSERVABLE_
			   this->tracer.Halo[ig].vmax=assigned_mass ;
#endif
			   this->tracer.Halo[ig].observed=true;
			   this->dm_properties_bins[index_bins].used_mass[i_mass_halo_label] = true; //mark this mass as already assigned
			   flag = true;
			   counter_masses_0++;
			   
			 }
		     }  // close while
#else
		   
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
			   if(mass_assigned<this->Mass_threshold_multi_scale_1) // just to be sure
			     {
			       this->tracer.Halo[ig].mass = mass_assigned;
			       this->tracer.Halo[ig].observed=true;
			       counter_masses_0++;
			       flag=true;
			     }
#elif  defined _USE_VMAX_AS_OBSERVABLE_
			   real_prec fraction_mass= xr*log10(this->tracer_ref.VMAXBmax[i_mass_halo_index]/this->tracer_ref.VMAXBmin[i_mass_halo_index]);
			   real_prec lmass_halo = log10(this->tracer_ref.VMAXBmin[i_mass_halo_index])+fraction_mass ;
			   real_prec mass_assigned=pow(10,lmass_halo);
			   if(mass_assigned<this->Mass_threshold_multi_scale_1) // just to be sure
			     {
			       this->tracer.Halo[ig].vmax = mass_assigned;
			       this->tracer.Halo[ig].observed=true;
			       counter_masses_0++;
			       flag=true;
			     }
#endif
			 }
		       
		     }
#endif
		 }
	       
#ifdef _USE_MASS_AS_OBSERVABLE_
               So.message_screen_flush("Number of masses assigned   = ",static_cast<int>(counter_masses_0));
#elif defined _USE_VMAX_AS_OBSERVABLE_
               So.message_screen_flush("Number of Vmax assigned   = ",static_cast<int>(counter_masses_0));
#endif
	       
             }// end loop over tracers
           So.DONE();
           N_masses_0=counter_masses_0;  // update this number. Asked below

           cout<<endl;
         } // end if N_masses_0 > 0

       // **********************************************END Level 0 *********************************************************
#endif // end for use__multiscale_approach


   }  // close if(true==initial_assignment): This implies that mass assignment as post-processing (i.e, using Vmax assigned info) is done only at a particle-bases
      
   else //this follows initial_assignment=false, to assign MASSES based on the information of  v max already assigned
     {
       N_masses_0=this->tracer.NOBJS;
       So.message_screen("Assigning masses at particle level. Expected = ", this->tracer.NOBJS);
       
#ifndef test_vmax_mass
       vector<ULONG>number_in_theta_ref(LENGHTdm,0);
       // container to track the number of available masses in each theta-bin
#pragma omp parallel for
       for(ULONG i=0;i<LENGHTdm; ++i)
         number_in_theta_ref[i]=this->dm_properties_bins[i].masses_bin_properties.size();
#endif

       ULONG counter_masses_0=0;
       ULONG counter_masses_or=0;

      ULONG Ntot=N_VMAX_BINS;
#ifdef _add_dm_density_
      Ntot*=this->NX;
#endif

#pragma omp parallel for
      for(ULONG ig=0; ig < this->tracer.NOBJS;++ig)
        {
          this->tracer.Halo[ig].observed=false;
          this->tracer.Halo[ig].mass=0;
        }



      for(ULONG ig=0; ig < this->tracer.NOBJS;++ig)
        {

#if !defined test_vmax_mass || defined (_add_dm_density_)
          ULONG id=this->tracer.Halo[ig].GridID;
          real_prec xdm = static_cast<real_prec>(this->delta_X[id]);
          int I_X  = get_bin(xdm,this->s_mins.prop1,this->NX,this->s_deltas.prop1,this->bin_accumulate_borders);
#endif


#ifndef test_vmax_mass
          ULONG id=this->tracer.Halo[ig].GridID;
          ULONG index_bins = cell_info_tr[id].Theta_bin;
          ULONG N_masses_left =  number_in_theta_ref[index_bins];
          ULONG N_masses_in_bin = this->dm_properties_bins[index_bins].masses_bin_properties.size();
          if(N_masses_left>0)
            {
              bool flag=false;
              while(flag == false) // este while obliga a elegir una masa -no elegida antes- de entre las disponibles
                {
                  int i_mass_halo_label= gsl_rng_uniform_int(rn,N_masses_in_bin);
                  bool used_mass = this->dm_properties_bins[index_bins].used_mass[i_mass_halo_label]; //true or false if the mass was already chosen or not
                  if(false==used_mass)// if the mass has not been used before, then assign that mass to the current particle i
                    {
                      this->tracer.Halo[ig].mass=dm_properties_bins[index_bins].masses_bin_properties[i_mass_halo_label] ;
                      this->tracer.Halo[ig].observed=true;
                      this->dm_properties_bins[index_bins].used_mass[i_mass_halo_label] = true; //mark this mass as already assigned
                      flag = true;
                      number_in_theta_ref[index_bins]--;
                      counter_masses_0++;
                    }
                }
            }
          else
            {
#endif

#ifdef test_vmax_mass
              ULONG index_bins=get_bin(log10(this->tracer.Halo[ig].vmax), log10(this->params._VMAXmin()),N_VMAX_BINS,log10(this->params._VMAXmax()/this->params._VMAXmin())/static_cast<double>(N_VMAX_BINS), true);
#endif
              real_prec aux_h=0;

              ULONG index_dm;
#ifdef _add_dm_density_
              index_dm=index_2d(index_bins,I_X,this->NX);
#else
              index_dm=index_bins;
#endif


              for(int iy=0;iy<this->params._NMASSbins();++iy)
#ifdef test_vmax_mass
                aux_h+=this->ABUNDANCE_normalized[index_2d(iy,index_dm,Ntot)];
#else
                aux_h+=this->ABUNDANCE_normalized[index_2d(iy,index_bins,LENGHTdm)];
#endif
              if(aux_h>0)
                {
                  real_prec prob=-10.0;
                  real_prec ran=10.0;
                  int i_mass_halo_index;
                  while(prob<ran)
                    {
                      i_mass_halo_index= gsl_rng_uniform_int(rn,this->params._NMASSbins());
#ifdef test_vmax_mass
                      ULONG index_or=index_2d(i_mass_halo_index,index_dm,Ntot);
#else
                       ULONG index_or=index_2d(i_mass_halo_index,index_bins,LENGHTdm);
#endif

                       prob = this->ABUNDANCE_normalized[index_or];
                       ran = gsl_rng_uniform(rn);
                     }
                   real_prec xr = static_cast<real_prec>(gsl_rng_uniform(rn));
                   real_prec fraction_mass= xr*log10(this->tracer_ref.MBmax[i_mass_halo_index]/this->tracer_ref.MBmin[i_mass_halo_index]);
                   real_prec lmass_halo = log10(this->tracer_ref.MBmin[i_mass_halo_index])+fraction_mass ;
                   this->tracer.Halo[ig].mass = pow(10,lmass_halo);
                   this->tracer.Halo[ig].observed=true;
                   counter_masses_0++;
                   counter_masses_or++;
                   // here we assign mass according to the P(M|Vmax,ð) relation obtained from the reference
                 }
#ifndef _MASS_ASSIGNMENT_TO_REFERENCE_
              else{

                      real_prec prob=-10.0;
                      real_prec ran=10.0;
                      real_prec mass_tracer;
                      while(prob<ran)
                        {
                           real_prec xr = static_cast<real_prec>(gsl_rng_uniform(rn));
                           real_prec fraction_mass = xr*log10(this->mass_max/this->mass_min);
                           mass_tracer=pow(10,log10(this->mass_min)+fraction_mass);
                           real_prec aux_mass= (mass_tracer <= this->tracer_ref.MBmin[0]? this->tracer_ref.MBmin[0]: mass_tracer);
                           aux_mass= (mass_tracer >= this->tracer_ref.MBmin[nbins_mf-1]? this->tracer_ref.MBin[nbins_mf-1]: mass_tracer);
                           if(aux_mass<this->tracer_ref.MBmin[0] || aux_mass >= this->tracer_ref.MBmin[nbins_mf-1])
                             prob=0.0;
                           else
                             {
                               ran  = static_cast<real_prec>(gsl_rng_uniform(rn));
                               prob = static_cast<real_prec>(gsl_spline_eval (spline, aux_mass, acc));
                             }
                           mass_tracer=aux_mass;
                         }
                        this->tracer.Halo[ig].mass=mass_tracer;
                        counter_masses_0++;
                  }
#endif



#ifndef test_vmax_mass
            }
#endif

#ifdef _USE_MASS_AS_OBSERVABLE_
           So.message_screen_flush("Number of masses assigned   = ",static_cast<int>(counter_masses_0));
#elif defined _USE_VMAX_AS_OBSERVABLE_
            So.message_screen_flush("Number of masses assigned   = ",static_cast<int>(counter_masses_0));
#endif
         } // end loop over tracers
       cout<<endl;
       So.message_screen("Number of masses assigned (using prob-dist)   = ",static_cast<int>(counter_masses_or));

       So.DONE();
       cout<<endl;

#ifndef test_vmax_mass
       number_in_theta_ref.clear();
       number_in_theta_ref.shrink_to_fit();
#endif

       N_masses_0=counter_masses_0;
     }




   So.message_screen("Freeing memmory acá");
   this->ABUNDANCE_normalized.clear(); this->ABUNDANCE_normalized.shrink_to_fit();
   So.DONE();


#ifndef test_vmax_mass
   this->dm_properties_bins.clear();
   this->dm_properties_bins.shrink_to_fit();
#endif


#ifdef _USE_TRACERS_IN_CELLS_
   MOCK_DEN_FIELD.clear();
   MOCK_DEN_FIELD.shrink_to_fit();
#endif
   So.DONE();


#ifdef _USE_MULTISCALE_MASS_ASSIGNMENT_ 
    ULONG Assigned_masses_mf=N_masses_0+N_masses_1+N_masses_2+N_masses_3+N_masses_4;
#else
   ULONG Assigned_masses_mf=counter_fmf;
#endif

   if(true==initial_assignment)  //this applies for vmax assignment when the global vmax function  is needed
     {
       if(Assigned_masses_mf  < this->tracer.NOBJS)
         {
           So.message_screen("Number of tracers with no property assigned = ", this->tracer.NOBJS-Assigned_masses_mf);
           cout<<endl;
           So.message_screen("Assigning mass with global *property* function");

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
                       real_prec fraction_mass = xr*log10(this->mass_max/this->mass_min);
                       mass_tracer=pow(10,log10(this->mass_min)+fraction_mass);
                       real_prec aux_mass= (mass_tracer <= this->tracer_ref.VMAXBmin[0]? this->tracer_ref.VMAXBmin[0]: mass_tracer);
                       aux_mass= (mass_tracer >= this->tracer_ref.VMAXBmin[nbins_mf-1]? this->tracer.VMAXBin[nbins_mf-1]: mass_tracer);
#ifdef _MASS_ASSIGNMENT_TO_REFERENCE_
                       // When assigned to mocks, this section is complementary to that assigning below the threshold
                       if(aux_mass<this->tracer_ref.VMAXBmin[0] || aux_mass>=this->tracer_ref.VMAXBmin[nbins_mf-1])
#else
                       if(aux_mass<this->tracer_ref.VMAXBmin[0] || aux_mass>=this->Mass_threshold_multi_scale_1)
#endif
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
                 }
           }
           So.DONE();
         }
    }
#ifndef test_vmax_mass
   else  // if false==initial_assignment, so this assings mass
    { // global mass assignment once vmax was done
      if(Assigned_masses_mf  < this->tracer.NOBJS)
       {

         So.message_screen("Number of tracers with no property assigned = ", this->tracer.NOBJS-Assigned_masses_mf);
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
                     real_prec fraction_mass = xr*log10(this->mass_max/this->mass_min);
                     mass_tracer=pow(10,log10(this->mass_min)+fraction_mass);
                     real_prec aux_mass= (mass_tracer <= this->mass_min ? this->mass_min : mass_tracer);
                     aux_mass= (mass_tracer >= this->mass_max ? this->mass_max: mass_tracer);
                     if(aux_mass<  this->mass_min  || aux_mass>=this->mass_max)
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
//             So.message_screen_flush("Number of masses assigned   = ",static_cast<int>(counter_m));
         }
         So.DONE();
       }
   }
#endif


   gsl_rng_free (rn);





   NTHREADS=omp_get_max_threads();
   omp_set_num_threads(NTHREADS);


#ifdef _USE_VMAX_AS_OBSERVABLE_
   string fname_mass_function_Y = this->Output_directory+"tracer_mock_abundance.txt";

   if(true==initial_assignment)
     {
       this->tracer.params.i_vmax_g=5; //this allows the tracer to get the vmax function
       this->tracer.params.i_mass_g=-5; //this allows the tracer to get the vmax function
       this->tracer.get_property_function(fname_mass_function_Y);

       real_prec residuals=0;
       int count_b=0;
       for(int i=0;i<this->params._NMASSbins_mf() ;++i)
         if(this->tracer.vmax_function[i]>0)
           {
             count_b++;
             residuals+= fabs(this->tracer_ref.vmax_function[i]/this->tracer.vmax_function[i] -1.0 );
           }
       residuals/=static_cast<real_prec>(count_b)/100.0;
       So.message_screen("Residuals from  vmax-function = ", residuals, "%");
       cout<<endl;
     }
   else
     {
       this->tracer.params.i_mass_g=4; //this allow the tracer to get the mass function
       this->tracer.get_property_function(fname_mass_function_Y);
       real_prec residuals=0;
       int count_b=0;
       for(int i=0;i<this->params._NMASSbins_mf() ;++i)
         if(this->tracer.mass_function[i]>0)
           {
             count_b++;
             residuals+= fabs(this->tracer_ref.mass_function[i]/this->tracer.mass_function[i]-1.0);
           }
       residuals/=static_cast<real_prec>(count_b)/100.0;
       So.message_screen("Residuals from mass-function = ", residuals, "%");
       cout<<endl;
     }
#elif defined _USE_MASS_AS_OBSERVABLE_
   real_prec residuals=0;
   for(int i=0;i<this->params._NMASSbins_mf() ;++i)
     if(this->tracer.mass_function[i]>0)
       residuals+= fabs(this->tracer_ref.mass_function[i]/this->tracer.mass_function[i]);
   residuals/=static_cast<real_prec>(this->params._NMASSbins_mf())/100.0;
   So.message_screen("Residuals from mass-function = ", residuals, "%");
   cout<<endl;
#endif


}

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
       pname= "_realization"+to_string(this->step - (this->N_iterations_Kernel)+this->N_dm_initial-1);
     string fname=this->Output_directory+"MOCK_TR"+pname+"_CAT_z"+to_string(this->redshift)+".txt";

     vector<real_prec>x_cart;
     vector<real_prec>y_cart;
     vector<real_prec>z_cart;

     real_prec delta=this->Lbox/static_cast<real_prec>(this->Nft);
     for(ULONG i=0;i<this->delta_Y_new.size();++i)
       {
         int Ngal=this->delta_Y_new[i];
         int counter=0;
         real_prec x_min, y_min, z_min;
         index2coords(this->Nft, i, &x_min, &y_min, &z_min);
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

#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_INVARIANT_TIDAL_FIELD_III_)
      this->cwclass.do_CWC(this->delta_X_ini);
#endif

#endif

#if (defined (_USE_NABLA2DELTA_) || defined (_USE_S2_)  || defined (_USE_S3_) || defined (_USE_S2DELTA_)) && (!defined (_USE_CWC_))
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

         string file=this->Output_directory+"2dbias"+prop_X+"_"+prop_Y+"_Nft"+to_string(this->Nft)+"_z"+to_string(this->redshift)+"_LambdaTH"+to_string(lambdath)+"_CW"+to_string(this->tstruct);

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
             So.message_screen("Transforming to overdensities for CWT",this->tstruct);


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
         So.message_screen("Getting X and Y bins used for the bias information");
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
               this->new_nbins_y = this->NY;
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
               this->new_nbins_y = this->NY;
               this->Ymin=this->ldelta_Y_min;
               if(type_Y=="linear")
                 this->Ymax=100;
               else if(type_Y=="log")
                 this->Ymax=this->ldelta_Y_max;
             }
         }
         So.DONE();



         // ******************************************************************************

         this->DELTAX=(this->Xmax-this->Xmin)/static_cast<real_prec>(this->NX);
         this->DELTAY=(this->Ymax-this->Ymin)/static_cast<real_prec>(this->new_nbins_y);
         cout<<DELTAX<<endl;
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
         for(ULONG i = 0; i < new_nbins_y; ++i)
           Y_bins[i]= (iMAS_Y == 0  && false==this->Convert_Density_to_Delta_Y  && type_Y=="linear") ? static_cast<real_prec>(i) :  this->Ymin+(i+0.5)*(this->Ymax-this->Ymin)/static_cast<real_prec>(new_nbins_y);

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

         So.message_screen("Computing BIAS(X,Y, CWT, ...) ");

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
                               //int IX  = (iMAS_X == 0  && false==this->Convert_Density_to_Delta_X && type_X=="linear") ? static_cast<int>(this->delta_X[ig]) : static_cast<int>(floor((delta_X[ig] - this->Xmin)/this->DELTAX));
                               int IX= get_bin(delta_X[ig], this->Xmin, this->NX,this->DELTAX, this->bin_accumulate_borders);
                               int IY= get_bin(delta_Y[ig], this->Ymin, this->NY,this->DELTAY, this->bin_accumulate_borders);

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

                 So.message_screen("Writting the joint distribution P(X,Y)");
                 vector<vector<real_prec> > Vaux;
                 Vaux.resize(new_nbins_x);
                 for(int i=0;i<new_nbins_x;++i)
                   for(int j=0;j<new_nbins_y;++j)
                     Vaux[i].push_back(aux_n[index_2d(i,j,new_nbins_y)]);

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
                 aux_v.clear(); aux_v.shrink_to_fit();
                 aux_n.clear(); aux_n.shrink_to_fit();
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
                                       int IX=get_bin(delta_X[i],this->Xmin,this->new_nbins_x,DELTAX,this->bin_accumulate_borders);
                                       DELTAXY[IX].push_back(delta_Y[i]);// IX ES EL BIN DE X: allocate all values of Y in the bins of dm.
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

   }



   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
#ifdef MOCK_MODE
   void Bam::collapse_randoms()
   {
     cout<<endl;
     this->So.message_screen("**Collapsing randoms towards the DM particles**");
     cout<<endl;


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

     this->So.message_screen("Separating DM and random:");
     ULONG N_dms=0;
     ULONG N_rand=0;

     vector<real_prec> x_dm_pos;
     vector<real_prec> y_dm_pos;
     vector<real_prec> z_dm_pos;
     vector<real_prec> mass_dm;

     vector<ULONG> dm_id;

     vector<real_prec> x_random_pos;
     vector<real_prec> y_random_pos;
     vector<real_prec> z_random_pos;
     vector<real_prec> mass_random;

     vector<ULONG> ran_id;

     for(ULONG i=0; i< Ntracers; ++i)
       {
         ULONG id=grid_ID(&box_collps, this->tracer.Halo[i].coord1,this->tracer.Halo[i].coord2,this->tracer.Halo[i].coord3);
         if(this->tracer.Halo[i].identity>0)
           {
             x_dm_pos.push_back(this->tracer.Halo[i].coord1);
             y_dm_pos.push_back(this->tracer.Halo[i].coord2);
             z_dm_pos.push_back(this->tracer.Halo[i].coord3);
             mass_dm.push_back(this->tracer.Halo[i].mass);
             dm_id.push_back(id);
             dm_count[id]++;
             N_dms++;
           }
         else
           {
             x_random_pos.push_back(this->tracer.Halo[i].coord1);
             y_random_pos.push_back(this->tracer.Halo[i].coord2);
             z_random_pos.push_back(this->tracer.Halo[i].coord3);
             mass_random.push_back(this->tracer.Halo[i].mass);
             ran_id.push_back(id);
             N_rand++;
           }
       }

     this->So.DONE();

     int max_count=get_max(dm_count); // get the maximum number of dm particles in one cell.
     So.message_screen("Maximum number of DM particles found in one cell (low resolution) = ", max_count);

     // Initialize this to NO_NUM, a negative number, for not all elements of the container
     // dm_index_cell will be filled and initializing with 0 makes confusion, for 0 is a valid/used entry


#define NO_NUM -999
     vector<vector<int> > dm_index_cell(box_collps.NGRID, vector<int>(max_count,NO_NUM));

     dm_count.clear();
     dm_count.resize(box_collps.NGRID,0);
     So.message_screen("Getting ids in cells");

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

     So.message_screen("Identifying closest DM particles for randoms");


#pragma omp parallel for
     for(ULONG i=0; i< N_rand; ++i)//loop over randoms
       {
         int count =0;
         int id_ran=ran_id[i];
         int N_dm_cell=dm_count[id_ran];
         vector<int>i_r_to_dm_dist;
         vector<int>n_index_tot;
         if(N_dm_cell>0)// if there are dm particles in the cell where this random is located, then
           {
             for(int jc = 0; jc < N_dm_cell; ++jc)//loop over the dm in that cell where the random is
               {
                 int jdm=dm_index_cell[id_ran][jc]; //indice entre (0,Ndm-1) que tiene cada particula de dm dentro de la celda id_ran
                 if(jdm!=NO_NUM)
                   {
                     real_prec dist_dm_r=sqrt(pow(x_random_pos[i]-x_dm_pos[jdm],2)+pow(y_random_pos[i]-y_dm_pos[jdm],2)+pow(z_random_pos[i]-z_dm_pos[jdm],2)); //distance between rand and teach central
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

     this->So.message_screen("Collapsing randoms now:");
     this->So.message_screen("Using fraction of distance to closest dm particle = ", this->Distance_fraction);

     this->So.message_screen("Computing new position of randoms");
     int counter=0;

     real_prec density=this->Cosmo.mean_matter_density(this->params._redshift(), (void *)&this->s_cosmo_pars);


     for(ULONG i=0; i<Ntracers; ++i)
       {
         if(this->tracer.Halo[i].identity<0)  //randoms
           {
             int index_dm_closer_a=dm_index_closer_tot[counter]; //total index (i.e, from 0 to N_dms) of the closest dm particle to the ramdom i. Usefil to retrieve coordinates
             real_prec xdm = x_dm_pos[index_dm_closer_a];
             real_prec ydm = y_dm_pos[index_dm_closer_a];
             real_prec zdm = z_dm_pos[index_dm_closer_a];
             real_prec dm_mass = mass_dm[index_dm_closer_a];

             real_prec radius_dm = pow(3.0*dm_mass/(4.*M_PI*density), 1./3.);

             // redefine ran coords to the ref sistem of its closest dm particle:
             real_prec new_x = x_random_pos[counter]-xdm;
             real_prec new_y = y_random_pos[counter]-ydm;
             real_prec new_z = z_random_pos[counter]-zdm;
             real_prec ran_mass = mass_random[counter];

             real_prec radius_ran = pow(3.0*ran_mass/(4.*M_PI*density), 1./3.);

             // get the distance bewteen randoms and their closest DM particle:
             real_prec dist_random_to_dm = sqrt(pow(new_x,2)+pow(new_y,2)+pow(new_z,2));

             // sum of the radii of the two halos
             real_prec RR=radius_dm+radius_ran;
             //          cout<<dm_mass<<"  "<<radius_dm<<"      "<<ran_mass<<"  "<<radius_ran<<"  "<<radius_dm+radius_ran<<"  "<<dist_random_to_dm<<endl;

             //get the angular coordinates:
             real_prec theta=acos(new_z/dist_random_to_dm);
             real_prec phi=atan2(new_y,new_x);

             real_prec new_distance = this->Distance_fraction*dist_random_to_dm;

             // The separation between halos must be greater or equal than the sum of their radii
             real_prec distance_with_exclusion=max(new_distance,RR);

             // transfor to cartesian given a reduced distance and return ot origin of coords:
             this->tracer.Halo[i].coord1=distance_with_exclusion*sin(theta)*cos(phi)+xdm;
             this->tracer.Halo[i].coord2=distance_with_exclusion*sin(theta)*sin(phi)+ydm;
             this->tracer.Halo[i].coord3=distance_with_exclusion*cos(theta)+zdm;
             counter++;
           }
       }

     this->So.DONE();

   }



   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   void Bam::collapse_randoms_isolated()
   {
     cout<<endl;
     this->So.message_screen("**Collapsing randoms towards the DM particles**");
     cout<<endl;


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

     this->So.message_screen("Separating DM and random:");
     ULONG N_dms=0;
     ULONG N_rand=0;

     // struct s_cell_info{
     //   vector<real_prec> id_within_cell;
     // };

     // vector<s_cell_info> dm_inf(box.NGRID);

     int identificator=7;

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
     So.message_screen("Number of tracers associated to DM particles =", N_dms);
     So.message_screen("(", 100.0*static_cast<real_prec>(N_dms)/(static_cast<real_prec>(N_dms)+static_cast<real_prec>(N_rand)), "%)");

     So.message_screen("Number of tracers associated to random particles =", N_rand);
     So.message_screen("(", 100.0*static_cast<real_prec>(N_rand)/(static_cast<real_prec>(N_dms)+static_cast<real_prec>(N_rand)), "%)");

     int max_count=get_max(dm_count); // get the maximum number of dm particles in one cell.
     So.message_screen("Maximum number of DM particles found in one cell (low resolution) = ", max_count);

     // Initialize this to NO_NUM, a negative number, for not all elements of the container
     // dm_index_cell will be filled and initializing with 0 makes confusion, for 0 is a valid/used entry



#define NO_NUM -999
     vector<vector<int> > dm_index_cell(box.NGRID, vector<int>(max_count,NO_NUM));

     dm_count.clear();
     dm_count.resize(box.NGRID,0);
     So.message_screen("Getting ids in cells");

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

     So.message_screen("Identifying closest DM particles for randoms");


#pragma omp parallel for
     for(ULONG i=0; i< N_rand; ++i)//loop over randoms
       {
         int count =0;
         int id_ran=ran_id[i];
         int N_dm_cell=dm_count[id_ran];
         vector<int>i_r_to_dm_dist;
         vector<int>n_index_tot;
         if(N_dm_cell>0)// if there are dm particles in the cell where this random is located, then
           {
             for(int jc = 0; jc < N_dm_cell; ++jc)//loop over the dm in that cell where the random is
               {
                 int jdm=dm_index_cell[id_ran][jc]; //indice entre (0,Ndm-1) que tiene cada particula de dm dentro de la celda id_ran
                 if(jdm!=NO_NUM)
                   {
                     real_prec dist_dm_r=sqrt(pow(x_random_pos[i]-x_dm_pos[jdm],2)+pow(y_random_pos[i]-y_dm_pos[jdm],2)+pow(z_random_pos[i]-z_dm_pos[jdm],2)); //distance between rand and teach central
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

     this->So.message_screen("Collapsing randoms now:");
     this->So.message_screen("Using fraction of distance to closest dm particle = ", this->Distance_fraction);


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
             real_prec dist_random_to_dm=sqrt(pow(new_x,2)+pow(new_y,2)+pow(new_z,2));
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

     this->So.message_screen("Writing new catalog");
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


#ifdef MOCK_MODE
     cout<<endl;
#ifdef _GET_BAM_REALIZATIONS_
     So.message_screen("GENERATING NEW TRACER DENSITY FIELD FROM BIAS AND KERNEL");
#else
     So.message_screen("CALIBRATING BIAS AND KERNEL FROM REFERNECE SIMULATION");
#endif
     cout<<endl;
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

#ifdef _USE_SAT_FRACTION_
     this->tracer.get_density_field_grid(_SAT_FRACTION_, file_Y_sat_frac);
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
     gsl_rng **gBaseRand;
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

#endif



     // If we need to do the calibration,
#ifdef _DO_BAM_CALIBRATION_

     //  we might want to use patchy either to creat the DM, or assign position to partices in mock, or both
#if defined(_USE_PATCHY_) || defined (_GET_BAM_CAT_)

     time_t start_patchy;
     time(&start_patchy);
#endif



#ifdef _USE_PATCHY_


#ifdef OMPPARRAN
     this->patchy.get_dm_field(gBaseRand);
     // Run Patchy!
#else
     this->patchy.get_dm_field(gBaseRand);
     // Run Patchy!
#endif
     this->So.message_screen("Patchy has created DMDF using ALPT");
     So.message_time2(start_patchy);
     cout<<endl;

#endif


#if defined (_USE_PATCHY_) & !defined(_READ_VELOCITIES_)
     // Read the file names generated in Patchy
     file_X =this->patchy.fnameDM+".dat";
     file_Vx=this->patchy.fnameVX+".dat";
     file_Vy=this->patchy.fnameVY+".dat";
     file_Vz=this->patchy.fnameVZ+".dat";
#endif



     this->So.message_screen("BAM:");

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

#if defined(_SEVERAL_REAL_BIAS_) || defined(_SEVERAL_REAL_CAL_)

     this->So.message_screen("Several calibrations from DM fields. Paths to DM fields are hard-coded in line", __LINE__);
     //Used if the Joint B is to be computed for a number of realizations of DM and DMH
     //OR if the calibration from differnet IC is to be performed


     for(int id=this->N_dm_initial ; id<this->N_dm_realizations; ++id)
       {
         //this->Output_directory="../Output_Minerva_R"+to_string(id)+"/";

         if(id<10)
           file_X="/scratch/marcos/DMforBAM/DM_dens_ALPT_000"+to_string(id)+"/densDMALPTrS20.0TETCICz1.000G500V1500.0S5.dat";
         else if(id>=10 && id<100)
           file_X="/scratch/marcos/DMforBAM/DM_dens_ALPT_00"+to_string(id)+"/densDMALPTrS20.0TETCICz1.000G500V1500.0S5.dat";
         else
           file_X="/scratch/marcos/DMforBAM/DM_dens_ALPT_0"+to_string(id)+"/densDMALPTrS20.0TETCICz1.000G500V1500.0S5.dat";

         // If we measure bias to make plots, use the cici VERSION OF THE CATS:
#ifdef BIAS_MODE
         file_Y="/scratch/balaguera/data/Numerics/HADRON-package/classlin/Minerva/TR_DENSITY_MAS1_Nft500_SmoothingScale0_Real"+to_string(id)+"_MCUT0_z1_LambdaTh0_CW0_CWclassMINERVA.dat";
#else
         file_Y="/scratch/balaguera/data/Numerics/HADRON-package/classlin/Minerva/TR_DENSITY_MAS0_Nft500_SmoothingScale0_Real"+to_string(id)+"_MCUT0_z1_LambdaTh0_CW0_CWclassMINERVA.dat";
#endif // end ifdef BIAS_MODE

#ifdef _USE_CWC_
         //this->Output_directory="../Output_Minerva_R"+to_string(id)+"_MCUT0/";
         this->Output_directory="../Output_Minerva_R"+to_string(id)+"_NEWCALIBRATION/";
#else
         //      this->Output_directory="../Output_Minerva_R"+to_string(id)+"_MCUT1_NOCLASS/";
         this->Output_directory="../Output_Minerva_R"+to_string(id)+"_NEWCALIBRATION_NOCLASS/";
#endif
         this->params.Output_directory=this->Output_directory;
#endif // end ifdef _SEVERAL_REAL_BIAS ||


         So.message_warning("Be aweare of this ifdef, I sometimes *comment it* when I need to *read* a reference, line", __LINE__);

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
         this->set_Fourier_vectors();

         // Intialize arrays for power spectrum of input fields
         So.message_warning("Be aware of this ifdef: I sometimnes comment it when I need to *read* a reference, line", __LINE__);
#if defined (_DO_BAM_CALIBRATION_)  || defined (BIAS_MODE)
         this->analyze_input(); // analyze the input references, requested also to create mock for limits
#endif
         // Get numbers from the input density fields and their power spectrum.
         // In this function, if requested from parameter file, an iterative process is performed in order to make the approx method
         // DM field match the reference DM power (if N_iterations_dm >0)


#ifdef _ONLY_PATCHY_
         exit(0);
#endif

         // *******************************************************************************************************
         // *******************************************************************************************************
         // load parameters for the cwclass
#if defined (_USE_CWC_) || defined (_USE_MASS_KNOTS_) || defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_INVARIANT_SHEAR_VFIELD_I_) || defined (_USE_INVARIANT_SHEAR_VFIELD_II_) || defined (_USE_INVARIANT_SHEAR_VFIELD_III_) || defined (_USE_NABLA2DELTA_) || defined (_USE_S2DELTA_) || defined (_USE_S3_)  || defined (_USE_S2_) ||  defined (_USE_DELTA3_)
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
#if defined (_USE_CWC_V_) || defined (_USE_INVARIANT_SHEAR_VFIELD_I_) || defined (_USE_INVARIANT_SHEAR_VFIELD_II_) || defined (_USE_INVARIANT_SHEAR_VFIELD_III_)
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

               time_t start_end;
               time_t start_aux;
               time(&start_end);
               start_aux=start_end;

#ifdef _DO_BAM_CALIBRATION_

               // Loop over the number if iterations demanded to perform the calibration of the Kernel and the Bias


               for(int i=init; i<=this->N_iterations_Kernel ;++i)
                 {
                   time_t start;
                   time(&start);
                   cout<<endl;
                   cout<<CYAN<<"**************************************************************************"<<RESET<<endl;
                   cout<<CYAN<<"**************************************************************************"<<RESET<<endl;
                   this->step=i;
#if defined _USE_CWC_ || defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_INVARIANT_TIDAL_FIELD_III_)
                   this->cwclass.step=i;
#endif
                   // Calculations performed during the iterative process.

                   if(i==0)
                     So.message_screen("BAM raw mapping" , i, start, start_aux);
                   else
                     So.message_screen("Iteration ", i, start,start_aux);
                   cout<<CYAN<<"**************************************************************************"<<RESET<<endl;

                   // ****************************************************************************
                   this->get_BAM_DM();

                   // step i) Do the CWC if requested
                   // Get the ratio T, update kernel K, convolve K with original DM field.
                   // The kernel is computed with the Power aspectrum or the mass weighted power spectrum, according to the preproc def _USE_MASS_WEIGHTED_KERNEL_

                   // ****************************************************************************
                   this->get_BIAS(this->Name_Property_Y);
                   // Step ii)
                   // Compute the halo bias from Numnber counts reference and DM from step i)
                   // ****************************************************************************

                   this->get_mock_grid(_COUNTS_);
                   // Step iii)
                   // Generate the Halo number density field by sampling the DM from step i) using the information of the bias from step ii)
                   // argument false indicates that the DM used in the one of the reference (either original Nbody or approximated)
                   // Also gets the power spectrum of the mock


                   // Since the info of the mass distribution and the sat fraction is not used for the clustering analysis to calibrate the bias,
                   // we can compute them in the last step of the iteration, with the DM already transformed with the Bam kernel.
                   if(i==this->N_iterations_Kernel)
                     {

#ifdef _USE_SAT_FRACTION_
                       this->get_mock_grid(_SAT_FRACTION_);
#endif
                     }

                   // ****************************************************************************
                   start_aux=start;
                 }
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
                   real_prec sigma = sqrt(1./(2.*M_PI*deltak*pow(this->kvec[i],2)*pow(this->Lbox,3)))*(this->Power_NEW[i]); // ESTE PARA EL KERNEL DE LA DM
                   residuals+=pow(this->Power_REF[i]-this->Power_NEW[i],2)/(sigma*sigma);
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

     time(&start_end);
     start_aux=start_end;
     // Sample the bias measured from the reference into other realzations of the Approximated density field.

     // Loop over the new DM fields

     int i=init;
#ifndef _ONLY_POST_PROC_
     for(i=init; i<=N_steps ;++i)
       {
#endif
         time_t start;
         time(&start);
         this->step=init;

         int i_realization=i-this->N_iterations_Kernel+this->N_dm_initial-1;

         this->step=i;
         cout<<CYAN<<"***********************************************************************"<<RESET<<endl;
         cout<<CYAN<<"***********************************************************************"<<RESET<<endl;

         So.message_screen("Creating mock, realization", i_realization, start, start_aux);

         // Get the new dm Df from ALPT and get its properties ({theta}). This will be used for the conditional mass function
         this->get_new_DM_field();

#ifndef _ONLY_POST_PROC_
         //#ifndef _test_mass_assign_
         // Get mock number counts
         this->get_mock_grid(_COUNTS_);
#endif
         //#endif

         //  For every new density filed, we need to compute the mass function N(reference|new DM) such that we can apply it to the new halo sample
         // built from the new dm field

#ifndef _ONLY_POST_PROC_
#ifdef _USE_SAT_FRACTION_
         // Get satellite fraction in cells
         this->get_mock_grid(_SAT_FRACTION_);
#endif
#endif

         // I) Using the information of the bias and kernel, apply these two to a new DM field
         // in order to generate a new halo number denisty field. Measure Power spectrum pof the mock.


         // **************************************************************************** CREAT CATALOG *********************



#ifdef _GET_BAM_CAT_

         //#ifndef _ONLY_POST_PROC_     //comment if the cat is done and only masses are to be assigned

#ifndef _MASS_ASSIGNMENT_TO_REFERENCE_
         this->fnameMOCK=this->Output_directory+"MOCK_TR_realization1_MASY"+to_string(this->iMAS_Y)+"_Nft"+to_string(this->Nft)+"_z"+to_string(this->redshift);
#endif



         //#ifndef _ONLY_POST_PROC_     //comment if the cat is done and only masses are to be assigned
         this->makecat(this->patchy.stradd,this->fnameMOCK,gBaseRand,i_realization);
         //#endif


         //#endif




         //AT THIS STAGE WE NEED TO ASSIGN MASSES TO THE CATALOG JUST CREATED
         //TO DO SO, BEFORE THE LOOP WE HAVE MEASURED THE MASS DISTRIBUTION OF OBJECTS IN CELLS
         // ACCORDING TO THE SAME PROPERTIES OF THE DM FIELD PLUS THE BINS IN HALO MASS, WHICH
         // WILL EFFFECTIVELY ACT AS MASS FUNCTION AS A FUNTION OF DM PROPERTIES
#endif //end get_bam_cat



         start_aux=start;

#ifndef _ONLY_POST_PROC_
       }
#endif


#endif //end of GET_BAM_REALIZATIONS

#endif //end of MOCK_MODE

   }


   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ********************************************************************************************************************************************************************************
   // ***********************************************************************************************************************************// ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
#ifdef MOCK_MODE
   void Bam::dark_matter_to_halos_analytical(const vector<real_prec>& in, vector<real_prec>&out)
   {
     So.message_screen("PERFOMRMING BAM TEST");
     gsl_rng_env_setup();
     gsl_rng_default_seed=75;
     const gsl_rng_type *  T= gsl_rng_ranlux;
     gsl_rng * r = gsl_rng_alloc (T);

#pragma omp parallel for
     for(ULONG i=0 ; i< this->NGRID ; ++i)
       out[i]=0;

     get_overdens(in, out);

     real_prec norm=0;
     for(ULONG i=0 ; i< this->NGRID ; ++i)
       norm+=static_cast<real_prec>(bias_test(out[i], 1.8, 5.0, 0.8));
     norm/=static_cast<real_prec>(this->NGRID);

     for(ULONG i=0 ; i< this->NGRID ; ++i)
       {
         real_prec lambda = (9.0*0.0334821/norm)*static_cast<real_prec>(bias_test(out[i], 1.8, 5.0, 0.8));
         out[i]= gsl_ran_poisson(r,lambda);   //static_cast<float>(lambda)+ gsl_ran_poisson(r,lambda));
       }

     this->params.Name_survey="TR_REF_TEST";
     PowerSpectrumF dPSFo(this->params);
     dPSFo.compute_power_spectrum_grid(out); //ojo, el argumento debe ser DENSIDAD
     dPSFo.write_power_and_modes();

     this->File.write_array("TR_REF_TEST", out);
     this->N_objects_Y=get_nobjects(out);
     So.message_screen("TEST done");

   }
#endif


   // ********************************************************************************************************************************************************************************
   // ***********************************************************************************************************************************// ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************
   // ***************************************************************************************************************************************************************************************

#ifdef MOCK_MODE
void Bam::makecat(string stradd,string fnameMOCK,gsl_rng ** gBaseRand,int ir)
   {
     // In this function we asign coordinates from DM + random particles based on the mock number counts
     // Masses are also assigned and the collapse of randoms towards dm (to correct small scale clusterin) is also performed after that
     // Catalog is then written
     cout<<endl;
     So.message_screen("***********************************************************************");
     So.message_screen("********************GENERATING MOCK CATALOG****************************");
     So.message_screen("***********************************************************************");
     So.message_screen("***********************************************************************");
     cout<<endl;

#ifndef _MASS_ASSIGNMENT_TO_REFERENCE_
     cout<<endl;
     this->So.message_screen("*****Assigning position and velocities to tracers using DM particles****");
     cout<<endl;

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
     // NOte that hwre the numberof dm particles is taht odf the grid

     vector<real_prec> posx(N_dm,0),posy(N_dm,0),posz(N_dm,0);
     vector<ULONG>index(N_dm,0);

     // Filenames of the bninary files containing the positions of the dark matter particles
     this->File.read_array(this->patchy.fnamePOSX+".dat",posx);
     this->File.read_array(this->patchy.fnamePOSY+".dat",posy);
     this->File.read_array(this->patchy.fnamePOSZ+".dat",posz);

     vector<real_prec> velx(N_dm,0),vely(N_dm,0),velz(N_dm,0);
     this->File.read_array(this->patchy.fnameVX+".dat",velx);
     this->File.read_array(this->patchy.fnameVY+".dat",vely);
     this->File.read_array(this->patchy.fnameVZ+".dat",velz);

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

     vector<real_prec> MOCK_DEN_FIELD(this->NGRID,0);
     this->File.read_array(fnameMOCK+".dat",MOCK_DEN_FIELD);


     vector<int> aux_cont(this->NGRID,0);
#pragma omp parallel for
     for(ULONG i=0;i<N_dm;++i)
       {
         real_prec x = static_cast<real_prec>(posx[i]);
         real_prec y = static_cast<real_prec>(posy[i]);
         real_prec z = static_cast<real_prec>(posz[i]);
         ULONG ind=grid_ID(&box, x,y,z);
         index[i]=ind;
#pragma omp atomic update
         aux_cont[ind]++;  //count the number of dm particle in a cell
       }

     vector<bool>dm_used(this->NGRID,false);
     vector<int>dm_cases(this->NGRID,0);

     //ACA DEBEMOS CONTEMPLAR 4 CASOS:

     // I) EN UNA CELDA HAY MAS TRACERS QUE DM
     // II) EN UNA CELDA HAY DM QUE TRACERS
     // III) En una celda hay Tracers pero no hay DM
     // IV) Celdas vacias, deben permanecer vacias


     ULONG Nobjects_mock=get_nobjects(MOCK_DEN_FIELD);
     ULONG Ntracers = Nobjects_mock;  //reduntant, but still
     So.message_screen("Total number of tracers = ", Ntracers);
     cout<<endl;

     ULONG empty_cells_original=0;
#pragma omp parallel for reduction(+:empty_cells_original)
     for(ULONG id=0;id<this->NGRID;++id)
       if(MOCK_DEN_FIELD[id]==0)
         {
           empty_cells_original++;
           dm_cases[id]=4;
           dm_used[id]=false;
         }

     //   So.message_screen("Original number of empty cells ", empty_cells_original);

     So.message_screen("Retrieving positions of dm particles");

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
         ULONG id   =  index[i];    //identify the cell where this particle lives
#ifdef _COLLAPSE_RANDOMS_AUX_
         ULONG id_collapse = grid_ID(&box_collps, x,y,z);
#endif
         if(MOCK_DEN_FIELD[id]>0 &&  aux_cont[id]>0)  // si hay tracers en esta celda y dm tambien
           {
             if(MOCK_DEN_FIELD[id]> aux_cont[id]) //si  hay mas o igual número tracers que dm (y hay dm), tomar todas las dm que hay en cada celda. Al haber mas tracers que dm, vendrá la necesidad de tener randoms
               {
                 real_prec vx=this->patchy.linInterp(x,y,z,velx);
                 real_prec vy=this->patchy.linInterp(x,y,z,vely);
                 real_prec vz=this->patchy.linInterp(x,y,z,velz);
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
                     real_prec vx=this->patchy.linInterp(x,y,z,velx);
                     real_prec vy=this->patchy.linInterp(x,y,z,vely);
                     real_prec vz=this->patchy.linInterp(x,y,z,velz);
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

#pragma omp parallel for
     for(ULONG id=0;id<this->NGRID;++id)
       if(MOCK_DEN_FIELD[id]>0)
         if(1==dm_cases[id] || 2==dm_cases[id])
           Nrandom_tracers[id]=MOCK_DEN_FIELD[id]-cell_inf_dm[id].posx_p.size();


     ULONG Ndm_used_I=0;
     ULONG Nrandoms1=0;
#pragma omp parallel for reduction(+:Ndm_used_I, Nrandoms1)
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
#pragma omp parallel for
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
#pragma omp parallel for reduction(+:Nrandoms3)
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
#pragma omp parallel for reduction(+:ncells_check)
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

         So.message_screen("Generating coordinates for random tracers");
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
                         real_prec rx= static_cast<real_prec>(gsl_rng_uniform(gBaseRand[jthread]));
                         real_prec ry= static_cast<real_prec>(gsl_rng_uniform(gBaseRand[jthread]));
                         real_prec rz= static_cast<real_prec>(gsl_rng_uniform(gBaseRand[jthread]));

                         real_prec x = d1*static_cast<real_prec>(i) + d1*rx;
                         real_prec y = d2*static_cast<real_prec>(j) + d2*ry;
                         real_prec z = d3*static_cast<real_prec>(k) + d3*rz;
                         cell_inf_ran[id].posx_p[ir]=x;
                         cell_inf_ran[id].posy_p[ir]=y;
                         cell_inf_ran[id].posz_p[ir]=z;
                         cell_inf_ran[id].velx_p[ir]=this->patchy.linInterp(x,y,z,velx);
                         cell_inf_ran[id].vely_p[ir]=this->patchy.linInterp(x,y,z,vely);
                         cell_inf_ran[id].velz_p[ir]=this->patchy.linInterp(x,y,z,velz);
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
     So.message_screen("Maximum number of DM particles found in one cell (low resolution) = ", max_count);

     // Initialize this to NO_NUM, a negative number, for not all elements of the container
     // dm_index_cell will be filled and initializing with 0 makes confusion, for 0 is a valid/used entry


#define NO_NUM -999
     vector<vector<int> > dm_index_cell(box_collps.NGRID, vector<int>(max_count,NO_NUM));

     dm_count.clear();
     dm_count.resize(box_collps.NGRID,0);
     So.message_screen("Getting ids in cells");

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

     So.message_screen("Identifying closest DM particles for randoms");


#pragma omp parallel for
     for(ULONG i=0; i< N_rand; ++i)//loop over randoms
       {
         int count =0;
         int id_ran=ran_id[i];
         int N_dm_cell=dm_count[id_ran];
         vector<int>i_r_to_dm_dist;
         vector<int>n_index_tot;
         if(N_dm_cell>0)// if there are dm particles in the cell where this random is located, then
           {
             for(int jc = 0; jc < N_dm_cell; ++jc)//loop over the dm in that cell where the random is
               {
                 int jdm=dm_index_cell[id_ran][jc]; //indice entre (0,Ndm-1) que tiene cada particula de dm dentro de la celda id_ran
                 if(jdm!=NO_NUM)
                   {
                     real_prec dist_dm_r=sqrt(pow(x_random_pos[i]-x_dm_pos[jdm],2)+pow(y_random_pos[i]-y_dm_pos[jdm],2)+pow(z_random_pos[i]-z_dm_pos[jdm],2)); //distance between rand and teach central
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

#endif




     // Before write to file we assign masses. This replaces the call of the function in bamrunner
     // Here we pass the coordinates to a prop vctor and pass it as vector& to the mass assignment function,
     // to then pass it to collapse randoms

     ULONG Nobjects_mock_rc=0;
     this->tracer.NOBJS=Nobjects_mock;
     this->tracer.Halo.resize(this->tracer.NOBJS);

     this->So.message_screen("Filling structure with dm positions:");

     //#pragma omp parallel for reduction (+:Nobjects_mock_rc)
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
                 this->tracer.Halo[Nobjects_mock_rc].identity = 1.0;
                 this->tracer.Halo[Nobjects_mock_rc].GridID = cell_inf_dm[id].id_p[in];
                 Nobjects_mock_rc++;
               }
       }
     So.DONE();

     So.message_screen("Number of tracers associated to DM particles =", Nobjects_mock_rc);
     So.message_screen("(", 100.0*static_cast<real_prec>(Nobjects_mock_rc)/(static_cast<real_prec>(Nobjects_mock_rc)+static_cast<real_prec>(Nrandoms1)), "%)");

     So.message_screen("Freeing memory");
     cell_inf_dm.clear(); cell_inf_dm.shrink_to_fit();
     So.DONE();


     ULONG N_dms_b= Nobjects_mock_rc;

     if(Nrandoms1>0)
       {
         this->So.message_screen("Adding random particles:");
         for(ULONG id=0;id<this->NGRID;++id)
           {
             if(MOCK_DEN_FIELD[id]>0)
               if(true==random_used[id])
                 for(int in=0;in<cell_inf_ran[id].posx_p.size();++in)
                   {
                     this->tracer.Halo[Nobjects_mock_rc].coord1=cell_inf_ran[id].posx_p[in];
                     this->tracer.Halo[Nobjects_mock_rc].coord2=cell_inf_ran[id].posy_p[in];
                     this->tracer.Halo[Nobjects_mock_rc].coord3=cell_inf_ran[id].posz_p[in];
                     this->tracer.Halo[Nobjects_mock_rc].vel1=cell_inf_ran[id].velx_p[in];
                     this->tracer.Halo[Nobjects_mock_rc].vel2=cell_inf_ran[id].vely_p[in];
                     this->tracer.Halo[Nobjects_mock_rc].vel3=cell_inf_ran[id].velz_p[in];
                     this->tracer.Halo[Nobjects_mock_rc].identity=-1.0;
                     this->tracer.Halo[Nobjects_mock_rc].GridID=cell_inf_ran[id].id_p[in];
                     Nobjects_mock_rc++;
                   }
           }
         So.DONE();

         So.message_screen("Number of tracers associated to random particles =", Nrandoms1);
         So.message_screen("(", 100.0*static_cast<real_prec>(Nrandoms1)/(static_cast<real_prec>(N_dms_b)+static_cast<real_prec>(Nrandoms1)), "%)");

         So.message_screen("Freeing memory");
         cell_inf_ran.clear(); cell_inf_ran.shrink_to_fit();
         random_used.clear();
         random_used.shrink_to_fit();
         So.DONE();
       }



#ifdef _COLLAPSE_RANDOMS_AUX_

     // THE GrdID kept in the sturcture Tracer is the *ORIGINAL*, not the one computed after the collapse of randoms

     cout<<endl;
     this->So.message_screen("Collapsing randoms.");
     this->So.message_screen("Using fraction of distance to closest dm particle = ", this->Distance_fraction);
     this->So.message_screen("Computing new position of randoms:");
     int counter=0;
     for(ULONG i=0; i<Ntracers; ++i)
       {
         if(this->tracer.Halo[i].identity<0)  //randoms
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
             real_prec dist_random_to_dm=sqrt(pow(new_x,2)+pow(new_y,2)+pow(new_z,2));
             //get the angular coordinates:
             real_prec theta=acos(new_z/dist_random_to_dm);
             real_prec phi=atan2(new_y,new_x);
             real_prec new_distance = this->Distance_fraction*dist_random_to_dm;
             // transfor to cartesian given a reduced distance and return ot origin of coords:
             this->tracer.Halo[i].coord1=new_distance*sin(theta)*cos(phi)+xdm;
             this->tracer.Halo[i].coord2=new_distance*sin(theta)*sin(phi)+ydm;
             this->tracer.Halo[i].coord3=new_distance*cos(theta)+zdm;
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

#endif




     // {
     //   cout<<endl;

     //   this->So.message_screen("****POWER SPECTRUM***PROBLEMS WITH THE CIC, DO NOT AGREE WITH THE CIC I DO FRM OUTPUT CATALOG. TRUST THE ONE BELOW, NOT THIS ONE HERE*");

     //   this->tracer.box.masskernel=1; // using CIC
     //   this->delta_Y_new.resize(this->NGRID,0);
     //   So.message_screen("Interpolating tracer catalog on a grid using MAS", this->tracer.box.masskernel);
     //   MAS_NEW(&(this->tracer.box),this->tracer.Halo, _COUNTS_, this->delta_Y_new);
     //   /// Get mock power spectrum from mock catalog
     //   this->get_power_spectrum("TR_MOCK_CATb");
     //   this->delta_Y_new.clear();
     //   this->delta_Y_new.shrink_to_fit();
     // }


#endif




     //  NOw we can call the assign_mass function


#ifdef _MASS_ASSIGNMENT_TO_REFERENCE_
         string outputFileName=this->Output_directory+"CAT_realization"+to_string(ir)+"_"+this->patchy.stradd+string(".txt_pos_from_ref");
#else
         string outputFileName=this->Output_directory+"CAT_realization"+to_string(ir)+"_"+this->patchy.stradd+string(".txt");
#endif



#ifdef _ASSIGN_PROPERTY_

     // Get info from the reference
     this->get_X_function();

     this->tracer.type_of_object="TRACER_MOCK_ONLY_COORDS"; // THis name is important.


#ifdef _MASS_ASSIGNMENT_TO_REFERENCE_
     cout<<endl;

     So.message_warning("***NOTE: The reference is currently being used as mock in order to re-assign masses anc do tests, see line", __LINE__);
     this->params.file_catalogue=this->Output_directory+"newcat.txt";
     //    this->params.file_catalogue=this->params._file_catalogue();
     this->params.i_coord1_g=0;
     this->params.i_coord2_g=1;
     this->params.i_coord3_g=2;
     this->params.i_vmax_g=-1;
     this->params.i_mass_g=-1;
     this->tracer.set_params_catalog(this->params);



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


#ifdef _USE_MASS_ASSIGNMENT_READING_REF_MASSES_
     this->assign_tracer_property_new_new(true);
#else
     this->assign_tracer_mass_new();
#endif
     this->So.DONE();


     this->tracer.type_of_object="TRACER_MOCK"; // Once masses are assigned, change the identity of this->tracer. This name is important.


   // Collapse randoms:  this applies if collapse_randoms_aux is undef
#ifdef _COLLAPSE_RANDOMS_
         this->collapse_randoms();
#endif


         this->patchy.fnameTRACERCAT=outputFileName;

         /*
           So.message_screen("================================================================");
           So.message_screen("Doing test writting same input masses M->pow(10,m*(1+N(0,sigma))");
           So.message_screen("================================================================");
           // We perturn the reference masses by a log-normal distributed perturbation to assess the impact on the mass-weighted power spectrum
           // So far with sigma = 0.2 effets are sizable on small scales
           const gsl_rng_type *  T;
           gsl_rng * r ;
           gsl_rng_env_setup();
           gsl_rng_default_seed=52511;
           T = gsl_rng_ranlux;
           r = gsl_rng_alloc (T);
           for(ULONG i = 0; i< this->tracer.NOBJS; ++i)
           this->tracer.Halo[i].mass=this->tracer_ref.Halo[i].mass*(1+gsl_ran_gaussian(r,0.1));
         */




        Params params_aux=params;
         params_aux.file_catalogue=outputFileName;
         params_aux.weight_with_mass=false;
#ifdef _MASS_WEIGHT_POWER_
      params_aux.weight_with_mass=true;
      params_aux.file_power="MOCK_mass_weight";
      params_aux.SN_correction=false;
#else
      params_aux.file_power="MOCK";
      params_aux.SN_correction=true;
#endif
      params_aux.Name_survey="BAM_CAT_from_ref";

      PowerSpectrumF Power(params_aux);
      Power.compute_power_spectrum(false,this->tracer.Halo);

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
        count_p=0;
         for(int i=0;i<this->Nft/2/this->ndel_data;++i)
         if(this->Power_REF_MW[i]>0)
          {
            count_p++;
            residuals+= fabs(this->Power_NEW_MW[i]/this->Power_REF_MW[i]-1.0);
          }
         residuals/=static_cast<real_prec>(count_p)/100.0;
         So.message_screen("At Niquyst frequency" );
         So.message_screen("Residuals from power spectrum  = ", residuals, "%");


         cout<<endl;


#ifdef _ASSIGN_MASS_POST_
// Here we can now complement the property assignment, using the inormation of vmax already assigned

         So.message_screen("********************************************************************************************");
         So.message_screen("Assigning halo masses using VMAX information::");
         So.message_screen("**************************************");
         cout<<endl;
         cout<<endl;
         this->params.i_vmax_g=8;
         this->tracer.type_of_object="TRACER_MOCK"; // THis name is important.
         this->get_X_function_complement();
         this->assign_tracer_property_new_new(false);
#endif

#endif // end if ASSIGN_PROPERTY


#ifdef _GET_DIST_MIN_SEP_MOCK_
    this->tracer_ref.get_distribution_min_separations(this->tracer.type_of_object, this->ncells_info);
#endif



         real_prec conversion_factor=1.0;
#ifdef _VEL_UNITS_MPC_PER_h_
         conversion_factor=cgs_km/(this->s_cosmo_info.Hubble_parameter); //this transforms to km/s with cgs_km and divide by 100h km/s /(Mpc/h) to leave units in Mpc/h
#endif



#define _write_
#ifdef _write_
             this->So.message_screen("Writting to file", outputFileName);
             ofstream outStream;
             outStream.open(outputFileName.c_str());
             assert(outStream.is_open());
             outStream.precision(_PREC_OUTPUT_);
             outStream.setf(ios::showpoint);
             outStream.setf(ios::scientific);
             for(ULONG i = 0; i< this->tracer.NOBJS; ++i)
#ifdef _ASSIGN_PROPERTY_
#ifdef _USE_MASS_AS_OBSERVABLE_
               outStream<<this->tracer.Halo[i].coord1<<" "<<this->tracer.Halo[i].coord2<<" "<<this->tracer.Halo[i].coord3<<" "<<this->tracer.Halo[i].vel1*conversion_factor<<" "<<this->tracer.Halo[i].vel2*conversion_factor<<" "<<this->tracer.Halo[i].vel3*conversion_factor<<" "<<pow(this->tracer.Halo[i].mass, 1./exponent_mass_tracer)<<endl;
#elif defined _USE_VMAX_AS_OBSERVABLE_
#ifdef _ASSIGN_MASS_POST_
                 outStream<<this->tracer.Halo[i].coord1<<" "<<this->tracer.Halo[i].coord2<<" "<<this->tracer.Halo[i].coord3<<" "<<this->tracer.Halo[i].vel1*conversion_factor<<" "<<this->tracer.Halo[i].vel2*conversion_factor<<" "<<this->tracer.Halo[i].vel3*conversion_factor<<" "<<this->tracer.Halo[i].mass<<" "<<this->tracer.Halo[i].vmax<<endl;
#endif
#else
                 outStream<<this->tracer.Halo[i].coord1<<" "<<this->tracer.Halo[i].coord2<<" "<<this->tracer.Halo[i].coord3<<" "<<this->tracer.Halo[i].vel1*conversion_factor<<" "<<this->tracer.Halo[i].vel2*conversion_factor<<" "<<this->tracer.Halo[i].vel3*conversion_factor<<"  "<<this->tracer.Halo[i].vmax<<endl;
#endif
#else
                 outStream<<this->tracer.Halo[i].coord1<<" "<<this->tracer.Halo[i].coord2<<" "<<this->tracer.Halo[i].coord3<<" "<<this->tracer.Halo[i].vel1*conversion_factor<<" "<<this->tracer.Halo[i].vel2*conversion_factor<<" "<<this->tracer.Halo[i].vel3*conversion_factor<<endl;
#endif

             outStream.close();
             So.DONE();
#endif


             So.message_screen("Freeing memmory from tracer");
             this->tracer.Halo.clear();
             this->tracer.Halo.shrink_to_fit();
             So.DONE();





}
#endif






// ********************************************************************************************************************************************************************************
// ********************************************************************************************************************************************************************************
// ********************************************************************************************************************************************************************************
// ********************************************************************************************************************************************************************************
// ***********************************************************************************************************************************// ***************************************************************************************************************************************************************************************
// ***************************************************************************************************************************************************************************************
// ***************************************************************************************************************************************************************************************
#ifdef MOCK_MODE

void Bam::correct_for_exclusion(ULONG LENGHTdm){

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
        int N_masses_in_bin = number_in_theta_mock[ih]; //Available Masses in bin of THETA, regardless of the separation between objects

        for(int j=0; j< this->tracer.masses_in_cells_min_sep[ih].mass_to_swap.size();++j)//loopover the masses to be swaped
          {

           ULONG old_ig=this->tracer.masses_in_cells_min_sep[ih].mass_to_swap[j]; // galaxy index IDg [0,NOBJS) of the mass requested to be swaped
           real_prec old_mass=this->tracer.Halo[old_ig].mass; //mass requested to be swaped

           bool baux=false;
           while(baux==false)
            {
              int i_mass_halo_label= gsl_rng_uniform_int(rn,N_masses_in_bin); // pick-up randomly some other mass in the same bin-theta
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
