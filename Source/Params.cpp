#include "../Headers/Params.h"
#include "../Headers/cosmo_parameters.h"

void Params::init_pars()
{

  this->NX = 200;
  this->NY = 100;
  this->n_sknot_massbin = 200;
  this->n_vknot_massbin = 200;
  this->NMASSbins = 1;
  this->NMASSbins_power = 1;
  this->iMAS_X = 0;
  this->iMAS_X_REF_PDF = 0;
  this->iMAS_X_NEW = 0;
  this->iMAS_Y = 0;
  this->lambdath = 0.0;
  this->lambdath_v = 0.0;
  this->realization = 1;
  this->redshift = 0.0;
  this->Initial_Redshift_DELTA=0;
  this->smscale = 0.0;
  this->delta_Y_min = 0.;
  this->delta_Y_max = 10;
  this->delta_X_min = 0;
  this->delta_X_max = 10.;
  this->ldelta_Y_min = 0;
  this->ldelta_Y_max = 10;
  this->ldelta_X_min = 0.;
  this->ldelta_X_max = 10.;
  this->Input_Directory_Y = "../Input";
  this->Name_Catalog_Y = "UNKNOWN_NAME";
  this->Name_Catalog_Y_HR = "UNKNOWN_NAME";
  this->Name_Catalog_Y_MWEIGHTED = "UNKNOWN_NAME";
  this->Input_Directory_X = "../Input";
  this->Input_Directory_X_REF = "../Input";
  this->Input_Directory_X_NEW = "../Input";
  this->Input_Directory_BIAS_KERNEL = "../Input";
  this->Input_Directory_BIAS_KERNEL_TWO = "../Input";
  this->XNAME = "XNAME";
  this->YNAME = "YNAME";
  this->Name_Catalog_X = "UNKNOWN_NAME";
  this->Name_VelFieldx_X = "UNKNOWN_NAME";
  this->Name_VelFieldy_X = "UNKNOWN_NAME";
  this->Name_VelFieldz_X = "UNKNOWN_NAME";
  this->Name_Catalog_X_REF_PDF = "UNKNOWN_NAME";
  this->Name_Catalog_X_NEW = "UNKNOWN_NAME";
  this->extra_info = "";
  this->Name_Property_X = "X_PREPERTY";
  this->Name_Property_Y = "Y_PROPERTY";
  this->Output_directory = "../output";
  this->Quantity= "QUANTITY";
  this->Scale_X= "linear";
  this->Scale_Y= "linear";
  this->N_iterations_Kernel = 1;
  this->N_dm_realizations = 0;
  this->N_dm_initial = 1;
  this->N_iterations_dm = 1;
  this->Comp_conditional_PDF = false;
  this->Apply_Rankordering = false;
  this->Apply_Rankordering_ab_initio = false;
  this->vel_units_g = "kmps";
  this->IC_index=0;
  this->Redefine_limits= false;
  this->Write_Scatter_Plot= false;
  this->Write_PDF_number_counts= false;
  this->Convert_Density_to_Delta_X= false;
  this->Convert_Density_to_Delta_Y= false;
  this->Comp_joint_PDF = false;
  this->write_files_for_histograms = false;
  this->n_cwt=ONE;
  this->n_cwv=ONE;
  this->M_exclusion=1e13;
  this->statistics = "fkp";
  this->Input_dir_cat = "../Input/";
  this->Input_dir_cat_TWO = "../Input/";
  this->dir_output = "../Output/";
  this->file_catalogue = "cat.dat";
  this->input_type = "catalog";
  this->input_type_two = "catalog";
  this->delta_grid_file = "delta_file";
  this->delta_grid_file2 = "delta_file";
  this->delta_grid_file3 = "delta_file";
  this->delta_grid_file4 = "delta_file";
  this->ngal_delta = 1;
  this->measure_cross = true;
  this->measure_cross_from_1 = 1;
  this->measure_cross_from_2 = 2;
  this->sys_of_coord_g = 1;
  this->redshift_space_coords_g =false;
  this->i_coord1_g = 0;
  this->i_coord2_g = 1;
  this->i_coord3_g = 2;
  this->i_mean_density_g = 3;
  this->i_weight1_g = 4;
  this->i_weight2_g = 5;
  this->i_weight3_g = 6;
  this->i_weight4_g = 7;
  this->i_rs_g = 6;
  this->i_virial_g = 7;
  this->i_coord1_dm = 0;
  this->i_coord2_dm = 1;
  this->i_coord3_dm = 2;
  this->i_mass_dm = 7;
  this->type_of_object="TRACER";

  this->use_weight1_g = true;
  this->use_weight2_g = true;
  this->use_weight3_g = true;
  this->use_weight4_g = true;
  this->i_property1_g = 0;
  this->i_property2_g = 1;
  this->i_property3_g = 2;
  this->i_property4_g = 3;
  this->i_property5_g = 4;
  this->i_property6_g = 5;
  this->new_Lbox = true;
  this->n_catalogues = 1;
  this->angles_units_g = "D";
  this->use_random_catalog = false;
  this->use_random_catalog_cl = false;
  this->sys_of_coord_r = 1;
  this->i_coord1_r = 0;
  this->i_coord2_r = 1;
  this->i_coord3_r = 2;
  this->i_mean_density_r = 3;
  this->i_weight1_r = 4;
  this->i_weight2_r = 5;
  this->i_weight3_r = 6;
  this->i_weight4_r = 7;

  this->i_property1_r = 0;
  this->i_property2_r = 1;
  this->i_property3_r = 2;
  this->i_property4_r = 3;
  this->i_property5_r = 4;
  this->i_property6_r = 5;

  this->i_mask_pixel= 0;
  this->i_mask_alpha= 1;
  this->i_mask_delta= 2;
  this->i_mask_flag= 3;

  this->use_weight1_r = true;
  this->use_weight2_r = true;
  this->use_weight3_r = true;
  this->use_weight4_r = true;
  this->angles_units_r = "D";
  this->file_random = "random";
  this->Name_survey = "survey";


  // parameters for the power spectrum
  this->new_los = false;
  this->Nft = 64;
  this->Nft_low = 64;
  this->Nft_random_collapse = 32;
  this->Distance_fraction=1.0;
  this->Lbox = 100;
  this->Lbox_low = 100;
  this->mass_assignment_scheme = "NGP";
  this->MAS_correction = false;
  this->type_of_binning = "linear";
  this->DeltaKmin = 0;
  this->N_log_bins = 10;
  this->ndel_data = 1;
  this->ndel_window = 1;
  this->N_mu_bins = 10;
  this->FKP_weight = false;
  this->Pest = 1000.0;
  this->SN_correction = false;
  this->FKP_error_bars = false;
  this->FKP_error_bars_exact = false;
  this->nbar_tabulated = false;
  this->constant_depth = false;
  this->N_z_bins = 3;
  this->redshift_min_sample = 0;
  this->redshift_max_sample = 1;
  N_dndz_bins = 2;
  new_N_dndz_bins = 2;
  area_survey = 1000;
  Healpix_resolution = 4;
  file_power = "power";
  file_power_log = "power_log";

  file_window = "window";
  file_dndz = "dndz";
  file_power2d = "power2d";
  file_power2d_mk = "power2d_mk";
  file_bispectrum = "bispectrum";

  //output files for Yamamoto estimator of moments
  kmax_y_ds = 0.2;

  //Parameters for the bispectrum
  kmax_bk = 0.2;
  kmin_bk = 0.1;
  use_fundamental_mode_as_kmin_bk = false;

  this->iteration_ini = 0;

  this->ic_alias_corrected = false;
  this->ic_input_type = DENSITY;

  this->masskernel=0;
  this->masskernel_vel=0;

  // cosmological parameters
  this->om_matter = COSMOPARS::Om_matter;
  this->om_radiation = COSMOPARS::Om_radiation;
  this->om_baryons = COSMOPARS::Om_baryons;
  this->om_cdm = this->om_matter-this->om_baryons;
  this->om_vac = COSMOPARS::Om_vac;
  this->om_k   = COSMOPARS::Om_k;
  this->Hubble = COSMOPARS::Hubble;
  this->hubble = COSMOPARS::hubble;
  this->spectral_index = COSMOPARS::n_s;
  this->w_eos = COSMOPARS::w_eos;
  this->N_eff =  COSMOPARS::N_eff;
  this->sigma8 = COSMOPARS::sigma8;
  this->Tcmb = COSMOPARS::Tcmb;
  this->use_wiggles = true;
  this->RR = COSMOPARS::RR;
  this->get_distribution_min_separations=false;

}


// *************************************************************************************************************
// *************************************************************************************************************
// *************************************************************************************************************
// *************************************************************************************************************
// *************************************************************************************************************
// *************************************************************************************************************
// *************************************************************************************************************
// *************************************************************************************************************
// *************************************************************************************************************
// *************************************************************************************************************
// *************************************************************************************************************
// *************************************************************************************************************



void Params::read_pars(string file)
{

  ifstream fin_parameters (file.c_str());


  if (fin_parameters.fail()) { 
    cerr <<"Error in opening the parameters file "<<file<<std::endl;
    cerr <<"Cosmicatlass stops here."<<std::endl; 
    cout<<std::endl;
    exit(1); 
  }

  string line_in_file;
  string par_name;
  string equality;
  string par_value;
  while (getline(fin_parameters,line_in_file))
    {
      // ignore lines starting with hashtag
      if (line_in_file[0] != '#' && line_in_file.empty()==0)
	{

	  // read parameter name
	  stringstream line_string (line_in_file);
	  line_string << line_in_file;
	  line_string >> par_name;


	  // check that second word is "="
	  line_string >> equality;

	  if (equality != "=")
	    {
	      std::cerr << "Error in parameters file at " << par_name << " = " << "???" << std::endl;
	      std::cerr << "Please check the parameter file "<<file<< std::endl; 
        std::cerr << "Cosmicatlass stops here."<<file<< std::endl; 
	      exit(1);
	    }

	  // read parameter value
	  line_string >> par_value;


	  if (par_value.empty())
	    {
	      cerr << "Value of " << par_name << " not specified in " <<file << endl;
	      cerr << "Assuming a default value for " << par_name << endl;
	      continue;
	    }

	  // Parameters for BAM
	  if (par_name == "NX")this->NX = atoi(par_value.c_str()); //this->NX.read_from_params=true;
	  else if (par_name == "NY")this->NY = atoi(par_value.c_str());
	  else if (par_name == "NY_MASS")this->NY_MASS = atoi(par_value.c_str());
	  else if (par_name == "NY_SAT_FRAC")this->NY_SAT_FRAC = atoi(par_value.c_str());
	  else if (par_name == "N_SKNOT_MASSBIN")this->n_sknot_massbin = atoi(par_value.c_str());
	  else if (par_name == "N_VKNOT_MASSBIN")this->n_vknot_massbin = atoi(par_value.c_str());
	  else if (par_name == "iMAS_X")this->iMAS_X = atoi(par_value.c_str());
	  else if (par_name == "iMAS_X_REF_PDF")this->iMAS_X_REF_PDF = atoi(par_value.c_str());
	  else if (par_name == "iMAS_X_NEW")this->iMAS_X_NEW = atoi(par_value.c_str());
	  else if (par_name == "iMAS_Y")this->iMAS_Y = atoi(par_value.c_str());
    else if (par_name == "Unitsim_plabel")this->Unitsim_plabel = atoi(par_value.c_str());
	  else if (par_name == "Number_of_chunks_new_dm")this->Number_of_chunks_new_dm = atoi(par_value.c_str());
	  else if (par_name == "lambdath")this->lambdath = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "lambdath_v")this->lambdath_v = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "Realization")this->realization = atoi(par_value.c_str());
	  else if (par_name == "Redshift")this->redshift = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "Initial_Redshift_DELTA")this->Initial_Redshift_DELTA = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "smscale")this->smscale = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "delta_Y_min")this->delta_Y_min = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "delta_Y_max")this->delta_Y_max = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "delta_X_min")this->delta_X_min = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "delta_X_max")this->delta_X_max = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "ldelta_Y_min")this->ldelta_Y_min = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "ldelta_Y_max")this->ldelta_Y_max = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "ldelta_X_min")this->ldelta_X_min = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "ldelta_X_max")this->ldelta_X_max =  static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "Input_Directory_Y")this->Input_Directory_Y = par_value;
	  else if (par_name == "Name_Catalog_Y")this->Name_Catalog_Y = par_value;
	  else if (par_name == "Name_Catalog_Y_HR")this->Name_Catalog_Y_HR = par_value;
	  else if (par_name == "Name_Catalog_Y_MWEIGHTED")this->Name_Catalog_Y_MWEIGHTED = par_value;
	  else if (par_name == "Input_Directory_X")this->Input_Directory_X = par_value;
	  else if (par_name == "Input_Directory_BIAS_KERNEL")this->Input_Directory_BIAS_KERNEL = par_value;
	  else if (par_name == "Input_Directory_BIAS_KERNEL_TWO")this->Input_Directory_BIAS_KERNEL_TWO = par_value;
	  else if (par_name == "Input_Directory_X_REF")this->Input_Directory_X_REF = par_value;
    else if (par_name == "Input_Directory_X_REF_TWO")this->Input_Directory_X_REF_TWO = par_value;
	  else if (par_name == "Input_Directory_X_NEW")this->Input_Directory_X_NEW = par_value;
	  else if (par_name == "XNAME")this->XNAME = par_value;
	  else if (par_name == "YNAME")this->YNAME = par_value;
	  else if (par_name == "Name_Catalog_X")this->Name_Catalog_X = par_value;
	  else if (par_name == "Name_VelFieldx_X")this->Name_VelFieldx_X = par_value;
	  else if (par_name == "Name_VelFieldy_X")this->Name_VelFieldy_X = par_value;
	  else if (par_name == "Name_VelFieldz_X")this->Name_VelFieldz_X = par_value;
	  else if (par_name == "Name_redshift_mask")this->Name_redshift_mask = par_value;
	  else if (par_name == "Name_binary_mask")this->Name_binary_mask = par_value;
	  else if (par_name == "IC_index")this->IC_index = atoi(par_value.c_str());
	  else if (par_name == "Name_Catalog_X_REF_PDF")this->Name_Catalog_X_REF_PDF = par_value;
	  else if (par_name == "Name_Catalog_X_NEW")this->Name_Catalog_X_NEW = par_value;
	  else if (par_name == "Name_Catalog_X_NEW_TWO")this->Name_Catalog_X_NEW_TWO = par_value;
	  else if (par_name == "Name_Property_X")this->Name_Property_X = par_value;
	  else if (par_name == "Name_Property_Y")this->Name_Property_Y = par_value;
	  else if (par_name == "Output_directory")this->Output_directory = par_value;
	  else if (par_name == "Quantity")this->Quantity= par_value;
	  else if (par_name == "Scale_X")this->Scale_X= par_value;
	  else if (par_name == "Scale_Y")this->Scale_Y= par_value;
	  else if (par_name == "N_iterations_Kernel")this->N_iterations_Kernel = atoi(par_value.c_str());
	  else if (par_name == "N_dm_realizations")this->N_dm_realizations = atoi(par_value.c_str());
	  else if (par_name == "N_dm_initial")this->N_dm_initial = atoi(par_value.c_str());
	  else if (par_name == "N_iterations_dm")this->N_iterations_dm = atoi(par_value.c_str());
	  else if (par_name == "iteration_ini")this->iteration_ini = atoi(par_value.c_str());
	  else if (par_name == "EXTRA_INFO")this->extra_info = par_value.c_str();
	  else if (par_name == "MASS_units") this->MASS_units = atof(par_value.c_str());
	  else if (par_name == "LOGMASSmin") this->LOGMASSmin = atof(par_value.c_str());
	  else if (par_name == "VMAXmax") this->VMAXmax = atof(par_value.c_str());
	  else if (par_name == "VMAXmin") this->VMAXmin = atof(par_value.c_str());
	  else if (par_name == "RSmax") this->RSmax = atof(par_value.c_str());
	  else if (par_name == "RSmin") this->RSmin = atof(par_value.c_str());
	  else if (par_name == "SPINmax") this->SPINmax = atof(par_value.c_str());
	  else if (par_name == "SPINmin") this->SPINmin = atof(par_value.c_str());
	  else if (par_name == "LOGMASSmax") this->LOGMASSmax = atof(par_value.c_str());
	  else if (par_name == "NMASSbins")this->NMASSbins = atoi(par_value.c_str());
	  else if (par_name == "NMASSbins_power")this->NMASSbins_power = atoi(par_value.c_str());
	  else if (par_name == "NMASSbins_mf") this->NMASSbins_mf = atoi(par_value.c_str());
	  else if (par_name == "file_bin_x_coord" ) this->file_bin_x_coord = par_value.c_str();
	  else if (par_name == "file_bin_y_coord" ) this->file_bin_y_coord = par_value.c_str();
	  else if (par_name == "file_bin_z_coord" ) this->file_bin_z_coord = par_value.c_str();
	  else if (par_name == "N_lines_binary") this->N_lines_binary = atoi(par_value.c_str());
	  else if (par_name == "Input_dir_cat") this->Input_dir_cat = par_value;
    else if (par_name == "Input_dir_cat_TWO") this->Input_dir_cat_TWO = par_value;
	  else if (par_name == "Type_of_object") this->type_of_object = par_value;
	  else if (par_name == "Prop_threshold_multi_scale_1") this->Prop_threshold_multi_scale_1 = atof(par_value.c_str());
	  else if (par_name == "Prop_threshold_multi_scale_2") this->Prop_threshold_multi_scale_2 = atof(par_value.c_str());
	  else if (par_name == "Prop_threshold_multi_scale_3") this->Prop_threshold_multi_scale_3 = atof(par_value.c_str());
	  else if (par_name == "Prop_threshold_multi_scale_4") this->Prop_threshold_multi_scale_4 = atof(par_value.c_str());
	  else if (par_name == "Tolerance_factor_l1") this->Tolerance_factor_l1 = atof(par_value.c_str());
	  else if (par_name == "Tolerance_factor_l2") this->Tolerance_factor_l2 = atof(par_value.c_str());
	  else if (par_name == "Tolerance_factor_l3") this->Tolerance_factor_l3 = atof(par_value.c_str());
	  else if (par_name == "Tolerance_factor_l4") this->Tolerance_factor_l4 = atof(par_value.c_str());
	  else if (par_name == "Nft_low_l1") this->Nft_low_l1 = static_cast<ULONG>(atoi(par_value.c_str()));
	  else if (par_name == "Nft_low_l2") this->Nft_low_l2 = static_cast<ULONG>(atoi(par_value.c_str()));
	  else if (par_name == "Nft_low_l3") this->Nft_low_l3 = static_cast<ULONG>(atoi(par_value.c_str()));
	  else if (par_name == "Nft_low_l4") this->Nft_low_l4 = static_cast<ULONG>(atoi(par_value.c_str()));
	  else if (par_name == "Comp_conditional_PDF")
	    {
	      if (par_value=="true")this->Comp_conditional_PDF = true;
	      else this->Comp_conditional_PDF = false;
	    }
	  else if (par_name == "Normalize_initial_redshift")
	    {
	      if (par_value=="true")this->Normalize_initial_redshift = true;
	      else this->Normalize_initial_redshift = false;
	    }
	  else if (par_name == "Apply_Rankordering")
	    {
	      if (par_value=="true")this->Apply_Rankordering = true;
	      else this->Apply_Rankordering = false;
	    }
	  else if (par_name == "Apply_Rankordering_ab_initio")
	    {
	      if (par_value=="true")this->Apply_Rankordering_ab_initio = true;
	      else this->Apply_Rankordering_ab_initio = false;
	    }
	  else if (par_name == "Redefine_limits")
	    {
	      if (par_value=="true")this->Redefine_limits= true;
	      else this->Redefine_limits = false;
	    }
	  else if (par_name == "Write_Scatter_Plot")
	    {
	      if (par_value=="true")this->Write_Scatter_Plot= true;
	      else this->Write_Scatter_Plot = false;
	    }
	  else if (par_name == "Write_PDF_number_counts")
	    {
	      if (par_value=="true")this->Write_PDF_number_counts= true;
	      else this->Write_PDF_number_counts = false;
	    }
	  else if (par_name == "Convert_Density_to_Delta_X")
	    {
	      if (par_value=="true")this->Convert_Density_to_Delta_X= true;
	      else this->Convert_Density_to_Delta_X = false;
	    }
	  else if (par_name == "Convert_Density_to_Delta_Y")
	    {
	      if (par_value=="true")this->Convert_Density_to_Delta_Y= true;
	      else this->Convert_Density_to_Delta_Y = false;
	    }
	  else if (par_name == "Comp_joint_PDF")
	    {
	      if (par_value=="true")this->Comp_joint_PDF = true;
	      else this->Comp_joint_PDF = false;
	    }
	  else if (par_name == "write_files_for_histograms")
	    {
	      if (par_value=="true")this->write_files_for_histograms = true;
	      else this->write_files_for_histograms = false;
	    }

#ifdef _USE_CWC_
	  else  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 't' && line_in_file[1] == '_'))
	    {
	      stringstream line_string (line_in_file);
	      line_string << line_in_file;
	      line_string >> par_name;   // read the first character, the name of the parameter
	      line_string >> equality; 	  // check that second character is "="
	      if (equality != "=")
      		{
		  cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
		  cerr << "Using a default value for " << par_name << endl; exit(1);
      		}
	      string par_value0;
	      int ending=0;
	      int n_m_lims=0;
	      do
		{
		  line_string >> par_value0;  //read value
		  n_m_lims++;
		  if(par_value0=="END")
		    ending=1;
		  this->cwt_used.push_back(atof(par_value0.c_str()));
		}while(ending!=1);
	      this->cwt_used.pop_back(); //delete the last element, which is END
	    }
#endif
#ifdef _USE_CWC_V_
          else  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'v' && line_in_file[1] == '_'))
            {
              stringstream line_string (line_in_file);
              line_string << line_in_file;
              line_string >> par_name;   // read the first character, the name of the parameter
              line_string >> equality; 	  // check that second character is "="
              if (equality != "=")
                {
                  cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
                  cerr << "Using a default value for " << par_name << endl; exit(1);
                }
              string par_value0;
              int ending=0;
              int n_m_lims=0;
              do{
                line_string >> par_value0;  //read value
                n_m_lims++;
                if(par_value0=="END")
                  ending=1;
                this->cwv_used.push_back(atof(par_value0.c_str()));
              }while(ending!=1);
              this->cwv_used.pop_back(); //delete the last element, which is END
            }
#endif
	  else  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'x' && line_in_file[1] == '_'))
	    {
	      stringstream line_string (line_in_file);
	      line_string << line_in_file;
	      line_string >> par_name;   // read the first character, the name of the parameter
	      line_string >> equality; 	  // check that second character is "="
	      if (equality != "=") {
	       	cerr << "Error in parameters file "<<file<<" at " << par_name << " = " << "???" << endl;
		      cerr << "Please check parameter file "<< std::endl; 
      		exit(1);
	      }
	      string par_value0;
	      int ending=0;
	      do{
    		line_string >> par_value0;  //read value
  	   	if(par_value0=="END")ending=1;
    	 	this->output_at_iteration.push_back(atof(par_value0.c_str()));
	      }while(ending!=1);
	      this->output_at_iteration.pop_back(); //delete the last element, which is END
	    }
#ifdef _USE_TWO_REFS_MOCKS_

#ifdef _SLICS_

	  else  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'l' && line_in_file[1] == '_'))
	    {
	      stringstream line_string (line_in_file);
	      line_string << line_in_file;
	      line_string >> par_name;   // read the first character, the name of the parameter
	      line_string >> equality;    // check that second character is "="
	      if (equality != "=") {
	       	cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
		      cerr << "Using a default value for " << par_name << endl; exit(1);
	      }
	      string par_value0;
	      int ending=0;
	      do{
		line_string >> par_value0;  //read value
		if(par_value0=="END")ending=1;
		this->list_new_dm_fields.push_back(atoi(par_value0.c_str()));
	      }while(ending!=1);
	      this->list_new_dm_fields.pop_back(); //delete the last element, which is END
	    }
	  else  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'g' && line_in_file[1] == '_'))
	    {
	      stringstream line_string (line_in_file);
	      line_string << line_in_file;
	      line_string >> par_name;   // read the first character, the name of the parameter
	      line_string >> equality;    // check that second character is "="
	      if (equality != "=") {
    	   	cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
		      cerr << "Using a default value for " << par_name << endl; exit(1);
	      }
	      string par_value0;
	      int ending=0;
	      do{
		      line_string >> par_value0;  //read value
		      if(par_value0=="END")ending=1;
		        this->list_bias_references.push_back(atoi(par_value0.c_str()));
	      }while(ending!=1);
	      this->list_bias_references.pop_back(); //delete the last element, which is END
	    }

	  // These arrays are inizialized here and updated in teh Bam.pp file 
	  this->Number_of_references=this->list_bias_references.size();
	  this->files_bias_references.resize(this->Number_of_references);
	  for(int i=0;i<this->Number_of_references;++i)
	    this->files_bias_references[i]="../Output_SLICS/CAL_R"+to_string(this->list_bias_references[i])+"new_IKWEB/Bam_Bias.dat";
	  this->files_kernel_references.resize(this->Number_of_references);
	  for(int i=0;i<this->Number_of_references;++i)
	    this->files_kernel_references[i]="../Output_SLICS/CAL_R"+to_string(this->list_bias_references[i])+"new_IKWEB/Bam_Kernel.dat";
	  this->files_tracer_references.resize(this->Number_of_references);
	  for(int i=0;i<this->Number_of_references;++i)
	    this->files_tracer_references[i]=this->Input_dir_cat+"Halos_SLICS_LOS"+to_string(this->list_bias_references[i])+".txt";
	  this->files_tracer_field_references.resize(this->Number_of_references);
	  for(int i=0;i<this->Number_of_references;++i)
	    this->files_tracer_field_references[i]=this->Input_Directory_Y+"SLICS_HALOS_LOS"+to_string(this->list_bias_references[i])+"_Nres192_MAS0.dat";
	  this->files_dm_references.resize(this->Number_of_references);
	  for(int i=0;i<this->Number_of_references;++i)
	    this->files_dm_references[i]=this->Input_Directory_X_NEW+"densDMALPTrS20.0TETCICz1.041G192V505.0S"+to_string(this->list_bias_references[i])+".dat";
	  this->Number_of_new_mocks=this->list_new_dm_fields.size();
	  this->files_new_dm_fields.resize(this->Number_of_new_mocks);
	  for(int i=0;i<this->Number_of_new_mocks;++i)
	    this->files_new_dm_fields[i]=this->Input_Directory_X_NEW+"densDMALPTrS20.0TETCICz1.041G192V505.0S"+to_string(this->list_new_dm_fields[i])+".dat";

#elif defined _UNITSIM_
      this->Number_of_references=2;
      this->files_bias_references.resize(this->Number_of_references);
      this->files_bias_references[0]=this->Input_Directory_BIAS_KERNEL+"Bam_Bias.dat";
      this->files_bias_references[1]=this->Input_Directory_BIAS_KERNEL_TWO+"Bam_Bias.dat";
      this->files_kernel_references.resize(this->Number_of_references);
      this->files_kernel_references[0]=this->Input_Directory_BIAS_KERNEL+"Bam_Kernel.dat";
      this->files_kernel_references[1]=this->Input_Directory_BIAS_KERNEL_TWO+"Bam_Kernel.dat";
      this->files_tracer_references.resize(this->Number_of_references);
      this->files_tracer_references[0]=this->Input_dir_cat+this->file_catalogue;
      this->files_tracer_references[1]=this->Input_dir_cat_TWO+this->file_catalogue;
      this->files_dm_references.resize(this->Number_of_references);
      this->files_dm_references[0]=this->Input_Directory_X_REF+this->Name_Catalog_X;
      this->files_dm_references[1]=this->Input_Directory_X_REF_TWO+this->Name_Catalog_X;// FOR Unitsim, the files will have the smae name regardless normal or Inv
      this->Number_of_new_mocks=1;
      this->list_new_dm_fields.resize(this->Number_of_new_mocks);
      this->files_new_dm_fields.resize(this->Number_of_new_mocks);
      this->files_new_dm_fields[0]=this->Input_Directory_X_NEW+this->Name_Catalog_X_NEW; 
#endif


#endif




#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_
    if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'u' && line_in_file[1] == '_'))
      {
        stringstream line_string (line_in_file);
        line_string << line_in_file;
        line_string >> par_name;   // read the first character, the name of the parameter
        line_string >> equality;    // check that second character is "="
        if (equality != "=") {
          cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
          cerr << "Using a default value for " << par_name << endl; exit(1);
        }
        string par_value0;
        int ending=0;
        do{
          line_string >> par_value0;  //read value
          if(par_value0=="END")ending=1;
           this->list_Props_Threshold_MultiLevels.push_back(atof(par_value0.c_str()));
        }while(ending!=1);
        this->list_Props_Threshold_MultiLevels.pop_back(); //delete the last element, which is END
      }
    else  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'q' && line_in_file[1] == '_'))
      {
        stringstream line_string (line_in_file);
        line_string << line_in_file;
        line_string >> par_name;   // read the first character, the name of the parameter
        line_string >> equality;    // check that second character is "="
        if (equality != "=") {
          cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
          cerr << "Using a default value for " << par_name << endl; exit(1);
        }
        string par_value0;
        int ending=0;
        do{
          line_string >> par_value0;  //read value
          if(par_value0=="END")ending=1;
            this->list_Nft_MultiLevels.push_back(atoi(par_value0.c_str()));
        }while(ending!=1);
        this->list_Nft_MultiLevels.pop_back(); //delete the last element, which is END
      }
    else  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 's' && line_in_file[1] == '_'))
      {
        stringstream line_string (line_in_file);
        line_string << line_in_file;
        line_string >> par_name;   // read the first character, the name of the parameter
        line_string >> equality;    // check that second character is "="
        if (equality != "=") {
          cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
          cerr << "Using a default value for " << par_name << endl; exit(1);
        }
        string par_value0;
        int ending=0;
        do{
          line_string >> par_value0;  //read value
          if(par_value0=="END")ending=1;
          this->list_Props_Tolerance_MultiLevels.push_back(atof(par_value0.c_str()));
        }while(ending!=1);
        this->list_Props_Tolerance_MultiLevels.pop_back(); //delete the last element, which is END
      }

    // These arrays are inizialized here and updated in teh Bam.pp file 
    this->Number_of_MultiLevels = this->list_Nft_MultiLevels.size();
    this->list_Ntracers_MultiLevels.resize(this->Number_of_MultiLevels,0);
#endif

	  // global and I/O parameters for Power spectrum
	  if (par_name == "Statistics") this->statistics = par_value;
	  else if (par_name == "Output_directory") this->dir_output = par_value;
	  else if (par_name == "Catalogue_file") this->file_catalogue = par_value;
	  else if (par_name == "Input_type") this->input_type = par_value;
	  else if (par_name == "Input_type_two") this->input_type = par_value;
	  else if (par_name == "delta_grid_file") this->delta_grid_file = par_value;
	  else if (par_name == "delta_grid_file2") this->delta_grid_file2 = par_value;
	  else if (par_name == "delta_grid_file3") this->delta_grid_file3 = par_value;
	  else if (par_name == "delta_grid_file4") this->delta_grid_file4 = par_value;
	  else if (par_name == "ngal_delta") this->ngal_delta = atof(par_value.c_str());
	  else if (par_name == "measure_cross")
	    {
	      if(par_value == "true")this->measure_cross = true;
	      else if(par_value == "false")this->measure_cross = false;
	    }
	  else if (par_name == "redshift_space_coords_g")
	    {
	      if(par_value == "true")this->redshift_space_coords_g  = true;
	      else if(par_value == "false")this->redshift_space_coords_g  = false;
	    }
	  else if (par_name == "measure_cross_from_1")this->measure_cross_from_1 = atoi(par_value.c_str());
	  else if (par_name == "measure_cross_from_2")this->measure_cross_from_2 = atoi(par_value.c_str());
	  else if (par_name == "sys_of_coord_g") this->sys_of_coord_g = atoi(par_value.c_str());
	  else if (par_name == "i_coord1_g") this->i_coord1_g = atoi(par_value.c_str());
	  else if (par_name == "i_coord2_g") this->i_coord2_g = atoi(par_value.c_str());
	  else if (par_name == "i_coord3_g") this->i_coord3_g = atoi(par_value.c_str());
	  else if (par_name == "i_mean_density_g") this->i_mean_density_g = atoi(par_value.c_str());
	  else if (par_name == "i_weight1_g") this->i_weight1_g = atoi(par_value.c_str());
	  else if (par_name == "i_weight2_g") this->i_weight2_g = atoi(par_value.c_str());
	  else if (par_name == "i_weight3_g") this->i_weight3_g = atoi(par_value.c_str());
	  else if (par_name == "i_weight4_g") this->i_weight4_g = atoi(par_value.c_str());
	  else if (par_name == "vel_units_g") this->vel_units_g = par_value;
	  else if (par_name == "i_v1_g") this->i_v1_g = atoi(par_value.c_str());
	  else if (par_name == "i_v2_g") this->i_v2_g = atoi(par_value.c_str());
	  else if (par_name == "i_v3_g") this->i_v3_g = atoi(par_value.c_str());
	  else if (par_name == "i_mass_g") this->i_mass_g = atoi(par_value.c_str());
	  else if (par_name == "i_vmax_g") this->i_vmax_g = atoi(par_value.c_str());
	  else if (par_name == "i_sf_g") this->i_sf_g = atoi(par_value.c_str());
	  else if (par_name == "i_rs_g") this->i_rs_g = atoi(par_value.c_str());
	  else if (par_name == "i_virial_g") this->i_virial_g = atoi(par_value.c_str());
	  else if (par_name == "i_spin_g") this->i_spin_g = atoi(par_value.c_str());
	  else if (par_name == "sys_of_coord_dm") this->sys_of_coord_dm = atoi(par_value.c_str());
	  else if (par_name == "i_coord1_dm") this->i_coord1_dm = atoi(par_value.c_str());
	  else if (par_name == "i_coord2_dm") this->i_coord2_dm = atoi(par_value.c_str());
	  else if (par_name == "i_coord3_dm") this->i_coord3_dm = atoi(par_value.c_str());
	  else if (par_name == "i_v1_dm") this->i_v1_dm = atoi(par_value.c_str());
	  else if (par_name == "i_v2_dm") this->i_v2_dm = atoi(par_value.c_str());
	  else if (par_name == "i_v3_dm") this->i_v3_dm = atoi(par_value.c_str());
	  else if (par_name == "weight_vel_with_mass")
	    {
	      if (par_value=="true")this->weight_vel_with_mass = true;
	      else this->weight_vel_with_mass = false;
	    }
	  else if (par_name == "weight_with_mass")
	    {
	      if (par_value=="true")this->weight_with_mass = true;
	      else this->weight_with_mass = false;
	    }
	  else if (par_name == "use_weight1_g")
	    {
	      if(par_value == "true")this->use_weight1_g = true;
	      else if(par_value == "false")this->use_weight1_g = false;
	    }
	  else if (par_name == "use_weight2_g")
	    {
	      if(par_value == "true")this->use_weight2_g = true;
	      else if(par_value == "false")this->use_weight2_g = false;
	    }
	  else if (par_name == "use_weight3_g")
	    {
	      if(par_value == "true")this->use_weight3_g = true;
	      else if(par_value == "false")this->use_weight3_g = false;
	    }
	  else if (par_name == "use_weight4_g") {
	    if(par_value == "true")this->use_weight4_g = true;
	    else if(par_value == "false")this->use_weight4_g = false;
	  }
	  else if (par_name == "i_property1_g") this->i_property1_r = atoi(par_value.c_str());
	  else if (par_name == "i_property2_g") this->i_property2_r = atoi(par_value.c_str());
	  else if (par_name == "i_property3_g") this->i_property3_r = atoi(par_value.c_str());
	  else if (par_name == "i_property4_g") this->i_property4_r = atoi(par_value.c_str());
	  else if (par_name == "i_property5_g") this->i_property5_r = atoi(par_value.c_str());
	  else if (par_name == "i_property6_g") this->i_property6_r = atoi(par_value.c_str());
	  else if (par_name == "new_Lbox"){
	    if (par_value == "true"){this->new_Lbox = true;}
	    else if (par_value == "false"){this->new_Lbox = false;}
	  }
	  else if (par_name == "n_catalogues")this->n_catalogues = atoi(par_value.c_str());
	  else if (par_name == "angles_units_g")
	    {
	      if (par_value=="D") this->angles_units_g = "D";
	      else if (par_value=="R") this->angles_units_g = "R";
	      else {cout<<"Unknown angle units"<<std::endl; exit(1) ;}
	    }
	  else if (par_name == "use_random_catalog")
	    {
	      if (par_value=="false") this->use_random_catalog = false;
	      else if (par_value=="true") this->use_random_catalog = true;
	    }
	  else if (par_name == "use_random_catalog_cl")
	    {
	      if (par_value=="no") this->use_random_catalog_cl = false;
	      else if (par_value=="yes") this->use_random_catalog_cl = true;
	    }
	  else if (par_name == "sys_of_coord_r") this->sys_of_coord_r = atoi(par_value.c_str());
	  else if (par_name == "i_coord1_r") this->i_coord1_r = atoi(par_value.c_str());
	  else if (par_name == "i_coord2_r") this->i_coord2_r = atoi(par_value.c_str());
	  else if (par_name == "i_coord3_r") this->i_coord3_r = atoi(par_value.c_str());
	  else if (par_name == "i_mean_density_r") this->i_mean_density_r = atoi(par_value.c_str());
	  else if (par_name == "i_weight1_r") this->i_weight1_r = atoi(par_value.c_str());
	  else if (par_name == "i_weight2_r") this->i_weight2_r = atoi(par_value.c_str());
	  else if (par_name == "i_weight3_r") this->i_weight3_r = atoi(par_value.c_str());
	  else if (par_name == "i_weight4_r") this->i_weight4_r = atoi(par_value.c_str());
	  else if (par_name == "i_property1_r") this->i_property1_r = atoi(par_value.c_str());
	  else if (par_name == "i_property2_r") this->i_property2_r = atoi(par_value.c_str());
	  else if (par_name == "i_property3_r") this->i_property3_r = atoi(par_value.c_str());
	  else if (par_name == "i_property4_r") this->i_property4_r = atoi(par_value.c_str());
	  else if (par_name == "i_property5_r") this->i_property5_r = atoi(par_value.c_str());
	  else if (par_name == "i_property6_r") this->i_property6_r = atoi(par_value.c_str());
	  else if (par_name == "i_mass_r") this->i_mass_r = atoi(par_value.c_str());
	  else if (par_name == "get_distribution_min_separations")
	    {
	      if (par_value=="true")this->get_distribution_min_separations = true;
	      else this->get_distribution_min_separations = false;
	    }
	  else if (par_name == "i_mask_pixel") this->i_mask_pixel= atoi(par_value.c_str());
	  else if (par_name == "i_mask_alpha") this->i_mask_alpha= atoi(par_value.c_str());
	  else if (par_name == "i_mask_delta") this->i_mask_delta= atoi(par_value.c_str());
	  else if (par_name == "i_mask_flag") this->i_mask_flag= atoi(par_value.c_str());
	  else if (par_name == "use_weight1_r")
	    {
	      if(par_value == "true")this->use_weight1_r = true;
	      else if (par_value == "false")this->use_weight1_r = false;
	    }

	  else if (par_name == "use_weight2_r")
	    {
	      if(par_value == "true")this->use_weight2_r = true;
	      else if (par_value == "false")use_weight2_r = false;
	    }
	  else if (par_name == "use_weight3_r") {
	    if(par_value == "true")use_weight3_r = true;
	    else if (par_value == "false")use_weight3_r = false;
	  }
	  else if (par_name == "use_weight4_r")
	    {
	      if(par_value == "true")use_weight4_r = true;
	      else if (par_value == "false")use_weight4_r = false;
	    }
	  else if (par_name == "angles_units_r")
	    {
	      if (par_value == "D") angles_units_r = "D";
	      else if (par_value == "R") angles_units_r = "R";
	      else { cout<<"Unknown angle units"<<std::endl; exit(1); }
	    }
	  else if (par_name == "Random_file") file_random = par_value;
	  else if (par_name == "Name_survey") Name_survey = par_value;
	  else if (par_name == "new_los")
	    {
	      if (par_value == "false") new_los = false;
	      else if (par_value == "true") new_los = true;
	    }
	  else if (par_name == "Nft") this->Nft = atoi(par_value.c_str());
	  else if (par_name == "Nft_low") this->Nft_low = atoi(par_value.c_str());
	  else if (par_name == "Nft_HR") this->Nft_HR = atoi(par_value.c_str());
	  else if (par_name == "Nft_random_collapse") this->Nft_random_collapse = atoi(par_value.c_str());
	  else if (par_name == "Lbox") Lbox = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "Distance_fraction") this->Distance_fraction = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "velbias_random") this->velbias_random = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "velbias_dm") this->velbias_dm = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "Lbox_low") this->Lbox_low = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "mass_assignment_scheme")
	    {
	      if (par_value=="NGP") this->mass_assignment_scheme = "NGP";
	      else if (par_value=="CIC") this->mass_assignment_scheme = "CIC";
	      else if (par_value=="TSC") this->mass_assignment_scheme = "TSC";
	      else if (par_value=="PCS") this->mass_assignment_scheme = "PCS";
	      else{cout<<"Unknown mass assigment scheme"<<std::endl;exit(1);}
	    }
	  else if (par_name == "MAS_correction")
	    {
	      if (par_value=="true") MAS_correction = true;
	      else if (par_value=="false") MAS_correction = false;
	    }
	  else if (par_name == "type_of_binning")
	    {
	      if (par_value=="linear") type_of_binning = "linear";
	      else if (par_value=="log") type_of_binning = "log";
	      else{cout<<"Unknown type of binning"<<std::endl; exit(1);}
	    }
	  else if (par_name == "DeltaKmin") this->DeltaKmin = atof(par_value.c_str());
	  else if (par_name == "N_log_bins") this->N_log_bins = atoi(par_value.c_str());
	  else if (par_name == "ndel_data") this->ndel_data = atoi(par_value.c_str());
	  else if (par_name == "ndel_window") this->ndel_window = atoi(par_value.c_str());
	  else if (par_name == "N_mu_bins") this->N_mu_bins = atoi(par_value.c_str());
	  else if (par_name == "FKP_weight"){
	    if (par_value=="false") this->FKP_weight = false;
	    else if (par_value=="true") this->FKP_weight = true;
	  }
	  else if (par_name == "Pest") this->Pest = atoi(par_value.c_str());
	  else if (par_name == "SN_correction")
	    {
	      if (par_value=="false") this->SN_correction = false;
	      else if (par_value=="true") this->SN_correction = true;
	    }
	  else if (par_name == "FKP_error_bars")
	    {
	      if (par_value=="false") this->FKP_error_bars = false;
	      else if (par_value=="true") this->FKP_error_bars = true;
	    }
	  else if (par_name == "FKP_error_bars_exact")
	    {
	      if (par_value=="false")this->FKP_error_bars_exact = false;
	      else if (par_value=="true") this->FKP_error_bars_exact = true;
	    }
	  else if (par_name == "nbar_tabulated"){
	    if (par_value=="false") nbar_tabulated = false;
	    else if (par_value=="true") nbar_tabulated = true;
	  }
	  else if (par_name == "constant_depth"){
	    if (par_value=="no") constant_depth = false;
	    else if (par_value=="yes") constant_depth = true;
	  }
	  else if (par_name == "N_z_bins") N_z_bins = atoi(par_value.c_str());
	  else if (par_name == "redshift_min_sample") redshift_min_sample = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "redshift_max_sample") redshift_max_sample = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "N_dndz_bins") N_dndz_bins = atoi(par_value.c_str());
	  else if (par_name == "new_N_dndz_bins") new_N_dndz_bins = atoi(par_value.c_str());
	  else if (par_name == "area_survey") area_survey = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "Healpix_resolution") Healpix_resolution = atoi(par_value.c_str());
	  else if (par_name == "file_power") file_power = par_value;
	  else if (par_name == "file_power_log") file_power_log = par_value;

	  else if (par_name == "file_window") file_window = par_value;
	  else if (par_name == "file_dndz") file_dndz = par_value;
	  else if (par_name == "file_power2d") file_power2d = par_value;
	  else if (par_name == "file_power2d_mk") file_power2d_mk = par_value;
	  else if (par_name == "file_bispectrum") file_bispectrum = par_value;
	  //output files for Yamamoto estimator of moments
	  else if (par_name == "kmax_y_ds")kmax_y_ds = static_cast<real_prec>(atof(par_value.c_str()));

	  //Parameters for the bispectrum
	  else if (par_name == "kmax_bk")kmax_bk = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "kmin_bk")kmin_bk = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "use_fundamental_mode_as_kmin_bk"){
	    if (par_value=="no") use_fundamental_mode_as_kmin_bk = false;
	    else if (par_value=="yes") use_fundamental_mode_as_kmin_bk = true;
	  }
#ifndef _USE_COSMO_PARS_
	  else if (par_name == "om_matter") this->om_matter = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "om_radiation") this->om_radiation = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "om_baryons") this->om_baryons = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "om_vac") this->om_vac = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "om_k") this->om_k = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "Hubble") this->Hubble = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "hubble") this->hubble = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "spectral_index") this->spectral_index = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "w_eos") this->w_eos = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "N_eff") this->N_eff = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "sigma8") this->sigma8 = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "Tcmb") this->Tcmb = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "RR") this->RR = static_cast<real_prec>(atof(par_value.c_str()));
#endif

	  //here now read the aprams from patchy
	  else if (par_name == "inputmode") this->inputmode = static_cast<int>(atof(par_value.c_str()));
	  else if (par_name == "seed") this->seed = static_cast<int>(atof(par_value.c_str()));
	  else if (par_name == "seed_ref") this->seed_ref = static_cast<int>(atof(par_value.c_str()));
	  else if (par_name == "runsim"){
	    if (par_value=="false") this->runsim = false;
	    else if (par_value=="true") this->runsim = true;
	  }
	  else if (par_name == "runv"){
	    if (par_value=="false") this->runv = false;
	    else if (par_value=="true") this->runv = true;
	  }
	  else if (par_name == "diffcosmorz"){
	    if (par_value=="false") this->diffcosmorz = false;
	    else if (par_value=="true") this->diffcosmorz = true;
	  }
          else if (par_name == "ic_alias_corrected"){
            if (par_value=="false") this->ic_alias_corrected = false;
            else if (par_value=="true") this->ic_alias_corrected = true;
          }
	  else if (par_name == "Initial_Redshift_ic_power_file") this->Initial_Redshift_ic_power_file = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "") this->ic_power_file = par_value;
	  else if (par_name == "ic_power_file") this->ic_power_file = par_value;
	  else if (par_name == "ic_WN_file") this->ic_WN_file = par_value;
	  else if (par_name == "ic_file") this->ic_file = par_value;
	  else if (par_name == "ic_input_type") this->ic_input_type = par_value;
	  else if (par_name == "ic_WN_dir") this->ic_WN_dir = par_value;
	  else if (par_name == "dir") this->dir = par_value;
	  else if (par_name == "lognden"){
	    if (par_value=="false") this->lognden = false;
	    else if (par_value=="true") this->lognden = true;
	  }
	  else if (par_name == "use_ic_file"){
	    if (par_value=="false") this->use_ic_file = false;
	    else if (par_value=="true") this->use_ic_file = true;
	  }
	  else if (par_name == "transf"){
	    if (par_value=="false") this->transf = false;
	    else if (par_value=="true") this->transf = true;
	  }

	  else if (par_name == "readPS"){
	    if (par_value=="false") this->readPS = false;
	    else if (par_value=="true") this->readPS = true;
	  }
	  else if (par_name == "use_vel_kernel"){
	    if (par_value=="false") this->use_vel_kernel = false;
	    else if (par_value=="true") this->use_vel_kernel = true;
	  }
	  else if (par_name == "Nchunk") this->Nchunk = static_cast<int>(atof(par_value.c_str()));
    else if (par_name == "vkernel_exponent") this->vkernel_exponent = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "sfmodel") this->sfmodel = (int)atof(par_value.c_str());
	  else if (par_name == "masskernel") this->masskernel = (int)atof(par_value.c_str());
	  else if (par_name == "masskernel_vel") this->masskernel_vel = (int)atof(par_value.c_str());
	  else if (par_name == "biasE") this->biasE = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "biasepsilon") this->biasepsilon = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "biasepsilon2") this->biasepsilon2 = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "biasrhoexp") this->biasrhoexp = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "biasrhoexp2") this->biasrhoexp2 = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "biasone") this->biasone = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "biassign") this->biassign = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "biassign2") this->biassign2 = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "devpois") this->devpois = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "deltathH") this->deltathH = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "deltath") this->deltath = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "Nmean") this->Nmean = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "cst") this->cst = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "cpsi") this->cpsi = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "biasL") this->biasL = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "sfac") this->sfac = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "ep") this->ep = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "slength") this->slength = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "slengthv") this->slengthv = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "velbias") this->velbias = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "vslength") this->vslength = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "cdeltas2") this->cdeltas2 = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "xllc") this->xllc = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "yllc") this->yllc = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "zllc") this->zllc = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "xobs") this->xobs = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "yobs") this->yobs = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "zobs") this->zobs = static_cast<real_prec>(atof(par_value.c_str()));
	}

    }


    int ml_aux=0;
#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_
#ifdef _USE_MULTISCALE_LEVEL_1_
    ml_aux++;
#endif
#ifdef _USE_MULTISCALE_LEVEL_2_
    ml_aux++;
#endif
#ifdef _USE_MULTISCALE_LEVEL_3_
    ml_aux++;
#endif
#ifdef _USE_MULTISCALE_LEVEL_4_
    ml_aux++;
#endif
    this->Number_of_MultiLevels =ml_aux;
#endif


#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_
if(this->Number_of_MultiLevels > this->list_Props_Tolerance_MultiLevels.size() || this->Number_of_MultiLevels < this->list_Props_Tolerance_MultiLevels.size()){
   cerr << "Error in parameters for Multilevel. Not the same amount of entries" << endl;
   cerr << "Please check input parameter file." << endl;
   cerr<<RED<<this->Number_of_MultiLevels<<"  "<<this->list_Props_Tolerance_MultiLevels.size()<<endl;
   exit(1);
  }
#endif


#ifdef _USE_TWO_REFS_MOCKS_
  if(Number_of_references<Number_of_new_mocks){
    cout<<RED<<"Warning: NUmber of references smaller than number of mocks to build simulstaneously"<<RESET<<std::endl;
    exit(0);
  }
#endif

#ifdef _USE_COSMO_PARS_
  this->om_matter = COSMOPARS::Om_matter;
  this->om_radiation = COSMOPARS::Om_radiation;
  this->om_baryons = COSMOPARS::Om_baryons;
  this->om_cdm = this->om_matter-this->om_baryons;
  this->om_vac = COSMOPARS::Om_vac;
  this->om_k   = COSMOPARS::Om_k;
  this->Hubble = COSMOPARS::Hubble;
  this->hubble = COSMOPARS::hubble;
  this->n_s = COSMOPARS::n_s;
  this->w_eos = COSMOPARS::w_eos;
  this->N_eff =  COSMOPARS::N_eff;
  this->sigma8 = COSMOPARS::sigma8;
  this->Tcmb = COSMOPARS::Tcmb;
  this->use_wiggles = true;
  this->RR = COSMOPARS::RR;
  this->alpha_s=COSMOPARS::alpha_s;
  this->Delta_SO=COSMOPARS::Delta_SO;
#endif

  this->s_cosmo_pars.cosmological_redshift=this->redshift;
  this->s_cosmo_pars.Hubble=this->Hubble;
  this->s_cosmo_pars.hubble=this->hubble;
  this->s_cosmo_pars.Om_matter =this->om_matter;
  this->s_cosmo_pars.Om_cdm =this->om_cdm;
  this->s_cosmo_pars.Om_baryons =this->om_baryons;
  this->s_cosmo_pars.Om_radiation =this->om_radiation;
  this->s_cosmo_pars.Om_vac =this->om_vac;
  this->s_cosmo_pars.Om_k =this->om_k;
  this->s_cosmo_pars.n_s =this->spectral_index;
  this->s_cosmo_pars.w_eos =this->w_eos;
  this->s_cosmo_pars.N_eff =this->N_eff;
  this->s_cosmo_pars.sigma8 =this->sigma8;
  this->s_cosmo_pars.f_baryon =this->om_baryons/this->om_matter;
  this->s_cosmo_pars.use_wiggles =this->use_wiggles;
  this->s_cosmo_pars.RR =this->RR;



  this->d1=this->Lbox/static_cast<real_prec>(Nft);		/* grid spacing x-direction */
  this->d2=this->Lbox/static_cast<real_prec>(Nft);		/* grid spacing y-direction */
  this->d3=this->Lbox/static_cast<real_prec>(Nft);		/* grid spacing z-direction */

#ifdef _USE_MULTISCALE_LEVEL_4_
  this->d1_low=this->Lbox/static_cast<real_prec>(Nft_low_l4);		/* grid spacing x-direction */
  this->d2_low=this->Lbox/static_cast<real_prec>(Nft_low_l4);		/* grid spacing y-direction */
  this->d3_low=this->Lbox/static_cast<real_prec>(Nft_low_l4);		/* grid spacing z-direction */
#endif

#ifdef _USE_X_ROOT_MASS_
  this->LOGMASSmax*=exponent_mass_tracer;
  this->LOGMASSmin*=exponent_mass_tracer;
#endif

#ifndef _USE_CWC_
#ifdef _USE_MASS_KNOTS_
  this->cwt_used.push_back(1);
#elif !defined  _USE_MASS_KNOTS_
  this->cwt_used.push_back(0);
#endif
#endif
  this->n_cwt=this->cwt_used.size();

#ifdef _USE_CWC_
  if(this->n_cwt==1)
    {
      if(this->cwt_used[0]==0)
    	{
	     cout<<RED<<"Warning: _USE_CWC_ is defined but not specified in parameter file"<<RESET<<std::endl;
	  exit(1);
	}
    }
#endif

#ifndef _USE_CW_V
#ifdef _USE_VEL_KNOTS_V_
  this->cwv_used.push_back(1);
#elif !defined  _USE_VEL_KNOTS_V_
  this->cwv_used.push_back(0);
#endif
#endif

  this->n_cwv=this->cwv_used.size();

#ifdef _USE_CWC_V_
  if(this->n_cwv==1)
    {
      if(this->cwv_used[0]==0){
          cout<<RED<<"Warning: _USE_CWC_V_ is defined but not specified in parameter file"<<RESET<<std::endl;
          exit(0);
        }
    }
#endif

#ifndef _USE_V_DISP_VEL
  this->n_vknot_massbin=1;
#endif

#ifndef _USE_MASS_KNOTS_
  this->n_sknot_massbin=1;
#endif

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

void Params::derived_pars(){
  this->NGRID = static_cast<ULONG>(this->Nft*this->Nft*this->Nft);
  this->NGRID_h = static_cast<ULONG>(this->Nft*this->Nft*(this->Nft/2+1));  // h is for half
  this->delta_x = this->Lbox/ (static_cast<double>(this->Nft)); 
  this->delta_y = this->Lbox/ (static_cast<double>(this->Nft)); 
  this->delta_z = this->Lbox/ (static_cast<double>(this->Nft)); 
  this->deltak_x = 2.*M_PI/static_cast<double>(this->Lbox);
  this->deltak_y = 2.*M_PI/static_cast<double>(this->Lbox);
  this->deltak_z = 2.*M_PI/static_cast<double>(this->Lbox);
  this->deltak_0 = sqrt(pow(this->deltak_x,2)+pow(this->deltak_y,2)+pow(this->deltak_z,2))/sqrt(3.0); 
  this->kmin     = this->DeltaKmin*this->deltak_0;
  this->kmax     = sqrt(pow(0.5*this->Nft*this->deltak_x,2)+pow(0.5*this->Nft*this->deltak_y,2)+pow(0.5*this->Nft*this->deltak_z,2))/sqrt(3.0);         
  this->DeltaK_data   = this->ndel_data*this->deltak_0; 
  this->DeltaK_window = this->ndel_window*this->deltak_0; 
  this->Nnp_data      = (this->type_of_binning=="log" ? this->N_log_bins : this->Nft/this->ndel_data/2); 
  this->Nnp_window    = (this->type_of_binning=="log" ? this->N_log_bins : this->Nft/this->ndel_window/ 2); 
  this->Deltal        = log10(this->kmax/this->kmin)/static_cast<double>(this->N_log_bins); 
  this->Deltamu       = 2.0/(static_cast<real_prec>(this->N_mu_bins));

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


void Params::show_params()
{

  std::cout<<BOLDYELLOW;
  std::cout<<"***************************************************************************"<<std::endl;
  std::cout<<"***************************************************************************"<<std::endl;
  std::cout<<"BAM                                                                       *"<<std::endl;
  std::cout<<"Bias Assignment Method for galaxy/halo mock catalogs                      *"<<std::endl;
  std::cout<<"Input values of parameters in parameter file                              *"<<std::endl;
  std::cout<<"***************************************************************************"<<std::endl;
  std::cout<<"***************************************************************************"<<RESET<<std::endl;

  std::cout<<"Redshift                   = "<<this->redshift<<std::endl;
  std::cout<<"Lbox                       = "<<this->Lbox<<std::endl;
  std::cout<<"NX                         = "<<this->NX<<std::endl;
  std::cout<<"NY                         = "<<this->NY<<std::endl;
  std::cout<<"Output_directory           = "<<this->Output_directory<<std::endl;
  std::cout<<"Input_Directory_X          = "<<this->Input_Directory_X<<std::endl;
  std::cout<<"Input_Directory_X_REF      = "<<this->Input_Directory_X_REF<<std::endl;
  std::cout<<"Input_Directory_X_REF_TWO  = "<<this->Input_Directory_X_REF_TWO<<std::endl;
  std::cout<<"Input_Directory_X_NEW      = "<<this->Input_Directory_X_NEW<<std::endl;
  std::cout<<"XNAME                      = "<<this->XNAME<<std::endl;
  std::cout<<"Name_Catalog_X             = "<<this->Name_Catalog_X<<std::endl;
  std::cout<<"Name_redshift_mask         = "<<this->Name_redshift_mask<<std::endl;
  std::cout<<"Name_Catalog_X_REF_PDF     = "<<this->Name_Catalog_X_REF_PDF<<std::endl;
  std::cout<<"Name_Catalog_X_NEW         = "<<this->Name_Catalog_X_NEW<<std::endl;
  std::cout<<"Name_Property_X            = "<<this->Name_Property_X<<std::endl;
  std::cout<<"iMAS_X                     = "<<this->iMAS_X<<std::endl;
  std::cout<<"iMAS_X_REF_PDF             = "<<this->iMAS_X_REF_PDF<<std::endl;
  std::cout<<"iMAS_X_NEW                 = "<<this->iMAS_X_NEW<<std::endl;
  std::cout<<"Input_Directory_Y          = "<<this->Input_Directory_Y<<std::endl;
  std::cout<<"YNAME                      = "<<this->YNAME<<std::endl;
  std::cout<<"Name_Catalog_Y             = "<<this->Name_Catalog_Y<<std::endl;
  std::cout<<"Name_Catalog_Y_MWEIGHTED   = "<<this->Name_Catalog_Y_MWEIGHTED<<std::endl;
  std::cout<<"Name_Catalog_Y_HR          = "<<this->Name_Catalog_Y_HR<<std::endl;
  std::cout<<"Name_Property_Y            = "<<this->Name_Property_Y<<std::endl;
  std::cout<<"iMAS_Y                     = "<<this->iMAS_Y<<std::endl;
  std::cout<<"delta_Y_max                = "<<this->delta_Y_max<<std::endl;
  std::cout<<"delta_Y_min                = "<<this->delta_Y_min<<std::endl;
  std::cout<<"delta_X_max                = "<<this->delta_X_max<<std::endl;
  std::cout<<"delta_X_min                = "<<this->delta_X_min<<std::endl;
  std::cout<<"ldelta_Y_max               = "<<this->ldelta_Y_max<<std::endl;
  std::cout<<"ldelta_Y_min               = "<<this->ldelta_Y_min<<std::endl;
  std::cout<<"ldelta_X_max               = "<<this->ldelta_X_max<<std::endl;
  std::cout<<"ldelta_X_min               = "<<this->ldelta_X_min<<std::endl;
  std::cout<<"Quantity                   = "<<this->Quantity<<std::endl;
  std::cout<<"NMASSbins                  = "<<this->NMASSbins<<std::endl;
  std::cout<<"smscale                    = "<<this->smscale<<std::endl;
  std::cout<<"realization                = "<<this->realization<<std::endl;
  std::cout<<"Nft                        = "<<this->Nft<<std::endl;
  std::cout<<"ndel_data                  = "<<this->ndel_data<<std::endl;
  std::cout<<"write_files_for_histograms = "<<this->write_files_for_histograms<<std::endl;
  std::cout<<"Redefine_limits            = "<<this->Redefine_limits<<std::endl;
  std::cout<<"Convert_Density_to_Delta_X = "<<this->Convert_Density_to_Delta_X<<std::endl;
  std::cout<<"Convert_Density_to_Delta_Y = "<<this->Convert_Density_to_Delta_Y<<std::endl;
  std::cout<<"lambdath                   = "<<this->lambdath<<std::endl;
  std::cout<<"Write_PDF_number_counts    = "<<this->Write_PDF_number_counts<<std::endl;
  std::cout<<"Scale_X                    = "<<this->Scale_X<<std::endl;
  std::cout<<"Scale_Y                    = "<<this->Scale_Y<<std::endl;
  std::cout<<"n_sknot_massbin            = "<<this->n_sknot_massbin<<std::endl;
  std::cout<<"n_vknot_massbin            = "<<this->n_vknot_massbin<<std::endl;
  std::cout<<"N_iterations_Kernel        = "<<this->N_iterations_Kernel<<std::endl;
  std::cout<<"N_iterations_dm            = "<<this->N_iterations_dm<<std::endl;
  std::cout<<"N_dm_realizations          = "<<this->N_dm_realizations<<std::endl;
  std::cout<<"n_cwt                      = "<<this->n_cwt<<std::endl;
  std::cout<<"n_cwv                      = "<<this->n_cwv<<std::endl;
  std::cout<<"ndel_data                  = "<<this->ndel_data<<std::endl;
  std::cout<<"binning type Fourier       = "<<this->type_of_binning<<std::endl;
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


