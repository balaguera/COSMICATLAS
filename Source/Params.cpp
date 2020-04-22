#include "../Headers/Params.h"
#include "../Headers/cosmo_parameters.h"

void Params::init_pars()
{
  
  this->NX = 200;
  this->NY = 100;
  this->n_sknot_massbin = 200;
  this->n_vknot_massbin = 200;
  this->NMASSbins = 1;
  this->iMAS_X = 0;
  this->iMAS_X_REF_PDF = 0;
  this->iMAS_X_NEW = 0;
  this->iMAS_Y = 0;
  this->lambdath = 0.0;
  this->lambdath_v = 0.0;
  this->realization = 1;
  this->redshift = 0.0;
  this->Redshift_initial;
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
  this->XNAME = "XNAME";
  this->YNAME = "YNAME";
  this->Name_Catalog_X = "UNKNOWN_NAME";
  this->Name_VelFieldx_X = "UNKNOWN_NAME";
  this->Name_VelFieldy_X = "UNKNOWN_NAME";
  this->Name_VelFieldz_X = "UNKNOWN_NAME";
  this->Name_Catalog_X_REF_PDF = "UNKNOWN_NAME";
  this->Name_Catalog_X_NEW = "UNKNOWN_NAME";
  this->extra_info = "";
  this->Name_Property_X = "X_PROEPRTY";
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

  this->Redefine_limits= false;
  this->Write_Scatter_Plot= false;
  this->Write_PDF_number_counts= false;
  this->Convert_Density_to_Delta_X= false;
  this->Convert_Density_to_Delta_Y= false;
  this->Comp_joint_PDF = false;
  this->write_files_for_histograms = false;
  //this->cwt_used.push_back(0);
  this->n_cwt=ONE;
  this->n_cwv=ONE;
  //  this->output_at_iteration.push_back(0);
  
  // global and I/O parameters for Power spectrum
  this->statistics = "fkp";
  this->dir_input = "../Input/";
  this->dir_output = "../Output/";
  this->file_catalogue = "cat.dat";
  this->input_type = "catalog";
  this->delta_grid_file = "delta_file";
  this->delta_grid_file2 = "delta_file";
  this->delta_grid_file3 = "delta_file";
  this->delta_grid_file4 = "delta_file";
  this->ngal_delta = 1;
  this->measure_cross = true;
  this->measure_cross_from_1 = 1;
  this->measure_cross_from_2 = 2;
  this->sys_of_coord_g = 1;
  this->i_coord1_g = 0;
  this->i_coord2_g = 1;
  this->i_coord3_g = 2;
  this->i_mean_density_g = 3;
  this->i_weight1_g = 4;
  this->i_weight2_g = 5;
  this->i_weight3_g = 6;
  this->i_weight4_g = 7;
  this->i_coord1_dm = 0;
  this->i_coord2_dm = 1;
  this->i_coord3_dm = 2;
  this->i_mass_dm = 7;
  

  this->use_weight1_g = true;
  this->use_weight2_g = true;
  this->use_weight3_g = true;
  this->use_weight4_g = true;
  this->i_property1_r = 0;
  this->i_property2_r = 1;
  this->i_property3_r = 2;
  this->i_property4_r = 3;
  this->i_property5_r = 4;
  this->i_property6_r = 5;
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
  this->k_bin_step = 0.5;
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

 

  // cosmological parameters
  this->om_matter = Cosmo_parameters_PLANCK::Om_matter;
  this->om_radiation = Cosmo_parameters_PLANCK::Om_radiation;
  this->om_baryons = Cosmo_parameters_PLANCK::Om_baryons;
  this->om_cdm = this->om_matter-this->om_baryons;
  this->om_vac = Cosmo_parameters_PLANCK::Om_vac;
  this->om_k   = Cosmo_parameters_PLANCK::Om_k;
  this->Hubble = Cosmo_parameters_PLANCK::Hubble;
  this->hubble = Cosmo_parameters_PLANCK::hubble;
  this->spectral_index = Cosmo_parameters_PLANCK::n_s;
  this->w_eos = Cosmo_parameters_PLANCK::w_eos;
  this->N_eff =  Cosmo_parameters_PLANCK::N_eff;
  this->sigma8 = Cosmo_parameters_PLANCK::sigma8;
  this->Tcmb = Cosmo_parameters_PLANCK::Tcmb;
  this->use_wiggles = true;
  this->RR = Cosmo_parameters_PLANCK::RR;
  
  
}


// *************************************************************************************************************
// *************************************************************************************************************
// *************************************************************************************************************



void Params::read_pars(string file)
{
  
  ifstream fin_parameters (file.c_str());

  
  if (!fin_parameters) { cerr <<"Error in opening the parameters file "<<file<<"!"<<endl; exit(1); }
  
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
	      cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
	      cerr << "Using a default value for " << par_name << endl; exit(1);
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
	  if (par_name == "NX")this->NX = atoi(par_value.c_str());
	  else if (par_name == "NY")this->NY = atoi(par_value.c_str());
	  else if (par_name == "NY_MASS")this->NY_MASS = atoi(par_value.c_str());
	  else if (par_name == "NY_SAT_FRAC")this->NY_SAT_FRAC = atoi(par_value.c_str());
	  else if (par_name == "N_SKNOT_MASSBIN")this->n_sknot_massbin = atoi(par_value.c_str());
          else if (par_name == "N_VKNOT_MASSBIN")this->n_vknot_massbin = atoi(par_value.c_str());
          else if (par_name == "NMASSbins")this->NMASSbins = atoi(par_value.c_str());
	  else if (par_name == "iMAS_X")this->iMAS_X = atoi(par_value.c_str());
	  else if (par_name == "iMAS_X_REF_PDF")this->iMAS_X_REF_PDF = atoi(par_value.c_str());
	  else if (par_name == "iMAS_X_NEW")this->iMAS_X_NEW = atoi(par_value.c_str());
	  else if (par_name == "iMAS_Y")this->iMAS_Y = atoi(par_value.c_str());
	  else if (par_name == "lambdath")this->lambdath = static_cast<real_prec>(atof(par_value.c_str()));
          else if (par_name == "lambdath_v")this->lambdath_v = static_cast<real_prec>(atof(par_value.c_str()));
          else if (par_name == "Realization")this->realization = atoi(par_value.c_str());
	  else if (par_name == "Redshift")this->redshift = static_cast<real_prec>(atof(par_value.c_str()));
	  else if (par_name == "Redshift_initial")this->Redshift_initial = static_cast<real_prec>(atof(par_value.c_str()));
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
	  else if (par_name == "Input_Directory_X_REF")this->Input_Directory_X_REF = par_value;
	  else if (par_name == "XNAME")this->XNAME = par_value;
	  else if (par_name == "YNAME")this->YNAME = par_value;
	  else if (par_name == "Name_Catalog_X")this->Name_Catalog_X = par_value;
          else if (par_name == "Name_VelFieldx_X")this->Name_VelFieldx_X = par_value;
          else if (par_name == "Name_VelFieldy_X")this->Name_VelFieldy_X = par_value;
          else if (par_name == "Name_VelFieldz_X")this->Name_VelFieldz_X = par_value;
	  else if (par_name == "Name_redshift_mask")this->Name_redshift_mask = par_value;
	  else if (par_name == "Name_binary_mask")this->Name_binary_mask = par_value;
	  

          else if (par_name == "Name_Catalog_X_REF_PDF")this->Name_Catalog_X_REF_PDF = par_value;
	  else if (par_name == "Name_Catalog_X_NEW")this->Name_Catalog_X_NEW = par_value;
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
          else if (par_name == "LOGMASSmax") this->LOGMASSmax = atof(par_value.c_str());
          else if (par_name == "NMASSbins") this->NMASSbins = atoi(par_value.c_str());
	  else if (par_name == "NMASSbins_mf") this->NMASSbins_mf = atoi(par_value.c_str());


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
	      do{
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
		cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
		cerr << "Using a default value for " << par_name << endl; exit(1);
	      }
	      string par_value0;
	      int ending=0;
	      int n_m_lims=0;
	      do{
		line_string >> par_value0;  //read value
		n_m_lims++;
		if(par_value0=="END")ending=1;
		this->output_at_iteration.push_back(atof(par_value0.c_str()));
	      }while(ending!=1);
	      this->output_at_iteration.pop_back(); //delete the last element, which is END
	    }
	  
	  
	  // global and I/O parameters for Power spectrum
	  if (par_name == "Statistics") this->statistics = par_value;
	  else if (par_name == "Input_directory") this->dir_input = par_value;
	  else if (par_name == "Output_directory") this->dir_output = par_value;
	  else if (par_name == "Catalogue_file") this->file_catalogue = par_value;
	  else if (par_name == "Input_type") this->input_type = par_value;
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
	  else if (par_name == "i_v1_g") this->i_v1_g = atoi(par_value.c_str()); 
	  else if (par_name == "i_v2_g") this->i_v2_g = atoi(par_value.c_str());
	  else if (par_name == "i_v3_g") this->i_v3_g = atoi(par_value.c_str());
	  else if (par_name == "i_mass_g") this->i_mass_g = atoi(par_value.c_str());
          else if (par_name == "i_vmax_g") this->i_vmax_g = atoi(par_value.c_str());
          else if (par_name == "i_sf_g") this->i_sf_g = atoi(par_value.c_str());

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
	      else {cout<<"Unknown angle units"<<endl; exit(1) ;}
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
	      else { cout<<"Unknown angle units"<<endl; exit(1); }
	    }
	  else if (par_name == "Random_file") file_random = par_value;
	  else if (par_name == "Name_survey") Name_survey = par_value;
	  
	  
      // parameters for the power spectrum 
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
      else if (par_name == "Lbox_low") this->Lbox_low = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "mass_assignment_scheme")
	{
          if (par_value=="NGP") this->mass_assignment_scheme = "NGP";
          else if (par_value=="CIC") this->mass_assignment_scheme = "CIC";
          else if (par_value=="TSC") this->mass_assignment_scheme = "TSC";
          else if (par_value=="PCS") this->mass_assignment_scheme = "PCS";
	  else{cout<<"Unknown mass assigment scheme"<<endl;exit(1);}
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
	  else{cout<<"Unknown type of binning"<<endl; exit(1);}
	}
      else if (par_name == "k_bin_step") this->k_bin_step = atof(par_value.c_str());
      else if (par_name == "N_log_bins") this->N_log_bins = atoi(par_value.c_str());
      else if (par_name == "ndel_data") this->ndel_data = atoi(par_value.c_str());
      else if (par_name == "ndel_window") this->ndel_window = atoi(par_value.c_str());
      else if (par_name == "N_mu_bins") this->N_mu_bins = atoi(par_value.c_str());
      else if (par_name == "FKP_weight"){
        if (par_value=="false") this->FKP_weight = false;
        else if (par_value=="true") this->FKP_weight = true;
      }
      else if (par_name == "Pest") this->Pest = atoi(par_value.c_str());
      else if (par_name == "SN_correction"){
        if (par_value=="false") this->SN_correction = false;
        else if (par_value=="true") this->SN_correction = true;
      }
      else if (par_name == "FKP_error_bars"){
        if (par_value=="false") this->FKP_error_bars = false;
        else if (par_value=="true") this->FKP_error_bars = true;
      }
      else if (par_name == "FKP_error_bars_exact"){
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
	  
	  // cosmological parameters. NOt yet in input_delta_minerva.ini file, so read from initialziation in params
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
      else if (par_name == "ic_power_file") this->ic_power_file = par_value;
      else if (par_name == "ic_WN_file") this->ic_WN_file = par_value;
      else if (par_name == "ic_file") this->ic_file = par_value;
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
      else if (par_name == "Nchunk") this->Nchunk = (int)atof(par_value.c_str());
      else if (par_name == "sfmodel") this->sfmodel = (int)atof(par_value.c_str());
      else if (par_name == "masskernel") this->masskernel = (int)atof(par_value.c_str());
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
    
  this->d1=this->Lbox/static_cast<real_prec>(Nft);		/* grid spacing x-direction */
  this->d2=this->Lbox/static_cast<real_prec>(Nft);		/* grid spacing y-direction */
  this->d3=this->Lbox/static_cast<real_prec>(Nft);		/* grid spacing z-direction */

#ifdef _USE_MULTISCALE_LEVEL_4_
  this->d1_low=this->Lbox/static_cast<real_prec>(NFT_LOW_4);		/* grid spacing x-direction */
  this->d2_low=this->Lbox/static_cast<real_prec>(NFT_LOW_4);		/* grid spacing y-direction */
  this->d3_low=this->Lbox/static_cast<real_prec>(NFT_LOW_4);		/* grid spacing z-direction */
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
	  cout<<RED<<"Warning: _USE_CWC_ is defined but not specified in parameter file"<<RESET<<endl;
	  exit(0);
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
      if(this->cwv_used[0]==0)
        {
          cout<<RED<<"Warning: _USE_CWC_V_ is defined but not specified in parameter file"<<RESET<<endl;
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



#ifdef MOCK_MODE
//  this->Scale_Y = "linear";
//  this->Scale_X = "log";
//  this->Convert_Density_to_Delta_X= true;
  //this->Convert_Density_to_Delta_Y= false;
#endif
  
  
#ifdef BIAS_MODE
/*
  if(this->Name_Property_X=="DENSITY")
    this->Convert_Density_to_Delta_X= true;
  else
    this->Convert_Density_to_Delta_X= false;
*/

/*
  if(this->Name_Property_Y=="DENSITY")
     this->Convert_Density_to_Delta_Y= true;
   else
     this->Convert_Density_to_Delta_Y= false
*/;
  
  
//   this->Scale_Y = "log";
//   this->Scale_X = "log";

#endif
  


  if(this->mass_assignment_scheme == "PCS")
   {
#ifndef _USE_FOURTH_ORDER_
      cout<<RED<<"Please define _FOURTH_ORDER_ in def.h file"<<RESET<<endl;
      exit(0);
#else
      cout<<RED<<"_FOURTH_ORDER_ requested and defined"<<RESET<<endl;
#endif
  }

  
}


