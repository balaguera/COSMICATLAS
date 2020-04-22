// ================================================================================================
// This file describes the class Parameters, used to set all the parameters adopted in the analysis
// ================================================================================================

#ifndef __PARAMETERS__
#define __PARAMETERS__

#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
using namespace std;



class ParametersCosmolib
{
  
 private :
  

  real_prec k_min_integration;
  real_prec k_max_integration;

  real_prec om_matter;
  real_prec om_cdm;
  real_prec om_radiation;
  real_prec om_baryons;
  real_prec om_vac;
  real_prec om_k;
  real_prec f_baryon;
  real_prec Hubble;
  real_prec hubble;
  real_prec w_eos;
  real_prec N_eff;
  real_prec sigma8;
  real_prec A_s;
  real_prec n_s;
  real_prec alpha_s;
  real_prec Tcmb;
  real_prec RR;
  real_prec M_reference;
  bool use_wiggles;
  bool fixed_redshift;
  
  real_prec redshift;
  real_prec redshift_min;
  real_prec redshift_max;

  string mass_function_fit;
  
  real_prec M_min_effective;
  real_prec M_max_effective;
  real_prec A_gas;
  real_prec B_gas;
  real_prec mstar;
  
  real_prec Delta_SO;
  int hod_model;
  real_prec muno_hod;
  real_prec alpha_hod;
  real_prec mmin_hod;
  real_prec scatter_hod;

  real_prec sigma_ln;
  real_prec sigma_red;
  real_prec missing_flux;

  string density_profile;


  real_prec M_min_mf;
  real_prec M_max_mf;
  string scale_mf;
  int n_points_mf;
  string mass_function_output_file;
  string halo_mass_bias_fit;
  string halo_mass_bias_output_file;
  string effective_halo_mass_bias_output_file;
  string effective_halo_mean_number_density_output_file;
  bool compute_output_linear_power_spectrum;
  bool compute_output_non_linear_power_spectrum;
  string scale_ps;
  real_prec k_min_ps;
  real_prec k_max_ps;
  real_prec kstar;
  real_prec GAL_BIAS;
  real_prec Amc;
  int n_points_ps;
  string linear_matter_ps_output_file;
  string non_linear_matter_ps_halo_fit_output_file;
  string non_linear_matter_ps_pt_output_file;
  real_prec coef_concentration_amp;
  real_prec coef_concentration;


  string galaxy_power_spectrum_halo_model_output_file;
  string galaxy_correlation_function_halo_model_output_file;

  string scale_cf;
  real_prec r_min_cf;
  real_prec r_max_cf;
  int n_points_cf;
  string linear_matter_cf_output_file;
  string non_linear_matter_cf_halo_fit_output_file;

  bool compute_output_linear_correlation_function;
  bool compute_output_non_linear_correlation_function;
  bool compute_density_profile;
  real_prec r_min_dp;
  real_prec r_max_dp;
  string scale_dp_r;
  int n_points_dp_r;
  string density_profile_r_output_file;
  real_prec k_min_dp;
  real_prec k_max_dp;
  string scale_dp_k;
  int n_points_dp_k;
  string density_profile_k_output_file ;




 public:

  // Default constructor
  ParametersCosmolib () {}
  // Constructor: the input is the parameter file
  ParametersCosmolib (string &);
  // Destructor
  ~ParametersCosmolib () {}

  // function to get private variables
  real_prec _k_max_integration () {return k_max_integration;}
  real_prec _k_min_integration () {return k_min_integration;}
  real_prec _om_matter () {return om_matter;}
  real_prec _om_cdm () {return om_cdm;}
  real_prec _om_radiation () {return om_radiation;}
  real_prec _om_baryons () {return om_baryons;}
  real_prec _om_vac () {return om_vac;}
  real_prec _om_k () {return om_k;}
  real_prec _f_baryon () {return f_baryon;}
  real_prec _Hubble () {return Hubble;}
  real_prec _hubble () {return hubble;}
  real_prec _n_s () {return n_s;}
  real_prec _w_eos () {return w_eos;}
  real_prec _N_eff () {return N_eff;}
  real_prec _sigma8 () {return sigma8;}
  real_prec _A_s () {return A_s;}
  real_prec _alpha_s () {return alpha_s;}
  real_prec _Tcmb() {return Tcmb;}
  real_prec _RR() {return RR;}
  real_prec _Delta_SO () {return Delta_SO;}
  bool _use_wiggles() {return use_wiggles;}
  bool _fixed_redshift() {return fixed_redshift;}
  real_prec _redshift() {return redshift;}
  real_prec _redshift_min() {return redshift_min;}
  real_prec _A_gas() {return A_gas;}
  real_prec _B_gas() {return B_gas;}
  real_prec _mstar() {return mstar;}
  real_prec _kstar() {return kstar;}
  real_prec _GAL_BIAS() {return GAL_BIAS;}
  real_prec _Amc() {return Amc;}

  real_prec _coef_concentration_amp() {return coef_concentration_amp;}
  real_prec _coef_concentration() {return coef_concentration;}


  string _mass_function_fit() {return mass_function_fit;}

  real_prec _sigma_ln() {return sigma_ln;}
  real_prec _sigma_red() {return sigma_red;}
  real_prec _missing_flux() {return missing_flux;}
  string _density_profile(){return density_profile;}
  real_prec _muno_hod() {return muno_hod;}
  real_prec _mmin_hod() {return mmin_hod;}
  real_prec _scatter_hod() {return scatter_hod;}
  int _hod_model() {return hod_model;}
  real_prec _alpha_hod() {return alpha_hod;}


  string _galaxy_power_spectrum_halo_model_output_file(){return galaxy_power_spectrum_halo_model_output_file;}
  string _galaxy_correlation_function_halo_model_output_file(){return galaxy_correlation_function_halo_model_output_file;}
  real_prec _M_min_effective() {return M_min_effective;}
  real_prec _M_max_effective() {return M_max_effective;}
  real_prec _M_min_mf() {return M_min_mf;}
  real_prec _M_max_mf() {return M_max_mf;}
  string _scale_mf() {return scale_mf ;}
  int _n_points_mf() {return n_points_mf ;}
  string _mass_function_output_file() {return mass_function_output_file ;}
  string _halo_mass_bias_fit() {return halo_mass_bias_fit ;}
  string _halo_mass_bias_output_file() {return halo_mass_bias_output_file ;}
  string _effective_halo_mass_bias_output_file() {return effective_halo_mass_bias_output_file ;}

  string _effective_halo_mean_number_density_output_file() {return effective_halo_mean_number_density_output_file ;}
  bool _compute_output_linear_power_spectrum() {return compute_output_linear_power_spectrum ;}
  bool _compute_output_non_linear_power_spectrum() {return compute_output_non_linear_power_spectrum  ;}
  string _scale_ps() {return scale_ps ;}
  real_prec _k_min_ps() {return k_min_ps ;}
  real_prec _k_max_ps() {return k_max_ps ;}
  int _n_points_ps() {return n_points_ps ;}
  string _linear_matter_ps_output_file() {return linear_matter_ps_output_file  ;}
  string _non_linear_matter_ps_halo_fit_output_file() {return non_linear_matter_ps_halo_fit_output_file ;}
  string _non_linear_matter_ps_pt_output_file() {return non_linear_matter_ps_pt_output_file ;}

  string _scale_cf() {return scale_cf ;}
  real_prec _r_min_cf() {return r_min_cf ;}
  real_prec _r_max_cf() {return r_max_cf ;}
  int _n_points_cf() {return n_points_cf;}
  string _linear_matter_cf_output_file() {return linear_matter_cf_output_file ;}
  string _non_linear_matter_cf_halo_fit_output_file() {return non_linear_matter_cf_halo_fit_output_file ;}

  bool _compute_density_profile() {return compute_density_profile ;}
  real_prec _r_min_dp() {return r_min_dp ;}
  real_prec _r_max_dp() {return r_max_dp ;}
  string _scale_dp_r() {return scale_dp_r ;}
  int _n_points_dp_r() {return n_points_dp_r ;}
  string _density_profile_r_output_file() {return density_profile_r_output_file ;}
  real_prec _k_min_dp() {return k_min_dp ;}
  real_prec _k_max_dp() {return k_max_dp ;}
  string _scale_dp_k() {return scale_dp_k ;}
  int _n_points_dp_k() {return n_points_dp_k ;}
  string _density_profile_k_output_file(){return density_profile_k_output_file ;}
  bool _compute_output_linear_correlation_function(){return compute_output_linear_correlation_function;}
  bool _compute_output_non_linear_correlation_function(){return compute_output_non_linear_correlation_function;}


};

#endif
