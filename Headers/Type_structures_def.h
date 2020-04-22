#ifndef _STRUCTURES_
#define _STRUCTURES_


# include "def.h"
# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include <vector>
using namespace std;

// *******************************************************************
// *******************************************************************

/**
* @struct<s_CosmoInfo>
* @brief The s_CosmoInfo struct
* @details Auxiliary structure containg cosmological quantities derived from cosmological parameter in class::Cosmology
*/

struct s_CosmoInfo
{
  
  /**
   *@brief Critical energy density of the Universe
  */
  real_prec critical_density;
  /**
    *@brief Density contrast for spherical collapse */
  real_prec density_contrast_top_hat;  
  /**
    *@brief Hubble parameter at input cosmological redshift */
  real_prec Hubble_parameter;    
  /**
    *@brief Comoving distance at input cosmological redshift */
  real_prec comoving_distance;   
  /**
    *@brief Comoving angular diameter distance */
  real_prec comoving_angular_diameter_distance; 
  /**
    *@brief Mean matter density of the Universe at input cosmological redshift */
  real_prec mean_matter_density;  
  /**
    *@brief Age of the Universe at input cosmological redshift */
  real_prec age_universe;       
  /**
    *@brief Comoving sound horizon at input cosmological redshift*/
  real_prec comoving_sound_horizon;  
  /**
    *@brief Growth factor (valid for LCDM) at input cosmological redshift*/
  real_prec growth_factor; 
  /**
    *@brief Auxiliary growth used by 2LPT */
  real_prec D2;     
  /**
    *@brief Growth index at input cosmological redshift*/
  real_prec growth_index;  
  /**
    *@brief Halo dynamical time*/
  real_prec halo_dynamical_time;
  /**
    *@brief Matter energy density parameter at input cosmological redshift*/
  real_prec omega_matter;   
  /**
    *@brief Distance modulos at input cosmological redshift*/
  real_prec Distance_Modulus;  
  /**
    *@brief Cosmological scale factor at input cosmological redshift*/
  real_prec scale_factor;      
};


// *******************************************************************
// *******************************************************************
// *******************************************************************
//Structure containing the cosmological parameters

/**
* @struct<s_CosmologicalParameters>
* @brief The s_CosmoInfo struct
* @details Auxiliary structure containg cosmological parameters
*/

struct s_CosmologicalParameters
{
  real_prec cosmological_redshift;
  /**
    *@brief  Energy density of total matter in units of the the critical density */
  real_prec Om_matter;     
  /**
*@brief  Energy density of cold_dark_matter in units of the the critical density */
  real_prec Om_cdm;      
  /**
*@brief  Energy density of radiation in units of the the critical density */
  real_prec Om_radiation;  
  /**
*@brief  Energy density of baryons in units of the the critical density */
  real_prec Om_baryons;    
  /**
*@brief  Energy density of dark energy in units of the the critical density */
  real_prec Om_vac;        
/**
*@brief  Energy density of curvature in units of the the critical density */
  real_prec Om_k;          
  /**
*@brief  Baryon fraction*/
  real_prec f_baryon;
  /**
*@brief  Hubble parameter in units of h km / s / Mpc */
  real_prec Hubble;        
  /**
*@brief  dimensionless Hubble parameter*/
  real_prec hubble;        
  /**
*@brief  Equation of state of dark energy */
  real_prec w_eos;         
  /**
*@brief  Effective number of relativistic particles*/
  real_prec N_eff;         
  /**
*@brief  RMS of matter fluctuations at R = 8 Mpc/h*/
  real_prec sigma8;        
  /**
*@brief   Amplitude of primordial ower spectrum*/
  real_prec A_s;            
  /**
*@brief  Primordial spectral index*/
  real_prec n_s;           
  /**
*@brief  primordial spectral index */
  real_prec alpha_s;           
/**
 *@brief  defined from patchy*/
  real_prec D1;   
  /**
*@brief  defined from patchy*/
  real_prec D2;   
  /**
*@brief  CMB temperature*/
  real_prec Tcmb;  
  /**
*@brief  Something related to a magnitude, used for vmax stuff*/
  real_prec Mabs;   
  /**
*@brief  Something related to a magnitude, used for vmax stuff*/
  real_prec mlim;   
  /**
*@brief  Comoving scale to compute the rms of mass fluctuation. Usually 8mpc/h */
  real_prec RR;           
 /**
*@brief  Spherical overdensity*/
  real_prec Delta_SO;     
/**
*@brief  Use baryonic wiggles in the linear P(k) */
  bool use_wiggles;        
/**
 *@brief   constant factor in the K-correction */
  real_prec K_index_a;  
/**
*@brief   constant factor in the K-correction */
  real_prec K_index_b;  
/**
*@brief   constant factor in the K-correction */
  real_prec K_index_c;  
/**
 *@brief   constant factor in the K-correction */
  real_prec K_index_d;  
/**
*@brief   constant factor in the K-correction */
  real_prec K_index_e;  
/**
*@brief   constant factor in the K-correction */
  real_prec K_index_f;  
/**
*@brief   constant factor in the K-correction*/
  real_prec K_index_g;  
/**
*@brief   constant factor in the K-correction*/
  real_prec K_index_h;  
/**
*@brief   constant factor in the K-correction*/
  real_prec K_index_i;  
  /**
*@brief   constant factor in the K-correction*/
  real_prec K_index_j;  

/**
*@brief   constant factor in the e-correction */
  real_prec e_index_a;  
/**
*@brief   constant factor in the e-correction */
  real_prec e_index_b;  
/**
*@brief   constant factor in the e-correction */
  real_prec e_index_c;  
/**
*@brief   constant factor in the e-correction */
  real_prec e_index_d;  
/**
*@brief   constant factor in the e-correction */
  real_prec e_index_e;  

  /**
*@brief   constant factor in the e-correction */
  real_prec e_index_zstar;  
/**
 *@brief  Beta factor used to mimic RSD*/
  real_prec beta_rsd;   

  /**
*@brief   constant factor in the k-correction */
  real_prec k_index;  
/**
*@brief   related to the density evolution */
  real_prec d_index; 
  /**
 *@brief  Galaxy bias, used in Cl */
  real_prec GAL_BIAS; 
  real_prec alpha_BIAS;
  
  
  // aca van variables que estan asociadas a integracion en k, o
  // que van servir para poner más parámetros y asi pasar una sola
  // estructura a las rutinas de integracion
  real_prec kmin_int;
  
  real_prec kmax_int;
  string mass_function_fit;
  string halo_mass_bias_fit;
  string density_profile;


  /**
 *@brief  values of integrals used in the HaloFit for intergaration wrt to z of P(k,z) */
  real_prec h4; 
  real_prec h2;
  
  real_prec M_max_mf;
  real_prec M_min_mf;
  real_prec kmin_ps;
  real_prec kmax_ps;

  /**
*@brief  use it for m*/
  real_prec aux_var1; 
/**
*@brief  use it for m*/
  real_prec aux_var2; 
  /**
*@brief  use it for z*/
  real_prec aux_var3; 
  /**
*@brief  use it for k*/
  real_prec aux_var4; 
  

  /**
*@brief  Non linaer mass scale*/
  real_prec Mnl;
  /**
*@brief  Non linear wavenumber*/
  real_prec knl;
  /**
*@brief  Non linear scale*/
  real_prec rnl;
  
  /**
*@brief  Non linear wavenumber from Halo-fit*/
  real_prec knl_hf;
  /**
 *@brief  Non linear scale from Halo-fit*/
  real_prec rnl_hf;
  /**
 *@brief  Non linear wavenumber from Halo-fit*/
  real_prec kstar;
  
  real_prec M_min_effective;
  real_prec M_max_effective;

  // Q-model related parameter
  real_prec A_PS;
  real_prec Q_PS;

  real_prec Amc; //Amplitude of the non linear correctionto P(k) in PT

  // Redshift dependent quantities
  real_prec critical_density;
  real_prec density_contrast_top_hat;
  real_prec mean_matter_density;
  real_prec growth_factor;
  real_prec pk_normalization;


  // Integration with respect to Mass
  int n_points_mass;
  int n_points_mass_integration;


/**
*@brief  Vector used to speed up calculations by means of interpolation. Must come in gsl_real precision if used in interpolations */
  vector<gsl_real> v_mass; 
  /**
*@brief  Container for sigma_mass */
  vector<gsl_real> v_sigma_mass; 
  /**
*@brief  Container for abundance (mass) */
  vector<gsl_real> v_mass_function; 
  /**
 *@brief  Container for bias as a functin of mass */
  vector<gsl_real> v_halo_mass_bias; 
  /**
*@brief  Container for wavenumbers in P(k) */
  vector<gsl_real> v_k_ps;  
/**
*@brief  Container for linear P(k) */
  vector<gsl_real> v_lin_power_spectrum; 
  /**
*@brief  Container for non-linear P(k) */
  vector<gsl_real> v_nl_power_spectrum; 
  /**
*@brief  Container for density profile in Fourier space */
  vector<gsl_real> v_density_profile_k; 
  /**
*@brief  Container for density profile in configuration space */
  vector<real_prec> v_density_profile_r; 
  /**
*@brief  Container for galaxy P(k) 1-halo term satelite-satellite */
  vector<gsl_real> v_galaxy_power_spectrum_1h_ss; 
  /**
*@brief  Container for galaxy P(k) 1-halo term satelite-central */
  vector<gsl_real> v_galaxy_power_spectrum_1h_sc; 
  /**
*@brief  Container for galaxy P(k) 2-halo term */
  vector<gsl_real> v_galaxy_power_spectrum_2h;
  
  bool use_K_correction;
  bool use_e_correction;

  // Vector used to speed up cosmo-calculations by means of interpolation
  // Must come in gsl_real precision
  /**
*@brief  Container for comoving distance r(z) */
  vector<gsl_real>rv;  
  /**
*@brief  Container for growth factor g(z) */
  vector<gsl_real>gv; 
  /**
 *@brief  Container for  redshift z */
  vector<gsl_real>zv; 
  /**
 *@brief  Container for transversal separation distance tr(z) */
  vector<gsl_real>trv; 

  /**
*@brief  Selection for HOD model*/
  int hod_model;  
  /**
*@brief  HOD parameter*/
  real_prec mmin_hod;  
  /**
*@brief  HOD parameter*/
  real_prec alpha_hod; 
  /**
*@brief  HOD parameter*/
  real_prec scatter_hod;
  /**
*@brief  HOD parameter*/
  real_prec muno_hod;
  /**
 *@brief  HOD parameter*/
  real_prec coef_concentration;
  /**
 *@brief  HOD parameter*/
  real_prec coef_concentration_amp;
};



// *******************************************************************
// *******************************************************************
// *******************************************************************
/**
 *@brief
 * @brief The s_hods_par struct
 * @details Structure containing parameers of simple model of HOD
 */
struct s_hods_par{
  int hod_model;
  real_prec mmin;
  real_prec alpha;
  real_prec scatter;
  real_prec muno;
};

// *******************************************************************
// *******************************************************************
// *******************************************************************
/**
 *@brief
* @struct<s_astrophysical_parameters>
* @brief The s_astrophysical_parameters struct
* @details Structure containing parameters of gas-mass scaling relation in clusters
 */
struct s_astrophysical_parameters{
  real_prec A_gas;
  real_prec B_gas;
  real_prec mstar;
  real_prec sigma_red;
  real_prec sigma_ln;
  real_prec missing_flux;
};


// *******************************************************************
// *******************************************************************
// *******************************************************************
/**
* @brief
* @struct<matrices>
* @brief The matrices struct
* @details Structure with 2d cointainers with c onariance matrices, used in FB analysis
*/

struct matrices{
/**
 *@brief  Measure
*/
 vector<real_prec>Cmed;

/**
*@brief  Measure fluctuation
*/
  vector<real_prec>Dmed;

/**
*@brief  Matrix R
*/

  vector<vector<real_prec> >R;
  /**
 *@brief  Matrix V
*/
  vector<vector<real_prec> >V;

  /**
*@brief  Noise matrix
*/
  vector<vector<real_prec> >N;

  /**
 *@brief  Inverse of covariance matrix
*/
  vector<vector<real_prec> >iCov;

  /**
*@brief  Covariance matrix
*/
  vector<vector<real_prec> >nCov;

  /**
  *@brief Something else
  */

  vector<vector<real_prec> >Step;
  real_prec det_matrix;
};



// *******************************************************************
// *******************************************************************
// *******************************************************************
/**
*@brief
* @struct<params_clth>
 * @brief The params_clth struct
 * @details  Auxiliary strucure for the class:ClFunctions computation of theoretical angular power spectrum
 */
struct params_clth{
/**
*@brief   comoving disance
*/
    real_prec r;
/**
 *@brief  Wave number
*/
    real_prec k;
 /**
*@brief  Angular multipole
*/
    int l;
  string wtype;
  /**
 *@brief  */
  real_prec zmax_bin;
  /**
 *@brief  */
  real_prec zmin_bin;
  /**
 *@brief  */
  real_prec rmax_bin;
  /**
 *@brief  */
  real_prec rmin_bin;
  /**
*@brief  */
  real_prec k_min_integration;
  /**
*@brief  */
  real_prec k_max_integration;
  /**
*@brief  */
  real_prec sigma_p_errors; 
/**
 *@brief  */
  string pdf_zerrors; 
  /**
*@brief  */
  real_prec zaux; 
  /**
*@brief  */
  vector<real_prec>pk; 
  /**
 *@brief  */
  vector<real_prec>kv; 
  /**
*@brief  */
  vector<real_prec>Fkernel; 
  /**
*@brief  */
   vector<real_prec>rv;/** *@brief  */
  vector<real_prec>zv; 
  /**
 *@brief  */
  vector<real_prec>gv; 
  /**
*@brief  */
  vector<real_prec>bias_zv; 
  /**
 *@brief  */
  vector<real_prec>gfv; 
  /**
 *@brief  */
  vector<real_prec>Hv; 
  /**
 *@brief  */
  vector<real_prec>dn_photo; 
  /**
 *@brief  */
  vector<real_prec>dn_spect; 
  /**
 *@brief  */
  vector<real_prec>dz; 
  /**
*@brief  */
  vector<real_prec>klnv; 
  /**
 *@brief  */
  vector<real_prec>sigma_mass; 
  /**
*@brief  */
  vector<real_prec>M_nl; 
  /**
*@brief  */
  vector<real_prec>HF_i4; 
  /**
 *@brief  */
  vector<real_prec>HF_i2; 
  /**
*@brief  */
  vector< vector<gsl_real> > sBessel; 
  /**
 *@brief  */
  vector<gsl_real> sxBessel; 
  /**
 *@brief  */
  string redshift; 
  /**
*@brief  */
  real_prec pk_normalization; 
  /**
 *@brief  */
  real_prec n_s;  
  /**
 *@brief  */
  bool use_non_linear_pk; 
  /**
 *@brief  */
  string type_of_nl_power; 
};


// *******************************************************************
// *******************************************************************
// *******************************************************************
/**
*@brief
* @struct<s_aux>
  @brief Template structure s_aux
  @details This structure is created in order to avoid missconfusing of information among threads when using OMP
*/
template<typename T>
struct s_aux{
/**
*@brief  */
  s_CosmologicalParameters *scp_a; /** *@brief  */
/**
*@brief  */
  params_clth *s_clth; 
  /**
 *@brief  */
  real_prec raux; 
  /**
*@brief  */
  int laux; 
  /**
*@brief  */
  real_prec kaux; 
  /**
*@brief  */
  real_prec zaux; 
/**
 *@brief  */
  vector<gsl_real>v_aux; 
  /**
 *@brief  */
  vector<gsl_real>XX; 
  /**
*@brief  */
  vector<gsl_real>WW; 
  /**
 *@brief  */
  vector<gsl_real>XX_mu; 
  /**
 *@brief  */
  vector<gsl_real>WW_mu; 
  /**
 *@brief  */
  vector<gsl_real>XX_z; 
  /**
*@brief  */
  vector<gsl_real>WW_z; 
  T Ps; 
};


// *******************************************************************
// *******************************************************************
// *******************************************************************
/**
 *@brief
* @struct<A1>
* @brief The A1 struct
* @details Structure used in the passage from Cl to Ps to spped up calculations
*/
struct A1{
  /**
*@brief  */
  s_CosmologicalParameters *s_cp; 
  /**
 *@brief  */
  vector<gsl_real>MASS; 
  /**
 *@brief  */
  vector<gsl_real>MASS_FUNCTION; 
  /**
*@brief  */
  vector<gsl_real>MASS_BIAS; 
  /**
*@brief  */
  real_prec aux_k; 
  /**
*@brief  */
  real_prec aux_z; 
  /**
*@brief  */
  real_prec aux_m; 
};


// *******************************************************************
// *******************************************************************
// *******************************************************************


/**
*@brief
* @struct<s_Halo>
 * @brief The s_Halo struct
 * @details Structure to allocate properties of catalogs
 * @details The full catalog can be contained in a container of this type
  *@code

   vector<s_Halo> halo(NOBJS);
   for(int i=0;i<NOBJECTS;++i)
     halo[i].coord1 = x; // x-coordinate of a catalog with NOBJECS
 *@endcode
 */
struct s_Halo
{
  /**
*@brief  Coordinate 1 (X,r or z)   */
  real_prec coord1; 
  /**
 *@brief   Coordinate 2 (Y,phi or phi) */
  real_prec coord2;
  /**
*@brief   Coordinate 1 (Z or theta, )  */
  real_prec coord3;
  /**
*@brief  Velocity component in  coord1  */
  real_prec vel1;
  /**
 *@brief   Velocity component in  coord2 */
  real_prec vel2;
  /**
 *@brief   Velocity component in  coord3 */
  real_prec vel3;
  /**
 *@brief  Tracer Mass  */
  real_prec mass;
  /**
 *@brief  Tracer Vmax   */
  real_prec vmax;
  /**
 *@brief   ID (in the mesh) where this tracer is located given the nomial resolution  */
  ULONG GridID;
  /**
 *@brief  ID (in the mesh) where this tracer is located given the resolution with l1  */
  ULONG GridID_l1;
  /**
 *@brief   ID (in the mesh) where this tracer is located given the resolution with l2  */
  ULONG GridID_l2;
  /**
*@brief   ID (in the mesh) where this tracer is located given the resolution with l3  */
  ULONG GridID_l3;
  /**
 *@brief   ID (in the mesh) where this tracer is located given the resolution with l4  */
  ULONG GridID_l4;
  /**
 *@brief   Number of substructures (if the tracer is parent halo) */
  int number_sub_structures;
  /**
 *@brief  Identify whether this is a random object (0) or a real tracer (1)  */
  int identity;
  /**
*@brief  Identify whether this has been observed (or used in different processes inside the code)  */
  bool observed;
  /**
 *@brief   Weight (if any) */
  real_prec weight1;
  /**
 *@brief   Weight (if any) */
  real_prec weight2;
  /**
*@brief   Weight (if any) */
  real_prec weight3;
  /**
*@brief   Weight (if any) */
  real_prec weight4;
  /**
 *@brief   Mean density (uf any) evaluated at the position of the tracer */
  real_prec mean_density;
  /**
*@brief  PropID */
  ULONG PropID;
};

// *******************************************************************
// *******************************************************************
// *******************************************************************
// *******************************************************************
// *******************************************************************
/**
*@brief
* @struct<experiments>
* @brief The experiments struct
* @details Containers used in the MCMC analysis of Cl. Used in class::HGAP
*/
struct experiments{
  vector<vector<real_prec> > acc_par1;
  vector<real_prec>  weight_par1;
  vector<vector<real_prec> > acc_par2;
  vector<real_prec>  weight_par2;
  vector<vector<real_prec> > acc_par3;
  vector<real_prec>  weight_par3;
  vector<vector<real_prec> > acc_par4;
  vector<real_prec>  weight_par4;
  vector<vector<real_prec> > acc_par5;
  vector<real_prec>  weight_par5;
  vector<vector<real_prec> > acc_par6;
  vector<real_prec>  weight_par6;
};




// *******************************************************************
// *******************************************************************
// *******************************************************************
//Structure containing the properties of the catalogs
/**
*@brief
* @struct<s_data_structure>
 * @brief The s_data_structure struct
 * #details Auxiliary structure used in the measrement of power spectrum class::PowerSpectrum
 */
struct s_data_structure{
  /**
 *@brief  Structure with tracer properties */
  vector<s_Halo> properties;
  /**
*@brief  Number of columns in input catalog*/
  int n_columns;
  /**
 *@brief  Coordinate system use din the input tracer cat*/
  int system_of_coordinates; 
  /**
 *@brief  This specifies whether this is the "random" or the "data" catalogue */
  string catalog;         
  /**
 *@brief  Mean density of catalog*/
  real_prec mean_density;
  /**
 *@brief  Ask whether the mean density is tabulated in the input catalog*/
  bool nbar_tabulated;
  /**
 *@brief  Ask whether the redshift distribution has to be computed from the input catalog */
  bool compute_dndz;
  /**
 *@brief  Container for redshift used if nbar hast to be computed here*/
  vector<gsl_real> zz_v;
  /**
*@brief  Container for redshift distribution*/
  vector<gsl_real> dndz_v;
  /**
*@brief  2D-Container for redshift distribution in different HEALPIX pixels*/
  vector< vector<gsl_real> > dndz_matrix;
};

// *******************************************************************
// *******************************************************************

/**
*@brief
* @struct<s_parameters_estimator>
 * @brief The s_parameters_estimator struct
   @details Auxiliary structure used in the measurement of power spectrum, fed in class::PowerSpectrum, used in class::FftwFunctions
 */
struct s_parameters_estimator{
/**
*@brief  Use random catalog
*/
  bool use_random_catalog;

  /**
*@brief  Number of selected real objects from which the P(k) will be measured
*/
  int number_of_objects;

  /**
 *@brief  number of selected random objects from which the P(k) will be measured
*/
  int number_of_randoms;

  /**
 *@brief  weighted number of selected real objects from which the P(k) will be measured
*/
   real_prec w_number_of_objects;

  /**
*@brief weighted number of selected random objects from which the P(k) will be measured
*/
  real_prec w_number_of_randoms;  

  /**
*@brief   w_number_of_objects divided by w_number_of_randoms
*/
  real_prec alpha;             

  /**
 *@brief  Sum of the squared of the weights in the real cataloge
*/
  real_prec S_g;

  /**
*@brief  Sum of the squared of the weights in the random cataloge
*/
  real_prec S_r;               

  /**
 *@brief  Noramlization of power spectrum
 */
  real_prec normalization;

  /**
*@brief  Normalization of window function
*/
  real_prec normalization_window;

  /**
 *@brief  Shot noise
*/
  real_prec shot_noise;         

  /**
 *@brief  Shot noise for window
*/
  real_prec shot_noise_window;  
};


// *******************************************************************
// *******************************************************************
/**
*@brief
* @struct<s_parameters_bis_estimator>
* @brief The s_parameters_bis_estimator struct
*/
struct s_parameters_bis_estimator{
  bool use_random_catalog;
  int number_of_objects;  //number of selected real objects from which the P(k) will be measured
  int number_of_randoms;  //number of selected random objects from which the P(k) will be measured
  real_prec w_number_of_objects;  //weighted number of selected real objects from which the P(k) will be measured
  real_prec w_number_of_randoms;  //weighted number of selected random objects from which the P(k) will be measured
  real_prec alpha;                //  w_number_of_objects divided by w_number_of_randoms
  real_prec S_g1;                  // Sum of the squared of the weights in the real cataloge
  real_prec S_g2;                  // Sum of the squared of the weights in the real cataloge
  real_prec S_r1;                  // Sum of the squared of the weights in the random cataloge
  real_prec S_r2;                  // Sum of the squared of the weights in the random cataloge
  real_prec normalization;  //Sum of the product of mean density times weights squared, used for the normalization
  real_prec shot_noise1;           //Shot noise of the power spectrum
  real_prec shot_noise2;           //Shot noise of the power spectrum

  real_prec normal_power;
  real_prec shot_noise_power;
};



// *******************************************************************
// *******************************************************************
/**
*@brief
* @struct<s_parameters_box>
* @brief The s_parameters_box struct
*/
struct s_parameters_box{
  bool use_random_catalog;
  string mas;       //Mass Assignment scheme
  string ave;        //Type of average
  real_prec k_bin_step;
  bool use_MAS_correction; 
  bool FKP_weight;
  bool FKP_error_bars;
  bool FKP_error_bars_exact;
  bool use_SN_correction;
  real_prec Pest;
  bool nbar_tabulated;
  bool compute_dndz;
  bool constant_depth;
  vector<gsl_real>zz;
  vector<gsl_real>rc;
  int n_dndz;
  int new_n_dndz;
  real_prec redshift_min_sample;
  real_prec redshift_max_sample;
  real_prec area_survey;
  int sys_of_coord_r;
  int i_coord1_r;
  int i_coord2_r;
  int i_coord3_r;
  int i_weight1_r;
  int i_weight2_r;
  int i_weight3_r;
  int i_weight4_r;
  bool use_weight1_r;
  bool use_weight2_r;
  bool use_weight3_r;
  bool use_weight4_r;
  int i_mean_density_r;
  string angles_units_r;
  int i_coord1_g;
  int i_coord2_g;
  int i_coord3_g;
  int i_mass_g;
  int i_weight1_g;
  int i_weight2_g;
  int i_weight3_g;
  int i_weight4_g;
  bool use_weight1_g;
  bool use_weight2_g;
  bool use_weight3_g;
  bool use_weight4_g;
  int i_mean_density_g;
  string angles_units_g;
  real_prec area_pixel;
  long npixels;
  long nside;
  string file_dndz;
  bool new_los;
  real_prec kmin_bk;
  real_prec kmax_bk;
  bool use_fundamental_mode_as_kmin_bk;
  bool measure_cross;
};

// *******************************************************************
// *******************************************************************
/**
*@brief
* @struct<s_parameters_box_yam>
* @brief The s_parameters_box_yam struct
*/
struct s_parameters_box_yam{
  int n1;           //Number of grid cells in the X-direction
  int n2;           //Number of grid cells in the Y-direction
  int n3;           //Number of grid cells in the Z-direction
  real_prec lx;        //LENGHT OF X SIDE OF FOURIER BOX, IN Mpc/h
  real_prec ly;        //LENGHT OF Y SIDE OF FOURIER BOX, IN Mpc/h
  real_prec lz;        //LENGHT OF Z SIDE OF FOURIER BOX, IN Mpc/h
  real_prec xmin;      //Minimim value of the X cordinate in Mpc/h
  real_prec ymin;      //Minimim value of the Y cordinate in Mpc/h
  real_prec zmin;      //Minimim value of the Z cordinate in Mpc/h
  bool use_random_catalog;
  string mas;       //Mass Assignment scheme
  int nnp_data;     //Number of modes to write output of P(k)
  int nnp_window;   //Number of modes to write output of W(k)
  real_prec Deltal;    //Log nin size in log-10
  real_prec DeltaK_data;//Linear bin size for P(k)
  real_prec DeltaK_window; //Linear bin size for W(k)
  string ave;        //Type of average
  real_prec kmin;       //Minimum wavenumber for log average
  real_prec kmax;
  bool use_MAS_correction; 
  bool FKP_weight;
  bool FKP_error_bars;
  bool FKP_error_bars_exact;
  bool use_SN_correction;
  int cartesian_los;
  real_prec ra_los;
  real_prec dec_los;
  real_prec Pest;
  int Nmu;
  bool nbar_tabulated;
  bool compute_dndz;
  bool constant_depth;
  vector<gsl_real>zz;
  vector<gsl_real>rc;
  int n_dndz;
  int new_n_dndz;
  real_prec redshift_min_sample;
  real_prec redshift_max_sample;
  real_prec area_survey;
  int sys_of_coord_r;
  int i_coord1_r;
  int i_coord2_r;
  int i_coord3_r;
  int i_weight1_r;
  int i_weight2_r;
  int i_weight3_r;
  int i_weight4_r;
  bool use_weight1_r;
  bool use_weight2_r;
  bool use_weight3_r;
  bool use_weight4_r;
  int i_mean_density_r;
  string angles_units_r;
  int i_coord1_g;
  int i_coord2_g;
  int i_coord3_g;
  int i_weight1_g;
  int i_weight2_g;
  int i_weight3_g;
  int i_weight4_g;
  bool use_weight1_g;
  bool use_weight2_g;
  bool use_weight3_g;
  bool use_weight4_g;
  int i_mean_density_g;
  string angles_units_g;
  real_prec area_pixel;
  long npixels;
  long nside;
  string file_dndz;
  bool new_los;
};



// *******************************************************************
/**
*@brief
* @struct<pic_storage>
 * @brief The pic_storage struct
 */
struct pic_storage {
  real_prec weight;
  int idx;
  //  int pad;            // padding would lower a bit the cache misses
};


// *******************************************************************
// In the file Coordinates.h we defiend the struct s_galaxy_operationsF
// whis is the same acept for the last element, which is vector<>
// used in the FFT measurements (as written by Luca)
// These should be merged once I transform all readsintgs to the
// one giving a vector.

/**
 *@brief
* @struct<s_galaxy_operation>
* @brief The s_galaxy_operations struct
*/

struct s_galaxy_operations{
  int sys_coord;
  int i_coord1;
  int i_coord2;
  int i_coord3;
  string angles_units;
  real_prec redshift_max;
  real_prec redshift_min;
  int nz;
  vector< vector<real_prec> > properties;
};


// *******************************************************************
/**
*@brief
* @struct<params>
 * @brief The params struct
  * @details IS THIS TO BE DEPRECATED???
 */
struct params
{
  int NX;
  int NY;
  int Nlambdath;
  int ndel_data;
  string Output_directory;
  string Input_Directory_Y;
  string Name_Catalog_Y;
  string Input_Directory_X;
  string Input_Directory_X_REF;
  string XNAME;
  string Name_Catalog_X;
  string Name_Catalog_X_REF_PDF;
  string Name_Property_X;
  string new_Name_Property_X;
  string YNAME;
  string Name_Property_Y;
  string new_Name_Property_Y;
  string Type_of_File_X;
  string Type_of_File_Y;
  int iMAS_X;
  int iMAS_X_REF_PDF;
  int iMAS_Y;
  real_prec delta_Y_max;
  real_prec delta_Y_min;
  real_prec delta_X_max;
  real_prec delta_X_min;
  real_prec ldelta_Y_max;
  real_prec ldelta_Y_min;
  real_prec ldelta_X_max;
  real_prec ldelta_X_min;
  string Quantity;
  int NCW_comb;
  int NMASSbins;
  real_prec redshift;
  real_prec smscale;
  int realization;
  int Nft;
  bool get_SKNOTS;
  bool Comp_conditional_PDF;
  bool Comp_joint_PDF;
  bool write_files_for_histograms;
  bool Redefine_limits;
  bool Convert_Density_to_Delta_X;
  bool Convert_Density_to_Delta_Y;
  real_prec lambdath;
  bool Write_Scatter_Plot;
  bool Write_PDF_number_counts;  
  string Scale_Y;
  string Scale_X;
  int n_sknot_massbin;
  bool get_mock;
  real_prec Lbox;
  bool CWC_in;
  int N_iterations;
  int N_iterations_dm;
  bool Apply_Rankordering;
  vector<int> cwt_used;
  int n_cwt;
};


// *******************************************************************
/**
*@brief
* @struct<s_params_box_mas>
* @brief Auxiliary structure used to pass box arguments to different functions
 */
struct s_params_box_mas{
  int masskernel;
  int Nft;
  ULONG NGRID;
  real_prec Lbox;
  real_prec d1;
  real_prec d2;
  real_prec d3;
  real_prec min1;
  real_prec min2;
  real_prec min3;
};


// *******************************************************************

/**
 *@brief
* @struct<s_Galaxy_struct>
 * @brief The s_Galaxy struct
 */
struct s_Galaxy
{
  string type_of_object;
  ULONG NOBJS;
  vector<real_prec> xgal;
  vector<real_prec> ygal;
  vector<real_prec> zgal;
  vector<real_prec> Pxgal; //P meant for vector properties such as velocity or momentum
  vector<real_prec> Pygal;
  vector<real_prec> Pzgal;
  vector<real_prec> VPgal;
  vector<real_prec> Mass;  // other scalar property such as the mass
  vector<real_prec> property;  // other scalar
  vector<real_prec> weight;

};

// *******************************************************************

/**
 *@brief
* @struct<s_minimums>
 * @brief The s_minimums struct
 * @details Auxiliary structure to allocate minimum values of properties used to characterize the bias in class::Bam

 */
struct s_minimums{
   real_prec prop0;  // Tracer property, counts
   real_prec prop0_mass;  // Tracer property, mass
   real_prec prop0_sf;  // Tracer property, satellite_fraction
   real_prec prop1;   // ð
   real_prec prop2;   // CWC
   real_prec prop3;   // Mk
   real_prec prop4;  // Inv1 or ð²
   real_prec prop5;  //InvII ot ð³
   real_prec prop6;  //Tidal Anis or s²
   real_prec prop7;  //InvI Shear or Nabla²ð
   real_prec prop8;  //InvII shear or s²ð
   real_prec prop9;  //InvIII Shear or s³
   real_prec prop10;
};

// *******************************************************************

/**
 *@brief
* @struct<s_maximums>
 * @brief The s_maximums struct
 * @details Auxiliary structure to allocate  maximum values of properties used to characterize the bias in class::Bam
 */
struct s_maximums{
  real_prec prop0;
  real_prec prop0_mass;  // Tracer property, mass
  real_prec prop0_sf;  // Tracer property, satellite_fraction
  real_prec prop1;
  real_prec prop2;
  real_prec prop3;
  real_prec prop4;
  real_prec prop5;
  real_prec prop6;
  real_prec prop7;
  real_prec prop8;
  real_prec prop9;
  real_prec prop10;
};

// *******************************************************************

/**
*@brief
* @struct<s_Deltas>
 * @details Auxiliary structure to allocate  bin sizes in the different properties used to characterize the bias in class::Bam
 */
struct s_Deltas{
 /**
 *@brief   Bin size for property0 */
  real_prec prop0;
    /**
 *@brief   Bin size for property0 using tracer mass */
  real_prec prop0_mass;  // Tracer property, mass
  /**
*@brief   Bin size for property0 satellite fraction */
  real_prec prop0_sf;  // Tracer property, satellite_fraction
  /**
 *@brief   Bin size for property1 */
  real_prec prop1;
  /**
 *@brief   Bin size for property2 */
  real_prec prop2;
  /**
*@brief   Bin size for property3 */
  real_prec prop3;
  /**
*@brief   Bin size for property4 */
  real_prec prop4;
  /**
 *@brief   Bin size for property5 */
  real_prec prop5;
  /**
*@brief   Bin size for property6 */
  real_prec prop6;
  /**
*@brief   Bin size for property7 */
  real_prec prop7;
  /**
*@brief   Bin size for property8 */
  real_prec prop8;
  /**
*@brief   Bin size for property9 */
  real_prec prop9;
  /**
 *@brief   Bin size for property10 */
  real_prec prop10;
};


// *******************************************************************

/**
 *@brief
 * @struct<s_mass_members>
 * @details Strucure to allocate set of masses, coordinates and their ID (in the mesh) for tracers within a given bin of the multidimendional Theta (DM) properties
 * @details Used in the assignation of halo properties.
 * @details Example: Private variable in class::Bam
 * @code
    vector<s_mass_members> dm_properties_bins;
 *@endcode
 */
struct s_mass_members{
    /**
*@brief  Masses (or any other property) in that particular theta_bin*/
    vector<real_prec> masses_bin_properties;
    /**
 *@brief  Flag to pinpoint properties used*/
    vector<bool> used_mass;
    /**
*@brief  x-coordinates of the objects in that theta-bin */
    vector<real_prec> x_coord_bin_properties;
    /**
*@brief  y-coordinates of the objects in that theta-bin */
    vector<real_prec> y_coord_bin_properties;
    /**
*@brief  z-coordinates of the objects in that theta-bin */
    vector<real_prec> z_coord_bin_properties;
    /**
*@brief  Grid ID of the objects in that theta-bin */
    vector<ULONG> GridID_bin_properties;
    /**
 *@brief  Object ID of the objects in that theta-bin */
    vector<ULONG> GalID_bin_properties;

};

// *******************************************************************

/**
*@brief
 * @struct<s_nearest_cells struct>.
   @details This is used in NumericalMethods::get_neighbour_cells()
   @details Structure containing the closest cell neirbours
   @code
   vector<s_nearest_cells>nearest_cells_to_cell(NGRID);
   @endcode
*/

struct s_nearest_cells{
 /**
 *@brief  Container with dimension Number_of_neighbours*/
  vector<ULONG> close_cell;
  /**
*@brief  Tracks info of boundary conditions in the x-direction */
  vector<int> bc_x;
  /**
 *@brief  Tracks info of boundary conditions in the y-direction */
  vector<int> bc_y;
  /**
 *@brief  Tracks info of boundary conditions in the z-direction */
  vector<int> bc_z;
};

// *******************************************************************

/**
*@brief
 * @struct<s_cell_info struct>
   @details Structure for each cell of the mesh, containing the coordinates, mass, other properties, tracer-index and the position of the cell in the Theta-histgrams.
   Used in class::Catalog
   @code
  vector<s_cell_info> cell_info_tr(NGRID);
  for(ULONG i=0;i<NOBJS ;++i)
    {
      ULONG ID=halo[i].GridID;
      cell_info_tr[ID].posx_p.push_back(halo[i].coord1);  // and so on
    }
 @endcode
 @details The size of the containers of this structure is the number of tracers in each cell
 @code
  vector<int>Ncounts(NGRID,0);
  for(ULONG i=0;i<NGRID ;++i)
     int Ncounts[i] = cell_info_tr[i].posx.size();
 @endcode
 */
struct s_cell_info{
  /**
 *@brief  Container to allocate the coordinates1 of objects in a cell*/
  vector<real_prec> posx_p;
  /**
*@brief  Container to allocate the coordinates2 of objects in a cell*/
  vector<real_prec> posy_p;
  /**
*@brief  Container to allocate the coordinates3 of objects in a cell*/
  vector<real_prec> posz_p;
  /**
 *@brief  Container to allocate the mass of objects in a cell*/
  vector<real_prec> mass;
  /**
*@brief  Container to allocate an extra-property of objects in a cell*/
  vector<real_prec> property;
  /**
*@brief  Container to allocate the id of each object in a cell*/
  vector<ULONG> gal_index;
  /**
*@brief  Position in the multidimensional DM property Theta of the current cell*/
  ULONG Theta_bin;
};

// *******************************************************************

/**
*@brief
 * @struct<s_dist_in_dmbins struct>
 * @details Strucure meant to allocate the list of masses in pairs that, within one bin of the dm properties,
   are separated by a distnace between m_min and MAXIMUM_DISTANCE_EXCLUSION
 */
struct s_dist_in_dmbins{
  /**
*@brief  Mass (or other prop) of one element of pair */
  vector<real_prec> M1;
 /**
*@brief  Mass (or other prop) of second element of pair */
  vector<real_prec> M2;
  vector<ULONG> mass_to_swap;
};
  
#endif
  
