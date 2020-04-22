/**
 *  @file Params.h
 *
 *  @brief The class Parameters. Reads the parameter 
 *
 */

#ifndef __PARAMS__
#define __PARAMS__

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include "def.h"
#include "NumericalMethods.h"
#include "Type_structures_def.h"

using namespace std;



class Params 
{
  
 private :
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string statistics;
  

  real_prec Redshift_initial;
  
  //////////////////////////////////////////////////////////
  /**
   * @brief input directory where the tracers, dm and random catalogues are stored
   **/
  string dir_input;

  //////////////////////////////////////////////////////////
  /**
   * @brief input directory where the tracers, dm and random catalogues are stored
   **/
  string ic_WN_dir;

  string ic_file;

  bool use_ic_file;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  unsigned long ngal_delta;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  string delta_grid_file;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/

  string delta_grid_file2;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/

  string delta_grid_file3;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  string delta_grid_file4;
  //////////////////////////////////////////////////////////
  /**
   * @brief ask if cross is desidered
  **/
  int measure_cross_from_1;
  //////////////////////////////////////////////////////////
  /**
   * @brief ask if cross is desidered
  **/

  int measure_cross_from_2;


    /**
   * @brief ask if cross is desidered
  **/

  string Name_redshift_mask;

  string Name_binary_mask;

  
  //////////////////////////////////////////////////////////
  /**
   * @brief Type if input file, options
   * are "cat" meaning catalog, and "grid_delta" meaning
   * that the input is alread the delta in the grid
   **/
  string file_random;

  //////////////////////////////////////////////////////////
  /** 
   *  @note coordinate system of the object catalogue
   *
   *  0 &rarr; positions are given in Cartesian coordinates (X, Y, Z)
   *
   *  1 &rarr; positions are given in equatorial coordinates (RA, Dec, r)
   *
   *  2 &rarr; positions are given in pseudo-equatorial coordinates (RA, Dec, z)
   *
   *  3 &rarr; positions are given in pseudo-equatorial coordinates (RA, Dec, z)
   *  
   *  @note with option 2, the code uses the set of cosmological
   *  parameters to transform redshift z to comoving distance; with
   *  option 3, the redshift is directly used as radial coordinate
   */
  int sys_of_coord_g;

  /**
   * @brief  column with the infoamtio0n of the numbr of sub_structures of the tracer
  **/
  int i_sf_g;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  int i_mass_dm;

  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/

  int sys_of_coord_dm;
  //////////////////////////////////////////////////////////

  /// the column where the first object coordinate (according to the system of coordinates of the catalog) is written
  /**
   * @brief 
  **/
  int i_coord1_dm;

  //////////////////////////////////////////////////////////
  /// the column where the second object coordinate (according to the system of coordinates of the catalog) is written
  /**
   * @brief 
   **/
  int i_coord2_dm;
  

  bool weight_vel_with_mass;


  real_prec MASS_units;
  real_prec LOGMASSmin;
  real_prec LOGMASSmax;  

  real_prec VMAXmin;
  real_prec VMAXmax;

  //////////////////////////////////////////////////////////
  /// the column where the third object coordinate (according to the system of coordinates of the catalog) is written
  /**
   * @brief 
  **/
  int i_coord3_dm;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  int i_v1_dm;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  int i_v2_dm;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  int i_v3_dm;
  
  //////////////////////////////////////////////////////////
  /// the column where the first object weight is written
  /**
   * @brief 
   **/
  int i_weight1_g;
  
  //////////////////////////////////////////////////////////
  /// the column where the second object weight is written
  /**
   * @brief 
   **/
  int i_weight2_g;

  //////////////////////////////////////////////////////////
  /// the column where the third object weight is written
  int i_weight3_g;

  //////////////////////////////////////////////////////////
  /// the column where the fourth object weight is written
  int i_weight4_g;

  //////////////////////////////////////////////////////////
  /// the column where the object mean number density is written
  int i_mean_density_g;

  //////////////////////////////////////////////////////////
  /// the column where the first object property is written
  int i_property1_g;

  //////////////////////////////////////////////////////////
  /// the column where the second object property is written
  int i_property2_g;

  //////////////////////////////////////////////////////////
  /// the column where the third object property is written
  int i_property3_g;

  //////////////////////////////////////////////////////////
  /// the column where the fourth object property is written
  int i_property4_g;

  //////////////////////////////////////////////////////////
  /// the column where the fifth object property is written
  int i_property5_g;

  //////////////////////////////////////////////////////////
  /// the column where the sixth object property is written
  int i_property6_g;

  //////////////////////////////////////////////////////////
  /// true &rarr; use the first weight; false &rarr; do not use the first weight
  bool use_weight1_g;

  //////////////////////////////////////////////////////////
  /// true &rarr; use the second weight; false &rarr; do not use the second weight
  bool use_weight2_g;

  //////////////////////////////////////////////////////////
  /// true &rarr; use the third weight; false &rarr; do not use the third weight
  bool use_weight3_g;

  //////////////////////////////////////////////////////////
  /// true &rarr; use the fourth weight; false &rarr; do not use the fourth weight
  bool use_weight4_g;
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief units of angles in the object catalogue: D &rarr; degrees; R &rarr; radians
  */
  string angles_units_g;

  //////////////////////////////////////////////////////////
  /**
   *  @brief units of angles in the object catalogue: D &rarr; degrees; R &rarr; radians
  */
  string angles_units_r;

  
  //////////////////////////////////////////////////////////
  /**
   *  @brief coordinate system of the object catalogue
   *  0 &rarr; positions are given in Cartesian coordinates (X, Y, Z)
   *  1 &rarr; positions are given in equatorial coordinates (RA, Dec, r)
   *  2 &rarr; positions are given in pseudo-equatorial coordinates (RA, Dec, z)
   *  3 &rarr; positions are given in pseudo-equatorial coordinates (RA, Dec, z)
   *  @note with option 2, the code uses the set of cosmological
   *  parameters to transform redshift z to comoving distance; with
   *  option 3, the redshift is directly used as radial coordinate
   */
  int sys_of_coord_r;

  //////////////////////////////////////////////////////////
  /// the column where the first random coordinate (according to the system of coordinates of the catalog) is written
  int i_coord1_r;

  //////////////////////////////////////////////////////////
  /// the column where the second random coordinate (according to the system of coordinates of the catalog) is written
  int i_coord2_r;

  //////////////////////////////////////////////////////////
  /// the column where the third random coordinate (according to the system of coordinates of the catalog) is written
  int i_coord3_r;

  //////////////////////////////////////////////////////////
  /// the column where the first random weight is written
  int i_weight1_r;

  //////////////////////////////////////////////////////////
  /// the column where the second random weight is written
  int i_weight2_r;

  //////////////////////////////////////////////////////////
  /// the column where the third random weight is written
  int i_weight3_r;

  //////////////////////////////////////////////////////////
  /// the column where the fourth random weight is written
  int i_weight4_r;

  //////////////////////////////////////////////////////////
  /// the column where the random mean number density is written
  int i_mean_density_r;

  //////////////////////////////////////////////////////////
  /// the column where the first random property is written
  int i_property1_r;

  //////////////////////////////////////////////////////////
  /// the column where the second random property is written
  int i_property2_r;

  //////////////////////////////////////////////////////////
  /// the column where the third random property is written
  int i_property3_r;

  //////////////////////////////////////////////////////////
  /// the column where the fourth random property is written
  int i_property4_r;

  //////////////////////////////////////////////////////////
  /// the column where the fifth random property is written
  int i_property5_r;

  //////////////////////////////////////////////////////////
  /// the column where the sixth random property is written
  int i_property6_r;

  //////////////////////////////////////////////////////////
  /// true &rarr; use the first weight; false &rarr; do not use the first weight
  bool use_weight1_r;

  //////////////////////////////////////////////////////////
  /// true &rarr; use the second weight; false &rarr; do not use the second weight
  bool use_weight2_r;

  //////////////////////////////////////////////////////////
  /// true &rarr; use the third weight; false &rarr; do not use the third weight
  bool use_weight3_r;

  //////////////////////////////////////////////////////////
  /// true &rarr; use the fourth weight; false &rarr; do not use the fourth weight
  bool use_weight4_r;


  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  /**
   *  @name FKP power spectrum
   */
  ///@{
  //////////////////////////////////////////////////////////  
  /// Number of grid cells /per dimension for the Discrete Fourier Trasnform


  int Nft_low;

  int Nft_HR;
  //////////////////////////////////////////////////////////
  ///  Lenght of the Fourier box in configuration space (in Mpc/h)
  real_prec Lbox;
  real_prec Lbox_low;

  //////////////////////////////////////////////////////////
  /// true &rarr; compute a new Lbox; false &rarr; the code uses Lbox
  bool new_Lbox;
  
  //////////////////////////////////////////////////////////
  /// the power spectrum is computed with the help of a random catalogue (if use_random_catalog = true), or from a simulation with known mean (if use_random_catalog = false) 
  bool use_random_catalog;
  
  
  //////////////////////////////////////////////////////////
  /// the angular power spectrum is computed with the help of a random catalogue (if use_random_catalog = true), or from a simulation with known mean (if use_random_catalog = false) 
  bool use_random_catalog_cl;
  
  


  //////////////////////////////////////////////////////////
  ///Type of binning in Fourier space: linear/log
  string type_of_binning; 

  //////////////////////////////////////////////////////////
  /// Type of linear binning in Fourier space: 0.5 or 1.
  /// for 0.5, the bins are k_i= (i+0.5)*Delta. In this case the center of the
  /// first bin is at the mid point between the zero mode and the fundamental mode
  /// if Delta=delta.
  /// for 1, the bins are k_i= (i+1)*Delta. In this case the center of the
  /// the first bin is at the fundamental mode, if Delta=delta.
  real_prec k_bin_step;


  //////////////////////////////////////////////////////////
  /// Number of log-spaced bins in Fourier space
  int N_log_bins; 

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int ndel_data;

 
  //////////////////////////////////////////////////////////
  /// Ratio between the shell-width and the fundamental mode for window
  int ndel_window; 

  //////////////////////////////////////////////////////////
  /// Number of mu-bins for P(k,mu)
  int N_mu_bins; 
 
  //////////////////////////////////////////////////////////
  ///Use FKP weights (yes/no)
  bool FKP_weight; 

  //////////////////////////////////////////////////////////
  /// Estimated power for FKP weights
  real_prec Pest;

  
  //////////////////////////////////////////////////////////
  /// Compute FKP error bars? (yes/no)
  bool FKP_error_bars; 
  
  //////////////////////////////////////////////////////////
  /// Compute error bars following FKP exact formula(yes/no)
  /// If this is no, and the previous is yes
  /// the code uses the Veff approximation for the variance.
  bool FKP_error_bars_exact; 
  
  //////////////////////////////////////////////////////////
  /// Is the mean number density tabulated?
  bool nbar_tabulated; 
  
  //////////////////////////////////////////////////////////
  /// Has the sample a constant depth?
  bool constant_depth; 
  
  //////////////////////////////////////////////////////////
  ///Number of redshift bins to measure dNdz
  int N_z_bins; 
  
  //////////////////////////////////////////////////////////
  /// Minimum redshift of the sample. Used when dNdz ahs to be measured 
  real_prec redshift_min_sample;

  //////////////////////////////////////////////////////////
  /// Maximum redshift of the sample. Used when dNdz ahs to be measured 
  real_prec redshift_max_sample;
  
  //////////////////////////////////////////////////////////
  ///  Number of dNdz bins to measure 
  int N_dndz_bins; 

  //////////////////////////////////////////////////////////
  /// Number of redshift bins withoin which the measuerd dNdz will be smoothed
  int new_N_dndz_bins; 

  //////////////////////////////////////////////////////////
  /// Area of the survey.
  real_prec area_survey;

  //////////////////////////////////////////////////////////
  /// Resolution Healpix for pixelization. Used when no nbar is tabulated and dNdz is to be computed from a non-constant depth sample
  int Healpix_resolution; 

  //////////////////////////////////////////////////////////
  /// output file for the redshift distribution
  string file_dndz;   

  //////////////////////////////////////////////////////////
  /// Redefine a new line of sight
  bool new_los;
  

  //////////////////////////////////////////////////////////
  /// output log file for the FKP power spectrum
  string file_power_log;

  //////////////////////////////////////////////////////////
  /// output file for the power spectrum of thw window function
  string file_window;
  //////////////////////////////////////////////////////////
  /// output file for the 2d power spectrum in cartesian coordinates
  string file_power2d; 

  //////////////////////////////////////////////////////////
  /// output file for the 2d power spectrum in polar coordinates
  string file_power2d_mk;

  //////////////////////////////////////////////////////////
  /// Maximum k value for the direct sum approach to Yamamoto-*Blake
  real_prec kmax_y_ds;


  //////////////////////////////////////////////////////////
  /// Maximum k value for the direct sum approach to Yamamoto-*Blake
  int n_catalogues;

  
  ///@}
  
  

  /**
   *  @name bispectrum
   */
  ///@{
  
  //////////////////////////////////////////////////////////
  /// name of the output file used to store the bispectrum
  string file_bispectrum;

  //////////////////////////////////////////////////////////
  /// These parameters is used to define the shells in k-space
  bool use_fundamental_mode_as_kmin_bk;

  
  /////////////////////////////////////////////////////////////
  /// Minimum k-value for constructing k-bins
  real_prec kmin_bk;

  /////////////////////////////////////////////////////////////
  /// Maximum k-value for constructing k-bins
  real_prec kmax_bk;


  ///@}
  


  /**
   *  @name angular power spectrum
   */
  ///@{



  string input_file_mask;  
  //////////////////////////////////////////////////////////
  /// the column where the pixel of the mask is written
  int i_mask_pixel;

  //////////////////////////////////////////////////////////
  /// the column where the RA of the pixel of the mask is written
  int i_mask_alpha;

  //////////////////////////////////////////////////////////
  /// the column where the DEC of the pixel of the mask is written
  int i_mask_delta;

  //////////////////////////////////////////////////////////
  /// the column where the FLAG of the mask is written
  int i_mask_flag;








  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  /**
   *  @name cosmological parameters
   */
  ///@{

  //////////////////////////////////////////////////////////
  /// &Omega;<SUB>M</SUB>: the density of baryons, cold dark matter and massive neutrinos (in units of the critical density) at z=0
  real_prec om_matter;


  //////////////////////////////////////////////////////////
  /// &Omega;<SUB>M</SUB>: the density of baryons, cold dark matter and massive neutrinos (in units of the critical density) at z=0
  real_prec om_cdm;

  //////////////////////////////////////////////////////////
  /// &Omega;<SUB>rad</SUB>: the radiation density at z=0
  real_prec om_radiation;

  ///////////////////////////////////////////////////y///////
  /// &Omega;<SUB>b</SUB>: the baryon density at z=0
  real_prec om_baryons;

  //////////////////////////////////////////////////////////
  /// &Omega;<SUB>DE</SUB>: the dark energy density at z=0
  real_prec om_vac;

  //////////////////////////////////////////////////////////
  /// &Omega;<SUB>k</SUB>: the density of curvature energy
  real_prec om_k;

  //////////////////////////////////////////////////////////
  /// H<SUB>0</SUB>: the Hubble constant at z=0 [km/sec/Mpc] 
  real_prec Hubble;

  //////////////////////////////////////////////////////////
  /// \e h: the Hubble parameter, H<SUB>0</SUB>/100
  real_prec hubble;

  //////////////////////////////////////////////////////////
  /// n<SUB>spec</SUB>: the primordial spectral index
  real_prec spectral_index;

  //////////////////////////////////////////////////////////
  /// w<SUB>0</SUB>: the parameter of the dark energy equation of state
  real_prec w_eos;

  //////////////////////////////////////////////////////////
  /// N<SUB>eff</SUB>: the effective number (for QED + non-instantaneous decoupling)
  real_prec N_eff;

  //////////////////////////////////////////////////////////
  /// &sigma;<SUB>8</SUB>: the power spectrum normalization
  real_prec sigma8;

  //////////////////////////////////////////////////////////
  /// T<SUB>CMB</SUB>: the present day CMB temperature [K]
  real_prec Tcmb;
  
  real_prec A_s;
  real_prec n_s;
  real_prec alpha_s;
  
  real_prec RR;
  real_prec M_reference;
  bool use_wiggles;
  bool fixed_redshift;
  


  
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  ///@}


 /**
   *  @name BAM parameters
   */
  ///@{



  
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string par_file;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  ULONG NGRID;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NX;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NY;

  int NY_MASS;

  int NY_SAT_FRAC;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int Nlambdath;

    //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Input_Directory_Y;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Name_Catalog_Y;

  string Name_Catalog_Y_HR;


  string Name_Catalog_Y_MWEIGHTED;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Input_Directory_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Input_Directory_X_REF;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string XNAME;
  //////////////////////////////////////////////////////////
  /*
   * @brief INput density field for dark matter
   */
  string Name_Catalog_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief INput density field for dark matter
   */
  string Name_VelFieldx_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief INput density field for dark matter
   */
  string Name_VelFieldy_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief INput density field for dark matter
   */
  string Name_VelFieldz_X;
  //////////////////////////////////////////////////////////

  /*
   * @brief
   */
  string Name_Catalog_X_REF_PDF;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Name_Catalog_X_NEW;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Name_Property_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string new_Name_Property_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string YNAME;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Name_Property_Y;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string new_Name_Property_Y;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string extra_info;

  //////////////////////////////////////////////////////////

  /*
   * @brief
   */
  int iMAS_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int iMAS_X_REF_PDF;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int iMAS_X_NEW;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int iMAS_Y;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  real_prec delta_Y_max;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  real_prec delta_Y_min;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  real_prec delta_X_max;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  real_prec delta_X_min;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  real_prec ldelta_Y_max;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  real_prec ldelta_Y_min;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  real_prec ldelta_X_max;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  real_prec ldelta_X_min;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Quantity;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NMASSbins;


  //////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  real_prec redshift;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  real_prec smscale;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int realization;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Comp_conditional_PDF;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Comp_joint_PDF;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool write_files_for_histograms;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Redefine_limits;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Convert_Density_to_Delta_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Convert_Density_to_Delta_Y;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the TWEB classification
   */
  real_prec lambdath;


  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec lambdath_v;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Write_Scatter_Plot;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Write_PDF_number_counts;  
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Scale_Y;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Scale_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int n_sknot_massbin;


  int n_vknot_massbin;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int N_iterations_Kernel;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int iteration_ini;
  
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int N_dm_realizations;


  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int N_dm_initial;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int N_iterations_dm;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Apply_Rankordering;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  vector<int> cwt_used;
  vector<int> cwv_used;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  vector<int> output_at_iteration;


  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int n_cwt;

  int n_cwv;

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////
  // PATCHY PARS
  /**
   * @brief 
  **/
  string dataFileName;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  int inputmode;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  int seed;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  int seed_ref;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  bool runsim;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  bool runv;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  bool diffcosmorz;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  string ic_power_file;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  string ic_WN_file;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  int sfmodel;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  string fastpmpath;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  bool lognden;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  string dir;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  bool transf;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  bool readPS;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  real_prec slength;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  real_prec slengthv;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  real_prec vslength;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  real_prec velbias;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  int Nchunk;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  real_prec biasE;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  real_prec biasepsilon; 
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  real_prec biasrhoexp;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  real_prec biasone;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  real_prec biassign;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  real_prec biassign2;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  real_prec biasepsilon2;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  real_prec biasrhoexp2;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  real_prec devpois;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  real_prec deltathH;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  real_prec Nmean;
  real_prec deltath;
  real_prec cs2;
  real_prec cs3;
  real_prec cst;
  real_prec cpsi;
  real_prec cdeltas2;
  real_prec biasL;
  real_prec sfac;
  real_prec ep;
  real_prec xllc;
  real_prec yllc;
  real_prec zllc;
  real_prec xobs;
  real_prec yobs;
  real_prec zobs;
  bool Normalize_initial_redshift;
  

  real_prec d1;
  real_prec d2;
  real_prec d3;


  real_prec d1_low;
  real_prec d2_low;
  real_prec d3_low;


  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************

  
 public:
  /**
   *  @brief default constructor
   *  @return object of class Parameters
   */
  Params(){}

  //////////////////////////////////////////////////////////
  
  /**
   *  @brief constructor 
   *  @param parameters_file parameter file
   *  @return object of class Parameters
   */
  
 Params(string _par_file):par_file (_par_file)
  {
    this->init_pars();
    this->read_pars(par_file);
    this->NGRID=static_cast<ULONG>(this->Nft*this->Nft*this->Nft);
  }
  //////////////////////////////////////////////////////////
  
  /**
   *  @brief default destructor
   *  @return none
   */
  ~Params(){}


  //////////////////////////////////////////////////////////
  /**
   * @brief output directory
   **/
  string dir_output;

  //////////////////////////////////////////////////////////
  /**
   * @brief output directory
   **/

  int Nft;

  //////////////////////////////////////////////////////////
  /// output file for the FKP power spectrum
  string file_power;

  //////////////////////////////////////////////////////////

  /// the column where the first object coordinate (according to the system of coordinates of the catalog) is written
  /**
   * @brief
  **/
  int i_coord1_g;

  //////////////////////////////////////////////////////////
  /// the column where the second object coordinate (according to the system of coordinates of the catalog) is written
  /**
   * @brief
   **/
  int i_coord2_g;

  //////////////////////////////////////////////////////////
  /// the column where the third object coordinate (according to the system of coordinates of the catalog) is written
  /**
   * @brief
  **/
  int i_coord3_g;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  **/
  int i_v1_g;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  **/
  int i_v2_g;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  **/
  int i_v3_g;

  //////////////////////////////////////////////////////////
  /**
   * @brief  column with the mass of the tracer
  **/
  int i_mass_g;
  int i_vmax_g;

  bool weight_with_mass;

  //////////////////////////////////////////////////////////
  /**
   * @brief Number of cells (per dimension) used in the process of finding closest neighbours in the context of collapsing random placed particles towards the closest dark amtter particles in Patchy.
   **/

  int Nft_random_collapse;
  //////////////////////////////////////////////////////////
  /**
   * @brief Fraction of distance to clsoest dm particle used to collapse randoms
   **/
  real_prec Distance_fraction;
  
  //////////////////////////////////////////////////////////
  /*
   * @brief Read input parameter file and allocate var as private variables of Params class.
   */
  string mass_assignment_scheme; 

  string file_catalogue;
  
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool MAS_correction;
  //////////////////////////////////////////////////////////
  ///  Use Poisson shot-noise correction (yes/no)
  bool SN_correction; 

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Name_survey;
  
  //////////////////////////////////////////////////////////
  /*
   * @brief Initialize all the parameteres defined as private variables of the class Params.
   */
  void init_pars();
  
  //////////////////////////////////////////////////////////
  /*
   * @brief Read input parameter file and allocate var as private variables of Params class.
   */
  void read_pars(string );

  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  int _NX(){return this->NX;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  int _NY(){return this->NY;}

  int _NY_MASS(){return this->NY_MASS;}
  
  int _NY_SAT_FRAC(){return this->NY_SAT_FRAC;}
    ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  int _Nlambdath(){return this->Nlambdath;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  string _Output_directory(){return this->Output_directory;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  string _Input_Directory_Y(){return this->Input_Directory_Y;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  string _Name_Catalog_Y(){return this->Name_Catalog_Y;}

  string _Name_Catalog_Y_HR(){return this->Name_Catalog_Y_HR;}


  string _Name_Catalog_Y_MWEIGHTED(){return this->Name_Catalog_Y_MWEIGHTED;}

  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  string _Input_Directory_X(){return this->Input_Directory_X;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  string _Input_Directory_X_REF(){return this->Input_Directory_X_REF;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  string _XNAME(){return this->XNAME;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  string _Name_Catalog_X(){return this->Name_Catalog_X;}

  //////////////////////////////////////////////////////////
   /**
    *  @brief get the value of the private member
    *  @return
    */
  string _Name_VelFieldx_X(){return this->Name_VelFieldx_X;}
   //////////////////////////////////////////////////////////
    /**
     *  @brief get the value of the private member
     *  @return
     */
  string _Name_VelFieldy_X(){return this->Name_VelFieldy_X;}
    //////////////////////////////////////////////////////////
     /**
      *  @brief get the value of the private member
      *  @return
      */
  string _Name_VelFieldz_X(){return this->Name_VelFieldz_X;}





  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  string _Name_Catalog_X_REF_PDF(){return this->Name_Catalog_X_REF_PDF;}

 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  string _Name_Catalog_X_NEW(){return this->Name_Catalog_X_NEW;}

  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  string _Name_Property_X(){return this->Name_Property_X;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  string _YNAME(){return this->YNAME;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  string _Name_Property_Y(){return this->Name_Property_Y;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  int _iMAS_X(){return this->iMAS_X;}

 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  int _iMAS_X_NEW(){return this->iMAS_X_NEW;}

  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  int _iMAS_X_REF_PDF(){return this->iMAS_X_REF_PDF;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  int _iMAS_Y(){return this->iMAS_Y;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  real_prec _delta_Y_max(){return this->delta_Y_max;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  real_prec _delta_Y_min(){return this->delta_Y_min;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  real_prec _delta_X_max(){return this->delta_X_max;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  real_prec _delta_X_min(){return this->delta_X_min;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  real_prec _ldelta_Y_max(){return this->ldelta_Y_max;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  real_prec _ldelta_Y_min(){return this->ldelta_Y_min;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  real_prec _ldelta_X_max(){return this->ldelta_X_max;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  real_prec _ldelta_X_min(){return this->ldelta_X_min;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  string _Quantity(){return this->Quantity;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  int _NMASSbins(){return this->NMASSbins;}

  int _NMASSbins_mf(){return this->NMASSbins_mf;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  real_prec _redshift(){return this->redshift;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  real_prec _smscale(){return this->smscale;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  int _realization(){return this->realization;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  bool _Comp_conditional_PDF(){return this->Comp_conditional_PDF;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  bool _Comp_joint_PDF(){return this->Comp_joint_PDF;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  bool _write_files_for_histograms(){return this->write_files_for_histograms;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  bool _Redefine_limits(){return this->Redefine_limits;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  bool _Convert_Density_to_Delta_X(){return this->Convert_Density_to_Delta_X;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  bool _Convert_Density_to_Delta_Y(){return this->Convert_Density_to_Delta_Y;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  real_prec _lambdath(){return this->lambdath;}


  real_prec _lambdath_v(){return this->lambdath_v;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  bool _Write_Scatter_Plot(){return this->Write_Scatter_Plot;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  bool _Write_PDF_number_counts(){return this->Write_PDF_number_counts;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  string _Scale_X(){return this->Scale_X;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  string _Scale_Y(){return this->Scale_Y;}


   ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  int _N_dm_initial(){return this->N_dm_initial;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  int _n_sknot_massbin(){return this-> n_sknot_massbin;}


  int _n_vknot_massbin(){return this-> n_vknot_massbin;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  int _N_iterations_Kernel(){return this->N_iterations_Kernel;}

  
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  int _N_dm_realizations(){return this->N_dm_realizations;}



  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  int _iteration_ini(){return this->iteration_ini;}

  
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  int _N_iterations_dm(){return this->N_iterations_dm;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  bool _Apply_Rankordering(){return this->Apply_Rankordering;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  vector<int> _cwt_used(){return this->cwt_used;}

  vector<int> _cwv_used(){return this->cwv_used;}
 ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  vector<int> _output_at_iteration(){return this->output_at_iteration;}

  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  @return 
   */
  int _n_cwt(){return this->n_cwt;}

  int _n_cwv(){return this->n_cwv;}

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //   Parameters of COSMOLIB    ///////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  
  // function to get private variables
  /**
   *  @brief get the value of the private member statistics
   *  @return statistics
   */
  string _statistics () {return statistics;}
  //////////////////////////////////////////////////////////  

  /**
   *  @brief get the value of the private member dir_input
   *  @return dir_input
   */
  string _dir_input () {return dir_input;}
  //////////////////////////////////////////////////////////
  /**
   * @brief ask if cross is desidered
  **/
  bool measure_cross;

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member dir_output
   *  @return dir_output
   */
  string _dir_output(){return this->dir_output;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member dir_output
   *  @return dir_output
   */
  string _input_type () {return input_type;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  unsigned long _ngal_delta() {return ngal_delta;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _delta_grid_file () {return delta_grid_file;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _delta_grid_file2 () {return delta_grid_file2;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _delta_grid_file3 () {return delta_grid_file3;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _delta_grid_file4 () {return delta_grid_file4;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  bool _measure_cross () {return measure_cross;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _measure_cross_from_1 () {return measure_cross_from_1;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _measure_cross_from_2 () {return measure_cross_from_2;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief set the value of the private member dir_input
   *  @param _dir_input input directory where the object and random
   *  catalogues are stored
   *  @return none
   */
  void set_dir_input (string _dir_input) {dir_input = _dir_input;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief set the value of the private member dir_output
   *  @param _dir_output output directory
   *  @return none
   */
  void set_dir_output (string _dir_output) {dir_output = _dir_output;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member file_catalogue
   *  @return file_catalogue
   */
  string _file_catalogue () {return file_catalogue;}
    
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member file_random
   *  @return file_random
   */
  string _file_random () {return file_random;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member new_Lbox
   *  @return new_Lbox
   */
  bool _new_Lbox () {return new_Lbox;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member sys_of_coord_g
   *  @return sys_of_coord_g
   */
  int _sys_of_coord_g () {return sys_of_coord_g;}

    //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member sys_of_coord_g
   *  @return sys_of_coord_dm
   */
  int _sys_of_coord_dm () {return sys_of_coord_dm;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord1_g
   *  @return i_coord1_g
   */
  int _i_coord1_g () {return i_coord1_g;}

  //////////////////////////////////////////////////////////  
  /**
   *  @brief get the value of the private member i_coord2_g
   *  @return i_coord2_g
   */
  int _i_coord2_g () {return i_coord2_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_g
   *  @return i_coord3_g
   */
  int _i_coord3_g () {return this->i_coord3_g;}

    //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord1_g
   *  @return i_coord1_g
   */
  int _i_v1_g () {return i_v1_g;}

  //////////////////////////////////////////////////////////  
  /**
   *  @brief get the value of the private member i_coord2_g
   *  @return i_coord2_g
   */
  int _i_v2_g () {return i_v2_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_g
   *  @return i_coord3_g
   */
  int _i_v3_g () {return this->i_v3_g;}


  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_g
   *  @return i_coord3_g
   */
  int _i_mass_g () {return this->i_mass_g;}

  int _i_vmax_g () {return this->i_vmax_g;}


  int _i_sf_g () {return this->i_sf_g;}

    //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_g
   *  @return i_coord3_g
   */
  int _i_mass_dm () {return this->i_mass_dm;}

  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight1_g
   *  @return i_weight1_g
   */
  int _i_weight1_g () {return this->i_weight1_g;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight2_g
   *  @return i_weight2_g
   */
  int _i_weight2_g () {return this->i_weight2_g;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight3_g
   *  @return i_weight3_g
   */
  int _i_weight3_g () {return this->i_weight3_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight4_g
   *  @return i_weight4_g
   */
  int _i_weight4_g () {return i_weight4_g;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight1_g
   *  @return use_weight1_g
   */
  bool _use_weight1_g () {return this->use_weight1_g;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight2_g
   *  @return use_weight2_g
   */
  bool _use_weight2_g () {return this->use_weight2_g;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight3_g
   *  @return use_weight3_g
   */
  bool _use_weight3_g () {return this->use_weight3_g;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight4_g
   *  @return use_weight4_g
   */
  bool _use_weight4_g () {return this->use_weight4_g;}


  bool _weight_vel_with_mass () {return this->weight_vel_with_mass;}




  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord1_g
   *  @return i_coord1_g
   */
  int _i_coord1_dm () {return i_coord1_dm;}

  //////////////////////////////////////////////////////////  
  /**
   *  @brief get the value of the private member i_coord2_g
   *  @return i_coord2_g
   */
  int _i_coord2_dm () {return i_coord2_dm;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_g
   *  @return i_coord3_g
   */
  int _i_coord3_dm () {return this->i_coord3_dm;}

    //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord1_g
   *  @return i_coord1_g
   */
  int _i_v1_dm () {return i_v1_dm;}

  //////////////////////////////////////////////////////////  
  /**
   *  @brief get the value of the private member i_coord2_g
   *  @return i_coord2_g
   */
  int _i_v2_dm () {return i_v2_dm;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_g
   *  @return i_coord3_g
   */
  int _i_v3_dm() {return this->i_v3_dm;}

  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member Name_survey
   *  @return Name_survey
   */
  string _Name_survey () {return this->Name_survey;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_mean_density_g
   *  @return i_mean_density_g
   */
  int _i_mean_density_g () {return i_mean_density_g;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member angles_units_g
   *  @return angles_units_g
   */
  string _angles_units_g () {return angles_units_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_random_catalog
   *  @return use_random_catalog
   */
  bool _use_random_catalog () {return use_random_catalog;}


  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_random_catalog
   *  @return use_random_catalog
   */
  bool _use_random_catalog_cl () {return use_random_catalog_cl;}
  

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member new_los
   *  @return new_los
   */
  bool _new_los () {return new_los;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member sys_of_coord_r
   *  @return sys_of_coord_r
   */
  int _sys_of_coord_r () {return sys_of_coord_r;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord1_r
   *  @return i_coord1_r
   */
  int _i_coord1_r () {return i_coord1_r;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord2_r
   *  @return i_coord2_r
   */
  int _i_coord2_r () {return i_coord2_r;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_r
   *  @return i_coord3_r
   */
  int _i_coord3_r () {return i_coord3_r;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight1_r
   *  @return i_weight1_r
   */
  int _i_weight1_r () {return i_weight1_r;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight2_r
   *  @return i_weight2_r
   */
  int _i_weight2_r () {return i_weight2_r;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight3_r
   *  @return i_weight3_r
   */
  int _i_weight3_r () {return i_weight3_r;}

   //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight4_r
   *  @return i_weight4_r
   */
  int _i_weight4_r () {return i_weight4_r;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property1_r
   *  @return i_property1_r
   */
  int _i_property1_r () {return i_property1_r;}

  //////////////////////////////////////////////////////////  
  /**
   *  @brief get the value of the private member i_property2_r
   *  @return i_property2_r
   */
  int _i_property2_r () {return i_property2_r;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property3_r
   *  @return i_property3_r
   */
  int _i_property3_r () {return i_property3_r;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property4_r
   *  @return i_property4_r
   */
  int _i_property4_r () {return i_property4_r;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property5_r
   *  @return i_property5_r
   */
  int _i_property5_r () {return i_property5_r;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property6_r
   *  @return i_property6_r
   */
  int _i_property6_r () {return i_property6_r;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property1_g
   *  @return i_property1_g
   */
  int _i_property1_g () {return i_property1_g;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property2_g
   *  @return i_property2_g
   */
  int _i_property2_g () {return i_property2_g;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property3_g
   *  @return i_property3_g
   */
  int _i_property3_g () {return i_property3_g;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property4_g
   *  @return i_property4_g
   */
  int _i_property4_g () {return i_property4_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property5_g
   *  @return i_property5_g
   */
  int _i_property5_g () {return i_property5_g;}
  //////////////////////////////////////////////////////////  
  /**
   *  @brief get the value of the private member i_property6_g
   *  @return i_property6_g
   */
  int _i_property6_g () {return i_property6_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight1_r
   *  @return use_weight1_r
   */
  bool _use_weight1_r () {return use_weight1_r;}

  bool _weight_with_mass () {return this->weight_with_mass;}

  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member use_weight2_r
   *  @return use_weight2_r
   */
  bool _use_weight2_r () {return use_weight2_r;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member use_weight3_r
   *  @return use_weight3_r
   */
  bool _use_weight3_r () {return use_weight3_r;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member use_weight4_r
   *  @return use_weight4_r
   */
  bool _use_weight4_r () {return use_weight4_r;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member i_mean_density_r
   *  @return i_mean_density_r
   */
  int _i_mean_density_r () {return i_mean_density_r;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member angles_units_r
   *  @return angles_units_r
   */
  string _angles_units_r () {return angles_units_r;}

  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member n_catalogues
   *  @return n_catalogues
   */
  int _n_catalogues () {return n_catalogues;}

  
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member Nft
   *  @return Nft
   */
  int _Nft () {return this->Nft;}


  int _Nft_low () {return this->Nft_low;}

  int _Nft_HR () {return this->Nft_HR;}

  int _Nft_random_collapse () {return this->Nft_random_collapse;}

  real_prec _Distance_fraction() {return this->Distance_fraction;}
  
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member Nft
   *  @return Nft
   */
  ULONG _NGRID(){return this->NGRID;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member Lbox
   *  @return Lbox
   */
  real_prec _Lbox () {return this->Lbox;}

  real_prec _Lbox_low () {return this->Lbox_low;}


  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member mass_assignment_scheme
   *  @return mass_assignment_scheme
   */
  string _mass_assignment_scheme () {return mass_assignment_scheme;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member type_of_binning
   *  @return type_of_binning
   */
  string _type_of_binning () {return type_of_binning;}


  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member type_of_binning
   *  @return type_of_binning
   */
  string _extra_info () {return this->extra_info;}

  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member k_bin_step
   *  @return k_bin_step
   */
  real_prec _k_bin_step () {return k_bin_step;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member N_log_bins
   *  @return N_log_bins
   */
  int _N_log_bins () {return N_log_bins;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member ndel_data 
   *  @return ndel_data 
   */
  int _ndel_data () {return this->ndel_data;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member ndel_window
   *  @return ndel_window
   */
  int _ndel_window () {return this->ndel_window;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member N_mu_bins
   *  @return N_mu_bins
   */
  int _N_mu_bins () {return N_mu_bins;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member MAS_correction 
   *  @return MAS_correction 
   */
  bool _MAS_correction () {return MAS_correction;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member FKP_weight
   *  @return FKP_weight
   */
  bool _FKP_weight () {return FKP_weight;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member SN_correction
   *  @return SN_correction
   */
  bool _SN_correction () {return SN_correction;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member FKP_error_bars
   *  @return FKP_error_bars
   */
  bool _FKP_error_bars () {return FKP_error_bars;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member FKP_error_bars_exact
   *  @return FKP_error_bars_exact
   */
  bool _FKP_error_bars_exact () {return FKP_error_bars_exact;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member Pest
   *  @return Pest
   */
  real_prec _Pest () {return Pest;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member nbar_tabulated
   *  @return nbar_tabulated
   */
  bool _nbar_tabulated () {return nbar_tabulated;}
  ///////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member constant_depth
   *  @return constant_depth 
   */
  bool _constant_depth () {return constant_depth;}
  ////////////////////////////////////////////////////////// 
  
  /**
   *  @brief get the value of the private member N_z_bins
   *  @return N_z_bins
   */
  int _N_z_bins () {return N_z_bins;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member redshift_min_sample
   *  @return redshift_min_sample
   */
  real_prec _redshift_min_sample () {return redshift_min_sample;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member redshift_max_sample
   *  @return redshift_max_sample
   */
  real_prec _redshift_max_sample () {return redshift_max_sample;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member N_dndz_bins
   *  @return N_dndz_bins
   */
  int _N_dndz_bins () {return N_dndz_bins;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member new_N_dndz_bins
   *  @return new_N_dndz_bins 
   */
  int _new_N_dndz_bins () {return new_N_dndz_bins;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member area_survey
   *  @return area_survey 
   */
  real_prec _area_survey () {return area_survey;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member Healpix_resolution
   *  @return Healpix_resolution 
   */
  int _Healpix_resolution () {return Healpix_resolution;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member file_dndz
   *  @return file_dndz
   */
  string _file_dndz () {return file_dndz;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member file_power
   *  @return file_power 
   */
  string _file_power () {return file_power;}
  
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member file_power_log
   *  @return file_power_log 
   */
  string _file_power_log () {return file_power_log;}
  ////////////////////////////////////////////////////////// 

  /**
   *  @brief get the value of the private member file_window
   *  @return file_window 
   */
  string _file_window () {return file_window;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member file_power2d
   *  @return file_power2d 
   */
  string _file_power2d () {return file_power2d;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member file_power2d_mk
   *  @return file_power2d_mk
   */
  string _file_power2d_mk () {return file_power2d_mk;}

  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member kmin_bk;
   *  @return kmin_bk
   */
  real_prec _kmin_bk () {return kmin_bk;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member kmin_bk;
   *  @return kmax_bk
   */
  real_prec _kmax_bk () {return kmax_bk;}


  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member use_fundamental_mode_as_kmin_bk
   *  @return Nshells_bk
   */
  bool _use_fundamental_mode_as_kmin_bk () {return use_fundamental_mode_as_kmin_bk;}

 

  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member file_bispectrum
   *  @return file_bispectrum
   */
  string _file_bispectrum () {return file_bispectrum;}



  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member kmax_y_ds
   *  @return kmax_y_ds
   */
  real_prec _kmax_y_ds () {return kmax_y_ds;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member input_file_mask
   *  @return input_file_mask
   */
  string _input_file_mask(){return input_file_mask;}
  ////////////////////////////////////////////////////////// 
  /**  
   *  @brief get the value of the the private member i_mask_pixel
   *  @return i_mask_pixel
   */
  int _i_mask_pixel(){return i_mask_pixel;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member i_mask_alpha
   *  @return  i_mask_alpha
   */
  int _i_mask_alpha(){return i_mask_delta;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member i_mask_delta
   *  @return  i_mask_delta
   */
  int _i_mask_delta(){return i_mask_delta;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member i_mask_flag
   *  @return  i_mask_flag
   */
  int _i_mask_flag(){return i_mask_flag;}


  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member 
   *  @return  
   */
  string _Name_redshift_mask(){return this->Name_redshift_mask;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member 
   *  @return  
   */
  string _Name_binary_mask(){return this->Name_binary_mask;}

  

  //////////////////////////////////////////////////////////
  /**
   * @brief Type if input file, options
   * are "cat" meaning catalog, and "grid_delta" meaning
   * that the input is alread the delta in the grid
  **/
  string input_type;

 
   ////////////////////////////////////////////////////////// 
  ////////////////////////////////////////////////////////// 
  ////////////////////////////////////////////////////////// 
  ////////////////////////////////////////////////////////// 
  ////////////////////////////////////////////////////////// 
  ////////////////////////////////////////////////////////// 
  ////COSMOLOGICAL PARAMETERS/////////////////////////////// 
  ////////////////////////////////////////////////////////// 
  ////////////////////////////////////////////////////////// 
  ////////////////////////////////////////////////////////// 



  /**
   *  @brief get the value of the private member &Omega;<SUB>M</SUB>
   *  @return &Omega;<SUB>M</SUB>
   */
  real_prec _om_matter () { return om_matter; }

  /**
   *  @brief get the value of the private member &Omega;<SUB>M</SUB>
   *  @return &Omega;<SUB>M</SUB>
   */
  real_prec _om_cdm () { return om_cdm; }
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member &Omega;<SUB>rad</SUB>: the radiation density at z=0
   *  @return &Omega;<SUB>rad</SUB>
   */
  real_prec _om_radiation () { return om_radiation; }
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member &Omega;<SUB>b</SUB>: the baryon density at z=0
   *  @return &Omega;<SUB>b</SUB>
   */
  real_prec _om_baryons () { return om_baryons; }
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member &Omega;<SUB>DE</SUB>: the dark energy density at z=0
   *  @return &Omega;<SUB>DE</SUB>
   */
  real_prec _om_vac () { return om_vac; }
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member &Omega;<SUB>k</SUB>: the density of curvature energy
   *  @return &Omega;<SUB>k</SUB>
   */
  real_prec _om_k () { return om_k; }
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member H<SUB>0</SUB>: the Hubble constant at z=0 [km/sec/Mpc] 
   *  @return H<SUB>0</SUB>
   */
  real_prec _Hubble () { return Hubble; }
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member \e h: the Hubble parameter, H<SUB>0</SUB>/100
   *  @return \e h
   */
  real_prec _hubble () { return hubble; }
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member n<SUB>spec</SUB>: the primordial spectral index
   *  @return n<SUB>spec</SUB>
   */
  real_prec _spectral_index () { return spectral_index; }
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member w<SUB>0</SUB>: the parameter of the dark energy equation of state
   *  @return w<SUB>0</SUB>
   */
  real_prec _w_eos () { return w_eos; }
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member N<SUB>eff</SUB>: the effective number (for QED + non-instantaneous decoupling)
   *  @return N<SUB>eff</SUB> 
   */
  real_prec _N_eff () { return N_eff; }
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member &sigma;<SUB>8</SUB>: the power spectrum normalization
   *  @return &sigma;<SUB>8</SUB>
   */
  real_prec _sigma8 () { return sigma8; }

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member T<SUB>CMB</SUB>: the present day CMB temperature [K]
   *  @return T<SUB>CMB</SUB>
   */
  real_prec _Tcmb () { return Tcmb; }


  bool _use_wiggles () { return use_wiggles; }

  bool _RR () { return RR; }
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /*
   * @brief This is made as public for it has to be modified under iterative processes
   */
  string Output_directory;

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NMASSbins_mf;


  //////////////////////////////////////////////////////////
  /**
   * @brief
  **/
  int masskernel;

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  // Params useful for Patchy. OJO QUE EN PARAMS.cpp aun no estan leyendo muchos de estos
  int _inputmode(){return this->inputmode;}
  string _dataFileName(){return this->dataFileName;}
  string _ic_power_file(){return this->ic_power_file;}
  string _ic_WN_file(){return this->ic_WN_file;}
  string _ic_file(){return this->ic_file;}
  bool _use_ic_file(){return this->use_ic_file;}
  string _ic_WN_dir(){return this->ic_WN_dir;}
  string _dir(){return this->dir;}
  string _fastpmpath(){return this->fastpmpath;}
  int _sfmodel(){return this->sfmodel;}
  bool _transf(){return this->transf;}
  bool _Normalize_initial_redshift(){return this->Normalize_initial_redshift;}
  bool _runsim(){return runsim;}
  bool _runv(){return runv;}
  bool _lognden(){return lognden;}
  bool _diffcosmorz(){return diffcosmorz;}
  bool _readPS(){return readPS;}
  real_prec _Nmean(){return Nmean;}
  real_prec _cs2(){return cs2;}
  real_prec _cdeltas2(){return cdeltas2;}
  real_prec _cpsi(){return cpsi;}
  real_prec _cst(){return cst;}
  real_prec _cs3(){return cs3;}
  real_prec _sfac(){return sfac;}
  real_prec _ep(){return ep;}
  real_prec _devpois(){return devpois;}
  //real_prec devpois2=params.find<real_prec>("devpois2");
  real_prec _deltathH(){return deltathH;}
  //real_prec deltathH2=params.find<real_prec>("deltathH2");
  ULONG _N1(){return this->Nft;}
  ULONG _Nchunk(){return this->Nchunk;}
  real_prec _L1(){return this->Lbox;}
  real_prec _biasL(){return this->biasL;}
  real_prec _biasE(){return this->biasE;}
  real_prec _biasrhoexp(){return this->biasrhoexp;}
  real_prec _biasepsilon(){return this->biasepsilon;}
  real_prec _biasrhoexp2(){return this->biasrhoexp2;}
  real_prec _biasepsilon2(){return this->biasepsilon2;}
  real_prec _biassign(){return this->biassign;}
  real_prec _biassign2(){return this->biassign2;}
  real_prec _biasone(){return this->biasone;}
  real_prec _velbias(){return this->velbias;}
  real_prec _xllc(){return this->xllc;}
  real_prec _yllc(){return this->yllc;}
  real_prec _zllc(){return this->zllc;}
  real_prec _xobs(){return this->xobs;}
  real_prec _yobs(){return this->yobs;}
  real_prec _zobs(){return this->zobs;}
  ULONG _seed(){return this->seed;}
  ULONG _seed_ref(){return this->seed_ref;}
  int _masskernel(){return this->masskernel;}
  real_prec _slength(){return this->slength;}
  real_prec _slengthv(){return this->slengthv;}
  real_prec _vslength(){return this->vslength;}
  real_prec _deltath(){return this->deltath;}

  real_prec _d1(){return this->d1;}
  real_prec _d2(){return this->d2;}
  real_prec _d3(){return this->d3;}

  real_prec _d1_low(){return this->d1_low;}
  real_prec _d2_low(){return this->d2_low;}
  real_prec _d3_low(){return this->d3_low;}


  real_prec _LOGMASSmin(){return this->LOGMASSmin;}
  real_prec _LOGMASSmax(){return this->LOGMASSmax;}

  real_prec _VMAXmin(){return this->VMAXmin;}
  real_prec _VMAXmax(){return this->VMAXmax;}

  real_prec _MASS_units(){return this->MASS_units;}
  real_prec _Redshift_initial(){return this->Redshift_initial;}

  
};

#endif


