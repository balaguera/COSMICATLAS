// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
/**
 *  @file Params.h
 *  @brief The class Parameters. Reads the parameter
 */
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************

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
#include "NumericalMethods.h"  // def.h is included in the NumericalMethods.h file
#include "Type_structures_def.h"

class Params
{



 private :
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string statistics;

  //////////////////////////////////////////////////////////
  /**
   *  @name 
   */
  real_prec Initial_Redshift_DELTA;

  //////////////////////////////////////////////////////////
  /**
   * @brief input directory where the tracers, dm and random catalogues are stored
   **/
  string Input_dir_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief input directory where the tracers, dm and random catalogues are stored
   **/
  string Input_dir_cat_TWO;

  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  string ic_WN_dir;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/

  string ic_file;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/

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
  //////////////////////////////////////////////////////////
   /**
   * @brief 
  **/
  string Name_redshift_mask;
  //////////////////////////////////////////////////////////
   /**
   * @brief 
  **/
  string Name_binary_mask;
  //////////////////////////////////////////////////////////
   /**
   * @brief 
  **/
  string type_of_object;
  //////////////////////////////////////////////////////////
   /**
   * @brief 
  **/
  int IC_index;
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
  int i_mass_r;

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

  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  int Number_of_chunks_new_dm;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  bool weight_vel_with_mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec MASS_units;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec LOGMASSmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec LOGMASSmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/

  real_prec VMAXmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec VMAXmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec RSmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec RSmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec SPINmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec SPINmax;

  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  int Number_of_MultiLevels;
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
  /**
   * @brief
   **/
  int i_weight1_g;

  //////////////////////////////////////////////////////////
  /**
   * @brief The column where the second object weight is written

   **/
  int i_weight2_g;

  //////////////////////////////////////////////////////////
  /**
   * @brief The column where the third object weight is written
   **/
  int i_weight3_g;

  //////////////////////////////////////////////////////////
  /**
   * @brief The column where the third object weight is written
   **/
  int i_weight4_g;

  //////////////////////////////////////////////////////////
  /**
   * @brief The column where the halo mean densityis written
   **/
  int i_mean_density_g;

  //////////////////////////////////////////////////////////
  /// the column where the object T/|W| ratio  (0.5=virialization)
  /**
   * @brief The column where the third object weight is written
   **/
  int i_virial_g;
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
   * @brief
   **/

  real_prec M_exclusion;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec Deltal;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec Deltamu;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  int Number_of_references;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/

  int Number_of_new_mocks;
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
  bool use_weight3_r;

  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  bool use_weight4_r;
  //////////////////////////////////////////////////////////

  bool get_distribution_min_separations;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */

  real_prec Prop_threshold_multi_scale_1;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec Prop_threshold_multi_scale_2;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec Prop_threshold_multi_scale_3;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec Prop_threshold_multi_scale_4;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec Prop_threshold_multi_scale_5;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  ULONG Nft_low_l1;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  ULONG Nft_low_l2;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  ULONG Nft_low_l3;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  ULONG Nft_low_l4;


  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec Tolerance_factor_l1;

  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec Tolerance_factor_l2;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec Tolerance_factor_l3;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec Tolerance_factor_l4;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec Tolerance_factor_l5;

  //////////////////////////////////////////////////////////

  /**
   *  @name FKP power spectrum
   */
  ///@{
  //////////////////////////////////////////////////////////
  /// Number of grid cells /per dimension for the Discrete Fourier Trasnform


  ULONG Nft_low;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/

  ULONG Nft_HR;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec Lbox_low;

  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec vkernel_exponent;
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
  /// Number of log-spaced bins in Fourier space
  int N_log_bins;

  //////////////////////////////////////////////////////////
  /**
   * @brief Size of bins in Foureir space
   * @details The size of k-bins is ndel times the fundamental mode
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
  /**
   * @brief
   **/
   bool use_vel_kernel;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  int n_catalogues;

  //////////////////////////////////////////////////////////
  /**
   * @brief name of the output file used to store the bispectrum
   **/
  string file_bispectrum;

  //////////////////////////////////////////////////////////
  /// These parameters is used to define the shells in k-space
  /**
   * @brief
   **/
  bool use_fundamental_mode_as_kmin_bk;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec kmin_bk;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec delta_x;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec delta_y;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec delta_z;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec deltak_x;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec deltak_y;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec deltak_z;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec deltak_0;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec kmin;

  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec kmax;

  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec DeltaK_data;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec DeltaK_window;


  /////////////////////////////////////////////////////////////
  /**
   * @brief   Maximum k-value for constructing k-bins

   **/
  real_prec kmax_bk;
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
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec n_s;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec alpha_s;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec RR;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec M_reference;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec Delta_SO;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  bool use_wiggles;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  bool fixed_redshift;

  //////////////////////////////////////////////////////////
  ///Type of binning in Fourier space: linear/log
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  string file_bin_x_coord;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  string file_bin_y_coord;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  string file_bin_z_coord;

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
  ULONG NGRID_h;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  ULONG N_lines_binary;

  //////////////////////////////////////////////////////////

  /**
   * @brief Number of bins in Y to construct BIAS
   * If the Y- density field is NGP, this values is transformed
   * to the maximum number of TR in one cell.
   */
  int NY;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of bins in X to construct BIAS
   */
  int NX;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of bins in Y to construct BIAS
   * If the Y- density field is NGP, this values is transformed
   * to the maximum number of TR in one cell.
   */
  int NY_MASS;
  //////////////////////////////////////////////////////////
  /*
      * @brief Number of bins in Y to construct BIAS
   * If the Y- density field is NGP, this values is transformed
   * to the maximum number of TR in one cell.

   */
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
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/

  string Name_Catalog_Y_HR;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/

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
  string Input_Directory_X_REF_TWO;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Input_Directory_BIAS_KERNEL;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Input_Directory_BIAS_KERNEL_TWO;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Input_Directory_X_NEW;
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
  string Name_Catalog_X_NEW_TWO;
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
  /**
   * @brief Maximum value of the TR overdensity
   * @details If requiested from parameter filw, this value is overloaded
   * by computing it from the input fields
   */
  real_prec delta_Y_max;
  //////////////////////////////////////////////////////////

  /**
   * @brief Minimum value of the TR overdensity
   * @details If requiested from parameter filw, this value is overloaded
   * by computing it from the input fields
   */
  real_prec delta_Y_min;
  //////////////////////////////////////////////////////////
  /**
   * @brief Maximum value of the DM overdensity
   * @details If requiested from parameter filw, this value is overloaded
   * by computing it from the input fields
   */
  real_prec delta_X_max;

 //////////////////////////////////////////////////////////
  /*
    * @brief Minimum value of the DM overdensity
   * @details If requiested from parameter filw, this value is overloaded
   * by computing it from the input fields
   */
  real_prec delta_X_min;


  /////////////////////////////////////////////////////////
  /**
   * @brief Maximum value of the log(1+overdensity) for TR
   * @detail If requiested from parameter filw, this value is overloaded
   * by computing it from the input fields
   */
  real_prec ldelta_Y_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Minimum value of the log(1+overdensity) for TR
   * @detail If requiested from parameter filw, this value is overloaded
   * by computing it from the input fields
   */
  real_prec ldelta_Y_min;

 //////////////////////////////////////////////////////////
  /*
   * @brief
   * @brief Maximum value of the log(1+overdensity) for DM
   * @detail If requiested from parameter filw, this value is overloaded
   * by computing it from the input fields
   */
  real_prec ldelta_X_max;
  //////////////////////////////////////////////////////////
  /*
  * @brief Minimum value of the log(1+overdensity) for DM
   * @detail If requiested from parameter filw, this value is overloaded
   * by computing it from the input fields
   */
  real_prec ldelta_X_min;





  //////////////////////////////////////////////////////////
  /*
    * @brief Identification of the input density field (density, delta)
   */
  string Quantity;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NMASSbins;
 //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NMASSbins_power;
   //////////////////////////////////////////////////////////
    /**
   * @brief Cosmological redshift
   * @brief Read from parameter file
   */
  real_prec redshift;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  real_prec smscale;
  //////////////////////////////////////////////////////////
  /*
   * @brief Identification for a realization
   * @brief Read from parameter file
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
  /**
   * @brief Container with the LOS of new DM fields used in the case in which seveal LOS are built simultaneously
   */
  vector<int> list_new_dm_fields;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container with the LOS of new DM fields used in the case in which seveal LOS are built simultaneously
   */
  vector<int> list_Nft_MultiLevels;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container with the LOS of new DM fields used in the case in which seveal LOS are built simultaneously
   */
  vector<ULONG> list_Ntracers_MultiLevels;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container with the LOS of new DM fields used in the case in which seveal LOS are built simultaneously
   */
  vector<real_prec> list_Props_Threshold_MultiLevels;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container with the LOS of new DM fields used in the case in which seveal LOS are built simultaneously
   */
  vector<real_prec> list_Props_Tolerance_MultiLevels;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  vector<string> files_new_dm_fields;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container with the LOS of new DM fields used in the case in which seveal LOS are built simultaneously
   */
  vector<string> files_dm_references;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  vector<int> list_bias_references;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  vector<string> files_bias_references;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  vector<string> files_kernel_references;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  vector<string> files_tracer_references;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  vector<string> files_tracer_field_references;


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
   * @brief Specify whether the histograms are done in log or linear scale for Y
   * @brief Read from parameter file
   */
  string Scale_Y;
  //////////////////////////////////////////////////////////
  /*
   * @brief Specify whether the histograms are done in log or linear scale for X
  */
  string Scale_X;
  //////////////////////////////////////////////////////////
  /*
     * @brief Number of bins in the Mk info to build BIAS.
   */
  int n_sknot_massbin;
  //////////////////////////////////////////////////////////
  /*
   * @brief Number of bins in the Veldisp info to build BIAS.
   */
  int n_vknot_massbin;

  //////////////////////////////////////////////////////////
  /**
   * @brief Number of BAM iterations
   * @brief Read from parameter file
   */
  int N_iterations_Kernel;
 //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int iteration_ini;

  //////////////////////////////////////////////////////////
  /*
  * @brief Number of DM that will be created by the approximated gravity solver
   * @ and on which the Kernel will act to create a halo mock density field.
   */
  int N_dm_realizations;
  //////////////////////////////////////////////////////////
  /*
   * @brief Initial label of the DM realizations
   */
  int N_dm_initial;
  //////////////////////////////////////////////////////////
  /*
   * @brief Number of iterations for the pre-processing of DM
   * Parameter set in the parameter file.
   * Default 0
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
  bool Apply_Rankordering_ab_initio;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  vector<int> cwt_used;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  vector<int> cwv_used;

  //////////////////////////////////////////////////////////
  /*
   * @brief  Container specifying iterations at which some outputs are produced
   */
  vector<int> output_at_iteration;

  //////////////////////////////////////////////////////////
  /**
   * @brief Number of CWT used, computed as the size of the
   * Bam::cwt_used container
   */

  int n_cwt;
  //////////////////////////////////////////////////////////
  /*
   * @brief Number of CWT based on the shear of velocity field used
   */

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
  bool ic_alias_corrected;


  //////////////////////////////////////////////////////////
  /**
   * @brief
  **/
  string ic_input_type;

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
  real_prec velbias_dm;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  **/
  real_prec velbias_random;
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
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec deltath;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec cs2;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec cs3;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec cst;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec cpsi;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec cdeltas2;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec biasL;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec sfac;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec ep;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec xllc;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec yllc;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec zllc;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec xobs;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec yobs;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec zobs;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  bool Normalize_initial_redshift;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec Initial_Redshift_ic_power_file;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec d1;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec d2;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec d3;

  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */

  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec d1_low;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec d2_low;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec d3_low;
  //////////////////////////////////////////////////////////
  /*
   * @brief Number if cells per dimention in the mesh. 
   */
  ULONG Nft;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */
  bool SN_correction;

  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  string Output_directory;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  int NMASSbins_mf;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  

  string vel_units_g;

  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  

  int masskernel;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  

  int masskernel_vel;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  string dir_output;
  //////////////////////////////////////////////////////////
  /**
   * @brief Size of the box in Mpc/h
   * @brief Read from parameter file
   */
  real_prec Lbox;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  real_prec DeltaKmin;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  string file_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  bool redshift_space_coords_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  int i_coord1_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  int i_coord2_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  int i_coord3_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  int i_v1_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  int i_v2_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  int i_v3_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  int i_mass_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  int i_vmax_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  int i_rs_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  int i_spin_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  bool weight_with_mass;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  ULONG Nft_random_collapse;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  real_prec Distance_fraction;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  string mass_assignment_scheme;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  string file_catalogue;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  bool MAS_correction;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  string Name_survey;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  ULONG Nnp_data;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  ULONG Nnp_window;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  int Unitsim_plabel;
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
  Params(){this->init_pars();}
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
    this->derived_pars();
  }
  //////////////////////////////////////////////////////////

  /**
   *  @brief default destructor
   *  @return none
   */
  ~Params(){}


// These variables smust be made private, 
  // andnif thy need to be modified, de define a method, e.g.
  // void set_Nft(int new_Nft){this->Nft=new_Nft;}

    //////////////////////////////////////////////////////////
  /**
   * @brief Structure containing cosmological parameters
   **/

  s_CosmologicalParameters s_cosmo_pars;
  
    //////////////////////////////////////////////////////////
  /**
   * @brief Initialize all the parameteres defined as private variables of the class Params.
   **/
  void init_pars();

  //////////////////////////////////////////////////////////
  /**
   * @brief Read input parameter file and allocate var as private variables of Params class.
  **/
  void read_pars(string );

  void derived_pars();
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
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */

  int _NY_MASS(){return this->NY_MASS;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */

  int _NY_SAT_FRAC(){return this->NY_SAT_FRAC;}
    //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _Nlambdath(){return this->Nlambdath;}
 //////////////////////////////////////////////////////////
  /**
   *  @brief get/set the value of the private member
   *  @return
   */
  string _Output_directory(){return this->Output_directory;}
  void set_Output_directory(string new_Output_directory){this->Output_directory=new_Output_directory;}
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
  string _Input_Directory_BIAS_KERNEL(){return this->Input_Directory_BIAS_KERNEL;}
 //////////////////////////////////////////////////////////

  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Input_Directory_BIAS_KERNEL_TWO(){return this->Input_Directory_BIAS_KERNEL_TWO;}
 //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _Number_of_references(){return this->Number_of_references;}
 //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  
  int _Number_of_new_mocks(){return this->Number_of_new_mocks;}

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
  string _Input_Directory_X_REF_TWO(){return this->Input_Directory_X_REF_TWO;}
  //////////////////////////////////////////////////////////
   /**
    *  @brief get the value of the private member
    *  @return
    */
   string _Input_Directory_X_NEW(){return this->Input_Directory_X_NEW;}
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
  string _Name_Catalog_X_NEW_TWO(){return this->Name_Catalog_X_NEW_TWO;}
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
//////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
 
  int _NMASSbins_mf(){return this->NMASSbins_mf;}
  void set_NMASSbins_mf(int new_NMASSbins_mf){this->NMASSbins_mf=new_NMASSbins_mf;}
//////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
 
  int _NMASSbins_power(){return this->NMASSbins_power;}

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
  void set_lambdath(real_prec new_lambdath){this->lambdath=new_lambdath;}

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
 //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
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
  void set_N_iterations_Kernel(int new_N_iterations_Kernel){this->N_iterations_Kernel=new_N_iterations_Kernel;}
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
   bool _Apply_Rankordering_ab_initio(){return this->Apply_Rankordering_ab_initio;}
 //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  vector<int> _cwt_used(){return this->cwt_used;}
 //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
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
   *  @brief get the value of the private member Input_dir_cat
   *  @return Input_dir_cat
   */
  string _Input_dir_cat () {return this->Input_dir_cat;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member Input_dir_cat
   *  @return Input_dir_cat
   */
  string _Input_dir_cat_TWO () {return this->Input_dir_cat_TWO;}
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
  void set_dir_output(string new_dir_output){this->dir_output=new_dir_output;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member dir_output
   *  @return dir_output
   */
  string _input_type () {return input_type;}
  void set_input_type(string new_input_type){input_type=new_input_type;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member dir_output
   *  @return dir_output
   */

  string _input_type_two () {return input_type_two;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  ULONG _ngal_delta() {return ngal_delta;}
  void set_ngal_delta(ULONG new_ngal_delta) {this->ngal_delta=new_ngal_delta;}
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
   *  @brief set the value of the private member Input_dir_cat
   *  @param _Input_dir_cat input directory where the object and random
   *  catalogues are stored
   *  @return none
   */
  void set_Input_dir_cat (string _Input_dir_cat) {Input_dir_cat = _Input_dir_cat;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member file_catalogue
   *  @return file_catalogue
   */
  string _file_catalogue () {return file_catalogue;}
  void set_file_catalogue (string new_file_catalogue) {file_catalogue=new_file_catalogue;}

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
  int _i_coord1_g () {return this->i_coord1_g;}
  void set_i_coord1_g(int new_i_coord1_g){this->i_coord1_g=new_i_coord1_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord2_g
   *  @return i_coord2_g
   */
  int _i_coord2_g () {return i_coord2_g;}
  void set_i_coord2_g(int new_i_coord2_g){this->i_coord2_g=new_i_coord2_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_g
   *  @return i_coord3_g
   */
  int _i_coord3_g () {return this->i_coord3_g;}
  void set_i_coord3_g(int new_i_coord3_g){this->i_coord3_g=new_i_coord3_g;}
    //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord1_g
   *  @return i_coord1_g
   */
  int _i_v1_g () {return i_v1_g;}
  void set_i_v1_g(int new_i_v1_g){this->i_v1_g=new_i_v1_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord2_g
   *  @return i_coord2_g
   */
  int _i_v2_g () {return i_v2_g;}
  void set_i_v2_g(int new_i_v2_g){this->i_v2_g=new_i_v2_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_g
   *  @return i_coord3_g
   */
  int _i_v3_g () {return this->i_v3_g;}
  void set_i_v3_g(int new_i_v3_g){this->i_v3_g=new_i_v3_g;}


  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_g
   *  @return i_coord3_g
   */
  int _i_mass_g () {return this->i_mass_g;}
  void set_i_mass_g (int new_i_mass_g) {this->i_mass_g=new_i_mass_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _i_vmax_g () {return this->i_vmax_g;}
  void set_i_vmax_g (int new_i_vmax_g) {this->i_vmax_g=new_i_vmax_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _i_sf_g () {return this->i_sf_g;}
  void set_i_sf_g (int new_i_sf_g) {this->i_sf_g=new_i_sf_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _i_rs_g () {return this->i_rs_g;}
  void set_i_rs_g (int new_i_rs_g) {this->i_rs_g=new_i_rs_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _i_virial_g () {return this->i_virial_g;}
  void set_virial_g (int new_i_virial_g) {this->i_virial_g=new_i_virial_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */

  int _i_spin_g () {return this->i_spin_g;}
  void set_i_spin_g (int new_i_spin_g) {this->i_spin_g=new_i_spin_g;}
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

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
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
  void set_Name_survey(string new_Name_survey){this->Name_survey=new_Name_survey;}
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
  void set_use_random_catalog (bool new_use_random_catalog) {this->use_random_catalog = new_use_random_catalog ;}


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
   *  @brief get the value of the private member i_property5_g
   *  @return i_property5_g
   */
  int _i_mass_r () {return i_mass_r;}

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
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight1_r
   *  @return use_weight1_r
   */

  bool _weight_with_mass () {return this->weight_with_mass;}
  void set_weight_with_mass(bool new_weight_with_mass){this->weight_with_mass=new_weight_with_mass;}
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
   *  @brief get / set the value of the private member Nft
   *  @return Nft
   */
  ULONG _Nft () {return this->Nft;}

  void set_Nft(ULONG new_Nft){
    this->Nft=new_Nft;
    this->NGRID=new_Nft*new_Nft*new_Nft;
    this->NGRID_h=new_Nft*new_Nft*(new_Nft/2+1);
}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG _Nft_low () {return this->Nft_low;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG _Nft_HR () {return this->Nft_HR;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG _Nft_random_collapse () {return this->Nft_random_collapse;}
  void set_Nft_random_collapse (int new_Nft_random_collapse) {this->Nft_random_collapse=new_Nft_random_collapse;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _Distance_fraction() {return this->Distance_fraction;}
  void set_Distance_fraction(real_prec new_Distance_fraction){this->Distance_fraction=new_Distance_fraction;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member Nft
   *  @return Nft
   */
  ULONG _NGRID(){return this->NGRID;}
  void set_NGRID(ULONG new_NGRID){this->NGRID=new_NGRID;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member Nft
   *  @return Nft
   */
  ULONG _NGRID_h(){return this->NGRID_h;}
  void set_NGRID_h(ULONG new_NGRID_h){this->NGRID_h=new_NGRID_h;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member Lbox
   *  @return Lbox
   */
  real_prec _Lbox () {return this->Lbox;}
  void set_Lbox(real_prec new_Lbox){this->Lbox=new_Lbox;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _Lbox_low () {return this->Lbox_low;}
  void set_Lbox_low(real_prec new_Lbox_low){this->Lbox_low=new_Lbox_low;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member mass_assignment_scheme
   *  @return mass_assignment_scheme
   */
  string _mass_assignment_scheme () {return this->mass_assignment_scheme;}
  void set_mass_assignment_scheme (string new_mass_assignment_scheme) {this->mass_assignment_scheme=new_mass_assignment_scheme;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member type_of_binning
   *  @return type_of_binning
   */
  string _type_of_binning () {return this->type_of_binning;}
  void set_type_of_binning(string new_type_of_binning){this->type_of_binning=new_type_of_binning;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member type_of_binning
   *  @return type_of_binning
   */

  real_prec _vkernel_exponent(){return this->vkernel_exponent;}

  //////////////////////////////////////////////////////////
  /**
   * @brief Identifies dark matter "DM" or tracer "TR"
   */
  void set_type_of_object(string new_type_of_object){this->type_of_object=new_type_of_object;}


  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member type_of_binning
   *  @return type_of_binning
   */
  string _extra_info () {return this->extra_info;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member N_log_bins
   *  @return N_log_bins
   */
  int _N_log_bins () {return this->N_log_bins;}
  void set_N_log_bins(int new_N_log_bins){this->N_log_bins=new_N_log_bins;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member ndel_data
   *  @return ndel_data
   */
  int _ndel_data () {return this->ndel_data;}
  void set_ndel_data(int new_ndel_data){this->ndel_data=new_ndel_data;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member ndel_window
   *  @return ndel_window
   */
  int _ndel_window () {return this->ndel_window;}
  void set_ndel_window(int new_ndel_window){this->ndel_data=new_ndel_window;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member N_mu_bins
   *  @return N_mu_bins
   */
  int _N_mu_bins () {return this->N_mu_bins;}
  void set_N_mu_bins(int new_N_mu_bins){this->N_mu_bins=new_N_mu_bins;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member MAS_correction
   *  @return MAS_correction
   */
  bool _MAS_correction () {return this->MAS_correction;}
  void set_MAS_correction(bool new_MAS_correction){this->MAS_correction=new_MAS_correction;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member FKP_weight
   *  @return FKP_weight
   */
  bool _FKP_weight () {return FKP_weight;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get/set the value of the private member SN_correction
   *  @return SN_correction
   */
  bool _SN_correction () {return SN_correction;}
  void set_SN_correction (bool new_sn_correction) {this->SN_correction=new_sn_correction;}
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
  string _file_power () {return this->file_power;}
  void set_file_power(string new_file_power){this->file_power=new_file_power;}
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
  void set_kmin_bk(real_prec new_kmin_bk){kmin_bk=new_kmin_bk;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member kmin_bk;
   *  @return kmax_bk
   */
  real_prec _kmax_bk () {return kmax_bk;}
  void set_kmax_bk(real_prec new_kmax_bk){kmax_bk=new_kmax_bk;}

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
  void set_kmax_y_ds (real_prec new_kmax_y_ds ){kmax_y_ds=new_kmax_y_ds;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _DeltaKmin () {return this->DeltaKmin;}
  void set_DeltaKmin(real_prec new_DeltaKmin){this->DeltaKmin=new_DeltaKmin;}
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
  /**
   *  @brief 
   *  @return
   */

  string input_type_two;

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

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  bool _use_wiggles () { return use_wiggles; }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  bool _redshift_space_coords_g () { return this->redshift_space_coords_g; }
  void  set_redshift_space_coords_g (bool new_redshift_space_coords_g) { this->redshift_space_coords_g=new_redshift_space_coords_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  string _vel_units_g () { return this->vel_units_g; }
  void set_vel_units_g (string new_vel_units_g) { this->vel_units_g=new_vel_units_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  bool _RR () { return RR; }
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  // Params useful for Patchy. OJO QUE EN PARAMS.cpp aun no estan leyendo muchos de estos
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _inputmode(){return this->inputmode;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _dataFileName(){return this->dataFileName;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _ic_power_file(){return this->ic_power_file;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _ic_WN_file(){return this->ic_WN_file;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _ic_file(){return this->ic_file;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _ic_input_type(){return this->ic_input_type;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _use_ic_file(){return this->use_ic_file;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _ic_alias_corrected(){return this->ic_alias_corrected;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _ic_WN_dir(){return this->ic_WN_dir;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _dir(){return this->dir;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _fastpmpath(){return this->fastpmpath;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _sfmodel(){return this->sfmodel;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _transf(){return this->transf;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _Normalize_initial_redshift(){return this->Normalize_initial_redshift;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _runsim(){return runsim;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _runv(){return runv;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _lognden(){return lognden;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _diffcosmorz(){return diffcosmorz;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _readPS(){return readPS;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _use_vel_kernel(){return this->use_vel_kernel;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _Nmean(){return Nmean;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _cs2(){return cs2;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _cdeltas2(){return cdeltas2;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _cpsi(){return cpsi;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _cst(){return cst;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _cs3(){return cs3;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _sfac(){return sfac;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _ep(){return ep;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _devpois(){return devpois;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _deltathH(){return deltathH;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG _N1(){return this->Nft;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG _Nchunk(){return this->Nchunk;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _L1(){return this->Lbox;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _biasL(){return this->biasL;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _biasE(){return this->biasE;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _biasrhoexp(){return this->biasrhoexp;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _biasepsilon(){return this->biasepsilon;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _biasrhoexp2(){return this->biasrhoexp2;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _biasepsilon2(){return this->biasepsilon2;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _biassign(){return this->biassign;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _biassign2(){return this->biassign2;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _biasone(){return this->biasone;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _velbias(){return this->velbias;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _velbias_dm(){return this->velbias_dm;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _velbias_random(){return this->velbias_random;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _xllc(){return this->xllc;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _yllc(){return this->yllc;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _zllc(){return this->zllc;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _xobs(){return this->xobs;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _yobs(){return this->yobs;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _zobs(){return this->zobs;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _seed(){return this->seed;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _seed_ref(){return this->seed_ref;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _masskernel(){return this->masskernel;}
  void set_masskernel(int new_masskernel){this->masskernel=new_masskernel;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _masskernel_vel(){return this->masskernel_vel;}
  void set_masskernel_vel(int new_masskernel_vel){this->masskernel_vel=new_masskernel_vel;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _slength(){return this->slength;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _slengthv(){return this->slengthv;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _vslength(){return this->vslength;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _deltath(){return this->deltath;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _file_bin_x_coord(){return this->file_bin_x_coord;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _file_bin_y_coord(){return this->file_bin_y_coord;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _file_bin_z_coord(){return this->file_bin_z_coord;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG _N_lines_binary(){return  this->N_lines_binary;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _get_distribution_min_separations(){return this->get_distribution_min_separations;} 
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _Number_of_chunks_new_dm(){return  this->Number_of_chunks_new_dm;}
 //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d1(){return this->d1;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d2(){return this->d2;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d3(){return this->d3;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d1_low(){return this->d1_low;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d2_low(){return this->d2_low;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d3_low(){return this->d3_low;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _M_exclusion(){return this->M_exclusion;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _LOGMASSmin(){return this->LOGMASSmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _LOGMASSmax(){return this->LOGMASSmax;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _VMAXmin(){return this->VMAXmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _VMAXmax(){return this->VMAXmax;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _RSmin(){return this->RSmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _RSmax(){return this->RSmax;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _SPINmin(){return this->SPINmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _SPINmax(){return this->SPINmax;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _MASS_units(){return this->MASS_units;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _Initial_Redshift_DELTA(){return this->Initial_Redshift_DELTA;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _Initial_Redshift_ic_power_file(){return this->Initial_Redshift_ic_power_file;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _Prop_threshold_multi_scale_1(){return this->Prop_threshold_multi_scale_1;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _Prop_threshold_multi_scale_2(){return this->Prop_threshold_multi_scale_2;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _Prop_threshold_multi_scale_3(){return this->Prop_threshold_multi_scale_3;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */  
  real_prec _Prop_threshold_multi_scale_4(){return this->Prop_threshold_multi_scale_4;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */  
  real_prec _Prop_threshold_multi_scale_5(){return this->Prop_threshold_multi_scale_5;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _Tolerance_factor_l1(){return this->Tolerance_factor_l1;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _Tolerance_factor_l2(){return this->Tolerance_factor_l2;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _Tolerance_factor_l3(){return this->Tolerance_factor_l3;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _Tolerance_factor_l4(){return this->Tolerance_factor_l4;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _Tolerance_factor_l5(){return this->Tolerance_factor_l5;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG _Nft_low_l1(){return this->Nft_low_l1;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG _Nft_low_l2(){return this->Nft_low_l2;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG _Nft_low_l3(){return this->Nft_low_l3;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG _Nft_low_l4(){return this->Nft_low_l4;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */

  int _IC_index(){return this->IC_index;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */

  string _files_new_dm_fields(int i){
    return this->files_new_dm_fields[i];
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _type_of_object(){
    return this->type_of_object;
  }

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */


  string _files_dm_references(int i){
    return this->files_dm_references[i];
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */

  string _files_tracer_references(int i){
    return this->files_tracer_references[i];
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _files_tracer_field_references(int i){
    return this->files_tracer_field_references[i];
  }

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _list_new_dm_fields(int i){
    return this->list_new_dm_fields[i];
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
    string _files_bias_references(int i){
    return this->files_bias_references[i];
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _files_kernel_references(int i){return this->files_kernel_references[i];}

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _Number_of_MultiLevels(){return this->Number_of_MultiLevels;}


  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG get_Nft_MultiLevels(int i){return this->list_Nft_MultiLevels[i];}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG get_Ntracers_MultiLevels(int i){return this->list_Ntracers_MultiLevels[i];}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG get_d_ml(int i){return this->Lbox/static_cast<real_prec>(this->list_Nft_MultiLevels[i]);}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
 real_prec get_PropThreshold_MultiLevels(int i){return this->list_Props_Threshold_MultiLevels[i];}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec get_Props_Tolerance_MultiLevels(int i){return this->list_Props_Tolerance_MultiLevels[i];}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  void set_Props_Tolerance_MultiLevels(int i, real_prec newTa){this->list_Props_Tolerance_MultiLevels[i]=newTa;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  void set_Ntracers_MultiLevels(int i, ULONG new_Nt){this->list_Ntracers_MultiLevels[i]=new_Nt;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  void show_params();


  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  // Derived paramerters
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
   real_prec _d_delta_x(){return this->delta_x;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
   real_prec _d_delta_y(){return this->delta_y;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
   real_prec _d_delta_z(){return this->delta_z;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
   real_prec _d_deltak_x(){return this->deltak_x;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
   real_prec _d_deltak_y(){return this->deltak_y;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
   real_prec _d_deltak_z(){return this->deltak_z;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
   real_prec _d_deltak_0(){return this->deltak_0;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
   real_prec _d_kmin(){return this->kmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
   real_prec _d_kmax(){return this->kmax;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
   real_prec _d_DeltaK_data(){return this->DeltaK_data;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
   real_prec _d_DeltaK_window(){return this->DeltaK_window;}
  
  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
   ULONG _d_Nnp_window(){return this->Nnp_window;}
  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
   ULONG _d_Nnp_data(){return this->Nnp_data;}

  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
   ULONG _d_Deltal(){return this->Deltal;}

  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
   ULONG _d_Deltamu(){return this->Deltamu;}

  ///////////////////////////////////////////////////////// 

};



#endif
