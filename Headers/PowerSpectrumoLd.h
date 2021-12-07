 /**
 * @class <PowerSpectrumF>
 * @brief Header file for the class PowerSpectrum::
 * @file PowerSpectrumF.h
 * @author Andres Balaguera-Antolínez
 * @version 1.0
 * @date  2020
 * @details This file defines the interface of the class PowerSpectrum, used to
 *  obtain the measurements of Power spectrum (3D, Angular).
 */


#ifndef __POWERSPECTRUM__
#define __POWERSPECTRUM__

#include "Params.h"
#include "FileOutput.h"
#include "Catalog.h"
#include "DnDz.h" 
#include "CoordinateSystem.h" 
#include "FftwFunctions.h"
#include "ScreenOutput.h"
//#include "cosmo_parameters.h"

using namespace std;


class PowerSpectrumF {


 private :

  //////////////////////////////////////////////////////////
  /** 
   * @brief Object of type CosmologicalFunctions 
   */
  Cosmology c_Cf;
  
  //////////////////////////////////////////////////////////
  /** 
   * @brief Object of type FileManager 
   */
  FileOutput c_Fm;

  //////////////////////////////////////////////////////////
  /** 
   * @brief Object of type OpStatistics
   */
  DnDz c_Op;
  
  //////////////////////////////////////////////////////////
  /** 
   * @brief Object of type ScreenOutput
   */
  ScreenOutput So;
  
  //////////////////////////////////////////////////////////
  /**
   * @brief Object of type GalaxyOperations
   */
  CoordinateSystem c_gal;
  
  //////////////////////////////////////////////////////////
  /**
   *  @name private variables for power spectrum
   */
  bool measure_cross;
  //////////////////////////////////////////////////////////
  /**
   *  @name private variables for power spectrum
   */
  int measure_cross_from_1;
  //////////////////////////////////////////////////////////
  /**
   *  @name private variables for power spectrum
   */
  int measure_cross_from_2;

  //////////////////////////////////////////////////////////
  /// Vector containing the k identifying spherical shells for the estimates of power spectrum
  //////////////////////////////////////////////////////////
  /**
   *  @name private variables for power spectrum
   */
  vector<real_prec> kvector_data;
  
  //////////////////////////////////////////////////////////
  /**
   *  @name Vector containing the k identifying spherical shells for the estimates pf the window function
   */
  vector<real_prec> kvector_window;
  
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector containing the k identifying spherical shells for the estimates of power spectrum* in 2d cart 
   */
  vector<real_prec> kvector_data2d;
  
  //////////////////////////////////////////////////////////
  /// Vector containing the k identifying spherical shells for the estimates of power spectrum* in 2d polar coordinates 
  vector<real_prec> muvector;
  
  //////////////////////////////////////////////////////////
  /// Vector containing the estimates of the monpole in spherical shells
  vector<real_prec> pk0;
  
  //////////////////////////////////////////////////////////
  /// Vector containing the estimates of the quadrupole in spherical shells
  vector<real_prec> pk2;
  
  //////////////////////////////////////////////////////////
  /// Vector containing the estimates of the hexadecapole in spherical  shells
  vector<real_prec> pk4;
  
  //////////////////////////////////////////////////////////
  /// Vector containing the estimates of the power of the window function in spherical shells
  vector<real_prec> pk_w;
  
  //////////////////////////////////////////////////////////
  /// Number of modes for the power spectrum..
  vector<int> modes_g; 

  //////////////////////////////////////////////////////////
  /// Vector containing the the estimates of power spectrum* in 2d cart 
  vector < vector<real_prec> > pkk;

  //////////////////////////////////////////////////////////
  /// Vector containing the the estimates of power spectrum* in 2d polar 
  vector < vector<real_prec> > pmk;
    
  //////////////////////////////////////////////////////////
  /// Gaussian variance of the estimator a la fkp. Is the same for Yamamoto (ds)
  vector<real_prec> sigma_fkp;
  
  //////////////////////////////////////////////////////////
  /// Gaussian variance for the quadrupole of Yamamoto ds
  vector<real_prec> sigma_y_l2;
  
  //////////////////////////////////////////////////////////
  /// Gaussian variance for hexadecapole of Yamamoto ds
  vector<real_prec> sigma_y_l4;
  
  //////////////////////////////////////////////////////////
  /// Vector of wavenumbers for Bispectrum 
  vector<real_prec> kvector_data_b;
  
  //////////////////////////////////////////////////////////
  /// Vector allocating estimates of Bispectrum*/
  vector<real_prec> bispectrum;


  //////////////////////////////////////////////////////////
  /// Vector allocating estimates of Bispectrum*/
  vector<real_prec> sn_bispectrum;

  
  //////////////////////////////////////////////////////////
  /// Vector allocating number of triangles */
  vector<int> modes_tri;
  //////////////////////////////////////////////////////////
  /// Output file for the angular power spectrum

  string input_file_mask;  

  //////////////////////////////////////////////////////////
  /// Column in the mask file where the number of the pixel is written
  int i_mask_pixel;

  //////////////////////////////////////////////////////////
  /// Column in the mask file where the RA of the pixel is written
  int i_mask_alpha;

  //////////////////////////////////////////////////////////
  /// Column in the mask file where the Dec of the pixel is written
  int i_mask_delta;

  //////////////////////////////////////////////////////////
  /// Column in the mask file where the flag (0/1) defining the mask is written
  int i_mask_flag;


  //////////////////////////////////////////////////////////
  /// Auxiliary string for output files.
  string rest_file;

 
  ///@}



  /**
   *  @name private variables for the 3d power spectrum
   */
  ///@{
  
  //////////////////////////////////////////////////////////
  /// Healpix related parameter used when an estimate of the mean number density on a non-constant depth survye is computed.
  long nside;
  
  //////////////////////////////////////////////////////////
  /// name of the survey
  string Name_survey;
  

  
  //////////////////////////////////////////////////////////
  /// Number of random objects in catalogue
  long N_random;

  //////////////////////////////////////////////////////////
  /// Statistics to be measured. Specified in the input parameter file.
  string statistics;
  //////////////////////////////////////////////////////////
  /// name of the object catalogue
  string input_type;

  //////////////////////////////////////////////////////////
  /// name of the object catalogue
  unsigned long ngal_delta;

  //////////////////////////////////////////////////////////
  /// name of the object catalogue
  string delta_grid_file2;
  //////////////////////////////////////////////////////////
  /// name of the object catalogue
  string delta_grid_file3;
  //////////////////////////////////////////////////////////
  /// name of the object catalogue
  string delta_grid_file4;

  //////////////////////////////////////////////////////////
  /// name of the random catalogue
  string file_random; 

  //////////////////////////////////////////////////////////
  /// name of output file for the dNdz 
  string file_dndz;   

  //////////////////////////////////////////////////////////
  /// name of output file for FKP power spectrum
  string file_power;  


  //////////////////////////////////////////////////////////
  /// name of output file for FKP power spectrum
  string file_MCF;  

  
   //////////////////////////////////////////////////////////
  /// name of log_file for the power spectrum
  string file_power_log;  

  //////////////////////////////////////////////////////////
  /// name of output file for 2d FKP power spectrum in cart. coordinates
  string file_power2d; 

  //////////////////////////////////////////////////////////
  /// name of output file for 2d FKP power spectrum in polar coordinates
  string file_power2d_mk;

  //////////////////////////////////////////////////////////
  /// name of output file for FKP power spectrum of the window function
  string file_window; 

  //////////////////////////////////////////////////////////
  /// name of output file for FKP bispectrum
  string file_bispectrum;

  //////////////////////////////////////////////////////////
   /// true &rarr; sample with constant depth. false &rarr; sample with varyiong depth. Demands pixelization in order to get estimates of nbar.
  bool constant_depth;   
  
  //////////////////////////////////////////////////////////
  /// Ctrue &rarr; ompute FKP error bars false
  bool FKP_error_bars;
  
  //////////////////////////////////////////////////////////
  /// Compute error bars following FKP exact formula(yes/no)
  /// If this is no, and the previous is yes
  /// the code uses the Veff approximation for the variance.
  bool FKP_error_bars_exact;

  //////////////////////////////////////////////////////////
  ///  true &rarr; use FKP weights; false &rarr; set weights to 1
  bool FKP_weight; 

  //////////////////////////////////////////////////////////
  /// Estimated power for FKP weights
  real_prec Pest;

  //////////////////////////////////////////////////////////
  ///  Lenght of the Fourier box in configuration space (in Mpc/h)
  real_prec Lbox;
  
  //////////////////////////////////////////////////////////
  /// Maximum wavenumber for the Yamamoto-Blake direct sum estimation.
  real_prec kmax_y_ds;


  //////////////////////////////////////////////////////////
  bool weight_with_mass;


  //////////////////////////////////////////////////////////
  /// Correct for the MAS: yes/no
  bool MAS_correction;

  //////////////////////////////////////////////////////////
  /// Select MAS: NGP=nearest grid point. CIC=cloud in cell. TSC= triangular shape cloud
  string mass_assignment_scheme; 

  //////////////////////////////////////////////////////////
  /// Estimate of the mean number density
  real_prec mean_density;

  //////////////////////////////////////////////////////////
  /// Number of mu-bins for P(k,mu)
  int N_mu_bins;
  
  //////////////////////////////////////////////////////////        
  /// Number of grid cells for the DFT mesh (per dimension)
  int Nft;           

  //////////////////////////////////////////////////////////
  /// Number of log-spaced bins in Fourier space
  int N_log_bins;           

  //////////////////////////////////////////////////////////
  ///  true &rarr; mean number density tabulated; false &rarr; compute it
  bool nbar_tabulated;


  //////////////////////////////////////////////////////////
  /// Ratio between the shell-width and the fundamental mode for window
  int ndel_window;   

  //////////////////////////////////////////////////////////
  /// true &rarr; Redefine a new line of sight; false &rarr; rotate sample such that the z-direction points towards the baricenter. 
  bool new_los;

  //////////////////////////////////////////////////////////
  /// Area of the survey.
  real_prec area_survey;


  //////////////////////////////////////////////////////////
  /// Type of k-binning 
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
  /// true &rarr; compute a new Lbox; false &rarr; the code uses Lbox
  bool new_Lbox;

  //////////////////////////////////////////////////////////
  /// the power spectrum is computed with the help of a random catalogue (if use_random_catalog = true), or from a simulation with known mean (if use_random_catalog = false) 
  bool use_random_catalog;



  //////////////////////////////////////////////////////////
  /// the power spectrum is computed with the help of a random catalogue (if use_random_catalog = true), or from a simulation with known mean (if use_random_catalog = false) 
  bool use_random_catalog_cl;



  /// coordinate system of the object catalogue
  /** options: 
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


  //////////////////////////////////////////////////////////
  /// the column where the first object coordinate (according to the system of coordinates of the catalog) is written
  int i_coord1_g;  

  //////////////////////////////////////////////////////////  
/// the column where the second object coordinate (according to the system of coordinates of the catalog) is written
  int i_coord2_g;  

  //////////////////////////////////////////////////////////
  /// the column where the third object coordinate (according to the system of coordinates of the catalog) is written
  int i_coord3_g;

  int i_mass_g;

  //////////////////////////////////////////////////////////
  /// the column where the object mean number density is written
  int i_mean_density_g;

  //////////////////////////////////////////////////////////
  /// the column where the first object weight is written
  int i_weight1_g;  

  //////////////////////////////////////////////////////////
  /// the column where the second object weight is written
  int i_weight2_g;  

  //////////////////////////////////////////////////////////
  /// the column where the third object weight is written
  int i_weight3_g;

  //////////////////////////////////////////////////////////
  /// the column where the fourth object weight is written
  int i_weight4_g;

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
  /// units of angles in the object catalogue: D &rarr; degrees; R &rarr; radians
  string angles_units_g;    



  int Nbins_r;
  real_prec rmin;
  real_prec rmax;
  string r_bin_type;
  
  //////////////////////////////////////////////////////////
  /// coordinate system of the random catalogue
  /** options: 
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
  /// the column where the random mean number density is written
  int i_mean_density_r; 

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
  /// units of angles in the random catalogue: D &rarr; degrees; R &rarr; radians
  string angles_units_r;

  //////////////////////////////////////////////////////////
  ///Number of redshift bins to measure dNdz
  int N_z_bins;

  //////////////////////////////////////////////////////////
  /// Minimum redshift of the sample. Used when dN dz ahs to be measured 
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
  /// Resolution Healpix for pixelization. Used when no nbar is tabulated and dNdz is to be computed from a non-constant depth sample
  int Healpix_resolution;

  ///@}
  
 /**
 *  @name bispectrum parameters
   */
  ///@{

  //////////////////////////////////////////////////////////
  /// These parameters is used to define the shells in k-space
  bool use_fundamental_mode_as_kmin_bk;

  
  /////////////////////////////////////////////////////////////
  /// Minimum k-value for constructing k-bins
  real_prec kmin_bk;

  /////////////////////////////////////////////////////////////
  /// Maximum k-value for constructing k-bins
  real_prec kmax_bk;

  //////////////////////////////////////////////////////////


 /**
 *  @name cosmological parameters
   */
  ///@{
  //////////////////////////////////////////////////////////
  /// &Omega;<SUB>M</SUB>: the density of baryons, cold dark matter and massive neutrinos (in units of the critical density) at z=0
  real_prec om_matter;

  //////////////////////////////////////////////////////////
  /// &Omega;<SUB>rad</SUB>: the radiation density at z=0
  real_prec om_radiation;

  //////////////////////////////////////////////////////////
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
  /// &sigma;<SUB>8</SUB>: the normalization of power spectrum
  real_prec sigma8;

  //////////////////////////////////////////////////////////
  /// T<SUB>CMB</SUB>: the present day CMB temperature [K]
  real_prec Tcmb;
  ///@}

  
  /* //////////////////////////////////////////////////////////   */
  /* // Real catalogue */
  /* Catalogue * real; */
  
  /* ////////////////////////////////////////////////////////// */
  /* // Synthetic catalogue */
  /* Catalogue * random; */




  
 public :

  //////////////////////////////////////////////////////////
 /**
   *  @brief Default constructor
   *  @return object of PowerSpectrum
   */
  PowerSpectrumF ():measure_cross(false),measure_cross_from_1(false),measure_cross_from_2(false),i_mask_alpha(0),i_mask_delta(0),i_mask_flag(0),i_mask_pixel(0),nside(0),N_random(0),Pest(0),Lbox(0),kmax_bk(0),kmax_y_ds(0),weight_with_mass(false),MAS_correction(false),mean_density(0),Nft(0),N_log_bins(0),ndel_data(0),ndel_window(0),ngal_delta(0),new_Lbox(0),new_los(false),area_survey(0),i_coord1_g(0),i_coord1_r(0),i_coord2_g(0),i_coord2_r(0),i_coord3_r(0),i_coord3_g(0),i_mean_density_g(0),i_mean_density_r(0),i_weight1_g(0),i_weight1_r(0),i_weight2_g(0),i_weight2_r(0),i_weight3_g(0),i_weight3_r(0)
,i_weight4_g(0),i_weight4_r(0),i_mass_g(0),i_property1_g(0),i_property1_r(0),i_property2_g(0),i_property2_r(0),i_property3_g(0),i_property3_r(0),i_property4_g(0),i_property4_r(0),i_property5_g(0),i_property5_r(0),Nbins_r(0),rmin(0),rmax(0){
      time_t time_bam;
      time(&time_bam);
      this->So.initial_time=time_bam;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief Constructor 
   *  @param params 
   */
  PowerSpectrumF (Params &_params): measure_cross(false),measure_cross_from_1(false),measure_cross_from_2(false),i_mask_alpha(0),i_mask_delta(0),i_mask_flag(0),i_mask_pixel(0),nside(0),N_random(0),Pest(0),Lbox(0),kmax_bk(0),kmax_y_ds(0),weight_with_mass(false),MAS_correction(false),mean_density(0),Nft(0),N_log_bins(0),ndel_data(0),ndel_window(0),ngal_delta(0),new_Lbox(0),new_los(false),area_survey(0),i_coord1_g(0),i_coord1_r(0),i_coord2_g(0),i_coord2_r(0),i_coord3_r(0),i_coord3_g(0),i_mean_density_g(0),i_mean_density_r(0),i_weight1_g(0),i_weight1_r(0),i_weight2_g(0),i_weight2_r(0),i_weight3_g(0),i_weight3_r(0)
  ,i_weight4_g(0),i_weight4_r(0),i_mass_g(0),i_property1_g(0),i_property1_r(0),i_property2_g(0),i_property2_r(0),i_property3_g(0),i_property3_r(0),i_property4_g(0),i_property4_r(0),i_property5_g(0),i_property5_r(0),Nbins_r(0),rmin(0),rmax(0)
  {
     this->params=_params;
     this->set_parameters_power (params);
     this->set_cosmological_parameters (params);
     time_t time_bam;
     time(&time_bam);
     this->So.initial_time=time_bam;
  }



  //////////////////////////////////////////////////////////
  /**
   *  @brief Constructor
   *  @param params
   */
  string delta_grid_file;

  
  //////////////////////////////////////////////////////////
  /**
   *  @brief Default destructor
   *  @return none
   */
  ~PowerSpectrumF () {}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief add the random and galaxy catalogs
   */
  /* void add_catalogues (Catalogue &_real, Catalogue &_random) */
  /* {  */
  /*   real = new Catalogue;  */
  /*   *real = _real;  */
  /*   random = new Catalogue;  */
  /*   *random = _random;  */
  /* } */
  //////////////////////////////////////////////////////////
#if defined (_USE_MASS_CUTS_PK_)  || defined (_USE_ALL_PK_)
  void add_catalogues(real_prec mcuts);
#elif defined (_USE_MASS_BINS_PK_)
  void add_catalogues(real_prec m_min, real_prec m_max);
#endif
  //////////////////////////////////////////////////////////
  /// Vector containing the galaxy catalog
  vector< real_prec > galaxy_catalog;

  //////////////////////////////////////////////////////////
  /// Vector containing the random catalog
  vector< real_prec > random_catalog;
  /// the number of real_precs for each random poin (to calculate displaements)
  int rc_n_columns;

  /// the number of real_precs for each galaxy (to calculate displaements)
  int gc_n_columns;


  //////////////////////////////////////////////////////////
  /// name of the object catalogue
  string file_data;


  //////////////////////////////////////////////////////////
  /// Ratio between the shell-width and the fundamental mode
  int ndel_data;     

  
  //////////////////////////////////////////////////////////
  /// Number of galaxies in catalogue
  long N_galaxy;



  //////////////////////////////////////////////////////////
  /// true &rarr; use Poisson shot-noise correction; false &rarr; set shot-noise to zero
  bool SN_correction; 

  real_prec shot_noise;

  //////////////////////////////////////////////////////////
  int get_catalog_index(int idx, int catalog)
  {
    if(catalog == 0)
      // return index of element idx in galaxy catalog 
      return gc_n_columns * idx;
    else
      // return index of element idx in random catalog 
      return rc_n_columns * idx;
  }

  real_prec* get_catalog_pointer(int idx, int catalog)
  {
    if(catalog == 0)
      return &this->galaxy_catalog[gc_n_columns * idx];
    else
      return &this->random_catalog[rc_n_columns * idx];
  }
  
  //////////////////////////////////////////////////////////
  /**
   * @brief Compute the Power_Spectrum and/or the Bispectrum
   * @details This function generates the estiametes of power spectrum ()
   *  (FKP, Yamamoto and their multipole decomposition) and the Bispectrum
   *  (using FKP).
   * @arg verbose use dto write on screen
*/
#ifdef _USE_POWER_IN_BAM_
  void compute_power_spectrum (bool verbose, vector<s_Halo>&tracer);
#else
  void compute_power_spectrum (bool verbose, bool mcuts);
#endif

  void compute_power_spectrum_grid();
  //////////////////////////////////////////////////////////
  /**
     
   */
  void compute_marked_correlation_function ();


  //////////////////////////////////////////////////////////
  /**
   * @brief Compute the Power_Spectrum and/or the Bispectrum
   * @details This function generates the estiametes of power spectrum ()
   *  (FKP, Yamamoto and their multipole decomposition) and the Bispectrum
   *  (using FKP).
   */
  void compute_power_spectrum_grid(const vector<real_prec>&);

  //////////////////////////////////////////////////////////
  /**
   * @brief Compute the Marked Power_Spectrum
   */

  void compute_marked_power_spectrum_grid(const vector<real_prec>&, const vector<real_prec>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief Compute the cross power spectrum between fields X and Y
    * @param dm bool set to true if X is dark matter so no shot-noise is computed for that component
 */
  void compute_cross_power_spectrum_grid(bool dm, vector<real_prec>&X,vector<real_prec>&Y);
  void compute_cross_power_spectrum_grid(bool dm, string file_X,string file_Y);



  //////////////////////////////////////////////////////////
  /** 
   * @brief Write output of power spectrum and/or bispectrum 
   * @details computed in the function compute_power_spectrum()
   */
  void write_power_spectrum ();


  void write_power_and_modes();

  void write_power_and_modes(string);
  //////////////////////////////////////////////////////////
  /** 
   * @brief Write output of power spectrum and/or bispectrum 
   * @details computed in the function compute_power_spectrum()
   */
  void write_power_spectrum_grid(string);

  
  //////////////////////////////////////////////////////////
  /**
   * @brief Compute the angular power spectrum in bins of redshift
   * @details using the Peebles estimator and the Healpix package
   */
  void compute_cl_power_spectrum ();
  
  //////////////////////////////////////////////////////////
 /**
   * @brief Write output of angular power spectrum in redshift bins
   */
  void write_cl_power_spectrum ();
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief Set the values of the private variables 
   */
  
  void set_parameters_power (Params &params)
  {
    this->statistics=params._statistics();
    this->Nft = params._Nft();
    this->Lbox = params._Lbox();
    this->input_type=params._input_type();
    this->delta_grid_file=params._delta_grid_file();
    this->delta_grid_file2=params._delta_grid_file2();
    this->delta_grid_file3=params._delta_grid_file3();
    this->delta_grid_file4=params._delta_grid_file4();
    this->measure_cross=params._measure_cross();
    this->measure_cross_from_1=params._measure_cross_from_1();
    this->measure_cross_from_2=params._measure_cross_from_2();
    this->ngal_delta=params._ngal_delta();
    
    this->use_random_catalog = params._use_random_catalog();
    this->kmax_y_ds= params._kmax_y_ds();
    this->new_los = params._new_los();
    this->mass_assignment_scheme = params._mass_assignment_scheme();
    this->type_of_binning = params._type_of_binning();
    this->k_bin_step = params._k_bin_step();
    this->new_Lbox=params._new_Lbox();
    this->N_log_bins = params._N_log_bins();                    
    this->ndel_data = params._ndel_data();             
    this->ndel_window = params._ndel_window();           
    this->N_mu_bins = params._N_mu_bins();                   
    this->MAS_correction = params._MAS_correction();   
    this->FKP_weight = params._FKP_weight();            
    this->SN_correction = params._SN_correction();         
    this->FKP_error_bars = params._FKP_error_bars();        
    this->FKP_error_bars_exact = params._FKP_error_bars_exact();  
    this->Pest = params._Pest();                  
    this->nbar_tabulated = params._nbar_tabulated();       
    this->constant_depth = params._constant_depth();        
    this->sys_of_coord_g = params._sys_of_coord_g(); 
    this->i_coord1_g = params._i_coord1_g();     
    this->i_coord2_g = params._i_coord2_g();     
    this->i_coord3_g = params._i_coord3_g();
    this->i_mass_g = params._i_mass_g();     
    this->i_mean_density_g = params._i_mean_density_g();

    
    this->angles_units_g = params._angles_units_g();      
    
    this->i_property1_g = params._i_property1_g();
    this->i_property2_g = params._i_property2_g();
    this->i_property3_g = params._i_property3_g();
    this->i_property4_g = params._i_property4_g();
    this->i_property5_g = params._i_property5_g();
    this->i_property6_g = params._i_property6_g();
    
    this->i_property1_r = params._i_property1_r();
    this->i_property2_r = params._i_property2_r();
    this->i_property3_r = params._i_property3_r();
    this->i_property4_r = params._i_property4_r();
    this->i_property5_r = params._i_property5_r();
    this->i_property6_r = params._i_property6_r();
    
    this->sys_of_coord_r = params._sys_of_coord_r();  
    this->i_coord1_r = params._i_coord1_r();         
    this->i_coord2_r = params._i_coord2_r();          
    this->i_coord3_r = params._i_coord3_r();         
    
    this->i_weight1_r = params._i_weight1_r();         
    this->i_weight2_r = params._i_weight2_r();          
    this->i_weight3_r = params._i_weight3_r();         
    this->i_weight4_r = params._i_weight4_r();         
    this->use_weight1_r = params._use_weight1_r();         
    this->use_weight2_r = params._use_weight2_r();          
    this->use_weight3_r = params._use_weight3_r();         
    this->use_weight4_r = params._use_weight4_r();         
    
    this->i_weight1_g = params._i_weight1_g();         
    this->i_weight2_g = params._i_weight2_g();          
    this->i_weight3_g = params._i_weight3_g();         
    this->i_weight4_g = params._i_weight4_g();         
    this->use_weight1_g = params._use_weight1_g();         
    this->use_weight2_g = params._use_weight2_g();          
    this->use_weight3_g = params._use_weight3_g();         
    this->use_weight4_g = params._use_weight4_g();
    this->weight_with_mass = params._weight_with_mass();
    this->Name_survey=params._Name_survey();
    
    this->i_mean_density_r = params._i_mean_density_r();   
    this->angles_units_r = params._angles_units_r();     
    this->N_z_bins = params._N_z_bins();             
    this->redshift_min_sample = params._redshift_min_sample();
    this->redshift_max_sample = params._redshift_max_sample();
    this->N_dndz_bins = params._N_dndz_bins();        
    this->new_N_dndz_bins = params._new_N_dndz_bins();         
    this->area_survey = params._area_survey();        
    this->Healpix_resolution = params._Healpix_resolution();       
    this->set_output_filenames(params);
    this->kmin_bk=params._kmin_bk();
    this->kmax_bk=params._kmax_bk();
    this->use_fundamental_mode_as_kmin_bk=params._use_fundamental_mode_as_kmin_bk();
    
    this->input_file_mask=params._input_file_mask(); 
    this->i_mask_flag=params._i_mask_flag();
    this->i_mask_alpha=params._i_mask_alpha();
    this->i_mask_delta=params._i_mask_delta();
    this->i_mask_pixel=params._i_mask_pixel();
  }
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   *  @brief Set the values of the private variables that regard filenames - with many catalogues
   */
  void set_output_filenames (Params &params, int icatalogue)
  {

    string newdir, suffix;



    if (params._n_catalogues()>1)
      {
	stringstream ss; ss << icatalogue;
	newdir = params._dir_output()+"Catalogue"+ss.str()+"/";
	suffix = "."+ss.str();
      }
    else
      {
	newdir="";
	suffix="";
      }


#ifdef _USE_REDSHIFT_BINS_
        Name_survey+="_zmin_"+to_string(redshift_min_sample) +"_zmax_"+to_string(redshift_max_sample);
#endif


    if(statistics=="Pk_y_ds")
      {
        file_power = params._dir_output()+newdir+statistics+"_"+Name_survey+"_kmax"+to_string(static_cast<float>(kmax_y_ds))+"_"+params._file_power()+".txt";
      }
    else
      file_power = params._dir_output()+newdir+statistics+"_"+Name_survey+"_Nft"+to_string(Nft)+"_"+mass_assignment_scheme+"_"+params._file_power()+".txt";
      

    file_power_log      = params._dir_output()+newdir+statistics+"_"+Name_survey+"_Nft"+to_string(Nft)+"_"+mass_assignment_scheme+"_"+params._file_power_log()+".log";

    file_dndz           = params._dir_output()+newdir+"dndz_"+Name_survey+params._file_dndz()+".txt";
    file_power2d        = params._dir_output()+newdir+statistics+"_"+Name_survey+"_Nft"+to_string(Nft)+"_"+mass_assignment_scheme+"_"+params._file_power2d()+".txt";
    file_power2d_mk     = params._dir_output()+newdir+statistics+"_"+Name_survey+"_Nft"+to_string(Nft)+"_"+mass_assignment_scheme+"_"+params._file_power2d_mk()+".txt";
    file_window         = params._dir_output()+newdir+statistics+"_"+Name_survey+"_Nft"+to_string(Nft)+"_"+mass_assignment_scheme+"_"+params._file_window();
    file_bispectrum     = params._dir_output()+newdir+statistics+"_"+Name_survey+"_Nft"+to_string(Nft)+"_"+mass_assignment_scheme+"_"+params._file_bispectrum()+".txt";
    file_data = params._file_catalogue()+suffix;
    file_random = params._file_random()+suffix;
    
  }
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////
  /**
   *  @brief Set the values of the private variables that regard filenames - with one catalogue
   */
  void set_output_filenames (Params &params)
  {
    
    file_MCF     = params._dir_output()+statistics+"_"+Name_survey+".txt";

#ifdef _USE_REDSHIFT_BINS_
        Name_survey+="_zmin_"+to_string(redshift_min_sample) +"_zmax_"+to_string(redshift_max_sample);
#endif


#ifdef _REDSHIFT_SPACE_
    file_power  = params._dir_output()+statistics+"_"+Name_survey+"_Nft"+to_string(Nft)+"_"+mass_assignment_scheme+"_"+params._file_power()+"_RSS.txt";
#else
    file_power  = params._dir_output()+statistics+"_"+Name_survey+"_Nft"+to_string(Nft)+"_"+mass_assignment_scheme+"_"+params._file_power()+".txt";
#endif


    file_power_log      = params._dir_output()+statistics+"_"+Name_survey+"_Nft"+to_string(Nft)+"_"+mass_assignment_scheme+"_"+params._file_power_log()+".log";

    
    file_dndz           = params._dir_output()+"dndz_"+Name_survey+params._file_dndz()+".txt";
    file_power2d        = params._dir_output()+statistics+"_"+Name_survey+"_Nft"+to_string(Nft)+"_"+mass_assignment_scheme+"_"+params._file_power2d()+".txt";
    file_power2d_mk     = params._dir_output()+statistics+"_"+Name_survey+"_Nft"+to_string(Nft)+"_"+mass_assignment_scheme+"_"+params._file_power2d_mk()+".txt";
    file_window         = params._dir_output()+statistics+"_"+Name_survey+"_Nft"+to_string(Nft)+"_"+mass_assignment_scheme+"_"+params._file_window();
    file_bispectrum     = params._dir_output()+statistics+"_"+Name_survey+"_Nft"+to_string(Nft)+"_"+mass_assignment_scheme+"_"+params._file_bispectrum()+".txt";

    //file_power_fb       = file_power_cl; // ???
    file_data = params._file_catalogue();
    file_random = params._file_random();

  }


  
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   *  @brief Set the values of the cosmological parameter private variable
   */
  void set_cosmological_parameters (Params &params) 
  {



    om_matter = params._om_matter();
    om_radiation = params._om_radiation();
    om_baryons = params._om_baryons();
    om_vac = params._om_vac();
    om_k = params._om_k();
    Hubble = params._Hubble();
    hubble = params._hubble();
    spectral_index =  params._spectral_index();
    w_eos  = params._w_eos();
    N_eff = params._N_eff();
    sigma8 = params._sigma8();
    Tcmb = params._Tcmb();
 }


  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member 
   *  PowerSpectrum::kvector_data
   *  @param i index of the vector kvector_data
   *  @return PowerSpectrum::kvector_data[i]
   */
  real_prec _kvector_data (int i) { return kvector_data[i]; }

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the size of the private member
   *  PowerSpectrum::kvector_data
   *  @return PowerSpectrum::kvector_data.size()
   */
  int _kvector_data_size () { return kvector_data.size(); }

  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  PowerSpectrum::pk0
   *  @param i index of the vector pk0
   *  @return PowerSpectrum::pk0[i]
   */
  real_prec _pk0 (int i) { return pk0[i]; }

  real_prec _nmodes_k (int i) { return modes_g[i]; }
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member 
   *  PowerSpectrum::pk2
   *  @param i index of the vector pk2
   *  @return PowerSpectrum::pk2[i]
   */
  real_prec _pk2 (int i) { return pk2[i]; }

  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  PowerSpectrum::pk4
   *  @param i index of the vector pk4
   *  @return PowerSpectrum::pk4[i]
   */
  real_prec _pk4 (int i) { return pk4[i]; }

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  PowerSpectrum::pk4
   *  @param i index of the vector pk4
   *  @return PowerSpectrum::pk4[i]
   */
  real_prec _sigma_fkp (int i) { return sigma_fkp[i]; }
  //////////////////////////////////////////////////////////

  Params params;

  Catalog tracer_cat;

  Catalog random_cat;
  //////////////////////////////////////////////////////////

  real_prec var_prop;

};

#endif
