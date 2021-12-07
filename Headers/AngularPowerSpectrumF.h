#ifndef __ANGULAR_POWER_SPECTRUM_F__
#define __ANGULAR_POWER_SPECTRUM_F__



#define nbess 200
#define NUM_ITER_MAP2ALM 0   /* Number of iterations for the map2alm_iter HealPix routine*/
#define EPSILON_MK 1.0
// ######################################################################
// ######################################################################

# include <iostream>
# include <fstream>
# include <vector>
# include <algorithm>
# include <math.h>
# include <sstream>
# include <iomanip>
# include <string>
# include "NumericalMethods.h"
# include "ScreenOutput.h"
# include "CosmologicalFunctions.h"
# include "FileOutput.h"
# include "Galaxy.h"
using namespace std;
#define no_zs -999.0






// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// Note that when we do averages over the m-modes,
// we do the following:
// (Sum_m Jlm) / (2l+1) from -l to +l is splitted
// in Jl0/(2l+1)+ sum_(m=1) Jlm/(l+0.5)
// Since we put all in the same loop, we make theinf (m==) to divide things by 2
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


/**
 *  @class AngularPowerF.h
 *  @brief The class AngularPowerF_
 *
 */


class AngularPowerF{
  
 protected:

 private:
  //////////////////////////////////////////////////////////
  /**
  * @brief Number of columns in galaxy cat
  */
  int n_columns_gal;
  //////////////////////////////////////////////////////////
  /**
  * @brief Number of columns in random
  */
  int n_columns_ran;

  //////////////////////////////////////////////////////////
  /**
   * @brief Number of columns in the mask
   */
  int n_columns_mask;
  
  //////////////////////////////////////////////////////////
  /**
   * @brief Object of type pointing, used in HealPix operations
  */
  pointing point;
  //////////////////////////////////////////////////////////
  /**
   * @brief RMS of galaxies in pixels
   */
  healpix_real rms_nran;
  //////////////////////////////////////////////////////////
  /**
   * @brief Mean number of random objects in pixels
   */
  healpix_real mean_number_randoms_pix;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of rings, depends on nside
   */
  int Nrings;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of redshift bins, set default zero
   */
  int nzbins;
  //////////////////////////////////////////////////////////
  /**
   * @brief Index identifying redhsift bin
   */
  int IZ;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector containing the minimum redshifts in the redshift bins
   */
  vector<healpix_real> z_min;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector containing the maximum redshifts in the redshift bins
   */
  vector<healpix_real> z_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the Wigner 3J symbols, used
   * to compute the mixing matrix Rll in it-s approximated
   * version,
   *  in the VLS
   */

  vector<vector<vector<healpix_real> > > Wigner3J;
  //////////////////////////////////////////////////////////






public:
  
  //////////////////////////////////////////////////////////
  /**
   * @brief Constructor
   */
  AngularPowerF(){
     std::cout.precision(6);
     std::cout.setf(ios::showpoint);
     std::cout.setf(ios::scientific);
  }

  AngularPowerF(string par_file){
     read_pars(par_file);
     std::cout.precision(6);
     std::cout.setf(ios::showpoint);
     std::cout.setf(ios::scientific);
  }


  //////////////////////////////////////////////////////////
  /**
   * @brief Default destructor
   */
  ~AngularPowerF(){}



  //////////////////////////////////////////////////////////
  /**
  *  @brief Galaxy catalog
  */
  vector<real_prec>   prop;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Random catalog
  */
  vector<real_prec>  prop_r;



  //////////////////////////////////////////////////////////
  /**
  *  @brief Read input cats
  */

  void read_input_cats(string cat, string file);


  //////////////////////////////////////////////////////////
  /**
  *  @brief Get redshift bins
  */
  void get_zbins_same_ngal(int, int, healpix_real, healpix_real,vector< vector<healpix_real> >&);
  //////////////////////////////////////////////////////////

  /**
  *    @brief Inpput/Output
  */
  FileOutput Fmi;
  //////////////////////////////////////////////////////////
  /**
   * @brief Fmd
  */
  FileOutput Fmd;
  //////////////////////////////////////////////////////////
  /**
   * @brief Random generator
  */
  gsl_rng * r;
  //////////////////////////////////////////////////////////
  /**
  * @brief  Type of code running. Do not touch it.
  */
  string code;
  //////////////////////////////////////////////////////////
  /**
   * @brief Statistics to measure, basically Cl
  */
  string statistics;
  //////////////////////////////////////////////////////////
  /**
   * @brief Name of the input catalog from wich the
   * Statistics will be measured. Overwritten by param file
  */
  string name_catalog;
  //////////////////////////////////////////////////////////
  /**
   * @brief Output directory where the Cl will be written
   * Overwritten by param file
  */
  string name_output_dir;
  //////////////////////////////////////////////////////////
  /**
   * @brief Input directory where the input cat is located
   * Overwritten by param file
  */
  string name_input_dir;
  //////////////////////////////////////////////////////////
  /**
   * @brief Name of catalog
   * Overwritten by param file
  */
  string input_file;
  //////////////////////////////////////////////////////////
  /**
   * @brief Type of file for the galaxy catalog, options are
   * ascii and fits
   * Overwritten by param file
  */
  string file_type;

  //////////////////////////////////////////////////////////
  /**
   * @brief Type of file for the Healpix mask, options are
   * ascii and fits
   * Overwritten by param file
  */
  string file_type_mask;

  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_alpha;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_delta;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_z;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_M;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_w;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_alpha_ran;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_delta_ran;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_z_ran;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_M_ran;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_w_ran;
  //////////////////////////////////////////////////////////
  /**


  //////////////////////////////////////////////////////////
  /**
   * @brief  We can use a random catalog.
   * In that case the map built from this catalog will
   * play the role of mask.
   * Overwritten by param file
  */
  bool use_random_cat;

  //////////////////////////////////////////////////////////
  /**
   * @brief This parameters
   *  will help us to avoid the calculation of Jlm or Ilm
   *  Options are full_sky /  masked_sky
   *  (def masked_sky)
  */
  string sky;

  //////////////////////////////////////////////////////////
/**
 * @brief Choose the coordinate system between galactic and
 * equatorial. This is only useful in case we split between
 * south and  north. When hemis=all, by default, we use galactic
 *  coordinates according to te input mask.
 *  The code internaly selects north / south if galacitic
 *  but demands an input mask if we split in equatorial
*/
  string coord;

  //////////////////////////////////////////////////////////
  /**
   * @brief Input file for the mask in HealPix format
  */
  string input_file_mask;
  //////////////////////////////////////////////////////////
  /**
   * @brief Input file for the north equatorial mask
  */
  string input_file_mask_north_equatorial;
  //////////////////////////////////////////////////////////
  /**
   * @brief Input file for the southern equatorial mask
  */
  string input_file_mask_south_equatorial;
  //////////////////////////////////////////////////////////
  /**
   * @brief Input file for the fulls sky mask
  */
  string input_file_mask_fs;

  //////////////////////////////////////////////////////////
  /**
   * @brief no (yes) if you (do not) want to generate fits files
   * containing the overdensity map. Note that if these are already
   * created, the code will stop. Delete fits files first.
  */
  bool generate_fits_files;

  //////////////////////////////////////////////////////////
  /**
   * @brief Number of REDSHIFT BINS z_bins, used if selection=fls
  */
  int n_z_bins;
  //////////////////////////////////////////////////////////
    /** @brief Select the type of redshift bins. If the variable
     *  define_z_bins is set to "delta",
     *  the code takes zmin and zmax (either given here or
     * taken from the catalog,  as has been specified in the
     * variable min_max_from_cat)  and generate n_z_bins in
     * redshift with the same width. If set  to "number",
     * the code take the same redshift range and  generate
     * n_z_bins each containing the same number of galaxies
*/
  string define_z_bins;
  //////////////////////////////////////////////////////////
  /**
   * @brief Determines hether to use Healpix or direct sum over galaxies
  */
  string sampling;
  //////////////////////////////////////////////////////////
  /**
   * @brief The kind of Peeble -lile estimator to use
  */
  string type_of_P_estimator;
  //////////////////////////////////////////////////////////
  /**
   * @brief Type of selection. This is fixed
  */
  string selection;
  //////////////////////////////////////////////////////////
  /**
   * @brief Hemisphere
  */
  string hemis;
  //////////////////////////////////////////////////////////
  /**
   * @brief Define the type of redshift
  */
  string redshift;
  //////////////////////////////////////////////////////////
  /**
   * @brief Type of Magnitude related property, off
  */
  string property_i_M_type;
  //////////////////////////////////////////////////////////
  /**
   * @brief TYpe of redshift bins
  */
  string select_z_bins;
  //////////////////////////////////////////////////////////
  /**
   * @brief Name of random ascii file
  */
  string input_file_random;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  string output_file_raw ;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  string output_file;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  string output_file_Jlm;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  string shot_noise_correction;
  ///////////////////////////////////////////////////////////
  /**
   * @brief
  */
  string output_file_window;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  string bin_type;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  string fits_map;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  bool compute_jlm;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int Lmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int Lmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int nside;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int N_L_bins;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int ngal;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  bool compute_mixing_matrix;
  //////////////////////////////////////////////////////////

  /**
   * @brief
  */
  bool mixing_matrix_exact=false;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  bool use_weight;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  bool use_Mk_min_max_from_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  bool use_z_min_max_from_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  bool compute_property_weighted_cl;
  //////////////////////////////////////////////////////////
  /**
   * @brief
*/
  bool use_external_pk_file;
  //////////////////////////////////////////////////////////
   /**
    * @brief
  */
  int i_mask_pixel;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int i_mask_alpha;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int i_mask_delta;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int i_mask_flag;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int n_M_bins;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  healpix_real MKmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  healpix_real MKmax;

  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  healpix_real zmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  healpix_real zmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int nran;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  healpix_real alpha;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  healpix_real zmin_bin;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  healpix_real zmax_bin;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int number_of_realizations;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  healpix_real zmax_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  healpix_real zmin_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  healpix_real MKmax_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  healpix_real MKmin_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  healpix_real area_pixel;
  //////////////////////////////////////////////////////////
  /**
   * @brief Expected number of galaxies in one pixel
  */
  healpix_real mean_number_galaxies_pix;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  healpix_real mean_number_galaxies; //expected number of galaxies. Ngal/Area
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  healpix_real rms_ngal;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  healpix_real shot_noise;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  healpix_real area_survey;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  healpix_real sky_fraction;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<healpix_real> theta_new;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<healpix_real> phi_new;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<healpix_real> theta_new_pix;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<healpix_real> Mean_ngal_pix; // mean surface density of galaxies per pixel
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<healpix_real> Mean_ngal; //mean surface density of galaxies in sample
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<healpix_real> Shot_Noise;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<healpix_real> pixmask;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<vector<healpix_real> > Jlm;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector to the Alms in Magnitude or Redshift bins
   */
  vector<vector<vector<complex<healpix_real> > > > Blm;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector to allocate zmin and zmax of nbins in z.
   * The first 0-0 refers to the fill z-interval. Used only for
   * FLS z-interval. Used only for FLS
   */
  vector< vector<healpix_real> >zminmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector containing the Window function
  */
  vector<healpix_real> Wl;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector containing the l-modes
  */
  vector<int> lvec;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<healpix_real> Clvec;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<healpix_real> lbin;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<healpix_real> lbin_min;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<healpix_real> lbin_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<healpix_real> Clbin;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<healpix_real> Clbin_meas;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<healpix_real> eClbin;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of angular modes
  */
  vector<int> nmodes;
  //////////////////////////////////////////////////////////
  /**
   * @brief ANgular power of the mask in lbins
  */
  vector<healpix_real> Wlbin;
  //////////////////////////////////////////////////////////
  /**
   * @brief Mixing matrix Rll
  */
  vector<vector<healpix_real> > R;
  //////////////////////////////////////////////////////////
  /**
   * @brief MIxing matrix is bins, l-lbin
  */
  vector<vector<healpix_real> > Rll_bin;
  //////////////////////////////////////////////////////////
  /**
   * @brief NUmber of used galaxies, i.e, not masked
  */
  int ngal_used;
  //////////////////////////////////////////////////////////
  /**
   * @brief NUmber of used random objects, i.e, not masked
  */
  int nran_used;
  //////////////////////////////////////////////////////////
  /**
   * @brief NUmber of observed pixels from the mask
  */
  int n_observed_pixels;
  //////////////////////////////////////////////////////////
  /**
   * @brief NUmber of pixels in the mask
  */
  int n_total_pixels;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int n_pixels;
  //////////////////////////////////////////////////////////
  /**
   * @brief Get the parameters from input file
  */
  void read_pars(string &);

  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  void set_vectors();
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  void set_ilm_jlm();
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  void set_mean_ngal_pix(int);

  //////////////////////////////////////////////////////////
  /**
   * @brief   Index used to pindown the redshift label
  */
  void set_IZ(int);
  //////////////////////////////////////////////////////////
  /**
   * @brief   Define the z-bins
   */
  void set_zbins();
  //////////////////////////////////////////////////////////
  /**
   * @brief   Define the type of L-bins
  */
  void set_Lbins();
  //////////////////////////////////////////////////////////
  /**
   * @brief  Function relating Healpix index and ring
  */
  int npix_ring(int, int);
  //////////////////////////////////////////////////////////
  /**
   * @brief  COMpute and set Healpix related quantities from the mask
  */
  void set_healpix_pars();
  //////////////////////////////////////////////////////////
  /**
   * @brief  Read the mask from input file
  */
  void get_mask(string);
  //////////////////////////////////////////////////////////
  /**
   * @brief   Get the map from the cat
  */
  void get_map(char,string, Healpix_Map<healpix_real>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief  Array for the ALm in redshift bins
  */
  void set_BLMZ(int);
  //////////////////////////////////////////////////////////
  /**
   * @brief  Get parameters associated to the CL estimator
  */
  void get_pars_cl_estimator();
  //////////////////////////////////////////////////////////
  /**
   * @brief   Get parmeters asociated to the mask
  */
  void get_pars_mask();
  //////////////////////////////////////////////////////////
  /**
   * @brief Get the mixing matrix
  */
  void get_mixing_matrix();
  //////////////////////////////////////////////////////////
  /**
   * @brief Write the parameters in screen
  */
  void write_pars_cl_estimator();
  //////////////////////////////////////////////////////////
  /**
   * @brief COmputes the Wigner symbols
  */
  void W3J();
  //////////////////////////////////////////////////////////
  /**
   * @brief   Get the Map from the catalog
  */
  void Cat2Map(char, Healpix_Map<healpix_real>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief  Get the map from the catalog,epressed in terms of deltas
  */
  void Cat2Map_delta(char, Healpix_Map<healpix_real>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  void Cat2Map_noz(char, Healpix_Map<healpix_real>&);

   //////////////////////////////////////////////////////////
  /**
   * @brief  Get the ALms from the map, our way
  */
  void Map2Alm(Healpix_Map<healpix_real>, Alm<xcomplex <healpix_real> >&,  Alm<xcomplex <healpix_real> >&, arr<arr<healpix_real> >& );
  //////////////////////////////////////////////////////////
  /**
   * @brief   Get the Alm from the map
  */
  void Map2Alm(Healpix_Map<healpix_real>, Alm<xcomplex <healpix_real> >&);
  //////////////////////////////////////////////////////////
  /**
   * @brief  Get the map from the Alms , using Healpix
  */
  void Alm2Map(Alm<xcomplex <healpix_real> >&, Healpix_Map<healpix_real>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief  Get Alm from the random map
  */
  void Map2Alm_ran(Healpix_Map<healpix_real>,Healpix_Map<healpix_real>, Alm<xcomplex <healpix_real> >&,  Alm<xcomplex <healpix_real> >&);
  //////////////////////////////////////////////////////////
  /**
   * @brief   Get the function Jlm used n the D-like estimator
  */
  void Map2Jlm(Healpix_Map<healpix_real>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  void Map2Ilm_Jlm(Alm<xcomplex <healpix_real> >&);
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  void Alm2Cl(Alm<xcomplex <healpix_real> >&, vector<int>&, vector<healpix_real>&);
  /**
   * @brief   .
  */
  //////////////////////////////////////////////////////////
  /**
   * @brief   COnvert Alm to Cl when partial sky coverage is prensent
  */
  void Alm2Cl(string est, Alm<xcomplex <healpix_real> >&, Alm<xcomplex <healpix_real> >&, arr<arr<healpix_real> >&, vector<int>&, vector<healpix_real>&, vector<healpix_real>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief  COnvert Alm to Cl for full sky
  */
  void Alm2Cl(int, int, Alm<xcomplex <healpix_real> >&,  vector<healpix_real> &);  // the one is
  //////////////////////////////////////////////////////////
  /**
   * @brief   Get Gaussian errors on the estimateos of Cl
  */
  void get_eCl(int, int, vector<healpix_real>,vector<healpix_real>,vector<healpix_real>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief  Get mixing matrix
  */
  void Mixing_matrix();
  //////////////////////////////////////////////////////////
  /**
   * @brief ANgular power in bins
  */
  void Cl_bins(vector<healpix_real>, vector<healpix_real>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief Gaussian errors in bins
  */
  void eCl_bins(vector<healpix_real>, vector<healpix_real>&);
   //////////////////////////////////////////////////////////
  /**
   * @brief Get the Mixing matix in bins
  */
  void Rll_bins();
  //////////////////////////////////////////////////////////
  /**
   * @brief Main member function doing all the job
  */
  void get_my_power();
  //////////////////////////////////////////////////////////
  /**
  *  @brief  Obtain minimum value of the i-th column of gal properties
  *          The column ic is passed explicitely
  */
  healpix_real get_min(char, int ic);
  //////////////////////////////////////////////////////////
   /**
   *  @brief  Obtain maximum value of the i-th column of gal properties
   *          The column is passed explicitely
   */
  healpix_real get_max(char, int);
  //////////////////////////////////////////////////////////
  /**
   * @brief Main member function doing all the job
  */
  void get_Cl_gal();






};

#endif
