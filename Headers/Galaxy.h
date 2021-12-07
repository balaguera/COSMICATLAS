/**
 *  @file Galaxy.h
 *
 *  @brief The class Galaxy
 *
 *  This file defines the interface of the class Galaxy.  Note: Our
 */


#ifndef __GALAXY__
#define __GALAXY__


//#define _USE_HEALPIX_
//#define _USE_MASK_
//#define _USE_mK_
//#define _USE_COLOR_
//#define _USE_JK_
//#define _USE_LF_COLE_   //meant to speed calculations of Mags within iterative procedures


#define _USE_MASS_

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

using namespace std;

#define fac M_PI/180.0
#define no_zs -999.0


struct s_bins_info{
    int index;
    real_prec min;
    real_prec max;
    int Nbins;
    real_prec delta;
    string type;
    void get_delta()
    {
     this->delta=(this->max-this->min)/static_cast<real_prec>(this->Nbins);
    };
    void show(){
      cout<<"Index = "<<this->index<<endl;
      cout<<"Min = "<<this->min<<endl;
      cout<<"Max = "<<this->max<<endl;
      cout<<"Nbins = "<<this->Nbins<<endl;
      cout<<"Delta = "<<this->delta<<endl;
      cout<<"Type of binning = "<<this->type<<endl;
    }
};



// ######################################################################
// ######################################################################


// Allocation of the position of different galaxy properties
// in the input catalog
struct gal_parameters{
  int i_ra;
  int i_dec;
  int i_lgal;
  int i_bgal;
  int i_zs;  
  int i_zp;
  int i_mJ;  
  int i_mH;
  int i_mK;
  int i_EBV;
  int i_lsd;
  int i_Kcsb;
  int i_mask_flag;
  real_prec z_min;
  real_prec z_max;
  real_prec mK_min;
  real_prec mK_max;
  real_prec MK_min;
  real_prec MK_max;
  real_prec color_min;
  real_prec color_max;
  int N_Mag_bins;
  int N_mag_bins;
  int N_z_bins;
  int N_color_bins;
  int N_iterations;
  real_prec area_in_deg;
  string ztype;

};

// ######################################################################  
// ######################################################################  
// ######################################################################  
// ######################################################################  

class GALAXY{


    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////

private:

  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  pointing point;

  //////////////////////////////////////////////////////////
  /**
  *  @brief Structure containing the properties
  */

  ofstream vdina;

  //////////////////////////////////////////////////////////
  /**
  *  @brief Get the absolute magnitude from
  */
  void get_MK(real_prec, real_prec, real_prec &Mks);

  void get_MK(real_prec, real_prec, real_prec, real_prec, real_prec &Mks);

  //////////////////////////////////////////////////////////
  /**
  *  @brief Auxiliary integer used to count objects
  */
  ULONG count;
  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */

  real_prec zzns;
  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */
  string output_file_dndz;


  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */
  string output_file_dndmk;

  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */
  string output_file_dndz_smooth;
  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */

  string output_file_nbar_smooth;
  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */

  string output_file_dndz_lf;
  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */

  string output_file_LF;
  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */
  string output_file_NM;
  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */

  string output_file_new_cat;

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////




 public:


  //////////////////////////////////////////////////////////
  /**
  *  @brief Constructor
  */
  GALAXY(){}
  //////////////////////////////////////////////////////////
  /**
  *  @brief Destructor
  */
  ~GALAXY(){}
  //////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */
  string LF_estimator;

  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */

  bool Get_new_catalog;


  bool Redefine_MKminmax;

  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */
  bool Measure_LF;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector to allocate the comoving distance
  */
  vector<gsl_real>rv;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector to allocate the transverse comoving distance
  */
  vector<gsl_real>trv;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector to allocate the redshift
  */
  vector<gsl_real>zv;
  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */

  vector<gsl_real>zv_n;

  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector to allocate the growth factor normalized to unity at z=0
  */
  vector<gsl_real>gv;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector to allocate the galaxy bias as a function of z
  */
  vector<gsl_real>bias_zv;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector to allocate the growth index
  */
  vector<gsl_real>gfv;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector to allocate the Distance modulus
  */
  vector<gsl_real>Dm;
  //////////////////////////////////////////////////////////
  // Vectors for redshift distributions
  /**
  *  @brief Vector to allocate the redshift for redshift distributions
  */
  vector<gsl_real> zn;
  vector<gsl_real> zn_low_res;

  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector to allocate the min redshift in arrays for for redshift distributions
  */
  vector<gsl_real> zn_min;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector to allocate the max redshift in arrays for for redshift distributions
  */
  vector<gsl_real> zn_max;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector to allocate the Hubble function
  */
  vector<gsl_real>Hv;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector to allocate Mags, and hist positions of gals, used when Vmax_dc is implemented
  */

  vector<int> index_z_Gal;
  vector<int> index_Mk_Gal;
  vector<int> index_Mk_low_Gal;
  vector<int> Gal_inside_intervals;

  vector<int> Gal_inside_intervals_z_m_c; //to check if galaxies are in the z, mK and color intervals only once and then change 3 ifs for only 1

  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector to allocate the Hubble function
  */
  vector<gsl_real>v_z_Kcorr;
  vector<gsl_real>v_Kcorr;
  vector<gsl_real>v_color_correction_Kcorr;



  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector to allocate the Maximum volume as a function of the z and M
  */
  vector<gsl_real> Vmax_v;


  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector to allocate the Maximum Area as a function of  K-magnitude
  */
  vector<gsl_real> Amax_v;

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector to allocate the Maximum Area as a function of  K-magnitude
  */
  string type_of_K_correction;

  real_prec kcorr_index;

  /**
  *  @brief Vector to allocate Magnitudes used in the interpolation to get Vmax
  */
  vector<gsl_real> MK_n;

  //////////////////////////////////////////////////////////
  /**
  *  @brief Object to the Cosmological Functions class
  */
  Cosmology Cf;

  ScreenOutput So;

  FileOutput File;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Structure to allocate cosmological parameters
  */
  s_CosmologicalParameters scp;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Galaxy catalog
  */
  vector<real_prec> prop;

  vector<int > gal_mask;

  //////////////////////////////////////////////////////////
  /**
  *  @brief Random catalog
  */
  vector<real_prec>  prop_r;

  vector<int > ran_mask;

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector to allocate the volume of an spherical shell
  */
  vector<gsl_real> Vshell;
  vector<gsl_real> Vshell_lowres;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<real_prec> nbarq;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<vector<real_prec> >S;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector for the redshift distribution
  */
  vector<gsl_real> dNdz;

  vector<real_prec> dNdM;

  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector for the redshift distribution
  */
  vector<gsl_real> dNdmk;

  vector<gsl_real> v_mk;

  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector for the redshift distribution
  */
  vector<real_prec> dNdz_low_res;

  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<real_prec> edNdz; // Poisson error in the dndz
  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */
  vector<real_prec> edNdz_low_res; // Poisson error in the dndz
  //////////////////////////////////////////////////////////
  /**
  *  @brief dndz from the luminosity function
  */
  vector<real_prec> dNdz_lf;
  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */
  vector<real_prec> dNdz_lf_low_res;

  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<gsl_real> nbar;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<gsl_real> zn_new;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<gsl_real> dNdz_new;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<gsl_real> nbar_new;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<gsl_real> prob;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector for the Luminosity function
  */
  vector<real_prec > Pm;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector for the Luminosity function
  */
  vector<real_prec > Pm_low_res;

  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector for the distribution of colors
  */
  vector<real_prec > Pc;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector for the joint distribution of color-Magnitude
  */
  vector<vector<real_prec> >PcMk;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector for the joint distribution of color-Magnitude
  */

  void get_PcMk();


  //////////////////////////////////////////////////////////
  /**
  *  @brief In this function we compute the joint probability P(X,Y) fron which the conditional P(Y|X) is computed
  * to be uysed later in construction of, e.g., randoms.
  */

  void get_P_X_Y(s_bins_info *bx, s_bins_info *by,string type);


  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector for the joint distribution of color-Magnitude
  */
  void get_PcMk_Vmax();
  //////////////////////////////////////////////////////////
  /**
  *  @brief Vector for the joint distribution of color-Magnitude
  */

  void get_Pzc_from_PcMk();

  //////////////////////////////////////////////////////////
    /**
  *  @brief Vector for the joint distribution of redshift-Magnitude
  */
  vector<vector<real_prec> >PzM;


  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */
  void get_PzMK();

  //////////////////////////////////////////////////////////
    /**
  *  @brief Vector for the joint distribution of redshift-Magnitude
  */
  vector<vector<real_prec> >Pzm;


  //////////////////////////////////////////////////////////
    /**
  *  @brief Vector for the joint distribution of redshift-Magnitude
  */
  vector<vector<real_prec> >PKc;

  //////////////////////////////////////////////////////////
    /**
  *  @brief Vector for the joint distribution of redshift-Magnitude
  */
  vector<real_prec> PXY;
  vector<real_prec> PXY_lowres;
  vector<real_prec> NXY;
  vector<real_prec> nbarXY;
  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */

  void get_PzK();

  //////////////////////////////////////////////////////////
    /**
  *  @brief Vector for the joint distribution of redshift-Magnitude
  */
  vector<vector<real_prec> >Pzc;


  vector<vector<real_prec> >Pzc_from_PzMk;
  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */

  void get_Pzc(bool dot);


  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */

  void get_PKc();


  //////////////////////////////////cM////////////////////////
  /**
  *  @brief Vector for the number of galaxies in a bin of Magnitude
  */
  vector<real_prec> number_m; //CAMBIAR A int
  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */

  vector<real_prec> number_m_low_res;  //CAMBIAR A int

  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<real_prec> number_z;  //CAMBIAR A int
  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */

  vector<real_prec> number_z_low_res;//CAMBIAR A int

  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<real_prec> number_c; //CAMBIAR A int
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<vector<real_prec> >number_cm;  //CAMBIAR A int
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<real_prec>v_Magnitude;
  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */

  vector<real_prec>v_Magnitude_low_res;

  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<real_prec>v_color;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<real_prec> Delta_it;
  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */

  vector<real_prec> Delta_it_low_res;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec alpha;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  string name_input_file_cat;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  string name_file_Kcorr;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  string name_file_mask;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  string name_file_mask_north_gal;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  string name_file_mask_south_gal;

  //////////////////////////////////////////////////////////
  /**
  *  @brief Column with the apparent magnitude in the J band
  */
  int i_mJ;



  //////////////////////////////////////////////////////////
  /**
  *  @brief Column with the apparent magnitude in the J band
  */
  int i_Kcsb;


  //////////////////////////////////////////////////////////
  /**
  *  @brief Column with the apparent magnitude in the K band
  */
  int i_mK;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Column with the apparent magnitude in the H band
  */
  int i_mH;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Column with the right ascencion
  */
  int i_ra;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Column with the declination
  */
  int i_dec;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Column with
  */
  int i_lgal;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Column with
  */
  int i_mass;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Column with
  */
  int i_bgal;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Column with
  */
  int i_EBV;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Column with
  */
  int i_lsd;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Column with
  */
  int i_zs;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Column with
  */
  int i_zp;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Column with
  */
  int i_mask_pixel;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Column with
  */
  int i_lpix;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Column with
  */
  int i_bpix;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Column with
  */
  int i_mask_flag;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Healpix resolution
  */
  int Healpix_resolution;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Column with
  */
  int N_iterations;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  int i_z_Kcorr;
  int i_flux_factor_Kcorr;
  int i_color_correction_Kcorr;
  //////////////////////////////////////////////////////////
  /**
  *  @brief type of redshift, s or p
  */
  string redshift_type;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  int N_bin_z;
  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */
  int N_bin_z_low_res;

  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec z_min;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec z_max;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec z_min_low_res;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec z_max_low_res;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  int N_bin_mag;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec mK_min;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec mK_max;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  int N_bin_color;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec color_min;

  real_prec color_max;

  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  int N_bin_lmass;
  int N_bin_lmass_lowres;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec lmass_min;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec lmass_max;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  int N_bin_Mag;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  int N_bin_Mag_low_res;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec MK_min;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec MK_max;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  int i_magnitude_limit;

  int i_z_limit;

  unsigned long NGAL;
  unsigned int NCOLS;

  unsigned long NRAN;
  unsigned int NCOLS_R;


  unsigned long NMASK;
  unsigned int NCOLS_MASK;
  
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<real_prec> mK_limits;
  vector<real_prec> z_limits;

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //   Definition of cosmological parameters
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  real_prec Om_matter;
  real_prec Om_cdm;
  real_prec Om_radiation ;
  real_prec Om_baryons;
  real_prec Om_vac;
  real_prec Om_k;
  real_prec Hubble;
  real_prec hubble;
  real_prec n_s;
  real_prec alpha_s;
  real_prec w_eos;
  real_prec N_eff;
  real_prec sigma8;
  real_prec Tcmb;
  real_prec GAL_BIAS;
  real_prec alpha_BIAS;
  real_prec kstar;
  bool use_K_correction;
  bool use_e_correction;
  real_prec RR;
  bool use_wiggles;

  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec total_area;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec z_min_vls;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec deltaz;
  real_prec deltaz_low_res;

  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec area_pixel;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  long n_pixels;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  int nside;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  int N_z_cosmo;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  int nM;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  int imk;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  string hemisphere;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Column added to propwith 1 naad 0 depending whether or not the alaxy is observed
  */
  int i_gal_mask;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Column added to propwith 1 naad 0 depending whether or not the alaxy is observed
  */
  int i_ran_mask;

  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  int iz_vls;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  int i_z;

  //////////////////////////////////////////////////////////

  /**
  *  @brief
  */
  bool use_galaxies_with_zp_AND_zs;

  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec area_in_deg;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<real_prec>  mask;

  vector<vector<real_prec> > mask_old;


  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<real_prec> lf;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void set_catalog(string);
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void set_vec();
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void set_pars(int);
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_random_cat(string, long); // to create random


  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void set_file_names(); // to create random


  string name_cat;


  int observed_pixels_in_mask;
  
  //////////////////////////////////////////////////////////
  /**
  *  @brief Get the index from the mask loaded in the class, Very SLOW
  */
  void get_ipix_from_mask(int, Healpix_Map<real_prec>, long &ipx ); // to create random



  //////////////////////////////////////////////////////////
  /**
  *  @brief  Obtain minimum value of the i-th column of gal properties
  *          The column ic is passed explicitely
  */
  real_prec get_mina(int ic);
  //////////////////////////////////////////////////////////
   /**
   *  @brief Obtain minimum value of the i-th column of random properties
   *          The column ic is passed explicitely
   */
  real_prec get_min_r(int ic);
  //////////////////////////////////////////////////////////
   /**
   *  @brief  Obtain maximum value of the i-th column of gal properties
   *          The column is passed explicitely
   */
  real_prec get_maxa(int);
  //////////////////////////////////////////////////////////
   /**
   *  @brief
   */
  real_prec get_max_r(int);
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_zbins_same_ngal(int, int, real_prec, real_prec,vector< vector<real_prec> >&);
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_mask();

  void get_mask_fits(string);
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_dndz();
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_dndmk();

  //////////////////////////////////////////////////////////
    /**
  *  @brief
  */
  int Ngal_dndz;
  //////////////////////////////////////////////////////////
  /**
  *  @brief get the distribution of the property x
  */
  void get_dndx(string x);
  //////////////////////////////////////////////////////////
  /**
  *  @brief get the joint distribution in the properties x,y
  */
  void get_dndx(string x, string y);

  //////////////////////////////////////////////////////////
  /**
  *  @brief get the Number of objects with property y greater than a value y, in bins of property x
  */
  void get_dndx_cumulative(string x, string y_cum);

  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_dndz_bcg(string);
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_nbar();
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_nbar_bcg();

  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_z_pdf();
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_M_pdf();

  //////////////////////////////////////////////////////////
  /**
  *  @brief  Get cosmological functions
  */
  void get_cosmo(int,bool);

  void read_input_cats(string);

  //////////////////////////////////////////////////////////
  /**
  *  @brief Get Number counts using the Amax estimator
  */
  void get_NS();

  //////////////////////////////////////////////////////////
  /**
  *  @brief Get luminosity function
  */
  void get_Vmax();
  //////////////////////////////////////////////////////////
  /**
  *  @brief Get luminosity function
  */
  void get_LF_Vmax(ULONG);
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_LF_Cole();

  void get_LF_Cole2();

  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_dndz_LF();
  void get_dndz_LF_bins();


  //////////////////////////////////////////////////////////
  /**
  *  @brief Generate a Healpix Mask masking north and south according
  * to the coordinate system passed as argument.
  */
  void get_new_mask(string coord);
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_new_cat();

  void read_get_Kcorr();

  void get_alpha_rsd();

  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_new_cat_bcg(string,string);
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void bspline(int, int);
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void set_nbar(vector<real_prec>,vector<real_prec>,vector<real_prec>);
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_jacknife_mask(); // specfic for 2mpz
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void read_pars(string &file);

  bool Measure_M_pdf;

  int index_mask;

  ULONG Ngal_expected;
  
  string output_dir;

  real_prec dec_max;
  real_prec dec_min;
  real_prec phi_max;
  real_prec phi_min;


  gal_parameters gp;


};


#endif
