//#define TIME
#undef TIME
#ifndef _COSMOLIB_
#define _COSMOLIB_

# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>

# include "NumericalMethods.h"
# include "Parameters_CosmoLib.h"
# include "CosmologicalFunctions.h"
# include "FileOutput.h"
# include "ScreenOutput.h"
# include "PowerSpectrumTH.h"
# include "CorrelationFunction.h"
# include "Astrophysics.h"
# include "ScalingRelations.h"
# include "Statistics.h"
# include "BiasFunctions.h"
# include "DensityProfiles.h"
# include "HOD.h"
# include "Marks.h"
# include "GalaxyOperations.h"
# include "AngularPowerSpectrum.h"
# include "McmcFunctions.h"
# include "Galaxy.h"

using namespace std;
using namespace Constants;


class CosmoLib{
 private:
    //////////////////////////////////////////////////////////
    /**
     *  @brief get the value of the the private member type_of_lbinning
     *  @return type_of_lbinning
     */

    Cosmology Cf;
    //////////////////////////////////////////////////////////
    /**
     *  @brief get the value of the the private member type_of_lbinning
     *  @return type_of_lbinning
     */
    FileOutput Fm;
    //////////////////////////////////////////////////////////
    /**
     *  @brief get the value of the the private member type_of_lbinning
     *  @return type_of_lbinning
     */
    Statistics Cs;
    //////////////////////////////////////////////////////////
    /**
     *  @brief get the value of the the private member type_of_lbinning
     *  @return type_of_lbinning
     */
    PowerSpectrum Ps;

    ScreenOutput So;

public:
    //////////////////////////////////////////////////////////
    /**
     *  @brief get the value of the the private member type_of_lbinning
     *  @return type_of_lbinning
     */

  CosmoLib(){}  // Default constructor
  //////////////////////////////////////////////////////////
  /**
   *  @brief Constructor passed the impot parameter file to
   *  @return 
   */
  CosmoLib(string par_file){
    ParametersCosmolib param (par_file);
    this->params=param;
    Statistics Csa(this->params._M_min_mf(), params._M_max_mf(), params._n_points_mf());
    this->Cs=Csa;
    PowerSpectrum Psa(params._k_min_integration() ,params._k_max_integration() ,params._n_points_dp_k(),10, params._M_min_mf(),params._M_max_mf(),params._n_points_mf());
    this->Ps=Psa;
  }


  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */

   ParametersCosmolib params;
   //////////////////////////////////////////////////////////
   /**
    *  @brief get the value of the the private member type_of_lbinning
    *  @return type_of_lbinning
    */
   s_CosmologicalParameters s_cosmo_par;
   //////////////////////////////////////////////////////////
   /**
    *  @brief get the value of the the private member type_of_lbinning
    *  @return type_of_lbinning
    */

  ~CosmoLib(){}


   //////////////////////////////////////////////////////////
   /**
    *  @brief get the value of the the private member type_of_lbinning
    *  @return type_of_lbinning
    */
  vector<gsl_real> v_mass;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_sigma_mass;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_mass_function;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_halo_mass_bias;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_effective_halo_mass_bias;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_effective_halo_mean_number_density;

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */

  vector<gsl_real> v_k_ps;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_nl_power_spectrum;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_nl_power_spectrum_pt;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */

  vector<gsl_real> v_l_power_spectrum;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */

  vector<gsl_real>pk_aux;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real>kk_aux;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */

  vector<real_prec> v_r_cf;

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<real_prec> v_nl_correlation_function;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<real_prec> v_l_correlation_function;

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */

  vector<real_prec> v_r_dp;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */

  vector<real_prec> v_density_profile_r;

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */

  vector<gsl_real> v_k_dp;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */

  vector<gsl_real> v_density_profile_k;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */


  vector<gsl_real> v_galaxy_power_1h_ss;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_galaxy_power_1h_sc;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_galaxy_power_2h;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_galaxy_matter_bias;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_galaxy_power_spectrum;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */


  vector<real_prec> v_galaxy_correlation_1h_ss;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<real_prec> v_galaxy_correlation_1h_sc;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<real_prec> v_galaxy_correlation_2h;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<real_prec> v_galaxy_correlation;


  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */

  void feed_cosmo_struct();

  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */

  void get_cosmolib();
  //////////////////////////////////////////////////////////
  /**
   *  @brief This does the same as the constructor above, but as an object of the type CosmoLib
   *  @return
   */
  void set_par_file(string);


};

#endif 





