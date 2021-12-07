// **************************************************************************************************
// **************************************************************************************************
/**
 * @class<Cosmology>
 * @brief Header file for the class Cosmology::
 * @file CosmologicalFunctions.h
 * @title Functions to compute cosmological dependent quantities
 * @author Andres Balaguera-Antolínez
 * @version 1.0
 * @date  2012-2020
 * @details The conversio of void * structure to the type s_CosmologicalParametes is
 * @code
     struct s_CosmologicalParameters * s_cp = (struct s_CosmologicalParameters *)p;
 * @endcode


 */
// **************************************************************************************************
// **************************************************************************************************



#ifndef __COSMOLOGY__
#define __COSMOLOGY__

# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include "Params.h"
# include "Type_structures_def.h"
# include "Constants.h"
using namespace std;
using namespace Constants;





#define ERROR_LIMIT_NR 0.01
// This demands 1% convergence precision for the Newton-Rhapson algorithm


class Cosmology{
  
 private:
  //////////////////////////////////////////////////////////
  /** 
   * @private
   * @brief Object of type CosmologicalFunctions 
   */
  static gsl_real gint(gsl_real, void *);
  //////////////////////////////////////////////////////////
  /** 
   * @private
   * @brief Object of type CosmologicalFunctions
   */
  static gsl_real i_rs(gsl_real, void *);
  //////////////////////////////////////////////////////////
  /** 
   * @private
   * @brief Inverse of the Hubble function H(z) used for integration
   */
  static gsl_real Hinv(gsl_real, void *);
  //////////////////////////////////////////////////////////
  /** 
   * @private
   * @brief Inverse of the Hubble function H(a) used for integration
   */
  static gsl_real Hinva(gsl_real, void *);
  //////////////////////////////////////////////////////////
  /** 
   * @private
   * @brief GSL function implemented in root finder (Newton-Rhapson)
   */
  static gsl_real Froot(gsl_real , void*);
  //////////////////////////////////////////////////////////
  /** 
   * @private
   * @brief GSL function implemented in root finder (Newton-Rhapson)
   */
  static gsl_real dFroot(gsl_real , void *);
  //////////////////////////////////////////////////////////
  /** 
   * @private
   * @brief GSL function implemented in root finder (Newton-Rhapson)
   */
  static void F_dF(gsl_real , void *, gsl_real *, gsl_real *);
    //////////////////////////////////////////////////////////
  /** 
   * @private
   * @brief Number of points to gsl integration
   */
  int NP;
  //////////////////////////////////////////////////////////
  /** 
   * @private
   * @brief Object of type CosmologicalFunctions
   */
  gsl_integration_glfixed_table *wf;


  
 public:
  
  //////////////////////////////////////////////////////////
  /**
   * @brief Default constructor
   */
  Cosmology(){}  // Default constructor
  //////////////////////////////////////////////////////////
  /**
   * @brief Default destructor
   */
  ~Cosmology(){}
  //////////////////////////////////////////////////////////

  /**
   * @brief Constructor
   * @param Np Number of points for the Gauss-Legendre integration of cosmological functions.
   * @return wf gsl_table for integration
   */
  Cosmology(int np){
    NP=np;
    wf=gsl_integration_glfixed_table_alloc(NP);
  } 

  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param p_structure void pointer to s_CosmologicalParameters structure
   * @return Hubble function at the input redshift

   */
  real_prec Hubble_function(real_prec z, void * p_structure);
    //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param p_structure void pointer to s_CosmologicalParameters structure
   * @return Comoving distance at the input redshift
   */
  real_prec comoving_distance(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   * @return Proper_angular_diameter_distance at the input redshift
   */
  real_prec proper_angular_diameter_distance(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   * @return Comoving_angular_diameter_distance at the input redshift
   */
  real_prec comoving_angular_diameter_distance(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   * @return Proper_angular_diameter_distance to be interpolated
   */
  real_prec inter_proper_angular_diameter_distance(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   * @return Luminosity Distance
   */
  real_prec luminosity_distance(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   * @return Luminosity Distance  to be interpolated
   */
  real_prec inter_luminosity_distance(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec transverse_comoving_distance(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec inter_transverse_comoving_distance(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec derivative_transverse_comoving_distance(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec mean_matter_density(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec age_universe(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec comoving_sound_horizon(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec growth_factor(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec growth_index(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec growth_index2(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec halo_dynamical_time(real_prec, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec omega_matter(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec omega_radiation(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec omega_curvature(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec omega_dark_energy(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec Distance_Modulus(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec inter_Distance_Modulus(real_prec z, void *p_structure);
    //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec zmax(void *);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec zmax_old(void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec zmax_old(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   * @return Comoving distance at input redshift
   */
  real_prec rr(real_prec , real_prec , void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   * @return Top-hot Density contrast
   */
  real_prec density_contrast_top_hat(real_prec z, void* p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   * @return Critical density at input redshift
   */
  real_prec critical_density(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   * @return K-correction at input redshift
   */
  real_prec K_correction(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /**
   * @public
   */
  real_prec K_correction(real_prec,real_prec,real_prec, void *p_structure);

  //////////////////////////////////////////////////////////
  /** 
   * @brief Object of type CosmologicalFunctions 
   */
  real_prec dK_correction_dz(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   * @return e-correction at input redshift
   */
  real_prec e_correction(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @param *p_structure void pointer to s_CosmologicalParameters structure
   */
  real_prec de_correction_dz(real_prec z, void *p_structure);
  //////////////////////////////////////////////////////////
  /** 
   * @brief Object of type CosmologicalFunctions 
   */
  void free_gsl_table();
  //////////////////////////////////////////////////////////
  /** 
   * @brief Object of type CosmologicalFunctions 
   */
  void check_cosmo_pars(s_CosmologicalParameters *sp);
  //////////////////////////////////////////////////////////
  /** 
   * @brief Object of type CosmologicalFunctions 
   */
  real_prec gsl_cosmo_integration(gsl_real (*function)(gsl_real, void *) ,void *p, gsl_real LowLimit,gsl_real UpLimit);
  //////////////////////////////////////////////////////////
  /** 
   * @brief Object of type CosmologicalFunctions 
   */
  void Comoving_distance_tabulated (real_prec, real_prec, void *, vector<gsl_real>&, vector<gsl_real>&);


};

#endif 
