/**
* @class
* @file DensityProfiles.h
* @brief Header file for the class DensityProfiles::
* @title Functions related to the generation of density profiles for dark matter haloes
* @details Bias Assignment method for mock catalogs
* @author Andres Balaguera-Antolínez
* @version 1.0
* @date    2020
* @details: This is an example of a main function to call bam. A file called cosmicatlas.cpp
*/


#ifndef _DensityProfiles_
#define _DensityProfiles_

# include "NumericalMethods.h"

using namespace std;


class DensityProfiles{

  static gsl_real idensity_rc(gsl_real, void *);
  static gsl_real idensity_fourier_c(gsl_real, void *);
  static gsl_real idensity_fourier(gsl_real, void *);
  void einasto_parameters(real_prec, real_prec, real_prec *, real_prec *, real_prec *,real_prec *,real_prec *, void *);
  void nfw_parameters(real_prec , real_prec, real_prec *, real_prec *,real_prec *,real_prec *,real_prec *, void*);
 public:
  DensityProfiles(){}
  ~DensityProfiles(){}
  real_prec density_r(real_prec,real_prec, real_prec, void *);
  real_prec density_fourier(real_prec,real_prec, real_prec, void *);

  real_prec mass_concentration_dis(real_prec,real_prec,real_prec, void*);
  real_prec density_rc(real_prec,real_prec, real_prec,  void *);
  real_prec density_k(real_prec,real_prec, real_prec,  s_CosmologicalParameters *);
  real_prec density_kc(real_prec,real_prec, real_prec, void *);
};





#endif
