#ifndef __POWER_SPECTRUMTH__
#define __POWER_SPECTRUMTH__

// CLASS TO COMPUTE LINEAR AND NON LINEAR MATTER POWER SPECTRUM
// BASED ON THE EISENSTEIN AND HU FITTING FORMULAE
// AND THE HALO-FIT FROM SMITH ET AL.


# include <iostream>
# include <math.h>
# include <cmath>
# include <cctype>
# include <stdio.h>
# include <fstream>
# include <omp.h>
# include "NumericalMethods.h"
# include "CosmologicalFunctions.h"
# include "ScreenOutput.h"


using namespace std;

class PowerSpectrum{
 private:
 //////////////////////////////////////////////////////////

  ScreenOutput So;

  /**
  * @brief
 */
  static gsl_real isigma(gsl_real ,void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real dsigma_dR(gsl_real, void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real fun(gsl_real ,void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real fun_nw(gsl_real ,void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real fun_aux_halo_fit(gsl_real , void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real fun_aux_halo_fit2(gsl_real , void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real fun_aux_halo_fit4(gsl_real , void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */

  static gsl_real fun_aux_halo_fit_dw(gsl_real , void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real fun_aux_halo_fit2_dw(gsl_real , void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real fun_aux_halo_fit4_dw(gsl_real , void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real Power_Spectrum_i(gsl_real, void*);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real iP1loop(gsl_real, void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real iFkernel(gsl_real, void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */

  static gsl_real i_Galaxy_power_spectrum_h1_ss(gsl_real, void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real i_Galaxy_power_spectrum_h1_sc(gsl_real, void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real i_Galaxy_matter_bias(gsl_real, void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real i_mean_galaxy_number_density(gsl_real, void*);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */


 public:

  //////////////////////////////////////////////////////////
  /**
   * @brief
  */

  PowerSpectrum(){}
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */

  PowerSpectrum(real_prec k1, real_prec k2, int nk, int nm){
    compute_int_table_k_mu(k1,k2, nk, nm);
  }
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */

  PowerSpectrum(real_prec M1, real_prec M2, int Np){
     compute_int_table_mass(M1,M2, Np);
  }
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec omhh; /* Omega_matter*h^2 */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec obhh;   /* Omega_baryon*h^2 */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec theta_cmb;  /* Tcmb in units of 2.7 K */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec z_equality; /* Redshift of matter-radiation equality, really 1+z */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec k_equality; /* Scale of equality, in Mpc^-1 */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec z_drag;   /* Redshift of drag epoch */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec R_drag;   /* Photon-baryon ratio at drag epoch */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec R_equality; /* Photon-baryon ratio at equality epoch */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec sound_horizon;  /* Sound horizon at drag epoch, in Mpc */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec k_silk;   /* Silk damping scale, in Mpc^-1 */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec alpha_c;  /* CDM suppression */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec beta_c;   /* CDM log shift */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec alpha_b;  /* Baryon suppression */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec beta_b;   /* Baryon envelope shift */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec beta_node;  /* Sound horizon shift */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec k_peak;   /* Fit to wavenumber of first peak, in Mpc^-1 */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec sound_horizon_fit;  /* Fit to sound horizon, in Mpc */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec alpha_gamma;  /* Gamma suppression in approximate TF */

  //////////////////////////////////////////////////////////
  /**
   * @brief
  */

  PowerSpectrum(real_prec k1, real_prec k2, int nk, int nmu,real_prec M1, real_prec M2,real_prec Nm){
    compute_int_table_k_mu(k1,k2, nk, nmu);
    compute_int_table_mass(M1,M2, Nm);
  }

  //////////////////////////////////////////////////////////
  /**
   * @brief
  */



  ~PowerSpectrum(){}
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */

  real_prec Linear_Matter_Power_Spectrum(s_CosmologicalParameters *, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Linear_Matter_Power_Spectrum_z(s_CosmologicalParameters *, real_prec, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Linear_Matter_Power_Spectrum_interpolated(s_CosmologicalParameters *, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Q_Model_Matter_Power_Spectrum(s_CosmologicalParameters *scp,real_prec k);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Primordial_Matter_Power_Spectrum(s_CosmologicalParameters *, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Linear_Matter_Power_Spectrum_NW(s_CosmologicalParameters *, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Linear_Matter_Power_Spectrum_DW(s_CosmologicalParameters *, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Linear_Matter_Power_Spectrum_G_DW(s_CosmologicalParameters *, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Linear_Matter_Power_Spectrum_G_NW(s_CosmologicalParameters *, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Non_Linear_Matter_Power_Spectrum_Halo_Fit(s_CosmologicalParameters *, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Non_Linear_Matter_Power_Spectrum_Halo_Fit_z(s_CosmologicalParameters *, real_prec, real_prec, real_prec, real_prec, real_prec, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
   real_prec Non_Linear_Matter_Power_Spectrum_PT(s_CosmologicalParameters *scp, real_prec k);
   //////////////////////////////////////////////////////////
   /**
    * @brief
   */
  real_prec P1loop(s_CosmologicalParameters *scp, real_prec k);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec sigma_masa(real_prec, real_prec, s_CosmologicalParameters *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec bias(real_prec, real_prec, s_CosmologicalParameters *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec mass_function_D(real_prec, real_prec, s_CosmologicalParameters *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void mass_function_M_Z(vector<real_prec>, s_CosmologicalParameters *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<vector<real_prec> > MASS_FUNCTION_M_Z;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<vector<real_prec> > MASS_BIAS_M_Z;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec mean_galaxy_number_density(real_prec redshift, s_CosmologicalParameters *scp);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec mean_galaxy_number_density(s_CosmologicalParameters *scp);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Non_Linear_Matter_Power_Spectrum_Halo_Fit_DW(s_CosmologicalParameters *, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Galaxy_power_spectrum_h1_ss(s_CosmologicalParameters *, real_prec, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Galaxy_power_spectrum_h1_sc(s_CosmologicalParameters *, real_prec, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Galaxy_matter_bias(s_CosmologicalParameters *, real_prec, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void TFset_parameters(real_prec omega0hh , real_prec f_baryon, real_prec Tcmb);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec TFfit_onek(real_prec , real_prec &, real_prec &);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void TFfit_hmpc(real_prec, real_prec, real_prec , real_prec,
                  int, real_prec *, real_prec *, real_prec *, real_prec *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec TFsound_horizon_fit(real_prec, real_prec , real_prec );
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec TFk_peak(real_prec , real_prec , real_prec );
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec TFnowiggles(real_prec, real_prec, real_prec,
                     real_prec , real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec TFzerobaryon(real_prec, real_prec , real_prec, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec window(real_prec,real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec windowg(real_prec,real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void halo_fit(real_prec, void*, real_prec*, real_prec*, real_prec*);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void halo_fit_integrals(void*, real_prec*, real_prec *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void halo_fit_z(real_prec , void *,  real_prec, real_prec,real_prec, real_prec, real_prec, real_prec *,real_prec *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void normalization(void *, real_prec &);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec normalization(void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void nl_scales_halo_fit(void *, real_prec *, real_prec*, vector<real_prec>&, vector<real_prec>&, bool);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void kstar_integral(void *, real_prec*);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void hf_aux(real_prec, real_prec,real_prec, real_prec, real_prec *, real_prec *, real_prec *, real_prec *,real_prec *,real_prec *,real_prec *, real_prec *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void compute_int_table_k_mu(real_prec, real_prec, int,  int);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void compute_int_table_mass(real_prec, real_prec, int);
  // These vectors are allocated after the computation of the linear power spectrum,
  // and are meant to speed up the calculation of the non-linear power spectrum using halo fit.
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> v_Pk_linear;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> v_kk;

  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  gsl_integration_glfixed_table *wfd;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  gsl_integration_glfixed_table *wf;

  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> v_mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> v_sigma_mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> v_mass_function;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> v_halo_mass_bias;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> XX_Mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> WW_Mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> XX;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> WW;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> XX_mu;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> WW_mu;


  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> kvector_external;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> power_external;

  bool use_external_power;

};

#endif  
  



