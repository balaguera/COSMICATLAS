
/**
* @namespace<Cosmo_parameters_test,Cosmo_parameters_test>
* @file cosmo_parametrs.h
* @title Cosmological parameters from different LSS analysis
* @author Andres Balaguera-Antolínez
* @version 1.0
* @date    2020
* @details: This is an example of a main function to call bam. A file called cosmicatlas.cpp
*/

#ifndef _cosmo_pars_
#define _cosmo_pars_

namespace Cosmo_parameters_test
{
  double Omega0=0.23;
  double Neff=3.046;
  double Omegarad=(1.3157e-5)*(1.+0.2271*Neff);
  double Omegavac=0.71;
  double Omegabaryon=0.0441;
  double Omegak=0.001;//1.-Omega0-Omegavac-Omegarad;
  double w=-1;
  // hubble constant H/100km/s/Mpc
  double hubble = 0.73;  
  double Hubble=100.0;
  // fraction of baryonic energy density to cdm?
  double f_baryon = Omegabaryon/Omega0; 
  double Tcmb = 2.725;              
  double spectral_index = 0.96;   
  double sigma8= 0.83;            
  double RR=8.0;

}


namespace Cosmo_parameters_PLANCK
{
  double Om_matter=0.3089;
  double N_eff=3.046;
  double Om_radiation=0;//(1.3157e-5)*(1.+0.2271*N_eff);
  double Om_baryons=0.044;
  double Om_k=0.00;//1.-Omega0-Omegavac-Omegarad;
  double w_eos=-1;
  double Om_vac=1.-Om_matter-Om_k-Om_radiation;

  // hubble constant H/100km/s/Mpc
  double hubble = 0.6774;
  double Hubble=100.0;

  // fraction of baryonic energy density to matter
  double f_baryon = Om_baryons/Om_matter;
  double Tcmb = 2.725;
  double n_s = 0.9667;
  double alpha_s = 0.;
  double sigma8= 0.8147;
  double RR=8.0;
}



namespace Cosmo_parameters_Minerva
{
  double Om_matter=0.285;
  double N_eff=3.046;
  double Om_radiation=0;//(1.3157e-5)*(1.+0.2271*N_eff);
  double Om_vac=0.715;
  double Om_baryons=0.044;
  double Om_k=0.00;//1.-Omega0-Omegavac-Omegarad;
  double w_eos=-1;
  // hubble constant H/100km/s/Mpc
  double hubble = 0.695;
  double Hubble=100.0;

  // fraction of baryonic energy density to matter
  double f_baryon = Om_baryons/Om_matter;
  double Tcmb = 2.725;
  double n_s = 0.9632;
  double alpha_s = 0.;
  double sigma8= 0.828;
  double RR=8.0;

}




#endif

