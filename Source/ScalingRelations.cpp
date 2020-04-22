# include "../Headers/ScalingRelations.h"
# include "../Headers/CosmologicalFunctions.h"



/*MASS-TEMPERATURE RELATION CALIBRATED(?), MASSES IN UNITS OF 10^14 Ms/h: RETURNS TEMPERATURE IN keV*/
double SCALING_RELATIONS::M2T(double m, void *p){
  return 0.595*pow(m,(double)0.34);
}

/*Converting mass in 10^{14}/h to the NATURAL LOGARITHM of 
  luminosities (bolometric or band) in units of 10^{44}erg/s/h^2  */
double SCALING_RELATIONS::RB_M2L(double m, void *p){/*Reipricht & Bohringer 2002 */
  return log(0.1169*pow(m,(double)1.462));
}

double SCALING_RELATIONS::FED_M2L(double m, void *p){/*Fedeli et al*/
  return log(0.0862*pow(m,(double)1.554));
}
double SCALING_RELATIONS::STA_M2L(double m, void *p){/*Measured, for Omega_matter=0.24 and sigma_ln M=0.39 leading to the */
  return log(0.1139*pow(m,(double)1.46));
}

double SCALING_RELATIONS::MANTZ_BOL_M2L(double m, void *p){
  /*Measured (bolometric) by Mantz et al 2009, also based on XLF of REFLEX data. The factors hubble^2 and 0.1
    come to transform from the units of Mantz ([M]=[10^15 Ms], [L]=10^44 erg/s) to our units ([M]=[10^14 Ms/h], [L]=10^44 erg/s/h^2) */
  
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p; 
  
  return log(pow(s_cp->hubble,2)*1.26*pow(0.1*m/s_cp->hubble,(double)1.59));
}

double SCALING_RELATIONS::MANTZ_BAND_M2L(double m, void *p){/*Measured in the ROSAT energy band by Mantz et al 2009, also based on XLF REFLEX data */
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p; 
  return log(pow(s_cp->hubble,2)*0.82*pow(0.1*m/s_cp->hubble,(double)1.29));   
}

double SCALING_RELATIONS::MOCKS_M2L(double m, void *p){ /*Calibrated for the reflex2 mocks to follow the reflex2 luminosity function. */
  // double a_ml       = -1.10;
  // double b_ml       =  1.70;
  // double c_ml       = -0.23;
  double a_ml       = -0.64;
  double b_ml       =  1.27;
  double c_ml       =  0.0;
  return  log(pow(10,a_ml+b_ml*log10(m)+c_ml*pow(log10(m),2)));
}
