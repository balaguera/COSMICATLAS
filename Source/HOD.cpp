# include "../Headers/NumericalMethods.h"
# include "../Headers/Constants.h"
# include "../Headers/Type_structures_def.h"
# include "../Headers/HOD.h"

double HOD::CENTRAL(double M, s_CosmologicalParameters *scp){
  double f;

  double mmin=scp->mmin_hod;
  double sm=scp->scatter_hod;
  switch(scp->hod_model){ 
    case(1): f=(M<mmin? 0.0 : 1.0);break;
    case(2): f=exp(-mmin/M);      break;
    case(3): f=0.5*(1.+gsl_sf_erf(log10(M/mmin)/sm));   break;
  }
  return f;
}  

double HOD::SATELLITE(double M, s_CosmologicalParameters *scp){
  double f;
  double mmin=scp->mmin_hod;
  double muno=scp->muno_hod;
  double al=scp->alpha_hod;
  switch(scp->hod_model){  
    case(1): f=pow(M/muno,al)*CENTRAL(M,scp); break;
    case(2): f=pow(M/muno,al)*CENTRAL(M,scp); break;
    case(3): f=(M<mmin? 0.0 : pow((M-mmin)/muno,al))   ;break;
  }
  return f;
}

