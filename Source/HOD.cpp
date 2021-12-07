# include "../Headers/HOD.h"
# include "../Headers/Constants.h"
# include "../Headers/Type_structures_def.h"

real_prec HOD::CENTRAL(real_prec M, s_CosmologicalParameters *scp){
  real_prec f;
  real_prec mmin=scp->mmin_hod;
  real_prec sm=scp->scatter_hod;
  switch(scp->hod_model){
    case(1): f=(M<mmin? 0.0 : 1.0);break;
    case(2): f=exp(-mmin/M);      break;
    case(3): f=0.5*(1.+gsl_sf_erf(log10(M/mmin)/sm));   break;
  }
  return f;
}

real_prec HOD::SATELLITE(real_prec M, s_CosmologicalParameters *scp){
  real_prec f;
  real_prec kappa=0.137;
  real_prec mmin=scp->mmin_hod;
  real_prec muno=scp->muno_hod;
  real_prec al=scp->alpha_hod;
  switch(scp->hod_model){
    case(1): f=pow(M/muno,al)*CENTRAL(M,scp); break;
    case(2): f=pow(M/muno,al)*CENTRAL(M,scp); break;
    case(3): f=(M<mmin? 0.0 : pow((M-kappa*mmin)/muno,al))   ;break;
  }
  return f;
}
