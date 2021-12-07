# include "../Headers/CosmologicalFunctions.h"
# include "../Headers/PowerSpectrumTH.h"
# include "../Headers/DensityProfiles.h"
# include "../Headers/HOD.h"
# include "../Headers/Marks.h"
# include "../Headers/BiasFunctions.h"

// *******************************************************************
// We will use this function to compute the HaloFit power spectrum or any other function
// that makes use of the linear matter power spectrum. Therefore, we compute Plin first
// and then we just interpolated


void PowerSpectrum::compute_int_table_k_mu(real_prec k_min_integration, real_prec k_max_integration, int nss_k, int nss_mu){
  wfd =  gsl_integration_glfixed_table_alloc (nss_k);
  this->WW.resize(nss_k);
  this->XX.resize(nss_k);
  gsl_get_GL_weights(static_cast<gsl_real>(log10(k_min_integration)),static_cast<gsl_real>(log10(k_max_integration)),this->wfd,this->XX,this->WW);

  // INtegrate with respect to mu
  wf =  gsl_integration_glfixed_table_alloc (nss_mu);
  this->WW_mu.resize(nss_mu);
  this->XX_mu.resize(nss_mu);
  gsl_get_GL_weights(-1.0, 1.0,this->wf,this->XX_mu,this->WW_mu);
}


// *******************************************************************

void PowerSpectrum::compute_int_table_mass(real_prec M_min_integration, real_prec M_max_integration, int nss_k){
  wfd =  gsl_integration_glfixed_table_alloc (nss_k);
  this->WW_Mass.resize(nss_k);
  this->XX_Mass.resize(nss_k);
  gsl_get_GL_weights(static_cast<gsl_real>(1.1*log10(M_min_integration)),static_cast<gsl_real>(0.9*log10(M_max_integration)),this->wfd,this->XX_Mass,this->WW_Mass);

}

// *******************************************************************



real_prec PowerSpectrum::Linear_Matter_Power_Spectrum_interpolated(s_CosmologicalParameters *scp, real_prec k){        /*this k comes in h/Mpc */
   return gsl_inter_new(scp->v_k_ps, scp->v_lin_power_spectrum, k);
}

// *******************************************************************


real_prec PowerSpectrum::Non_Linear_Matter_Power_Spectrum_PT(s_CosmologicalParameters *scp, real_prec k){
    real_prec g   = scp->growth_factor;

    // The growth factor is factorized out and multiplied at the end. Hence we pass here 1.0 in the third argument
    real_prec P1  = exp(-0.5*pow(k/scp->kstar,2))*this->Linear_Matter_Power_Spectrum_z(scp,k,1.0);
    real_prec P2 = scp->Amc*this->P1loop(scp, k);
     return pow(g,2)*P1+ pow(g,4)*P2;
}




// *******************************************************************
// P1looop evaluated at z=0
real_prec PowerSpectrum::P1loop(s_CosmologicalParameters *scp, real_prec k){
    s_aux<PowerSpectrum> ssa;
    ssa.kaux=k;
    ssa.scp_a=scp;
    ssa.WW=this->WW_mu;
    ssa.XX=this->XX_mu;
    return gsl_integration2(iP1loop, (void *)&ssa,this->XX, this->WW)/(2.*M_PI*M_PI);
}

// *******************************************************************
gsl_real PowerSpectrum::iP1loop(gsl_real lq, void *p){
   struct s_aux<PowerSpectrum> * saux = (struct s_aux<PowerSpectrum> *)p;
   struct s_CosmologicalParameters * scpa= saux->scp_a;
   PowerSpectrum Ps;
   real_prec q=pow(10,lq);

   // Need to redefine here the aux structure to pass it at iFkernel
   s_aux<PowerSpectrum> ssaa;
   ssaa.kaux=saux->kaux; //k
   ssaa.raux=q; //q
   ssaa.scp_a=scpa;  //sca

   // INtegral with respect to mu = k*q
   real_prec IF=gsl_integration2(iFkernel, (void *)&ssaa, saux->XX, saux->WW);

   // The Plin is evaluated at z = 0  and extrapolated with the growth factor
   return static_cast<gsl_real>((log(10.0)*q)*q*q*Ps.Linear_Matter_Power_Spectrum_z(scpa, q, 1.0)*IF);
}

// *******************************************************************
// *******************************************************************
gsl_real PowerSpectrum::iFkernel(gsl_real mu, void *p){
    struct s_aux<PowerSpectrum> * saux = (struct s_aux<PowerSpectrum> *)p;
    struct s_CosmologicalParameters * scpa= saux->scp_a;
    PowerSpectrum Ps;
    real_prec k = saux->kaux;
    real_prec q = saux->raux;
    real_prec kmq=sqrt(q*q+k*k-2.0*k*q*mu);
/*
    // Terms com,ing from int !F(k-q,q)!^2 d phi
    real_prec A=17./21.0+((mu*k-q)/(2.*q))*(1+pow(q/kmq,2))+(4.*mu*mu*pow(k-mu*q,2))/(21.0*kmq*kmq) - (4*q*mu*(k-mu*q)*(1-mu*mu))/(21.0*kmq*kmq);
    real_prec B= (4.0*q*q*(1.0-mu*mu)/(21*kmq*kmq));
    real_prec C=-(4.0*q*mu*(1-mu*mu)*(k-mu*q))/(21.0*kmq*kmq);
    real_prec D=-(4.0*q*(k-q*mu)*pow(fabs(1-mu*mu), 1.5))/(21.0*kmq*kmq);
    //real_prec F2 = (M_PI/8.)*(16*A*A+28.*A*B+13.*B*B+5.0*C*C+2.*C*D+D*D);
*/
    real_prec F2 = (5./7.)+0.5*((mu*k-q)/q)*(1.+pow(q/kmq,2))+(2./7.)*pow(mu*k-q,2)/(kmq*kmq);
    F2 = pow(F2,2);


//    real_prec F2=pow( (17./21.)+((mu*k-q)/(2.*q))*(1+pow(q/kmq,2))+(2./7.)*(pow(mu*k-q,2)/(kmq*kmq)-1./3.) , 2 );

    // The Plin is evaluated at z = 0  and extrapolated with the growth factor
    return static_cast<gsl_real>(F2*Ps.Linear_Matter_Power_Spectrum_z(scpa,kmq, 1.0));
  }




// *******************************************************************
// *******************************************************************

real_prec PowerSpectrum::Linear_Matter_Power_Spectrum(s_CosmologicalParameters *scp,real_prec k){        /*this k comes in h/Mpc */
  // *********************************************
  // Linear matter power spectrum 
  // *********************************************
     real_prec tf_thisk, baryon_piece, cdm_piece;
     TFset_parameters((scp->Om_matter)*scp->hubble*scp->hubble, scp->f_baryon, scp->Tcmb);
     tf_thisk = TFfit_onek(k*scp->hubble, baryon_piece, cdm_piece);
     real_prec Tfunc = true==scp->use_wiggles ? pow(fabs(tf_thisk),2.): pow(fabs(cdm_piece),2.) ;
     real_prec power=scp->pk_normalization*Tfunc*Primordial_Matter_Power_Spectrum(scp,k)*pow(scp->growth_factor,2);
     return power;
  }

// *******************************************************************
// *******************************************************************
real_prec PowerSpectrum::Linear_Matter_Power_Spectrum_z(s_CosmologicalParameters *scp,real_prec k, real_prec g){        /*this k comes in h/Mpc */
  // *********************************************
  // Linear matter power spectrum
  // *********************************************
  real_prec tf_thisk, baryon_piece, cdm_piece;
  TFset_parameters((scp->Om_matter)*scp->hubble*scp->hubble, scp->f_baryon, scp->Tcmb);
  tf_thisk = TFfit_onek(k*scp->hubble, baryon_piece, cdm_piece);
  real_prec Tfunc = true==scp->use_wiggles ? pow(fabs(tf_thisk),2.): pow(fabs(cdm_piece),2.) ;
  return scp->pk_normalization*Tfunc*this->Primordial_Matter_Power_Spectrum(scp,k)*pow(g,2);
}


// *******************************************************************
// *******************************************************************

real_prec PowerSpectrum::Q_Model_Matter_Power_Spectrum(s_CosmologicalParameters *scp,real_prec k){        /*this k comes in h/Mpc */

  // *********************************************
  // Linear matter power spectrum 
  // *********************************************

  real_prec ans;
  real_prec tf_thisk, baryon_piece, cdm_piece;
  real_prec growth_factor=scp->growth_factor;
  real_prec normalization=scp->pk_normalization;
  real_prec alpha_s=scp->alpha_s;
  real_prec A_ps=scp->A_PS;
  real_prec Q_ps=scp->Q_PS;
  TFset_parameters((scp->Om_matter)*scp->hubble*scp->hubble, scp->f_baryon, scp->Tcmb);
  tf_thisk = TFfit_onek(k*scp->hubble, baryon_piece, cdm_piece);
  if(true==scp->use_wiggles)ans=Primordial_Matter_Power_Spectrum(scp,k)*pow(fabs(tf_thisk),2.); 
  else if(!scp->use_wiggles)ans=Primordial_Matter_Power_Spectrum(scp,k)*pow(fabs(cdm_piece),2.); 
  return normalization*ans*(1+k*Q_ps)/(1+A_ps*k)*pow(growth_factor,2);
}




// *******************************************************************
// *******************************************************************
real_prec PowerSpectrum::Primordial_Matter_Power_Spectrum(s_CosmologicalParameters *scp,real_prec k){        /*this k comes in h/Mpc */

  // *********************************************
  // Primordial matter power spectrum 
  // *********************************************
  real_prec ans;
  real_prec tf_thisk, baryon_piece, cdm_piece;
  real_prec A_s=scp->A_s;
  real_prec alpha_s=scp->alpha_s;
  //return A_s*exp( (scp->n_s-1.)*log(k/0.05)+alpha_s*pow(log(k/0.05),2));
  //  return A_s*pow(k/0.05,scp->n_s);
  return pow(k,scp->n_s);
}

// *******************************************************************
// *******************************************************************

real_prec PowerSpectrum::Linear_Matter_Power_Spectrum_NW(s_CosmologicalParameters *scp, real_prec k){        /*this k comes in h/Mpc */

  // *********************************************
  // Linear matter power spectrum 
  // *********************************************

  real_prec nm,nm_w, ans;
  real_prec tf_thisk, baryon_piece, cdm_piece;
  real_prec growth_factor=scp->growth_factor;
  real_prec alpha_s=scp->alpha_s;
  real_prec normalization=scp->pk_normalization;
  TFset_parameters((scp->Om_matter)*scp->hubble*scp->hubble, scp->f_baryon, scp->Tcmb);
  tf_thisk =  TFfit_onek(k*scp->hubble, baryon_piece, cdm_piece);
  ans= Primordial_Matter_Power_Spectrum(scp,k)*pow(fabs(cdm_piece),2.); 
  return normalization*ans*pow(growth_factor,2);
}



// *******************************************************************
// *******************************************************************

real_prec  PowerSpectrum ::Linear_Matter_Power_Spectrum_DW(s_CosmologicalParameters *scp ,real_prec k){        /*this k comes in h/Mpc */
  real_prec G=exp(-0.5*pow(k/(scp->kstar),2));
  return Linear_Matter_Power_Spectrum(scp,k)*G+Linear_Matter_Power_Spectrum_NW(scp,k)*(1.-G);
}


// *******************************************************************
// *******************************************************************


real_prec  PowerSpectrum ::Linear_Matter_Power_Spectrum_G_NW(s_CosmologicalParameters *scp, real_prec k){        /*this k comes in h/Mpc */
  return Linear_Matter_Power_Spectrum_NW(scp, k)*exp(-0.5*pow(k/(scp->kstar),2));
}


// *******************************************************************
// *******************************************************************

void PowerSpectrum ::normalization(void *p, real_prec &nm){
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  if(true==scp->use_wiggles)nm=pow(scp->sigma8,2)/gsl_integration(fun,p, -7,7);
  else nm=pow(scp->sigma8,2)/gsl_integration(fun_nw,p, -7,7);
}



real_prec PowerSpectrum ::normalization(void *p){
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  real_prec nm;
  if(scp->use_wiggles==true)
    nm=pow(scp->sigma8,2)/gsl_integration(fun,p, -7,7);
  else
    nm=pow(scp->sigma8,2)/gsl_integration(fun_nw,p, -7,7);
  return nm;
}


// *******************************************************************
// *******************************************************************
gsl_real PowerSpectrum ::fun(gsl_real lk, void *p){
  real_prec k=pow(10,lk );
  PowerSpectrum  Ps;
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  real_prec tf_thisk, baryon_piece, cdm_piece;
  Ps.TFset_parameters((scp->Om_matter)*pow(scp->hubble, 2), scp->f_baryon, scp->Tcmb);
  tf_thisk =  Ps.TFfit_onek(k*scp->hubble, baryon_piece, cdm_piece);
  return static_cast<gsl_real>((log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*Ps.Primordial_Matter_Power_Spectrum(scp,k)*pow(tf_thisk,2)*pow(Ps.window(k,scp->RR),2));
}

// *******************************************************************
// *******************************************************************
gsl_real PowerSpectrum ::fun_nw(gsl_real lk, void *p){
 
  PowerSpectrum  Ps;
  real_prec k=pow(10,lk);
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  Ps.TFset_parameters((scp->Om_matter)*scp->hubble*scp->hubble, scp->f_baryon, scp->Tcmb);
  real_prec tfnw=Ps.TFnowiggles(scp->Om_cdm, scp->f_baryon,scp->hubble, scp->Tcmb, k);
  return static_cast<gsl_real>((log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*Ps.Primordial_Matter_Power_Spectrum(scp,k)*pow(tfnw,2)*pow(Ps.window(k,scp->RR),2));//*pow(scp->growth_factor,2);
}

// *******************************************************************
// *******************************************************************

// Compute Kastar at z=0
void  PowerSpectrum::kstar_integral(void *p, real_prec *ksta){
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  *ksta=pow((real_prec)((1./(6.*M_PI*M_PI))*gsl_integration(Power_Spectrum_i,(void *)scp,log10(scp->kmin_int),log10(scp->kmax_int))),(real_prec)-0.5);
}

// *******************************************************************
// *******************************************************************

gsl_real PowerSpectrum::Power_Spectrum_i(gsl_real lk, void *p){        /*this k comes in h/Mpc */
  PowerSpectrum Ps;
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  real_prec k=pow(10,lk);
  return static_cast<gsl_real>((log(10.0)*k)*Ps.Linear_Matter_Power_Spectrum_z(scp,k,1.0));
}

// *******************************************************************
// *******************************************************************
// WINDOW FUNCTION IN K SPACE
// *******************************************************************

real_prec PowerSpectrum ::window(real_prec k,real_prec R){
  /*Top hat window function, normalized suich that \int \dtx W(\xv) = 1*/
  real_prec y=k*R;
  return 3.*((sin(y)/(pow(y,3)))-(cos(y)/(pow(y,2))));
}
// *******************************************************************
real_prec PowerSpectrum ::windowg(real_prec k,real_prec R){return exp(-0.5*pow(k*R,2));}
// *******************************************************************

// *******************************************************************
// *******************************************************************
// *******************************************************************
// *******************************************************************

// NON-LINEAR MATTER POWER SPECTRUM WITH HALO_FIT
// *******************************************************************


real_prec PowerSpectrum::Non_Linear_Matter_Power_Spectrum_Halo_Fit(s_CosmologicalParameters *scp, real_prec k){
  real_prec ql, p, p_dw;
  halo_fit(k,(void *)scp,&ql,&p,&p_dw);
  return p;
}


// *******************************************************************
void PowerSpectrum::halo_fit(real_prec k, void *p,real_prec *ql,real_prec *pp,real_prec *pp_dw){
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  real_prec w_eos=scp->w_eos;
  real_prec fac   = (1./(2.0*M_PI*M_PI))*pow(k,3);
  real_prec dh,dh_dw, dql, dql_dw;
  Cosmology cCf;
  real_prec Omz=cCf.omega_matter(scp->cosmological_redshift,p);
  real_prec Omv=cCf.omega_dark_energy(scp->cosmological_redshift,p);

  real_prec y     = k/(scp->knl_hf);
  real_prec f1a    = pow(Omz,-0.0732);  /// Omega matter aca es al redhift, OJO
  real_prec f2a    = pow(Omz,-0.1423);
  real_prec f3a    = pow(Omz,+0.0725);

  real_prec f1b    = pow(Omz,-0.0307);  /// Omega matter aca es al redhift, OJO
  real_prec f2b    = pow(Omz,-0.0585);
  real_prec f3b    = pow(Omz,+0.0743);

  real_prec frac =Omv/(1.-Omz);
  real_prec f1=frac*f1b+(1-frac)*f1a;
  real_prec f2=frac*f2b+(1-frac)*f2a;
  real_prec f3=frac*f3b+(1-frac)*f3a;


  real_prec kmin  = log10(scp->kmin_int);
  real_prec kmax  = log10(scp->kmax_int);

  /*for the wiggled power spectrum*/
  real_prec rnl_hf=scp->rnl_hf;
  real_prec a, b, c, alpha, beta, gama, nu, mu;
  scp->aux_var3=rnl_hf;

  real_prec integration_aux4=gsl_integration(fun_aux_halo_fit4,(void *)scp, kmin,kmax);
  real_prec integration_aux2=gsl_integration(fun_aux_halo_fit2,(void *)scp, kmin,kmax);

  real_prec cc    = 4.0*pow(rnl_hf,2)*integration_aux4+4*pow(rnl_hf,4)*pow(integration_aux2,2);

  real_prec neff  =-3.0+2.0*pow(rnl_hf,2.)*integration_aux2;


  hf_aux(neff,cc,w_eos, Omv,&a,&b,&c,&alpha,&beta,&gama,&mu,&nu);
  real_prec dl  = fac*this->Linear_Matter_Power_Spectrum(scp,k);
  dql   = dl*(pow(1.+dl, beta)/(1.+alpha*dl))*exp(-y/4.-y*y/8.);
  real_prec dhp   = a*pow(y,3.0*f1)/(1.+b*pow(y,f2)+pow(c*f3*y,3-gama));
  dh    = dhp/(1.+(mu/y)+nu*pow(y,-2));

  *ql=dql/fac;
  *pp=(dh+dql)/fac;
  *pp_dw=(dh_dw+dql_dw)/fac;
}


// *******************************************************************
real_prec PowerSpectrum::Non_Linear_Matter_Power_Spectrum_Halo_Fit_z(s_CosmologicalParameters *scp, real_prec k, real_prec z, real_prec g, real_prec h4, real_prec h2, real_prec kln){
  real_prec ql, p;
  halo_fit_z(k,(void *)scp, z,g,h4, h2, kln, &ql,&p);
  return p;
}
// *******************************************************************
void PowerSpectrum::halo_fit_z(real_prec k, void *p, real_prec z, real_prec gf, real_prec h4, real_prec h2, real_prec knl, real_prec *ql,real_prec *pp){
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;

  real_prec fac   = (1./(2.0*M_PI*M_PI))*pow(k,3);

  Cosmology cCf;
  real_prec Omz=cCf.omega_matter(z,p);
  real_prec Omv=cCf.omega_dark_energy(z,p);

  real_prec rnl_hf= 1./knl;
  real_prec y     = k/knl;

  real_prec f1a    = pow(Omz,-0.0732);  /// Omega matter aca es al redhift, OJO
  real_prec f2a    = pow(Omz,-0.1423);
  real_prec f3a    = pow(Omz,+0.0725);

  real_prec f1b    = pow(Omz,-0.0307);  /// Omega matter aca es al redhift, OJO
  real_prec f2b    = pow(Omz,-0.0585);
  real_prec f3b    = pow(Omz,+0.0743);

  real_prec frac =  Omv/(1.-Omz); // =1 appears in the original paper of Takahashi. The original paper by R SMith suggest interpolation. CLASS uses frac like this:
  real_prec f1=frac*f1b+(1.-frac)*f1a;
  real_prec f2=frac*f2b+(1.-frac)*f2a;
  real_prec f3=frac*f3b+(1.-frac)*f3a;

  // Las cantidades h2 y h4 ya tienen incorporado el gf**2 respectivo, que viene desde el cl_model.

  real_prec cc    = 4.0*pow(rnl_hf,2)*h4+4.0*pow(rnl_hf,4)*pow(h2,2);
  real_prec neff  =-3.0+2.*pow(rnl_hf,2.)*h2;

  real_prec a, b, c, alpha, beta, gama, mu, nu;
  hf_aux(neff,cc,scp->w_eos, Omv,&a,&b,&c,&alpha,&beta,&gama,&mu,&nu);

  // Two-halo term Delta Q
  real_prec dl    = fac*this->Linear_Matter_Power_Spectrum_z(scp,k,gf);
  real_prec dql   = dl*(pow(1.+dl, beta)/(1.+alpha*dl))*exp(-y/4.0-y*y/8.0);

  // One halo term Delta H
  real_prec dhp   = a*pow(y,3.0*f1)/(1.+b*pow(y,f2)+pow(c*f3*y, 3.-gama));
  real_prec dh    = dhp/(1.+nu*pow(y,-2));  // mu = 0 for the fit of Takahashi

    *ql=dh/fac;
  *pp=(dh+dql)/fac;
//   *pp=(dh)/fac;
}

// *****************************************
// *****************************************

void PowerSpectrum::halo_fit_integrals(void *p, real_prec *inte4, real_prec *inte2 ){
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  real_prec kmin  = log10(scp->kmin_int);
  real_prec kmax  = log10(scp->kmax_int);

  int Ni = 100;
  real_prec integ4, integ2;
  if (true==scp->use_wiggles){
   integ4=gsl_integration3(Ni,fun_aux_halo_fit4,(void *)scp, kmin,kmax);
   integ2=gsl_integration3(Ni,fun_aux_halo_fit2,(void *)scp, kmin,kmax);
  }
  else
    {
      /*for the de-wigled power spectrum*/
      integ4=gsl_integration3(Ni,fun_aux_halo_fit4_dw,(void *)scp,kmin,kmax);
      integ2=gsl_integration3(Ni,fun_aux_halo_fit2_dw,(void *)scp, kmin,kmax);
    }
  *inte4=integ4;
  *inte2=integ2;
  
}



// *******************************************************************************
// *******************************************************************************

void PowerSpectrum::hf_aux(real_prec index, real_prec cc,real_prec weos, real_prec Omv,real_prec *a, real_prec *b, real_prec *c, real_prec *alpha,real_prec *beta,real_prec *gama,real_prec *mu, real_prec *nu){
  /*Funciones auxiliares del halo fit*/
  *a   = pow(10,1.5222+2.8553*index+2.3706*index*index+0.9903*index*index*index+0.2250*pow(index,4)-0.6038*cc+0.1749*Omv*(1+weos));
  *b   = pow(10,-0.5642+0.5864*index+0.5716*index*index-1.5474*cc+0.2279*Omv*(1+weos));
  *c   = pow(10,0.3698+2.0404*index+0.8161*index*index+0.5869*cc);
  *alpha= abs(6.0835+1.3373*index-0.1959*index*index-5.5274*cc);
  *beta = 2.0379-0.7354*index+0.3157*index*index+1.2490*pow(index,3)+0.3980*pow(index,4)-0.1682*cc;
  *gama = 0.1971-0.0843*index+0.8460*cc;
  *mu   = 0.0;         //pow(10,-3.5442+0.1908*index);
  *nu   = pow(10,5.2105+3.6902*index);
  return ;
}

// *******************************************************************************
// *******************************************************************************
void PowerSpectrum::nl_scales_halo_fit(void *p, real_prec *knl_hf,real_prec *rnl_hf, vector<real_prec>&rr, vector<real_prec>&sums, bool silence){
  if(silence)
    So.message_screen("Computing non linear scales for halo fit using Halo Fit by S03 and revised by Takahashi et al 2012");
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  int nr=rr.size();
  real_prec rinic;
  rinic = scp->Om_cdm<0.07? 1e-10 : 1e-3  ;
  real_prec rfinal=2e2;

  // We have an issue here. If we use pragma omp for, we get the numbers, but not that fast.
  // If we use pragma omp parallel for, the code executes jobs in parallel through different threads
  // but the information encoded in the structure is messed up, so no good answer. Issue Open.
  // Solution: I created another structure (struct s_aux, in Type_def*h) with two members, the cosmological_parameters structure
  // and a real_prec. Then I define an object of this type withon the loop so I get sure that every index
  // has its own strucuture defined, as if it were private.

  fill(sums.begin(), sums.end(), 0);

  omp_set_num_threads(1);
  vector<gsl_real> sums_aux (sums.size(),0);
  vector<gsl_real> rr_aux (sums.size(),0);

#ifdef SINGLE_PREC
  for(int i=0;i<sums.size();++i)
       sums_aux[i]=static_cast<gsl_real>(sums[i]);
  for(int i=0;i<sums.size();++i)
       rr_aux[i]=static_cast<gsl_real>(rr[i]);
#else
  sums_aux=sums;
#endif
//  #pragma omp parallel for
  for(int i=0;i<rr.size();++i){
    real_prec Rr= pow(10, log10(rinic)+i*log10(rfinal/rinic)/((real_prec)nr-1.));
    s_aux<PowerSpectrum>ssa;
    ssa.raux=Rr;
    ssa.scp_a=scp;
    rr_aux[nr-1-i]=Rr;
    // Compute sigma**2 (R,z=0)
    sums_aux[nr-1-i]= gsl_integration3(400,fun_aux_halo_fit,(void *)&ssa,log10(scp->kmin_int),log10(scp->kmax_int));
  }
   cout<<scp->kmin_int<<"  "<<scp->kmax_int<<endl;

  *rnl_hf=gsl_inter_new(sums_aux,rr_aux,num_1);
  *knl_hf=1./(*rnl_hf);

  silence=true;
  if(silence){
    s_aux<PowerSpectrum> ssb;
    ssb.raux=*rnl_hf;
    ssb.scp_a=scp;
    real_prec ss_check=gsl_integration3(400,fun_aux_halo_fit,(void *)&ssb,log10(scp->kmin_int),log10(scp->kmax_int));
    So.message_screen("Check: Sigma(r) at r = Rnl",ss_check);
    So.message_screen("Absolute error",100.0*abs(1- ss_check)," %");
  }


#ifdef SINGLE_PREC
  for(int i=0;i<sums.size();++i)
       sums[i]=static_cast<real_prec>(sums_aux[i]);
  for(int i=0;i<sums.size();++i)
       rr[i]=static_cast<real_prec>(rr_aux[i]);
#else
  sums=sums_aux;
  rr=rr_aux;
#endif


  return;
 }

  // *******************************************************************************
  // *******************************************************************************

 gsl_real PowerSpectrum::fun_aux_halo_fit(gsl_real lk, void *p){
   PowerSpectrum Ps;
   struct s_aux<PowerSpectrum> * saa = (struct s_aux<PowerSpectrum> *)p;
   s_CosmologicalParameters * scp = saa->scp_a ;
   real_prec k=pow(10,lk);
   real_prec r = saa->raux;
   real_prec power=Ps.Linear_Matter_Power_Spectrum(scp,k);
   real_prec ans=static_cast<gsl_real>((log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*power*exp(-pow(k*r,2)));
   return ans ;
 }

  // *******************************************************************************
  // *******************************************************************************

  gsl_real PowerSpectrum::fun_aux_halo_fit_dw(gsl_real lk, void *p){
    PowerSpectrum Ps;
    struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
    real_prec k=pow(10,lk);
    real_prec r = scp->aux_var3;
    return  static_cast<gsl_real>((log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*Ps.Linear_Matter_Power_Spectrum_DW(scp,k)*exp(-pow(k*r,2)));
  }

  // *******************************************************************************
  // *******************************************************************************

  gsl_real PowerSpectrum::fun_aux_halo_fit2(gsl_real lk, void *p){
    PowerSpectrum Ps;
    struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
    real_prec k=pow(10,lk);
    real_prec r = scp->aux_var3;
    return (log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*Ps.Linear_Matter_Power_Spectrum(scp,k)*exp(-pow(k*r,2))*pow(k,2);
  }

  // *******************************************************************************
  // *******************************************************************************

  gsl_real PowerSpectrum::fun_aux_halo_fit2_dw(gsl_real lk, void *p){
    PowerSpectrum Ps;
    struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
    real_prec k=pow(10,lk);
    real_prec r = scp->aux_var3;
    return  (log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*Ps.Linear_Matter_Power_Spectrum_DW(scp,k)*exp(-pow(k*r,2))*pow(k,2);
  }

  // *******************************************************************************
  // *******************************************************************************

  gsl_real PowerSpectrum::fun_aux_halo_fit4(gsl_real lk, void *p){
    PowerSpectrum Ps;
    struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
    real_prec k=pow(10,lk);
    real_prec r = scp->aux_var3;
    return (log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*Ps.Linear_Matter_Power_Spectrum(scp,k)*exp(-pow(k*r,2))*pow(k,2)*(1.-pow(r*k,2));
  }

  // *******************************************************************************
  // *******************************************************************************

  gsl_real PowerSpectrum::fun_aux_halo_fit4_dw(gsl_real lk, void *p){
    PowerSpectrum Ps;
    struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
    real_prec k=pow(10,lk);
    real_prec r = scp->aux_var3;
    real_prec ans=(log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*Ps.Linear_Matter_Power_Spectrum_DW(scp,k)*exp(-pow(k*r,2))*k*k*(1.-pow(r*k,2));
    return  static_cast<gsl_real>(ans);
  }

  // *****************************************************************************
  // *****************************************************************************
  // *****************************************************************************
  // GALAXY POWER SPECTRUM USING HALO MODEL
  // *****************************************************************************


  // *****************************************************************************
  // *****************************************************************************

  real_prec PowerSpectrum::Galaxy_power_spectrum_h1_ss(s_CosmologicalParameters *scp, real_prec k, real_prec z){

    A1 sA1;
    sA1.MASS=this->v_mass;
    sA1.MASS_FUNCTION=this->v_mass_function;
    sA1.MASS_BIAS=this->v_halo_mass_bias;
    sA1.s_cp=scp;
    sA1.aux_k=k;
    sA1.aux_z=z;
    return gsl_integration2(i_Galaxy_power_spectrum_h1_ss, (void *)&sA1,this->XX_Mass, this->WW_Mass);
  }

  // *****************************************************************************
  // *****************************************************************************

 gsl_real PowerSpectrum::i_Galaxy_power_spectrum_h1_ss(gsl_real m, void *p){
    struct A1 * sA1= (struct A1 *)p;
    s_CosmologicalParameters *scp = sA1->s_cp;
    HOD Shod;
    DensityProfiles Dp;
    real_prec M=pow(10,m);
    real_prec z=sA1->aux_z;
    real_prec k=sA1->aux_k;
    real_prec jacobian=(log(10.0)*M);
    real_prec uden=Dp.density_k(k,M,z, scp);
    return static_cast<gsl_real>(jacobian*gsl_inter_new(sA1->MASS, sA1->MASS_FUNCTION,m)*pow(Shod.SATELLITE(M, scp)*uden,2));

  }

  // *****************************************************************************
  // ****************************************************** ***********************
  real_prec PowerSpectrum::Galaxy_power_spectrum_h1_sc(s_CosmologicalParameters *scp, real_prec k, real_prec z){
    scp->aux_var4=k;
    scp->aux_var2=z;
    A1 sA1;
    sA1.MASS=this->v_mass;
    sA1.MASS_FUNCTION=this->v_mass_function;
    sA1.s_cp=scp;
    sA1.aux_k=k;
    sA1.aux_z=z;
    return gsl_integration2(i_Galaxy_power_spectrum_h1_sc, (void *)&sA1,this->XX_Mass, this->WW_Mass);
  }

  // *****************************************************************************
  gsl_real PowerSpectrum::i_Galaxy_power_spectrum_h1_sc(gsl_real m, void *p)
  {
    struct A1 * sA1= (struct A1 *)p;
    s_CosmologicalParameters *scp = sA1->s_cp;
    HOD Shod;
    DensityProfiles Dp;
//    MARKS Smark;
    real_prec M=pow(10,m);
    real_prec z=sA1->aux_z;
    real_prec k=sA1->aux_k;
    real_prec jacobian=(log(10.0)*M);
    real_prec uden=Dp.density_k(k,M,z,scp);
    return static_cast<gsl_real>(jacobian*2.0*Shod.CENTRAL(M, scp)*Shod.SATELLITE(M,scp)*gsl_inter_new(sA1->MASS, sA1->MASS_FUNCTION,m)*uden);
  }

  // *****************************************************************************
  // *****************************************************************************
  real_prec PowerSpectrum::Galaxy_matter_bias(s_CosmologicalParameters *scp, real_prec k, real_prec z){
    A1 sA1;
    sA1.MASS=this->v_mass;
    sA1.MASS_FUNCTION = this->v_mass_function;
    sA1.MASS_BIAS = this->v_halo_mass_bias;
    sA1.s_cp=scp;
    sA1.aux_k=k;
    sA1.aux_z=z;
    return gsl_integration2(i_Galaxy_matter_bias, (void *)&sA1,this->XX_Mass, this->WW_Mass);
  }

  // *****************************************************************************
  gsl_real PowerSpectrum::i_Galaxy_matter_bias(gsl_real m, void *p){
    struct A1 * sA1= (struct A1 *)p;
    s_CosmologicalParameters *scp = sA1->s_cp;
    HOD Shod;
    DensityProfiles Dp;
    real_prec M=pow(10,m);
    real_prec k=sA1->aux_k;
    real_prec z=sA1->aux_z;
    real_prec jacobian=(log(10.0)*M);
    real_prec uk=Dp.density_k(k,M,z,scp);
    return static_cast<gsl_real>(jacobian*(Shod.CENTRAL(M,scp)+Shod.SATELLITE(M,scp)*uk)*gsl_inter_new(sA1->MASS, sA1->MASS_BIAS,m)*gsl_inter_new(sA1->MASS, sA1->MASS_FUNCTION,m));
  }

  // *****************************************************************************

  real_prec PowerSpectrum::mean_galaxy_number_density(real_prec redshift, s_CosmologicalParameters *scp){
    scp->aux_var3=redshift;
    //Integrate from the value mmin_hod
    return gsl_integration2(i_mean_galaxy_number_density,(void *)scp,this->XX_Mass,this->WW_Mass);
  }

  // *****************************************************************************
  real_prec PowerSpectrum::mean_galaxy_number_density(s_CosmologicalParameters *scp){
    A1 sA1;
    sA1.MASS=this->v_mass;
    sA1.MASS_FUNCTION=this->v_mass_function;
    sA1.s_cp=scp;
    //Integrate from the value mmin_hod
    return gsl_integration2(i_mean_galaxy_number_density,(void *)&sA1,this->XX_Mass,this->WW_Mass);
  }

  // **********************************************************************************
  gsl_real PowerSpectrum::i_mean_galaxy_number_density(gsl_real m, void *p){
    struct A1 * sA1= (struct A1 *)p;
    s_CosmologicalParameters *scp = sA1->s_cp;
    HOD Shod;
    real_prec M=pow(10,m);
    real_prec jacobian=log(10)*M;
    return static_cast<gsl_real>(jacobian*gsl_inter_new(sA1->MASS,sA1->MASS_FUNCTION,m)*(Shod.SATELLITE(M,scp)+Shod.CENTRAL(M,scp)));
  }


  // *****************************************************************************
  // *****************************************************************************
  // *****************************************************************************
  // *****************************************************************************


  /* ------------------------ FITTING FORMULAE ROUTINES ----------------- */

  /* There are two routines here.  TFset_parameters() sets all the scalar
  parameters, while TFfit_onek() calculates the transfer function for a
  given wavenumber k.  TFfit_onek() may be called many times after a single
  call to TFset_parameters() */

  /* Global variables -- We've left many of the intermediate results as
  global variables in case you wish to access them, e.g. by declaring
  them as extern variables in your main program. */
  /* Note that all internal scales are in Mpc, without any Hubble constants! */


  /* Convenience from Numerical Recipes in C, 2nd edition */
  static real_prec sqrarg;
  #define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
  static real_prec cubearg;
  #define CUBE(a) ((cubearg=(a)) == 0.0 ? 0.0 : cubearg*cubearg*cubearg)
  static real_prec pow4arg;
  #define POW4(a) ((pow4arg=(a)) == 0.0 ? 0.0 : pow4arg*pow4arg*pow4arg*pow4arg)


  void PowerSpectrum ::TFset_parameters(real_prec omega0hh, real_prec f_baryon, real_prec Tcmb)
  /* Set all the scalars quantities for Eisenstein & Hu 1997 fitting formula */
  /* Input: omega0hh -- The density of CDM and baryons, in units of critical dens,
          multiplied by the square of the Hubble constant, in units
          of 100 km/s/Mpc */
  /* 	  f_baryon -- The fraction of baryons to CDM */

  /*        Tcmb -- The temperature of the CMB in Kelvin.  Tcmb<=0 forces use
              of the COBE value of  2.728 K. */
  /* Output: Nothing, but set many global variables used in TFfit_onek().
  You can access them yourself, if you want. */
  /* Note: Units are always Mpc, never h^-1 Mpc. */
  {

      real_prec z_drag_b1, z_drag_b2;
      real_prec alpha_c_a1, alpha_c_a2, beta_c_b1, beta_c_b2, alpha_b_G, y;

      if (f_baryon<=0.0 || omega0hh<=0.0) {
//         f_baryon=fabs(f_baryon);
 //        omega0hh=fabs(omega0hh);
         fprintf(stderr, "TFset_parameters(): Illegal input.");
         cout<<f_baryon<<"  "<<omega0hh<<endl;
       exit(1);
      }
      this->omhh = omega0hh;
      this->obhh = omhh*f_baryon;
      if (Tcmb<=0.0) Tcmb=2.728;	/* COBE FIRAS */
      this->theta_cmb = Tcmb/2.7;

      this->z_equality = 2.50e4*omhh/POW4(theta_cmb);  /* Really 1+z */
      this->k_equality = 0.0746*omhh/SQR(theta_cmb);

      z_drag_b1 = 0.313*pow(this->omhh,-0.419)*(1+0.607*pow(this->omhh,0.674));
      z_drag_b2 = 0.238*pow(this->omhh,0.223);
      this->z_drag = 1291.0*pow(omhh,0.251)/(1+0.659*pow(this->omhh,0.828))*(1.0+z_drag_b1*pow(this->obhh,z_drag_b2));

      this->R_drag = 31.5*obhh/POW4(theta_cmb)*(1000.0/(1.+z_drag));
      this->R_equality = 31.5*obhh/POW4(theta_cmb)*(1000./z_equality);

      this->sound_horizon = 2./3./k_equality*sqrt(6./R_equality)*log((sqrt(1+R_drag)+sqrt(R_drag+R_equality))/(1+sqrt(R_equality)));

      this->k_silk = 1.6*pow(obhh,0.52)*pow(omhh,0.73)*(1+pow(10.4*omhh,-0.95));

      alpha_c_a1 = pow(46.9*omhh,0.670)*(1+pow(32.1*omhh,-0.532));
      alpha_c_a2 = pow(12.0*omhh,0.424)*(1+pow(45.0*omhh,-0.582));
      this->alpha_c = pow(alpha_c_a1,-f_baryon)*pow(alpha_c_a2,-CUBE(f_baryon));

      beta_c_b1 = 0.944/(1.+pow(458.0*this->omhh,-0.708));
      beta_c_b2 = pow(0.395*omhh, -0.0266);
      this->beta_c = 1.0/(1+beta_c_b1*(pow(1-f_baryon, beta_c_b2)-1));

      y = this->z_equality/(1.+this->z_drag);
      alpha_b_G = y*(-6.*sqrt(1+y)+(2.+3.*y)*log((sqrt(1+y)+1)/(sqrt(1+y)-1)));
      alpha_b = 2.07*k_equality*sound_horizon*pow(1+R_drag,-0.75)*alpha_b_G;

      this->beta_node = 8.41*pow(omhh, 0.435);
      this->beta_b = 0.5+f_baryon+(3.-2.*f_baryon)*sqrt(pow(17.2*omhh,2.0)+1);

      this->k_peak = 2.5*3.14159*(1+0.217*omhh)/sound_horizon;
      this->sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1+10.0*pow(obhh,0.75));

      this->alpha_gamma = 1-0.328*log(431.0*omhh)*f_baryon + 0.38*log(22.3*omhh)*SQR(f_baryon);

      return;
  }


  // *****************************************************************************
  // *****************************************************************************

  real_prec PowerSpectrum ::TFfit_onek(real_prec k, real_prec &tf_baryon, real_prec &tf_cdm)
  /* Input: k -- Wavenumber at which to calculate transfer function, in Mpc^-1.
        *tf_baryon, *tf_cdm -- Input value not used; replaced on output if
                  the input was not NULL. */
  /* Output: Returns the value of the full transfer function fitting formula.
          This is the form given in Section 3 of Eisenstein & Hu (1997).
        *tf_baryon -- The baryonic contribution to the full fit.
        *tf_cdm -- The CDM contribution to the full fit. */
  /* Notes: Units are Mpc, not h^-1 Mpc. */
  {
      real_prec T_c_ln_beta, T_c_ln_nobeta, T_c_C_alpha, T_c_C_noalpha;
      real_prec q, xx, xx_tilde, q_eff;
      real_prec T_c_f, T_c, s_tilde, T_b_T0, T_b, f_baryon, T_full;
      real_prec T_0_L0, T_0_C0, T_0, gamma_eff;
      real_prec T_nowiggles_L0, T_nowiggles_C0, T_nowiggles;

      k = fabs(k);	/* Just define negative k as positive */
      if (k==0.0) 
      {
        if (&tf_baryon!=NULL) tf_baryon = 1.0;
        if (&tf_cdm!=NULL) tf_cdm = 1.0;
        return 1.0;
      }
      else
      {

      
      q = k/13.41/k_equality;
      xx = k*sound_horizon;

      T_c_ln_beta = log(2.718282+1.8*this->beta_c*q);
      T_c_ln_nobeta = log(2.718282+1.8*q);
      T_c_C_alpha = 14.2/this->alpha_c + 386.0/(1+69.9*pow(q,1.08));
      T_c_C_noalpha = 14.2 + 386.0/(1+69.9*pow(q,1.08));

      T_c_f = 1.0/(1.0+POW4(xx/5.4));
      T_c = T_c_f*T_c_ln_beta/(T_c_ln_beta+T_c_C_noalpha*SQR(q)) +
          (1-T_c_f)*T_c_ln_beta/(T_c_ln_beta+T_c_C_alpha*SQR(q));

      s_tilde = this->sound_horizon*pow(1+CUBE(beta_node/xx),-1./3.);
      xx_tilde = k*s_tilde;

      T_b_T0 = T_c_ln_nobeta/(T_c_ln_nobeta+T_c_C_noalpha*SQR(q));
      T_b = sin(xx_tilde)/(xx_tilde)*(T_b_T0/(1+SQR(xx/5.2))+alpha_b*exp(-pow(k/k_silk,1.4)))/(1.+CUBE(beta_b/xx));

      f_baryon = obhh/omhh;
      T_full = f_baryon*T_b+ (1.0-f_baryon)*T_c;


      /* Now to store these transfer functions */
      if (&tf_baryon!=NULL) tf_baryon = T_b;
      if (&tf_cdm!=NULL) tf_cdm = T_c;
      return T_full;
    }
  }

  /* ======================= Approximate forms =========================== */

  real_prec PowerSpectrum ::TFsound_horizon_fit(real_prec omega0, real_prec f_baryon, real_prec hubble)
  /* Input: omega0 -- CDM density, in units of critical density
        f_baryon -- Baryon fraction, the ratio of baryon to CDM density.
        hubble -- Hubble constant, in units of 100 km/s/Mpc
  /* Output: The approximate value of the sound horizon, in h^-1 Mpc. */
  /* Note: If you prefer to have the answer in  units of Mpc, use hubble -> 1
  and omega0 -> omega0*hubble^2. */
  {
      real_prec omhh, sound_horizon_fit_mpc;
      omhh = omega0*hubble*hubble;
      sound_horizon_fit_mpc = 44.5*log(9.83/omhh)/sqrt(1.+10.0*pow(omhh*f_baryon,0.75));
      return sound_horizon_fit_mpc*hubble;
  }

  // *****************************************************************************
  // *****************************************************************************

  real_prec PowerSpectrum ::TFk_peak(real_prec omega0, real_prec f_baryon, real_prec hubble)
  /* Input: omega0 -- CDM density, in units of critical density
        f_baryon -- Baryon fraction, the ratio of baryon to CDM density.
        hubble -- Hubble constant, in units of 100 km/s/Mpc
  /* Output: The approximate location of the first baryonic peak, in h Mpc^-1 */
  /* Note: If you prefer to have the answer in  units of Mpc^-1, use hubble -> 1
  and omega0 -> omega0*hubble^2. */
  {
      real_prec omhh, k_peak_mpc;
      omhh = omega0*hubble*hubble;
      k_peak_mpc = 2.5*3.14159*(1+0.217*omhh)/TFsound_horizon_fit(omhh,f_baryon,1.0);
      return k_peak_mpc/hubble;
  }


  // *****************************************************************************
  // *****************************************************************************

  real_prec PowerSpectrum ::TFnowiggles(real_prec omega0, real_prec f_baryon, real_prec hubble,
          real_prec Tcmb, real_prec k_hmpc)
  /* Input: omega0 -- CDM density, in units of critical density
        f_baryon -- Baryon fraction, the ratio of baryon to CDM density.
        hubble -- Hubble constant, in units of 100 km/s/Mpc
        Tcmb -- Temperature of the CMB in Kelvin; Tcmb<=0 forces use of
              COBE FIRAS value of 2.728 K
        k_hmpc -- Wavenumber in units of (h Mpc^-1). */
  /* Output: The value of an approximate transfer function that captures the
  non-oscillatory part of a partial baryon transfer function.  In other words,
  the baryon oscillations are left out, but the suppression of power below
  the sound horizon is included. See equations (30) and (31).  */
  /* Note: If you prefer to use wavenumbers in units of Mpc^-1, use hubble -> 1
  and omega0 -> omega0*hubble^2. */
  {
      real_prec k, omhh, theta_cmb, k_equality, q, xx, alpha_gamma, gamma_eff;
      real_prec q_eff, T_nowiggles_L0, T_nowiggles_C0;

      k = k_hmpc*hubble;	/* Convert to Mpc^-1 */
      omhh = omega0*hubble*hubble;
      if (Tcmb<=0.0) Tcmb=2.728;	/* COBE FIRAS */
      theta_cmb = Tcmb/2.7;

      k_equality = 0.0746*omhh/SQR(theta_cmb);
      q = k/13.41/k_equality;
      xx = k*TFsound_horizon_fit(omhh, f_baryon, 1.0);

      alpha_gamma = 1-0.328*log(431.0*omhh)*f_baryon + 0.38*log(22.3*omhh)*
          SQR(f_baryon);
      gamma_eff = omhh*(alpha_gamma+(1-alpha_gamma)/(1+POW4(0.43*xx)));
      q_eff = q*omhh/gamma_eff;

      T_nowiggles_L0 = log(2.0*2.718282+1.8*q_eff);
      T_nowiggles_C0 = 14.2 + 731.0/(1+62.5*q_eff);
      return T_nowiggles_L0/(T_nowiggles_L0+T_nowiggles_C0*SQR(q_eff));
  }


  // *******************************************************************
  // *******************************************************************
  // ======================= Zero Baryon Formula ======================

  real_prec PowerSpectrum ::TFzerobaryon(real_prec omega0, real_prec hubble, real_prec Tcmb, real_prec k_hmpc)
  /* Input: omega0 -- CDM density, in units of critical density
        hubble -- Hubble constant, in units of 100 km/s/Mpc
        Tcmb -- Temperature of the CMB in Kelvin; Tcmb<=0 forces use of
              COBE FIRAS value of 2.728 K
        k_hmpc -- Wavenumber in units of (h Mpc^-1). */
  /* Output: The value of the transfer function for a zero-baryon universe. */
  /* Note: If you prefer to use wavenumbers in units of Mpc^-1, use hubble -> 1
  and omega0 -> omega0*hubble^2. */
  {
      real_prec k, omhh, theta_cmb, k_equality, q, T_0_L0, T_0_C0;

      k = k_hmpc*hubble;	/* Convert to Mpc^-1 */
      omhh = omega0*hubble*hubble;
      if (Tcmb<=0.0) Tcmb=2.728;	/* COBE FIRAS */
      theta_cmb = Tcmb/2.7;

      k_equality = 0.0746*omhh/SQR(theta_cmb);
      q = k/13.41/k_equality;

      T_0_L0 = log(2.0*2.718282+1.8*q);
      T_0_C0 = 14.2 + 731.0/(1+62.5*q);
      return T_0_L0/(T_0_L0+T_0_C0*q*q);
  }
  // *******************************************************************
  // *******************************************************************
  // *******************************************************************
  // *******************************************************************
  // *******************************************************************
  // *******************************************************************
  // *******************************************************************
  // *******************************************************************
  // *******************************************************************
  // *******************************************************************
  // *******************************************************************
  // *******************************************************************
  // *******************************************************************
  // *******************************************************************
  // These functions are also in Statistics, but in order to use them in the
  // Cl analysis with omp, I had to copy them here


  real_prec PowerSpectrum::sigma_masa(real_prec m, real_prec z, s_CosmologicalParameters *scp){
    /* Sigma as a function of mass */
    scp->aux_var1 = m;  //log10M
    scp->aux_var2 = z;
    real_prec ans=gsl_integration2(isigma,(void *)scp,this->XX,this->WW);
    return sqrt(ans);
  }

  // ******************************************************************************

 gsl_real PowerSpectrum::isigma(gsl_real lk,void *p){  /* Integrand for sigma^2 */
    PowerSpectrum ps;
    Cosmology cf;
    real_prec k=pow(10,lk);
    s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
    real_prec m= s_cp->aux_var1;
    real_prec z= s_cp->aux_var2;
    real_prec M=pow(10,m);
    real_prec ans=(log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*ps.Linear_Matter_Power_Spectrum(s_cp,k)*pow(ps.window(k,cf.rr(M,z,p)),2);
    return ans;
  }

  // ******************************************************************************
  real_prec PowerSpectrum::mass_function_D(real_prec m, real_prec z,s_CosmologicalParameters  *scp)
  {
    Cosmology Cf;
    MASS_BIAS_FUNCTIONS mb;

    real_prec deltac = Cf.critical_density(z,scp);
    real_prec Smasa=gsl_inter_new(this->v_mass, this->v_sigma_mass,m);    //     sigma_masa(M,z,scp);
  //  real_prec Smasa=sigma_masa(M,z,scp);
    real_prec M=pow(10,m);
    real_prec nu=pow(deltac/Smasa,2);
    real_prec mean_matter_density=Cf.mean_matter_density(z,scp);
    real_prec com_density=(mean_matter_density)*pow(1+z,(real_prec)-3.0);

    A1 sA1;
    sA1.aux_m = m;
    sA1.s_cp=scp;
    sA1.aux_z = z;
    real_prec dsdr=gsl_integration2(dsigma_dR,(void *)&sA1, this->XX,this->WW);
    return (com_density/M)*(mb.mass_function(nu,z, (void *)scp))*pow(Smasa,-2)*(Cf.rr(M,z,(void *)scp)/(3.*M))*dsdr;
  }


// ******************************************************************************
  void PowerSpectrum::mass_function_M_Z(vector<real_prec>XZ, s_CosmologicalParameters  *scp){

   for(int i=0; i<this->v_mass.size();i++){
      this->v_mass[i]=(log10(scp->M_min_effective)+i*(log10(scp->M_max_effective)-log10(scp->M_min_effective))/this->v_mass.size());
   }

   for(int i=0; i<this->v_mass.size();i++){  //ESTO HACE EL SIGMA MASS DE HM SCALES INNECESARIO
      this->v_sigma_mass[i]=this->sigma_masa(this->v_mass[i],0.0,scp);
   }


   this->MASS_FUNCTION_M_Z.resize(XZ.size());
   for(int i=0;i<XZ.size();++i)this->MASS_FUNCTION_M_Z[i].resize(this->v_mass.size(),0);


   this->MASS_BIAS_M_Z.resize(XZ.size());
   for(int i=0;i<XZ.size();++i)this->MASS_BIAS_M_Z[i].resize(this->v_mass.size(),0);

   for(int i=0;i<XZ.size();++i)
    {
      real_prec gg=gsl_inter_new(scp->zv,scp->gv,XZ[i]);
      for(int j=0; j<this->v_mass.size();j++)
       {
         this->v_sigma_mass[j]*=gg;
         this->MASS_FUNCTION_M_Z[i][j]=mass_function_D(this->v_mass[j],XZ[i],scp);
         this->MASS_BIAS_M_Z[i][j]=bias(this->v_mass[j],XZ[i],scp);
       }
    }
  }



  // ******************************************************************************


  gsl_real PowerSpectrum::dsigma_dR(gsl_real lk, void *p){  /* integrand to calculate the derivative of sigma^2 with respect to  R(M)*/
    Cosmology Cf;
    PowerSpectrum ps;
    struct A1 * sA1= (struct A1 *)p;
    s_CosmologicalParameters *s_cp = sA1->s_cp;

    real_prec m= sA1->aux_m;
    real_prec M=pow(10,m);
    real_prec z= sA1->aux_z;
    real_prec k=pow(10,lk);
    real_prec r =Cf.rr(M,z,(void *)s_cp);
    real_prec y=r*k;
    return   (log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*ps.Linear_Matter_Power_Spectrum_z(s_cp,k, 1.0)*2.0*fabs(ps.window(k,r))*k*fabs(((3.*pow(y,-2))*sin(y)-9.*pow(y,-4)*(sin(y)-y*cos(y))));
  }


  real_prec PowerSpectrum::bias(real_prec m, real_prec z, s_CosmologicalParameters *scp){
    Cosmology Cf;
    MASS_BIAS_FUNCTIONS mb;
    real_prec Smasa=gsl_inter_new(this->v_mass, this->v_sigma_mass,m);    //     sigma_masa(M,z,scp);
    real_prec deltac = Cf.critical_density(z,scp);
    real_prec nu=pow(deltac/Smasa,2.);
    return mb.dm_h_bias(z,nu,(void *)scp);
  }

