#include "../Headers/CosmologicalFunctions.h"


// ***************************************************************************
// ***************************************************************************

gsl_real Cosmology::gint(gsl_real lscale_factor, void *p){   /*Used to obtain D(redshift)*/
  Cosmology cf;
  real_prec scale_factor=pow(10,lscale_factor);
  real_prec redshift=-1.0+1./scale_factor;
  return static_cast<gsl_real>((scale_factor*log(10.))*pow(scale_factor*cf.Hubble_function(redshift, p),-3));
}

// ***************************************************************************
// ***************************************************************************

gsl_real Cosmology::i_rs(gsl_real a, void *p){  /*Used to obtain rs(redshift)*/
  Cosmology cf;
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec redshift = 1./a -1.0;
  real_prec Om_photons= (s_cp->Om_radiation)/(1.+0.2271*s_cp->N_eff);
  real_prec Om_baryons= (s_cp->Om_baryons);
  real_prec Rg= (3.*Om_baryons)/(4.*Om_photons);
  return static_cast<gsl_real>(1./(a*a*cf.Hubble_function(redshift,p)*sqrt(1.+Rg*a)));
}
// ***************************************************************************
// ***************************************************************************

real_prec Cosmology::comoving_sound_horizon(real_prec redshift, void *p)
{
  return (Constants::speed_light/sqrt(3.))*static_cast<gsl_real>(gsl_integration(i_rs, p, 0, static_cast<gsl_real>(1./(1.+redshift))));
}

// ***************************************************************************
// ***************************************************************************

real_prec Cosmology::rr(real_prec M, real_prec z, void *p){
  Cosmology cf;
  /*este es r/a, i.e, rr is a comoving radius, 
    only used when the top-hat window function is needed
    In the end is redshift independent
    Output in Mpc/h
    here z is the mass in the units of this code, z=x*/
  return (1+z)*pow((real_prec)3.*M/(4.*M_PI*cf.mean_matter_density(z,p)),(real_prec)1./3.);
}

// ***************************************************************************
// ***************************************************************************

real_prec Cosmology::critical_density(real_prec z, void *p){
  Cosmology cf;
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec PI=acos(-1.0);
  return  (3./20.)*pow(12.*PI,2./3.)*(1.0-(5./936)*log(1.0+(1.-s_cp->Om_matter)/(s_cp->Om_matter*(1+z)))); /*parameterization of linear extrapol. overdensity*/
  
}

// ***************************************************************************
// ***************************************************************************

real_prec Cosmology::density_contrast_top_hat(real_prec z, void *p){
  // density contrast to be used in the fitting formulae of Tinker
  // so chech there whether this is relative to mean or background density
  return  (18*pow(M_PI,2)+82.0*omega_matter(z,p)-39.0*pow(omega_matter(z,p),2))/omega_matter(z,p); 
}


// ***************************************************************************
// ***************************************************************************


real_prec Cosmology::Hubble_function(real_prec redshift, void *p){
  Cosmology cf;
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  return s_cp->Hubble*sqrt(s_cp->Om_matter*pow(1+redshift,3)+s_cp->Om_radiation*pow(1+redshift,4)*(1+0.227*s_cp->N_eff)+s_cp->Om_k*pow(1+redshift,2)+s_cp->Om_vac*pow(1+redshift,3*(1.+s_cp->w_eos)));
}
// ***************************************************************************
// ***************************************************************************

gsl_real Cosmology::Hinv(gsl_real redshift, void *p){  /*Used to obtain r(redshift)*/
  Cosmology cf;
  return static_cast<gsl_real>(1./cf.Hubble_function(redshift, p));
} 
// ***************************************************************************

real_prec Cosmology::comoving_distance(real_prec redshift, void *p){
  Cosmology cf;
  real_prec ans=gsl_integration(Hinv, p, 0.0, redshift)*Constants::speed_light;
  return ans;
}

// ***************************************************************************
// ***************************************************************************


real_prec Cosmology::transverse_comoving_distance(real_prec redshift, void *p)
{
  Cosmology cf;
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec fac=s_cp->Hubble*sqrt(fabs(s_cp->Om_k));
  real_prec ans;

  if(s_cp->Om_k<0){
    ans= Constants::speed_light*sin(fac*cf.comoving_distance(redshift, p)/Constants::speed_light)/fac;
  }
  if(s_cp->Om_k==0){
    ans= cf.comoving_distance(redshift, p);
  }
  if(s_cp->Om_k>0){
    ans= Constants::speed_light*sinh(fac*cf.comoving_distance(redshift, p)/Constants::speed_light)/fac;
  }
  return ans;
}

// ***************************************************************************
// ***************************************************************************
real_prec Cosmology::inter_transverse_comoving_distance(real_prec redshift, void *p)
{
  //  Cosmology cf;
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec fac=s_cp->Hubble*sqrt(fabs(s_cp->Om_k));
  real_prec ans;
  real_prec cd=gsl_inter_new(s_cp->zv, s_cp->rv,redshift);
  if(s_cp->Om_k<0){
    ans= Constants::speed_light*sin(fac*cd/Constants::speed_light)/fac;
  }
  if(s_cp->Om_k==0){
    ans=cd;
  }
  if(s_cp->Om_k>0){
    ans= Constants::speed_light*sinh(fac*cd/Constants::speed_light)/fac;
  }
  return ans;
}


// ***************************************************************************
// ***************************************************************************

real_prec Cosmology::derivative_transverse_comoving_distance(real_prec redshift, void *p)
{
  Cosmology cf;
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;

  real_prec fac=s_cp->Hubble*sqrt(fabs(s_cp->Om_k));
  real_prec cd=gsl_inter_new(s_cp->zv, s_cp->rv, redshift);
  real_prec ans;
  if(s_cp->Om_k<0){
    ans= (Constants::speed_light/Hubble_function(redshift,p))*cos(fac*cd/Constants::speed_light)/fac;
  }
  if(s_cp->Om_k==0){
    ans= Constants::speed_light/Hubble_function(redshift, p);
  }
  if(s_cp->Om_k>0){
    ans= (Constants::speed_light/Hubble_function(redshift,p))*cosh(fac*cd/Constants::speed_light)/fac;
  }
  return ans;
}


// ***************************************************************************
// ***************************************************************************
real_prec Cosmology::comoving_angular_diameter_distance(real_prec redshift, void *p)
//This is proper
{
  return transverse_comoving_distance(redshift,p);
}


// ***************************************************************************
real_prec Cosmology::proper_angular_diameter_distance(real_prec redshift, void *p)
//This is proper
{
  return transverse_comoving_distance(redshift,p)/(1.+redshift);
}

real_prec Cosmology::inter_proper_angular_diameter_distance(real_prec redshift, void *p)
//This is proper
{
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  return gsl_inter_new(s_cp->zv, s_cp->trv, redshift)/(1.+redshift);
}

// ***************************************************************************
// ***************************************************************************

real_prec Cosmology::luminosity_distance(real_prec redshift, void *p)
{
  return transverse_comoving_distance(redshift,p)*(1.+redshift);
}

real_prec Cosmology::inter_luminosity_distance(real_prec redshift, void *p)
{
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  return gsl_inter_new(s_cp->zv, s_cp->trv, redshift)*(1.+redshift);
}

// ***************************************************************************
// ***************************************************************************


real_prec Cosmology::mean_matter_density(real_prec redshift, void *p)
/*Mean background density as a function of redshift in units (Ms/h)/(Mpc/h)^(-3) */
{
  Cosmology cf;
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec PI=acos(-1.0);
  return (3.*s_cp->Hubble*s_cp->Hubble/(8.*PI*Gravitational_constant))*s_cp->Om_matter*pow(1+redshift,3)*(Mpc_to_km/Solar_mass);  
}

// ***************************************************************************
// ***************************************************************************


real_prec Cosmology::age_universe(real_prec redshift, void *p)
/*Age of the universe in years / h*/
{
  return gsl_integration(Hinva,p,static_cast<gsl_real>(redshift),static_cast<gsl_real>(1e4))*(Constants::Mpc_to_km/Constants::years_to_sec);
}


// ***************************************************************************
gsl_real Cosmology::Hinva(gsl_real redshift, void *p){ /*Used to obtan Age(redshift)*/
   Cosmology cf;
   return 1./((1.+redshift)*cf.Hubble_function(static_cast<real_prec>(redshift),p));
 }


// ***************************************************************************
// ***************************************************************************
//  D(z) satisfying D(a)=a for a matter dom. time.
real_prec Cosmology::growth_factor(real_prec redshift, void *p){
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec a=1./(1.+redshift);
  real_prec ans1=gsl_integration(gint,p,-10,log10(a));
  return (5./2.)*(s_cp->Om_matter)*this->Hubble_function(redshift,p)*ans1*pow(s_cp->Hubble,2);
}

// ***************************************************************************
// **************************************************************************

real_prec Cosmology::growth_index(real_prec redshift, void *p){
  Cosmology cf;
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec om=s_cp->Om_matter*pow(1+redshift,3)/sqrt(s_cp->Om_matter*pow(1+redshift,3)+s_cp->Om_vac*pow(1+redshift,3*(1.+s_cp->w_eos)));
  return pow(om, (real_prec)0.55);
}

// ***************************************************************************
// ***************************************************************************

real_prec Cosmology::halo_dynamical_time(real_prec redshift, void *p){
  Cosmology cf;
  return (0.1/cf.Hubble_function(redshift, p))*(Mpc_to_km/years_to_sec);
}

// ***************************************************************************
// ***************************************************************************

real_prec Cosmology::omega_matter(real_prec redshift, void *p){
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  return s_cp->Om_matter*pow(1+redshift,3)*pow(Hubble_function(redshift, p)/s_cp->Hubble,-2);
}

real_prec Cosmology::omega_radiation(real_prec redshift, void *p){
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  return s_cp->Om_radiation*pow(1+redshift,4)*pow(Hubble_function(redshift, p)/s_cp->Hubble,-2);
}


real_prec Cosmology::omega_curvature(real_prec redshift, void *p){
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  return s_cp->Om_k*pow(1+redshift,2)*pow(Hubble_function(redshift, p)/s_cp->Hubble,-2);
}

real_prec Cosmology::omega_dark_energy(real_prec redshift, void *p){
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  return s_cp->Om_vac*pow(1+redshift,1+s_cp->w_eos)*pow(Hubble_function(redshift, p)/s_cp->Hubble,-2);
}

// ***************************************************************************
// ***************************************************************************

real_prec Cosmology:: Distance_Modulus(real_prec redshift, void *p){
  // ***************************
  // DISTANCE MODULUS
  // ***************************
  return 25.+5.0*log10(luminosity_distance(redshift, p));
}


real_prec Cosmology:: inter_Distance_Modulus(real_prec redshift, void *p){
  // ***************************
  // DISTANCE MODULUS
  // ***************************
  return 25.+5.0*log10(inter_luminosity_distance(redshift, p));
}


// ***************************************************************************
// ***************************************************************************
real_prec Cosmology:: K_correction(real_prec redshift, void *p){
  // ***************************
  // K-CORRECTION
  // ***************************
  //  return -6.0*log10(1+redshift);
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
//  return s_cp->k_index*log10(1+redshift);

    real_prec a0=0.01529;
    real_prec a1=3.2838;
    real_prec a2=-13.3196;
    real_prec a3=120.51;
    real_prec a4=-391.612;
    real_prec a5=395.23;
    return a0+a1*redshift+a2*pow(redshift,2)+a3*pow(redshift,3)+a4*pow(redshift,4)+a5*pow(redshift,5);
}

// ***************************************************************************
// ***************************************************************************

real_prec Cosmology:: K_correction(real_prec redshift, real_prec color1, real_prec color2,  void *p){
  // ***************************
  // K-CORRECTION
  // ***************************
  //  return -6.0*log10(1+redshift);
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec a0=s_cp->K_index_a;
  real_prec a1=s_cp->K_index_b;
  real_prec a2=s_cp->K_index_c;
  real_prec a3=s_cp->K_index_d;
  real_prec a4=s_cp->K_index_e;
  real_prec a5=s_cp->K_index_f;
  real_prec a6=s_cp->K_index_g;
  real_prec a7=s_cp->K_index_h;
  real_prec a8=s_cp->K_index_i;
  real_prec a9=s_cp->K_index_j;

//  real_prec lred=log10(1.0+redshift);
  real_prec lred=redshift;
  real_prec color_poly= a0+a1*color1 + a2*color1*color1 + a3*pow(color1,3)+a4*pow(color1,4);
  real_prec z_poly =    a5*lred + a6*lred*lred + a7*pow(lred,3) + a8*pow(lred,4)+a9*pow(lred,5);
  return color_poly*z_poly;
}


// ***************************************************************************
real_prec Cosmology:: dK_correction_dz(real_prec redshift, void *p){
  // ***************************
  // K-CORRECTION
  // ***************************
  //  return -6.0*log10(1+redshift);
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
//  return s_cp->k_index*log10(1+redshift);

    real_prec a1=3.2838;
    real_prec a2=-13.3196;
    real_prec a3=120.51;
    real_prec a4=-391.612;
    real_prec a5=395.23;
    return s_cp->use_e_correction ? a1+2.*a2*redshift+3.*a3*pow(redshift,2)+4.*a4*pow(redshift,3)+5.*a5*pow(redshift,4): 0. ;

}


// ***************************************************************************
// ***************************************************************************

real_prec Cosmology:: e_correction(real_prec redshift, void *p){
  // ***************************
  // e-CORRECTION
  // ***************************
  //return 3.04*redshift;
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  return s_cp->use_e_correction ?  s_cp->e_index_a*pow(redshift,2)+s_cp->e_index_b*redshift+s_cp->e_index_c  :   0.0 ;

}
real_prec Cosmology:: de_correction_dz(real_prec redshift, void *p){
  // ***************************
  // e-CORRECTION
  // ***************************
  //return 3.04*redshift;
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  return s_cp->use_e_correction ? 2.*s_cp->e_index_a*redshift+s_cp->e_index_b:  0. ;

}

// ***************************************************************************
// ***************************************************************************



real_prec Cosmology::zmax_old(void *p){
  // ****************************************************************************
  // Given the limiting magnitude of the sample and the absolute
  // magnitude of the object, find the zmax to which this object could have been
  // observed.
  // ****************************************************************************
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec z_ini=0.01;
  real_prec z=z_ini;

  vector<real_prec> zmm;
  int i=0;
  do{
    i++;
    if(z<0)std::cerr<<"Negative reds    hift present"<<endl;
    real_prec tcd = gsl_inter_new(s_cp->zv, s_cp->trv, z);
    real_prec F  = s_cp->Mabs-s_cp->mlim+25.0+5.0*log10(tcd*(1+z))+this->e_correction(z,p)+this->K_correction(z,p) ;
    real_prec dF=  (5.0/log(10.0))*( 1./(1.+z) + derivative_transverse_comoving_distance(z,p)/tcd)+ this->dK_correction_dz(z,p)  +  this->de_correction_dz(z,p);
    z -= F/dF;
    //z=fabs(z);
    cout<<s_cp->Mabs<<"  "<<z<<"  "<<F<<"  "<<dF<<endl;
    zmm.push_back(z);

  }while(fabs((zmm[zmm.size()-1]-zmm[zmm.size()-2])/zmm[zmm.size()-2])>ERROR_LIMIT_NR);
  return z;
}


real_prec Cosmology::zmax_old(real_prec z_ini, void *p){
  // ****************************************************************************
  // Given the limiting magnitude of the sample and the absolute
  // magnitude of the object, find the zmax to which this object could have been
  // observed.
  // ****************************************************************************
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec Mabs=s_cp->Mabs;
  real_prec mlim=s_cp->mlim;
  real_prec z=z_ini;

   vector<real_prec> zmm;
  int i=0;
  do{
    i++;
    real_prec tcd = gsl_inter_new(s_cp->zv, s_cp->trv, z);
    real_prec F  = Mabs-mlim+25.0+5.0*log10(tcd*(1+z))+K_correction(z,p)+e_correction(z,p);
    real_prec dF=  (5.0/log(10.0))*( 1./(1.+z) + derivative_transverse_comoving_distance(z,p)/tcd)+ this->dK_correction_dz(z,p)  +  this->de_correction_dz(z,p);
    z -= F/dF;
    zmm.push_back(z);
  }while(fabs(zmm[zmm.size()-1]-zmm[zmm.size()-2])>ERROR_LIMIT_NR);

return z;
}


// ***************************************************************************
// ***************************************************************************

real_prec Cosmology::zmax(void *p){
  // ****************************************************************************
  // Given the limiting magnitude of the sample and the absolute
  // magnitude of the object, find the zmax to which this object could have been
  // observed. 
  // applied.
  // ****************************************************************************

    
  int status;
  int iter = 0, max_iter = 50;
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  real_prec x0, x = 0.005;
  gsl_function_fdf FDF;
  
  FDF.f = &this->Froot;
  FDF.df = &this->dFroot;
  FDF.fdf = &this->F_dF;
  FDF.params = p;

  T = gsl_root_fdfsolver_secant;
  s = gsl_root_fdfsolver_alloc (T);
  gsl_root_fdfsolver_set (s, &FDF, x);

  do{
    iter++;
    status = gsl_root_fdfsolver_iterate (s);
    x0 = x;
    x = gsl_root_fdfsolver_root (s);
    status = gsl_root_test_delta (x, x0, 0, 1e-3);
  }
  while (status == GSL_CONTINUE && iter < max_iter);
  
  gsl_root_fdfsolver_free (s);
  return x;
}

// ***************************************************************************
// **************************************************************************

gsl_real Cosmology::Froot(gsl_real z, void *p){
  Cosmology cf;
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  gsl_real tcd=gsl_inter_new(s_cp->zv, s_cp->trv,z);
  gsl_real ans=s_cp->Mabs-s_cp->mlim+25.0+5.0*log10(tcd*(1+z))+cf.K_correction(z,p)+cf.e_correction(z,p);
  return ans;
}

// ***************************************************************************
// **************************************************************************

gsl_real Cosmology::dFroot(gsl_real z, void *p){
  Cosmology cf;
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec tcd=gsl_inter_new(s_cp->zv, s_cp->trv,z);
  real_prec ans=(5.0/log(10.0))*( 1./(1.+z) + cf.derivative_transverse_comoving_distance(z,p)/tcd)+cf.dK_correction_dz(z,p)  +  cf.de_correction_dz(z,p);
  return ans;
}

// ***************************************************************************
// ***************************************************************************

void Cosmology::F_dF(gsl_real z, void *p, gsl_real *y, gsl_real *dy){
  Cosmology cf;
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec tcd=static_cast<real_prec>(gsl_inter_new(s_cp->zv, s_cp->trv, z));
  *y=s_cp->Mabs-s_cp->mlim+25.0+5.0*log10(tcd*(1+z))+cf.K_correction(static_cast<real_prec>(z),p)+cf.e_correction(static_cast<real_prec>(z),p);
  *dy= (5.0/log(10.0))*( 1./(1.+z) + cf.derivative_transverse_comoving_distance(static_cast<real_prec>(z),p)/tcd)+ cf.dK_correction_dz(static_cast<real_prec>(z),p)  +  cf.de_correction_dz(static_cast<real_prec>(z),p);
}
// ***************************************************************************
// ***************************************************************************
void Cosmology::free_gsl_table(){
  gsl_integration_glfixed_table_free(this->wf);
}

// ***************************************************************************
// ***************************************************************************

real_prec Cosmology::gsl_cosmo_integration(gsl_real (*function)(gsl_real, void *) ,void *p,gsl_real LowLimit,gsl_real UpLimit){
  gsl_function F;
  F.params   = p;  
  F.function = function;
  return gsl_integration_glfixed(&F,LowLimit,UpLimit,this->wf);
}    
// ***************************************************************************
// ***************************************************************************
void Cosmology::check_cosmo_pars(s_CosmologicalParameters *scp){
   scp->Om_matter = scp->Om_cdm+scp->Om_baryons;
   scp->f_baryon  = scp->Om_baryons/scp->Om_matter;
   scp->Om_vac    = 1.- scp->Om_matter-scp->Om_radiation- scp->Om_k;
}

// ***************************************************************************
// ***************************************************************************

void Cosmology::Comoving_distance_tabulated(real_prec z_min, real_prec z_max, void *p, vector<gsl_real> & zz, vector<gsl_real> & rc)
{
  ////////////////////////////////////////////////////////
  // Tabulation of the comoving distance.
  ////////////////////////////////////////////////////////
  for (size_t i=0; i<zz.size(); i++) {
    zz[i]=static_cast<gsl_real>(z_min+i*(z_max-z_min)/(zz.size()-1.0));
    rc[i]=static_cast<gsl_real>(comoving_distance(zz[i],p));
  }  
}
