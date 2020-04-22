# include "../Headers/CosmologicalFunctions.h"
# include "../Headers/Astrophysics.h"
# include "../Headers/ScalingRelations.h"



double ASTROPHYSICS::mass2temp(double x, double z, void *p){
    ASTROPHYSICS ap;
    /////////////////////////////////////////////////////////////////
    /// Temperature in Kelvin as a function of circular velocity 
    /// (given in km/s). Mass in units of this code
    /////////////////////////////////////////////////////////////////

    return 0.5*mean_mol_weight*proton_mass*pow(vc(x,z,p),(double)2.0)/Boltzmann_constant; 
  }
  

double ASTROPHYSICS::cooling_function(double T, void*p){
  /*
    Cooling function as a function of the temperature in kelvins, taken from Robison and Silk 2000, after equation 5
    in units of kg m⁵ /s³ (these are J m³ /s)
    An isothermal gas density profile is assumed in order to give an analytial prescription
    x denotes temperature 
  */
  return (4e-36)*pow(T*Boltzmann_constant/kev,(double)-0.5);
}

  
double ASTROPHYSICS::filtering_mass(double z, void *p){
  /*Function implemented to introduce Reionization effects on the Halo Mass function, following Gnedin
    and the Fitting formula of Kravtstov et al 2004
    Filtering mass as a function of redshift in the Code units*/
  
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  double zr=7.0;
  double zo=8.0;
  double alpha_r=6.0;
  double fa;
  double a=1./(1.+z);
  double ar=1./(1.+zr);
  double ao=1./(1.+zo);
  double m_mol_w=0.59;
  double Jeans_Mass=2.5*sqrt(s_cp->Om_matter)*pow(m_mol_w,(double)-1.5)*(1.e11)/(M_reference);
  
  if(z>=zo){
    fa=(3.*a/((2.+alpha_r)*(5.+2.0*alpha_r)))*pow(a/ao,alpha_r);
  }
  if(z<zo && z>=zr ){
    fa=(3./a)*(pow(1.+zo,-2)*(pow(2.+alpha_r,(double)-1.0)-((2*pow(a*(1+zo),(double)-0.5))/(5.+2.*alpha_r)))+0.1*a*a-0.1*pow(1.+zo,(double)-2.0)*(5.0-4.*pow(a*(1.+zo),(double)-0.5)));
  }
  if(z<zr){
    fa=(3./a)*(pow((double)ao,(double)2)*(pow(2.+alpha_r,(double)-1.0)-((2*pow(a/ao,(double)-0.5))/(5.+2.*alpha_r)))+0.1*pow(ar,(double)2)*(5.0-4.*pow(a/ar,(double)-0.5))-0.1*pow(ao,(double)2)*(5.0-4.*pow(a/ao,(double)-0.5))+(1./3.)*(a*ar)-(1./3.)*pow(ar,(double)2)*(3.-2.*pow(a/ar,(double)-0.5)));
  }
  return Jeans_Mass*fa;
}


double ASTROPHYSICS::baryon_fraction(double x, double z, void *p){
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  double f_baryon=s_cp->Om_baryons/s_cp->Om_matter;
  double Mass=pow(10,x)*(M_reference);
  double ans= (f_baryon*(s_cp->Om_matter-s_cp->Om_baryons)/s_cp->Om_matter)/pow(1+0.26*filtering_mass(z,p)/Mass,(double)3.0);
  return ans;
}



double ASTROPHYSICS::cooling_radius(double x, double l, double z, void *p){
  /*
    Cooling radius determined from the assumptions on White and Frenk 1999, using equations 16 and 17
    in units of Mpc/h
    Up to now I'm not using this expression explicitely anywhere
  */
  
  Cosmology cf;
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;


  double time_scale=cf.halo_dynamical_time(z,p)*years_to_sec;  /*Halo dynamical time in seconds/h. This quantity is suggested by Crotton, see their discussion after eq. 4*/
  double fbary=baryon_fraction(x,z,p);
  double temp=0.5*mean_mol_weight*proton_mass*pow(vc(x,z,p),(double)2.0)/Boltzmann_constant;  /*en Kelvin*/
  double ans=(1.e-8)*(49.*fbary*s_cp->Om_baryons)*cooling_function(temp, p)/(192.*M_PI*Constants::Gravitational_constant*pow(proton_mass,(double)2)); 
  /*esta tiene unidades de km² / s : el factor 1e-8 viene del factor 1/(1000)⁵ que sale de pasar las uniades de metros a kilometros, ya que la function de cooling tiene unidades de metros*/
  //  ans=ans*((2./3.)/(Hubble*pow(1+(*red),3./2.)))/Mpc_to_km
  //  ans=sqrt(ans);
  
  ans=(ans/s_cp->hubble)*(time_scale); /*hago las unidades de tiempo en ans⁻¹ [sec/h]*/
  ans=sqrt(ans);               /*en unidades de km*/
  ans=ans*s_cp->hubble/Constants::Mpc_to_km;           /*en unidades de Mpc/h*/
  return ans;
}


double ASTROPHYSICS::vc(double x, double z, void *p){
  /*
    Circular velocity in a dark matter halo observed at the same time of collapse
    as a function of the comoving lagrangian radius rr(x)
    in a matter dominated universe. Taken from  Robinson and Silk 2000
    En unidades de km/s
  */
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  Cosmology cf;
  double ans=(1./sqrt(2))*(s_cp->Hubble)*cf.rr(x, z, p)*sqrt(s_cp->Om_matter)*sqrt(1+z)*pow((cf.density_contrast_top_hat(z,p))/(cf.omega_matter(z, p)),(double)1./6.);
  return ans; 
}



double ASTROPHYSICS::cooling_mass_rate(double x, double z, void *p){
  /*Cooling mass rate, equation 20 from White and Frenk 1991    in units of (Solar_masses/h) /(years/h)  */
  Cosmology cf;
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  //return (3./4.)*f_baryon*omegabaryon*Hubble*pow(1+(*red),3./2.)*cooling_radius(x)*pow(vc(x),2)*years_to_sec/Gravitational_constant/Solar_mass;
  double temp=mass2temp(x,z,p);
  /*Equation 8 of Robinson and Silk*/
  //  double time_scale=cf.age_universe();    /*In years/h, Original suggestion from White and Frenk*/
  double time_scale=cf.halo_dynamical_time(z,p); /*In years/h, Suggested by Crotton, see their discussion after eq. 4*/
  double fbary=baryon_fraction(x,z,p); 
  return 0.0045*s_cp->hubble*(s_cp->Om_matter*pow((double)fbary,(double)1.5)/sqrt(s_cp->Hubble*(Constants::years_to_sec/Constants::Mpc_to_km)*time_scale))*sqrt(cooling_function(temp,p)/1e-36)*pow(vc(x,z,p),2);
  /*el factor (years_to_sec/Mpc_to_km) convierte las unidades de Ho a km/s/MPc/h a Mpc/Year /Mpc/h*/
}

double ASTROPHYSICS::infall_mass_rate(double x, double z, void *p){
  //equation 7 of Robison and Silk, in Ms/h /(yrs/h)
  double fbary=baryon_fraction(x,z,p); 
  return 0.00024*Constants::sfr_efficiency*fbary*pow(vc(x,z,p),3);
}



double ASTROPHYSICS::star_formation_rate(double x, double l, double z, void *p){
  /*infall massa accretion rate in units Solar masses/h /years/h, equation 1 from White and Frenk 1990, 
    the term in the denominator is from their equation 23*/
  //  Minf=hubble*0.15*sfr_efficiency*f_baryon*pow(vc(x),3)/Gravitational_constant*years_to_sec/Solar_mass;
  //  return (vc(x)<min_circular_vel ? 0. : DMIN(infall_mass_rate(x),cooling_mass_rate(x))/(1.+0.02*pow(700./vc(x),2))); 
  return cooling_mass_rate(x,z,p); 
}



double ASTROPHYSICS::baryon_mass_virial_mass_distribution(double xxb, double xx, double z, void *p){
  double mf     =   pow(10,filtering_mass(z,p))*(M_reference);   /*Filtering mass*/
  double mv     =   pow(10,xx)*(M_reference);                   /*DMH mass*/
  double mg_mean=   mv*baryon_fraction(xxb,z,p);                /*Mean Mvir-Mg relation*/
  double mg     =   pow(10,xxb)*(M_reference);                  /*Gas mass*/
  double sigma_mass_reion=mf/(3.*mv);                       /*Intrinsic dispersion*/
  /*Log-normal distribution, Gnedin 2000:*/
  return (1./(sqrt(2.0*M_PI)*sigma_mass_reion))*exp(-0.5*pow(sigma_mass_reion,-2)*pow(log(mg/mg_mean)+0.5*pow(sigma_mass_reion,-2),2));
}


  
  
double ASTROPHYSICS::total_mass_nb_mass_relation(double x, double z, void *cp, void *ap){
  /*recibe la masa de simulaciones (nb)y de acuerdo a la fraccion de bariones
    como function de la masa total fg(M)g enera la relacion Mt(Mnb) mediante 
    el metodo Newton-Rhapson.
    solucionando Mt=Mnb*((1-fc)/(1-fg))  
    Retorna la masa total en unidades de Ms/h
    Lo que llamamos la masa de las simulaciones es la masa total Mnb=Mdm+fc*Mnb  = Mdm/(1-fc) donde fc
    es la fraccion cosmical de bariones. Esta masa sobre-estima la masa  real de un cluster con un modelo de fgas realista.
  */
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)cp;
  struct s_astrophysical_parameters * s_ap= (struct s_astrophysical_parameters *)ap;
  double f_baryon=s_cp->Om_baryons/s_cp->Om_matter;
  double xnb=x;
  double Mass_nb=pow(10,xnb)*(M_reference);
  double masita,Fg,fg,dFg;
  double alpha = (s_ap->A_gas)*(s_cp->f_baryon);
  double beta  = (s_ap->B_gas);
  masita=Mass_nb;
  for(int j=1;j<=20;j++){ /*Newton-Rhapson loop*/
    fg    = alpha*pow(masita/s_ap->mstar,beta);
    Fg    = masita-Mass_nb*(1-s_cp->f_baryon)/(1-fg);
    dFg   = 1.-(1-s_cp->f_baryon)*alpha*beta*pow(masita/s_ap->mstar,beta-1.0)*pow(1.-fg,-2.0)*Mass_nb/s_ap->mstar;
    masita-= Fg/dFg;
  }
  masita=(fg>=(s_cp->f_baryon)? Mass_nb: masita);
  return masita;
}


//double ASTROPHYSICS::Mobs_Mnb_distribution(double x_nb, double xt, double z, void *p){
  // double xx_nb=(M_reference)*pow(10,x_nb);
  // double xxt=(M_reference)*pow(10,xt);
  
  // ************************************************************************************
  // Commented ver si se deja o no
  // ************************************************************************************
  
  
  // /*Do MC realizations to characterize the PDF */
  // double pdf[3];
    // pdf[0]=0;
    // nb2m_mc(x_nb,pdf);
    // double xxt_mean=pdf[1];  /*mean of ln(Mtot(Mnb)) */
    // double ss      =pdf[2];
    // ss=sqrt(pow(scatter_Mobs_Mtot,2)+ss);

    // /*Use fit of sigma_mm from total_mass_nb*cpp and NR solver to avoid the MC realizations*/
    // // ASTROPHYSICS ap(x_nb,xt,z);
    // // double xxt_mean=log(ap.total_mass_nb_mass_relation());
    // // double ss=sqrt(pow(scatter_Mobs_Mtot,2)+pow(sigma_mm(x_nb),2));

    // // //Si queremos ignorar los efectos de fgas
    // // double xxt_mean=log(xx_nb);
    // // double ss=scatter_Mobs_Mtot;
    
    // double ans=(1./(sqrt(2.0*M_PI)*ss))*exp(-0.5*pow(ss,-2)*pow(log(xxt)-xxt_mean-log(bias_mass),2.0))/xxt; /*log-normal*/
    // //cout<<ss<<"  "<<xxt_mean<<endl;
    // return ans;
  // return 1.0;
  // }
 


double ASTROPHYSICS::mass_lum_distribution(double x, double l, double z, void *p, void *ap){
  /*PROBABILITY DISTRIBUTION FUNCTION FOR P(L|M)*/
  /*Log-normal distribution of mass-luminosities with intrinsec scatter. 
    Constant factors are not included for they calncel in the normalization
    x is the log10(mass/1e12)
    l is the natural log of the luminosity in units of  10^44 erg/sec h^-2 
    Nota importante: esta distribucion es log-normal, de modo que esta definida con p(ln L|M)d \ln L = (p(ln L|M)/L)dL */
  struct s_astrophysical_parameters * s_ap= (struct s_astrophysical_parameters *)ap;
  double y         = pow(10,x)*(M_reference)/(1e14); //ML RELATIONS RECIVE MASAS EN UNIDADES 1E14/h
  SCALING_RELATIONS ml;
  return      exp(-0.5*pow(l-ml.MOCKS_M2L(y,p),2)/(pow(s_ap->sigma_ln,2)))/exp(l);
}

double ASTROPHYSICS::mass_lum_distribution_errors(double x, double l, double z, void *cp, void *ap){
  /*this is the convolution of p(L|M) with another log-normal distribution with sigma for the flux errors
    Constant factors are not included for they cancel in the normalization
    x is the log10(mass/1e12)
    l is the natural log of the luminosity in units of  10^44 erg/sec h^-2 
    include missing flux correction in the factor missing flux */
  struct s_astrophysical_parameters * s_ap= (struct s_astrophysical_parameters *)ap;
  
  double y          =  pow(10,x)*(M_reference)/(1e14); //ML RELATIONS RECIVE MASAS EN UNIDADES 1E14/h
  SCALING_RELATIONS ml;
  l                 =  log(exp(l)/s_ap->missing_flux);
  return            exp(-0.5*pow(l-ml.MOCKS_M2L(y,cp),2)/pow(s_ap->sigma_red,2))/exp(l);
}



double ASTROPHYSICS::G(double x, double l, double z, void *p, void *ap){  
  /*Analytic integration of the function  mass_lum_distribution_errors for objects with lum greater than L. 
    No coloco los factores multiplicativos pues se cancelan */
  struct s_astrophysical_parameters * s_ap= (struct s_astrophysical_parameters *)ap;
  double factores     =  1./sqrt(2.*s_ap->sigma_red*s_ap->sigma_red);   
  double y=pow(10,x)*(M_reference)/(1e14);             //ML RELATIONS RECIVE MASAS EN UNIDADES 1E14/h
  SCALING_RELATIONS ml;
  return gsl_sf_erf(factores*(l-ml.MOCKS_M2L(y,p)));
}



double ASTROPHYSICS::G_Lmin_Lmax(double x, double l, double z, void *cp, void *ap){  
    /*Analytic integration of the function  mass_lum_distribution_errors for objects with lum greater than L. 
      and smaller than a fixed value Lmax.
      No coloco los factores multiplicativos pues se cancelan en las funciones de normalizacion*/
  struct s_astrophysical_parameters * s_ap= (struct s_astrophysical_parameters *)ap;
  double factores     =  1./sqrt(2.*s_ap->sigma_red*s_ap->sigma_red);   
  double y=pow(10,x)*(M_reference)/(1e14);             //ML RELATIONS RECIVE MASAS EN UNIDADES 1E14/h
  double lmin=l;
  double lmax=z;
  double interval=(exp(lmax)-exp(lmin));
  SCALING_RELATIONS ml;
  return 0.5*(gsl_sf_erf(factores*(lmax-ml.MOCKS_M2L(y,cp)))-gsl_sf_erf(factores*(lmin-ml.MOCKS_M2L(y,cp))))/interval;
  }



