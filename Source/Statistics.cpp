#ifndef __STATISTICS__
#define __STATISTICS__



# include "../Headers/PowerSpectrumTH.h"
# include "../Headers/Astrophysics.h"
# include "../Headers/ScalingRelations.h"
# include "../Headers/BiasFunctions.h"
# include "../Headers/Statistics.h"
# include "../Headers/HOD.h"
# include "../Headers/Marks.h"


void Statistics::compute_int_table_mass(real_prec M_min_integration, real_prec M_max_integration, int nss_k){
  wfd =  gsl_integration_glfixed_table_alloc (nss_k);
  this->WW_Mass.resize(nss_k);
  this->XX_Mass.resize(nss_k);
  gsl_get_GL_weights(static_cast<gsl_real>(log10(M_min_integration)),static_cast<gsl_real>(0.99*log10(M_max_integration)),this->wfd,this->XX_Mass,this->WW_Mass);

}


// aca tengo los mismos problemas que tenia en DENSITY_PROFILES
// con la cantidad de parametros que tengo.

// ******************************************************************************
// ******************************************************************************

real_prec  Statistics::As2sigma8(s_CosmologicalParameters *scp){
Cosmology cf;
  real_prec lkmin=log10(scp->kmin_int);
  real_prec lkmax=log10(scp->kmax_int);
  real_prec ans=sqrt(gsl_integration(iAs2sigma8,(void *)scp,lkmin,lkmax));
  return ans;
}

// ******************************************************************************
gsl_real Statistics::iAs2sigma8(gsl_real lk,void *p){  /* Integrand for sigma^2 */
  PowerSpectrum ps;
  Cosmology cf;
  real_prec k=pow(10,lk);
  s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  gsl_real ans=(log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*ps.Linear_Matter_Power_Spectrum(scp,k)*pow(ps.window(k,scp->RR),2);
  return ans;
}
// ******************************************************************************
// ******************************************************************************
real_prec Statistics::sigma_masa(real_prec m, real_prec z, s_CosmologicalParameters *scp){
  /* Sigma as a function of mass */

  real_prec lkmin=log10(1.09*scp->kmin_int);
  real_prec lkmax=log10(0.99*scp->kmax_int);
  scp->aux_var1 = m;  //log10M
  scp->aux_var2 = z;

  real_prec ans=gsl_integration3(60,isigma,(void *)scp,lkmin,lkmax);
  return sqrt(ans);
}

real_prec Statistics::sigma_masa(real_prec m, real_prec z, s_CosmologicalParameters *scp, vector<gsl_real>XX, vector<gsl_real>WW){
  /* Sigma as a function of mass */
  real_prec lkmin=log10(scp->kmin_int);
  real_prec lkmax=log10(scp->kmax_int);
  scp->aux_var1 = m;
  scp->aux_var2 = z;
  real_prec ans=gsl_integration2(isigma, (void *)scp, XX, WW);
  return sqrt(ans);
}


// ******************************************************************************
// ******************************************************************************

gsl_real Statistics::isigma(gsl_real lk,void *p){  /* Integrand for sigma^2 */
  PowerSpectrum ps;
  Cosmology cf;
  real_prec k=pow(10,lk);
  s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  ps.use_external_power=s_cp->use_external_power;
  real_prec m= s_cp->aux_var1;
  real_prec z= s_cp->aux_var2;
  real_prec M=pow(10,m);
  real_prec power;
  if(true==s_cp->use_external_power)
      power= pow(s_cp->growth_factor,2)*gsl_inter_new(s_cp->kvector_external, s_cp->power_external, k);
   else
      power=ps.Linear_Matter_Power_Spectrum(s_cp,k);
  real_prec ans=(log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*power*pow(ps.window(k,cf.rr(M,z,p)),2);
  return static_cast<gsl_real>(ans);
}

// // ******************************************************************************
// // ******************************************************************************
// real_prec Statistics::sigma_l(real_prec l, real_prec z, void *p){
//     /* Sigma as a function of X-ray luminosity for the REFLEX II sample, 
//        using R=Vmax**1/3, with Vmax as a function of Lx
//     */
//     struct two_pars tpar=  {l,z};
//     return sqrt(gsl_integration(isigma_l,&tpar,k_min,k_max));
//   }

// ******************************************************************************
// ******************************************************************************
real_prec Statistics::mass_function_D(real_prec m, real_prec z,s_CosmologicalParameters  *scp)
{
  Cosmology Cf;
  MASS_BIAS_FUNCTIONS mb;
  real_prec lkmin=log10(scp->kmin_int);
  real_prec lkmax=log10(scp->kmax_int);
 
  real_prec deltac=scp->critical_density;

  real_prec Smasa=gsl_inter_new(scp->v_mass, scp->v_sigma_mass,static_cast<gsl_real>(m));    //     sigma_masa(M,z,scp);
//  real_prec Smasa=sigma_masa(M,z,scp);
  real_prec M=pow(10,m);
  real_prec nu=pow(deltac/Smasa,2);
  real_prec com_density=(scp->mean_matter_density)*pow(1+z,(real_prec)-3.0);
  scp->aux_var1 = m;
  scp->aux_var2 = z;
  real_prec dsdr=gsl_integration3(60,dsigma_dR,(void *)scp, lkmin,lkmax);
  return (com_density/static_cast<double>(M))*(mb.mass_function(nu,z, (void *)scp))*pow(Smasa,-2)*(Cf.rr(M,z,(void *)scp)/(3.*static_cast<double>(M)))*dsdr;
}



// ******************************************************************************
real_prec Statistics::mass_function_D(real_prec m, real_prec z,s_CosmologicalParameters  *scp,vector<gsl_real>XX, vector<gsl_real>WW)
{
  Cosmology Cf;
  MASS_BIAS_FUNCTIONS mb;
  real_prec deltac=scp->critical_density;
  real_prec Smasa=gsl_inter_new(scp->v_mass, scp->v_sigma_mass,static_cast<gsl_real>(m));    //     sigma_masa(M,z,scp);
  real_prec M=pow(10,m);
  real_prec nu=pow(deltac/Smasa,2);
  real_prec com_density=(scp->mean_matter_density)*pow(1+z,(real_prec)-3.0);
  scp->aux_var1 = m;
  scp->aux_var2 = z;
  real_prec dsdr=gsl_integration2(dsigma_dR,(void *)scp,XX,WW);
  return (com_density/M)*(mb.mass_function(nu,z, (void *)scp))*pow(Smasa,-2)*(Cf.rr(M,z,(void *)scp)/(3.*M))*dsdr;
}



// // ******************************************************************************
// // ******************************************************************************



gsl_real Statistics::dsigma_dR(gsl_real lk, void *p){  /* integrand to calculate the derivative of sigma^2 with respect to  R(M)*/
  Cosmology Cf;
  PowerSpectrum ps;
  s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec m= s_cp->aux_var1;
  real_prec M=pow(10,m);
  real_prec z= s_cp->aux_var2;
  real_prec k=pow(10,lk);
  real_prec r =Cf.rr(M,z,p);
  real_prec y=r*k;
  real_prec power=1.0;
  if(true==s_cp->use_external_power)
      power= pow(s_cp->growth_factor,2)*gsl_inter_new(s_cp->kvector_external, s_cp->power_external, k);
   else
      power=ps.Linear_Matter_Power_Spectrum(s_cp,k);
  return  (log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*power*2.0*fabs(ps.window(k,r))*k*fabs(((3.*pow(y,-2))*sin(y)-9.*pow(y,-4)*(sin(y)-y*cos(y))));
}


// ******************************************************************************
// ******************************************************************************

real_prec Statistics::bias(real_prec m, real_prec z, s_CosmologicalParameters *scp){
  Cosmology Cf;
  MASS_BIAS_FUNCTIONS mb;  
  real_prec Smasa=gsl_inter_new(scp->v_mass, scp->v_sigma_mass,m);    //     sigma_masa(M,z,scp);
  real_prec nu=pow((scp->critical_density)/Smasa,2.);
  return mb.dm_h_bias(z,nu,(void *)scp);
}
// ******************************************************************************
// ******************************************************************************
// ******************************************************************************
// ******************************************************************************

real_prec Statistics::effective_halo_mass_bias(real_prec m, real_prec z, s_CosmologicalParameters *scp){
  // This computes the halo_mass bias of haloes with masses greater than M
  // note that since we are interpolating over quantities already computed, 
  // we do not need now explicitely the redshift
   scp->aux_var2=z;
  return gsl_integration(i_effective_halo_mass_bias,(void *)scp, m, log10(scp->M_max_effective))/gsl_integration(i_mass_function,(void *)scp, m, log10(scp->M_max_effective));
}

// ******************************************************************************
// ******************************************************************************

gsl_real Statistics::i_effective_halo_mass_bias(gsl_real m, void *p){
  s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec M=pow(10,m);
  real_prec jacobian= log(10.0)*M;
  vector<gsl_real> v_mass = s_cp->v_mass;
  vector<gsl_real> v_mass_function = s_cp->v_mass_function;
  vector<gsl_real> v_halo_mass_bias = s_cp->v_halo_mass_bias;
  return jacobian*gsl_inter_new(v_mass,v_mass_function,m)*gsl_inter_new(v_mass,v_halo_mass_bias,m);
}

// ******************************************************************************
// ******************************************************************************
real_prec Statistics::effective_halo_mean_number_density(real_prec m, real_prec z, s_CosmologicalParameters *scp){

  // This computes the mean number density of objects with masses greater than M

  scp->aux_var2=z;     
  // note that since we are interpolating over quantities already computed, 
  // we do not need now explicitely the redshift
  return gsl_integration(i_mass_function,(void *)scp, m, log10(scp->M_max_effective));
}

// ******************************************************************************
// ******************************************************************************

gsl_real Statistics::i_mass_function(gsl_real m, void *p){
  // I integrate wrt log10(M); foir this reason, this function expects m=lg10(M)
  // insted of M. If I integrate with respect to M, I should pass M to this function.
  // In that case, jacobian = 1
  real_prec M=pow(10,m);
  real_prec jacobian= log(10.0)*M;
  s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  vector<gsl_real> v_mass = s_cp->v_mass;
  vector<gsl_real> v_mass_function = s_cp->v_mass_function;
  return jacobian*gsl_inter_new(v_mass,v_mass_function,m);
}
// ******************************************************************************
// ******************************************************************************

// real_prec Statistics::cluster_mass_function(real_prec x, real_prec xt, real_prec z, s_CosmologicalParameters *scp){
//   return gsl_integration(i_cluster_mass_function_feedback,p,m_min,m_max);
// }

// real_prec Statistics::i_cluster_mass_function_feedback(real_prec m, void *p){

//   real_prec xt=(params->p1);
//   real_prec z=(params->p2);
//   ASTROPHYSICS ap;
//   return (log(10.0)*pow(10,m)*(M_reference))*(ap.Mobs_Mnb_distribution(m,0, p))*gsl_inter_pointers(mass_array_p,mass_function_p,nn,m);
// }







// real_prec Statistics::scale_dependent_bias(real_prec x, real_prec z, void *p){
//   real_prec r=x;
//   /*r en log(scale)*/
//   /*Tinker parametrization of scale-dependent bias in terms of the non-linear correlatin function*/
//   real_prec xi= gsl_inter_pointers(xRp,XI_NLp,nn,r);
//   //  real_prec f= ( r<*xRp ? 0 : pow(1.+1.17*xi,1.49)/pow(1.+0.69*xi,2.09));
//   real_prec f=  pow(1.+1.17*xi,1.49)/pow(1.+0.69*xi,2.09);
//   return f;
// }
// // ******************************************************************************
// // ******************************************************************************

// real_prec Statistics::delta_fof(real_prec x, real_prec z, void *p){
//   /*Calculo de Delta para halos FOF siguiendo los resultados de MOre et al 2011, en donde 
//     se asume un perfil NFW y una concentracion dada por la masa*/
//   real_prec al,fc,rv,c,rhos,rs,ans;
//   /*NOT ACCURATE. WARNING*/
//   //nfw_parameters(x,&fc,&rv,&mean_c,&rhos,&rs);
//   c=8.0;
//   fc   = log(1+(c))-((c)/(1.+(c)));
//   ans=3.*(0.652960)*pow(0.2,-3)*fc*(1+c)*pow(c,2);
//   return ans;
// }
// // ******************************************************************************
// // ******************************************************************************



// real_prec Statistics::mass_function_light_cone(){
//   real_prec ans;
//   real_prec rcmax,rcmin;
//   real_prec zmax=0.2;
//   real_prec zmin=0.001;
//   Cosmology cf;
//   real_prec mfz[nzpoints+1], zp[nzpoints+1];
//   for(int i=0;i<=nzpoints;i++)zp[i]=0;
//   for(int i=0;i<=nzpoints;i++)mfz[i]=0;  
//   mf_z_generator(x,1.1*zmax,zp,mfz);   /*THE FACTOR 1.1 PREVENTS NANS WHEN ONE NEEDS TO INTERPOLATE...*/
//   mfz_p=&mfz[0];
//   zp_p= &zp[0];
//   rcmax=cf.comoving_distance(zmax);
//   rcmin=cf.comoving_distance(zmin);
//   struct my_parameters par ={0,0,0};
//   return (3./(pow(rcmax,3)-pow(rcmin,3)))*gsl_integration(i_mass_function_m_z,&par,zmin,zmax);
// }
// // ******************************************************************************
// // ******************************************************************************

// real_prec Statistics::mass_function_prediction(real_prec x, real_prec z, void *p){  /*Mass function as a function of x=log10(M/masa_ns)*/
//   /*ACA LAS ENTRADAS SON MASAS; EL Z LO TOMAMOS EL PUNTERO*/
//   real_prec xmax=x;
//   real_prec xmin=z;
//   real_prec dxmax=(M_reference)*pow(10,xmax);
//   real_prec dxmin=(M_reference)*pow(10,xmin);
//   struct my_parameters par ={0,0,0};
//   return gsl_integration(imass_function_prediction,&par,xmin,xmax)/(dxmax-dxmin); // las masas estan en unidades de 10 a la 14  cuando mido n(m) de LBASICC
// }
// // ******************************************************************************
// // ******************************************************************************


// real_prec Statistics::lum_func_prediction(real_prec x, real_prec z, void *p){
//   real_prec l=x;
//   struct my_parameters par=  {l};
//   return gsl_integration(ilum_func_prediction,&par,m_min,m_max);
// }
// // ******************************************************************************
// // ******************************************************************************



// real_prec Statistics::occupancy_variance_numerator(real_prec x, real_prec z, void *p){   /*See Smith et al 2011*/
//   real_prec l=x;
//   struct my_parameters par=  {l};
//   return gsl_integration(i_occupancy_variance,&par,m_min,m_max);
// }

// // ******************************************************************************
// // ******************************************************************************


// real_prec Statistics::lum_func_reflex(real_prec x, real_prec z, void *p){
//   real_prec L=exp(x);                      /*L en unidades de 10^44 erg/s/h^2*/
//   return ncero*pow(L/lstar,alpha_lf+1)*Qexponential(qq,-L/lstar)/L; 

// // ******************************************************************************
// // ******************************************************************************
// real_prec Statistics::i_baryon_mass_function_reionization(real_prec x, real_prec z, void *p){
//   // struct my_parameters * params = (struct my_parameters *)p;
//   // real_prec xb=(params->p1);
//   // real_prec z =(params->p2);
//   ASTROPHYSICS ap;
//   real_prec ans=(log(10)*pow(10,x)*(M_reference))*ap.baryon_mass_virial_mass_distribution(x,z,p)*gsl_inter_pointers(mass_array_p,mass_function_p,nn,x);
//   return  ans;
// }




// real_prec Statistics::baryon_mass_function_reionization(real_prec x, real_prec l, real_prec z, void *p){
//   struct my_parameters par=  {x,z};
//   return gsl_integration(i_baryon_mass_function_reionization,&par,m_min,m_max);
// }



// ***********************************************************************************
// ***********************************************************************************


void Statistics::non_linear_scales(void *p, real_prec *knl,real_prec *Mnl,real_prec *rnl, real_prec *sigman){
  //cout<<"  "<<endl;
//  cout<<"Calculating non linear scales defined as the values of M such that Sigma(M)=1 "<<endl;
  Cosmology Cf;
  PowerSpectrum ps;
  s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec z=s_cp->aux_var1;
  vector<gsl_real> v_mass = s_cp->v_mass;
  vector<gsl_real> v_sigma_mass = s_cp->v_sigma_mass;


  real_prec normalization=s_cp->pk_normalization;
  real_prec den=s_cp->mean_matter_density*pow(1+z,-3);

  real_prec sig2,fac2,ms2,df2,dm2;
  /*Looking for mass scale where nu=1
    and scales where sigma(M)=1
    with Newton-Rhapson algorithm*/
  ms2=-2;
  //  if(ms1<m_min_interp){
  if(ms2<-4){
    cout<<"Potential error in function non_linear_scales:"<<endl;
    cout<<"minimum value of mass smaller than mn_it. CHECK!"<<endl;
  }
  int nr=7;
  for(int i=0;i<nr;i++){
    s_cp->aux_var1 = log10(pow(10,ms2)*M_reference);
    df2=-gsl_integration(dsigma_dR,(void *)s_cp,log10(s_cp->kmin_int),log10(s_cp->kmax_int));
    fac2=(4.*M_PI*den*pow(Cf.rr(pow(10,ms2)*M_reference,z,p),2))/(log(10.0)*(M_reference)*pow(10,ms2));
    sig2=gsl_inter_new(v_mass, v_sigma_mass, log10(pow(10,ms2)*M_reference));
    dm2=fac2*(pow(sig2,2)-1.0)/df2;
    ms2-=dm2;
  }
  *Mnl=pow(10,ms2)*(M_reference);
  *knl=pow((real_prec) 6*(pow(10,ms2)*(M_reference))/(M_PI*(Cf.mean_matter_density(z,p))),(real_prec)-1/3);
  *sigman=sigma_masa(log10(pow(10,ms2)*M_reference),z,s_cp); 
  *rnl=1./(*knl);
}
// ***********************************************************************************
// ***********************************************************************************
/*
real_prec Statistics::mean_galaxy_number_density(real_prec redshift, s_CosmologicalParameters *scp){
  scp->aux_var3=redshift;
  //Integrate from the value mmin_hod
  return gsl_integration3(30,i_mean_galaxy_number_density,(void *)scp,log10(scp->M_min_effective),log10(scp->M_max_effective));
}
*/
// ***********************************************************************************
real_prec Statistics::mean_galaxy_number_density(real_prec redshift, s_CosmologicalParameters *scp, vector<gsl_real>XX, vector<gsl_real>WW){
  scp->aux_var3=redshift;
  //Integrate from the value mmin_hod
  return gsl_integration2(i_mean_galaxy_number_density,(void *)scp,XX,WW);
}


real_prec Statistics::mean_galaxy_number_density(real_prec redshift, s_CosmologicalParameters *scp){
  scp->aux_var3=redshift;
  //Integrate from the value mmin_hod
  return gsl_integration2(i_mean_galaxy_number_density,(void *)scp,this->XX_Mass,this->WW_Mass);
}


// ***********************************************************************************
gsl_real  Statistics::i_mean_galaxy_number_density(gsl_real m, void *p){
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;  
  HOD Shod;
  MARKS Smark;
  real_prec M=pow(10,m);
  real_prec jacobian=log(10.)*M;
  return jacobian*gsl_inter_new(scp->v_mass,scp->v_mass_function,m)*(Shod.SATELLITE(M,scp)+Shod.CENTRAL(M,scp));
}
// ***********************************************************************************
// ***********************************************************************************
// ***********************************************************************************
// ***********************************************************************************

real_prec Statistics::mean_halo_number_density(s_CosmologicalParameters *scp){
  // Integrates the mass function wrt the mass to get the mean number density of objects
  return gsl_integration(i_mass_function,(void *)scp,log10(scp->M_min_effective),log10(scp->M_max_effective));

}

// ***********************************************************************************

#endif
