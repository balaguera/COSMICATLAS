# include "../Headers/CosmoLib.h"

using namespace std;
void message(string mess){
  cout<<BOLDRED<<mess<<RESET<<endl;
}

void welcome_message(){
 time_t rawtime;
 time (&rawtime);

 cout<<CYAN<<"***********************************************************************************"<<endl;
 cout<<"COSMOLIB: Some cosmology-related numbers "<<endl;
 cout<<"VERSION 1.2"<<endl;
 cout<<"\t"<<endl;
 cout<<"Starting time "<<ctime (&rawtime);
 cout<<"***********************************************************************************"<<RESET<<endl;
}



void welcome_read_chains(){
  time_t rawtime;
 time (&rawtime);
 cout<<"***********************************************************************************"<<endl;
 cout<<"Reading MCMChains"<<endl;
 cout<<"\t"<<endl;
 cout<<"Starting time "<<ctime (&rawtime);
 cout<<"***********************************************************************************"<<endl;
}


void message_screen(string s, real_prec v, string units)
{
  cout<<"Derived quantity"<<endl;
  cout<<CYAN;
  cout<<s<<" = "<<v<<" "<<units<<RESET<<endl;
}

void message_screen(string s, real_prec v)
{
  cout<<"Derived quantity"<<endl;
  cout<<CYAN<<s<<" = "<<v<<RESET<<endl;
}
void message_screen_ini(string s, real_prec v)
{
  cout<<"Input parameter"<<endl;
  cout<<CYAN<<s<<" = "<<v<<RESET<<endl;
}


void welcome_chains(){
 time_t rawtime;
 time (&rawtime);
 cout<<"***********************************************************************************"<<endl;
 cout<<"MCMChains"<<endl;
 cout<<"\t"<<endl;
 cout<<"Starting time "<<ctime (&rawtime);
 cout<<"***********************************************************************************"<<endl;
}


// **************************************************************************
// **************************************************************************

void done(string text){
  cout<<"...leaving "<<text<<endl;
  cout<<"\t"<<endl;
}

// **************************************************************************
// **************************************************************************

void enter(string text){
  cout<<"\t"<<endl;
  cout<<"Going to:  "<<text<<endl;
}

// **************************************************************************
// **************************************************************************


void ending_message(){
  cout<<"    "<<endl;
  time_t rawtime;
  time ( &rawtime );
  cout<<"Ending time:"<<ctime (&rawtime)<<endl;
  cout<<"***********************************************************************************"<<endl;
}

// **************************************************************************
// **************************************************************************
void CosmoLib::set_par_file(string par_file)
{
  ParametersCosmolib param (par_file);
  this->params=param;
  Statistics Csa(this->params._M_min_mf(), params._M_max_mf(), params._n_points_mf());
  this->Cs=Csa;
  PowerSpectrum Psa(params._k_min_integration() ,params._k_max_integration() ,params._n_points_dp_k(),10, params._M_min_mf(),params._M_max_mf(),params._n_points_mf());
  this->Ps=Psa;
}

// **************************************************************************
// **************************************************************************

void write_cosmo_parameters(void *p)
{
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  ScreenOutput So;
  So.message_screen("COSMOLOGICAL PARAMETERS");
  So.message_screen("Omega Matter", s_cp->Om_matter);
  So.message_screen("Omega Vacuum", s_cp->Om_vac);
  So.message_screen("Omega Baryons", s_cp->Om_baryons);
  So.message_screen("Omega curvature", s_cp->Om_k);
  So.message_screen("Hubble Ho", s_cp->Hubble);
  So.message_screen("Spectral index", s_cp->n_s);
  So.message_screen("Sigma 8", s_cp->sigma8);
  So.message_screen("Primordial amplitude", s_cp->A_s);
  So.message_screen("Running index", s_cp->alpha_s);
  So.message_screen("Mean CMB temperature", s_cp->Tcmb);
}

// **************************************************************************
// **************************************************************************

// **************************************************************************
// **************************************************************************


void comp_time(time_t start, long full, long step){

  real_prec fraction=100.0*((real_prec)(step+1))/((real_prec)full);
  time_t end;
  time(&end);
  real_prec lapse=difftime(end,start);
  if(lapse<=60)cout<<"\r "<<fraction<<" % completed. Time elapsed: "<<lapse<<"  secs \r";cout.flush();
  if(lapse>60)cout<<"\r "<<fraction<<" % completed. Time elapsed: "<<lapse/60.<<"  mins \r";cout.flush();
  if(lapse>3600)cout<<"\r "<<fraction<<" % completed. Time elapsed: "<<lapse/3600.<<"  hrs \r";cout.flush();
  if (fraction==25) cout <<"\r  ..25% completed \r";cout.flush();
  if (fraction==50) cout <<"\r  ..50% completed \r";cout.flush();
  if (fraction==75) cout <<"\r  ..75% completed \r";cout.flush();
  if (fraction==100) cout<<"\r ..100% completed \r";cout.flush();

}


void comp_time_MH(time_t start, long full, long step, real_prec H){
  std::cout<<RED;
  real_prec fraction=100.0*((real_prec)(step+1))/((real_prec)full);
  time_t end;
  time(&end);
  real_prec lapse=difftime(end,start);
  if(lapse<=60)cout<<"\r "<<fraction<<" % completed. Time elapsed: "<<lapse<<"  secs.  H(q,p) = "<<H<<"\r";cout.flush();
  if(lapse>60)cout<<"\r "<<fraction<<" % completed. Time elapsed: "<<lapse/60<<"  secs.  H(q,p) = "<<H<<"\r";cout.flush();
  if(lapse>3600)cout<<"\r "<<fraction<<" % completed. Time elapsed: "<<lapse/3600<<"  secs.  H(q,p) = "<<H<<"\r";cout.flush();
  std::cout<<RESET;


  if (fraction==25) cout <<"\r  ..25% completed \r";cout.flush();
  if (fraction==50) cout <<"\r  ..50% completed \r";cout.flush();
  if (fraction==75) cout <<"\r  ..75% completed \r";cout.flush();
  if (fraction==100) cout<<"\r ..100% completed \r";cout.flush();

}


// **********************************************************************************
// **********************************************************************************
// **********************************************************************************
// **********************************************************************************
// **********************************************************************************
void CosmoLib::feed_cosmo_struct(){

    s_cosmo_par.Om_radiation=params._om_radiation();
    s_cosmo_par.Om_cdm=params._om_cdm();
    s_cosmo_par.Om_baryons=params._om_baryons();
  //  s_cosmo_par.Om_vac=params._om_vac();
    s_cosmo_par.Om_k=params._om_k();
    s_cosmo_par.Hubble=params._Hubble();
    s_cosmo_par.hubble=params._hubble();
    s_cosmo_par.w_eos=params._w_eos();
    s_cosmo_par.N_eff=params._N_eff();
    s_cosmo_par.sigma8=params._sigma8();
    s_cosmo_par.A_s=params._A_s();
    s_cosmo_par.n_s=params._n_s();
    s_cosmo_par.alpha_s=params._alpha_s();
    s_cosmo_par.Tcmb=params._Tcmb();
    s_cosmo_par.Mabs=0.;
    s_cosmo_par.mlim=0.;
    s_cosmo_par.RR=params._RR();
    s_cosmo_par.GAL_BIAS=params._GAL_BIAS();
    s_cosmo_par.Amc=params._Amc();

    s_cosmo_par.kstar=params._kstar();
    s_cosmo_par.Delta_SO=params._Delta_SO();
    s_cosmo_par.use_wiggles=params._use_wiggles();
    s_cosmo_par.kmin_int=params._k_min_integration();
    s_cosmo_par.kmax_int=params._k_max_integration();
    s_cosmo_par.mass_function_fit=params._mass_function_fit();
    s_cosmo_par.halo_mass_bias_fit=params._halo_mass_bias_fit();
    s_cosmo_par.density_profile=params._density_profile();

    s_cosmo_par.coef_concentration_amp=params._coef_concentration_amp();
    s_cosmo_par.coef_concentration=params._coef_concentration();


    real_prec redshift=params._redshift();
    s_cosmo_par.cosmological_redshift=redshift;



    Cf.check_cosmo_pars(&s_cosmo_par);



}

// **********************************************************************************
// **********************************************************************************
// **********************************************************************************
// **********************************************************************************
// **********************************************************************************
// **********************************************************************************
// **********************************************************************************

void CosmoLib::get_cosmolib(){

  welcome_message();
  
  write_cosmo_parameters(&s_cosmo_par);

  
  real_prec redshift=params._redshift();
  ofstream hout;
  hout.open("hm_check.log");

  // *********************************************************

  real_prec critical_density=Cf.critical_density(redshift, (void *)&s_cosmo_par);
  real_prec density_contrast_top_hat=Cf.density_contrast_top_hat(redshift, (void *)&s_cosmo_par);
  real_prec Hubble_function=Cf.Hubble_function(redshift, (void *)&s_cosmo_par);
  real_prec comoving_distance=Cf.comoving_distance(redshift, (void *)&s_cosmo_par);
  real_prec comoving_angular_diameter_distance=Cf.comoving_angular_diameter_distance(redshift, (void *)&s_cosmo_par);
  real_prec mean_matter_density=Cf.mean_matter_density(redshift, (void *)&s_cosmo_par);
  real_prec age_universe=Cf.age_universe(redshift, (void *)&s_cosmo_par);
  real_prec comoving_sound_horizon=Cf.comoving_sound_horizon(redshift, (void *)&s_cosmo_par);
  real_prec growth_factor=Cf.growth_factor(redshift, (void *)&s_cosmo_par);
  real_prec growth_index=Cf.growth_index(redshift, (void *)&s_cosmo_par);
  real_prec halo_dynamical_time=Cf.halo_dynamical_time(redshift, (void *)&s_cosmo_par);
  real_prec omega_matter=Cf.omega_matter(redshift, (void *)&s_cosmo_par);
  real_prec Distance_Modulus=Cf.Distance_Modulus(redshift, (void *)&s_cosmo_par);
  real_prec pk_normalization;
  Ps.normalization((void *)&s_cosmo_par, pk_normalization);
  
  So.message_screen("Redshift =", redshift);
  So.message_screen("Omega matter at this redshift =", omega_matter);
  So.message_screen("Hubble parameter at this redshift =", Hubble_function);
  So.message_screen("Mean matter density =", mean_matter_density, "(Ms/h)/(Mpc/h)^(-3)");
  So.message_screen("Comoving distance to current redshift =", comoving_distance, "Mpc/h");
  So.message_screen("Comoving angular diameter distance to current resdhift =", comoving_angular_diameter_distance, "Mpc/h");
  So.message_screen("Comoving sound horizon =", comoving_sound_horizon, "Mpc/h");
  So.message_screen("Age of the Universe at current redshift =", age_universe/1e9, "Gyr/h");
  So.message_screen("Distance Modulus =", Distance_Modulus);
  So.message_screen("Growing mode at this redshift D(z) =", growth_factor);
  So.message_screen("Growth index at this redshift g(z) =", growth_index);
  So.message_screen("Halo-dynamical_time =", halo_dynamical_time/1e9, "Gyr/h");
  So.message_screen("Critical overdensity linearly extrapolated =", critical_density);
  So.message_screen("Top-hat density contrast at virial =",density_contrast_top_hat);
  So.message_screen("Normalization of matter power spectrum at current redshift =", pk_normalization);
  growth_factor/=Cf.growth_factor(0, (void *)&s_cosmo_par);
  // Normalize the growth factor to compute the processed linear matter power spectrum
 
  // Aca reacomodo algunos de estos factores en la estructura grande
  s_cosmo_par.critical_density=critical_density;
  s_cosmo_par.density_contrast_top_hat=density_contrast_top_hat;
  s_cosmo_par.mean_matter_density=mean_matter_density;

  s_cosmo_par.growth_factor=growth_factor;

  s_cosmo_par.pk_normalization=pk_normalization;

  // ***********************************************************************************************
  real_prec sigma_from_As=Cs.As2sigma8(&s_cosmo_par);
  
  // ***********************************************************************************************
  // Generating mass function
  // ***********************************************************************************************
  v_mass.resize(params._n_points_mf(),0);
  v_sigma_mass.resize(params._n_points_mf(),0);
  v_mass_function.resize(params._n_points_mf(),0);
  v_halo_mass_bias.resize(params._n_points_mf(),0);

  time_t start;
  time (&start);

  cout<<endl;
  So.message_screen("Computing Sigma(M)");
  
  for(int i=0; i<v_mass.size();i++)
    {
#ifdef TIME
      comp_time(start, v_mass.size(), i);
#endif
      
      if(params._scale_mf()=="linear")
	this->v_mass[i]=static_cast<gsl_real>(log10(params._M_min_mf()+i*(params._M_max_mf()-params._M_min_mf())/static_cast<double>(v_mass.size())));
      
      else if(params._scale_mf()=="log")
	this->v_mass[i]=static_cast<gsl_real>(log10(params._M_min_mf())+i*(log10(params._M_max_mf())-log10(params._M_min_mf()))/static_cast<double>(v_mass.size()));

      this->v_sigma_mass[i]=static_cast<gsl_real>(Cs.sigma_masa(v_mass[i],redshift,&s_cosmo_par));
    }



  
  So.DONE();
  
  s_cosmo_par.v_mass=this->v_mass;
  s_cosmo_par.v_sigma_mass=this->v_sigma_mass;

  So.message_screen("Computing halo mass-function");

  for(int i=0; i<v_mass.size();i++){
#ifdef TIME
    comp_time(start, v_mass.size(), i);
#endif    
    v_mass_function[i]=static_cast<gsl_real>(Cs.mass_function_D(static_cast<real_prec>(this->v_mass[i]),redshift,&s_cosmo_par));
    v_halo_mass_bias[i]=static_cast<gsl_real>(Cs.bias(static_cast<real_prec>(this->v_mass[i]),redshift,&s_cosmo_par));
  }


  Fm.write_to_file(params._mass_function_output_file(),v_mass,v_mass_function,v_sigma_mass);
  Fm.write_to_file(params._halo_mass_bias_output_file(),v_mass,v_halo_mass_bias);
  // ALlcate new comp[uted vectors in a structure to be interpolated later
  s_cosmo_par.M_max_mf=params._M_max_mf();
  s_cosmo_par.M_min_mf=params._M_min_mf();
  s_cosmo_par.v_mass=v_mass;

  s_cosmo_par.v_mass_function=v_mass_function;
  s_cosmo_par.v_halo_mass_bias=v_halo_mass_bias;


  // Once this has been created, we can compute integrals with respect to the mass.
  // Let us go for effective halo_mass bias as a function of a minimum mass
  
  v_effective_halo_mass_bias.resize(params._n_points_mf(),0);
  v_effective_halo_mean_number_density.resize(params._n_points_mf(),0);
  s_cosmo_par.M_min_effective=params._M_min_effective();
  s_cosmo_par.M_max_effective=params._M_max_effective();

  cout<<endl;
  So.message_screen("Computing effective halo mass-bias");

//#pragma omp parallel for
  for(int i=0; i<v_mass.size();i++){
#ifdef TIME
    comp_time(start, v_mass.size(), i);
#endif

    v_effective_halo_mass_bias[i]=Cs.effective_halo_mass_bias(v_mass[i],redshift,&s_cosmo_par);
    v_effective_halo_mean_number_density[i]=Cs.effective_halo_mean_number_density(v_mass[i],redshift,&s_cosmo_par);
  }  
  So.DONE();
  
  Fm.write_to_file(params._effective_halo_mass_bias_output_file(),v_mass, v_effective_halo_mass_bias);
  Fm.write_to_file(params._effective_halo_mean_number_density_output_file(),v_mass,v_effective_halo_mean_number_density);

  // ***********************************************************************************************
  // ***********************************************************************************************
  
  real_prec eff_bias=Cs.effective_halo_mass_bias( log10(params._M_min_effective()),redshift,&s_cosmo_par);
  So.message_screen("Effective halo-mass bias at the resolution (min) mhalo mass", eff_bias);
  
  // ***********************************************************************************************
  // ***********************************************************************************************
  // Calculamos las escalas que definen que es linal approx, usando sigma = 1
  // ***********************************************************************************************
  real_prec Mnl, rnl, sigman, knl;
  s_cosmo_par.aux_var1=params._redshift(); //ESTO LO PUEDO HACER MEJOR CON LAS NON SCALES DEL HALO FIT, VER CODIO cl_functions
  Cs.non_linear_scales((void *)&s_cosmo_par, &knl,&Mnl,&rnl,&sigman);
  s_cosmo_par.Mnl=Mnl;
  s_cosmo_par.knl=knl;
  s_cosmo_par.rnl=rnl;

  So.message_screen("Non linear mass scale Mnl",Mnl,"Ms/h");
  So.message_screen("Non linear wave number knl",knl,"h/Mpc");
  So.message_screen("Non linear scales rnl",rnl,"Mpc/h");
  So.message_screen("Sanity check: sigma(Mnl)",sigman);
  hout<<"Non linear scales : mass = "<<Mnl<<" Ms/h "<<endl;
  hout<<"Non linear scales : k    = "<<knl<<" h/Mpc "<<endl;
  hout<<"Non linear scales : r    = "<<rnl<<" Mpc/h "<<endl;
  hout<<"Sanity check: sigma(Mnl) = "<<sigman<<endl;

  // ***********************************************************************************************
  // ***********************************************************************************************
  // NOW WE CAN COMPUTE HALO FIT
  // ***********************************************************************************************



  v_k_ps.resize(params._n_points_ps(),0);
  v_nl_power_spectrum.resize(params._n_points_ps(),0);
  v_nl_power_spectrum_pt.resize(params._n_points_ps(),0);
  v_l_power_spectrum.resize(params._n_points_ps(),0);

  pk_aux.resize(v_k_ps.size(),0);
  kk_aux.resize(v_k_ps.size(),0);
  for(int i=0;i<v_k_ps.size();++i){
    if(params._scale_ps()=="linear")
      kk_aux[i]=(0.5*params._k_min_ps()+i*(2.*params._k_max_ps()-0.5*params._k_min_ps())/static_cast<double>(v_k_ps.size()));
    else if(params._scale_ps()=="log")
      kk_aux[i]=pow(10,(log10(0.5*params._k_min_ps())+i*(log10(2.*params._k_max_ps())-log10(0.5*params._k_min_ps()))/static_cast<double>(v_k_ps.size())));
  }
  s_cosmo_par.v_k_ps=kk_aux;
  kk_aux.clear();

  for(int i=0;i<v_k_ps.size();++i)pk_aux[i]=Ps.Linear_Matter_Power_Spectrum(&s_cosmo_par, kk_aux[i]);
  s_cosmo_par.v_lin_power_spectrum=pk_aux;
  pk_aux.clear();


  real_prec knl_hf, rnl_hf;
  vector<real_prec> rrs(100,0);
  vector<real_prec> sss(100,0);

  Ps.nl_scales_halo_fit((void *)&s_cosmo_par, &knl_hf, &rnl_hf, rrs, sss,true);

  s_cosmo_par.knl_hf=knl_hf;
  s_cosmo_par.rnl_hf=rnl_hf;
  real_prec kstar;
  Ps.kstar_integral((void *)&s_cosmo_par, &kstar);
  s_cosmo_par.kstar=kstar;


  hout<<"Non linear scales halo fit: k       = "<<knl_hf<<" h/Mpc "<<endl;
  hout<<"Non linear scales halo fit: r       = "<<rnl_hf<<" Mpc/h "<<endl;
  So.message_screen("Non linear scales halo fit k",knl_hf,"h/Mpc");
  So.message_screen("Non linear scales halo fit: r", rnl_hf," Mpc/h");
  //  hout<<"Sigma8 from As                      = "<<sigma_from_As<<endl;
  So.message_screen("k*",kstar,"h/Mpc");


  
  time (&start);
  cout<<endl;
  So.message_screen("Computing Non linear matter power spectrum");
  
  //#pragma omp parallel for
  for(int i=0; i<v_k_ps.size();i++){
#ifdef TIME
    comp_time(start, v_k_ps.size(), i);
#endif
    
    if(params._scale_ps()=="linear")v_k_ps[i]=(params._k_min_ps()+i*(params._k_max_ps()-params._k_min_ps())/v_k_ps.size());
    else if(params._scale_ps()=="log")v_k_ps[i]=pow(10,(log10(params._k_min_ps())+i*(log10(params._k_max_ps())-log10(params._k_min_ps()))/v_k_ps.size()));
    
    v_nl_power_spectrum[i]=Ps.Non_Linear_Matter_Power_Spectrum_Halo_Fit(&s_cosmo_par,v_k_ps[i]);
    v_nl_power_spectrum_pt[i]=Ps.Non_Linear_Matter_Power_Spectrum_PT(&s_cosmo_par,v_k_ps[i]);
    
    v_l_power_spectrum[i]=Ps.Linear_Matter_Power_Spectrum(&s_cosmo_par,v_k_ps[i]);
  }
  Fm.write_to_file(params._linear_matter_ps_output_file(),v_k_ps,v_l_power_spectrum);
  Fm.write_to_file(params._non_linear_matter_ps_halo_fit_output_file(),v_k_ps,v_nl_power_spectrum);
  Fm.write_to_file(params._non_linear_matter_ps_pt_output_file(),v_k_ps,v_nl_power_spectrum_pt);
  
  s_cosmo_par.v_k_ps=v_k_ps;
  s_cosmo_par.v_nl_power_spectrum=v_nl_power_spectrum;
  So.DONE();
  
  
  // ***********************************************************************************************
  // ***********************************************************************************************
  // CORRELATION FUNCTION
  // ***********************************************************************************************
  
  CorrelationFunctionTH SCf(1000, 1.);  // Check NumericalMethods to understand these inputs
  
  v_r_cf.resize(params._n_points_cf(),0);
  if(params._compute_output_linear_correlation_function())
    {
      So.message_screen("Computing matter correlation function");
      v_nl_correlation_function.resize(params._n_points_cf(),0);
      v_l_correlation_function.resize(params._n_points_cf(),0);
      time (&start);
      
      for(int i=0;i<v_r_cf.size();i++){
	//      comp_time(start, v_r_cf.size(), i);
	if(params._scale_cf()=="linear")v_r_cf[i]=(params._r_min_cf()+i*(params._r_max_cf()-params._r_min_cf())/v_r_cf.size());
	else if(params._scale_cf()=="log")v_r_cf[i]=pow(10,(log10(params._r_min_cf())+i*(log10(params._r_max_cf())-log10(params._r_min_cf()))/v_r_cf.size()));
	v_l_correlation_function[i]=SCf.Linear_Matter_Correlation_Function(&s_cosmo_par,v_r_cf[i]);
	v_nl_correlation_function[i]=SCf.Non_Linear_Matter_Correlation_Function_Halo_Fit(&s_cosmo_par,v_r_cf[i]);
      }
      Fm.write_to_file(params._linear_matter_cf_output_file(),v_r_cf,v_l_correlation_function);
      Fm.write_to_file(params._non_linear_matter_cf_halo_fit_output_file(),v_r_cf,v_nl_correlation_function);
      So.DONE();
  }
  
  // ***********************************************************************************************
  // ***********************************************************************************************
  // ***********************************************************************************************
  // *********************************************************************************************** 
 


  So.message_screen("Density profiles in configuration space");
  DensityProfiles Dp;
  v_r_dp.resize(params._n_points_dp_r(),0);
  
  v_density_profile_r.resize(params._n_points_dp_r(),0);
  for(int i=0;i<v_r_dp.size();i++){
#ifdef TIME
    comp_time(start, v_r_dp.size(), i);
#endif
    
    if(params._scale_dp_r()=="linear")v_r_dp[i]=(params._r_min_dp()+i*(params._r_max_dp()-params._r_min_dp())/v_r_dp.size());
    else if(params._scale_dp_r()=="log")v_r_dp[i]=pow(10,(log10(params._r_min_dp())+i*(log10(params._r_max_dp())-log10(params._r_min_dp()))/v_r_dp.size()));
    v_density_profile_r[i]=Dp.density_r(v_r_dp[i], log10(Mnl/M_reference), redshift, (void *)&s_cosmo_par);
  }

  So.DONE();
  Fm.write_to_file(params._density_profile_r_output_file(),v_r_dp,v_density_profile_r);

  
  So.message_screen("Density profiles in Fourier space");
  v_k_dp.resize(params._n_points_dp_k(),0);
  v_density_profile_k.resize(params._n_points_dp_k(),0);
  for(int i=0;i<v_k_dp.size();i++){
#ifdef TIME
    comp_time(start, v_k_dp.size(), i);
#endif

    if(params._scale_dp_k()=="linear")v_k_dp[i]=(params._k_min_dp()+i*(params._k_max_dp()-params._k_min_dp())/v_r_dp.size());
    else if(params._scale_dp_k()=="log")v_k_dp[i]=pow(10,(log10(params._k_min_dp())+i*(log10(params._k_max_dp())-log10(params._k_min_dp()))/v_k_dp.size()));
    v_density_profile_k[i]=Dp.density_k(v_k_dp[i], log10(Mnl/M_reference), redshift, &s_cosmo_par);
  }
  Fm.write_to_file(params._density_profile_k_output_file(),v_k_dp,v_density_profile_k);

  So.DONE();
  
  s_cosmo_par.v_density_profile_k=v_density_profile_k;
  s_cosmo_par.v_density_profile_r=v_density_profile_r;
  
  
  // ***********************************************************************************************
  // ***********************************************************************************************
  // POWER SPECTRUM IN THE HALO MODEL.
  // ***********************************************************************************************
  s_cosmo_par.hod_model=params._hod_model();
  s_cosmo_par.mmin_hod=params._mmin_hod();
  s_cosmo_par.alpha_hod=params._alpha_hod();
  s_cosmo_par.scatter_hod=params._scatter_hod();
  s_cosmo_par.muno_hod=params._muno_hod();

  
  // ***********************************************************************************************
  // 1 HALO TERM: SATELLITE-SATELLIT AND CENTRAL-SATELLITE.
  // 2 HALO TERM: SATELLITE1-SATELLIT2, CENTRAL 1-SATELLITE 2, CENTRAL-CENTRAL:
  // condensed in a single scale dependnet bias and the linear power spectrum
    // ***********************************************************************************************
  
  v_galaxy_power_1h_ss.resize(params._n_points_ps(),0);
  v_galaxy_power_1h_sc.resize(params._n_points_ps(),0);
  v_galaxy_power_2h.resize(params._n_points_ps(),0);
  v_galaxy_matter_bias.resize(params._n_points_ps(),0);
  v_galaxy_power_spectrum.resize(params._n_points_ps(),0);
  
  real_prec mean_gal_density=Cs.mean_galaxy_number_density(redshift, &s_cosmo_par);
  So.message_screen("Mean galaxy number density",mean_gal_density,"(h/Mpc)^3 ");
  real_prec mean_halo_density=Cs.mean_halo_number_density(&s_cosmo_par);
  So.message_screen("Mean halo number density",mean_halo_density,"(h/Mpc)^3");
  
  Ps.v_mass=v_mass;
  Ps.v_mass_function=v_mass_function;
  Ps.v_halo_mass_bias=v_halo_mass_bias;
  
  
  // ***********************************************************************************************
  So.message_screen("Computing galaxy power spectrum");
  for(int i=0;i<v_k_ps.size();i++)
    {
#ifdef TIME
      comp_time(start, v_k_ps.size(), i);
#endif
      v_galaxy_power_1h_ss[i]=Ps.Galaxy_power_spectrum_h1_ss(&s_cosmo_par, v_k_ps[i], redshift)/pow(mean_gal_density,2);
      v_galaxy_power_1h_sc[i]=Ps.Galaxy_power_spectrum_h1_sc(&s_cosmo_par, v_k_ps[i], redshift)/pow(mean_gal_density,2);
      v_galaxy_matter_bias[i]=Ps.Galaxy_matter_bias(&s_cosmo_par, v_k_ps[i], redshift)/mean_gal_density;
      v_galaxy_power_2h[i]=v_l_power_spectrum[i]*(v_galaxy_matter_bias[i],2);
      v_galaxy_power_spectrum[i]=v_galaxy_power_1h_ss[i]+v_galaxy_power_1h_sc[i]+v_galaxy_power_2h[i];
    }
  Fm.write_to_file(params._galaxy_power_spectrum_halo_model_output_file(),v_k_ps,v_galaxy_power_1h_ss,v_galaxy_power_1h_sc,v_galaxy_matter_bias,v_galaxy_power_2h,v_galaxy_power_spectrum);
  
  
  s_cosmo_par.v_galaxy_power_spectrum_1h_ss=this->v_galaxy_power_1h_ss;
  s_cosmo_par.v_galaxy_power_spectrum_1h_sc=this->v_galaxy_power_1h_sc;
  s_cosmo_par.v_galaxy_power_spectrum_2h=this->v_galaxy_power_2h;


  // ***********************************************************************************************
  // ***********************************************************************************************
  // GALAXY_CORRELATION FUNCTION HALO MODEL
  // Para calcular esto es mandatorio calcular el espectro

  // ***********************************************************************************************
  if(params._compute_output_non_linear_correlation_function())
    {
      v_galaxy_correlation_1h_ss.resize(params._n_points_cf(),0);
      v_galaxy_correlation_1h_sc.resize(params._n_points_cf(),0);
      v_galaxy_correlation_2h.resize(params._n_points_cf(),0);
      v_galaxy_correlation.resize(params._n_points_cf(),0);

      So.message_screen("Computing galaxy correlation function");
      for(int i=0;i<v_r_cf.size();i++){
#ifdef TIME
	comp_time(start, v_r_cf.size(), i);
#endif
	v_galaxy_correlation_1h_ss[i]=SCf.Galaxy_Correlation_Function_1h_ss(&s_cosmo_par, v_r_cf[i]);
	v_galaxy_correlation_1h_sc[i]=SCf.Galaxy_Correlation_Function_1h_sc(&s_cosmo_par, v_r_cf[i]);
	v_galaxy_correlation_2h[i]=SCf.Galaxy_Correlation_Function_2h(&s_cosmo_par, v_r_cf[i]);
	v_galaxy_correlation[i]=v_galaxy_correlation_2h[i]+v_galaxy_correlation_1h_ss[i]+v_galaxy_correlation_1h_sc[i];
      }
      So.DONE();
      Fm.write_to_file(params._galaxy_correlation_function_halo_model_output_file(),
		       v_r_cf,
		       v_galaxy_correlation_1h_ss,
		       v_galaxy_correlation_1h_sc,
		       v_galaxy_correlation_2h,
		       v_galaxy_correlation
		       );
      
    }
  
  
  // ***********************************************************************************************
  // ***********************************************************************************************
  ending_message();
  
}


