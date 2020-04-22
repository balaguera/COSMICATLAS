
/** @file Screen Output
 *
 *  @brief 
 */



# include "../Headers/ScreenOutput.h"

// **************************************************************************
// **************************************************************************

void ScreenOutput::message(string mess){cout<<mess<<endl;}

void ScreenOutput::welcome_message(){
 time_t rawtime;
 time (&rawtime);
 cout<<YELLOW<<endl;
 cout<<"************************************************************************************"<<endl;
 cout<<"************************************************************************************"<<endl;
 cout<<"POWER SPECTRUM OF COSMOLOGICAL MASS TRACERS USING FKP (+FFTW) ESTIMATOR            *"<<endl;
 cout<<"VERSION 1.1                                                                        *"<<endl;
 cout<<"For more documentation see README/README.pdf                                       *"<<endl;
 cout<<"************************************************************************************"<<endl;
 cout<<"Starting time and date "<<ctime (&rawtime);
 cout<<"************************************************************************************"<<endl;
 cout<<RESET<<endl;

}


// **************************************************************************
// **************************************************************************
void ScreenOutput::welcome_message_cl(){
 time_t rawtime;
 time (&rawtime);
 cout<<YELLOW<<endl;
 cout<<"************************************************************************************"<<endl;
 cout<<"************************************************************************************"<<endl;
 cout<<"ANGULAR POWER SPECTRUM OF COSMOLOGICAL MASS TRACERS                                *"<<endl; 
 cout<<"USING HEALPIX DECOMPOSITION AND PEEBLES ESTIMATOR                                  * "<<endl;
 cout<<"VERSION 1.0                                                                        *"<<endl ;
 cout<<"For more documentation see README/README.pdf                                       *"<<endl;
 cout<<"************************************************************************************"<<endl;
 cout<<"Starting time and date"<<ctime (&rawtime);
 cout<<"************************************************************************************"<<endl;
 cout<<RESET<<endl;

}

// **************************************************************************
// **************************************************************************
void ScreenOutput::welcome_message_fb(){
 time_t rawtime;
 time (&rawtime);
 cout<<YELLOW<<endl;
 cout<<"************************************************************************************"<<endl;
 cout<<"************************************************************************************"<<endl;
 cout<<"3D POWER SPECTRUM OF COSMOLOGICAL MASS TRACERS                                     *"<<endl; 
 cout<<"USING THE FOURIER BESSEL DECOMPOSITION                                             * "<<endl;
 cout<<"VERSION 1.0                                                                        *"<<endl ;
 cout<<"For more documentation see README/README.pdf                                       *"<<endl;
 cout<<"************************************************************************************"<<endl;
 cout<<"Starting time "<<ctime (&rawtime);
 cout<<"************************************************************************************"<<endl;
 cout<<RESET<<endl;

}
// **************************************************************************
// **************************************************************************


void ScreenOutput::welcome_message_yama(){
  time_t rawtime;
  time (&rawtime);
  cout<<YELLOW<<endl;
  cout<<"***********************************************************************************"<<endl;
  cout<<"***********************************************************************************"<<endl;
  cout<<"MULTIPOLE DECOMPOSITION OF THE THREE DIMENSIONAL POWER SPECTRUM OF COSMOLOGICAL   *"<<endl;
  cout<<"MASS TRACERS USING YAMAMOTO ESTIMATOR (FFTW-based)                                *"<<endl;
  cout<<"VERSION 1.1                                                                       *"<<endl;
  cout<<"For more documentation see README/README.pdf                                      *"<<endl;
  cout<<"***********************************************************************************"<<endl;
  cout<<"Starting time "<<ctime (&rawtime);
  cout<<"***********************************************************************************"<<endl;
  cout<<RESET<<endl;

}

// **************************************************************************
// **************************************************************************

void ScreenOutput::welcome_message_bispectrum(){
 cout<<YELLOW<<endl;
 time_t rawtime;
 time (&rawtime);
 cout<<"************************************************************************************"<<endl;
 cout<<"************************************************************************************"<<endl;
 cout<<"BISPECTRUM OF COSMOLOGICAL MASS TRACERS USING FKP ESTIMATOR                        *"<<endl; 
 cout<<"VERSION 1.1                                                                        *"<<endl;
 cout<<"Documentation and warnings in readme.ps                                            *"<<endl;
 cout<<"************************************************************************************"<<endl;
 cout<<"Starting time "<<ctime (&rawtime);
 cout<<"************************************************************************************"<<endl;
 cout<<RESET;

}


void ScreenOutput::welcome_message_bispectrum_fast()
{
  time_t rawtime;
  time (&rawtime);
  cout<<YELLOW<<endl;
  cout<<"************************************************************************************"<<endl;
  cout<<"************************************************************************************"<<endl;
  cout<<"BISPECTRUM OF COSMOLOGICAL MASS TRACERS USING FKP ESTIMATOR                        *"<<endl; 
  cout<<"AND INVERSE FFTW TRICK                                                             *"<<endl; 
  cout<<"VERSION 1.0                                                                        *"<<endl;
  cout<<"Documentation and warnings in readme.ps                                            *"<<endl;
  cout<<"************************************************************************************"<<endl;
  cout<<"Starting time "<<ctime (&rawtime);
  cout<<"************************************************************************************"<<endl;
  cout<<RESET;
}

// **************************************************************************
// **************************************************************************

void ScreenOutput::message_interrupt()
{
  cout<<YELLOW<<endl;
  cout<<"***********************************************************************************"<<endl;
  cout<<"Code interrupted at"<<endl; 
  time_t rawtime;
  time ( &rawtime );
  cout<<ctime (&rawtime);
  cout<<"***********************************************************************************"<<endl;
  cout<<"   "<<endl;
}

// **************************************************************************
// **************************************************************************

void ScreenOutput::done(string text)
{
  cout<<"...leaving "<<text<<endl;
  cout<<"\t"<<endl;
}
// **************************************************************************
// **************************************************************************
void ScreenOutput::enter(string text)
{
  cout<<"\t"<<endl;
  cout<<"Going to:  "<<text<<endl;
}
// **************************************************************************
// **************************************************************************
void ScreenOutput::ending_message()
{

  cout<<YELLOW<<endl;
  time_t rawtime;
  time ( &rawtime );
  cout<<"Ending date:"<<ctime (&rawtime)<<endl;
  cout<<RESET<<endl;
}
// **************************************************************************
// **************************************************************************

void ScreenOutput::write_fftw_parameters(void *p)
{
  struct s_parameters_box * s_cp= (struct s_parameters_box *)p;
  ScreenOutput So;
  cout<<BLUE;
  cout<<endl;
  std::cout<<CYAN<<"Selected options :"<<RESET<<std::endl;
  cout<<YELLOW<<"MAS:"<<CYAN<<(s_cp->mas)<<endl;
  if(s_cp->use_MAS_correction) cout<<YELLOW<<"MAS correction"<<CYAN<<" enabled"<<endl;
  else cout<<YELLOW<<"MAS correction"<<CYAN<<" disabled"<<endl;
  if(s_cp->FKP_weight)cout<<YELLOW<<"Using FKP weights with Pest = "<<CYAN<<s_cp->Pest<<endl;
  else cout<<YELLOW<<"Using weights = "<<CYAN<<" 1"<<endl;
  if(s_cp->FKP_error_bars)cout<<YELLOW<<"Computing FKP error bars"<<endl;
  else cout<<YELLOW<<"Estimate without error bars."<<endl;
  if(s_cp->use_SN_correction)cout<<YELLOW<<"Shot-noise correction:"<<CYAN<<" enabled"<<endl;
  else cout<<YELLOW<<"Shot-noise correction"<<CYAN<<" disabled  "<<RESET<<endl;
  cout<<RESET;
}

// **************************************************************************
// **************************************************************************

void ScreenOutput::write_parameters_estimator(void *p)
{

  struct s_parameters_estimator * s_cp= (struct s_parameters_estimator *)p;
  cout<<YELLOW<<"PARAMETERS of the FKP estimator: "<<RESET<<endl;
  cout<<YELLOW<<"Number of objects = "<<CYAN<<s_cp->number_of_objects<<endl;
  if(s_cp->use_random_catalog)cout<<CYAN<<"Number of random objects= "<<s_cp->number_of_randoms<<endl;
  cout<<YELLOW<<"Weighted number of objects = "<<CYAN<<s_cp->w_number_of_objects<<endl;
  if(s_cp->use_random_catalog)cout<<"Weighted number of random objects = "<<CYAN<<s_cp->w_number_of_randoms<<endl;
  cout<<YELLOW<<"alpha = "<<CYAN<<s_cp->alpha<<endl;
  cout<<YELLOW<<"Normalization = "<<CYAN<<s_cp->normalization<<endl;
  cout<<YELLOW<<"Shot_noise (power) = "<<CYAN<<s_cp->shot_noise<<"  (Mpc h^-1)^3"<<endl;
  cout<<YELLOW<<"Shot_noise (window) = "<<CYAN<<s_cp->shot_noise_window<<endl;
  if(s_cp->use_random_catalog)cout<<YELLOW<<"Mean number density (from weights) = "<<CYAN<<(s_cp->number_of_randoms-s_cp->w_number_of_randoms)/(20000.0*s_cp->w_number_of_randoms)<<" (Mpc h^-1)^(-3)"<<endl; 
  cout<<RESET;
}

// **************************************************************************
// **************************************************************************


void ScreenOutput::write_parameters_b_estimator(void *p)
{

  struct s_parameters_bis_estimator * s_cp= (struct s_parameters_bis_estimator *)p;
  cout<<YELLOW;
  cout<<"PARAMETERS bispectrum estimator                                                   *"<<endl;
  cout<<YELLOW<<"Number of objects = "<<CYAN<<s_cp->number_of_objects<<endl;
  if(s_cp->use_random_catalog)cout<<"Number of random  = "<<CYAN<<s_cp->number_of_randoms<<endl;
  cout<<"Weighted number of objects = "<<CYAN<<s_cp->w_number_of_objects<<endl;
  if(s_cp->use_random_catalog)cout<<"Weighted number of random  = "<<CYAN<<s_cp->w_number_of_randoms<<endl;
  cout<<YELLOW<<"alpha = "<<s_cp->alpha<<endl;
  cout<<YELLOW<<"Normalization bispectrum= "<<CYAN<<s_cp->normalization<<"  (Mpc h^-1)^6"<<endl;
  cout<<YELLOW<<"Normalization power spectrum = "<<CYAN<<s_cp->normal_power<<"  (Mpc h^-1)^3"<<endl;
  cout<<YELLOW<<"Shot_noise S1 B= "<<CYAN<<s_cp->shot_noise1<<"  (Mpc h^-1)^-3"<<endl;
  cout<<YELLOW<<"Shot_noise S2 B= "<<CYAN<<s_cp->shot_noise2<<"  (Mpc h^-1)^-3"<<endl;
  cout<<YELLOW<<"Shot_noise  P(k) = "<<CYAN<<s_cp->shot_noise_power<<"  (Mpc h^-1)^-3"<<endl;
  cout<<"***********************************************************************************"<<endl;
  cout<<endl;
  cout<<RESET;
}

// **************************************************************************
// **************************************************************************

void ScreenOutput::write_cosmo_parameters(void *p, void *pa)
{
  cout<<endl;  
  this->message_screen("**Cosmological parameters");
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  struct s_CosmoInfo * s_ci= (struct s_CosmoInfo *)pa;
  this->message_screen("Input redshift z =", s_cp->cosmological_redshift);
  this->message_screen("Omega matter =", s_cp->Om_matter);
  this->message_screen("Omega vac =", s_cp->Om_vac);
  this->message_screen("Omega CDM", s_cp->Om_cdm);
  this->message_screen("Omega baryons =", s_cp->Om_baryons);
  this->message_screen("Omega curv =", s_cp->Om_k);
  this->message_screen("Hubble par =", s_cp->Hubble);
  this->message_screen("hubble par =", s_cp->hubble);
  this->message_screen("Spectral index =", s_cp->n_s);
  this->message_screen("Sigma8 =", s_cp->sigma8);
  this->message_screen("BAO wiggles = ",s_cp->use_wiggles);
  cout<<endl;
  this->message_screen("**Cosmological information");  
  this->message_screen("Scale factor a =", 1./(s_cp->cosmological_redshift+1.));
  this->message_screen("Growth factor D1(z) (D1(z=0)=1) =",s_ci->growth_factor);
  this->message_screen("Growth factor D2(z) =",s_ci->D2);
  this->message_screen("Growth index f(z) =",s_ci->growth_index);
  this->message_screen("H(z) =",s_ci->Hubble_parameter, "(km/s)/(Mpc/h)");
  this->message_screen("Comoving distance r(z) =",s_ci->comoving_distance, "Mpc/h");
  this->message_screen("Comoving AD distance d(z) =", s_ci->comoving_angular_diameter_distance, "Mpc/h");
  this->message_screen("Comoving sound horizon rs(z) =",s_ci->comoving_sound_horizon, "Mpc/h");
  this->message_screen("Mean mass density =",s_ci->mean_matter_density, "Ms/(Mpc/h)^3");
  this->message_screen("Age of the Universe =",s_ci->age_universe, "Years/h");
  this->message_screen("Omega matter (z) =",s_ci->omega_matter);
#ifdef _USE_PATCHY_
  this->message_screen("Normalization of Linear Matter P(k,z=0) =",s_cp->pk_normalization);  
#endif
}

void ScreenOutput::write_cosmo_parameters(void *p)
{
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  cout<<BOLDCYAN<<"Input Cosmological parameters"<<RESET<<endl;
  this->message_screen("Omega matter =", s_cp->Om_matter);
  this->message_screen("Omega vac =", s_cp->Om_vac);
  this->message_screen("Omega baryons =", s_cp->Om_baryons);
  this->message_screen("Omega curv =", s_cp->Om_k);
  this->message_screen("Hubble par =", s_cp->Hubble);
  this->message_screen("Spectral index =", s_cp->n_s);
  this->message_screen("Sigma8 =", s_cp->sigma8);
}

// **************************************************************************
// **************************************************************************

void ScreenOutput::error_ncolumns(string fname){
  cout<<RED;
  cout<<fname<<" with less than 4 columns"<<endl;
  char yn;
  ScreenOutput So;
  cout<<"Continue (y/n)?"<<endl; cin>>yn;
  if(yn !='y' && yn !='n')
    {
      cout<<"Please answer  yes (y) or not (n)"<<endl;
      cout<<"Continue ?"<<endl; cin>>yn;
    }
  else{
    if(yn !='y'){
    }
    So.message_interrupt();
  }
  cout<<RESET;
}

// **************************************************************************
// **************************************************************************

void ScreenOutput::comp_time(time_t start, unsigned long full, unsigned long step){
  // THIS IS TIME CONSUMING
  double fraction=100.0*((double)(step))/((double)full);
  time_t end;
  time(&end);
  double lapse=difftime(end,start);
  if(lapse<=60)cout<<"\r "<<fraction<<" % completed. Time elapsed: "<<lapse<<"  secs \r";cout.flush();
  if(lapse>60)cout<<"\r "<<fraction<<" % completed. Time elapsed: "<<lapse/60.<<"  mins \r";cout.flush();
  if(lapse>3600)cout<<"\r "<<fraction<<" % completed. Time elapsed: "<<lapse/3600.<<"  hrs \r";cout.flush();
  if (fraction==25.0) cout <<"..25%  completed"<<endl;
  if (fraction==50.0) cout <<"..50%  completed"<<endl;
  if (fraction==75.0) cout <<"..75%  completed"<<endl;
  if (fraction==100.0)cout <<"..100% completed"<<endl;  
}
// **************************************************************************
// **************************************************************************

//##################################################################################
//##################################################################################
//##################################################################################
void ScreenOutput::message_warning(string ss)
{
  std::cout<<RED<<"Warning"<<RESET<<std::endl;
  std::cout<<GREEN<<ss<<RESET<<std::endl;
}

void ScreenOutput::message_warning(string ss, ULONG line)
{
  std::cout<<RED<<"Warning"<<RESET<<std::endl;
  std::cout<<GREEN<<ss<<"  "<<line<<RESET<<std::endl;
}

//##################################################################################

int ScreenOutput::message_error(string ss)
{
  std::cerr<<CYAN<<ss<<RESET<<std::endl;
  std::cerr<<RED<<"...Exiting..."<<endl;
  exit (0);
}
//##################################################################################
int ScreenOutput:: message_error(string ss, double a)
{
  std::cerr<<CYAN<<ss<<" "<<a<<RESET<<std::endl;
  std::cerr<<RED<<"...Exiting..."<<endl;
  exit (0);
}
//##################################################################################

void ScreenOutput::message_error(string ss, double a, string sa)
{
  std::cerr<<CYAN<<ss<<""<<a<<""<<sa<<RESET<<std::endl;
  std::cerr<<RED<<"...Exiting..."<<endl;
  exit (0);
}


//##################################################################################
void ScreenOutput::message_screen(string ss)
{
  std::clog<<YELLOW<<ss<<RESET<<std::endl;
}

void ScreenOutput::message_screen_flush(string ss, real_prec s2, string sa, real_prec s3)
{
  std::cout<<"\r"<<YELLOW<<ss<<RESET<<" "<<s2<<" "<<sa<<"  "<<s3<<RESET; cout.flush();
}

void ScreenOutput::message_screen_flush(string ss, int s2)
{
  std::cout<<"\r"<<YELLOW<<ss<<RESET<<" "<<s2<<RESET;cout.flush();
}

//##################################################################################
void ScreenOutput::message_screen(string ss, double s2)
{
  //  cout.precision(12);
  std::cout<<YELLOW<<ss<<RESET<<" "<<CYAN<<s2<<RESET<<std::endl;
}

void ScreenOutput::message_screen(string ss, int i, string sa,double s2)
{
  std::cout<<YELLOW<<ss<<RESET<<" "<<i<<" "<<sa<<"  "<<s2<<RESET<<std::endl;
}

//##################################################################################
void ScreenOutput::message_screen(string ss, int s2)
{
  std::clog<<YELLOW<<ss<<RESET<<" "<<CYAN<<s2<<RESET<<std::endl;
}

void ScreenOutput::message_screen(string ss, ULONG s2)
{
  std::clog<<YELLOW<<ss<<RESET<<" "<<CYAN<<s2<<RESET<<std::endl;
}

//##################################################################################
void ScreenOutput::message_screen(string ss, string s2){
  std::clog<<YELLOW<<ss<<RESET<<" "<<CYAN<<s2<<RESET<<std::endl;
}



//##################################################################################
void ScreenOutput::message_screen(string ss, double a, string s2, string f){
  std::clog<<YELLOW<<ss<<RESET<<" "<<CYAN<<a<<" "<<s2<<"  "<<f<<RESET<<std::endl;
}


//##################################################################################
void ScreenOutput::message_screen(string ss, double d2, string s2)
{
  std::clog<<YELLOW<<ss<<RESET<<" "<<CYAN<<d2<<" "<<s2<<RESET<<std::endl;
}


//##################################################################################
void ScreenOutput::message_screen(string ss, double d2, time_t time)
{
  double lapse=difftime(time,this->initial_time);
  if(lapse<60)
    std::clog<<RED<<ss<<" "<<d2<<" (T_init + "<<lapse<<" secs)"<<RESET<<std::endl;
  else if(lapse>60)
    std::clog<<RED<<ss<<" "<<d2<<" (T_init + "<<lapse/60<<" mins)"<<RESET<<std::endl;
  else if(lapse>3600)
    std::clog<<RED<<ss<<" "<<d2<<" (T_init + "<<lapse/3600<<" hrs)"<<RESET<<std::endl;
}
//##################################################################################

void ScreenOutput::message_screen(string ss, double d2, time_t time, time_t time2)
{
  double lapse=difftime(time,this->initial_time);
  double lapse2=difftime(time,time2);

  if(lapse<60)
    std::clog<<BOLDCYAN<<ss<<" "<<d2<<" (T_init + "<<lapse<<" secs, time spent in last step "<<lapse2<<" s)"<<RESET<<std::endl;
  else if(lapse>60)
    std::clog<<BOLDCYAN<<ss<<" "<<d2<<" (T_init + "<<lapse/60<<" mins, time spent in last step "<<lapse2<<" s)"<<RESET<<std::endl;
  else if(lapse>3600)
    std::clog<<BOLDCYAN<<ss<<" "<<d2<<" (T_init + "<<lapse/3600<<" hrs, time spent in last step "<<lapse2<<" s)"<<RESET<<std::endl;
}


//##################################################################################

void ScreenOutput::message_screen(string ss, string s2, string s3, double d3)
{
  std::clog<<GREEN<<ss<<CYAN<<" "<<s2<<" "<<BLUE<<s3<<"  "<<CYAN<<d3<<RESET<<std::endl;
}

//##################################################################################
void ScreenOutput::message_time(time_t start_all)
{
  time_t end;
  time(&end);
  cout<<CYAN<<"**************************************************************************"<<RESET<<endl;
  cout<<CYAN<<"**************************************************************************"<<RESET<<endl;
  cout<<CYAN<<"**************************************************************************"<<RESET<<endl;
  double lapse=difftime(end,start_all);
  cout<<"BAM ends on "<<__DATE__<<" at "<<__TIME__<<"                             *"<<endl;
  cout<<"Time elapsed: "<<lapse<<"   secs "<<RESET<<endl;
  cout<<CYAN<<"**************************************************************************"<<RESET<<endl;
  cout<<CYAN<<"**************************************************************************"<<RESET<<endl;
  cout<<CYAN<<"**************************************************************************"<<RESET<<endl;
}

//##################################################################################
void ScreenOutput::message_time_mock(time_t start_all)
{
  time_t end;
  time(&end);
  double lapse=difftime(end,start_all);
  cout<<YELLOW<<"Mock DF generated in : "<<lapse<<"   secs "<<RESET<<endl;
}



void ScreenOutput::message_time2(time_t start_all)
{
  time_t end;
  time(&end);
  double lapse=difftime(end,start_all);
  cout<<YELLOW<<"Time of calculation = "<<CYAN<<lapse<<" secs "<<RESET<<endl;
}

//##################################################################################
void ScreenOutput::usage(string s)
{
  cout<<BOLDGREEN;  
  cout<<"*****************************************************************"<<endl;
  cout<<"*****************************************************************"<<endl;
  cout<<"*****************************************************************"<<endl;
  cout<<"*****************************************************************"<<endl;
  cout<<"COSMICATLASS                                                    *"<<endl;
  cout<<"CosmologicalCATalogs for LArge Scale Structure                  *"<<endl;
  cout<<"Usage:"<<s<<" [-option] [argument]                              *"<<endl;
  cout<<"Options: -h current help, no arguement                          *"<<endl;
  cout<<"         -a for information on the author, no argument          *"<<endl;
  cout<<"         -c parameter_file.ini: to execute BAM                  *"<<endl;
  cout<<"         -m parameter_file.ini: to mesure power                 *"<<endl;
  cout<<"         -i parameter_file.ini: to show input pars              *"<<endl;
  cout<<"         -p parameter_file.ini: to execute Patchy               *"<<endl;
  cout<<"Consult ../Headers/def.h for pre-procesor directives:           *"<<endl;
  cout<<YELLOW;
  cout<<"*****************************************************************"<<endl;
  cout<<"*****************************************************************"<<endl;
  cout<<"*****************************************************************"<<endl;
  cout<<"*****************************************************************"<<RESET<<endl;
  cout<<RESET<<endl;                                       
}
//##################################################################################

void ScreenOutput::author(){
  cout<<BOLDGREEN;  
  cout<<"************************************************************************"<<endl;
  cout<<"************************************************************************"<<endl;
  cout<<BOLDCYAN;
  cout<<" ****  **  |*** *      * **  ****    *   ***** *       *    **** "<<endl;
  cout<<"|     *  *  *   * *  * * ** *       * *    *   *      * *     *  "<<endl;
  cout<<"|     *  *   *| *  **  * ** *      *****   *   *     *****     * "<<endl;
  cout<<" ****  **  **** *      * **  **** *     *  *   **** *     * **** "<<RESET<<endl;
  cout<<BOLDGREEN;  
  cout<<"***********************************************************************"<<endl;
  cout<<"***********************************************************************"<<endl;
  cout<<"COSMICATLASS                                                          *"<<endl;
  cout<<"CosmologicalCATalogs for LArge Scale Structure                        *"<<endl;
  cout<<"based on the the Bias Assignment Method and PATCHY code               *"<<endl;
  cout<<"Andres Balaguera Antolinez (balaguera@iac.es)                         *"<<endl;
  cout<<"Francisco-Shu Kitaura (fkitaura@iac.es)                               *"<<endl;
  cout<<"Instituto de Astrofisica de Canarias                                  *"<<endl;
  cout<<"Code developed for the paper https://arxiv.org/abs/1806.05870         *"<<endl;
  cout<<"FoF adapted from HADRON (Cheng Zhao)                                  *"<<endl;
  cout<<"Cosmic web classification adapted from CLASSLIN (F. Kitaura)          *"<<endl;
  cout<<"Power Spectrum from COSMOLib (A. Balaguera)                           *"<<endl;
  cout<<"***********************************************************************"<<endl;
  cout<<"***********************************************************************"<<endl;
  cout<<"***********************************************************************"<<endl;
  cout<<"***********************************************************************"<<endl;
  cout<<RESET<<endl;
}
//##################################################################################

void  ScreenOutput::message(time_t start_all)
{
  cout<<CYAN; 
  cout<<"*****************************************************************"<<endl;
  cout<<"*****************************************************************"<<endl;
  cout<<"*****************************************************************"<<endl;
  cout<<"*****************************************************************"<<endl;
  cout<<"*                        COSMICATLASS                           *"<<endl;
  cout<<"*       Cosmological CATalogs for LArge Scale Structure         *"<<endl;
  cout<<"*                       IAC 2017-2019                           *"<<endl;
  cout<<"*****************************************************************"<<endl;
  cout<<"*****************************************************************"<<endl;
  cout<<"*****************************************************************"<<endl;
  cout<<"*****************************************************************"<<endl;
  cout<<"                    "<<endl;
  cout<<"*****************************************************************"<<endl;
  cout<<"Compilation began on "<<__DATE__<<" at "<<__TIME__<<"                    *"<<endl;
#ifdef _USE_OMP_
  cout<<"Using OMP with "<<omp_get_max_threads()<<" threads                                       *"<<endl;
#endif
  cout<<"*****************************************************************"<<endl;
  cout<<RESET<<endl; 
}


void ScreenOutput::message_output_file(string s, ULONG l)
{
  cout<<YELLOW<<"Writting output file "<<CYAN<<s<<" with "<<l<<" lines"<<RESET<<endl;
}


void ScreenOutput::message_output_file(string s, int nlines, int ncols)
{
  cout<<YELLOW<<"Writting output file "<<CYAN<<s<<" with "<<nlines<<" lines and "<<ncols<<" columns"<<RESET<<endl;
}


void ScreenOutput::DONE()
{
  cout<<BOLDGREEN<<"                                               ["<<BOLDBLUE<<"DONE"<<BOLDGREEN<<"]"<<RESET<<endl;
  cout<<endl;
}

  
