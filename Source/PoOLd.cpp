/** @file PowerSpectrumF.cpp
 *
 * @brief This file contains a public PowerSpectrum-class member function
 * @details Generates thee estimates of power spectrum
 * @author Andres Balaguera Antolinez
 */



#undef MINERVA
//#define MINERVA

#include "../Headers/PowerSpectrumF.h"

// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
/*
#ifdef _USE_ALL_PK_
void PowerSpectrumF::add_catalogues()
 {
    this->tracer_cat.set_params_catalog(this->params);
    this->tracer_cat.type_of_object="TRACfER";
    this->tracer_cat.read_catalog(this->file_data,0);
    this->gc_n_columns=this->tracer_cat.NCOLS;
    this->N_galaxy=this->tracer_cat.NOBJS;

    if(true==use_random_catalog)
      {
        this->tracer_cat.set_params_catalog(this->params);
        this->tracer_cat.type_of_object="RANDOM";
        this->random_cat.read_catalog(this->file_random,0);
        this->rc_n_columns=this->random_cat.NCOLS;
        this->N_random=this->random_cat.NOBJS;
      }
}
#endif
*/

// *************************************************************************************
// *************************************************************************************


#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
void PowerSpectrumF::add_catalogues(real_prec mcut)
 {
#elif defined (_USE_MASS_BINS_PK_)
void PowerSpectrumF::add_catalogues(real_prec m_min, real_prec m_max)
 {
#endif

    this->So.enter(__PRETTY_FUNCTION__);

    this->tracer_cat.set_params_catalog(this->params);
   this->tracer_cat.type_of_object="TRACER";

#if defined(_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
   this->tracer_cat.read_catalog_power(this->file_data, mcut);
#elif defined (_USE_MASS_BINS_PK_)
       this->tracer_cat.read_catalog_power(file_data,m_min,m_max);
#endif
    this->gc_n_columns=this->tracer_cat.NCOLS;
    this->N_galaxy=this->tracer_cat.NOBJS;

   if(true==use_random_catalog)
     {
       this->random_cat.set_params_catalog(this->params);
       this->random_cat.type_of_object="RANDOM";
#if defined(_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
      this->random_cat.read_catalog_power(this->file_random, mcut);
#elif defined (_USE_MASS_BINS_PK_)
       this->random_cat.read_catalog_power(file_random,m_min,m_max);
#endif
       this->rc_n_columns=this->random_cat.NCOLS;
       this->N_random=this->random_cat.NOBJS;
    }


}


// *************************************************************************************
// *************************************************************************************

void PowerSpectrumF::write_power_spectrum()
{
  // Write P(k) to file:

    this->So.enter(__PRETTY_FUNCTION__);


  if(statistics=="Pk_fkp")
    {
#ifdef _WRITE_MULTIPOLES_
      c_Fm.write_to_file(this->file_power,this->kvector_data,this->pk0,this->pk2,this->pk4,this->sigma_fkp,this->modes_g);
#else
      cout<<"writting power"<<endl;
      c_Fm.write_to_file(this->file_power,this->kvector_data,this->pk0,this->modes_g);
#endif

      // Write P(kperp, kpar) to file:
#ifdef _WRITE_2DPOWER_
      c_Fm.write_to_file(this->file_power2d,this->kvector_data2d,this->kvector_data2d,this->pkk);
      //Write P(k, mu) to file:
      c_Fm.write_to_file(file_power2d_mk,this->muvector,this->kvector_data2d,this->pmk);
#endif
      // Write W(k):
      if(true==use_random_catalog)
	c_Fm.write_to_file(this->file_window,this->kvector_window,this->pk_w);

    }

  else if(statistics=="Pk_ys" || statistics=="Pk_yb"  || statistics=="Pk_ybc" || statistics=="Pk_ysc" || statistics=="Pk_y_ds")
    {
      c_Fm.write_to_file(file_power,this->kvector_data,this->pk0,this->pk2,this->pk4,this->sigma_fkp,this->modes_g);
    }
  else if(statistics=="Bk_fkp")
    c_Fm.write_to_file(file_bispectrum,kvector_data_b,bispectrum,sn_bispectrum, modes_tri);
}



// *************************************************************************************
// *************************************************************************************

void PowerSpectrumF::write_power_and_modes()
{
   this->So.enter(__PRETTY_FUNCTION__);
   c_Fm.write_to_file2(this->file_power,this->kvector_data,this->pk0,this->modes_g,true);
   if(this->use_random_catalog)
      c_Fm.write_to_file2(this->file_window,this->kvector_data,this->pk_w,this->modes_g,true);
}
// *************************************************************************************
void PowerSpectrumF::write_power_and_modes(string file)
{
    this->So.enter(__PRETTY_FUNCTION__);
  c_Fm.write_to_file2(file,this->kvector_data,this->pk0,this->modes_g, true);
}
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************

void PowerSpectrumF::write_power_spectrum_grid(string output_file)
{
  this->So.enter(__PRETTY_FUNCTION__);
#ifdef _FULL_VERBOSE_
    So.message_screen("Writing outputs");
#endif
    c_Fm.write_to_file(output_file,this->kvector_data,this->pk0,this->pk2,this->pk4,this->sigma_fkp,this->modes_g);


}
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************

void PowerSpectrumF::compute_marked_correlation_function()
{

  this->So.enter(__PRETTY_FUNCTION__);
  time_t start;
  time (&start);

  cout.precision(12);
  cout<<CYAN<<"Measuring marked correlation function. SO FAR INPUTS ARE IN CARTESSIAN COORDINATES"<<RESET<<endl;
  cout<<CYAN<<"The info of the weight is taken from column "<<i_weight1_g<<" in the input catalog "<<RESET<<endl;
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
 this->add_catalogues(0);
#else
  this->add_catalogues(0,1e20);
#endif

  if(this->gc_n_columns<3)cout<<RED<<" Catalog with no extra information"<<RESET<<endl;

  vector<int>count(Nbins_r,0);
  vector<real_prec>mcount(Nbins_r,0);
  vector<real_prec> rbin(Nbins_r,0);

  real_prec Deltar=r_bin_type == "LIN"  ?   (rmax-rmin)/((real_prec)Nbins_r):  (log10(rmax/rmin))/((real_prec)Nbins_r);

  int NTHREADS =  omp_get_max_threads();
  cout<<CYAN<<"Measuring Marked Correlation function using "<<NTHREADS<<" threads ..."<<RESET<<endl;

  if(r_bin_type=="LIN")for (int i=0;i<rbin.size();++i)rbin[i]=rmin+(i+0.5)*Deltar;
  else for (int i=0;i<rbin.size();++i)rbin[i]=pow(10, log10(rmin)+(i+0.5)*Deltar);

  real_prec mean_mark=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:mean_mark)
#endif
  for (int i=0;i<this->N_galaxy;++i)
    mean_mark+=this->tracer_cat.Halo[i].weight1;

  mean_mark=(static_cast<double>(mean_mark))/(static_cast<double>(this->N_galaxy));
  cout<<CYAN<<"Mean mark = "<<BLUE<<mean_mark<<RESET<<endl;


#pragma omp parallel num_threads(NTHREADS)
  {

    vector<int> count_priv(Nbins_r,0);
    vector<real_prec> mcount_priv(Nbins_r,0);

#pragma omp for nowait
    for (int i=0;i<this->N_galaxy;++i)
      {
        real_prec x=this->tracer_cat.Halo[i].coord1;
        real_prec y=this->tracer_cat.Halo[i].coord2;
        real_prec z=this->tracer_cat.Halo[i].coord3;
        real_prec weight=this->tracer_cat.Halo[i].weight1;
	for (int j=i+1;j<this->N_galaxy;++j)
	  {
            real_prec xp=this->tracer_cat.Halo[j].coord1;
            real_prec yp=this->tracer_cat.Halo[j].coord1;
            real_prec zp=this->tracer_cat.Halo[j].coord1;
            real_prec weightp=this->tracer_cat.Halo[j].weight1;

            real_prec r=sqrt(pow(xp-x,2)+pow(yp-y,2)+pow(zp-z,2));
	    if(r<rmax && r>=rmin)
	      {
		int ind = r_bin_type == "LIN" ? (int)floor((r-rmin)/Deltar)  :  (int)floor((log10(r/rmin))/Deltar);
		count_priv[ind]++;
		mcount_priv[ind]+=(weight*weightp);
	      }
	  }
      }

#pragma critical
    {
      for(int i=0;i<mcount.size();++i)mcount[i]+= mcount_priv[i];
      for(int i=0;i<count.size() ;++i) count[i]+= count_priv[i];
    }
  }

  this->tracer_cat.Halo.clear();
  this->tracer_cat.Halo.shrink_to_fit();

  So.DONE();
  cout<<CYAN<<" Done"<<RESET<<endl;

  for(int i=0;i<mcount.size();++i)mcount[i]=(mcount[i]/( pow(mean_mark,2)*(real_prec)count[i]));
  for(int i=0;i<mcount.size();++i)cout<<rbin[i]<<"  "<<mcount[i]<<"  "<<count[i]<<endl;

  string ov_name= file_MCF;
  ofstream mcor; mcor.open(ov_name.c_str());
  mcor.precision(6);
  for(int i=0;i<mcount.size();++i)mcor<<rbin[i]<<"\t"<<mcount[i]<<"\t"<<count[i]<<endl;
  mcor.close();

  cout<<CYAN<<"Output in file "<<GREEN<<ov_name<<RESET<<endl;
  cout<<CYAN<<"Done"<<endl;

}


// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
void PowerSpectrumF::compute_cross_power_spectrum_grid(bool dm, vector<real_prec>&X,vector<real_prec>&Y)
{

    this->So.enter(__PRETTY_FUNCTION__);

    int NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);

    FftwFunctions c_Ff(this->params);
    s_parameters_box s_p_box;
    s_p_box.mas=this->mass_assignment_scheme;
    s_p_box.ave=this->type_of_binning;
    s_p_box.k_bin_step=this->k_bin_step;
    s_p_box.use_MAS_correction=this->MAS_correction;
    s_p_box.FKP_weight=this->FKP_weight;
    s_p_box.FKP_error_bars=this->FKP_error_bars;
    s_p_box.FKP_error_bars_exact= this->FKP_error_bars_exact;
    s_p_box.use_SN_correction=this->SN_correction;
    s_p_box.Pest=this->Pest;
    c_Ff.Lside=this->Lbox;
    c_Ff.NT = this->Nft*this->Nft*this->Nft;
    c_Ff.Nnp_data = this->Nft/2/this->ndel_data;
    real_prec ngal_new=0;

    c_Ff.data_g.clear();
    c_Ff.data_g.shrink_to_fit();
    c_Ff.data_g.resize(c_Ff.NT,0);

    if (input_type=="delta_grid")
      {
#pragma omp parallel
        for(ULONG i=0;i<c_Ff.NT;++i)
          c_Ff.data_g[i]=static_cast<real_prec>(Y[i]);
        ngal_new=this->ngal_delta;
      }
    else if (input_type=="density_grid")
        {
#pragma omp parallel for reduction(+:ngal_new)
          for(ULONG i=0;i<c_Ff.NT;++i)
            ngal_new+=static_cast<real_prec>(Y[i]);
          real_prec nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(c_Ff.NT);

#pragma omp parallel for
          for(ULONG i=0;i<c_Ff.NT;++i)
            c_Ff.data_g[i]=(static_cast<real_prec>(Y[i])/static_cast<real_prec>(nmean))-1.;
          this->ngal_delta=ngal_new;
        }
    if(true==this->SN_correction)
      c_Ff.shot_noise=static_cast<real_prec>(pow(c_Ff.Lside,3))/static_cast<real_prec>(ngal_new);
    real_prec factor=pow(static_cast<real_prec>(c_Ff.Lside),1.5)/static_cast<real_prec>(c_Ff.data_g.size());
    c_Ff.normal_power=pow(factor,-2);


    c_Ff.data_gp.clear();
    c_Ff.data_gp.shrink_to_fit();
    c_Ff.data_gp.resize(c_Ff.NT,0);
    ngal_new=0;
    if (input_type=="delta_grid")
#ifdef _USE_OMP_
#pragma omp parallel
#endif
        for(ULONG i=0;i<c_Ff.NT;++i)
          c_Ff.data_gp[i]=static_cast<real_prec>(X[i]);

    else if (input_type=="density_grid")
        {
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ngal_new)
#endif
          for(ULONG i=0;i<c_Ff.NT;++i)
            ngal_new+=static_cast<real_prec>(X[i]);
          real_prec nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(c_Ff.NT);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
           for(ULONG i=0;i<c_Ff.NT;++i)
            c_Ff.data_gp[i]=(static_cast<real_prec>(X[i])/static_cast<real_prec>(nmean))-1.;
        }

      c_Ff.shot_noise2=0;
      factor=pow(static_cast<real_prec>(c_Ff.Lside),1.5)/static_cast<real_prec>(c_Ff.data_g.size());
      c_Ff.normal_power_two=pow(factor,-2);


      if(false==dm && true==this->SN_correction)
          c_Ff.shot_noise2=static_cast<real_prec>(pow(c_Ff.Lside,3))/static_cast<real_prec>(ngal_new);


         c_Ff.set_pars(this->Lbox,kmax_y_ds, ndel_data, ndel_window, N_log_bins,N_mu_bins,MAS_correction,mass_assignment_scheme,kmin_bk,kmax_bk);
         c_Ff.set_bins(type_of_binning);
#ifdef _FULL_VERBOSE_
         So.write_fftw_parameters((void *)&s_p_box);
#endif
         //          c_Ff.write_fftw_parameters();
          c_Ff.fftw_vectors(use_random_catalog);



          kvector_data.clear();
          kvector_data.shrink_to_fit();
          for(int i=0;i<c_Ff.Nnp_data;i++)
            kvector_data.push_back(c_Ff.DeltaK_data*(i+k_bin_step));


          this->pk0.clear();
          this->pk0.shrink_to_fit();
          this->pk0.resize(c_Ff.Nnp_data,0); //Monopole
          this->modes_g.clear();
          this->modes_g.resize(c_Ff.Nnp_data,0); //Monopole
          c_Ff.cross_power_spectrum_fkp(&s_p_box, this->pk0,this->modes_g);

}



// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
void PowerSpectrumF::compute_cross_power_spectrum_grid(bool dm,string file_X, string file_Y)
{
    this->So.enter(__PRETTY_FUNCTION__);

    int NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);
    FileOutput File;
    FftwFunctions c_Ff(this->params);
    s_parameters_box s_p_box;
    s_p_box.mas=this->mass_assignment_scheme;
    s_p_box.ave=this->type_of_binning;
    s_p_box.k_bin_step=this->k_bin_step;
    s_p_box.use_MAS_correction=this->MAS_correction;
    s_p_box.FKP_weight=this->FKP_weight;
    s_p_box.FKP_error_bars=this->FKP_error_bars;
    s_p_box.FKP_error_bars_exact= this->FKP_error_bars_exact;
    s_p_box.use_SN_correction=this->SN_correction;
    s_p_box.Pest=this->Pest;
    c_Ff.Lside=this->Lbox;
    c_Ff.NT = this->Nft*this->Nft*this->Nft;
    c_Ff.Nnp_data = this->Nft/2/this->ndel_data;
    real_prec ngal_new=0;


    vector<float>X(c_Ff.NT,0);
    File.read_array_t<float>(file_X, X);

    vector<real_prec>Y(c_Ff.NT,0);
    File.read_array_t<PrecType_Y>(file_Y, Y);


    c_Ff.data_g.clear();
    c_Ff.data_g.shrink_to_fit();
    c_Ff.data_g.resize(c_Ff.NT,0);



    if (input_type=="delta_grid")
      {
#pragma omp parallel
        for(ULONG i=0;i<c_Ff.NT;++i)
          c_Ff.data_g[i]=static_cast<real_prec>(Y[i]);
        ngal_new=this->ngal_delta;
      }
    else if (input_type=="density_grid")
      {
#pragma omp parallel for reduction(+:ngal_new)
	for(ULONG i=0;i<c_Ff.NT;++i)
	  ngal_new+=static_cast<real_prec>(Y[i]);
	real_prec nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(c_Ff.NT);

#pragma omp parallel for
	for(ULONG i=0;i<c_Ff.NT;++i)
	  c_Ff.data_g[i]=(static_cast<real_prec>(Y[i])/static_cast<real_prec>(nmean))-1.;
	this->ngal_delta=ngal_new;
      }


    if(true==this->SN_correction)
      c_Ff.shot_noise=static_cast<real_prec>(pow(c_Ff.Lside,3))/static_cast<real_prec>(ngal_new);

    real_prec factor=pow(static_cast<real_prec>(c_Ff.Lside),1.5)/static_cast<real_prec>(c_Ff.data_g.size());
    c_Ff.normal_power=pow(factor,-2);

    Y.clear(); Y.shrink_to_fit();



    c_Ff.data_gp.clear();
    c_Ff.data_gp.shrink_to_fit();
    c_Ff.data_gp.resize(c_Ff.NT,0);
    ngal_new=0;

    for(ULONG i=0;i<c_Ff.NT;++i)
      c_Ff.data_gp[i]=static_cast<real_prec>(X[i]);


    if (input_type=="delta_grid")
#ifdef _USE_OMP_
#pragma omp parallel
#endif
        for(ULONG i=0;i<c_Ff.NT;++i)
          c_Ff.data_gp[i]=static_cast<real_prec>(X[i]);

    else if (input_type=="density_grid")
      {
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ngal_new)
#endif
	for(ULONG i=0;i<c_Ff.NT;++i)
	  ngal_new+=static_cast<real_prec>(X[i]);
	real_prec nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(c_Ff.NT);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
	for(ULONG i=0;i<c_Ff.NT;++i)
	  c_Ff.data_gp[i]=(static_cast<real_prec>(X[i])/static_cast<real_prec>(nmean))-1.;
      }

    c_Ff.shot_noise2=0;
    factor=pow(static_cast<real_prec>(c_Ff.Lside),1.5)/static_cast<real_prec>(c_Ff.data_g.size());
    c_Ff.normal_power_two=pow(factor,-2);

    X.clear(); X.shrink_to_fit();



    if(false==dm && true==this->SN_correction)
      c_Ff.shot_noise2=static_cast<real_prec>(pow(c_Ff.Lside,3))/static_cast<real_prec>(ngal_new);


    c_Ff.set_pars(this->Lbox,kmax_y_ds, ndel_data, ndel_window, N_log_bins,N_mu_bins,MAS_correction,mass_assignment_scheme,kmin_bk,kmax_bk);
    c_Ff.set_bins(type_of_binning);
#ifdef _FULL_VERBOSE_
    So.write_fftw_parameters((void *)&s_p_box);
#endif
    c_Ff.write_fftw_parameters();
    c_Ff.fftw_vectors(use_random_catalog);



    kvector_data.clear();
    kvector_data.shrink_to_fit();
    for(int i=0;i<c_Ff.Nnp_data;i++)
      kvector_data.push_back(c_Ff.DeltaK_data*(i+k_bin_step));

    this->pk0.clear();
    this->pk0.shrink_to_fit();
    this->pk0.resize(c_Ff.Nnp_data,0); //Monopole
    this->modes_g.clear();
    this->modes_g.resize(c_Ff.Nnp_data,0); //Monopole


      c_Ff.cross_power_spectrum_fkp(&s_p_box, this->pk0,this->modes_g);
      write_power_and_modes();
      cout<<"ebony"<<endl;
}

// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
#ifdef _USE_POWER_IN_BAM_
void PowerSpectrumF::compute_power_spectrum(bool verbose, vector<s_Halo>& tracer_cat){    
#else
void PowerSpectrumF::compute_power_spectrum(bool verbose, bool mcut){
#endif

  So.enter(__PRETTY_FUNCTION__);


#ifdef _USE_OMP_
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

  #ifdef _MASS_WEIGHT_POWER_
   if(true==this->params._weight_with_mass())
     So.message_screen("Measuring the Mass-power spectrum");
#elif defined(_USE_MASS_CUTS_PK_)
   So.message_screen("Measuring power spectrum in mass cuts");
#endif

  string file_pow=this->file_power;

#ifdef _USE_MASS_CUTS_PK_
  vector<real_prec>mass_cuts;
#elif defined (_USE_MASS_BINS_PK_)
  vector<real_prec>mass_bin_min;
  vector<real_prec>mass_bin_max;
#endif


  int im=0;
#ifdef _USE_MASS_CUTS_PK_
  if(true==mcut)
    {
#ifdef _USE_MASS_AS_OBSERVABLE_POWER_
    So.message_warning("Mass-cuts are defined in PowerSpectrumF::compute_power_spectrum");
     mass_cuts.push_back(9e11);
     mass_cuts.push_back(3e12);
     mass_cuts.push_back(6e12);
     mass_cuts.push_back(1e13);
     mass_cuts.push_back(3e13);
     mass_cuts.push_back(6e13);
#elif defined _USE_VMAX_AS_OBSERVABLE_POWER_
     So.message_warning("VMAX-cuts are defined in PowerSpectrumF::compute_power_spectrum");
     mass_cuts.push_back(this->params._VMAXmin());
     mass_cuts.push_back(300);
     mass_cuts.push_back(400);
     mass_cuts.push_back(600);
     mass_cuts.push_back(700);
     mass_cuts.push_back(900);
#endif
    }
  else
    {
      real_prec mcut_aux=MINIMUM_PROP_CUT;
#ifdef _USE_MASS_AS_OBSERVABLE_POWER_
      So.message_screen("Power spectrum mesasured with one mass-cut at" , mcut_aux, "Ms/h");
#elif defined _USE_VMAX_AS_OBSERVABLE_POWER_
      So.message_screen("Power spectrum mesasured with one VMAX cut at" , mcut_aux, "km/s");
#endif
      if(true==this->params._weight_with_mass())
      {
        So.message_screen("Mass weighted power spectrum");
        mcut_aux=pow(10,this->params._LOGMASSmax());
      }
      mass_cuts.push_back(mcut_aux);
    }//closes else


  for(im=0; im< mass_cuts.size();++im)// meant for cuts
    {


#elif defined (_USE_MASS_BINS_PK_)
#ifdef _FULL_VERBOSE_
      So.message_warning("Property-bins are defined in PowerSpectrumF::compute_power_spectrum");
#endif
      int N_intervals=this->params._NMASSbins_power()+1;
      mass_bin_min.resize(N_intervals,0);
      mass_bin_max.resize(N_intervals,0);
#ifdef _USE_MASS_AS_OBSERVABLE_POWER_
      real_prec l_minp=this->params._LOGMASSmin();
      real_prec l_maxp=this->params._LOGMASSmax();
#elif defined _USE_VMAX_AS_OBSERVABLE_POWER_
      real_prec l_minp=log10(this->params._VMAXmin());
      real_prec l_maxp=log10(this->params._VMAXmax());
#elif defined _USE_RS_AS_OBSERVABLE_POWER_
      real_prec l_minp=log10(this->params._RSmin());
      real_prec l_maxp=log10(this->params._RSmax());
#elif defined _USE_SPIN_AS_OBSERVABLE_POWER_
      real_prec l_minp=log10(this->params._SPINmin());
      real_prec l_maxp=log10(this->params._SPINmax());
#endif

      mass_bin_min[0]=pow(10,l_minp);  // Let bin 0 be the full sample
      mass_bin_max[0]=pow(10,l_maxp);
      real_prec ldelta_prop=(l_maxp-l_minp)/static_cast<real_prec>(N_intervals-1);


  for(int im=1;im<N_intervals;++im)// works for bins
     {
       mass_bin_min[im]=pow(10,l_minp+(im-1)*ldelta_prop);
       mass_bin_max[im]=pow(10,l_minp+(im)*ldelta_prop);
#ifdef _FULL_VERBOSE_
       cout<<"Bin "<<im<<"  :"<<mass_bin_min[im]<<"  "<<mass_bin_max[im]<<endl;
#endif
     }

   for(im=0; im<N_intervals;++im) // meant for bins
     {
#endif



      //FftwFunctions c_Ff(statistics,Nft, measure_cross); //old declaration of class object
      FftwFunctions c_Ff(this->params);

      c_Ff.imcut=im;

#ifdef _USE_MASS_CUTS_PK_
      this->file_power=file_pow+"_masscut"+to_string(im);
#elif defined (_USE_MASS_BINS_PK_)
#ifdef _USE_MASS_AS_OBSERVABLE_POWER_
      this->file_power=file_pow+"_massbin"+to_string(im);
#elif defined (_USE_VMAX_AS_OBSERVABLE_POWER_)
      this->file_power=file_pow+"_vmaxbin"+to_string(im);
#elif defined (_USE_RS_AS_OBSERVABLE_POWER_)
      this->file_power=file_pow+"_rsbin"+to_string(im);
#elif defined (_USE_SPIN_AS_OBSERVABLE_POWER_)
      this->file_power=file_pow+"_spinbin"+to_string(im);
#endif
#endif


#ifdef _FULL_VERBOSE_
  if(true==verbose)
    {
      if(statistics=="Pk_fkp")So.welcome_message();
      if(statistics=="Bk_fkp")So.welcome_message_bispectrum();
      if(statistics=="Bk_fkp_fast")So.welcome_message_bispectrum_fast();
      if(statistics=="Pk_ys" || statistics=="Pk_yb" || statistics=="Pk_ybc" || statistics=="Pk_ysc" || statistics=="Pk_y_ds")So.welcome_message_yama();
    }
#endif
  // ***********************************************************
  // Load structure with  cosmological parameters

  s_CosmologicalParameters s_cosmo_par;
  s_cosmo_par.Om_matter=om_matter;
  s_cosmo_par.Om_baryons=om_baryons;
  s_cosmo_par.Om_radiation=om_radiation;
  s_cosmo_par.Om_k=om_k;
  s_cosmo_par.Om_vac=om_vac;
  s_cosmo_par.Hubble=Hubble;
  s_cosmo_par.hubble=hubble;
  s_cosmo_par.w_eos =w_eos;
  s_cosmo_par.N_eff=N_eff;
  s_cosmo_par.n_s=spectral_index;
  s_cosmo_par.sigma8=sigma8;
  s_cosmo_par.Tcmb=Tcmb;



  if(sys_of_coord_g==2)
    So.write_cosmo_parameters((void *)&s_cosmo_par);
  // ***********************************************************

  // Define some structures here
  s_parameters_box s_p_box;
  s_data_structure s_data_struct_r; // this is defiend here for the fkp error bars need the info from the randoms. So far at this point it is only defiend but not filled
  //
  // ***********************************************************
  // Add random and galaxy catalogues
  this->c_Fm.input_type=this->input_type;

//#if defined (_USE_MASS_CUTS_PK_) ||  defined (_USE_ALL_PK_)
  if(this->input_type!="catalog")
    So.message_warning("Warning: expecting catalog to perform mass cuts when density field is provided");

  real_prec mm=0;


#ifdef _USE_POWER_IN_BAM_

  this->use_random_catalog=false;
// Add random and galaxy catalogues
  this->c_Fm.input_type=this->input_type;
  this->N_galaxy=tracer_cat.size();

#else  // else for ifdef _USE_POWER_IN_BAM_


  if("catalog"==this->input_type)
#ifdef _USE_ALL_PK_
    if(this->params._i_mass_g()>0)mm=pow(10,this->params._LOGMASSmin());
    else mm=0;
#else
    mm=mass_cuts[im];
#endif


#ifdef _USE_ALL_PK_
  if("catalog"==this->input_type)
    this->add_catalogues(mm);
#elif defined (_USE_MASS_BINS_PK_)
  this->add_catalogues(mass_bin_min[im],mass_bin_max[im]);
#endif

#endif


  if(this->input_type=="catalog")
    {


      // Tabulate r-z relation if radial coordinate is redshift
      vector<gsl_real>vzz(N_z_bins,0);
      vector<gsl_real>vrc(N_z_bins,0);
      So.message_screen("Computing comovong distances");
      if(sys_of_coord_g==1 ||sys_of_coord_g==2 )
        c_Cf.Comoving_distance_tabulated(Z_MIN,Z_MAX, (void *)&s_cosmo_par,vzz,vrc);
      So.DONE();
      // **************************************************************************
	  // Derived parameters associated to
	  //  c_Ff.set_healpix_pars(Healpix_resolution);

	  // **************************************************************************
	  // If applies, give a first estimate of mean number density from a box
	  // in case we do not use random catalog                                                    *
	  //If a random catalog is used, set this numer to 1
	  //and use the nmbar tabulated in the catalogs or computed in this code
	  mean_density=1.0;
    if(false==use_random_catalog)
	    {
	      So.message_screen("Using particles in a box.");
        mean_density=static_cast<real_prec>(this->N_galaxy)/pow(Lbox,3);
	      So.message_screen("Mean Number density = ",mean_density," (Mpc/h)^(-3)");
	    }


	  // *************************************************************************
	  // If nbar is not tabulated, compute it from the random catalog
     bool compute_dndz=false;
     if(true==use_random_catalog && false==nbar_tabulated)
      compute_dndz=true;

	  // *************************************************************************
	  // Give a first estimate of the alpha-parameter                                            *
     real_prec alpha_0=1.0;
     if(true==use_random_catalog)
       alpha_0 = ((real_prec)this->N_galaxy)/((real_prec)this->N_random);


	  // This structure will be used in 1P_statistics
    s_dndz s_dndz_params;
    if(use_random_catalog==true && nbar_tabulated==false)
      s_dndz_params={
       this->random_cat.Halo,
       rc_n_columns,
       alpha_0,
       nbar_tabulated,
       compute_dndz,
       constant_depth,
       vzz,
       vrc,
       N_dndz_bins,
       new_N_dndz_bins,
       redshift_min_sample,
       redshift_max_sample,
       area_survey,
       sys_of_coord_r,
       i_coord1_r,
       i_coord2_r,
       i_coord3_r,
       c_Ff.area_pixel,
       c_Ff.npixels,
       c_Ff.nside,
       file_dndz
      };



	  // This should only be defined if we do not use randoms with nbar tabulated.
	  // IN FFTWFUNCTIONS ALL CALLS TO HEALPIX AND MAP-DEFINED VECTORS ARE COMMENTED
     vector<gsl_real> z_v;
     vector<gsl_real> dndz_v;
     vector< vector<gsl_real> > dndz_matrix;
	  // If nbar is not tabulated, Compute a smoothed version of dN/dz from randoms to get nbar  *
     if(true==use_random_catalog && false==nbar_tabulated)
	    {
	      z_v.resize(new_N_dndz_bins,0);
	      dndz_v.resize(z_v.size(),0);
	      dndz_matrix.resize(z_v.size());
	      for(int i=0;i<dndz_matrix.size(); ++i)dndz_matrix[i].resize(c_Ff.npixels,0);
	      c_Op.dndz((void *)&s_dndz_params,z_v, dndz_v, dndz_matrix);
	    }
	  // *************************************************************************
     So.message_screen("Setting data in structure for TRACERS");
     s_data_structure s_data_struct_g={
     this->tracer_cat.Halo,
     gc_n_columns,
     sys_of_coord_g,
	    "data",
	    mean_density,
	    nbar_tabulated,
	    compute_dndz,
	    z_v,
	    dndz_v,
	    dndz_matrix
	  };
    So.DONE();
#ifdef _FULL_VERBOSE_
    So.message_screen("Freeing memmory prop in ", __PRETTY_FUNCTION__);
#endif
   this->tracer_cat.Halo.clear();
   this->tracer_cat.Halo.shrink_to_fit();
   So.DONE();

	  // *************************************************************************
	  // Allocate structure with random catalogue and information of dN/dz                       *
	  // *************************************************************************
#ifdef _FULL_VERBOSE_
    So.message_screen("Setting data in structure for RANDOMS");
#endif
    s_data_structure s_data_struct_r={
    this->random_cat.Halo,
    rc_n_columns,
    sys_of_coord_r,
    "random",
    mean_density,
    nbar_tabulated,
    compute_dndz,
    z_v,
    dndz_v,
    dndz_matrix
	  };
   So.DONE();
#ifdef _FULL_VERBOSE_
   So.message_screen("Freeing memmory prop in ", __PRETTY_FUNCTION__);
#endif
   this->random_cat.Halo.clear();
   this->random_cat.Halo.shrink_to_fit();
    So.DONE();

	  // *****************************************************************************************
	  // Allocate a strucuture
	  // to be passed to the class member functions. We pass only those parameters
	  // that were not passed through the set_pars method.

          s_p_box={
	    use_random_catalog,
	    mass_assignment_scheme,
	    type_of_binning,
	    k_bin_step,
	    MAS_correction,
	    FKP_weight,
	    FKP_error_bars,
	    FKP_error_bars_exact,
	    SN_correction,
	    Pest,
	    nbar_tabulated,
	    compute_dndz,
	    constant_depth,
	    vzz,
	    vrc,
	    N_dndz_bins,
	    new_N_dndz_bins,
	    redshift_min_sample,
	    redshift_max_sample,
	    area_survey,
	    sys_of_coord_r,
	    i_coord1_r,
	    i_coord2_r,
	    i_coord3_r,
	    i_weight1_r,
	    i_weight2_r,
	    i_weight3_r,
	    i_weight4_r,
	    use_weight1_r,
	    use_weight2_r,
	    use_weight3_r,
	    use_weight4_r,
	    i_mean_density_r,
	    angles_units_r,
	    i_coord1_g,
	    i_coord2_g,
	    i_coord3_g,
	    i_mass_g,
	    i_weight1_g,
	    i_weight2_g,
	    i_weight3_g,
	    i_weight4_g,
	    use_weight1_g,
	    use_weight2_g,
	    use_weight3_g,
	    use_weight4_g,
	    i_mean_density_g,
	    angles_units_g,
	    c_Ff.area_pixel,
	    c_Ff.npixels,
	    c_Ff.nside,
	    file_dndz,
	    new_los,
	    kmin_bk,
	    kmax_bk,
	    use_fundamental_mode_as_kmin_bk,
	    false
	  };


	  // *****************************************************************************************
	  // Transforming to cartessian coord.- and searching for box side lenght.
          if(this->sys_of_coord_g!=0)
            {
              So.message_screen("Transform Coordinates in tracer catalogue to cartesian coordinates");
              c_Ff.cart_coordinates(&s_p_box,&s_data_struct_g);
              So.DONE();
           }

          if(true==this->use_random_catalog)
	    {
              if(this->sys_of_coord_r!=0)
               {
                 So.message_screen("Transform Coordinates in random catalogue to cartesian coordinates");
                 c_Ff.cart_coordinates(&s_p_box,&s_data_struct_r);
                 So.DONE();
              }
	    }

	  // *****************************************************************************************
	  // Determine size of box for the density interpolation and the Fourier transform:
	  // New version: we do this in the same loop that convert to cartesian coords.
	  // We then save other loop over the galaxies!
	  // The side of the box computed from the catalog is a public variable of FFTW Lside_data.
	  // Given the fact that we might want to have it fixed from the parameter file
	  // we select the new Lside and pass it again as Lside
          real_prec Lside=Lbox;
          if(true==new_Lbox)
	    Lside=c_Ff.Lside_data;

	  // *****************************************************************************************
	  c_Ff.set_pars(Lside,kmax_y_ds, ndel_data, ndel_window, N_log_bins,N_mu_bins,MAS_correction,mass_assignment_scheme,kmin_bk,kmax_bk);

	  // **********************************************************************************
	  // Define binning
          c_Ff.set_bins(type_of_binning);

	  // *****************************************************************************************
	  // FFTW ARRAYS : must be *AFTER* set_pars.
          c_Ff.fftw_vectors(this->use_random_catalog);

	  // **********************************************************************************
	  // Estimate of the mean number density
	  //  cout<<"Maximum Nfft allowed by mean density = "<<(int)(2.*pow(N_galaxy, 1./3.))<<endl;
          mean_density=this->N_galaxy/pow(Lside,3);  //This is  raw estimate!!
          if(false==use_random_catalog)
              mean_density=this->N_galaxy/pow(Lside,3);
          else mean_density=1.0;

          if(true==verbose)
              {
#ifdef _FULL_VERBOSE_
                  So.write_fftw_parameters((void *)&s_p_box);
#endif
                  c_Ff.write_fftw_parameters();
             }
	  // ***********************************************************************************
	  // Build interpolated galaxy density field
	  if(statistics=="Pk_y_ds")
	    {
              cout<<endl;
	      So.message_screen("Creating galaxy density field on a Fourier grid...");
	      c_Ff.get_power_moments_fourier_grid_ds_yam(&s_p_box,&s_data_struct_g);
	      So.DONE();
	    }
	  else
	    {
              cout<<endl;
              So.message_screen("Interpolating galaxy density field on a grid");
#ifdef _USE_VECTORIZED_GRID_ASSIGNMENT_
              c_Ff.get_interpolated_density_field_old(&s_p_box,&s_data_struct_g);
#else
              c_Ff.get_interpolated_density_field(&s_p_box,&s_data_struct_g);

#endif

              So.DONE();
	    }
	  // ***********************************************************************************
	  // Build interpolated random density field
          if(true==use_random_catalog)
	    {
	      if(statistics=="Pk_y_ds")
		{
                  cout<<endl;
		  So.message_screen("Creating random density field on a Fourier grid");
		  c_Ff.get_power_moments_fourier_grid_ds_yam(&s_p_box,&s_data_struct_r);
		  So.DONE();
		}
	      else
		{
		  So.message_screen("Interpolating random density field on a grid");
#ifdef _USE_VECTORIZED_GRID_ASSIGNMENT_
                  c_Ff.get_interpolated_density_field_old(&s_p_box,&s_data_struct_r);
#else
                  c_Ff.get_interpolated_density_field(&s_p_box,&s_data_struct_r);
#endif

                  So.DONE();
		}

	    }
	  else
	    {
              real_prec vol=pow(Lside,3);
              c_Ff.raw_sampling(vol);
	    }

	  // **********************************************************************************
          c_Ff.get_parameters_estimator(use_random_catalog, verbose);
	  // **********************************************************************************
	  // Build fluctuation                                                                       *
          c_Ff.get_fluctuation(this->use_random_catalog);

	}


      else if (input_type=="delta_grid" || input_type=="density_grid"  ) // Else, we read the delta from this input file.
	{

	  So.message_screen("Starting density field on a grid");

	  s_p_box.mas=mass_assignment_scheme;
	  s_p_box.ave=type_of_binning;
	  s_p_box.k_bin_step=k_bin_step;
	  s_p_box.use_MAS_correction=MAS_correction;
	  s_p_box.FKP_weight=FKP_weight;
	  s_p_box.FKP_error_bars=FKP_error_bars;
	  s_p_box.FKP_error_bars_exact= FKP_error_bars_exact;
	  s_p_box.use_SN_correction=SN_correction;
	  s_p_box.Pest=Pest;
	  c_Ff.Lside=this->Lbox;
	  c_Ff.NT = Nft*Nft*Nft;

          real_prec ngal_new;
	  // *****************************************************************************************

	  bool measure_diff=false;



	  if(false==measure_cross && false==measure_diff) // If no crossed power, read per default thedelta_grid_file
	    {
	      fftw_array<float> dummy(c_Ff.NT);

	      switch(this->measure_cross_from_1)
		{ //choose the file to get the auto power from
		case(1):
		  this->c_Fm.read_array(this->delta_grid_file, dummy,c_Ff.NT);
		  break;
		case(2):
		  this->c_Fm.read_array(this->delta_grid_file2, dummy,c_Ff.NT);
		  break;
		case(3):
		  this->c_Fm.read_array(this->delta_grid_file3, dummy,c_Ff.NT);
		  break;
		case(4):
		  this->c_Fm.read_array(this->delta_grid_file4, dummy,c_Ff.NT);
		  break;
		}

	      c_Ff.data_g.clear();
	      c_Ff.data_g.shrink_to_fit();
	      c_Ff.data_g.resize(c_Ff.NT,0);

	      if (input_type=="delta_grid")
                {
#ifdef _USE_OMP_
#pragma omp parallel
#endif
		  for(ULONG i=0;i<c_Ff.NT;++i)
                    c_Ff.data_g[i]=static_cast<real_prec>(dummy[i]);

                  ngal_new=this->ngal_delta;
		}
              else if (input_type=="density_grid")
		  {
		    ngal_new=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ngal_new)
#endif
                    for(ULONG i=0;i<c_Ff.NT;++i)
                      ngal_new+=static_cast<real_prec>(dummy[i]);
                    real_prec nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(c_Ff.NT);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
		    for(ULONG i=0;i<c_Ff.NT;++i)
                      c_Ff.data_g[i]=(static_cast<real_prec>(dummy[i])/static_cast<real_prec>(nmean))-1.;
		    this->ngal_delta=ngal_new;
		    c_Ff.n_gal =ngal_new;
		  }

              real_prec factor=pow(static_cast<real_prec>(c_Ff.Lside),1.5)/static_cast<real_prec>(c_Ff.data_g.size());
              c_Ff.normal_power=pow(factor,-2);
              c_Ff.shot_noise=static_cast<real_prec>(pow(c_Ff.Lside,3))/static_cast<real_prec>(ngal_new);
	      c_Ff.set_pars(this->Lbox,kmax_y_ds, ndel_data, ndel_window, N_log_bins,N_mu_bins,MAS_correction,mass_assignment_scheme,kmin_bk,kmax_bk);
	      c_Ff.set_bins(type_of_binning);
              if(true==verbose)
                {
#ifdef _FULL_VERBOSE_
                 So.write_fftw_parameters((void *)&s_p_box);
#endif
                 c_Ff.write_fftw_parameters();
#ifdef _FULL_VERBOSE_
                 So.message_screen("Shot Noise =",c_Ff.shot_noise);
#endif
               }
              c_Ff.fftw_vectors(use_random_catalog);
	    }

	  else if(true==measure_diff)
	    {

	      fftw_array<float> dummy(c_Ff.NT);
	      fftw_array<float> dummy2(c_Ff.NT);

	      this->c_Fm.read_array(this->delta_grid_file, dummy,c_Ff.NT);
	      this->c_Fm.read_array(this->delta_grid_file2, dummy2,c_Ff.NT);

	      c_Ff.data_g.resize(c_Ff.NT,0);


	      ULONG ngal_new2=0;
	      ULONG ngal_new1=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ngal_new1)
#endif
	      for(ULONG i=0;i<c_Ff.NT;++i)
                ngal_new1+=static_cast<real_prec>(dummy[i]);


#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ngal_new2)
#endif
	      for(ULONG i=0;i<c_Ff.NT;++i)
                ngal_new2+=static_cast<real_prec>(dummy2[i]);


              real_prec nmean1=static_cast<real_prec>(ngal_new1)/static_cast<real_prec>(c_Ff.NT);
              real_prec nmean2=static_cast<real_prec>(ngal_new2)/static_cast<real_prec>(c_Ff.NT);



              c_Ff.shot_noise=static_cast<real_prec>(pow(c_Ff.Lside,3))/static_cast<real_prec>(ngal_new1);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
	      for(ULONG i=0;i<c_Ff.NT;++i)
                //	    c_Ff.data_g[i]=(static_cast<real_prec>(dummy[i])/static_cast<real_prec>(nmean1) - static_cast<real_prec>(dummy2[i])/static_cast<real_prec>(nmean2));///c_Ff.shot_noise;
                c_Ff.data_g[i]=(static_cast<real_prec>(dummy[i]) - static_cast<real_prec>(dummy2[i])/static_cast<real_prec>(nmean2));///c_Ff.shot_noise;

	      this->ngal_delta=ngal_new1;
	      c_Ff.n_gal =ngal_new1;


              real_prec factor=pow(static_cast<real_prec>(c_Ff.Lside),1.5)/static_cast<real_prec>(c_Ff.data_g.size());
	      c_Ff.normal_power=pow(factor,-2);
	      c_Ff.shot_noise=0;

	      c_Ff.set_pars(this->Lbox,kmax_y_ds, ndel_data, ndel_window, N_log_bins,N_mu_bins,MAS_correction,mass_assignment_scheme,kmin_bk,kmax_bk);
	      c_Ff.set_bins(type_of_binning);
              if(true==verbose)
                  {
#ifdef _FULL_VERBOSE_
                      So.write_fftw_parameters((void *)&s_p_box);
#endif
                      c_Ff.write_fftw_parameters();
                }
              c_Ff.fftw_vectors(use_random_catalog);

                  }

	  else  if(true==measure_cross) // If no crossed power, read per default thedelta_grid_file  // if we measure the cross, then
	    {
	      c_Ff.data_g.resize(c_Ff.NT,0);
	      fftw_array<float> dummy(c_Ff.NT);
	      switch(measure_cross_from_1){ //choose the file to get the auto power from
	      case(1):
		this->c_Fm.read_array(this->delta_grid_file, dummy,c_Ff.NT);
		break;
	      case(2):
		this->c_Fm.read_array(this->delta_grid_file2, dummy,c_Ff.NT);
		break;
	      case(3):
		this->c_Fm.read_array(this->delta_grid_file3, dummy,c_Ff.NT);
		break;
	      case(4):
		this->c_Fm.read_array(this->delta_grid_file4, dummy,c_Ff.NT);
		break;
	      }

	      if (input_type=="delta_grid")
		{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
                  for(ULONG i=0;i<c_Ff.NT;++i)c_Ff.data_g[i]=static_cast<real_prec>(dummy[i]);
		  ngal_new=this->ngal_delta;
		}
	      else if (input_type=="density_grid")
		{
		  ngal_new=0;
                  for(ULONG i=0;i<c_Ff.NT;++i)ngal_new+=static_cast<real_prec>(dummy[i]);
                  real_prec nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(c_Ff.NT);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
		  for(ULONG i=0;i<c_Ff.NT;++i)
                    c_Ff.data_g[i]=(static_cast<real_prec>(dummy[i])/static_cast<real_prec>(nmean))-1.0;
		}

              c_Ff.shot_noise=0;
              if(true==this->SN_correction)
                c_Ff.shot_noise=static_cast<real_prec>(pow(c_Ff.Lside,3))/static_cast<real_prec>(ngal_new);
              So.message_screen("Shot Noise =",c_Ff.shot_noise);

	      switch(measure_cross_from_2)
		{ //choose the file to get the auto power from
		case(1):
		  this->c_Fm.read_array(this->delta_grid_file, dummy,c_Ff.NT);
		  break;
		case(2):
		  this->c_Fm.read_array(this->delta_grid_file2, dummy,c_Ff.NT);
		  break;
		case(3):
		  this->c_Fm.read_array(this->delta_grid_file3, dummy, c_Ff.NT);
		  break;
		case(4):
		  this->c_Fm.read_array(this->delta_grid_file4, dummy, c_Ff.NT);
		  break;
		}

	      c_Ff.data_gp.resize(c_Ff.NT,0);

	      if (input_type=="delta_grid")
		{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
		  for(ULONG i=0;i<c_Ff.NT;++i)
                    c_Ff.data_gp[i]=static_cast<real_prec>(dummy[i]);


		  ngal_new=this->ngal_delta;
		}
	      else if (input_type=="density_grid")
		{
		  ngal_new=0;




//                  for(ULONG i=0;i<c_Ff.NT;++i)c_Ff.data_gp[i]=static_cast<real_prec>(dummy[i]);

#ifdef _KONV_
		  ULONG ntt=Nft*Nft*(Nft/2+1);
		  FileOutput File;
                  vector<real_prec>kernel(ntt,0);
		  File.read_array("/net/vaina/scratch/balaguera/data/Numerics/IACmocks/ANALYSIS/BAM/Output_Minerva_R220_III/Bam_Kernel.dat",kernel);
		  convolvek(Nft,c_Ff.data_gp, kernel,c_Ff.data_gp);
#endif

//                 for(ULONG i=0;i<c_Ff.NT;++i)
//                    ngal_new+=c_Ff.data_gp[i];

//                  real_prec nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(c_Ff.NT);
//#pragma omp parallel for
//                  for(ULONG i=0;i<c_Ff.NT;++i)
//                    c_Ff.data_gp[i]=(static_cast<real_prec>(c_Ff.data_gp[i])/static_cast<real_prec>(nmean))-1.0;

                this->So.message_warning("In PowerSpectrumF.cpp I have manually set that the field2 is already a delta. This can be solved if I define a new variable input_type for the second field in the cross power");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
                for(ULONG i=0;i<c_Ff.NT;++i)
                    c_Ff.data_gp[i]=static_cast<real_prec>(dummy[i]);


              }


	      // Here we find 1+delta in the grid
              real_prec factor=pow(static_cast<real_prec>(c_Ff.Lside),1.5)/static_cast<real_prec>(c_Ff.data_g.size());
	      c_Ff.normal_power=pow(factor, -2);

	      //IF IT IS DM
              c_Ff.shot_noise2=0;//static_cast<real_prec>(pow(c_Ff.Lside,3))/static_cast<real_prec>(ngal_new);

	      if(true==this->SN_correction)
		{
		  So.message_screen("Shot noise 1 =", c_Ff.shot_noise);
		  So.message_screen("Shot noise 2 =", c_Ff.shot_noise2);
		}
	      else
		So.message_screen("No SN correction");


	      c_Ff.set_pars(this->Lbox,kmax_y_ds, ndel_data, ndel_window, N_log_bins,N_mu_bins,MAS_correction,mass_assignment_scheme,kmin_bk,kmax_bk);
	      c_Ff.set_bins(type_of_binning);
              if(true==verbose)
                  {
#ifdef _FULL_VERBOSE_
	      So.write_fftw_parameters((void *)&s_p_box);
#endif
              c_Ff.write_fftw_parameters();
                }
              c_Ff.fftw_vectors(use_random_catalog);
	    }
	}

      // *****************************************************************************************
      // *****************************************************************************************
      // WELCOME TO FOURIER SPACE
      // *****************************************************************************************
      if(statistics=="Pk_fkp" || statistics=="Pk_ys" || statistics=="Pk_yb"  || statistics=="Pk_ybc"   || statistics=="Pk_ysc" ||  statistics=="Pk_y_ds")
	{

	  // *****************************************************************************************
	  // FFTW and estimates of power spectrum                                                    *
	  // *****************************************************************************************




          kvector_data.clear();
          kvector_data.shrink_to_fit();
          kvector_window.clear();
          kvector_window.shrink_to_fit();
          if("linear" == type_of_binning)
            {
              for(int i=0;i<c_Ff.Nnp_data;i++)
                kvector_data.push_back(c_Ff.DeltaK_data*(i+k_bin_step));

              for(int i=0;i<c_Ff.Nnp_window;i++)
                kvector_window.push_back(c_Ff.DeltaK_window*(i+k_bin_step));
            }

          else if("log"==type_of_binning)
            {
               for(int i=0;i<kvector_data.size();i++)
                    kvector_data.push_back(c_Ff.kmin*pow(10,(i-0.5)*c_Ff.Deltal));
               for(int i=0;i<kvector_window.size();i++)
                    kvector_window.push_back(c_Ff.kmin*pow(10,(i-0.5)*c_Ff.Deltal));
             }

          for(int i=0;i<c_Ff.Nnp_data;i++)
              kvector_data2d.push_back(c_Ff.DeltaK_data*(i+k_bin_step));
          for(int i=0;i<N_mu_bins;i++)
              muvector.push_back(-1.0+c_Ff.Deltamu*(i+0.5));


	  // *****************************************************************************
	  // Resize arrays for P(k), and 2d P(k). Compute and write to file

          this->pk0.clear();
          this->pk0.shrink_to_fit();
	  this->pk0.resize(c_Ff.Nnp_data,0); //Monopole
	  //#ifdef _WRITE_MULTIPOLES_
          this->pk2.resize(c_Ff.Nnp_data,0); //Quadrupole
	  this->pk4.resize(c_Ff.Nnp_data,0); //Hexadecapole
	  //#endif
          this->pk_w.clear();
          this->pk_w.shrink_to_fit();
          this->pk_w.resize(c_Ff.Nnp_window,0); //W(k)

          this->modes_g.clear();
          this->modes_g.shrink_to_fit();
          this->modes_g.resize(c_Ff.Nnp_data,0); //Needed in case we use the Veff for the variance

	  //#ifdef _WRITE_2DPOWER_
          this->pkk.resize(c_Ff.Nnp_data);
	  this->pmk.resize(N_mu_bins);
	  for(int i=0;i<c_Ff.Nnp_data;i++)this->pkk[i].resize(c_Ff.Nnp_data,0);
	  for(int i=0;i<N_mu_bins;i++)this->pmk[i].resize(c_Ff.Nnp_data,0);
	  //#endif

          this->sigma_fkp.clear();
          this->sigma_fkp.shrink_to_fit();
          this->sigma_fkp.resize(c_Ff.Nnp_data,0);

	  // ****************************************************************************
	  // Get power spectrum and more

	  if(statistics=="Pk_fkp")
	    {
              c_Ff.get_power_spectrum_fkp(&s_p_box, this->pk0,this->pk2,this->pk4,this->pk_w,this->pkk,this->pmk,this->modes_g);

              sigma_y_l2.resize(c_Ff.Nnp_data,0);
	      sigma_y_l4.resize(c_Ff.Nnp_data,0);
              if(true==FKP_error_bars)
		{
		  So.message("Computing FKP error bars");
		  c_Ff.get_fkp_error_bars(&s_p_box,&s_data_struct_r, kvector_data, this->pk0, this->modes_g, this->sigma_fkp);
		}
	    }
	  else if(statistics=="Pk_yb" || statistics=="Pk_ybc" || statistics=="Pk_ys" || statistics=="Pk_y_ds" || statistics=="Pk_ysc" )
	      c_Ff.get_power_spectrum_yamamoto(&s_p_box, this->pk0,this->pk2,this->pk4,this->modes_g);

	  //MISSINGN ERROR BARS FROM YAMAMOTO HERE.

	}

      // Estimates of Bispectrum. Using the DFT already done for P(k)
      else if(statistics=="Bk_fkp")
	{

	  if(type_of_binning=="linear")
	    for(int i=0;i<c_Ff.Nnp_data;i++)
	      kvector_data_b.push_back(c_Ff.DeltaK_data*(i+0.5)); //Oficcial binning
	  else
	    if(type_of_binning=="log"){
	      for(int i=0;i<c_Ff.Nnp_data;i++)
		kvector_data_b.push_back(c_Ff.kmin*pow(10,(i-0.5)*c_Ff.Deltal));
	    }
	  bispectrum.resize(Nft*Nft*Nft);
	  sn_bispectrum.resize(Nft*Nft*Nft);
	  modes_tri.resize(Nft*Nft*Nft);

	  c_Ff.get_bispectrum_fkp('d', &s_p_box, bispectrum, sn_bispectrum, modes_tri);
	  c_Fm.write_to_file(file_bispectrum,kvector_data_b,bispectrum,modes_tri);
	}

      // Estimates of Bispectrum for FKP using fast version
      else if(statistics=="Bk_fkp_fast")
	{
	  //for(int i=0;i<c_Ff.Nshells_bk;i++)kvector_data_b.push_back(c_Ff.DeltaK_data*(i+0.5)); //Oficcial binning

	  this->pk0.resize(c_Ff.Nnp_data,0);

	  for(int i=0;i<c_Ff.Nshells_bk;i++)
	    kvector_data_b.push_back(c_Ff.DeltaK_data*(i+1)); //Jennifer's binning

	  bispectrum.resize(c_Ff.Nshells_bk*c_Ff.Nshells_bk*c_Ff.Nshells_bk,0);
	  sn_bispectrum.resize(c_Ff.Nshells_bk*c_Ff.Nshells_bk*c_Ff.Nshells_bk,0);
	  modes_tri.resize(c_Ff.Nshells_bk*c_Ff.Nshells_bk*c_Ff.Nshells_bk,0);

	  c_Ff.get_power_spectrum_for_bispectrum(&s_p_box, this->pk0);
	  c_Ff.get_bispectrum_fkp_fast(&s_p_box,this->pk0,bispectrum,modes_tri,file_bispectrum);

	}


      // *********************************************************************
      // Write Log file.
      if(true==verbose)
          c_Ff.write_fftw_parameters((void *)&s_p_box,file_power_log);
      // *********************************************************************q

#ifndef _WRITE_MULTIPOLES_
      write_power_and_modes();
#else
      write_power_spectrum();
#endif


#if defined(_USE_MASS_CUTS_PK_) || defined (_USE_MASS_BINS_PK_)
    }
#endif

}





// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// Used specifically in case we supply grids 1+delta1 and 1+delta2 and we want to compute the auto and cross power
// threof
void PowerSpectrumF::compute_power_spectrum_grid(const vector<real_prec> &data_in)
{
    this->So.enter(__PRETTY_FUNCTION__);

#ifdef _USE_OMP_
   int NTHREADS=_NTHREADS_;
   omp_set_num_threads(NTHREADS);
#endif

   FftwFunctions c_Ff(this->statistics,this->Nft, false);
  s_parameters_box s_p_box;

  s_p_box.mas=this->mass_assignment_scheme;
  s_p_box.ave=type_of_binning;
  s_p_box.k_bin_step=k_bin_step;
  s_p_box.use_MAS_correction=this->MAS_correction;
  s_p_box.FKP_weight=this->FKP_weight;
  s_p_box.FKP_error_bars=FKP_error_bars;
  s_p_box.FKP_error_bars_exact= FKP_error_bars_exact;
  s_p_box.use_SN_correction=this->SN_correction;
  s_p_box.Pest=this->Pest;
  c_Ff.Lside=this->Lbox;
  c_Ff.NT = data_in.size();
  c_Ff.set_pars(this->Lbox,kmax_y_ds, ndel_data, ndel_window, N_log_bins,N_mu_bins,MAS_correction,mass_assignment_scheme,kmin_bk,kmax_bk);
  c_Ff.set_bins(type_of_binning);
  c_Ff.write_fftw_parameters();
  c_Ff.fftw_vectors(use_random_catalog);

  real_prec ngal_new,nmean;
  real_prec vol = static_cast<real_prec>(pow(this->Lbox,3));

  // if delta  = rho - mean, (with Ngal given from par file) use normal=Ngal²/Vol (i.e, constructing delta from a catalog)
  // if delta = rho/mean -1, use normal = Ncells²/Vol (i.e, givng delta from outside.

  if(this->input_type=="density_grid")
    {
      ngal_new=get_nobjects(data_in);
      this->ngal_delta=ngal_new;
      c_Ff.n_gal= ngal_new;
      nmean=static_cast<double>(ngal_new)/static_cast<double>(c_Ff.NT);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<c_Ff.NT;++i)
        c_Ff.data_g[i]=static_cast<real_prec>(data_in[i])-static_cast<real_prec>(nmean);
       c_Ff.normal_power=static_cast<double>((static_cast<double>(ngal_new)/static_cast<double>(vol))*static_cast<double>(ngal_new));
    }
  else if(this->input_type=="delta_grid")
      {
	      ngal_new=1.0;
        this->ngal_delta=ngal_new;
        c_Ff.n_gal= ngal_new;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      	for(ULONG i=0;i<c_Ff.NT;++i)
	        c_Ff.data_g[i]=data_in[i];
        c_Ff.normal_power=static_cast<real_prec>(c_Ff.NT)*static_cast<real_prec>(c_Ff.NT)/static_cast<real_prec>(vol);
      }

#ifdef _FULL_VERBOSE_
  if(ngal_new>1)
    So.message_screen("Number of tracers in input density field = ",ngal_new);
  So.message_screen("Normalization of power = ",c_Ff.normal_power);
#endif
  c_Ff.shot_noise=vol/static_cast<real_prec>(ngal_new);
  this->shot_noise=c_Ff.shot_noise;
#ifdef _FULL_VERBOSE_
  if(true==this->SN_correction)
    So.message_screen("Poisson Shot noise = ",c_Ff.shot_noise);
  else
    So.message_screen("Poisson Shot noise (not applied) = ",c_Ff.shot_noise);
#endif

  kvector_data.resize(c_Ff.Nnp_data, 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<c_Ff.Nnp_data;i++)
    this->kvector_data[i]=(i+k_bin_step)*c_Ff.DeltaK_data;

  this->pk0.resize(c_Ff.Nnp_data,0); //Monopole
  this->modes_g.resize(c_Ff.Nnp_data,0); //Needed in case we use the Veff for the variance
  c_Ff.get_power_spectrum_fkp(&s_p_box, this->pk0,this->modes_g);
  //  c_Ff.free_fftw_vectors();

}




// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
void PowerSpectrumF::compute_marked_power_spectrum_grid(const vector<real_prec> &data_in,const vector<real_prec> &data_in_MW)
{


#ifdef _USE_OMP_
    int NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif

  FftwFunctions c_Ff(this->statistics,this->Nft, false);
  s_parameters_box s_p_box;

  s_p_box.mas=this->mass_assignment_scheme;
  s_p_box.ave=type_of_binning;
  s_p_box.k_bin_step=k_bin_step;
  s_p_box.use_MAS_correction=MAS_correction;
  s_p_box.FKP_weight=FKP_weight;
  s_p_box.FKP_error_bars=FKP_error_bars;
  s_p_box.FKP_error_bars_exact= FKP_error_bars_exact;
  s_p_box.use_SN_correction=SN_correction;
  s_p_box.Pest=Pest;
  c_Ff.Lside=this->Lbox;
  c_Ff.NT = data_in.size();
  c_Ff.set_pars(this->Lbox,kmax_y_ds, ndel_data, ndel_window, N_log_bins,N_mu_bins,MAS_correction,mass_assignment_scheme,kmin_bk,kmax_bk);
  c_Ff.set_bins(type_of_binning);

  c_Ff.fftw_vectors(use_random_catalog);

  real_prec ngal_new=get_nobjects(data_in);

  real_prec nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(c_Ff.NT);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<c_Ff.NT;++i)
    c_Ff.data_g[i]=(static_cast<real_prec>(data_in_MW[i])-static_cast<real_prec>(data_in[i]))/static_cast<real_prec>(nmean);

  this->ngal_delta=ngal_new;
  c_Ff.n_gal= ngal_new;

  real_prec factor=pow(static_cast<real_prec>(c_Ff.Lside),1.5)/static_cast<real_prec>(c_Ff.data_g.size());
  c_Ff.normal_power=pow(factor,-2);
  c_Ff.shot_noise=(this->var_prop-1.0)*static_cast<real_prec>(pow(c_Ff.Lside,3))/static_cast<real_prec>(ngal_new);

  kvector_data.resize(c_Ff.Nnp_data, 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<c_Ff.Nnp_data;i++)
    kvector_data[i]=c_Ff.DeltaK_data*(i+k_bin_step);

  this->pk0.resize(c_Ff.Nnp_data,0); //Monopole
  this->modes_g.resize(c_Ff.Nnp_data,0); //Needed in case we use the Veff for the variance

  c_Ff.get_power_spectrum_fkp(&s_p_box, this->pk0,this->modes_g);
}


// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
void PowerSpectrumF::compute_power_spectrum_grid()
{

#ifdef _USE_OMP_
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

  int ir=1;
  FftwFunctions c_Ff(this->statistics,this->Nft, this->measure_cross);

  c_Ff.NT=this->Nft*this->Nft*this->Nft;
  vector<real_prec> data_in(c_Ff.NT,0);
  //      fftw_array<re> dummy(c_Ff.NT);
  this->c_Fm.read_array_t<real_prec>(this->delta_grid_file, data_in);


      //      this->c_Fm.read_array(this->delta_grid_file, dummy,c_Ff.NT);
// #ifdef _USE_OMP_
// #pragma omp parallel for
// #endif
//   for(ULONG i=0;i<c_Ff.NT;++i)
//     data_in[i]=  (data_in[i]>1e3 ? data_in[i]: 0);// static_cast<real_prec>(dummy[i]);


  s_parameters_box s_p_box;

  s_p_box.mas=this->mass_assignment_scheme;
  s_p_box.ave=type_of_binning;
  s_p_box.k_bin_step=k_bin_step;
  s_p_box.use_MAS_correction=MAS_correction;
  s_p_box.FKP_weight=FKP_weight;
  s_p_box.FKP_error_bars=FKP_error_bars;
  s_p_box.FKP_error_bars_exact= FKP_error_bars_exact;
  s_p_box.use_SN_correction=SN_correction;
  s_p_box.Pest=Pest;
  c_Ff.Lside=this->Lbox;
  c_Ff.NT = data_in.size();
  c_Ff.set_pars(this->Lbox,kmax_y_ds, ndel_data, ndel_window, N_log_bins,N_mu_bins,MAS_correction,mass_assignment_scheme,kmin_bk,kmax_bk);
  c_Ff.set_bins(type_of_binning);
  c_Ff.fftw_vectors(use_random_catalog);

#ifdef _FULL_VERBOSE_
  c_Ff.write_fftw_parameters();
#endif

  real_prec ngal_new=get_nobjects(data_in);
#ifdef _FULL_VERBOSE_
  So.message_screen("Number of objects =",ngal_new);
#endif

  real_prec nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(c_Ff.NT);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<c_Ff.NT;++i)
    c_Ff.data_g[i]=(static_cast<real_prec>(data_in[i])/static_cast<real_prec>(nmean))-1.;
  this->ngal_delta=ngal_new;
  c_Ff.n_gal =ngal_new;


  real_prec factor=pow(static_cast<real_prec>(c_Ff.Lside),1.5)/static_cast<real_prec>(c_Ff.data_g.size());
  c_Ff.normal_power=pow(factor,-2);

  c_Ff.shot_noise=0;
  if(true==this->params._SN_correction())
    c_Ff.shot_noise=static_cast<real_prec>(pow(c_Ff.Lside,3))/static_cast<real_prec>(ngal_new);
#ifdef _FULL_VERBOSE_
  So.message_screen("Shot Noise =",c_Ff.shot_noise);
#endif

  kvector_data.resize(c_Ff.Nnp_data, 0);

 if("linear" == type_of_binning)
   {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(int i=0;i<c_Ff.Nnp_data;i++)
       kvector_data[i]=c_Ff.DeltaK_data*(i+k_bin_step);

   }
 else if("log"==type_of_binning)
   {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(int i=0;i<c_Ff.Nnp_data;i++)
       kvector_data[i]=c_Ff.kmin*pow(10,(i-0.5)*c_Ff.Deltal);
   }


  this->pk0.resize(c_Ff.Nnp_data,0); //Monopole
  this->modes_g.resize(c_Ff.Nnp_data,0); //Needed in case we use the Veff for the variance
  c_Ff.get_power_spectrum_fkp(&s_p_box, this->pk0,this->modes_g);

  write_power_and_modes();

}

// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
