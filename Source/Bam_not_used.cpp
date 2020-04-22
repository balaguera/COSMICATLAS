// Not used functions in BAM

void Bam::move_particles_lpt()
{
  
  cout<<endl;
  this->So.message_screen("Moving tracers to fix small-scale clustering");
  cout<<endl;
  //i. Get overdensity field of tracers. The vector delta_Y_new is istill active under this preprocdir
  this->delta_Y_new.resize(this->NGRID,0);
  string pname= "_realization1";
  string fname=this->Output_directory+"MOCK_TR"+pname+"_"+"MASY"+to_string(this->iMAS_Y)+"_Nft"+to_string(this->Nft)+"_z"+dto_string(this->redshift);
  this->File.read_array_t<PrecType_Y>(fname+".dat", delta_Y_new);
  get_overdens(this->delta_Y_new, this->delta_Y_new);
  vector<real_prec> vaux1(this->NGRID,0);
  vector<real_prec> vaux2(this->NGRID,0);
  // Convolve the term factor * ð with a Gaussian kernel with width = cell size
  // The factor is introduced after the convolution
  real_prec Rscale=this->params._d1(); // Smoothing using the cell size
  kernelcomp(Lbox,Lbox,Lbox,this->params._d1(),this->params._d2(), this->params._d3(),Nft,Nft,Nft,Rscale,1,this->Output_directory);
  convcomp(Lbox,Lbox,Lbox,this->params._d1(),this->params._d2(),this->params._d3(), Nft, Nft,Nft, delta_Y_new, vaux1, 1,Rscale,this->Output_directory);
  this->So.DONE();
  real_prec factor=-0.005;
#pragma omp parallel for
  for(ULONG i=0;i<this->NGRID;i++)
    {
      real_prec psilin=-factor*vaux1[i];
      //      real_prec psilin=-factor*delta_Y_new[i];
      real_prec insidesqrt=static_cast<real_prec>(1.-2./3.*psilin);
      real_prec psisc=3.;
      
      if (insidesqrt > 0.)
	psisc=-static_cast<real_prec>(3.*(sqrt(insidesqrt)-1.));
      vaux1[i]=psisc;
    }
#define Psi delta_Y_new
  Rscale=20.0;
  kernelcomp(Lbox,Lbox,Lbox,this->params._d1(),this->params._d2(), this->params._d3(),Nft,Nft,Nft,Rscale,1,this->Output_directory);
  convcomp(Lbox,Lbox,Lbox,this->params._d1(),this->params._d2(),this->params._d3(), Nft, Nft,Nft, delta_Y_new, Psi, 1,Rscale,this->Output_directory);
  convcomp(Lbox,Lbox,Lbox,this->params._d1(),this->params._d2(),this->params._d3(), Nft, Nft,Nft, vaux1, vaux2, 1,Rscale,this->Output_directory);
#pragma omp parallel for
  for(ULONG i=0;i<this->NGRID;i++)
    Psi[i]=(-factor*Psi[i])+(vaux1[i]-vaux2[i]);   // factor *  delta\otimes K +  Psi^SC- K \otimes Psi^SC
  //get x-comp
  this->patchy.sfmodel=1;
  this->So.message_screen("Getting Vel_x");
  //this->patchy.Lag2Eul_comp(1.0, kth,1,true,vaux1,1);
  this->patchy.theta2velcomp(Psi,vaux1,false,false,1);
  this->So.DONE();
  //get y-comp
  this->So.message_screen("Getting Vel_y");
  //this->patchy.Lag2Eul_compB(kth,1,true,vaux2,2);
  this->patchy.theta2velcomp(Psi,vaux2,false,false,2);
  this->So.DONE();
  
  //get z-comp
  this->So.message_screen("Getting Vel_z");
  vector<real_prec> vaux3(this->NGRID,0);
  //this->patchy.Lag2Eul_compB(kth,1,true,vaux3,3);
  this->patchy.theta2velcomp(Psi,vaux3,false,false,3);
  this->So.DONE();
  
#undef Psi
  s_params_box_mas box;
  box.min1=this->params._xllc();
  box.min2=this->params._yllc();
  box.min3=this->params._zllc();
  box.d1=this->params._d1();
  box.d2=this->params._d2();
  box.d3=this->params._d3();
  box.Lbox=this->params._Lbox();
  box.Nft=this->params._Nft();
  box.NGRID=this->NGRID;
  
  // Displace particles
  vector<real_prec> prop;
  this->patchy.fnameTRACERCAT=this->Output_directory+"CAT_realization"+to_string(1)+"_"+this->patchy.stradd+string(".txt");
  
  ULONG Ntracers=File.read_file(this->patchy.fnameTRACERCAT,prop, omp_get_max_threads());
  int NCOLS=(static_cast<ULONG>(prop.size()/Ntracers));
  
  this->So.message_screen("Moving tracers");
  
  
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0; i< Ntracers; ++i)
    {
      real_prec x=prop[0+i*NCOLS];
      real_prec y=prop[1+i*NCOLS];
      real_prec z=prop[2+i*NCOLS];
      ULONG index=grid_ID(&box, x,y,z);
      prop[0+i*NCOLS]+=vaux1[index];
      prop[1+i*NCOLS]+=vaux2[index];
      prop[2+i*NCOLS]+=vaux3[index];
    }
  this->So.DONE();
  
  this->So.message_screen("Writing new catalog");
  ofstream trout;
  string newcat=this->patchy.fnameTRACERCAT+"NEW";
  trout.open(newcat);
  for(ULONG i=0; i<Ntracers;++i)
    {
      for(int j=0; j<NCOLS; ++j)
	trout<<prop[j+i*NCOLS]<<"\t";
      trout<<endl;
    }
  trout.close();
  this->So.DONE();
  
  this->patchy.fnameTRACERCAT=newcat;
  prop.clear();
  prop.shrink_to_fit();
  vaux1.clear();
  vaux1.shrink_to_fit();
  vaux2.clear();
  vaux2.shrink_to_fit();
  vaux3.clear();
  vaux3.shrink_to_fit();
  
}

// *****************************************************************************************************
// *****************************************************************************************************
// *****************************************************************************************************
// *****************************************************************************************************
// ********************************************************************************************************************************************************************************
// ********************************************************************************************************************************************************************************
// ********************************************************************************************************************************************************************************
// ********************************************************************************************************************************************************************************

void Bam::assign_tracer_mass()
{

  //#define verify_cells

#define _consistency_test_
  
  
  this->So.message_screen("**Assigning mass to particles**");
  cout<<endl;

  real_prec mass_min=pow(10,params._LOGMASSmin())*params._MASS_units();
  real_prec mass_max=pow(10,params._LOGMASSmax())*params._MASS_units();

  gsl_rng_env_setup();
  gsl_rng_default_seed=71232;
  const gsl_rng_type *  T= gsl_rng_ranlux;
  gsl_rng * r = gsl_rng_alloc (T);

  this->tracer.set_params_catalog(this->params); // this was only loaded under calibration and read_ref_cat
  



  s_params_box_mas box;
  box.min1=this->params._xllc();
  box.min2=this->params._yllc();
  box.min3=this->params._zllc();
  box.Lbox=this->params._Lbox();
  box.Nft=this->Nft;
  box.d1= box.Lbox/static_cast<real_prec>(box.Nft);		/* grid spacing x-direction */
  box.d2= box.d1;
  box.d3= box.d1;
  box.NGRID=this->NGRID;
  
  
  
  // Ojo, el NGRID aca puede ser el los, el que se usó ar la calibración
  // Load mass field of the current mock

#ifdef _ONLY_POST_PROC_

#ifndef _consistency_test_
  string fname=this->Output_directory+"MOCK_TR_MASS_realization1_MASY"+to_string(this->iMAS_Y)+"_Nft"+to_string(this->Nft)+"_z"+dto_string(this->redshift);
#else
  string fname=this->Output_directory+"TR_MASS_DENS_FIELD";
#endif
  this->patchy.fname_MOCK_MASS=fname;
  
#ifndef _consistency_test_  
  fname=this->Output_directory+"MOCK_TR_realization1_MASY"+to_string(this->iMAS_Y)+"_Nft"+to_string(this->Nft)+"_z"+dto_string(this->redshift);
#else
  fname=this->Output_directory+"TR_DENS_FIELD";
#endif
  this->patchy.fname_MOCK_NCOUNTS=fname;
#endif

  vector<real_prec> mass_field(this->NGRID,0);
#ifndef _consistency_test_
  this->File.read_array(this->patchy.fname_MOCK_MASS+".dat", mass_field);  // COMMEntED FOR TESTING. ENABLE WHEN READY
#endif
  

  vector<real_prec> v_ncounts(this->NGRID,0);
#ifndef _consistency_test_
  this->File.read_array(this->patchy.fname_MOCK_NCOUNTS+".dat", v_ncounts); // COMMEntED FOR TESTING. ENABLE WHEN READY
#endif
  
  // string fname2=this->Output_directory+"MOCK_TR_MASS_iteration30_MASY"+to_string(this->iMAS_Y)+"_Nft"+to_string(this->Nft)+"_z"+dto_string(this->redshift);
  // vector<real_prec> mass_field2(this->NGRID,0);
  // this->File.read_array(fname2+".dat", mass_field2);

  // fname=this->Output_directory+"MOCK_TR_iteration30_MASY"+to_string(this->iMAS_Y)+"_Nft"+to_string(this->Nft)+"_z"+dto_string(this->redshift);
  // this->patchy.fname_MOCK_NCOUNTS=fname;
  // vector<real_prec> v_ncounts2(this->NGRID,0);
  // this->File.read_array(fname+".dat", v_ncounts2);


  //ofstream tea;tea.open(this->Output_directory+"test.txt");
  // for (ULONG id=0;id<this->NGRID; ++id)
  //   if(v_ncounts[id]>0)
  //     {
  // 	if( static_cast<int>(floor(mass_field[id]/mass_min))<v_ncounts[id])
  // 	  cout<<mass_field[id]<<"  "<<v_ncounts[id]<<endl;
  //     }
  //     //tea<<floor(mass_field[id]/mass_min)<<"  "<<v_ncounts[id]<<endl;
  // // tea.close();
 



  vector<real_prec> prop;
#ifndef _consistency_test_  
  this->patchy.fnameTRACERCAT=this->Output_directory+"CAT_realization"+to_string(1)+"_"+this->patchy.stradd+string(".txt");
#else
  this->patchy.fnameTRACERCAT=this->Output_directory+"newcat.txt";
#endif
  //this->patchy.fnameTRACERCAT=this->params._file_catalogue();
  ULONG Ntracers=File.read_file(this->patchy.fnameTRACERCAT,prop, omp_get_max_threads());
  int NCOLS=static_cast<int>(prop.size()/Ntracers);
  
  //#define method1

  //#define method2
#ifdef method2
#define _use_mass_function_
#endif
  
  //#define method3
#ifdef method3
    //#define _use_mass_function_
#endif


#define method4
#ifdef method4
  //#define _use_mass_function_
#endif

  
#ifdef verify_cells    
  vector<int> counter_a(this->NGRID,0);
#endif
  
  vector<ULONG> index(Ntracers,0);
  


  // ************************** test with reference to check wich method works to give masses**********************  

  // ofstream newc; newc.open("newcat.txt");
  // ULONG nnn=0;
  // for(ULONG i=0; i< Ntracers; ++i)
  //   {
  //     if(prop[7+i*NCOLS]>=mass_min)
  // 	{
  // 	  real_prec x=prop[1+i*NCOLS];  
  // 	  real_prec y=prop[2+i*NCOLS];
  // 	  real_prec z=prop[3+i*NCOLS];
  // 	  //ULONG id=grid_ID(&box, x,y,z);
  // 	  //index[nnn]=id;
  // 	  //nnn++;
  // 	  newc<<x<<"\t"<<y<<"\t"<<z<<"  "<<prop[7+i*NCOLS]<<endl;
  // 	}
  //   }
  // newc.close();
  // exit(0);
  // ********************************************************************

  struct s_rinfo{
    vector<real_prec>mass_tracers;
  };

  vector<s_rinfo> mass_info(this->NGRID);

#pragma omp parallel for 
  for(ULONG i=0; i< Ntracers; ++i)
    {
      real_prec x=prop[0+i*NCOLS];
      real_prec y=prop[1+i*NCOLS];
      real_prec z=prop[2+i*NCOLS];
#ifdef _consistency_test_
      real_prec mass=prop[3+i*NCOLS];
#endif
      ULONG id=grid_ID(&box, x,y,z);
      index[i]=id;
#ifdef verify_cells    
#pragma omp atomic update
      counter_a[id]++;
#endif
#ifdef _consistency_test_
      v_ncounts[id]++;  // just for the test
      mass_field[id]+=mass;
      mass_info[id].mass_tracers.push_back(mass);


#endif
    }



  
//   this->tracer.Halo.resize(this->NGRID);
//   for(ULONG i=0;i<this->NGRID;++i)
//     this->tracer.Halo[i].mass=mass_field[i];

// #ifdef _consistency_test_
//   string auxf = this->Output_directory+"tracer_ref_Cell_mass_function.txt";
// #else
//   string auxf = this->Output_directory+"tracer_mock_Cell_mass_function.txt";
// #endif
//   this->tracer.NOBJS=NGRID;
//   this->tracer.get_mass_function(auxf);




  
  
  
#ifdef verify_cells    
  /*  // Hemos verificado que la masa es cero cuando no hay objetos.
      for(ULONG id =0 ; id<this->NGRID;++id)
      if(v_ncounts[id]!=0 && mass_field[id]>0)
      cout<<v_ncounts[id]<<"  "<<mass_field[id]<<endl;
  */
  
  vector<bool> off_cells(this->NGRID, false); // This is jut becuase the particles from patchy do not generat the same number counts as the number counts they are built from
  ULONG nfc=0;
  ULONG npfc=0;
  real_prec mfc=0;
  
  // Loop dedicated to veryfy that the ncounts of the cat are those of the ncounts generated by BM
  for(ULONG id =0 ; id<this->NGRID;++id)
    if(v_ncounts[id]!=counter_a[id])
      {
	off_cells[id]=true;
        npfc+=v_ncounts[id];
        nfc++;
        mfc+=mass_field[id];
      }
  if(nfc>0)
    {
      this->So.message_screen("Number of mis-identified cells =", nfc);
      this->So.message_screen("Number of mising particles =", npfc);
      this->So.message_screen("Mising mass =", mfc);
    }
#endif
  
  
  
  this->So.message_screen("Adding masses to catalog:");    
  cout<<endl;
  
  ofstream trout;
  string newcat=this->patchy.fnameTRACERCAT+"MASS";
  this->patchy.fnameTRACERCAT=newcat;
  trout.open(newcat.c_str());
  trout.precision(8);
  trout.setf(ios::showpoint);
  trout.setf(ios::scientific);
  
  
  // ------------------------------------------------------------------------------------------------------------
  
#ifdef method2
  
  this->tracer.Halo.resize(Ntracers);
  
#ifdef _use_mass_function_
  this->So.message_screen("Using reference mass function to complement assignment.");
  this->tracer.define_mass_bins();
  string massf_file=this->Output_directory+"tracer_ref_mass_function.txt";
  vector<real_prec>massf;
  int nbins_mf=this->File.read_file(massf_file,massf,1);
  vector<real_prec>mass_bins(nbins_mf,0);
  vector<real_prec>mfunc(nbins_mf,0);

  
  
#pragma omp parallel for
  for(int i=0; i< nbins_mf; ++i)
    {
      mass_bins[i]=massf[0+i*2];
      real_prec deltaM=this->tracer.MBmax[i]-this->tracer.MBmin[i];
      mfunc[i]=massf[1+i*2]*deltaM;// with n(M)*dM I get the right n(M) dividing by the max and sampling below
    }
  massf.clear();massf.shrink_to_fit();
  real_prec max_mf=static_cast<real_prec>(get_max(mfunc));

  vector<real_prec>mcount(nbins_mf,0);
  for(int i=0; i< nbins_mf; ++i)
    for(int j=i; j< nbins_mf; ++j)
      mcount[i]+=mfunc[j];

  ofstream sal;
  this->So.message_screen("Measuring cumulative MF of the reference:");
  sal.open(massf_file+"_cumulative"); 
  for(int i=0; i< nbins_mf; ++i)
    sal<<mass_bins[i]<<"\t"<<mcount[i]<<endl;
  sal.close();

  // Normalize to maximum
  for(int i=0; i< nbins_mf; ++i)
    mfunc[i]/=max_mf;
  


  mass_max=mass_bins[mass_bins.size()-1];
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, mass_bins.size());
  gsl_spline_init (spline, &mass_bins[0], &mfunc[0], mass_bins.size());

#endif
  
  this->So.message_screen("Assigning masses and writting to file ...");
  ULONG IDD=16751686;
  real_prec massid=mass_field[IDD];
  real_prec auxm=0;
  for(int i=0; i<mass_info[IDD].mass_tracers.size();++i)
    {
      auxm+=mass_info[IDD].mass_tracers[i];
      cout<<mass_info[IDD].mass_tracers[i]<<"   "<<auxm<<"  "<<massid<<endl;
  }

  massid=0;


  vector<int> cells_with_one(this->NGRID,-1);
  vector<bool> nomass(Ntracers,false);
  ULONG nomass_p=0;
  for(ULONG id=0; id< this->NGRID; ++id)
    {
      if(v_ncounts[id]==1)
        cells_with_one[id]=1;
      else if(v_ncounts[id]>1)
          cells_with_one[id]=0;
    }


  for(ULONG i=0; i< Ntracers; ++i)
    {
      ULONG id=index[i];
      real_prec mass_tracer=1e10;
      real_prec fraction_mass;
      
      
      if(1==cells_with_one[id])
	mass_tracer=mass_field[id];
      
      else if(0==cells_with_one[id])
	{
	  
	  if(v_ncounts[id]>1)
	    {
	      
#ifdef _use_mass_function_
	      real_prec prob=-10.0;
	      real_prec xrnew=10.0;
	      
	      while(prob<xrnew)
		{
#endif
		  real_prec xr = static_cast<real_prec>(gsl_rng_uniform(r));
		  fraction_mass = xr*(mass_field[id]-mass_min*(v_ncounts[id]-1)); // my last try
		  mass_tracer=mass_min+fraction_mass;
		  		  
#ifdef _use_mass_function_
		  real_prec aux_mass= (mass_tracer < mass_bins[0]? mass_bins[0]: mass_tracer);
		  if(aux_mass<mass_min || aux_mass>mass_max)
		    prob=0.0;
		  else
		    {
		      xrnew = static_cast<real_prec>(gsl_rng_uniform(r));
		      prob  = static_cast<real_prec>(gsl_spline_eval (spline, aux_mass, acc));
		    }
		}
#endif
	      
	      mass_field[id] -= mass_tracer;
	      v_ncounts[id] --;
	      
	    }
	  
	  else if(1==v_ncounts[id])  // siq ueda una, asígna lo que queda
	    {
	      mass_tracer= mass_field[id];   /// OJO ACA, QUE  ESTABA HACIENDO LA CUENTA DE ASIGNAR AL ULTIM MASA, MAL
	      mass_field[id]-=mass_tracer;
	      v_ncounts[id]=0; //once assigned, subtract
	    }
	  
	  
	  // if(id==IDD)
	  //   {
	  //     massid+=mass_tracer;
	  //     cout<<BLUE<<mass_tracer<<"  "<<massid<<RESET<<"   "<<mass_field[id]<<"   "<<v_ncounts[id]<<endl;
	  //   }
	}
      
      if(mass_tracer<mass_min)
	{
	  nomass_p++;
	  nomass[i]=true;
	}
      else
	{
	  this->tracer.Halo[i].mass=mass_tracer;
	  this->tracer.Halo[i].GridID = id;
	  
	}
    }
  this->So.DONE();
  

  this->So.message_screen("Particles without mass = ", nomass_p);
  this->So.message_screen("Assigning mass using mass fuction");
  real_prec frac=static_cast<real_prec>(nomass_p)/static_cast<real_prec>(Ntracers);
  
  for(ULONG i=0;i<Ntracers;++i)
    {
      if(nomass[i]==true)
	{
	  real_prec prob=-10.0;
	  real_prec xrnew=10.0;
	  real_prec mass_tracer;
	  while(prob<xrnew)
	    {
	      real_prec xr = static_cast<real_prec>(gsl_rng_uniform(r));
	      real_prec fraction_mass = xr*log10(mass_max/mass_min);
	      mass_tracer=pow(10,log10(mass_min)+fraction_mass);
	      
	      real_prec aux_mass= (mass_tracer < mass_bins[0]? mass_bins[0]: mass_tracer);
	      if(aux_mass<mass_min || aux_mass>mass_max)
		prob=0.0;
	      else
		{
		  xrnew = static_cast<real_prec>(gsl_rng_uniform(r));
		  prob  = static_cast<real_prec>(gsl_spline_eval (spline, aux_mass, acc))*frac;
		}
	    }
	  this->tracer.Halo[i].mass=mass_tracer;
	  this->tracer.Halo[i].GridID = index[i];
 
	  
	}

      
      /*    UNCOMMENT THIS WHEN TESTS ARE DONE AND WRITE TO FILE PROCEED
	    prop[6+i*NCOLS]=mass_tracer;
	    for(int j = 0; j<NCOLS; ++j)
		trout<<prop[j+i*NCOLS]<<"\t";
		trout<<endl;
      */
      
    }
  this->So.DONE();
  
  prop.clear();
  prop.shrink_to_fit();
  
  
  real_prec total_mass=0;
#pragma omp parallel for reduction(+:total_mass)
  for(ULONG id=0; id < this->NGRID; ++id)
    total_mass+=mass_field[id];  // should be zero
  So.message_screen("Mass left to be assigned =", total_mass);
  cout<<endl;
  
  
  
  trout.close();
  
  
  
  string fname_mass_function_Y = this->Output_directory+"tracer_mock_mass_function.txt";
  this->tracer.NOBJS=Ntracers;
  this->tracer.get_mass_function(fname_mass_function_Y);
  //    }
  mass_field.clear();mass_field.shrink_to_fit();
  index.clear();
  index.shrink_to_fit();
  prop.clear();
  prop.shrink_to_fit();
  
  
  
  
#endif
  
  // ------------------------------------------------------------------------------------------------------------
  
#ifdef method3
    
  this->So.message_screen("Assigning masses within cells");
  
  struct s_rinfo{
    vector<real_prec>mass_tracers;
  };
  
  vector<s_rinfo>  m_info(this->NGRID);
 
  
  vector<int>  nparts_assigned(this->NGRID,0);
  vector<real_prec>mass_left(this->NGRID,0);
  for(ULONG id=0; id< this->NGRID; ++id)
    {
      int nobs=v_ncounts[id];  // number of objects in the cell
      if(nobs>1)
        {
	  m_info[id].mass_tracers.resize(nobs,mass_min);  //assign the minimum mass to all particles in the cell, to start
	  mass_left[id]=mass_field[id]-static_cast<real_prec>(v_ncounts[id])*mass_min;
	}
      else if(nobs==1) // cells with only one particle
	m_info[id].mass_tracers.push_back(mass_field[id]);
    }
  
  
  So.DONE();
  
  
  // ************************ MASS FUNCTION***************************************
#ifdef _use_mass_function_
  
  this->So.message_screen("Using reference mass function to complement assignment.");
  this->tracer.define_mass_bins();
  string massf_file=this->Output_directory+"tracer_ref_mass_function.txt";
  vector<real_prec>massf;
  int nbins_mf=this->File.read_file(massf_file,massf,1);
  vector<real_prec>mass_bins(nbins_mf,0);
  vector<real_prec>mfunc(nbins_mf,0);

 #pragma omp parallel for
   for(int i=0; i< nbins_mf; ++i)
     {
       mass_bins[i]=massf[0+i*2];
       real_prec deltaM=this->tracer.MBmax[i]-this->tracer.MBmin[i];
       mfunc[i]=(massf[1+i*2])*deltaM;
     }
   massf.clear();massf.shrink_to_fit();

   real_prec max_mf=static_cast<real_prec>(get_max(mfunc));

 #pragma omp parallel for
   for(int i=0; i< nbins_mf; ++i)
     mfunc[i]/=max_mf;

   gsl_interp_accel *acc = gsl_interp_accel_alloc ();
   gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, mass_bins.size());
   gsl_spline_init (spline, &mass_bins[0], &mfunc[0], mass_bins.size());

#endif
   
   // ************************
   
   
   // Assign left mass to particles in cell
   for(ULONG id=0; id< this->NGRID; ++id)
     {

       int nobs=v_ncounts[id];  // number of objects in the cell
       if(nobs>1)
	 {
	   
	   real_prec normalization=mass_left[id]/mass_min; // Mi = Mmin(1+xi), such that sum_i xi = normalization. We randomly add masses in agreement with what we have left.

	   vector<real_prec> mass_aux (nobs,0);
	   
#ifdef _use_mass_function_
	   bool flag=false;
	   
	   real_prec prob=-10.;
	   while(flag==false)
	     {
#endif
	       vector<real_prec> randoms (nobs,0);
	       
               real_prec sran=0;
	       // get the sumn of the random numbers
	       for(int nr=0; nr< nobs; ++nr)
		 {
		   randoms[nr]=static_cast<real_prec>(gsl_rng_uniform(r)); // random between (0,1)
		   sran+=randoms[nr];
		 }
	       
	       // normalize the random numbers such that their sum is equal to normalization. The last one is just normalization - the sum of the normalized
	       for(int nr=0; nr< nobs ; ++nr)
		 {
		   real_prec deltar=randoms[nr]*(normalization/sran);
		   
                   mass_aux[nr]=m_info[id].mass_tracers[nr]*(1.0+deltar);
		   
                   real_prec mass_tr=mass_aux[nr];
#ifdef _use_mass_function_
		   if(mass_tr<mass_bins[0])
                     mass_tr=mass_bins[0];
		   
		   prob = static_cast<real_prec>(gsl_spline_eval (spline, mass_tr, acc));
		   flag=true;
		   if(prob<randoms[nr])
		     flag=false;
                   break;
#endif
                 }
	       
#ifdef _use_mass_function_
	     }
#endif
	   
           for(int nr=0; nr< nobs ; ++nr)
	     m_info[id].mass_tracers[nr]=mass_aux[nr];
	   
	 }
     }
   
   So.DONE();
   
   this->So.message_screen("Writting to file");  
   
   
   this->tracer.Halo.resize(Ntracers);
   
   vector<int>counter(this->NGRID,0);
   
   for(ULONG i=0; i<  Ntracers; ++i)
     {
       ULONG id=index[i];
       real_prec mass_tracer;
       //       cout<<i<<"  "<<id<<"  "<<m_info[id].mass_tracers.size()<<"  "<<counter[id]<<endl;
       mass_tracer=m_info[id].mass_tracers[counter[id]];
       mass_field[id]-=mass_tracer;
       counter[id]++;
       this->tracer.Halo[i].mass=mass_tracer;
       //this->tracer.Halo[i].GridID = id;
       
       // for(ULONG i = 0; i<Ntracers; ++i)
       // 	{
       // 	  for(int j = 0; j<NCOLS; ++j)
       // 	    trout<<prop[j+i*NCOLS]<<"\t";
       // 	  trout<<endl;
       // 	}
     }
   counter.clear();counter.shrink_to_fit();
   So.DONE();
   
   

   real_prec mass_left_p=0;
#pragma omp parallel for reduction(+:mass_left_p)
  for(ULONG id=0; id < this->NGRID; ++id)
    mass_left_p+=mass_field[id];  // should be zero
 
  So.message_screen("Mass left to be assigned", mass_left_p);
  cout<<endl; 
  
  this->So.DONE();
  index.clear();
  index.shrink_to_fit();
  prop.clear();
  prop.shrink_to_fit();

  
  
 //   real_prec total_mass=0;
// #pragma omp parallel for reduction(+:total_mass)
//   for(ULONG id=0; id < this->NGRID; ++id)
//     total_mass+=mass_field[id];  // should be zero
//   So.message_screen("Mass left to be assigned =", total_mass);
//   cout<<endl;
  
  mass_field.clear();mass_field.shrink_to_fit();
  trout.close();

  prop.clear();
  prop.shrink_to_fit();


  string fname_mass_function_Y = this->Output_directory+"tracer_mock_mass_function.txt";
  this->tracer.NOBJS=Ntracers;
  this->tracer.get_mass_function(fname_mass_function_Y);

#endif
  

  // ************************************************************************************************************************
  
#ifdef method4


#ifdef _use_mass_function_
  this->So.message_screen("Using reference mass function to complement assignment.");
  this->tracer.define_mass_bins();
  string massf_file=this->Output_directory+"tracer_ref_mass_function.txt";
  vector<real_prec>massf;
  int nbins_mf=this->File.read_file(massf_file,massf,1);
  vector<real_prec>mass_bins(nbins_mf,0);
  vector<real_prec>mfunc(nbins_mf,0);



#pragma omp parallel for
  for(int i=0; i< nbins_mf; ++i)
    {
      mass_bins[i]=massf[0+i*2];
      real_prec deltaM=this->tracer.MBmax[i]-this->tracer.MBmin[i];
      mfunc[i]=massf[1+i*2]*deltaM;// with n(M)*dM I get the right n(M) dividing by the max and sampling below
    }
  massf.clear();massf.shrink_to_fit();
  real_prec max_mf=static_cast<real_prec>(get_max(mfunc));

  vector<real_prec>mcount(nbins_mf,0);
  for(int i=0; i< nbins_mf; ++i)
    for(int j=i; j< nbins_mf; ++j)
      mcount[i]+=mfunc[j];

  ofstream sal;
  this->So.message_screen("Measuring cumulative MF of the reference:");
  sal.open(massf_file+"_cumulative");
  for(int i=0; i< nbins_mf; ++i)
    sal<<mass_bins[i]<<"\t"<<mcount[i]<<endl;
  sal.close();

  // Normalize to maximum
  for(int i=0; i< nbins_mf; ++i)
    mfunc[i]/=max_mf;


  mass_max=mass_bins[mass_bins.size()-1];
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, mass_bins.size());
  gsl_spline_init (spline, &mass_bins[0], &mfunc[0], mass_bins.size());

#endif


  

  this->So.message_screen("Assigning masses within cells");

  vector<s_rinfo>  new_mass_info(this->NGRID);


  ULONG IDD=16751686;
  /*
  real_prec massid=mass_field[IDD];
  real_prec auxm=0;
  for(int i=0; i<mass_info[IDD].mass_tracers.size();++i)
    {
      auxm+=mass_info[IDD].mass_tracers[i];
      cout<<mass_info[IDD].mass_tracers[i]<<"   "<<auxm<<"  "<<massid<<endl;
  }

*/

  real_prec total_mass=0;
#pragma omp parallel for reduction(+:total_mass)
  for(ULONG id=0; id < this->NGRID; ++id)
    total_mass+=mass_field[id];






  for(int ie=0;ie<=0;++ie)
     {
      real_prec expo=-1.0;//-0.8+(ie+0.5)*(-1.3+0.8)/20.;

  for(ULONG id=0; id< this->NGRID; ++id)
    {
      int nobs=v_ncounts[id];  // number of objects in the cell

      // If Neff>nobs, then we have less real particles than effective ones (i.e, Mid/Mmin), and mass can be distributed with some randomicity, respecting the limits mlim and Mid
      // If Neff<=nobs, then most likely all particles have the same mass, Mmin, plus some random selected from Mid-Mlim and distributed among all particles
      // if(nobs>1)
      // 	cout<<id<<"   Neff = "<<Neff<<"    Ntr = "<<nobs<<"    Mid = "<<mass_field[id]<<"    Lefts"<<mass_field[id]-static_cast<real_prec>(nobs)*mass_min<<endl;
      if(nobs>1)
	{
	  
	  int Neff=static_cast<int>(floor(mass_field[id]/mass_min)); //Total mass divided by the minimum = effective number of halos with Mlim
	  int alpha=static_cast<int>(floor(mass_field[id]/(nobs*mass_min))); //Total mass divided by the minimum = effective number of halos with Mlim
          real_prec epsilon=mass_field[id]-static_cast<real_prec>(nobs)*mass_min;


	  
	  bool flag=false;
	  while(false==flag)
	    {
	      new_mass_info[id].mass_tracers.resize(nobs,0);  //assign the minimum mass to all particles in the cell, to start
	      real_prec mass_aux=0;

              flag=true;

	      for(int i=0;i<nobs-1;++i)
		{
		  real_prec mass_tr;
                  real_prec epsilon_assigned;

#ifdef _use_mass_function_
		  real_prec prob=-10;
		  real_prec newxr=10;
		  while(prob<newxr)
		    {
#endif
		      real_prec xr;
		      //     if(alpha>1)
                      if(alpha==1)
                        {
                          epsilon_assigned=epsilon/(static_cast<real_prec>(nobs));
                          xr=gsl_rng_uniform(r);
                          mass_tr= mass_min+xr*epsilon_assigned;
                        }
                      else
                        {
                          real_prec mu = pow(static_cast<real_prec>(Neff)/static_cast<real_prec>(nobs), -expo);
                          xr=static_cast<real_prec>(gsl_ran_exponential(r,mu));
                          mass_tr= mass_min+xr*(mass_field[id] - mass_aux - (nobs-1-i)*mass_min);
                        }

#ifdef _use_mass_function_
		      if(mass_tr<mass_bins[0])mass_tr=mass_bins[0];
		      prob = static_cast<real_prec>(gsl_spline_eval (spline,mass_tr, acc));
		      newxr=static_cast<real_prec>(gsl_rng_uniform(r));
		    }
#endif

                  new_mass_info[id].mass_tracers[i] = mass_tr;    // assign the mass to a particle
                  mass_aux+=mass_tr;   // get the mass that has been assigned to this cell
                  }

              real_prec mass_last= mass_field[id]-mass_aux;//ramaining mass

              new_mass_info[id].mass_tracers[nobs-1]=mass_last; // Assign the remaining mass to the last particle within the cell.
	      
              mass_aux+=mass_last; // add that remaining mass to the total mass assigned
	      /*
		cout<<endl;
		if(id==IDD)
                for(int i=0;i<nobs;++i)
		cout<<id<<"   "<<i<<"   "<<alpha<<"  "<<new_mass_info[id].mass_tracers[i]<<"  "<<mass_aux<<"  "<<v_ncounts[id]<<"  "<<mass_field[id]<<endl;
	      */
	      
              if(mass_aux>mass_field[id])
                flag=false;
	      
              for(int i=0;i<nobs;++i)
		{
		  if(new_mass_info[id].mass_tracers[i]  < mass_min)
		    {
		      flag=false;
		      break;
		    }
		}
	      
	    }
	}
      else if(nobs==1)
	new_mass_info[id].mass_tracers.push_back(mass_field[id]);
    }
  
  
  So.DONE();
  this->So.message_screen("Writting to file");
  
  this->tracer.Halo.resize(Ntracers);
  vector<int>counter(this->NGRID,0);
  
  real_prec assigned_mass=0;
  for(ULONG i=0; i<  Ntracers; ++i)
    {
      ULONG id=index[i];
      real_prec mass_tracer;
      mass_tracer=new_mass_info[id].mass_tracers[counter[id]];
      assigned_mass+=mass_tracer;
      counter[id]++;
      this->tracer.Halo[i].mass=mass_tracer;
      //this->tracer.Halo[i].GridID = id;
      //prop[7+i*NCOLS]=mass_tracer;
      // for(ULONG i = 0; i<Ntracers; ++i)
      // 	{
      // 	  for(int j = 0; j<NCOLS; ++j)
      // 	    trout<<prop[j+i*NCOLS]<<"\t";
      // 	  trout<<endl;
      // 	}
    }
  counter.clear();counter.shrink_to_fit();
  So.DONE();
  trout.close();

  
  So.message_screen("Mass left to be assigned", total_mass-assigned_mass);
  cout<<endl;
  

  this->So.DONE();
  



  string fname_mass_function_Y = this->Output_directory+"tracer_mock_mass_function.txt";
  this->tracer.NOBJS=Ntracers;
  this->tracer.get_mass_function(fname_mass_function_Y);

 }

  mass_field.clear();mass_field.shrink_to_fit();

  index.clear();
  index.shrink_to_fit();
  prop.clear();
  prop.shrink_to_fit();


#endif






}


// ***************************************************************************************************************************************************************************************
// ***************************************************************************************************************************************************************************************
// ***************************************************************************************************************************************************************************************
// ***************************************************************************************************************************************************************************************
// ***************************************************************************************************************************************************************************************
// ***************************************************************************************************************************************************************************************
