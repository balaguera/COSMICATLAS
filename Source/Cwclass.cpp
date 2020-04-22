//##################################################################################
//##################################################################################
/** @file Cwclass.cpp
 *
 *  @brief Cosmic web classification
 *  @author: Andrés Balaguera-Antolínez, Francisco-Shu Kitaura, IAC, 2017-2019
 */

# include "../Headers/Cwclass.h"

int sign(real_prec x)
{

    int ans=1;
    if(x<0)
        ans=-1;
    return ans;

}

//##################################################################################
//##################################################################################
void Cwclass::set_params_cwclass(Params par)
{
  So.message_screen("Loading parameters for CWclass");
  this->Nft=par._Nft();
  this->lambdath=par._lambdath();
  this->lambdath_v=par._lambdath_v();
  this->n_sknot_massbin=par._n_sknot_massbin();
  this->n_vknot_massbin=par._n_vknot_massbin();
  this->Lbox=par._Lbox();
  this->cwt_used=par._cwt_used();
  this->cwv_used=par._cwv_used();
  this->n_cwt=par._n_cwt();
  this->n_cwv=par._n_cwv();
  this->NGRID=par._NGRID();
  this->NTT=static_cast<ULONG>(this->Nft*this->Nft*(this->Nft/2+1));
  this->params=par;
  So.DONE();
}





//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
void Cwclass::get_bias_terms(const vector<real_prec>&delta)
{

  So.message_screen("Computing bias terms");
  
  vector<real_prec>phi(this->NGRID,0);

  real_prec max_C1;
  real_prec min_C1;
  
  
#ifdef _USE_DELTA2_
  So.message_screen("Computing ð²");
  this->DELTA2.resize(this->NGRID,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->NGRID;++i)
    this->DELTA2[i]=delta[i]*delta[i];
  So.DONE();
#ifdef _WRITE_DELTA2_
  File.write_array(this->params._Output_directory()+"DELTA2", DELTA2);
#endif
  So.message_screen("DELTA2-term :");
  max_C1=get_max<real_prec>(DELTA2);
  So.message_screen("Maximum ð² =", max_C1);
  min_C1=get_min<real_prec>(DELTA2);
  So.message_screen("Minimum ð² =", min_C1);
  cout<<endl;
#endif



  
#ifdef _USE_DELTA3_
  So.message_screen("Computing ð³");
  this->DELTA3.resize(this->NGRID,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->NGRID;++i)
    DELTA3[i]=delta[i]*delta[i]*delta[i];
  So.DONE();
#ifdef _WRITE_DELTA3_
  File.write_array(this->params._Output_directory()+"DELTA3", DELTA3);
#endif
  So.message_screen("DELTA3-term :");
  max_C1=get_max<real_prec>(DELTA3);
  So.message_screen("Maximum ð³  =", max_C1);
  min_C1=get_min<real_prec>(DELTA3);
  So.message_screen("Minimum ð³ =", min_C1);
  cout<<endl;
#endif
  
  
  // Here we do not need the eigenvalues. Hence we do not resize these vectors
  So.message_screen("Computing gravitational potential");
  PoissonSolver(this->Lbox, this->Nft,delta,phi);
  So.DONE();
  
#ifdef _USE_S2_
  So.message_screen("Computing s²");
  this->S2.resize(this->NGRID,0);
  So.DONE();
#endif
  
#ifdef _USE_S3_
  So.message_screen("Computing s³");
  this->S3.resize(this->NGRID,0);
  So.DONE();
#endif
  
#ifdef _USE_NABLA2DELTA_
  So.message_screen("Computing NABLA2DELTA term");
  this->N2D.resize(this->NGRID,0);
#endif

#if (defined (_USE_NABLA2DELTA_) || defined (_USE_S2_) || defined (_USE_S3_) || defined (_USE_S2DELTA_))
  EigenValuesTweb_bias(this->Nft,this->Lbox,delta, phi,this->S2, this->S3, this->N2D);
  So.DONE();
#endif
  
  

  
  
#ifdef _USE_S2_
#ifdef _WRITE_S2_
  File.write_array(this->params._Output_directory()+"S2", S2);
#endif
  /*  // This will be written in Bam.pp, member function get_new_min_max_limits()
  So.message_screen("s²-term :");
  max_C1=get_max<real_prec>(S2);
  So.message_screen("Maximum s²  =", max_C1);
  min_C1=get_min<real_prec>(S2);
  So.message_screen("Minimum s² =", min_C1);
  cout<<endl;
  */
#endif
  
#ifdef _USE_S3_
#ifdef _WRITE_S3_
  File.write_array(this->params._Output_directory()+"S3", S3);
#endif
  /*
  So.message_screen("s³-term :");
  max_C1=get_max<real_prec>(S3);
  So.message_screen("Maximum s³  =", max_C1);
  min_C1=get_min<real_prec>(S3);
  So.message_screen("Minimum s³ =", min_C1);
  cout<<endl;
  */
#endif
  

  
#ifdef _USE_S2DELTA_
  this->S2DELTA.resize(this->NGRID,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->NGRID;++i)
    S2DELTA[i]=S2[i]*delta[i];
#ifdef _WRITE_S2DELTA_
  File.write_array(this->params._Output_directory()+"S2DELTA", S2DELTA);
#endif
  /*
  So.message_screen("s²ð-term :");
  max_C1=get_max<real_prec>(S2DELTA);
  So.message_screen("Maximum s²ð  =", max_C1);
  min_C1=get_min<real_prec>(S2DELTA);
  So.message_screen("Minimum s²ð =", min_C1);
  cout<<endl;
  */
#endif
  

#ifdef _USE_NABLA2DELTA_
#ifdef _WRITE_NABLA2DELTA_
  File.write_array(this->params._Output_directory()+"NABLA2DELTA", N2D);
#endif
  /*
    So.message_screen("Nabla²ð:");
    max_C1=get_max<real_prec>(N2D);
    So.message_screen("Maximum NABLA²ð  =", max_C1);
    min_C1=get_min<real_prec>(N2D);
    So.message_screen("Minimum Nabla²ð =", min_C1);
    cout<<endl;
  */
#endif
  
}


//##################################################################################
//##################################################################################
//##################################################################################


void Cwclass::do_CWC(const vector<real_prec>&delta)
{


  vector<real_prec>phi(this->NGRID,0);


  So.message_screen("Computing gravitational potential");
#ifdef _USE_ZERO_PADDING_POT_
  So.message_screen("using zero-padding");
#endif

  PoissonSolver(this->Lbox, this->Nft,delta,phi);
  So.DONE();
   
  this->lambda1.resize(this->NGRID,0);
  this->lambda2.resize(this->NGRID,0);
  this->lambda3.resize(this->NGRID,0);
  
  
#ifdef _USE_DELTA2_
    this->DELTA2.resize(this->NGRID,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<this->NGRID;++i)
      DELTA2[i]=delta[i]*delta[i];
#endif
    
#ifdef _USE_DELTA3_
    this->DELTA3.resize(this->NGRID,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<this->NGRID;++i)
      DELTA3[i]=delta[i]*delta[i]*delta[i];
#endif
    
    
    
    
    So.message_screen("Computing Eigenvalues of tidal field");
    // Here we do not compute S2 or N2D. Therefore we do not resize these containers.
    // #ifdef _USE_S2_
    //   this->S2.resize(this->NGRID,0);
    // #endif
    
    // #ifdef _USE_S3_
    //   this->S3.resize(this->NGRID,0);
    // #endif
    
    
    // #ifdef _USE_NABLA2DELTA_
    //       this->N2D.resize(this->NGRID,0);
    // #endif
    
    EigenValuesTweb(this->Nft,this->Lbox,delta, phi,this->lambda1,this->lambda2,this->lambda3, this->S2, this->S3, this->N2D);
    So.DONE();
    
  phi.clear();
  phi.shrink_to_fit();


#ifdef _USE_S2DELTA_
    this->S2DELTA.resize(this->NGRID,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<this->NGRID;++i)
      S2DELTA[i]=S2[i]*delta[i];
#endif


    real_prec max_C1;
    real_prec min_C1;


  /*
  max_C1=get_max<real_prec>(lambda1);
   So.message_screen("Maximum lambda1 =", max_C1);
   max_C1=get_min<real_prec>(lambda1);
   So.message_screen("Minimum lambda1 =", max_C1);
   max_C1=get_max<real_prec>(lambda2);
   So.message_screen("Maximum lambda2 =", max_C1);
   max_C1=get_min<real_prec>(lambda2);
   So.message_screen("Minimum lambda2 =", max_C1);
   max_C1=get_max<real_prec>(lambda3);
   So.message_screen("Maximum lambda3 =", max_C1);
   max_C1=get_min<real_prec>(lambda3);
   So.message_screen("Minimum lambda3 =", max_C1);
*/



#ifdef _USE_S2_
  So.message_screen("S2-term :");
  max_C1=get_max<real_prec>(S2);
  So.message_screen("Maximum s2  =", max_C1);
  max_C1=get_min<real_prec>(S2);
  So.message_screen("Minimum s2 =", max_C1);
  cout<<endl;
  if(this->step==0)
   File.write_array(this->params._Output_directory()+"S2", S2);

#endif

#ifdef _USE_S3_
  So.message_screen("S3-term :");
  max_C1=get_max<real_prec>(S3);
  So.message_screen("Maximum s3  =", max_C1);
  max_C1=get_min<real_prec>(S3);
  So.message_screen("Minimum s3 =", max_C1);
  cout<<endl;
  if(this->step==0)
   File.write_array(this->params._Output_directory()+"S3", S3);

#endif


#ifdef _USE_NABLA2DELTA_
  So.message_screen("Nabala² delta:");
  max_C1=get_max<real_prec>(N2D);
  So.message_screen("Maximum s2  =", max_C1);
  max_C1=get_min<real_prec>(N2D);
  So.message_screen("Minimum s2 =", max_C1);
  cout<<endl;
  if(this->step==0)
    File.write_array(this->params._Output_directory()+"NABLA2DELTA", N2D);

#endif






#ifdef _USE_INVARIANT_TIDAL_FIELD_I_
    So.message_screen("Invariant Tidal field I");
    this->Invariant_TF_I.resize(this->NGRID,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->NGRID;++index)
     {
      real_prec x=invariant_field_I(this->lambda1[index],this->lambda2[index],this->lambda3[index]);
#ifdef _USE_EXPONENT_INVARIANT_I_
      int sign_x=1;
#ifdef _USE_SIGN_INVARIANT_I_
      sign_x=sign(x);
#endif
      this->Invariant_TF_I[index]=sign_x*pow(abs(x), EXPONENT_INVARIANT_I);
#else
     this->Invariant_TF_I[index]=x;
#endif
     }


#ifdef _MAP_TO_INTERVAL_INV_I_
  {
  real_prec xmin=get_min(this->Invariant_TF_I);
  real_prec xmax=get_max(this->Invariant_TF_I);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->NGRID;++index)
      this->Invariant_TF_I[index]=(NEWMAX_INV_I-NEWMIN_INV_I)*(this->Invariant_TF_I[index]-xmin)/(xmax-xmin)+NEWMIN_INV_I;
}
#endif


  if(0==this->step)
#ifdef _USE_EXPONENT_INVARIANT_I_
      File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_I", this->Invariant_TF_I);
#else
      File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_I_original", this->Invariant_TF_I);
#endif

   if(this->step ==this->params._N_iterations_Kernel())
#ifdef _USE_EXPONENT_INVARIANT_I_
      File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_I_iteration"+to_string(this->step), this->Invariant_TF_I);
#else
      File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_I_original_iteration"+to_string(this->step), this->Invariant_TF_I);
#endif

  cout<<endl;
#endif



#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
    So.message_screen("Invariant Tidal field II");
    this->Invariant_TF_II.resize(this->NGRID,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->NGRID;++index)
     {
      real_prec x=invariant_field_II(this->lambda1[index],this->lambda2[index],this->lambda3[index]);
#ifdef _USE_EXPONENT_INVARIANT_II_
     int sign_x=1;
#ifdef _USE_SIGN_INVARIANT_II_
      sign_x=sign(x);
#endif
      this->Invariant_TF_II[index]=sign_x*pow(abs(x), EXPONENT_INVARIANT_II);
#else
     this->Invariant_TF_II[index]=x;
#endif
     }

#ifdef _MAP_TO_INTERVAL_INV_II_
 {
  real_prec xmin=get_min(this->Invariant_TF_II);
  real_prec xmax=get_max(this->Invariant_TF_II);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->NGRID;++index)
      this->Invariant_TF_II[index]=(NEWMAX_INV_II-NEWMIN_INV_II)*(this->Invariant_TF_II[index]-xmin)/(xmax-xmin)+NEWMIN_INV_II;
  }
#endif

  if(0==this->step)
#ifdef _USE_EXPONENT_INVARIANT_II_
      File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_II", this->Invariant_TF_II);
#else
      File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_II_original", this->Invariant_TF_II);
#endif

   if(this->step ==this->params._N_iterations_Kernel())
#ifdef _USE_EXPONENT_INVARIANT_II_
      File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_II_iteration"+to_string(this->step), this->Invariant_TF_II);
#else
      File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_II_original_iteration"+to_string(this->step), this->Invariant_TF_II);
#endif

  cout<<endl;
#endif





#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
  So.message_screen("Invariant Tidal field III");
  this->Invariant_TF_III.resize(this->NGRID,0);


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->NGRID;++index)
   {
     real_prec x=invariant_field_III(this->lambda1[index],this->lambda2[index],this->lambda3[index]);
#ifdef _USE_EXPONENT_INVARIANT_III_
     int sign_x=1;
#ifdef _USE_SIGN_INVARIANT_II_
     sign_x=sign(x);
#endif
     this->Invariant_TF_III[index]=sign_x*pow(abs(x), EXPONENT_INVARIANT_III);
#else
     this->Invariant_TF_III[index]=x;
#endif

  }


#ifdef _MAP_TO_INTERVAL_INV_III_
  {
  real_prec xmin=get_min(this->Invariant_TF_III);
  real_prec xmax=get_max(this->Invariant_TF_III);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->NGRID;++index)
      this->Invariant_TF_III[index]=(NEWMAX_INV_III-NEWMIN_INV_III)*(this->Invariant_TF_III[index]-xmin)/(xmax-xmin)+NEWMIN_INV_III;
}
#endif


  if(0==this->step)
#ifdef _USE_EXPONENT_INVARIANT_III_
      File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_III", this->Invariant_TF_II);
#else
      File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_III_original", this->Invariant_TF_III);
#endif

  if(this->step ==this->params._N_iterations_Kernel())
#ifdef _USE_EXPONENT_INVARIANT_III_
     File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_III_iteration"+to_string(this->step), this->Invariant_TF_III);
#else
     File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_III_original_iteration"+to_string(this->step), this->Invariant_TF_III);
#endif




#endif


#ifdef _USE_TIDAL_ANISOTROPY_
  So.message_screen("Invariant Tidal anisotropy");
  this->Tidal_Anisotropy.resize(this->NGRID,0);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->NGRID;++index)
    this->Tidal_Anisotropy[index]=tidal_anisotropy(this->lambda1[index],this->lambda2[index],this->lambda3[index]);
  max_C1=get_max<real_prec>(INV1);
  So.message_screen("Maximum Invariant =", max_C1);
  min_C1=get_min<real_prec>(INV1);
  So.message_screen("Minimum Invariant =", min_C1);
  cout<<endl;
  if(this->step==0)
   File.write_array(this->params._Output_directory()+"TIDAL_ANISOTROPY", INV);
#endif



  this->CWClass.clear();
#if defined (_USE_CWC_) || defined (_USE_MASS_KNOTS_)
  this->CWClass.resize(this->NGRID,0);
#endif


#ifdef _USE_CWC_
  ULONG nknots=0;
  ULONG nfilaments=0;
  ULONG nsheets=0;
  ULONG nvoids=0;
  ULONG nrest=0;
  
  So.message_screen("Extracting the Cosmic Web Classification:");
  So.message_screen("Lambda threshold = ", this->lambdath);

#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:nknots,nfilaments,nsheets,nvoids,nrest)
#endif
  for (ULONG index=0;index<this->NGRID;++index)
    {
      if (this->lambda1[index]>this->lambdath && this->lambda2[index]>this->lambdath && this->lambda3[index]>this->lambdath)
	{
	  this->CWClass[index]=I_KNOT;
	  nknots++;
	}
      
      else if (this->lambda1[index]<this->lambdath && lambda2[index]<this->lambdath && lambda3[index]<this->lambdath)
	{
	  this->CWClass[index]=I_VOID;
	  nvoids++;
	}
      
      else if ((lambda1[index]<this->lambdath && lambda2[index]<this->lambdath && lambda3[index]>this->lambdath) || (lambda1[index]<this->lambdath && lambda2[index]>this->lambdath && lambda3[index]<this->lambdath) || (lambda1[index]>this->lambdath && lambda2[index]<lambdath && lambda3[index]<this->lambdath))
	{
	  this->CWClass[index]=I_SHEET;
	  nsheets++;
	}
      
      else if ((lambda1[index]<this->lambdath && lambda2[index]>this->lambdath && lambda3[index]>this->lambdath) || (lambda1[index]>this->lambdath && lambda2[index]>this->lambdath && lambda3[index]<this->lambdath) || (lambda1[index]>this->lambdath && lambda2[index]<this->lambdath && lambda3[index]>this->lambdath))
	{
	  this->CWClass[index]=I_FILAMENT;
	  nfilaments++;
	}
      else
	nrest ++;
    }
  So.DONE();


  this->volume_knots=static_cast<real_prec>(nknots)/static_cast<real_prec>(this->NGRID); // Volume in knots / Total volume
  this->volume_filaments=static_cast<real_prec>(nfilaments)/static_cast<real_prec>(this->NGRID); // Volume in knots / Total volume
  this->volume_sheets=static_cast<real_prec>(nsheets)/static_cast<real_prec>(this->NGRID); // Volume in knots / Total volume
  this->volume_voids=static_cast<real_prec>(nvoids)/static_cast<real_prec>(this->NGRID); // Volume in knots / Total volume

//#ifdef _VERBOSE_
  So.message_screen("Summary of T-web classification, Lambda threshold = ", this->lambdath);
  So.message_screen("Knots (%) =", volume_knots*100.0);
  So.message_screen("Filaments (%) =", volume_filaments*100.0);
  So.message_screen("Sheets (%) =", volume_sheets*100.);
  So.message_screen("Voids (%) =", volume_voids*100.0);
//#endif


#ifdef _VERBOSE_
#ifndef _TEST_THRESHOLDS_RESIDUALS_
  int index= (this->step <=this->params._N_iterations_Kernel())  ?  this->step  : this->step - (this->params._N_iterations_Kernel())+1;
  string label_aux = this->step <= this->params._N_iterations_Kernel() ? "_iteration": "_realization";
  
  string file_cwc=this->params._Output_directory()+"CWC"+label_aux+to_string(index)+".txt";
  ofstream fcwc; fcwc.open(file_cwc.c_str());
  So.message_screen("Writing CWC fractions in file", file_cwc);
  fcwc<<static_cast<real_prec>(nknots)*100./static_cast<real_prec>(this->NGRID)<<endl;
  fcwc<<static_cast<real_prec>(nfilaments)*100./static_cast<real_prec>(this->NGRID)<<endl;
  fcwc<<static_cast<real_prec>(nsheets)*100./static_cast<real_prec>(this->NGRID)<<endl;
  fcwc<<static_cast<real_prec>(nvoids)*100./static_cast<real_prec>(this->NGRID)<<endl;
  fcwc.close();
  So.DONE();
#endif


  if(nrest>0)
    So.message_screen("Unclassified (%) =",static_cast<real_prec>(nrest)*100./static_cast<real_prec>(this->NGRID));
#endif

  // If we do not want to use the CWC but still use the infor from the knots, then
#elif !defined _USE_CWC_
#ifdef _USE_MASS_KNOTS_
  So.message_screen("Doing Cosmic Web Classification. Knots Only. Lambda threshold = ", this->lambdath);
#pragma omp parallel for
  for (ULONG index=0;index<this->NGRID;++index)
    if (lambda1[index]>this->lambdath && lambda2[index]>this->lambdath && lambda3[index]>this->lambdath)
      this->CWClass[index]=I_KNOT;
  So.DONE();
#endif
#endif

}




//##################################################################################
//##################################################################################
//##################################################################################


void Cwclass::do_CWC_V(vector<real_prec>&Vx,vector<real_prec>&Vy,vector<real_prec>&Vz)
{

    So.message_screen("V-WEB classification requested");

    // Transform the velocities to v / H(z), putting the units of H(z) with the cvel, such that the shear is dimensionless
/*
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0; i<this->NGRID;++i)
    {
       Vx[i]/=cvel;
       Vy[i]/=cvel;
       Vz[i]/=cvel;
    }

*/

  this->lambda1_vs.resize(this->NGRID,0);
  this->lambda2_vs.resize(this->NGRID,0);
  this->lambda3_vs.resize(this->NGRID,0);
  this->Divergence_VelField.resize(this->NGRID,0);
  So.message_screen("Computing Divergence and Eigenvalues of the shear of the velocity field");
  EigenValuesVweb(this->Nft,Lbox,Vx,Vy,Vz,this->Divergence_VelField, this->lambda1_vs,this->lambda2_vs,this->lambda3_vs);
  So.DONE();


#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
    So.message_screen("Invariant Shear Vfield I");
    this->Invariant_VS_I.resize(this->NGRID,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->NGRID;++index)
     {
      real_prec x=invariant_field_I(this->lambda1_vs[index],this->lambda2_vs[index],this->lambda3_vs[index]);
      int sign_x=1;
 #ifdef _USE_SIGN_INVARIANT_VS_I_
      sign_x=sign(x);
 #endif
      this->Invariant_VS_I[index]=sign_x*pow(abs(x), EXPONENT_INVARIANT_VS_I);
     }
  cout<<endl;
#ifdef _MAP_TO_INTERVAL_INV_SHEAR_I_
  {
  real_prec xmin=get_min(this->Invariant_VS_I);
  real_prec xmax=get_max(this->Invariant_VS_I);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->NGRID;++index)
      this->Invariant_VS_I[index]=(NEWMAX_INV_SHEAR_I-NEWMIN_INV_SHEAR_I)*(this->Invariant_VS_I[index]-xmin)/(xmax-xmin)+NEWMIN_INV_SHEAR_I;
}
#endif



#ifdef _USE_EXPONENT_INVARIANT_VS_I_
      File.write_array(this->params._Output_directory()+"INVARIANT_SHEAR_I", this->Invariant_VS_I);
#else
      File.write_array(this->params._Output_directory()+"INVARIANT_SHEAR_I_original", this->Invariant_VS_I);
#endif

#endif


#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
    So.message_screen("Invariant Shear Vfield II");
    this->Invariant_VS_II.resize(this->NGRID,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->NGRID;++index)
     {
      real_prec x=invariant_field_II(this->lambda1_vs[index],this->lambda2_vs[index],this->lambda3_vs[index]);
      int sign_x=1;
 #ifdef _USE_SIGN_INVARIANT_VS_II_
      sign_x=sign(x);
 #endif
      this->Invariant_VS_II[index]=sign_x*pow(abs(x), EXPONENT_INVARIANT_VS_II);
     }
  cout<<endl;
#ifdef _MAP_TO_INTERVAL_INV_SHEAR_II_
  {
  real_prec xmin=get_min(this->Invariant_VS_II);
  real_prec xmax=get_max(this->Invariant_VS_II);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->NGRID;++index)
      this->Invariant_VS_II[index]=(NEWMAX_INV_SHEAR_II-NEWMIN_INV_SHEAR_II)*(this->Invariant_VS_II[index]-xmin)/(xmax-xmin)+NEWMIN_INV_SHEAR_II;
}
#endif

#ifdef _USE_EXPONENT_INVARIANT_VS_II_
      File.write_array(this->params._Output_directory()+"INVARIANT_SHEAR_II", this->Invariant_VS_II);
#else
      File.write_array(this->params._Output_directory()+"INVARIANT_SHEAR_II_original", this->Invariant_VS_II);
#endif

#endif


#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
    So.message_screen("Invariant Shear Vfield III");
    this->Invariant_VS_III.resize(this->NGRID,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->NGRID;++index)
     {
      real_prec x=invariant_field_III(this->lambda1_vs[index],this->lambda2_vs[index],this->lambda3_vs[index]);
      int sign_x=1;
 #ifdef _USE_SIGN_INVARIANT_VS_III_
      sign_x=sign(x);
 #endif
      this->Invariant_VS_III[index]=sign_x*pow(abs(x), EXPONENT_INVARIANT_VS_III);
     }
  cout<<endl;
#ifdef _MAP_TO_INTERVAL_INV_SHEAR_III_
  {
  real_prec xmin=get_min(this->Invariant_VS_III);
  real_prec xmax=get_max(this->Invariant_VS_III);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->NGRID;++index)
      this->Invariant_VS_III[index]=(NEWMAX_INV_SHEAR_III-NEWMIN_INV_SHEAR_III)*(this->Invariant_VS_III[index]-xmin)/(xmax-xmin)+NEWMIN_INV_SHEAR_III;
}
#endif

  #ifdef _USE_EXPONENT_INVARIANT_VS_III_
      File.write_array(this->params._Output_directory()+"INVARIANT_SHEAR_III", this->Invariant_VS_III);
#else
      File.write_array(this->params._Output_directory()+"INVARIANT_SHEAR_III_original", this->Invariant_VS_III);
#endif
#endif





#ifdef _USE_CWC_V_
  this->CWClass_V.clear();
  this->CWClass_V.resize(this->NGRID,0);

  ULONG nknots=0;
  ULONG nfilaments=0;
  ULONG nsheets=0;
  ULONG nvoids=0;
  ULONG nrest=0;


  So.message_screen("Extracting the Cosmic Web Classification:");

#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:nknots,nfilaments,nsheets,nvoids,nrest)
#endif
  for (ULONG index=0;index<this->NGRID;++index)
    {
      if (this->lambda1_vs[index]>this->lambdath_v && this->lambda2_vs[index]>this->lambdath_v && this->lambda3_vs[index]>this->lambdath_v)
        {
          this->CWClass_V[index]=I_KNOT;
          nknots++;
        }

      else if (this->lambda1_vs[index]<this->lambdath_v && lambda2_vs[index]<this->lambdath_v && lambda3_vs[index]<this->lambdath_v)
        {
          this->CWClass_V[index]=I_VOID;
          nvoids++;
        }

      else if ((lambda1_vs[index]<this->lambdath_v && lambda2_vs[index]<this->lambdath_v && lambda3_vs[index]>this->lambdath_v) || (lambda1_vs[index]<this->lambdath_v && lambda2_vs[index]>this->lambdath_v && lambda3_vs[index]<this->lambdath_v) || (lambda1_vs[index]>this->lambdath_v && lambda2_vs[index]<lambdath_v && lambda3_vs[index]<this->lambdath_v))
        {
          this->CWClass_V[index]=I_SHEET;
          nsheets++;
        }

      else if ((lambda1_vs[index]<this->lambdath_v && lambda2_vs[index]>this->lambdath_v && lambda3_vs[index]>this->lambdath_v) || (lambda1_vs[index]>this->lambdath_v && lambda2_vs[index]>this->lambdath_v && lambda3_vs[index]<this->lambdath_v) || (lambda1_vs[index]>this->lambdath_v && lambda2_vs[index]<this->lambdath_v && lambda3_vs[index]>this->lambdath_v))
        {
          this->CWClass_V[index]=I_FILAMENT;
          nfilaments++;
        }
      else
        nrest ++;
    }
  So.DONE();
#ifdef _VERBOSE_
  So.message_screen("Summary of V-web classification, Lambda threshold = ", this->lambdath_v);
  So.message_screen("V_Knots (%) =",static_cast<real_prec>(nknots)*100./static_cast<real_prec>(this->NGRID));
  So.message_screen("V_Filaments (%) =",static_cast<real_prec>(nfilaments)*100./static_cast<real_prec>(this->NGRID));
  So.message_screen("V_Sheets (%) =",static_cast<real_prec>(nsheets)*100./static_cast<real_prec>(this->NGRID));
  So.message_screen("V_Voids (%) =",static_cast<real_prec>(nvoids)*100./static_cast<real_prec>(this->NGRID));
#endif


#ifndef _TEST_THRESHOLDS_RESIDUALS_
  int index= (this->step <=this->params._N_iterations_Kernel())  ?  this->step  : this->step - (this->params._N_iterations_Kernel())+1;
  string label_aux = this->step <= this->params._N_iterations_Kernel() ? "_iteration": "_realization";

  string file_cwc=this->params._Output_directory()+"CWC_V"+label_aux+to_string(index)+".txt";
  ofstream fcwc; fcwc.open(file_cwc.c_str());
  So.message_screen("Writing CWC-V fractions in file", file_cwc);
  fcwc<<static_cast<real_prec>(nknots)*100./static_cast<real_prec>(this->NGRID)<<endl;
  fcwc<<static_cast<real_prec>(nfilaments)*100./static_cast<real_prec>(this->NGRID)<<endl;
  fcwc<<static_cast<real_prec>(nsheets)*100./static_cast<real_prec>(this->NGRID)<<endl;
  fcwc<<static_cast<real_prec>(nvoids)*100./static_cast<real_prec>(this->NGRID)<<endl;
  fcwc.close();
  So.DONE();
#endif

  if(nrest>0)
    So.message_screen("Unclassified (%) =",static_cast<real_prec>(nrest)*100./static_cast<real_prec>(this->NGRID));

  // If we do not want to use the CWC but still use the infor from the knots, then
#elif !defined _USE_CWC_
#ifdef _USE_MASS_KNOTS_
  So.message_screen("Doing Cosmic Web Classification. Knots Only. Lambda threshold = ", this->lambdath);
#pragma omp parallel for
  for (ULONG index=0;index<this->NGRID;++index)
    if (lambda1_vs[index]>this->lambdath && lambda2_vs[index]>this->lambdath && lambda3_vs[index]>this->lambdath)
      this->CWClass_V[index]=I_KNOT;
  So.DONE();
#endif
#endif

}







//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
void Cwclass::get_Mk_collapsing_regions(vector<real_prec>&in, real_prec Nmean)
{
  /*
    Routine taken from HADRON
  */

    // If we are in _GET_REALIZATIONS_, the limits here must be fixed, since different realizations might have different min and max,
    
  cout<<endl;
  So.message_screen("Getting FoF from Knots");
  So.message_screen("mean number density =", Nmean);

  vector<real_prec>rho(in.size(),0);

  
  if(Nmean>0)  // If Nob is zero, then we are passing a density field. If not, convert delta to density with the value of Nobs passed
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<this->NGRID;++i)
        rho[i]= Nmean*(num_1+in[i]);
    }
  else
    rho=in;
  

    ULONG aux=0;
#pragma omp parallel for reduction(+:aux)
    for(ULONG i=0;i<this->NGRID;++i)
      if(this->CWClass[i]==I_KNOT)
	aux++;

  // Initialize to -1 to then identify whether a cell it is used or not
  // This array will be finally used to allocate the bin in wehich a mass knot is found.
  this->SKNOT_M_info.clear();
  this->SKNOT_M_info.shrink_to_fit();
  this->SKNOT_M_info.resize(this->NGRID, -1);
  
  // This vector contains the number of grids that belongs to a SKNOT. One expects
  // KNOT_MAX_SIZE superknots
  vector<ULONG>knot(KNOT_MAX_SIZE,0);  
  
  int k;
  real_prec knotmass;
  real_prec MKrange_max = -1e38;  // Maximum  and Min of Mass of DM in knots. Mass here is number of DM particles in the super knot
  real_prec MKrange_min = +1e38; // These two variables are updated when knot masses are computed

  for(ULONG i=0;i < this->NGRID; ++i)
    {
      if(this->CWClass[i]==I_KNOT && rho[i]>0)  // If cell is classified as knot, proceed. This must be the classification of the REFERENCE density  field. The >0 indicates that there might be some cells classified as knots but empty
	{
	  if(this->SKNOT_M_info[i] == -1) // if has not been already used in the fof, proceed
	    {
	      knot[0]=i;            // The first cell classified as knot. knot is always inizialized like this for every cell
	      this->SKNOT_M_info[i]=1;    //Mark the first knot as visited
	      int n=1; // Starts from 1. It gets added 1 if one of the neighbour cells is a knot. Percolation.
	      
              for(int j=0;j<n;++j)  // j runs over the neighbour cells of the current one, up to n, when n is increased inside the "if" when a neighbour knot is found
		{
		  if( (k= knot[j]-1)>=0 && k < this->NGRID)
                    if( (this->CWClass[k]==I_KNOT && rho[k]>0) && this->SKNOT_M_info[k]== -1)   // If density is greater than zero and cell not yet used
		      {
			knot[n] = k;   
			this->SKNOT_M_info[n] = 1; //mark this cell as already visited
			n++;
			if(n >= KNOT_MAX_SIZE)
			  So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
		      }
		  
		  if( (k = knot[j] + 1) >= 0 && k < this->NGRID)
                    if( (this->CWClass[k]==I_KNOT && rho[k]>0) && this->SKNOT_M_info[k] == -1)
		      {
			knot[n] = k;
			this->SKNOT_M_info[k] = 1;
			n++;
			if(n >= KNOT_MAX_SIZE)
			  So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
		      }
		  
		  if( (k = knot[j] - this->Nft) >= 0 && k < this->NGRID)
                    if( (this->CWClass[k]==I_KNOT && rho[k]>0) && this->SKNOT_M_info[k] == -1)
		      {
			knot[n] = k;
			this->SKNOT_M_info[k] = 1;
			n++;
			if(n >= KNOT_MAX_SIZE) 
			  So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
		      }
		  
		  if( (k = knot[j] + this->Nft) >= 0 && k < this->NGRID)
                    if((this->CWClass[k]==I_KNOT && rho[k]>0) && this->SKNOT_M_info[k] == -1)
		      {
			knot[n] = k;
			this->SKNOT_M_info[k] = 1;
			n++;
			if(n >= KNOT_MAX_SIZE) 
			  So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
		      }
		  
		  if( (k = knot[j] - this->Nft*this->Nft) >= 0 && k < this->NGRID)
                    if( (this->CWClass[k]==I_KNOT && rho[k]>0) && this->SKNOT_M_info[k] == -1)
		      {
			knot[n] = k;
			this->SKNOT_M_info[k] = 1;
			n++;
			if(n >= KNOT_MAX_SIZE)
			  So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
		      }
		  
		  if( (k = knot[j] + this->Nft*this->Nft) >= 0 && k < this->NGRID)
                    if((this->CWClass[k]==I_KNOT && rho[k]>0) && this->SKNOT_M_info[k] == -1)
		      {
			knot[n] = k;
			this->SKNOT_M_info[k] = 1;
			n++;
			if(n >= KNOT_MAX_SIZE)
			  So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
		      }
		}
	      
	      knotmass = 0;
              for(int j = 0; j < n; j++) // Compute the mass of the SK by adding the mass of its n parts. Note that if n=1, i.e, only
		knotmass += rho[knot[j]];

	      for(int j = 0; j < n; j++) //Assign to all cells in this superknot the mass of the sp they are in.
		this->SKNOT_M_info[knot[j]] = knotmass;
	      
	      // Adjust the interval for the bins in KMass to the current values min and max of the variable knotmass
	      if(MKrange_min > knotmass) MKrange_min = knotmass;
	      if(MKrange_max < knotmass) MKrange_max = knotmass;
	    }
	}
      else
	this->SKNOT_M_info[i] = 0;
    }

 knot.clear();
 knot.shrink_to_fit();



#ifdef _WRITE_MKNOTS_
 if(this->step==0)
  File.write_array(this->params._Output_directory()+"MASS_KNOTS", this->SKNOT_M_info);
#endif

 real_prec number_log;
#ifdef _DM_NEW_UNITS_
 number_log=0.0;
#else
 number_log=1.0;
#endif


 So.message_screen("Min and Max from FoF: ");
 So.message_screen("log Minimum Knot-mass (in units of log_10 Mass of DM particle) =", log10(MKrange_min+number_log));
 So.message_screen("log Maximim Knot-mass (in units of log_10 Mass of DM particle) =", log10(MKrange_max+number_log));

    // *** TYhis region is commnted if we want the current limits to be used in each step of the iteration
#ifndef _MODIFY_LIMITS_
  MKrange_max = MKMAX;
  MKrange_min = MKMIN;
  So.message_screen("log Used minimum Knot-mass (in units of log_10 Mass of DM particle) =", log10(MKrange_min+number_log));
  So.message_screen("log Used maximim Knot-mass (in units of log_10 Mass of DM particle) =", log10(MKrange_max+number_log));
  So.message_screen("");
#endif


  real_prec deltaMK= static_cast<real_prec>((log10(MKrange_max+number_log)-log10(MKrange_min+number_log))/static_cast<real_prec>(this->n_sknot_massbin));
  for(ULONG i = 0; i < this->NGRID; i++)
    if(this->SKNOT_M_info[i] == -1) // If this cell was not used, exit
      this->SKNOT_M_info[i] = 0;
    else
     // Get the bin in Knot mass. At this point SKNOT_M_info[i] encodes the mass of the Super knot where this cell is included
     this->SKNOT_M_info[i] = get_bin(log10(SKNOT_M_info[i]+number_log),log10(MKrange_min+number_log), this->n_sknot_massbin, deltaMK,true);


  So.DONE();


}


//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
void Cwclass::get_SigmaVel_collapsing_regions(const vector<real_prec>&in, vector<real_prec>&Vx, vector<real_prec>&Vy,vector<real_prec>&Vz, real_prec Nmean)
{
  /*
    Routine taken from HADRON
  */

    // If we are in _GET_REALIZATIONS_, the limits here must be fixed, since different realizations might have different min and max,

  cout<<endl;
  So.message_screen("Getting FoF from V-Knots");


  real_prec Hubble_function=this->s_cosmo_info.Hubble_parameter;
  real_prec cvel=Hubble_function/(cgs_km/cgs_Mpc);
  // Transform the velocities to v / H(z), putting the units of H(z) with the cvel, such that the shear is dimensionless

/*
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
for(int i=0; i<this->NGRID;++i)
  {
     Vx[i]/=cvel;
     Vy[i]/=cvel;
     Vz[i]/=cvel;
  }

  */


  vector<real_prec>rho(in.size(),0);


  //Now compute teh velocity dispersion
  So.message_screen("Computing Velocity dispersion");
  this->Dispersion_VelField.resize(this->NGRID,0);
  real_prec mean_vel_x=0;
  real_prec mean_vel_y=0;
  real_prec mean_vel_z=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:mean_vel_x,mean_vel_y,mean_vel_z)
#endif
  for(int i=0; i<this->NGRID;++i)
    {
      mean_vel_x+=Vx[i];
      mean_vel_y+=Vy[i];
      mean_vel_z+=Vz[i];
    }
  mean_vel_x/=static_cast<double>(this->NGRID);
  mean_vel_y/=static_cast<double>(this->NGRID);
  mean_vel_z/=static_cast<double>(this->NGRID);
  mean_vel_x=sqrt(pow(mean_vel_x,2)+pow(mean_vel_y,2)+pow(mean_vel_z,2));
  So.message_screen("Mean DM velocity = ", mean_vel_x, "Mpc/h");

  if(Nmean>0)  // If Nob is zero, then we are passing a density field. If not, convert delta to density with the value of Nobs passed
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<this->NGRID;++i)
        rho[i]= Nmean*(num_1+in[i]);
    }
  else
    rho=in;


    ULONG aux=0;
#pragma omp parallel for reduction(+:aux)
    for(ULONG i=0;i<this->NGRID;++i)
      if(this->CWClass_V[i]==I_KNOT)
        aux++;

  // Initialize to -1 to then identify whether it is used or not

  this->VDISP_KNOT_info.clear();
  this->VDISP_KNOT_info.resize(this->NGRID, -1);

  // This vector contains the number of grids that belongs to a SKNOT. One expects
  // KNOT_MAX_SIZE superknots
  vector<ULONG>knot(V_MAX_SIZE,0);

  int k;
  real_prec v_disp;
  real_prec VKrange_max = -1e38;  // Maximum  and Min of Mass of DM in knots. Mass here is number of DM particles in the super knot
  real_prec VKrange_min = +1e38; // These two variables are updated when knot masses are computed

  for(ULONG i=0;i < this->NGRID; ++i)
    {
      if(this->CWClass_V[i]==I_KNOT && rho[i]>0)  // If cell is classified as knot, proceed. This must be the classification of the REFERENCE density  field. The >0 indicates that there might be some cells classified as knots but empty
        {
          if(this->VDISP_KNOT_info[i] == -1) // if has not been already used in the fof, proceed
            {
              knot[0]=i;            // The first cell classified as knot. knot is always inizialized like this for every cell
              this->VDISP_KNOT_info[i]=1;    //Mark the first knot as visited
              int n=1; // Starts from 1. It gets added 1 if one of the neighbour cells is a knot. Percolation.

              for(int j=0;j<n;++j)  // j runs over the neighbour cells of the current one, up to n, when n is increased inside if a neirbour knot is found
                {
                  if( (k= knot[j]-1)>=0 && k < this->NGRID)
                    if( (this->CWClass_V[k]==I_KNOT && rho[k]>0) && this->VDISP_KNOT_info[k]== -1)   // If density is greater than zero and cell not yet used
                      {
                        knot[n] = k;
                        this->VDISP_KNOT_info[n] = 1; //mark this cell as already visited
                        n++;
                        if(n >= V_MAX_SIZE)
                          So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
                      }

                  if( (k = knot[j] + 1) >= 0 && k < this->NGRID)
                    if( (this->CWClass_V[k]==I_KNOT && rho[k]>0) && this->VDISP_KNOT_info[k] == -1)
                      {
                        knot[n] = k;
                        this->VDISP_KNOT_info[k] = 1;
                        n++;
                        if(n >= V_MAX_SIZE)
                          So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
                      }

                  if( (k = knot[j] - this->Nft) >= 0 && k < this->NGRID)
                    if( (this->CWClass_V[k]==I_KNOT && rho[k]>0) && this->VDISP_KNOT_info[k] == -1)
                      {
                        knot[n] = k;
                        this->VDISP_KNOT_info[k] = 1;
                        n++;
                        if(n >= V_MAX_SIZE)
                          So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
                      }

                  if( (k = knot[j] + this->Nft) >= 0 && k < this->NGRID)
                    if((this->CWClass_V[k]==I_KNOT && rho[k]>0) && this->VDISP_KNOT_info[k] == -1)
                      {
                        knot[n] = k;
                        this->VDISP_KNOT_info[k] = 1;
                        n++;
                        if(n >= V_MAX_SIZE)
                          So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
                      }

                  if( (k = knot[j] - this->Nft*this->Nft) >= 0 && k < this->NGRID)
                    if( (this->CWClass_V[k]==I_KNOT && rho[k]>0) && this->VDISP_KNOT_info[k] == -1)
                      {
                        knot[n] = k;
                        this->VDISP_KNOT_info[k] = 1;
                        n++;
                        if(n >= V_MAX_SIZE)
                          So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
                      }

                  if( (k = knot[j] + this->Nft*this->Nft) >= 0 && k < this->NGRID)
                    if((this->CWClass_V[k]==I_KNOT && rho[k]>0) && this->VDISP_KNOT_info[k] == -1)
                      {
                        knot[n] = k;
                        this->VDISP_KNOT_info[k] = 1;
                        n++;
                        if(n >= V_MAX_SIZE)
                          So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
                      }
                }

              v_disp = 0;

              for(int j = 0; j < n; j++) // Compute the mass*sigma² of the SK by adding the mass of its n parts
                 {
                   int index=knot[j]; // identify the index of the cell
                   real_prec vel_knots=sqrt(pow(Vx[index],2)+pow(Vy[index],2)+pow(Vz[index],2));
                   v_disp += rho[index]* pow(vel_knots - mean_vel_x,2)/static_cast<double>(n); // m sigma²
                }

              for(int j = 0; j < n; j++) //Assign to all cells in this superknot the mass of the sp they are in.
                this->VDISP_KNOT_info[knot[j]] = v_disp;

              // Adjust the interval for the bins in KMass to the current values min and max of the variable knotmass
              if(VKrange_min > v_disp) VKrange_min = v_disp;
              if(VKrange_max < v_disp) VKrange_max = v_disp;
            }
        }
      else
        this->VDISP_KNOT_info[i] = 0;
    }
 So.DONE();


 real_prec number_log;
#ifdef _DM_NEW_UNITS_
 number_log=0.0;
#else
 number_log=1.0;
#endif


 So.message_screen("Min and Max of Kinetic from FoF regions: ");
  So.message_screen("Minimum log(1+Kinetic) =",  log10(VKrange_min+number_log));
  So.message_screen("Maximim log(1+Kinetic) =",  log10(VKrange_max+number_log));

    // *** TYhis region is commnted if we want the current limits to be used in each step of the iteration
#ifndef _MODIFY_LIMITS_
  VKrange_max = VKMAX;
  VKrange_min = VKMIN;
  So.message_screen("log Used minimum Kinetic =", log10(VKrange_min+number_log));
  So.message_screen("log Used maximum Kinetic =", log10(VKrange_max+number_log));
  So.message_screen("");
#endif


  real_prec deltaMK= static_cast<real_prec>((log10(VKrange_max+number_log)-log10(VKrange_min+number_log))/static_cast<real_prec>(this->n_vknot_massbin));
  for(ULONG i = 0; i < this->NGRID; i++)
    {
      if(this->CWClass_V[i]==I_KNOT && rho[i]>0)
        {
          if(this->VDISP_KNOT_info[i] == -1) // If this cell was not used, exit
            {
              this->VDISP_KNOT_info[i] = 0;
              continue;
            }
          // GetA the bin in Knot mass. At this point SKNOT_M_info[i] encodes the mass of the Super knot where this cell is included
          if(log10(VDISP_KNOT_info[i]+number_log)<=log10(VKrange_max+number_log) && log10(VDISP_KNOT_info[i]+number_log)>=log10(VKrange_min+number_log) )
            {
              int j = static_cast<int>(floor((log10(VDISP_KNOT_info[i]+number_log)-log10(VKrange_min+number_log))/deltaMK));
              if(j == this->n_sknot_massbin) j--;
              // Assign the bin to the sknot_m_info to be used in the P(x,y,...)
              this->VDISP_KNOT_info[i] = j;
            }
        }
      else
        this->VDISP_KNOT_info[i] = 0;
    }


}

//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################


// This function helps to asses whether a cell should be counted (true) or not (false)
// according to the CW classification (member object) and the classifications
// requested in the input parameter file.
bool Cwclass::get_cell_classified(int sua, ULONG ig)
{
  bool res=false;
  if(this->cwt_used.size()==1 && this->cwt_used[0]==0)
    res=true;
  else
    {
      int classi=this->cwt_used[sua]; // This is the value passed fro the parameter file
      if(classi==0)
	res=true;
      
      else if(classi>0 && classi<5)  // This refers to the 4 typical classifications
	{
	  if(this->CWClass[ig]==classi)
	    res=true;
	}
      else         // This referes to combinations
	if(classi==234)
	  {
            if(this->CWClass[ig]==I_FILAMENT || this->CWClass[ig]==I_SHEET ||this->CWClass[ig]==I_VOID)
	      res=true;
	  }
	else if(classi==12)
	  {
            if(this->CWClass[ig]==I_KNOT || this->CWClass[ig]==I_FILAMENT)
	      res=true;
	  }
	else if(classi==23)
	  {
            if(this->CWClass[ig]==I_FILAMENT || this->CWClass[ig]==I_SHEET)
	      res=true;
	  }
	else if(classi==34)
	  {
            if(this->CWClass[ig]==I_SHEET || this->CWClass[ig]==I_VOID)
	      res= true;
	  }
	else if(classi==123)
	  {
            if(this->CWClass[ig]==I_KNOT || this->CWClass[ig]==I_FILAMENT || this->CWClass[ig]==I_SHEET)
	      res= true;
	  }
    }
  return res;
}

//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################

/** 
 * @brief This function provides the index identifying the CWT classification
 * as proposed in the parameter file, for each cell, and according
 * to the CWClass array
 */

int Cwclass::get_Tclassification(int ig)
{

  int ans=0;

  if(this->cwt_used.size()==1)
    ans=this->cwt_used[0];
  else
    {
      for(int ic=0;ic<this->cwt_used.size();++ic)// loop over the number of cwt required from the par file
	{
	  int cwc=this->cwt_used[ic];  // the cwt asked in par file. E.g, 1, 234
	  if(cwc<5)
	    {
	      if(this->CWClass[ig]==cwc)
		ans=ic;
	    }
	  else{
	    if(cwc==234)
	      {
                if(this->CWClass[ig]==I_FILAMENT || this->CWClass[ig]==I_SHEET || this->CWClass[ig]==I_VOID)
		  ans=ic;
	      }
	    else if(cwc==12)
	      {
                if(this->CWClass[ig]==I_KNOT || this->CWClass[ig]==I_FILAMENT)
		  ans=ic;
	      }
	    else if(cwc==23)
	      {
                if(this->CWClass[ig]==I_FILAMENT || this->CWClass[ig]==I_SHEET)
		  ans=ic;
	      }
	    else if(cwc==34)
	      {
                if(this->CWClass[ig]==I_SHEET || this->CWClass[ig]==I_VOID)
		  ans=ic;
	      }
	    
	    else if(cwc==123)
	      {
                if(this->CWClass[ig]==I_KNOT ||  this->CWClass[ig]==I_FILAMENT || this->CWClass[ig]==I_SHEET)
		  ans=ic;
	      }
	    
	  }
	}
    }
  
  return ans;
}


//##################################################################################

/**
 * @brief This function provides the index identifying the CWT classification
 * as proposed in the parameter file, for each cell, and according
 * to the CWClass array
 */

int Cwclass::get_Vclassification(int ig)
{

  int ans=0;

  if(this->cwv_used.size()==1)
    ans=this->cwv_used[0];
  else
    {
      for(int ic=0;ic<this->cwv_used.size();++ic)// loop over the number of cwt required from the par file
        {
          int cwc=this->cwv_used[ic];  // the cwt asked in par file. E.g, 1, 234
          if(cwc<5)
            {
              if(this->CWClass_V[ig]==cwc)
                ans=ic;
            }
          else{
            if(cwc==234)
              {
                if(this->CWClass_V[ig]==I_FILAMENT || this->CWClass_V[ig]==I_SHEET || this->CWClass_V[ig]==I_VOID)
                  ans=ic;
              }
            else if(cwc==12)
              {
                if(this->CWClass_V[ig]==I_KNOT || this->CWClass_V[ig]==I_FILAMENT)
                  ans=ic;
              }
            else if(cwc==23)
              {
                if(this->CWClass_V[ig]==I_FILAMENT || this->CWClass_V[ig]==I_SHEET)
                  ans=ic;
              }
            else if(cwc==34)
              {
                if(this->CWClass_V[ig]==I_SHEET || this->CWClass_V[ig]==I_VOID)
                  ans=ic;
              }

            else if(cwc==123)
              {
                if(this->CWClass_V[ig]==I_KNOT ||  this->CWClass_V[ig]==I_FILAMENT || this->CWClass_V[ig]==I_SHEET)
                  ans=ic;
              }

          }
        }
    }

  return ans;
}


//##################################################################################
//##################################################################################



