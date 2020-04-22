/** @file FftwFunctions.cpp
 *
 *  @brief Methods of the class FftwFunctions
 *
 *  This file contains the implementation of the methods of the class
 *  FftwFunctions, used to measure the 3D power spectrum (FKP and Yamamoto)
 *  based on the FFTW algrithm
*/

#include "../Headers/FftwFunctions.h"

#include <string.h>

#ifdef _USE_OMP_
#include <omp.h>
#endif

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// GENERAL STUFF
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


void FftwFunctions::set_healpix_pars(int Nres)
{
  //////////////////////////////////////////////////////////
  // passes the parameters associated to HEALPIX
  //////////////////////////////////////////////////////////
  nside=pow(2.,Nres);
  npixels=12.*nside*nside;
  area_pixel=4.*M_PI/npixels;
}
// *******************************************************************************************
void FftwFunctions::set_pars(real_prec lsi, real_prec kmy, int ndel_d, int ndel_w, int Nlb, int Nmb, bool masc, string mass, real_prec kmnb, real_prec kmxb)
{
  //////////////////////////////////////////////////////////
  // passes the parameters associated to the DFT
  //////////////////////////////////////////////////////////

  string statistics=this->params._statistics();

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  int nff;
  if(statistics!="Pk_y_ds")nff=Nft;
  // Get the number of grid-cells for the direct sum DFT as a function of kmax
  else if(statistics=="Pk_y_ds")
    {
      nff=(int)(kmy*lsi/M_PI); 
      if(nff%2!=0)nff++;
      this->Nft=nff;
      this->sgrid=nff;
    }
  
  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Get number of grid cells per dimension used for the bispectrum_fast given kmax
  int n_sgrid= static_cast<int>(kmxb*lsi/M_PI);
  
  if(statistics=="Bk_fkp_fast")
    {
      if(n_sgrid>this->Nft)
	{
          So.message_warning("Warning: kmax greater than Nyquist frequency");
          So.message_screen("Setting kmax to the Nyquist");
	  n_sgrid=this->Nft;
	  kmxb=Nft*M_PI/lsi;
	}
      if(n_sgrid%2!=0)n_sgrid++;
      this->sgrid=n_sgrid;
      this->kmax_bk = kmxb;
      this->kmin_bk = kmnb;         
      if(kmnb<2.*M_PI/lsi){
        So.message_warning("Warning: kmin smaller than fundamental mode");
        So.message_screen("Setting kmin as fundamental mode");
	this->kmin_bk = 2.*M_PI/lsi;     
      }
    }
  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(statistics=="Bk_fkp_fast")
    {
      // Get effective number of grid cells 
      // resulting from looping over half of the positive quadrant.
    // The result will be used to allocate memory for the arrays
    // used in the Bispectrum as performed by Jennifer.
      int new_sd=0;
      for(int i=0;i<this->sgrid;++i){
	for(int j=i;j<this->sgrid;++j){
	  for(int k=j;k<this->sgrid;++k){
	  if(i*i+j*j+k*k>0)new_sd++;
	  }
	}
      }
      this->new_sgrid=new_sd;
    }
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Initialize other private variables for the Pk
  this->Nft2=(int)(nff*nff);
  this->rNft = (real_prec)1.0 / this->Nft;
  this->NT=(int)(nff*nff*nff);
  this->NTT=(int)(nff*nff*(nff+2));
  this->kmax_y_ds=kmy;
  this->Lside =lsi;
  this->ndel_data = ndel_d;
  this->ndel_window = ndel_w;
  this->N_log_bins = Nlb;
  this->N_mu_bins = Nmb;
  this->MASS = mass;
  this->MAS_correction = masc;
  this->deltax = this->Lside / ((real_prec)this->Nft);
  this->deltay = this->deltax;
  this->deltaz = this->deltax;
  this->rdeltax = (real_prec)1.0/this->deltax;
  this->rdeltay = this->rdeltax;
  this->rdeltaz = this->rdeltax;
}


// *******************************************************************************************
// *******************************************************************************************
// *******************************************************************************************

#ifdef HEALPIX
inline void my_get_mean_density_interpolated(Healpix_Map<real_prec> &map, int new_n_dndz, real_prec redshift_min_sample, real_prec redshift_max_sample, real_prec ra, real_prec dec, real_prec redshift, vector< vector<real_prec> > &dndz_m, real_prec *nb)
{

  ///////////////////////////////////////////////////////////////
  // Function used when the mean number density is not tabulated
  ///////////////////////////////////////////////////////////////
  pointing point;
  
  real_prec fac=M_PI/180.0;
  point.phi=ra*fac;  
  point.theta=0.5*M_PI-dec*fac;
  long ipix=map.ang2pix(point);
  real_prec Delta_Z=(redshift_max_sample-redshift_min_sample)/dndz_m.size();
  int iz=(int)floor((float)((redshift-redshift_min_sample)/Delta_Z));
  *nb=dndz_m[iz][ipix];
}
#endif


// *******************************************************************************************
// *******************************************************************************************
// *******************************************************************************************
void FftwFunctions::free_fftw_vectors()
{
//  cout<<"Releasing memmory"<<endl;
  // This frees memmory. STILL MISSING VECTORS TO FREE MEMORY. HERE ONLY THOSE RELATED TO FKP
  fftw_free(this->data_out_g);
  fftw_free(this->data_out_r);
  // fftw_free(this->n);

  // this->data_g.clear();
  // this->data_g.shrink_to_fit();
  // this->data_r.clear();
  // this->data_r.shrink_to_fit();
  // this->SN.clear();
  // this->SN.shrink_to_fit();
  // this->Q.clear();
  // this->Q.shrink_to_fit();


}

// *******************************************************************************************
// *******************************************************************************************
// *******************************************************************************************

void FftwFunctions::fftw_vectors(bool use_random_catalog)
{
  //////////////////////////////////////////////////////////
  // This function initializes the input and 
  // output vectors for the FFTW.#
  // Must be called *AFTER* set_pars.            
  //////////////////////////////////////////////////////////


  // Array used in the FFTW
  // this->n= (int *)fftw_malloc(ic_rank*sizeof(float));
  // this->n[0]=this->Nft; 
  // this->n[1]=this->Nft;
  // this->n[2]=this->Nft;
  
  this->data_g.resize(this->NT);

#ifdef _MASS_WEIGHT_POWER_
    this->data_g_mw.resize(this->NT);
#endif


  if(true==this->measure_cross)
    data_gp.resize(this->NT);
  

  if(statistics=="Pk_ys" || statistics=="Pk_ybc"  || statistics=="Pk_ysc"){
    data_g_xx.resize(this->NT,0.0);
    data_g_yy.resize(this->NT,0.0);
    data_g_zz.resize(this->NT,0.0);
    data_g_xy.resize(this->NT,0.0);
    data_g_xz.resize(this->NT,0.0);
    data_g_yz.resize(this->NT,0.0);
  }

  if(statistics=="Pk_yb" || statistics=="Pk_ybc" ){
    data_g_xx.resize(this->NT,0.0);
    data_g_yy.resize(this->NT,0.0);
    data_g_zz.resize(this->NT,0.0);
    data_g_xy.resize(this->NT,0.0);
    data_g_xz.resize(this->NT,0.0);
    data_g_yz.resize(this->NT,0.0);    
    data_g_xxx.resize(this->NT,0.0);
    data_g_yyy.resize(this->NT,0.0);
    data_g_zzz.resize(this->NT,0.0);
    data_g_xxy.resize(this->NT,0.0);
    data_g_xxz.resize(this->NT,0.0);
    data_g_yyx.resize(this->NT,0.0);
    data_g_yyz.resize(this->NT,0.0);
    data_g_zzx.resize(this->NT,0.0);
    data_g_zzy.resize(this->NT,0.0);
    data_g_xyy.resize(this->NT,0.0);
    data_g_xzz.resize(this->NT,0.0);
    data_g_yzz.resize(this->NT,0.0);
    data_g_xyz.resize(this->NT,0.0);
    data_g_yxz.resize(this->NT,0.0);
    data_g_zxy.resize(this->NT,0.0);
  }
  
  if(statistics=="Pk_ybc" || statistics=="Pk_ysc"){
    cell_x.resize(this->NT,0.0);
    cell_y.resize(this->NT,0.0);
    cell_z.resize(this->NT,0.0);
  }


  this->data_out_g   =(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));

  if(true==this->measure_cross)
    data_out_gp   =(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
  

  if(statistics=="Pk_ys" || statistics=="Pk_ysc"){
    data_out_g_xx=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_yy=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_zz=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_xy=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_xz=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_yz=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
  }  
  
  if(statistics=="Pk_yb"  || statistics=="Pk_ybc"){
    data_out_g_xx=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_yy=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_zz=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_xy=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_xz=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_yz=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_xxx=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_yyy=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_zzz=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_xxy=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_xxz=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_yyx=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_yyz=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_zzx=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_zzy=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_xyy=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_xzz=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_yzz=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_xyz=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_yxz=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
    data_out_g_zxy=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
  }

  if(true==use_random_catalog)
     data_r.resize(this->NT);

  if(statistics=="Pk_ys" || statistics=="Pk_ysc"){
    data_r_xx.resize(this->NT,0.0);
    data_r_yy.resize(this->NT,0.0);
    data_r_zz.resize(this->NT,0.0);
    data_r_xy.resize(this->NT,0.0);
    data_r_xz.resize(this->NT,0.0);
    data_r_yz.resize(this->NT,0.0);
  }
  if(statistics=="Pk_yb" || statistics=="Pk_ybc" ){
    data_r_xx.resize(this->NT,0.0);
    data_r_yy.resize(this->NT,0.0);
    data_r_zz.resize(this->NT,0.0);
    data_r_xy.resize(this->NT,0.0);
    data_r_xz.resize(this->NT,0.0);
    data_r_yz.resize(this->NT,0.0);
    data_r_xxx.resize(this->NT,0.0);
    data_r_yyy.resize(this->NT,0.0);
    data_r_zzz.resize(this->NT,0.0);
    data_r_xxy.resize(this->NT,0.0);
    data_r_xxz.resize(this->NT,0.0);
    data_r_yyx.resize(this->NT,0.0);
    data_r_yyz.resize(this->NT,0.0);
    data_r_zzx.resize(this->NT,0.0);
    data_r_zzy.resize(this->NT,0.0);
    data_r_xyy.resize(this->NT,0.0);
    data_r_xzz.resize(this->NT,0.0);
    data_r_yzz.resize(this->NT,0.0);
    data_r_xyz.resize(this->NT,0.0);
    data_r_yxz.resize(this->NT,0.0);
    data_r_zxy.resize(this->NT,0.0);
  }

  data_out_r=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
  
  if(statistics=="Pk_y_ds"){
    data_g_out_y0.resize(this->NTT,0.0);
    data_g_out_y2.resize(this->NTT,0.0);
    data_g_out_y4.resize(this->NTT,0.0);
    data_r_out_y0.resize(this->NTT,0.0);
    data_r_out_y2.resize(this->NTT,0.0);
    data_r_out_y4.resize(this->NTT,0.0);
    SN_r_out_y2.resize(this->NTT,0.0);
    SN_r_out_y4.resize(this->NTT,0.0);
    SN_g_out_y2.resize(this->NTT,0.0);
    SN_g_out_y4.resize(this->NTT,0.0);
    data_g_y0.resize(this->NTT,0.0);
    data_g_y2.resize(this->NTT,0.0);
    data_g_y4.resize(this->NTT,0.0);
  }

  if(statistics=="Pk_fkp"){
    SN.resize(this->NT,0);
    Q.resize(this->NT,0);
  }
  
  
  if(statistics=="Bk_fkp_fast"){
    Array_corr.resize(this->NTT,1.0);
    Arraykx.resize(this->new_sgrid,0);
    Arrayky.resize(this->new_sgrid,0);
    Arraykz.resize(this->new_sgrid,0);
    Arraykk.resize(this->new_sgrid,0);
    ArrayID.resize(this->new_sgrid,0);
    VecArray.resize(this->new_sgrid,0);
    Kbin.resize(this->new_sgrid,0);
    Bmodes.resize(this->Nshells_bk,0);
    kkminID.resize(this->Nshells_bk,0);    
    kbins_bk.resize(this->Nshells_bk,0);
    Ngrids_bk.resize(this->Nshells_bk,0);

    
  }


}


// *******************************************************************************************
// *******************************************************************************************
// *******************************************************************************************

// original definition, accounts for potentially different dimensions for x,y and z axis
// This is a row-major order
#define ijk(i, j, k, nn1, nn2, nn3) ((int)((k)+(nn3)*( (j)+ (nn2)*(i))))

// column-major order, used when needed to compare against Florian's code.
//#define ijk(i, j, k, nn1, nn2, nn3) ((int)((i)+(nn1)*( (j)+ (nn2)*(k))))


// *****************************************************************************************
// *****************************************************************************************
// *****************************************************************************************

real_prec FftwFunctions::MAS(real_prec x)
{
  // seems it is not used 
  //////////////////////////////////////////////////////////
  // Mass assignment scheme. Three different MAS     
  // are available for interpolation of the density  
  // field. 
  //////////////////////////////////////////////////////////

  real_prec daux=fabs(x);
  real_prec ans;
  
  if(this->MASS=="NGP"){
    if(daux<0.5){
      ans=1.0;
    }
    if(daux==0.5){
      ans=0.5;
    }
    if(daux>0.5){
      ans=0.0;
    }
  }

  if(this->MASS=="CIC"){
    ans=(daux<1 ? 1.-daux: 0);
  }

  if(this->MASS=="TSC"){
    if(daux<0.5)ans=(0.75-x*x);
    if(daux>=0.5 && daux<1.5) ans= (0.5*pow(1.5-daux,2));
    if(daux>=1.5)ans=0.0;
  }

  if(this->MASS=="PCS"){
    if(daux>=0 && daux<1)ans=(1./6.)*(4-6*x*x+3.*pow(daux,3));
    else if(daux>=1 && daux<2)ans= (1./6.)*pow(2-daux,3);
    else ans=0.0;
  }
  return ans;
}


// *****************************************************************************************

inline real_prec FftwFunctions::MAS_NGP(real_prec x)
{
  //////////////////////////////////////////////////////////
  // NGP Mass assignment scheme. 
  //////////////////////////////////////////////////////////

  // takes the absolute value
  x=fabs(x);
  
  if(x<0.5) return 1.0;
  if(x>0.5) return 0.0;
  if(x==0.5) return 0.5;
}

// *****************************************************************************************

inline real_prec FftwFunctions::MAS_TSC(real_prec x)
{
  //////////////////////////////////////////////////////////
  // TSC Mass assignment scheme. 
  //////////////////////////////////////////////////////////

    // takes the absolute value
  x=fabs(x);

  if(x<0.5)
    return (0.75-x*x);
  if(x>=0.5 && x<1.5) {
    real_prec r = 1.5 - x;
    r *= r;
    return (0.5*r);
  }
  return 0;
}

// *****************************************************************************************
// *****************************************************************************************
inline real_prec FftwFunctions::MAS_CIC(real_prec x)
{
  //////////////////////////////////////////////////////////
  // CIC Mass assignment scheme. 
  //////////////////////////////////////////////////////////
  x=fabs(x);
  if(x<1)return 1.-x;
  else return 0.0;
}




// *****************************************************************************************
// *****************************************************************************************

inline real_prec FftwFunctions::MAS_PCS(real_prec x)
{
  //////////////////////////////////////////////////////////
  // PCS Mass assignment scheme.
  //////////////////////////////////////////////////////////
  real_prec ans;
  real_prec y=fabs(x);
  if(y>=0 && y<1)
    ans=(1./6.)*(4.0- 6.0*x*x + 3.*pow(y,3));
  else if(y>=1. && y<2.)
    ans= (1./6.)*pow(2.0-y, 3);
  else
    ans=0.;
  return ans;
}

// ************************************************************************************************
// ************************************************************************************************
//  TO BE DEPRECATED. PARAMETERS ARE DEFINED INSTEAD AT ONE OF THE CONSTRUCTOR
void FftwFunctions::set_bins(string type_of_binning)
{
  //////////////////////////////////////////////////////////
  // This function computes the parameters of the k-bins
  // e.g., width,, kmin, kmax, etc
  //////////////////////////////////////////////////////////

  this->binning = type_of_binning;
  this->deltak_x      = 2.*M_PI/this->Lside;
  this->deltak_y      = 2.*M_PI/this->Lside;
  this->deltak_z      = 2.*M_PI/this->Lside;
  this->deltak_0      = sqrt(pow(deltak_x,2)+pow(deltak_y,2)+pow(deltak_z,2))/sqrt(3.0); 
  this->DeltaK_data   = this->ndel_data*deltak_0; 
  this->DeltaK_window = this->ndel_window*deltak_0; 
  this->kmin          = 0.5*this->deltak_0;
  this->kmax          = sqrt(pow(0.5*this->Nft*deltak_x,2)+pow(0.5*this->Nft*deltak_y,2)+pow(0.5*this->Nft*deltak_z,2))/sqrt(3.0); 
  this->Deltal        = log10(this->kmax/this->kmin)/this->N_log_bins; 
  this->Deltamu       = 2.0/(static_cast<real_prec>(this->N_mu_bins));
  this->Nnp_data      = (type_of_binning=="log" ? this->N_log_bins : this->Nft/this->ndel_data/2); 
  this->Nnp_window    = (type_of_binning=="log" ? this->N_log_bins : this->Nft/this->ndel_window/2); 
  this->DeltaK_Bis    = this->DeltaK_data;// (this->kmax_bk-this->kmin_bk)/((real_prec)this->Nshells_bk)
  this->Nshells_bk =   static_cast<int>(this->kmax_bk/this->DeltaK_Bis); 

}

// *******************************************************************
// *******************************************************************
// *******************************************************************
// *******************************************************************

void FftwFunctions::cart_coordinates(s_parameters_box *s_box,s_data_structure *s_data,vector<s_Halo >&prop)
{
  
  /**
   * Transform the coordinates of the input catalogues to cartessian coordinates
   * returning the same input vector in which the three first columns correspond to
   * to the X,Y,Z coordinates and in the i_nbar (see parameter file)
   * the mean number density
   */
  
  cout.precision(12);
  
  
  
  gsl_spline *spline_zro, *spline_nbar;
  //  Healpix_Map<real_prec>map(log2(this->nside), RING);


  real_prec aXMAX=-1e6,aYMAX=-1e6, aZMAX=-1e6;
  real_prec aXMIN=+1e6,aYMIN=+1e6, aZMIN=+1e6;
  
  real_prec fac=M_PI/180.0;
  std::cout.precision(8);  
  vector <gsl_real>zzv =  s_data->zz_v;
  vector <gsl_real>dndz = s_data->dndz_v;
  // *****************************************************************************
  // Identify columns in the corresponding catalog                               *
  // *****************************************************************************
  int i_coord1     = (s_data->catalog=="data"? s_box->i_coord1_g : s_box->i_coord1_r);
  int i_coord2     = (s_data->catalog=="data"? s_box->i_coord2_g : s_box->i_coord2_r);
  int i_coord3     = (s_data->catalog=="data"? s_box->i_coord3_g : s_box->i_coord3_r);
  
  int i_mean_density= (s_data->catalog=="data"? s_box->i_mean_density_g : s_box->i_mean_density_r);
  string angles_units = (s_data->catalog=="data" ? s_box->angles_units_g : s_box->angles_units_r);
  
  // *****************************************************************************
  // Identify catalogue                                                          *
  // *****************************************************************************
  int n_columns = s_data->n_columns;
  int nlines = prop.size() / n_columns;

  vector<gsl_real> zz=s_box->zz;
  vector<gsl_real> rc=s_box->rc;

  int n_rc, n_zzv;
  
  // if needed initialize gsl structure to interpolate cmv dist

  if(s_box->use_random_catalog && !s_data->nbar_tabulated)
    {
      n_zzv = zzv.size();
      vector<gsl_real>x (n_zzv,0);
      vector<gsl_real>y (n_zzv,0);
      
      spline_nbar = gsl_spline_alloc (gsl_interp_linear, zzv.size());
      for(int i = 0; i < n_zzv; i++)
       {
         x[i] = zzv[i];
         y[i] = dndz[i];
       }
      gsl_spline_init (spline_nbar, &x[0], &y[0], n_zzv);
    }
  
  if(s_box->use_random_catalog)
    {
      n_rc = rc.size();
      vector<gsl_real>x (n_rc,0);
      vector<gsl_real>y (n_rc,0);
      for(int i = 0; i < n_rc; i++)
      {
        x[i] = rc[i];
        y[i] = zz[i];
      }
      spline_zro = gsl_spline_alloc (gsl_interp_linear, n_rc);
      switch(s_data->system_of_coordinates) {
      case(1):      
          gsl_spline_init (spline_zro, &x[0], &y[0], n_rc);
        break;
      case(2):
        gsl_spline_init (spline_zro, &y[0], &x[0], n_rc);
        break;
      }
    }


  time_t start;
  time (&start);
  
  


  int NTHREADS = omp_get_max_threads();


  So.message_screen("Transforming to cartessian coordinates");

    
#pragma omp parallel num_threads(NTHREADS)
  {
    int myID;
    real_prec *pprop;
    gsl_interp_accel *spline_acc_zro;
    gsl_interp_accel *spline_acc_nbar;

    spline_acc_zro = gsl_interp_accel_alloc ();
    if(s_box->use_random_catalog &&
       !s_data->nbar_tabulated) {
      spline_acc_nbar = gsl_interp_accel_alloc ();
    }
    
    myID = omp_get_thread_num();
    

    // notew for improving the speed of the code
    //HERE WE CAN MOVE THE SIWTCH OUTSIDE THE LOOPS AND MAKE THREE CASES


#pragma omp for schedule(static, NTHREADS) reduction(min:aXMIN, aYMIN, aZMIN) reduction(max:aXMAX, aYMAX, aZMAX) nowait
    for(int i=0;i<nlines;++i)
    {

      real_prec x,y,z,rr, zro, nbar, ra_s, dec_s;
      
//      pprop = &prop[i*n_columns];
      
      switch(s_data->system_of_coordinates){

      // Choose Cartesian coordinates:
      case(0):   
        x=prop[i].coord1; //pprop[i_coord1];
        y=prop[i].coord1; //pprop[i_coord2];
        z=prop[i].coord1; //pprop[i_coord3];
        //compute nbar:
        if(s_box->use_random_catalog){
          if(s_data->nbar_tabulated)nbar=prop[i].weight1;//pprop[i_mean_density];
	  // If nbar is not tabulated in the case (i.e., cartesian coords), then we do not generate an estimate, 
	  // for the redshift-axis is not necessary defined in the z-component. Then, if not tabulated, we set it to 1:
	  else{nbar=1;}
        }
	// If a random catalog is not used, then use the
	// mean number density obtained from the information of L and Ngal
	// passed through the structure s_data
        else{nbar=s_data->mean_density;}
     
        break;
      
	
        //Choose Equatorial coordinates	
      case(1): 
        rr=pprop[i_coord3];
        if(angles_units=="D"){
          ra_s=prop[i].coord1;//pprop[i_coord1];
          dec_s=prop[i].coord2;//pprop[i_coord2];
        }
        else{
          ra_s=prop[i].coord1/fac;
          dec_s=prop[i].coord2/fac;
        }

        // *******************************************************************
        // Transform to cartesian coord:
        Go.equatorial_to_cartesian(ra_s,dec_s,rr,&x, &y, &z);  
     
        // *******************************************************************
        // Compute the mean number density if not tabulated :          
        if(s_box->use_random_catalog){
          if(s_data->nbar_tabulated)nbar=prop[i].mean_density;
          else{
            if(s_box->constant_depth){
              zro= gsl_spline_eval(spline_zro, prop[i].coord3, spline_acc_zro);
              nbar= gsl_spline_eval(spline_nbar, zro, spline_acc_nbar); 
            }
            else{
	      // Interpolate to get redshift
              zro= gsl_spline_eval(spline_zro, prop[i].coord3, spline_acc_zro);
	      // Get ane stimate of the men number density:
	      #ifdef HEALPIX
	      //	      my_get_mean_density_interpolated(map, s_box->new_n_dndz,s_box->redshift_min_sample, s_box->redshift_max_sample,ra_s, dec_s, zro, dndz_m, &nbar);
	      #endif
            }
          }
        }
	// If a random catalog is not used, then use the
	// mean number density obtained from the information of L and Ngal
	// passed through the structure s_data
        else{nbar=s_data->mean_density; }     
        break;
     
      case(2):// Choose pseudo equatorial coordinates, and follow instructions as above
    
        if(angles_units=="D"){
          ra_s=prop[i].coord1;
          dec_s=prop[i].coord2;
        }
        else{
          ra_s=prop[i].coord1/fac;
          dec_s=prop[i].coord2/fac;
        }
     
        // ******************************************************************
        // Transform to comoving distance *
        rr=gsl_spline_eval(spline_zro, prop[i].coord3, spline_acc_zro);

       // *******************************************************************
        // Transform to cartesian coord: 
        Go.equatorial_to_cartesian(ra_s, dec_s, rr, &x, &y, &z);
     
        // **************************************************************************
        // Compute the mean number density if not tabulated                         *
     
        if(s_box->use_random_catalog){
          if(s_data->nbar_tabulated)
              nbar=prop[i].mean_density;
          else 
            if(s_box->constant_depth){
              zro=prop[i].coord3;
              nbar = gsl_spline_eval(spline_nbar, zro, spline_acc_nbar);
            }
            else{
              zro=pprop[i_coord3];
	      #ifdef HEALPIX
	      //              my_get_mean_density_interpolated(map, s_box->new_n_dndz,s_box->redshift_min_sample, s_box->redshift_max_sample,ra_s, dec_s, zro, dndz_m, &nbar);
	      #endif
            }
	}
	else{nbar=s_data->mean_density;}
        break;
	
     
     
      case(3):// Choose pseudo equatorial coordinates: use redshift as radial coordinate
        rr=pprop[i_coord3];
        if(angles_units=="D"){
          ra_s=prop[i].coord1;
          dec_s=prop[i].coord2;
        }
        else{
          ra_s=prop[i].coord1/fac;
          dec_s=prop[i].coord2/fac;
        }

        // **************************************************************************
        // Transform to cartesian coord. *
        // **************************************************************************
        Go.equatorial_to_cartesian(ra_s, dec_s, rr, &x, &y, &z);
        // **************************************************************************
        // Compute the mean number density if not tabulated                         *
        // **************************************************************************
        if(s_box->use_random_catalog){
          if(s_data->nbar_tabulated)
              nbar=prop[i].mean_density;
          else{
            if(s_box->constant_depth){
              zro=prop[i].coord3;
              nbar = gsl_spline_eval(spline_nbar, zro, spline_acc_nbar);            
            }
            else{
              zro=prop[i].coord3;
#ifdef HEALPIX
	      // my_get_mean_density_interpolated(map, s_box->new_n_dndz,s_box->redshift_min_sample, s_box->redshift_max_sample,ra_s, dec_s, zro, dndz_m, &nbar);
#endif
            }
          }
        }
        else{nbar=s_data->mean_density;}
        break;
      }
      
      if(s_data->system_of_coordinates != 0)
	{
          prop[i].coord1=x;
          prop[i].coord2=y;
          prop[i].coord3=z;
	}
      
      
      // This if statetment is important: 
      // If we do not have a fourth column in the input catalog
      // we will be calling a portion of memory not allocated
      // or allocated with other information. In other words, a mess ;-)
      if(s_box->use_random_catalog) pprop[i_mean_density]=nbar;

   
      // Identify min and max coordinates in order to set size of cube.
      if(x > aXMAX)
        aXMAX = x;
      if(y > aYMAX)
        aYMAX = y;
      if(z > aZMAX)
        aZMAX = z;
      if(x < aXMIN)
        aXMIN = x;
      if(y < aYMIN)
        aYMIN = y;
      if(z < aZMIN)
        aZMIN = z;    
    } // END parallel loop
    
  } // END parallel region



  if(s_data->catalog=="data"){
    this->Xmin=aXMIN;  
    this->Ymin=aYMIN;
    this->Zmin=aZMIN; 
    this->Xmax=aXMAX;  
    this->Ymax=aYMAX;
    this->Zmax=aZMAX;  
    
    // Determine dimension of box from the data
    // Get min and maxs:
    real_prec llx=fabs(aXMAX-aXMIN);
    real_prec lly=fabs(aYMAX-aYMIN);
    real_prec llz=fabs(aZMAX-aZMIN);
    llx=max(llx,lly);
    lly=max(lly,llz);
    llz=max(llx,lly);
    this->Lside_data=llz;
    
    
    So.message_screen("Computed box:");
    cout <<YELLOW<<" x [ " << aXMIN << " : " << aXMAX << " ] " << endl;
    cout <<YELLOW<< " y [ " << aYMIN << " : " << aYMAX << " ] " << endl;
    cout <<YELLOW<< " z [ " << aZMIN << " : " << aZMAX << " ] " << endl;
    cout <<YELLOW<< "Lside estimated = " <<CYAN<<this->Lside_data<<RESET<<endl;
    this->Xoffset = 0.5*(aXMAX+aXMIN);
    this->Yoffset = 0.5*(aYMAX+aYMIN);
    this->Zoffset = 0.5*(aZMAX+aZMIN);
    So.message_screen("Xoffset = ",this->Xoffset);
    So.message_screen("Yoffset = ",this->Yoffset);
    So.message_screen("Zoffset = ",this->Zoffset);
  } 

  
  std::cout<<YELLOW;
  time_t end; time (&end);
  real_prec diff = difftime(end,start);
  if (diff<60) std::cout <<"Lapse: "<<diff<<" secs"<<endl;
  else if (diff<3600) std::cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
  else std::cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
  std::cout<<RESET;
  return;
}


// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************

#ifdef _MASS_WEIGHT_POWER_
void FftwFunctions::grid_assignment(real_prec x,real_prec y,real_prec z,real_prec weight, real_prec weight_mass, vector<real_prec>& field, vector<real_prec>& field_mw)
#else
void FftwFunctions::grid_assignment(real_prec x,real_prec y,real_prec z,real_prec weight, vector<real_prec>& field)
#endif
{
  //////////////////////////////////////////////////////////
  // Interpolation of the galaxy overdensity field                    
  // We apply periodic bounday conditions to remap objects 
  // with coords. outside the range [0,Lside] within the box                           
  //////////////////////////////////////////////////////////  
  
  // ******************************************************
  // The second term places the mid point
  // in each component in the origin. The third
  // places the center of the sample (i.e, its x coord) in the center of the cube 
  // to place it inside the box, 
  x = x - this->Xoffset + 0.5*this->Lside;
  y = y - this->Yoffset + 0.5*this->Lside;
  z = z - this->Zoffset + 0.5*this->Lside;

  // ******************************************************
  // Boundary conditions:
  // ******************************************************

  while(x< 0){x+=this->Lside;}
  while(y< 0){y+=this->Lside;}
  while(z< 0){z+=this->Lside;}
  while(x>=this->Lside){x-=this->Lside;}
  while(y>=this->Lside){y-=this->Lside;}
  while(z>=this->Lside){z-=this->Lside;}

  
  // ***************************************************** 
  // ***************************************************** 
  // Determining the number of cell in each direction:   *
  // Output in the range [0,N-1]                         *
  // ***************************************************** 
  int xc=(int)floor((real_prec)(x*this->rdeltax));
  int yc=(int)floor((real_prec)(y*this->rdeltay));
  int zc=(int)floor((real_prec)(z*this->rdeltaz));


  // ******************************************************
  // This is just to be extremily sure that
  // we do not go beyond a cell with index greater
  // than the one we are setting. It should never happen
  // by construction.
  if(xc==this->Nft)xc=this->Nft-1;  
  if(yc==this->Nft)yc=this->Nft-1;
  if(zc==this->Nft)zc=this->Nft-1;
  // ******************************************************

  // ******************************************************
  // ******************************************************
  // Identify the cell to which this particle belongs     *
  // with its centre.                                     *     
  // ******************************************************
  real_prec xx  = this->deltax*((real_prec)xc+0.5);
  real_prec yy  = this->deltay*((real_prec)yc+0.5);
  real_prec zz  = this->deltaz*((real_prec)zc+0.5);
  
  // ******************************************************
  // ******************************************************
  // Center of each cell forwards                         *
  // ******************************************************
  real_prec xxf = this->deltax*((real_prec)xc+1.5);
  real_prec yyf = this->deltay*((real_prec)yc+1.5);
  real_prec zzf = this->deltaz*((real_prec)zc+1.5);


  // ******************************************************
  // ******************************************************
  // Center of each cell forwards twice                        *
  // ******************************************************
#ifdef _USE_FOURTH_ORDER_
  real_prec xxff = this->deltax*((real_prec)xc+2.5);
  real_prec yyff = this->deltay*((real_prec)yc+2.5);
  real_prec zzff = this->deltaz*((real_prec)zc+2.5);
#endif
  // ******************************************************
  // ******************************************************
  // Center of each cell backwards                        *
  // ******************************************************
  real_prec xxb = this->deltax*((real_prec)xc-0.5);
  real_prec yyb = this->deltay*((real_prec)yc-0.5);
  real_prec zzb = this->deltaz*((real_prec)zc-0.5);

  // ******************************************************
#ifdef _USE_FOURTH_ORDER_
  // ******************************************************
  // Center of each cell backwards twice                       *
  // ******************************************************
  real_prec xxbb = this->deltax*((real_prec)xc-1.5);
  real_prec yybb = this->deltay*((real_prec)yc-1.5);
  real_prec zzbb = this->deltaz*((real_prec)zc-1.5);
#endif

  // ******************************************************
  // ******************************************************
  // Particular positions (borders)
  int Xb=(xc==0 ? this->Nft: xc);
  int Yb=(yc==0 ? this->Nft: yc);
  int Zb=(zc==0 ? this->Nft: zc);
  
  int Xf=(xc==this->Nft-1 ? -1: xc);
  int Yf=(yc==this->Nft-1 ? -1: yc);
  int Zf=(zc==this->Nft-1 ? -1: zc);

#ifdef _USE_FOURTH_ORDER_
  int Xbb=(xc==0 ? this->Nft: xc); // if at 0, when subtracting 2, it goes to Nft-2
  int Ybb=(yc==0 ? this->Nft: yc);
  int Zbb=(zc==0 ? this->Nft: zc);

  if(xc!=0)
    Xbb=(xc==1 ? this->Nft+1: xc);  // if at 1, when subtracting 2, it goes to Nft-1
  if(yc!=0)
   Ybb=(yc==1 ? this->Nft+1: yc);
  if(zc!=0)
   Zbb=(zc==1 ? this->Nft+1: zc);

  int Xff=(xc==this->Nft-1 ? -1: xc);  //if at Nft-1, when adding 2, it goes to 1
  int Yff=(yc==this->Nft-1 ? -1: yc);
  int Zff=(zc==this->Nft-1 ? -1: zc);

  if(xc!=this->Nft-1)
     Xff=(xc==this->Nft-2 ? -2: xc); // if at Nft-2, if adding 2, it goes to zero
  if(yc!=this->Nft-1)
    Yff=(yc==this->Nft-2 ? -2: yc);
  if(zc!=this->Nft-1)
    Zff=(zc==this->Nft-2 ? -2: zc);
#endif

  // ******************************************************
  // ******************************************************
  // Some controls
  // ******************************************************
#ifdef _CHECK_PARTICLES_IN_BOX_

  if(Zb-1>Nft-1 || zc >Nft-1 || Zf+1>Nft-1)
    {
#pragma omp critical
      {
        cout<<"z"<<"  "<<z<<"  "<<Lside<<" " << z/Lside << endl;
        cout<< Zb-1 << " " << zc << " " << Zf+1 << " [" << Nft-1 << "]" << endl;
      }
    }
  
  if(Yb-1>Nft-1 || yc >Nft-1 || Yf+1>Nft-1)
    {
#pragma omp critical
      {
        cout<<"y"<<"  "<<y<<"  "<<Lside<<" " << y/Lside << endl;
        cout<< Yb-1 << " " << yc << " " << Yf+1 << " [" << Nft-1 << "]" << endl;
      }
    }

  if(Xb-1>Nft-1 || xc >Nft-1 || Xf+1>Nft-1)
    {
#pragma omp critical
      {
	cout<<"x"<<"  "<<x<<"  "<<Lside<<" " << x/Lside <<"  "<<this->Xoffset<<endl;
	cout<< Xb-1 << " " << xc << " " << Xf+1 << " [" << Nft-1 << "]" << endl;
      }
    }
  // ******************************************************
  // ******************************************************
#endif



#ifdef _USE_FOURTH_ORDER_
  int i_idx[MAX_MAS_DEG] = {Xbb-2, Xb-1, xc, Xf+1, Xff+2};
  int j_idx[MAX_MAS_DEG] = {Ybb-2, Yb-1, yc, Yf+1, Xff+2};
  int k_idx[MAX_MAS_DEG] = {Zbb-2, Zb-1, zc, Zf+1, Zff+2};
  real_prec MAS_xx[MAX_MAS_DEG] =
    {
      (this->*MAS_ptr)((xxbb - x)*rdeltax),
      (this->*MAS_ptr)((xxb  - x)*rdeltax),
      (this->*MAS_ptr)((xx   - x)*rdeltax),
      (this->*MAS_ptr)((xxf  - x)*rdeltax),
      (this->*MAS_ptr)((xxff - x)*rdeltax)
  };
  real_prec MAS_yy[MAX_MAS_DEG] =
    {
      (this->*MAS_ptr)((yybb - y)*rdeltay),
      (this->*MAS_ptr)((yyb  - y)*rdeltay),
      (this->*MAS_ptr)((yy   - y)*rdeltay),
      (this->*MAS_ptr)((yyf  - y)*rdeltay),
      (this->*MAS_ptr)((yyff - y)*rdeltay)
    };
  real_prec MAS_zz[MAX_MAS_DEG] =
    {
      (this->*MAS_ptr)((zzbb - z)*rdeltaz),
      (this->*MAS_ptr)((zzb  - z)*rdeltaz),
      (this->*MAS_ptr)((zz   - z)*rdeltaz),
      (this->*MAS_ptr)((zzf  - z)*rdeltaz),
      (this->*MAS_ptr)((zzff - z)*rdeltaz)
    };

#else

  vector<int>i_idx(MAX_MAS_DEG,0);
  i_idx[0]=(Xb-1);
  i_idx[1]=(xc);
  i_idx[2]=(Xf+1);

  vector<int>j_idx(MAX_MAS_DEG,0);
  j_idx[0]=(Yb-1);
  j_idx[1]=(yc);
  j_idx[2]=(Yf+1);

  vector<int>k_idx(MAX_MAS_DEG,0);
  k_idx[0]=(Zb-1);
  k_idx[1]=(zc);
  k_idx[2]=(Zf+1);


  vector<real_prec>MAS_xx;
  MAS_xx.push_back((this->*MAS_ptr)((xxb  - x)*rdeltax));
  MAS_xx.push_back((this->*MAS_ptr)((xx   - x)*rdeltax));
  MAS_xx.push_back((this->*MAS_ptr)((xxf  - x)*rdeltax));

  vector<real_prec>MAS_yy;
  MAS_yy.push_back((this->*MAS_ptr)((yyb  - y)*rdeltay));
  MAS_yy.push_back((this->*MAS_ptr)((yy   - y)*rdeltay));
  MAS_yy.push_back((this->*MAS_ptr)((yyf  - y)*rdeltay));

  vector<real_prec>MAS_zz;
  MAS_zz.push_back((this->*MAS_ptr)((zzb  - z)*rdeltaz));
  MAS_zz.push_back((this->*MAS_ptr)((zz   - z)*rdeltaz));
  MAS_zz.push_back((this->*MAS_ptr)((zzf  - z)*rdeltaz));
#endif

  for(int ih=0;ih<MAX_MAS_DEG;++ih)
    for(int jh=0;jh<MAX_MAS_DEG;++jh)
      for(int kh=0;kh<MAX_MAS_DEG;++kh)
        {
          real_prec weight_aux=weight*MAS_xx[ih]*MAS_yy[jh]*MAS_zz[kh];
          ULONG index_aux = index_3d(i_idx[ih], j_idx[jh], k_idx[kh],this->Nft,this->Nft);
          field[index_aux] +=weight_aux;
#ifdef _MASS_WEIGHT_POWER_
          field_mw[index_aux] +=weight_aux*weight_mass;
#endif
        }

}





// *************************************************************************************************
// *************************************************************************************************
// *************************************************************************************************
// *************************************************************************************************
#ifdef _USE_VECTORIZED_GRID_ASSIGNMENT_
void FftwFunctions::grid_assignment_old(  real_prec x,
                                      real_prec y,
                                      real_prec z,
                                      real_prec *weights,
                                      pic_storage **data,
                                      int offset, int N_storages
                                      )
{
  //////////////////////////////////////////////////////////
  // Interpolation of the galaxy overdensity field
  // We apply periodic bounday conditions to remap objects
  // with coords. outside the range [0,Lside] within the box
  //////////////////////////////////////////////////////////

  // ******************************************************
  // Part suggested by Florian code:
  // The second term places the mid point
  // in each component in the origin. The third
  // places the center of the sample (i.e, its x coord) in the center of the cube
  // to place it inside the box,
  x = x - this->Xoffset + 0.5*this->Lside;
  y = y - this->Yoffset + 0.5*this->Lside;
  z = z - this->Zoffset + 0.5*this->Lside;

  // ******************************************************
  // Boundary conditions:
  // ******************************************************

  while(x< 0){x+=this->Lside;}
  while(y< 0){y+=this->Lside;}
  while(z< 0){z+=this->Lside;}
  while(x>=this->Lside){x-=this->Lside;}
  while(y>=this->Lside){y-=this->Lside;}
  while(z>=this->Lside){z-=this->Lside;}


  // *****************************************************
  // *****************************************************
  // Determining the number of cell in each direction:   *
  // Output in the range [0,N-1]                         *
  // *****************************************************
  int xc=(int)floor((real_prec)(x*this->rdeltax));
  int yc=(int)floor((real_prec)(y*this->rdeltay));
  int zc=(int)floor((real_prec)(z*this->rdeltaz));


  // ******************************************************
  // This is just to be extremily sure that
  // we do not go beyond a cell with index greater
  // than the one we are setting. It should never happen
  // by construction.
  if(xc==this->Nft)xc=this->Nft-1;
  if(yc==this->Nft)yc=this->Nft-1;
  if(zc==this->Nft)zc=this->Nft-1;
  // ******************************************************

  // ******************************************************
  // ******************************************************
  // Identify the cell to which this particle belongs     *
  // with its centre.                                     *
  // ******************************************************
  real_prec xx  = this->deltax*((real_prec)xc+0.5);
  real_prec yy  = this->deltay*((real_prec)yc+0.5);
  real_prec zz  = this->deltaz*((real_prec)zc+0.5);

  // ******************************************************
  // ******************************************************
  // Center of each cell forwards                         *
  // ******************************************************
  real_prec xxf = this->deltax*((real_prec)xc+1.5);
  real_prec yyf = this->deltay*((real_prec)yc+1.5);
  real_prec zzf = this->deltaz*((real_prec)zc+1.5);


  // ******************************************************
  // ******************************************************
  // Center of each cell forwards twice                        *
  // ******************************************************

  real_prec xxff = this->deltax*((real_prec)xc+2.5);
  real_prec yyff = this->deltay*((real_prec)yc+2.5);
  real_prec zzff = this->deltaz*((real_prec)zc+2.5);

  // ******************************************************
  // ******************************************************
  // Center of each cell backwards                        *
  // ******************************************************
  real_prec xxb = this->deltax*((real_prec)xc-0.5);
  real_prec yyb = this->deltay*((real_prec)yc-0.5);
  real_prec zzb = this->deltaz*((real_prec)zc-0.5);

  // ******************************************************
  // ******************************************************
  // Center of each cell backwards twice                       *
  // ******************************************************
  real_prec xxbb = this->deltax*((real_prec)xc-1.5);
  real_prec yybb = this->deltay*((real_prec)yc-1.5);
  real_prec zzbb = this->deltaz*((real_prec)zc-1.5);


  // ******************************************************
  // ******************************************************
  // Particular positions (borders)
  int Xb=(xc==0 ? this->Nft: xc);
  int Yb=(yc==0 ? this->Nft: yc);
  int Zb=(zc==0 ? this->Nft: zc);

  int Xf=(xc==this->Nft-1 ? -1: xc);
  int Yf=(yc==this->Nft-1 ? -1: yc);
  int Zf=(zc==this->Nft-1 ? -1: zc);

  // ******************************************************
  // ******************************************************
  // Some controls
  // ******************************************************

  if(Zb-1>Nft-1 || zc >Nft-1 || Zf+1>Nft-1)
    {
#pragma omp critical
      {
        cout<<"z"<<"  "<<z<<"  "<<Lside<<" " << z/Lside << endl;
        cout<< Zb-1 << " " << zc << " " << Zf+1 << " [" << Nft-1 << "]" << endl;
      }
    }

  if(Yb-1>Nft-1 || yc >Nft-1 || Yf+1>Nft-1)
    {
#pragma omp critical
      {
        cout<<"y"<<"  "<<y<<"  "<<Lside<<" " << y/Lside << endl;
        cout<< Yb-1 << " " << yc << " " << Yf+1 << " [" << Nft-1 << "]" << endl;
      }
    }

  if(Xb-1>Nft-1 || xc >Nft-1 || Xf+1>Nft-1)
    {
#pragma omp critical
      {
        cout<<"x"<<"  "<<x<<"  "<<Lside<<" " << x/Lside <<"  "<<this->Xoffset<<endl;
        cout<< Xb-1 << " " << xc << " " << Xf+1 << " [" << Nft-1 << "]" << endl;
      }
    }
  // ******************************************************
  // ******************************************************

  int i_idx[MAX_MAS_DEG] = {Xb-2, Xb-1, xc, Xf+1, Xf+2};
  int j_idx[MAX_MAS_DEG] = {Yb-2, Yb-1, yc, Yf+1, Xf+2};
  int k_idx[MAX_MAS_DEG] = {Zb-2, Zb-1, zc, Zf+1, Zf+2};


  int grid_indexes[CHUNK];

  real_prec MAS_xx[MAX_MAS_DEG] =
    {
      (this->*MAS_ptr)((xxbb - x)*rdeltax),
      (this->*MAS_ptr)((xxb  - x)*rdeltax),
      (this->*MAS_ptr)((xx   - x)*rdeltax),
      (this->*MAS_ptr)((xxf  - x)*rdeltax),
      (this->*MAS_ptr)((xxff - x)*rdeltax)
    };
  real_prec MAS_yy[MAX_MAS_DEG] =
    {
      (this->*MAS_ptr)((yybb - y)*rdeltay),
      (this->*MAS_ptr)((yyb  - y)*rdeltay),
      (this->*MAS_ptr)((yy   - y)*rdeltay),
      (this->*MAS_ptr)((yyf  - y)*rdeltay),
      (this->*MAS_ptr)((yyff - y)*rdeltay)
    };
  real_prec MAS_zz[MAX_MAS_DEG] =
    {
      (this->*MAS_ptr)((zzbb - z)*rdeltaz),
      (this->*MAS_ptr)((zzb  - z)*rdeltaz),
      (this->*MAS_ptr)((zz   - z)*rdeltaz),
      (this->*MAS_ptr)((zzf  - z)*rdeltaz),
      (this->*MAS_ptr)((zzff - z)*rdeltaz)
    };

  real_prec grid_weights[CHUNK];



  // /////////////////////////////////////////////////////////////////
  for(int iii = 0, q = 0; iii < MAX_MAS_DEG; iii++)
    {
      int I = i_idx[iii] * this->Nft2;
#pragma ivdep
      for(int jjj = 0; jjj < MAX_MAS_DEG; jjj++)
        {
          int J = j_idx[jjj] * this->Nft;
          grid_indexes[q]   = I + J + k_idx[0];
          grid_indexes[q+1] = I + J + k_idx[1];
          grid_indexes[q+2] = I + J + k_idx[2];
          grid_indexes[q+3] = I + J + k_idx[3];
          grid_indexes[q+4] = I + J + k_idx[4];
          grid_indexes[q+5] = I + J + k_idx[5];
          q += 5;
        }
    }


  for(int iii = 0, q = 0; iii < MAX_MAS_DEG; iii++) {
    real_prec wi = MAS_xx[iii];
#pragma ivdep
    for(int jjj = 0; jjj < MAX_MAS_DEG; jjj++) {
      real_prec wj = wi * MAS_yy[jjj];
      grid_weights[q]   = wj * MAS_zz[0];
      grid_weights[q+1] = wj * MAS_zz[1];
      grid_weights[q+2] = wj * MAS_zz[2];
      grid_weights[q+3] = wj * MAS_zz[3];
      grid_weights[q+4] = wj * MAS_zz[4];
      grid_weights[q+5] = wj * MAS_zz[5];
      q += 5;
    }
  }
  // /////////////////////////////////////////////////////////////////

  // ******************************************************
  // ******************************************************
  // Filling grid cells                                   *
  // ******************************************************

  pic_storage *storage;

  for(int sss = 0; sss < N_storages; sss++)
    {
      storage = data[sss] + offset;
#pragma ivdep
      for(int iii = 0; iii < CHUNK; iii++)
       {
          storage[iii].weight = weights[sss] * grid_weights[iii];
          storage[iii].idx = grid_indexes[iii];
       }
   }



}
#endif


// ***************************************************************************************
// ***************************************************************************************
// ***************************************************************************************
// ***************************************************************************************
// ***************************************************************************************

real_prec FftwFunctions::correction_MAS(int i,
				     int j, 
				     int k
				     )
{
  //////////////////////////////////////////////////////////  
  // Correcting for the MAS mode-by-mode
  //////////////////////////////////////////////////////////  
  real_prec xx=i*M_PI/this->Nft;
  real_prec yy=j*M_PI/this->Nft;
  real_prec zz=k*M_PI/this->Nft;
  real_prec sx=(i==0? 1.00: sin(xx)/xx);
  real_prec sy=(j==0? 1.00: sin(yy)/yy);
  real_prec sz=(k==0? 1.00: sin(zz)/zz);
  
  real_prec res = sx*sy*sz;
  real_prec ans;
  
  switch((int)correction_MAS_exp)
    {
    case(0): ans=res; break;
    case(1): ans=res*res; break;
    case(2): ans=pow(res,3); break;
    case(3): ans=pow(res,4); break;
    }
  
  
  // for(int ii = 0; ii < correction_MAS_exp; ii++)
  //   res *= res;
    
  return ans;
}
// *************************************************************************************************
// *************************************************************************************************
// *************************************************************************************************
inline real_prec SN_correction_MAS(int exp,real_prec xx,real_prec yy,real_prec zz)
{
  //////////////////////////////////////////////////////////  
  // Correcting for the MAS mode-by-mode
  //////////////////////////////////////////////////////////  
  real_prec ans;
  
  switch(exp)
    {
    case(0):
      ans=1; break;
    case(1):
      ans=(1.0-(2./3.)*xx*xx)*(1.0-(2./3.)*yy*yy)*(1.0-(2./3.)*zz*zz); break;
    case(2):
      ans=(1.0-xx*xx+(2./15.)*xx*xx*xx*xx)*(1.0-yy*yy+(2./15.)*yy*yy*yy*yy)*(1.0-zz*zz+(2./15.)*zz*zz*zz*zz)  ; break;
    case(3):
      ans=1; break;  // This is one becuase it has not been yet computed
    }
  return ans;
}
// *************************************************************************************************
// *************************************************************************************************
// *************************************************************************************************


inline real_prec my_correction_MAS(int exp,
                                real_prec xx,
                                real_prec yy,
                                real_prec zz
                                )
{
  //////////////////////////////////////////////////////////
  // Correcting for the MAS mode-by-mode
  //////////////////////////////////////////////////////////

  real_prec res = xx*yy*zz;
  real_prec ans;

  switch(exp)
    {
    case(0): ans=res; break;
    case(1): ans=res*res; break;
    case(2): ans=(res*res)*res; break;
    case(3): ans=(res*res)*res*res; break;
    }

  return ans;
}

// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************


void FftwFunctions::get_interpolated_density_field(s_parameters_box *s_box, 
				   s_data_structure *s_data)
{ 
  /////////////////////////////////////////////////////////  
  // Sampling of the random or real catalog and 
  // interpolation of the object density field into a grid.
  // This method expects the coordinates of the input 
  // catalogue already transformed into cartessian, 
  // with the information of the mean number density written
  // in the corresponding slot as that specified in the
  // parameter file.
  //////////////////////////////////////////////////////////  

  int NTHREADS = 1;//omp_get_max_threads();
//  omp_set_num_threads(NTHREADS);


  // *****************************************************************************
  // Identify columns in the corresponding catalog                               *
  // *****************************************************************************
 // int i_coord1      = (s_data->catalog=="data"? s_box->i_coord1_g : s_box->i_coord1_r);
 // int i_coord2      = (s_data->catalog=="data"? s_box->i_coord2_g : s_box->i_coord2_r);
//  int i_coord3      = (s_data->catalog=="data"? s_box->i_coord3_g : s_box->i_coord3_r);
#ifdef _USE_WEIGHTS_IN_POWER_
  int i_weight1     = (s_data->catalog=="data"? s_box->i_weight1_g : s_box->i_weight1_r);
  int i_weight2     = (s_data->catalog=="data"? s_box->i_weight2_g : s_box->i_weight2_r);
  int i_weight3     = (s_data->catalog=="data"? s_box->i_weight3_g : s_box->i_weight3_r);
  int i_weight4     = (s_data->catalog=="data"? s_box->i_weight4_g : s_box->i_weight4_r);
  bool use_weight1  = (s_data->catalog=="data"? s_box->use_weight1_g : s_box->use_weight1_r);
  bool use_weight2  = (s_data->catalog=="data"? s_box->use_weight2_g : s_box->use_weight2_r);
  bool use_weight3  = (s_data->catalog=="data"? s_box->use_weight3_g : s_box->use_weight3_r);
  bool use_weight4  = (s_data->catalog=="data"? s_box->use_weight4_g : s_box->use_weight4_r);
#endif

  //  vector < real_prec >prop = s_data->properties;
  //int i_mean_density= (s_data->catalog=="data"? s_box->i_mean_density_g : s_box->i_mean_density_r);

  ULONG nlines= s_data->properties.size();//  prop prop.size() / n_columns;


  if(this->MASS=="NGP")
    {
      this->MAS_ptr = &FftwFunctions::MAS_NGP;
      correction_MAS_exp = 0;
    }
  else if(this->MASS=="CIC")
    {
      this->MAS_ptr = &FftwFunctions::MAS_CIC;
      correction_MAS_exp = 1;
    }
  else if(this->MASS=="TSC")
    {
      this->MAS_ptr = &FftwFunctions::MAS_TSC;
      correction_MAS_exp = 2;
    }
  else if(this->MASS=="PCS")
    {
      this->MAS_ptr = &FftwFunctions::MAS_PCS;
      correction_MAS_exp = 3;
    }



  // *****************************************************************************
  // Initialize variables used in the estimation of Pk / Bk
  int n_selected=0;
  real_prec S_r_power=0;
  real_prec normal_power=0;
  real_prec W_r=0;
  real_prec S_r1=0;
  real_prec S_r2=0;
  real_prec normal_bispectrum=0;


 vector<real_prec>field(this->NTT,0);

#ifdef _MASS_WEIGHT_POWER_
 vector<real_prec>field_mw(this->NTT,0);
#endif

#ifdef _USE_WEIGHTS_IN_POWER_
    vector<real_prec> ow(MAX_NUMBER_WEIGHTS,0);
#endif


    for (ULONG i=0;i< nlines ;++i)
      {
        real_prec x=s_data->properties[i].coord1;
        real_prec y=s_data->properties[i].coord2;
        real_prec z=s_data->properties[i].coord3;
#ifdef _USE_MASS_AS_OBSERVABLE_
        real_prec property=s_data->properties[i].mass;
#elif defined _USE_VMAX_AS_OBSERVABLE_
        real_prec property=s_data->properties[i].vmax;
#endif

#ifdef _SET_GLOBAL_MASS_CUT_
        if(property>=MINIMUM_PROP_CUT)
            {
#endif

#ifdef _MASS_WEIGHT_POWER_
        real_prec mass=s_data->properties[i].mass;
#endif
        real_prec nbar=s_data->mean_density;
	if(true==s_box->use_random_catalog)
          nbar=s_data->properties[i].mean_density;


        real_prec ptotal_weight = 1.0;
#ifdef _USE_WEIGHTS_IN_POWER_
        ow[0] = (use_weight1 && (i_weight1<n_columns))? s_data->properties[i].weight1 : 1.0;
        ow[1] = (use_weight2 && (i_weight2<n_columns))? s_data->properties[i].weight2 : 1.0;
        ow[2] = (use_weight3 && (i_weight3<n_columns))? s_data->properties[i].weight3 : 1.0;
        ow[3] = (use_weight4 && (i_weight4<n_columns))? s_data->properties[i].weight4 : 1.0;
        ptotal_weight = ow[0] * ow[1] * ow[2] * ow[3];
#endif
        real_prec we_fkp=1.0;
        if(true==s_box->FKP_weight)
           we_fkp=1.0/(1.0+s_box->Pest*nbar);

        ptotal_weight*=we_fkp;
	
	// ***********************************************************
	// Parameters for the Pk
	n_selected++;                                //number of selected objects
	W_r+=ptotal_weight;                          //weighted number of selected objects
	real_prec ptotal_weight2 = ptotal_weight * ptotal_weight;
	S_r_power+= ptotal_weight2;                   //sum of the squared of weight
	normal_power += nbar*ptotal_weight2;    //normalization of power spectrum
	
	// ***********************************************************
	// Parameters for the bispectrum: shot noise
#ifdef _GET_BISPECTRUM_NUMBERS_
        real_prec ptotal_weight3 = ptotal_weight2 * ptotal_weight;
	S_r1 += ptotal_weight3*nbar;           
	S_r2 += ptotal_weight3;                
	normal_bispectrum+= nbar * nbar * ptotal_weight3;
#endif
         // *****************************
	// Interpolation of the object density field:
#ifdef _MASS_WEIGHT_POWER_
        this->grid_assignment(x, y, z, ptotal_weight, mass, field, field_mw);
#else
        this->grid_assignment(x, y, z, ptotal_weight, field);
#endif

#ifdef _SET_GLOBAL_MASS_CUT_
            }
#endif

    }
  


  // When a random catalog is used, the
  // quantities below, computed from the random, 
  // are used to get the normalization and the shot-noise.
  // If a simulation with no random is used, 
  // the values below are computed in the function
  // raw sampling.
  
  if(s_data->catalog=="data")
    {
      this->data_g=field;
#ifdef _MASS_WEIGHT_POWER_
      this->data_g_mw=field_mw;
#endif
      this->n_gal=n_selected;
      this->w_g=W_r;
      this->normal_p=1.0; 
      this->normal_b=1.0;
      this->s_g=S_r_power;
      this->sg1=S_r1;
      this->sg2=S_r2;
    }
  else if(s_data->catalog=="random")
   {
      this->data_r=field;
      this->n_ran=n_selected;
      this->w_r=W_r;
      this->s_r=S_r_power;
      this->normal_p=normal_power;
      this->normal_b=normal_bispectrum;
      this->sr1=S_r1;
      this->sr2=S_r2;
    }

  field.clear();
  field.shrink_to_fit();
#ifdef _MASS_WEIGHT_POWER_
  field_mw.clear();
  field_mw.shrink_to_fit();
#endif


}

// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// *******************************************************************************************************
#ifdef _USE_VECTORIZED_GRID_ASSIGNMENT_
// This is the function written by L. Tornatore, using vectorization. It had implemented the "wrong" approach
// to te yamamoto estimator, and had some memmory problems that I could not find. I moved to the original version and dubbed this as "old"
void FftwFunctions::get_interpolated_density_field_old(s_parameters_box *s_box,
                                   s_data_structure *s_data)
{
  /////////////////////////////////////////////////////////
  // Sampling of the random or real catalog and
  // interpolation of the object density field into a grid.
  // This method expects the coordinates of the input
  // catalogue already transformed into cartessian,
  // with the information of the mean number density written
  // in the corresponding slot as that specified in the
  // parameter file.
  //////////////////////////////////////////////////////////

  // version en paralelo

#define Pk_ys 1
#define Pk_yb 2
#define Pk_ybc -1
#define Pk_ysc -2

  cout.precision(12);

  // *****************************************************************************
  // Identify columns in the corresponding catalog                               *
  // *****************************************************************************
  int i_coord1      = (s_data->catalog=="data"? s_box->i_coord1_g : s_box->i_coord1_r);
  int i_coord2      = (s_data->catalog=="data"? s_box->i_coord2_g : s_box->i_coord2_r);
  int i_coord3      = (s_data->catalog=="data"? s_box->i_coord3_g : s_box->i_coord3_r);
  int i_weight1     = (s_data->catalog=="data"? s_box->i_weight1_g : s_box->i_weight1_r);
  int i_weight2     = (s_data->catalog=="data"? s_box->i_weight2_g : s_box->i_weight2_r);
  int i_weight3     = (s_data->catalog=="data"? s_box->i_weight3_g : s_box->i_weight3_r);
  int i_weight4     = (s_data->catalog=="data"? s_box->i_weight4_g : s_box->i_weight4_r);
  bool use_weight1  = (s_data->catalog=="data"? s_box->use_weight1_g : s_box->use_weight1_r);
  bool use_weight2  = (s_data->catalog=="data"? s_box->use_weight2_g : s_box->use_weight2_r);
  bool use_weight3  = (s_data->catalog=="data"? s_box->use_weight3_g : s_box->use_weight3_r);
  bool use_weight4  = (s_data->catalog=="data"? s_box->use_weight4_g : s_box->use_weight4_r);

  vector < real_prec >prop = s_data->properties;
  int i_mean_density= (s_data->catalog=="data"? s_box->i_mean_density_g : s_box->i_mean_density_r);
  int n_columns= s_data->n_columns;
  ULONG nlines= prop.size() / n_columns;

  // *****************************************************************************
  // Initialize variables used in the estimation of Pk / Bk
  int n_selected=0;
  real_prec S_r_power=0;
  real_prec normal_power=0;
  real_prec W_r=0;
  real_prec S_r1=0;
  real_prec S_r2=0;
  real_prec normal_bispectrum=0;
  // *****************************************************************************
  int my_statistics;  // avoid string comparison in the cycle
  // *****************************************************************************
  // cout<<RED<<"Warning in get_interpolated(); using >1 threads give unaccurate results"<<endl;
  // cout<<"Set NTHREADS=1"<<RESET<<endl;
  int NTHREADS = 1;//omp_get_max_threads();
  // *****************************************************************************

  void **pdata;  // use void* instead of vector<real_prec>* to avoid use of [] operator and memcpy
  int pdata_nidx;
  pic_storage *wstorage[22];


  pdata_nidx = 1;
  if(statistics=="Pk_ys")
    {
      my_statistics = Pk_ys;
      pdata_nidx += 6;
    }
  else if(statistics == "Pk_yb")
    {
      my_statistics = Pk_yb;
      pdata_nidx += 21;
    }
  else
    my_statistics = 0;  // This contains FKP and the right implmentation of Yamamoto estimator.

  for(int ii = 0; ii < pdata_nidx; ii++)
    wstorage[ii] = new pic_storage [NTHREADS * CHUNK];

  pdata = new void* [pdata_nidx];
  if(s_data->catalog=="data")
    pdata[0] = (void*)&(this->data_g);
  else if(s_data->catalog=="random")
    pdata[0] = (void*)&(this->data_r);

  if(my_statistics > 0)
    {  //This is only done for Yamamoto-Sco, the wrong way
      if(s_data->catalog=="data")
        {
          pdata[1] = (void*)&(this->data_g_xx); // x*x
          pdata[4] = (void*)&(this->data_g_yy); // y*y
          pdata[6] = (void*)&(this->data_g_zz); // z*z
          pdata[2] = (void*)&(this->data_g_xy); // x*y
          pdata[3] = (void*)&(this->data_g_xz); // x*z
          pdata[5] = (void*)&(this->data_g_yz); // y*z
        }
      else if(s_data->catalog=="random")
        {
          pdata[1] = (void*)&(this->data_r_xx);
          pdata[4] = (void*)&(this->data_r_yy);
          pdata[6] = (void*)&(this->data_r_zz);
          pdata[2] = (void*)&(this->data_r_xy);
          pdata[3] = (void*)&(this->data_r_xz);
          pdata[5] = (void*)&(this->data_r_yz);
        }

      if(my_statistics == Pk_yb)
        {  //This is only done for Yamamoto-Bianchi, the wrong way
          if(s_data->catalog=="data")
            {
              pdata[7 ] = (void*)&(this->data_g_xxx); // xxxx
              pdata[17] = (void*)&(this->data_g_yyy); // yyyy
              pdata[21] = (void*)&(this->data_g_zzz); // zzzz
              pdata[8 ] = (void*)&(this->data_g_xxy); // xxxy
              pdata[9 ] = (void*)&(this->data_g_xxz); // xxxz
              pdata[13] = (void*)&(this->data_g_yyx); // yyyx
              pdata[18] = (void*)&(this->data_g_yyz); // yyyz
              pdata[16] = (void*)&(this->data_g_zzx); // zzzx
              pdata[20] = (void*)&(this->data_g_zzy); // zzzy
              pdata[10] = (void*)&(this->data_g_xyy); // xxyy
              pdata[12] = (void*)&(this->data_g_xzz); // xxzz
              pdata[19] = (void*)&(this->data_g_yzz); // yyzz
              pdata[11] = (void*)&(this->data_g_xyz); // xxyz
              pdata[14] = (void*)&(this->data_g_yxz); // yyxz
              pdata[15] = (void*)&(this->data_g_zxy); // zzxy
            }
          else if(s_data->catalog=="random")
            {
              pdata[7 ] = (void*)&(this->data_r_xxx);
              pdata[17] = (void*)&(this->data_r_yyy);
              pdata[21] = (void*)&(this->data_r_zzz);
              pdata[8 ] = (void*)&(this->data_r_xxy);
              pdata[9 ] = (void*)&(this->data_r_xxz);
              pdata[13] = (void*)&(this->data_r_yyx);
              pdata[18] = (void*)&(this->data_r_yyz);
              pdata[16] = (void*)&(this->data_r_zzx);
              pdata[20] = (void*)&(this->data_r_zzy);
              pdata[10] = (void*)&(this->data_r_xyy);
              pdata[12] = (void*)&(this->data_r_xzz);
              pdata[19] = (void*)&(this->data_r_yzz);
              pdata[11] = (void*)&(this->data_r_xyz);
              pdata[14] = (void*)&(this->data_r_yxz);
              pdata[15] = (void*)&(this->data_r_zxy);
            }
        }
    }

  if(this->MASS=="NGP")
    {
      MAS_ptr = &FftwFunctions::MAS_NGP;
      correction_MAS_exp = 0;
    }
  else if(this->MASS=="CIC")
    {
      MAS_ptr = &FftwFunctions::MAS_CIC;
      correction_MAS_exp = 1;
    }
  else if(this->MASS=="TSC")
    {
      MAS_ptr = &FftwFunctions::MAS_TSC;
      correction_MAS_exp = 2;
    }
  else if(this->MASS=="PCS")
    {
      MAS_ptr = &FftwFunctions::MAS_PCS;
      correction_MAS_exp = 3;
    }
  time_t start;
  time (&start);

  int alldata = NTHREADS*CHUNK;
  int ALLSTOP = 0;




#pragma omp parallel num_threads(NTHREADS) reduction(+:n_selected,normal_power,W_r,S_r_power,S_r1,S_r2,normal_bispectrum)
  {

    int myID = omp_get_thread_num();
    real_prec ow[MAX_NUMBER_WEIGHTS];

#pragma omp single
    So.message_screen("Using", NTHREADS, " threads");

    real_prec ptotal_weight=1.0;
    real_prec *pprop;
    real_prec *weights;

    weights = new real_prec[pdata_nidx];
    int GO = 1;
    int i = myID;

    while(ALLSTOP < NTHREADS)
      {
        if(GO)
          {

            pprop = &prop[i*n_columns];

            real_prec we;
            real_prec x=pprop[i_coord1];
            real_prec y=pprop[i_coord2];
            real_prec z=pprop[i_coord3];

            real_prec nbar=s_data->mean_density;
            if(true==s_box->use_random_catalog)
              nbar=pprop[i_mean_density];

            // distance squared, used in the wrong Yamamoto-related interpolations
            // and to be discarted once the correct one is implemented succesfully.
            real_prec rr=(x*x+y*y+z*z);
            real_prec irr = 1.0 / rr;
            real_prec r4=rr*rr;
            real_prec ir4 = 1.0 / r4;
            // **************************************************************************
            // Compute the weights for each galaxy

            ow[0] = (use_weight1 && (i_weight1<n_columns))? pprop[i_weight1] : 1.0;
            ow[1] = (use_weight2 && (i_weight2<n_columns))? pprop[i_weight2] : 1.0;
            ow[2] = (use_weight3 && (i_weight3<n_columns))? pprop[i_weight3] : 1.0;
            ow[3] = (use_weight4 && (i_weight4<n_columns))? pprop[i_weight4] : 1.0;

            ptotal_weight = ow[0] * ow[1] * ow[2] * ow[3];
            we=1.0;
            if(s_box->FKP_weight)
             {
               we=1.0/(1.0+s_box->Pest*nbar);
             }
            ptotal_weight*=we;

            // ***********************************************************
            // Parameters for the Pk
            n_selected++;                                //number of selected objects
            W_r+=ptotal_weight;                          //weighted number of selected objects
            real_prec ptotal_weight2 = ptotal_weight * ptotal_weight;
            S_r_power+= ptotal_weight2;                   //sum of the squared of weight
            normal_power += nbar*ptotal_weight2;    //normalization of power spectrum

            // ***********************************************************
            // Parameters for the bispectrum: shot noise
            real_prec ptotal_weight3 = ptotal_weight2 * ptotal_weight;
            S_r1 += ptotal_weight3*nbar;
            S_r2 += ptotal_weight3;
            normal_bispectrum+= nbar * nbar * ptotal_weight3;

            // *****************************
            // Interpolation of the object density field:

            weights[0] = ptotal_weight;

            // to be discarted once the correct one is implemented succesfully.
            real_prec carray[3] = {x, y, z};

            if(my_statistics > 0)
              { //This is Yamamoto, the wrong way (wither Scoc or Bianchi
                int q = 1;
                real_prec weight = ptotal_weight * irr;

                for(int iii = 0; iii < 3; iii++)
                  for(int jjj = iii; jjj < 3; jjj++)
                    weights[q++] = weight * carray[iii] * carray[jjj];

                if(my_statistics == Pk_yb)
                  {
                    weight = ptotal_weight * ir4;

                    for(int iii = 0; iii < 3; iii++)
                      {
                        real_prec wi = carray[iii];
                        for(int jjj = iii; jjj < 3; jjj++)
                          {
                            real_prec wj = wi * carray[jjj];
                            for(int kkk = jjj; kkk < 3; kkk++)
                              {
                                real_prec wk = wj * carray[kkk];
                                for(int lll = kkk; lll < 3; lll++)
                                  {
                                    weights[q++] = weight * wk * carray[lll];
                                  }
                              }
                          }
                      }
                  }
              }
            this->grid_assignment_old(x, y, z, weights, wstorage, myID * CHUNK, pdata_nidx);

          }
#pragma omp barrier

        // WARNING:
        // there is a problem here with large
        // valus of Nft. Found by Jennifer Pollack when running the code for the Bispectrum.
        for(int jj = myID; jj < pdata_nidx; jj+=NTHREADS)
          {
            // each task deal with a different grid, so we do not need to worry about
            // memory concurrency
            real_prec *pointer = *(real_prec**)pdata[jj];
            int last;

            last = (partition(wstorage[jj], wstorage[jj]+alldata, [](const pic_storage &A) -> bool {return (A.weight !=0);}) - wstorage[jj]);

            for(int ii = 0; ii < last; ii++)
              pointer[wstorage[jj][ii].idx] += wstorage[jj][ii].weight;

          }


        if(GO)
          {
            i += NTHREADS;
            if(i >= nlines)
              {
                GO = 0;
#pragma omp atomic update
                ALLSTOP += 1;
              }
          }



#pragma omp barrier
      }
  }  // END parallelized region.


  // When a random catalog is used, the
  // quantities below, computed from the random,
  // are used to get the normalization and the shot-nouise.
  // If a simulation with no random is used,
  // the values below are computed in the function
  // raw sampling.

  if(s_data->catalog=="data")
    {
      this->n_gal=n_selected;
      this->w_g=W_r;
      this->normal_p=1.0;
      this->normal_b=1.0;
      this->s_g=S_r_power;
      this->sg1=S_r1;
      this->sg2=S_r2;
    }
  else if(s_data->catalog=="random")
   {
      this->n_ran=n_selected;
      this->w_r=W_r;
      this->s_r=S_r_power;
      this->normal_p=normal_power;
      this->normal_b=normal_bispectrum;
      this->sr1=S_r1;
      this->sr2=S_r2;
    }

  std::cout<<YELLOW;
  time_t end; time (&end);
  real_prec diff = difftime(end,start);
  if (diff<60) std::cout <<"Lapse: "<<diff<<" secs"<<endl;
  else if (diff<3600) std::cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
  else std::cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
  std::cout<<RESET;
}
#endif



// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************

void FftwFunctions::raw_sampling(real_prec vol)
{
  //////////////////////////////////////////////////////////
  // Function used when no random catalogue is implemented.
  // A vector is filled with the mean number density of
  // the simulation. Such mean number density is computed
  // from the information if the size of the box
  // as given in the parameter file, together with the
  // number of objects.
  //////////////////////////////////////////////////////////

//  So.message_screen("Doing raw sampling of random grid");
//  omp_set_num_threads(omp_get_max_threads());


  /*
#pragma omp parallel for
  for(int i=0;i<this->NT;++i)
    data_r[i] = static_cast<real_prec>(this->n_gal)/static_cast<real_prec>(this->NT);
  So.DONE();
*/

  this->normal_p=static_cast<real_prec>(this->n_gal)*static_cast<real_prec>(this->n_gal)/vol;
  this->normal_b=static_cast<real_prec>(this->n_gal)*pow(this->n_gal/vol,2);
  this->w_r=static_cast<real_prec>(this->n_gal);
  this->s_r=vol/static_cast<real_prec>(this->n_gal);
  this->sr1=this->n_gal*(this->n_gal/vol);
  this->sr2=static_cast<real_prec>(this->n_gal);
  this->n_ran=0;
  this->alpha=1.0;
  this->normal_power=normal_p;
  this->shot_noise=this->s_r;
  this->shot_noise_window=0;
  this->normal_window=1.0 ;

}

// **************************************************************************************************
// **************************************************************************************************

void FftwFunctions::get_parameters_estimator(bool use_random, bool verbose)
{
  //////////////////////////////////////////////////////////  
  // This function computes the parameters associated to the
  // FKP estimator, e.g., normalization, shot noise, etc.
  //////////////////////////////////////////////////////////  

  if(use_random){
    this->alpha= (real_prec)w_g/(real_prec)w_r;
    this->normal_power= alpha* normal_p;
    //this->shot_noise= alpha*(alpha+1)*(s_r/ normal_power);    
    this->shot_noise= s_g/ normal_power + pow( alpha,2)*( s_r/ normal_power);
    this->normal_window= normal_power/pow(alpha,2);
    this->shot_noise_window=pow( alpha, 2)* s_r/ normal_power;
    this->normal_bispectrum= this->alpha* normal_b;     
    this->shot_noise_b1= sr1/( this->alpha* this->normal_bispectrum);
    this->shot_noise_b2=(1.- this->alpha* this->alpha)* sr2/( this->alpha* this->normal_bispectrum);

  }
  else{
    this->alpha=1.0; //This alfa should be zero when using a box. We set it to one correspondingly with that case, 
    // for all quantities are explicitely comnputed in that situation without assuming anything on this parameter
    this->normal= normal_power;
    this->shot_noise= s_r;
    this->normal_window=1.0;
    this->shot_noise_window=1.0; //care
    this->normal_bispectrum= normal_b;
    this->shot_noise_b1= this->sr1/(this->normal_b);
    this->shot_noise_b2= this->sr2/(normal_b);
  }



  
  if(statistics=="Pk_fkp" || statistics=="Pk_yb" || statistics=="Pk_ys"){
    s_parameters_estimator s_p_est={
      use_random,
      n_gal, 
      n_ran, 
      w_g,
      w_r,
      alpha,
      s_g,
      s_r,
      normal_power,
      normal_window,
      shot_noise,
      shot_noise_window
    };   
    if(true==verbose)
    So.write_parameters_estimator(&s_p_est);        /*Write parameters to the screen*/
    
  }
  else   if(statistics=="Bk_fkp" || statistics=="Bk_fkp_fast" || statistics=="Bk_ys_fast"){
    s_parameters_bis_estimator  s_p_est={
      use_random,
      n_gal, 
      n_ran, 
      w_g,
      w_r,
      alpha,
      sg1,
      sg2,
      sr1,
      sr2,
      normal_bispectrum,
      shot_noise_b1,
      shot_noise_b2,
      normal_power,
      shot_noise
    };   
    if(true==verbose)
    So.write_parameters_b_estimator(&s_p_est);        /*Write parameters to the screen*/
  }
}

// **************************************************************************************************
// **************************************************************************************************



void FftwFunctions::get_fluctuation(bool use_randoms)
{
  //////////////////////////////////////////////////////////  
  // Build the galaxy fluctuation by subtracting data 
  // and random catalogue with factor alpha
  //////////////////////////////////////////////////////////  


  // Define the coordinates of the cells within the grid
  // to be used in the Yamamoto estimator
  if(statistics=="Pk_ybc"  || statistics=="Pk_ysc"){
    for(int i=0;i<this->Nft;++i){
      for(int j=0;j<this->Nft;++j){
	for(int k=0;k<this->Nft;++k){
          int lp=index_3d(i ,j ,k, this->Nft,this->Nft);
	  cell_x[lp]=(i+0.5)*deltax + this->Xoffset - this->Lside*0.5;
	  cell_y[lp]=(j+0.5)*deltay + this->Yoffset - this->Lside*0.5;
	  cell_z[lp]=(k+0.5)*deltaz + this->Zoffset - this->Lside*0.5;
	}
      }
    }
    
  }
  

  
  if(statistics=="Pk_ds"){
    real_prec alpha = this->alpha;
    for(int i=0;i<this->NTT;++i){
      this->data_g_out_y0[i] -= alpha*this->data_r_out_y0[i];
      this->data_g_out_y2[i] -= alpha*this->data_r_out_y2[i];
      this->data_g_out_y4[i] -= alpha*this->data_r_out_y4[i];
    }
  }


  omp_set_num_threads(omp_get_max_threads());
#pragma omp parallel
  {
    
    real_prec alpha = this->alpha;

  if(true==use_randoms)
    {
#pragma omp for nowait
      for(int i=0;i<this->NT;++i)
        this->data_g[i]-=alpha*this->data_r[i];
    }
    else
    {

#ifdef _MASS_WEIGHT_POWER_
#pragma omp for nowait
     for(int i=0;i<this->NT;++i)
        this->data_g[i]=this->data_g_mw[i]-data_g[i];
#else
#pragma omp for nowait
     for(int i=0;i<this->NT;++i)
        this->data_g[i]-=static_cast<real_prec>(this->n_gal)/static_cast<real_prec>(this->NT);
#endif
    }



#pragma omp for nowait
    for(int i=0;i<this->NT;++i)
      {
	
	// Define the fluctuation F(r) interpolated 

        if(statistics=="Pk_ys")
	  {
	    // Seems this interpolation is correct
	    // Compute squared of the distance of the cell to the origin
            real_prec absr2=this->cell_x[i]*this->cell_x[i]+this->cell_y[i]*this->cell_y[i]+this->cell_z[i]*this->cell_z[i];
            //At this point this delta is already the fluctuation
            real_prec delta_cc=this->data_g[i];
	    
	    this->data_g_xx[i]=(cell_x[i]*cell_x[i]/absr2)*delta_cc;
	    this->data_g_yy[i]=(cell_y[i]*cell_y[i]/absr2)*delta_cc;
	    this->data_g_zz[i]=(cell_z[i]*cell_z[i]/absr2)*delta_cc;
	    this->data_g_xy[i]=(cell_x[i]*cell_y[i]/absr2)*delta_cc;
	    this->data_g_xz[i]=(cell_x[i]*cell_z[i]/absr2)*delta_cc;
	    this->data_g_yz[i]=(cell_y[i]*cell_z[i]/absr2)*delta_cc;
	  }
	
	
        else if(statistics=="Pk_yb")
	  {
	    // Seems this interpolation is correct
	    // Compute squared of the distance of the cell to the origin
            real_prec absr2=this->cell_x[i]*this->cell_x[i]+this->cell_y[i]*this->cell_y[i]+this->cell_z[i]*this->cell_z[i];
            real_prec iabsr4=1./(absr2*absr2);
	    //At this point this delta is already the fluctuation
            real_prec delta_cc=this->data_g[i];
	    
	    this->data_g_xx[i]=(cell_x[i]*cell_x[i]/absr2)*delta_cc;
	    this->data_g_yy[i]=(cell_y[i]*cell_y[i]/absr2)*delta_cc;
	    this->data_g_zz[i]=(cell_z[i]*cell_z[i]/absr2)*delta_cc;
	    this->data_g_xy[i]=(cell_x[i]*cell_y[i]/absr2)*delta_cc;
	    this->data_g_xz[i]=(cell_x[i]*cell_z[i]/absr2)*delta_cc;
	    this->data_g_yz[i]=(cell_y[i]*cell_z[i]/absr2)*delta_cc;
	    this->data_g_xxx[i]=iabsr4*(pow(cell_x[i],2)*cell_x[i]*cell_x[i])*delta_cc;
	    this->data_g_yyy[i]=iabsr4*(pow(cell_y[i],2)*cell_y[i]*cell_y[i])*delta_cc;
	    this->data_g_zzz[i]=iabsr4*(pow(cell_z[i],2)*cell_z[i]*cell_z[i])*delta_cc;
	    this->data_g_xxy[i]=iabsr4*(pow(cell_x[i],2)*cell_x[i]*cell_y[i])*delta_cc;
	    this->data_g_xxz[i]=iabsr4*(pow(cell_x[i],2)*cell_x[i]*cell_z[i])*delta_cc;
	    this->data_g_yyx[i]=iabsr4*(pow(cell_y[i],2)*cell_y[i]*cell_x[i])*delta_cc;
	    this->data_g_yyz[i]=iabsr4*(pow(cell_y[i],2)*cell_y[i]*cell_z[i])*delta_cc;
	    this->data_g_zzx[i]=iabsr4*(pow(cell_z[i],2)*cell_z[i]*cell_x[i])*delta_cc;
	    this->data_g_zzy[i]=iabsr4*(pow(cell_z[i],2)*cell_z[i]*cell_y[i])*delta_cc;
	    this->data_g_xyy[i]=iabsr4*(pow(cell_x[i],2)*cell_y[i]*cell_y[i])*delta_cc;
	    this->data_g_xzz[i]=iabsr4*(pow(cell_x[i],2)*cell_z[i]*cell_z[i])*delta_cc;
	    this->data_g_yzz[i]=iabsr4*(pow(cell_y[i],2)*cell_z[i]*cell_z[i])*delta_cc;
	    this->data_g_xyz[i]=iabsr4*(pow(cell_x[i],2)*cell_y[i]*cell_z[i])*delta_cc;
	    this->data_g_yxz[i]=iabsr4*(pow(cell_y[i],2)*cell_x[i]*cell_z[i])*delta_cc;
	    this->data_g_zxy[i]=iabsr4*(pow(cell_z[i],2)*cell_x[i]*cell_y[i])*delta_cc;
	  }
      }
  }

}
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************

void FftwFunctions::get_power_spectrum_fkp(s_parameters_box *s_box, vector<real_prec> &power_g0, vector<real_prec> &power_g2,vector<real_prec> &power_g4,vector<real_prec> &power_r,  vector< vector<real_prec> >&power2d_cart,vector< vector<real_prec> >&power2d_spher,vector<int> &mod_g)
{
  //////////////////////////////////////////////////////////  
  // Computes the outputs for the P(k),P(kp, kp), P(k,mu)
  // and W(k)(
  // Also generates the quadrupole and hexadecapole.
  // With this, the 2d() functions are no longer needed
  // and we only do onle loop over Fourier space to retrieve
  // all possible forms of the 3d power spectrum.
  //////////////////////////////////////////////////////////  
  
  time_t start;
  time (&start);
  
  do_fftw_r2c(this->Nft,this->data_g, this->data_out_g);
  
  if(true==this->measure_cross)
    do_fftw_r2c(this->Nft,this->data_gp, this->data_out_gp);
  
  if(true==s_box->use_random_catalog)
    do_fftw_r2c(this->Nft,this->data_r, this->data_out_r);
  
  time (&start); 
  this->power_spectrum_fkp(s_box,power_g0,power_g2,power_g4,power_r,power2d_cart,power2d_spher,mod_g); 
}


//##################################################################################
//##################################################################################
void FftwFunctions::get_power_spectrum_fkp(s_parameters_box *s_box, vector<real_prec> &power_g0,vector<int> &mod_g)
{
  do_fftw_r2c(this->Nft,this->data_g, this->data_out_g);
  power_spectrum_fkp(s_box,power_g0,mod_g); 
}


// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************

// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************

void FftwFunctions::get_power_spectrum_yamamoto(s_parameters_box *s_box, vector<real_prec> &power_g0, vector<real_prec> &power_g2,vector<real_prec> &power_g4,vector<int> &mod_g)
{
  //////////////////////////////////////////////////////////  
  // Generates the quadrupole and hexadecapole.
  // For Yamamotos
  //////////////////////////////////////////////////////////  
  
  
  if(this->statistics=="Pk_y_ds"){
    
    time_t start;
    time (&start);
    
    real_prec sn_aux =  s_box->use_SN_correction ? 1.0: 0.0;
    //    if(s_box->use_SN_correction)sn_aux=1.0; else sn_aux=0.;  
    
    real_prec r_normal_power = 1.0/this->normal_power;
    real_prec alpha2 = this->alpha * this->alpha;
    
    
#ifdef _USE_OMP_
    omp_set_num_threads(omp_get_max_threads());
#pragma omp parallel for
#endif
    for(int i=0;i<this->NTT;++i){
      
      // *******************************************************
      // Moment l=0: should be equal to the monopole from FKP                                  
      // Given that this->shot_noise ya viene normalizado, multiplico por normal_power
      // para cancelar y dividir todo el parentesis de nuevo por la normalizacion
      this->data_g_y0[i] = (norm(this->data_g_out_y0[i]) - sn_aux*this->shot_noise*this->normal_power) * r_normal_power;   
      
      // Moment l=2                                                                           
      this->data_g_y2[i]=(2.*2.0+1.)*(this->data_g_out_y2[i].real()*this->data_g_out_y0[i].real()+this->data_g_out_y2[i].imag()*this->data_g_out_y0[i].imag()
				      - sn_aux*(this->SN_g_out_y2[i]+alpha2*SN_r_out_y2[i])) * r_normal_power;
      
      // ******************************************************
      // Moment l=4
      this->data_g_y4[i]=(2.*4.0+1.)*(this->data_g_out_y4[i].real()*this->data_g_out_y0[i].real()+this->data_g_out_y4[i].imag()*this->data_g_out_y0[i].imag()
				      - sn_aux*(this->SN_g_out_y4[i]+ alpha2*this->SN_r_out_y4[i])) * r_normal_power;
      
    }
    
    // Shell averages:
    power_yam_1d_ds(s_box, power_g0, power_g2, power_g4,mod_g);
    
    
    
    std::cout<<RED;
    time_t end;
    time(&end);
    real_prec diff=difftime(end,start);
    if (diff<60) std::cout <<"Lapse:  "<<diff<<" secs"<<endl;
    else if (diff<3600) std::cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
    else std::cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
    std::cout<<RESET;
    
  }
  
  else{
    
    // Do FFTS
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads()); 

    So.message_screen("Evaluating DFT using ", omp_get_max_threads(), " threads");
    
    time_t start;
    time (&start);
      

    do_fftw_r2c(this->Nft,this->data_g, this->data_out_g);
    
    if(this->statistics=="Pk_ys"  || this->statistics=="Pk_ysc"){
      do_fftw_r2c(this->Nft,this->data_g_xx, this->data_out_g_xx);
      do_fftw_r2c(this->Nft,this->data_g_zz, this->data_out_g_zz);
      do_fftw_r2c(this->Nft,this->data_g_xy, this->data_out_g_xy);
      do_fftw_r2c(this->Nft,this->data_g_xz, this->data_out_g_xz);
      do_fftw_r2c(this->Nft,this->data_g_yz, this->data_out_g_yz);
    }   
    
    
    //     // If we want to implement the approach by Bianchi et al, 
    //     // we need to compute the following Fourier transforms
    
    else if(this->statistics=="Pk_yb" || this->statistics=="Pk_ybc")
      {
	do_fftw_r2c(this->Nft,this->data_g_xx, this->data_out_g_xx);
	do_fftw_r2c(this->Nft,this->data_g_zz, this->data_out_g_zz);
	do_fftw_r2c(this->Nft,this->data_g_xy, this->data_out_g_xy);
	do_fftw_r2c(this->Nft,this->data_g_xz, this->data_out_g_xz);
	do_fftw_r2c(this->Nft,this->data_g_yz, this->data_out_g_yz);
	do_fftw_r2c(this->Nft,this->data_g_xxx, this->data_out_g_xxx);
	do_fftw_r2c(this->Nft,this->data_g_yyy, this->data_out_g_yyy);
	do_fftw_r2c(this->Nft,this->data_g_zzz, this->data_out_g_zzz);
	do_fftw_r2c(this->Nft,this->data_g_xxy, this->data_out_g_xxy);
	do_fftw_r2c(this->Nft,this->data_g_xxz, this->data_out_g_xxz);
	do_fftw_r2c(this->Nft,this->data_g_yyx, this->data_out_g_yyx);
	do_fftw_r2c(this->Nft,this->data_g_yyz, this->data_out_g_yyz);
	do_fftw_r2c(this->Nft,this->data_g_zzx, this->data_out_g_zzx);
	do_fftw_r2c(this->Nft,this->data_g_zzy, this->data_out_g_zzy);
	do_fftw_r2c(this->Nft,this->data_g_xyy, this->data_out_g_xyy);
	do_fftw_r2c(this->Nft,this->data_g_xzz, this->data_out_g_xzz);
	do_fftw_r2c(this->Nft,this->data_g_yzz, this->data_out_g_yzz);
	do_fftw_r2c(this->Nft,this->data_g_xyz, this->data_out_g_xyz);
	do_fftw_r2c(this->Nft,this->data_g_yxz, this->data_out_g_yxz);
	do_fftw_r2c(this->Nft,this->data_g_zxy, this->data_out_g_zxy);
    }
  
    if(s_box->use_random_catalog)
      do_fftw_r2c(this->Nft,this->data_r, this->data_out_r);
    
    std::cout<<RED;
    time_t end;
    time(&end);
    real_prec diff=difftime(end,start);
    if (diff<60) std::cout <<"Lapse:  "<<diff<<" secs"<<endl;
    else if (diff<3600) std::cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
    else std::cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
    std::cout<<RESET;
    
    // ***************************************************************
    // Compute shell-average //
    time (&start);
    So.message_screen("Shell average");
    power_spectrum_yamamoto(s_box,power_g0,power_g2,power_g4,mod_g); 
    So.DONE();
    // ***************************************************************


    std::cout<<RED;
    time(&end);
    diff=difftime(end,start);
    if (diff<60) std::cout <<"Lapse:  "<<diff<<" secs"<<endl;
    else if (diff<3600) std::cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
    else std::cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
    std::cout<<RESET;
    
  }
  
  
  
}

// **************************************************************************************************
// **************************************************************************************************


void FftwFunctions::power_spectrum_fkp(s_parameters_box *s_box,
                                       vector<real_prec> & powerk_g0,
                                       vector<real_prec> & powerk_g2,
                                       vector<real_prec> & powerk_g4,
                                       vector<real_prec> & powerk_r,
                                       vector< vector<real_prec> >&power2d_cart,
                                       vector< vector<real_prec> >&power2d_spher,
				       vector<int> & mod_g)
{
  //////////////////////////////////////////////////////////  
  // Shell average in Fourier space. 
  // This function returns the spherical average estimate 
  // of the monopole, quadrupole, hexadecapole, and the 
  // window function of based on the FKP estiamtor.
  // Only the vector with the modes used for the actual P(K)
  // (instead of those used for W(k)- though the same
  // eventually-, since these will be usweful when 
  // computing the variance using the Veff approach.
  // Note that we do loops only through an octant of the volume sampled
  // (i.e, kz>0 modes) in Fourier space. The other
  // three octants are obtained by retreiving the indices
  // The account for the kz<0 region, we multuply by two.
  // Again, this is valid since we want coordinates squared.
  // For other uses (e.g., the bispectrum,) we must be explicit.
  //////////////////////////////////////////////////////////  
  
#ifdef _WRITE_2DPOWER_
  vector < vector<int> > mod_cart(power2d_cart.size(), vector<int>(power2d_cart[0].size(),0));
  vector < vector<int> > mod_spher(power2d_spher.size(), vector<int>(power2d_spher[0].size(),0));
#endif
  vector<int> mod_r(powerk_r.size(),0);
  
#pragma omp parallel for
  for(int i=0;i<mod_g.size();++i)
    {
      mod_g[i]=0;
      powerk_g0[i]=0;
#ifdef _WRITE_MULTIPOLES_
      powerk_g2[i]=0;
      powerk_g4[i]=0;
#endif
  }

  real_prec rDeltaK_data;
  real_prec rDeltaK_window;
  real_prec rkmin;
  
  int my_s_box_ave;
  real_prec k_bin_index=(s_box->k_bin_step >= 1.? 0.5:0.0);
  
  
#define s_box_linear 1
#define s_box_log 2
#define s_box_other 0
  
  my_s_box_ave = s_box_other;

  if(s_box->ave=="linear") {
     my_s_box_ave = s_box_linear;
    rDeltaK_data = 1.0 / this->DeltaK_data;
    rDeltaK_window = 1.0 / this->DeltaK_window;
  }
  else if(s_box->ave=="log")
   {
    my_s_box_ave = s_box_log;
    rkmin = 1.0/this->kmin;
   }

  vector<real_prec>yy_MAS_array(this->Nft/2,0);
  vector<real_prec>zz_MAS_array(this->Nft/2+1,0);
  vector<real_prec>sn_yy_MAS_array(this->Nft/2,0);
  vector<real_prec>sn_zz_MAS_array(this->Nft/2+1,0);


  int NTHREADS = omp_get_max_threads();
  
  real_prec M_PI_rNft = M_PI * this->rNft;
  
  time_t start;
  time (&start);


  real_prec SN_aux=0;
  if(true==s_box->use_SN_correction)
    SN_aux=this->shot_noise;

  real_prec SN_aux2=0;
  if(true==s_box->use_SN_correction)
    SN_aux2=this->shot_noise2;

#pragma omp parallel num_threads(NTHREADS)
  {
    real_prec i_per_fact, j_per_fact;
    real_prec i_deltak_x_2, j_deltak_y_2, k_deltak_z_2;
    real_prec xx_MAS, yy_MAS, zz_MAS;
    real_prec sn_xx_MAS, sn_yy_MAS, sn_zz_MAS;

    
#pragma omp sections
    {
#pragma omp section
      {
        yy_MAS_array[0] = 1.0;
        sn_yy_MAS_array[0] = 0.0;
        for(int j=1;j<this->Nft/2;++j) {
          yy_MAS = j * M_PI_rNft;
          yy_MAS_array[j] = sin(yy_MAS)/yy_MAS;
          sn_yy_MAS_array[j] = sin(j * M_PI_rNft);
	}
      }
#pragma omp section
      {
        zz_MAS_array[0] = 1.0;
        sn_zz_MAS_array[0] = 0.0;
        for(int k=1;k<=this->Nft/2;++k){  // loop over octant kz>0: this is all FFTW gives.
          zz_MAS = k * M_PI_rNft;
          zz_MAS_array[k] = (zz_MAS = sin(zz_MAS)/zz_MAS);
          sn_zz_MAS_array[k] = sin(k * M_PI_rNft);
         }
      }
    }
    
#pragma omp barrier
    
#pragma omp for nowait  
    for(int i=0;i<this->Nft/2;++i)
      {  // loop over octant kx>0. Half of what FFT gives
	if(i > 0)
	  {
	    xx_MAS = i * M_PI_rNft;
	    xx_MAS = sin(xx_MAS)/xx_MAS;
	  }
	else xx_MAS = 1.0;
	
        sn_xx_MAS = sin(i * M_PI_rNft);
	
	
	i_deltak_x_2 = i * i * this->deltak_x * this->deltak_x ;
	i_per_fact = i_deltak_x_2;
	if(my_s_box_ave == s_box_log)
	  i_deltak_x_2 *= rkmin * rkmin;
	
	for(int j=0;j<this->Nft/2;++j)
	  {  // loop over octant ky>0. Half of what FFT gives
	    yy_MAS = yy_MAS_array[j];
            sn_yy_MAS=sn_yy_MAS_array[j];
	    
	    j_deltak_y_2 = j * j * this->deltak_y * this->deltak_y;
	    j_per_fact = i_per_fact + j_deltak_y_2;
	    j_per_fact = sqrt(j_per_fact);
	    if(my_s_box_ave == s_box_log)
	      j_deltak_y_2 *= rkmin * rkmin;
	    
	    for(int k=0;k<=this->Nft/2;++k)
	      {  // loop over octant kz>0: this is all FFTW gives.
		zz_MAS = zz_MAS_array[k];
                sn_zz_MAS=sn_zz_MAS_array[k];

		k_deltak_z_2 = k * k * this->deltak_z * this->deltak_z;
		if(my_s_box_ave == s_box_log)
		  k_deltak_z_2 *= rkmin * rkmin;
		
                real_prec kv;
		int kmod_g, kmod_r;
		
		// Compute k-shell index 
		if(my_s_box_ave == s_box_linear)
		  {
		    kv=sqrt(i_deltak_x_2 + j_deltak_y_2 + k_deltak_z_2);
                    kmod_g=(int)floor(static_cast<float>(kv * rDeltaK_data)+k_bin_index);
                    kmod_r=(int)floor(static_cast<float>(kv * rDeltaK_window));
		  }
		else
		  {
		    if(my_s_box_ave == s_box_log)
		      {
			kv= i_deltak_x_2 + j_deltak_y_2 + k_deltak_z_2;
			if(kv!=0){
                          kmod_g=(int)floor(static_cast<float>(log10(sqrt(kv))/this->Deltal))+1;
			  kmod_r=kmod_g;
			}
			else{kmod_g=0;kmod_r=0;}
		      }
		  }
		
		
		
		// Compute correction for MAS:
                real_prec corr=(s_box->use_MAS_correction ? my_correction_MAS(correction_MAS_exp, xx_MAS, yy_MAS, zz_MAS): 1.0);
                real_prec icorr2 = 1.0/(corr*corr);
		
                real_prec sn_corr=SN_aux*SN_correction_MAS(correction_MAS_exp, sn_xx_MAS, sn_yy_MAS, sn_zz_MAS);



		// Compute index in c-order
                int lp=index_3d(i,j,k,this->Nft,this->Nft/2+1);
		
		
		// Compute thisgs related to 2d power spectrum and multipole 
		// decomposition for the first octant
		// Define the angle bewtween k and los:       
		// In the FKP we need to specify the LOS direction. Set to z by default.
		
#ifdef _WRITE_2DPOWER_
		// Compute mu = kz/k. 
                real_prec mu = k*this->deltak_z/kv;
                int i_par = static_cast<int>(floor((float)(this->deltak_z*k * rDeltaK_data)));
                int i_per = static_cast<int>(floor((float)(j_per_fact * rDeltaK_data)));
#endif

#ifdef _WRITE_MULTIPOLES_
                // Compute index in mu-direction
                real_prec i_mu  = static_cast<int>(floor((float)((mu-(-1.0))/Deltamu)));
		i_mu  = (i_mu==this->N_mu_bins ? 0: i_mu); 

		// Compute Legendre functions neeeded for the Multipole decomposition:
                real_prec leg2=0.5*(3.*mu*mu-1);
                real_prec leg4=(35.*pow(mu,4)-30.*mu*mu+3.)/8.;
#endif

		//  *****************************************
		// Now, add modes in the first octant
		// Note that k>0 from here
		// Note also that Leg are invariant under the
		// shift os sigh kz-> - kz, and therefore we do not 
		// need explicitely to account for the kz<0 region; 
		// We can simply multiply by two both in the power as
		// in the number of counted modes. This factor two cancells
		// in the end, and therefore we don't employ it.
		//  *****************************************
		//Exclude the zero-frequency onLY when computing the P: for the W keep it:
		if(lp!=0)
		  {
                    real_prec Pk, Pk1, Pk2;
		    if(true==measure_cross)
		      {
			Pk= (this->data_out_g[lp][0] * this->data_out_gp[lp][0] + this->data_out_g[lp][1] * this->data_out_gp[lp][1]) * icorr2;
			Pk1= (this->data_out_g[lp][0] * this->data_out_g[lp][0] + this->data_out_g[lp][1] * this->data_out_g[lp][1]) * icorr2;
			Pk2= (this->data_out_gp[lp][0] * this->data_out_gp[lp][0] + this->data_out_gp[lp][1] * this->data_out_gp[lp][1]) * icorr2;
			//Pk=Pk/sqrt((Pk1-SN_aux)*(Pk2-SN_aux2));
		      }
		    else
                      Pk= (this->data_out_g[lp][0] * this->data_out_g[lp][0]  + this->data_out_g[lp][1] * this->data_out_g[lp][1]  - sn_corr*normal_power ) * icorr2 ;
		    
		    
		    // 1d Spherical average
		    if(kmod_g<powerk_g0.size())
		      {
			if(true==measure_cross)
			  {
#pragma omp atomic update             
			    powerk_g0[kmod_g]  += Pk;
#pragma omp atomic update             
			    powerk_g2[kmod_g]  += Pk1;
#pragma omp atomic update       
			    powerk_g4[kmod_g]  += Pk2;
			  }
			else
			  {
#pragma omp atomic update             
			    powerk_g0[kmod_g]  += Pk;
#ifdef _WRITE_MULTIPOLES_
#pragma omp atomic update             
			    powerk_g2[kmod_g]  += leg2*Pk;
#pragma omp atomic update       
			    powerk_g4[kmod_g]  += leg4*Pk;
#endif
			  }
#pragma omp atomic update       
			mod_g[kmod_g]++ ;    
		      }

#ifdef _WRITE_2DPOWER_
		    // 2d in Spherical coordinates:
		    if(i_mu<power2d_spher.size() && kmod_g<power2d_spher[0].size())
		      {
#pragma omp atomic update              
			power2d_spher[i_mu][kmod_g]+= Pk;
#pragma omp atomic  update             
			mod_spher[i_mu][kmod_g]++;    
		      }
		    
		    // 2d in Cartesian Coordinates:
                    if(i_par<power2d_cart.size() && i_per<power2d_cart[0].size())
                      {
#pragma omp atomic  update                           
		      power2d_cart[i_par][i_per] += Pk;
#pragma omp atomic  update                           
		      mod_cart[i_par][i_per] ++ ;    
		    }
#endif
                }

                else
                   {  //Compute 1d spherical average for the window function only
                   if(s_box->use_random_catalog)
		    {
		      if(kmod_r<powerk_r.size()){
#pragma omp atomic                              
			powerk_r[kmod_r]  += (data_out_r[lp][0] * data_out_r[lp][0] + data_out_r[lp][1]*data_out_r[lp][1]) * icorr2;
#pragma omp atomic                            
			mod_r[kmod_r]++ ;    
		      }
		    }
		}
		
		
		
		//  *****************************************
		//  Now, add modes in other octants:
		
		//  *****************************************
		//  Add negative frequencies in y:
		//  kv is the same, for ky->-ky.
		//  The same applies below for kx.
		//  Therefore, kmod_g is also the same.
		//  We only compute lp again to go to the right octant
		//  *****************************************
		if(j>0  && k>0)
		  {    
		    lp=ijk(i,this->Nft-j,k,this->Nft,this->Nft,this->Nft/2+1);   // Go to the ky<0 octant
		    
                    real_prec Pk, Pk1,Pk2;
		    if(true==measure_cross)
		      {
			Pk= (this->data_out_g[lp][0] * this->data_out_gp[lp][0] + this->data_out_g[lp][1] * this->data_out_gp[lp][1]) * icorr2;
			Pk1= (this->data_out_g[lp][0] * this->data_out_g[lp][0] + this->data_out_g[lp][1] * this->data_out_g[lp][1]) * icorr2;
			Pk2= (this->data_out_gp[lp][0] * this->data_out_gp[lp][0] + this->data_out_gp[lp][1] * this->data_out_gp[lp][1]) * icorr2;
			//Pk=Pk/sqrt((Pk1-SN_aux)*(Pk2-SN_aux2));
		      }
		    else 
                      Pk= (this->data_out_g[lp][0] * this->data_out_g[lp][0] + this->data_out_g[lp][1] * this->data_out_g[lp][1]  - sn_corr*normal_power) * icorr2 ;
		    
		    
		    // 1d Spherical average
		    if(kmod_g<powerk_g0.size())
		      {
			if(true==measure_cross)
			  {
#pragma omp atomic  update                             
			    powerk_g0[kmod_g]  += Pk;
#pragma omp atomic  update                             
			    powerk_g2[kmod_g]  += Pk1;
#pragma omp atomic  update             
			    powerk_g4[kmod_g]  += Pk2;
			  }

			else
			  {
#pragma omp atomic  update                             
			    powerk_g0[kmod_g]  += Pk;
#ifdef _WRITE_MULTIPOLES_
#pragma omp atomic  update                             
			    powerk_g2[kmod_g]  += leg2*Pk;
#pragma omp atomic  update             
			    powerk_g4[kmod_g]  += leg4*Pk;
#endif
			  }
			
#pragma omp atomic update                              
			mod_g[kmod_g]++ ;    
		      }

#ifdef _WRITE_2DPOWER_
		    // 2d in Spherical coordinates:
		    if(i_mu<power2d_spher.size() && kmod_g<power2d_spher[0].size())
		      {
#pragma omp atomic update                           
			power2d_spher[i_mu][kmod_g]+= Pk;
#pragma omp atomic update                         
			mod_spher[i_mu][kmod_g]++;    
		      }
		    // 2d in Cartesian Coordinates:
                    if(i_par<power2d_cart.size() && i_per<power2d_cart[0].size())
                      {
#pragma omp atomic update                              
		      power2d_cart[i_par][i_per] += Pk;
#pragma omp atomic update                              
		      mod_cart[i_par][i_per] ++ ;    
		    }
#endif
		    //Compute 1d spherical average for the window function only
		    if(s_box->use_random_catalog){
		      if(kmod_r<powerk_r.size()){
#pragma omp atomic                                
			powerk_r[kmod_r]  += (data_out_r[lp][0]*data_out_r[lp][0] + data_out_r[lp][1]*data_out_r[lp][1]) * icorr2;
#pragma omp atomic                              
			mod_r[kmod_r]++ ;    
		      }
		    }
		  }
		
		
		
		// ******************************************
		// Add negative frequencies in x:
		// ******************************************
		if(i>0  && (j>0 || k>0)){
		  lp=ijk(this->Nft-i,j,k,this->Nft,this->Nft,this->Nft/2+1);
		  
                  real_prec Pk,Pk1,Pk2;
		  if(true==measure_cross)
		    {
		      Pk= (this->data_out_g[lp][0] * this->data_out_gp[lp][0] + this->data_out_g[lp][1] * this->data_out_gp[lp][1]) * icorr2;
		      Pk1= (this->data_out_g[lp][0] * this->data_out_g[lp][0] + this->data_out_g[lp][1] * this->data_out_g[lp][1]) * icorr2;
		      Pk2= (this->data_out_gp[lp][0] * this->data_out_gp[lp][0] + this->data_out_gp[lp][1] * this->data_out_gp[lp][1]) * icorr2;
		      //Pk=Pk/sqrt((Pk1-SN_aux)*(Pk2-SN_aux2));
		    }
		  else
		    if(false==measure_cross)
                      Pk= (this->data_out_g[lp][0] * this->data_out_g[lp][0] + this->data_out_g[lp][1] * this->data_out_g[lp][1]  - sn_corr*normal_power) * icorr2 ;
		  
		  
		  if(kmod_g<powerk_g0.size())
		    {
		      if(true==measure_cross)
			{
#pragma omp atomic update                                 
			  powerk_g0[kmod_g]  += Pk;
#pragma omp atomic update                              
			  powerk_g2[kmod_g]  += Pk1;
#pragma omp atomic update                              
			  powerk_g4[kmod_g]  += Pk2;
			}
		      
		      else{
#pragma omp atomic update                                 
			powerk_g0[kmod_g]  += Pk;
#ifdef _WRITE_MULTIPOLES_

#pragma omp atomic update                              
			powerk_g2[kmod_g]  += leg2*Pk;
#pragma omp atomic update                              
			powerk_g4[kmod_g]  += leg4*Pk;
#endif
		      }
		      
#pragma omp atomic update                              
		      mod_g[kmod_g]++;
		    }    

#ifdef _WRITE_2DPOWER_
                  // 2d in Spherical coordinates:
		  if(i_mu<power2d_spher.size() && kmod_g<power2d_spher[0].size())
		    {
#pragma omp atomic update                                
		      power2d_spher[i_mu][kmod_g]+= Pk;
#pragma omp atomic update                              
		      mod_spher[i_mu][kmod_g]++;    
		    }
		  // 2d in Cartesian Coordinates:
		  if(i_par<power2d_cart.size() && i_per<power2d_cart[0].size())
		    {
#pragma omp atomic update                                
		      power2d_cart[i_par][i_per] += Pk;
#pragma omp atomic update                              
		      mod_cart[i_par][i_per] ++ ;    
		    }
#endif
		  //Compute 1d spherical average for the window function only
		  if(s_box->use_random_catalog)
		    {
		      if(kmod_r<powerk_r.size())
			{
#pragma omp atomic                                
			  powerk_r[kmod_r]  += (data_out_r[lp][0]*data_out_r[lp][0]+ data_out_r[lp][1]*data_out_r[lp][1]) * icorr2;
#pragma omp atomic                              
			  mod_r[kmod_r]++ ;    
			}
		    }
		}
		
		// ******************************************
		//  Add negative frequencies in x and y:
		// ******************************************
		if(i>0  && j>0  && k>0){
		  lp=ijk(this->Nft-i,this->Nft-j,k,this->Nft,this->Nft,this->Nft/2+1);
		  
                  real_prec Pk,Pk1,Pk2;
		  if(true==measure_cross)
		    {
		      Pk= (this->data_out_g[lp][0] * this->data_out_gp[lp][0] + this->data_out_g[lp][1] * this->data_out_gp[lp][1]) * icorr2;
		      Pk1= (this->data_out_g[lp][0] * this->data_out_g[lp][0] + this->data_out_g[lp][1] * this->data_out_g[lp][1]) * icorr2;
		      Pk2= (this->data_out_gp[lp][0] * this->data_out_gp[lp][0] + this->data_out_gp[lp][1] * this->data_out_gp[lp][1]) * icorr2;
		      //Pk=Pk/sqrt((Pk1-SN_aux)*(Pk2-SN_aux2));
		    }
		  else if(false==measure_cross)
                    Pk= (this->data_out_g[lp][0] * this->data_out_g[lp][0] + this->data_out_g[lp][1] * this->data_out_g[lp][1] - sn_corr*normal_power) * icorr2  ;
		  
		  
		  if(kmod_g<powerk_g0.size())
		    {
		      if(true==measure_cross)
			{
#pragma omp atomic update                                  
			  powerk_g0[kmod_g] += Pk;
#pragma omp atomic update                                  
			  powerk_g2[kmod_g] += Pk1;
#pragma omp atomic update                                   
			  powerk_g4[kmod_g] += Pk2;
			  
			}
		      else{
#pragma omp atomic update                                  
		      powerk_g0[kmod_g] += Pk;
#ifdef _WRITE_MULTIPOLES_
#pragma omp atomic update                                  
			powerk_g2[kmod_g] += leg2*Pk;
#pragma omp atomic update                                   
			powerk_g4[kmod_g] += leg4*Pk;
#endif

                      }
#pragma omp atomic update
		      mod_g[kmod_g] ++ ;    
		    }
	    
#ifdef _WRITE_2DPOWER_
                  // 2d in Spherical coordinates:
	    if(i_mu<power2d_spher.size() && kmod_g<power2d_spher[0].size())
	      {
#pragma omp atomic                                
		power2d_spher[i_mu][kmod_g]+= Pk;
#pragma omp atomic                                
		mod_spher[i_mu][kmod_g]++;    
	      }
	    // 2d in Cartesian Coordinates:
	    if(i_par<power2d_cart.size() && i_per<power2d_cart[0].size())
	      {
#pragma omp atomic update                                
		power2d_cart[i_par][i_per] += Pk;
#pragma omp atomic update                                
		mod_cart[i_par][i_per] ++ ;    
	      }
#endif
	    //Compute 1d spherical average for the window function only
	    
	    if(s_box->use_random_catalog)
	      {
		if(kmod_r<powerk_r.size())
		  {
#pragma omp atomic update                                  
		    powerk_r[kmod_r] += (data_out_r[lp][0]*data_out_r[lp][0]+data_out_r[lp][1]*data_out_r[lp][1]) * icorr2;
#pragma omp atomic update                                  
		    mod_r[kmod_r] ++ ;    
		  }
	      }
		}
	      }
	  }
      } // END of parallel loop
    
  } // END of parallel regione  
  
  
  // Subtract shot noise and normalize only he 1d power spectrum.
 
  // For the monopole:
  if(true==measure_cross) // use the arrays g2 and g4 to get auto power needed for the correlation
    {
#pragma omp parallel for
      for(int i=0;i<powerk_g0.size();++i)
	powerk_g0[i]=(mod_g[i]== 0 ? 0 : powerk_g0[i]/(static_cast<real_prec>(mod_g[i])*normal_power)); // use this if corr is defined using shell averages of spectra
      
#pragma omp parallel for
      for(int i=0;i<powerk_g2.size();++i)
        powerk_g2[i]=(mod_g[i]== 0 ? 0 : powerk_g2[i]/(static_cast<real_prec>(mod_g[i])*normal_power));
      
      // For the hexaecapole:
#pragma omp parallel for
      for(int i=0;i<powerk_g4.size();++i)
        powerk_g4[i]=(mod_g[i]== 0 ? 0 : powerk_g4[i]/(static_cast<real_prec>(mod_g[i])*normal_power));
      
      //      Get the cross power normalized:// use this if corr is defined using shell averages of spectra
#pragma omp parallel for
      for(int i=0;i<powerk_g0.size();++i)
	powerk_g0[i]/=static_cast<real_prec>(sqrt(powerk_g2[i]*powerk_g4[i]));

      
    }

  else{
#pragma omp parallel for
    for(int i=0;i<powerk_g0.size();++i)
      powerk_g0[i]=(mod_g[i]== 0 ? 0 : powerk_g0[i]/(static_cast<real_prec>(mod_g[i])*normal_power));

#ifdef _WRITE_MULTIPOLES_
#pragma omp parallel for
    for(int i=0;i<powerk_g2.size();++i)
      powerk_g2[i]=(mod_g[i]== 0 ? 0 : 5.*powerk_g2[i]/((real_prec)mod_g[i]*normal_power));
    
    // For the hexaecapole:
#pragma omp parallel for
    for(int i=0;i<powerk_g4.size();++i)
      powerk_g4[i]=(mod_g[i]== 0 ? 0 : 9.*powerk_g4[i]/((real_prec)mod_g[i]*normal_power));
#endif

  }
  
  
  // For the quadrupole:
  
  // Also for the window function:
#pragma omp parallel for
  for(int i=0;i<powerk_r.size();++i)
    powerk_r[i]=(mod_r[i]== 0 ? 0 : powerk_r[i]/(static_cast<real_prec>(mod_r[i])*normal_window));
  
#ifdef _WRITE_2DPOWER_
  // And for the 2d estimates in cartesian coordinates
#pragma omp parallel for
  for(int i=0;i<power2d_cart.size();++i)
    for(int j=0;j<power2d_cart[0].size();++j)
      power2d_cart[i][j]=(mod_cart[i][j]==0? 0. : power2d_cart[i][j]/(( static_cast<real_prec>((mod_cart[i][j]))*normal_power)-SN_aux);
  
  // And for the 2d estimates in polar coordinates
#pragma omp parallel for
  for(int i=0;i<power2d_spher.size();++i)
    for(int j=0;j<power2d_spher[0].size();++j)
      power2d_spher[i][j]=(mod_spher[i][j]==0? 0. : power2d_spher[i][j]/((static_cast<real_prec>((mod_spher[i][j]))*normal_power)-SN_aux);
#endif
    
}



// ====================================================
// This function is an override of the previous function, just to compute the monopole
void FftwFunctions::power_spectrum_fkp(s_parameters_box *s_box,vector<real_prec> & powerk_g0,vector<int> & mod_g)
{
  
#pragma omp parallel for
  for(int i=0;i<mod_g.size();++i)
    {
      mod_g[i]=0;
      powerk_g0[i]=0;
    }
  
  real_prec rDeltaK_data;
  real_prec rDeltaK_window;
  real_prec rkmin;
  
  int my_s_box_ave;
  real_prec k_bin_index=(s_box->k_bin_step >=1.? 0.5:0.0);
  
#define s_box_linear 1
#define s_box_log 2
#define s_box_other 0
  
  my_s_box_ave = s_box_other;
  
  if(s_box->ave=="linear") {
     my_s_box_ave = s_box_linear;
    rDeltaK_data = 1.0 / this->DeltaK_data;
    rDeltaK_window = 1.0 / this->DeltaK_window;
  }
  else if(s_box->ave=="log")
   {
    my_s_box_ave = s_box_log;
    rkmin = 1.0/this->kmin;
   }


  vector<real_prec>yy_MAS_array(this->Nft/2,0);
  vector<real_prec>zz_MAS_array(this->Nft/2+1,0);
  vector<real_prec>sn_yy_MAS_array(this->Nft/2,0);
  vector<real_prec>sn_zz_MAS_array(this->Nft/2+1,0);

  
  int NTHREADS = omp_get_max_threads();
  
  real_prec M_PI_rNft = M_PI * this->rNft;
  
  time_t start;
  time (&start);
  


  real_prec SN_aux=0;
  if(true==s_box->use_SN_correction)
    SN_aux=this->shot_noise;

  

#pragma omp parallel num_threads(NTHREADS)
  {
    real_prec i_per_fact, j_per_fact;
    real_prec i_deltak_x_2, j_deltak_y_2, k_deltak_z_2;
    real_prec xx_MAS, yy_MAS, zz_MAS;
    real_prec sn_xx_MAS, sn_yy_MAS, sn_zz_MAS;
    
#pragma omp sections
    {
#pragma omp section
      {
        yy_MAS_array[0] = 1.0;
        sn_yy_MAS_array[0] = 0.0;
        for(int j=1;j<this->Nft/2;++j) {
          yy_MAS = j * M_PI_rNft;
          yy_MAS_array[j] = sin(yy_MAS)/yy_MAS;
          sn_yy_MAS_array[j] = sin(j * M_PI_rNft);
        }
      }
#pragma omp section
      {
        zz_MAS_array[0] = 1.0;
        sn_zz_MAS_array[0] = 0.0;
        for(int k=1;k<=this->Nft/2;++k){  // loop over octant kz>0: this is all FFTW gives.
          zz_MAS = k * M_PI_rNft;
          zz_MAS_array[k] = (zz_MAS = sin(zz_MAS)/zz_MAS);
          sn_zz_MAS_array[k] = sin(k * M_PI_rNft);
        }
      }
    }
    
    
    
#pragma omp barrier
    
#pragma omp for nowait  
    for(int i=0;i<this->Nft/2;++i)
      {  // loop over octant kx>0. Half of what FFT gives
	if(i > 0)
	  {
	    xx_MAS = i * M_PI_rNft;
	    xx_MAS = sin(xx_MAS)/xx_MAS;
	  }
	else xx_MAS = 1.0;
	
        sn_xx_MAS = sin(i * M_PI_rNft);
	
	
	i_deltak_x_2 = i * i * this->deltak_x * this->deltak_x ;
	i_per_fact = i_deltak_x_2;
	if(my_s_box_ave == s_box_log)
	  i_deltak_x_2 *= rkmin * rkmin;
	
	for(int j=0;j<this->Nft/2;++j)
	  {  // loop over octant ky>0. Half of what FFT gives
	    yy_MAS = yy_MAS_array[j];
            sn_yy_MAS=sn_yy_MAS_array[j];

	    
	    j_deltak_y_2 = j * j * this->deltak_y * this->deltak_y;
	    j_per_fact = i_per_fact + j_deltak_y_2;
	    j_per_fact = sqrt(j_per_fact);
	    if(my_s_box_ave == s_box_log)
	      j_deltak_y_2 *= rkmin * rkmin;
	    
	    for(int k=0;k<=this->Nft/2;++k)
	      {  // loop over octant kz>0: this is all FFTW gives.
		zz_MAS = zz_MAS_array[k];
		sn_zz_MAS=sn_zz_MAS_array[k];
		
		k_deltak_z_2 = k * k * this->deltak_z * this->deltak_z;
		if(my_s_box_ave == s_box_log)
		  k_deltak_z_2 *= rkmin * rkmin;
		
                real_prec kv;
		int kmod_g, kmod_r;
		
		// Compute k-shell index 
		if(my_s_box_ave == s_box_linear)
		  {
		    kv=sqrt(i_deltak_x_2 + j_deltak_y_2 + k_deltak_z_2);
		    kmod_g=(int)floor((float)(kv * rDeltaK_data)+k_bin_index);
		    kmod_r=(int)floor((float)(kv * rDeltaK_window));
		  }
		else
		  {
		    if(my_s_box_ave == s_box_log)
		      {
			kv= i_deltak_x_2 + j_deltak_y_2 + k_deltak_z_2;
			if(kv!=0)
			  {
			    kmod_g=(int)floor((float)( log10(sqrt(kv))/this->Deltal))+1;
			    kmod_r=kmod_g;
			  }
			else{
			  kmod_g=0;kmod_r=0;
			}
		      }
		  }
		
		
		
		// Compute correction for MAS:
                real_prec corr=(s_box->use_MAS_correction ? my_correction_MAS(correction_MAS_exp, xx_MAS, yy_MAS, zz_MAS): 1.0);
                real_prec icorr2 = 1.0/(corr*corr);
		


                real_prec sn_corr=SN_aux*SN_correction_MAS(correction_MAS_exp, sn_xx_MAS, sn_yy_MAS, sn_zz_MAS);

		// Compute index in c-order
                int lp=index_3d(i,j,k,this->Nft,this->Nft/2+1);
		
		//Exclude the zero-frequency onLY when computing the P: for the W keep it:
		if(lp!=0)
		  {
                    //real_prec Pk=gsl_inter_new(ka,pa,kv);
                    real_prec Pk;
                    Pk= (this->data_out_g[lp][0] * this->data_out_g[lp][0] + this->data_out_g[lp][1] * this->data_out_g[lp][1] - sn_corr*this->normal_power) * icorr2 ;
		    
		    // 1d Spherical average
		    if(kmod_g<powerk_g0.size())
		      {
#pragma omp atomic update             
			powerk_g0[kmod_g]  += Pk;
#pragma omp atomic update       
			mod_g[kmod_g]++ ;    
		      }
		  }
		
		
		if(j>0  && k>0)
		  {    
                    lp=index_3d(i,this->Nft-j,k,this->Nft,this->Nft/2+1);   // Go to the ky<0 octant
		    
                    real_prec Pk;
                    Pk= (this->data_out_g[lp][0] * this->data_out_g[lp][0] + this->data_out_g[lp][1] * this->data_out_g[lp][1]- sn_corr*this->normal_power)* icorr2 ;
		    // 1d Spherical average
		    if(kmod_g<powerk_g0.size())
		      {
#pragma omp atomic  update                             
			powerk_g0[kmod_g]  += Pk;
#pragma omp atomic update                              
			mod_g[kmod_g]++ ;    
		      }
		  }
		
		
		// ******************************************
		// Add negative frequencies in x:
		// ******************************************
		if(i>0  && (j>0 || k>0))
		  {
                    lp=index_3d(this->Nft-i,j,k,this->Nft,this->Nft/2+1);
		    
                    real_prec Pk;
                    Pk= (this->data_out_g[lp][0] * this->data_out_g[lp][0] + this->data_out_g[lp][1] * this->data_out_g[lp][1]- sn_corr*this->normal_power) * icorr2 ;
		    
		    
		    if(kmod_g<powerk_g0.size())
		      {
#pragma omp atomic update                                 
			powerk_g0[kmod_g]  += Pk;
#pragma omp atomic update                              
			mod_g[kmod_g]++;
		      }
		  }
		
		// ******************************************
		//  Add negative frequencies in x and y:
		// ******************************************
		if(i>0  && j>0  && k>0)
		  {
                    lp=index_3d(this->Nft-i,this->Nft-j,k,this->Nft,this->Nft/2+1);
		    
                    real_prec Pk;
                    Pk= (this->data_out_g[lp][0] * this->data_out_g[lp][0] + this->data_out_g[lp][1] * this->data_out_g[lp][1] - sn_corr*this->normal_power) * icorr2;
		    
		    if(kmod_g<powerk_g0.size())
		      {
#pragma omp atomic update                                  
			powerk_g0[kmod_g] += Pk;
#pragma omp atomic update                                  
			mod_g[kmod_g] ++ ;    
		      }
		  }
	      }
	  }
      } // END of parallel loop
    
  } // END of parallel regione  
  
#pragma omp parallel for
  for(int i=0;i<powerk_g0.size();++i)
    powerk_g0[i]=(mod_g[i]== 0 ? 0 : powerk_g0[i]/(static_cast<real_prec>(mod_g[i])*this->normal_power));
  
}


// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************

// // USE THIS FUNCTION TO GET PK USING THE ROW-MAJOR ORDER
// // 
void FftwFunctions::power_spectrum_yamamoto(s_parameters_box *s_box,
                                            vector<real_prec> & powerk_g0,
                                            vector<real_prec> & powerk_g2,
                                            vector<real_prec> & powerk_g4,
					    vector<int> & mod_g)
{
	      
 // Power spectrum and multipole decomposition
  // Using the Yamamoto estimator with the algorithm of Scocimarro
  // that permits the use of FFTW.
  // Loops over all Available modes in Fourier space 
  
  int nbink=mod_g.size();
  for(int i=0;i<nbink;++i)mod_g[i]=0;
  for(int i=0;i<nbink;++i)powerk_g0[i]=0;
  for(int i=0;i<nbink;++i)powerk_g2[i]=0;
  for(int i=0;i<nbink;++i)powerk_g4[i]=0;
  
  int my_s_box_ave;
  real_prec rkmin, rDeltaK_data, rDeltaK_window;
  
#define s_box_linear 1
#define s_box_log 2
#define s_box_other 0
  
  my_s_box_ave = s_box_other;
  if(s_box->ave=="linear")
    my_s_box_ave = s_box_linear;
  else if(s_box->ave=="log")
    my_s_box_ave = s_box_log;
  
  
  rDeltaK_data = 1.0 / this->DeltaK_data;    
  rDeltaK_window = 1.0 / this->DeltaK_window;
  rkmin = 1.0/this->kmin;
  
  int NTHREADS = omp_get_max_threads();
  real_prec M_PI_rNft = M_PI * this->rNft;
  
  real_prec *yy_MAS_array, *zz_MAS_array;
  real_prec *_ky_array, *_ky2_array;
  yy_MAS_array = new real_prec[this->Nft+1];
  zz_MAS_array = new real_prec[this->Nft+1];
  _ky_array = new real_prec[this->Nft+1];
  _ky2_array = new real_prec[this->Nft+1];
  
  
  time_t start;
  time (&start);
  
  
  
  
  
#pragma omp parallel num_threads(NTHREADS)
  {
    real_prec _kx, _ky, kz;
    real_prec _kx2, _ky2, kz2;
    real_prec _kx_kmin2, _ky_kmin2;
    
    real_prec xx_MAS, yy_MAS, zz_MAS;
    real_prec rxx_MAS, ryy_MAS, rzz_MAS;
    
    int q_i, q_j, q_k;
    
#pragma omp sections
    {
#pragma omp section      
      
      for(int j = 0; j <= this->Nft; j++)
	{
	  q_j= (j<=this->Nft/2? j: j-(this->Nft+1)); 
	  _ky_array[j] = (_ky = q_j*this->deltak_y);
	  _ky2_array[j] = (_ky2 = _ky * _ky);
	}
      
#pragma omp section      
      for(int j = 1; j<= this->Nft; j++)
	{
	  q_j= (j<=this->Nft/2? j: j-(this->Nft+1)); 
	  yy_MAS = q_j * M_PI_rNft;
	  yy_MAS_array[j] = (yy_MAS = sin(yy_MAS)/yy_MAS);
	}

#pragma omp section
      for(int k = 1; k <= this->Nft; k++)
	{
	  q_k= (k<=this->Nft/2? k: k-(this->Nft+1)); 
	  zz_MAS = q_k * M_PI_rNft;
	  zz_MAS_array[k] = (zz_MAS = sin(zz_MAS)/zz_MAS);
	}
    }
    
#pragma omp single    
    yy_MAS_array[0] = 1.0;
    zz_MAS_array[0] = 1.0;
    
    
#pragma omp for nowait

    for(int i=0;i<=this->Nft;++i)
      {
	//set coordinate in units of the fundamental mode      
	q_i= (i<=this->Nft/2? i: i-(this->Nft+1)); 
	_kx = q_i*this->deltak_x;
	_kx2 = _kx * _kx; 
	if(my_s_box_ave == s_box_log)
	  _kx_kmin2 = _kx2 * rkmin * rkmin;
	if(i > 0) {
	  xx_MAS = q_i * M_PI_rNft;
	  xx_MAS = sin(xx_MAS)/xx_MAS;
	}
	else xx_MAS = 1.0;
	
	
	for(int j=0;j<=this->Nft;++j)
	  {
	    //set coordinate in units of the fundamental mode:
	    _ky = _ky_array[j];
	    _ky2 = _ky2_array[j];
	    yy_MAS = yy_MAS_array[j];
	    
	    if(my_s_box_ave == s_box_log)
	      _ky_kmin2 = _ky2 * rkmin * rkmin;
	    
	    
	    for(int k=0;k<=this->Nft;++k)
	      {
		//set coordinate in units of the fundamental mode:
		q_k= (k<=this->Nft/2? k: k-(this->Nft+1)); 
		
		zz_MAS = zz_MAS_array[k];
		
		
		// Note k>=Nft/2+1 is the region in which kz<0 (q_k = k-(Nft-1))
		// Then, in order to apply Hermitian symmetry
		// delta(kx,ky,-kz)=delta(-kx,-ky,kz)*
		// when we cover this region we reverse the sign of kx and ky
		// Since kz is reversed (see the definiton of q_k) in this region, we also reverse it here 
		// to make it positive. 
		// Note that the terms assocaited to the multipole decomposition
		// are products of k-coordinates: these products do not change
		// when we apply Hermitian symetri, so there is no need
		// to remap coordinates.
		
		
		int kmod_g, kmod_r;
		
		
		
		// Find the indices assocaited to the original fftw array
		int o_i, o_j, o_k;
                real_prec factor,factor2;
		remap(this->Nft, this->Nft, this->Nft, i, j, k, &o_i, &o_j, &o_k, &factor);
		factor2 = factor * factor;
		// Recall that the value of factor is +1 if kz>0 and
		// -1 if kz<0. We could multiply the imaginary part of the array
		// by this factor such that Hermitian symmetry follows.
		// However, note that factor comes squared, 
		// then will be always one.
		// The coordinates however need to be reversed in the kz<0 region, 
		// again, to apply HS. This because we need explicitely need the coordinates
		// in this estimator, not only the modulos of the wavevector.
		// This factor in the end appears on power of 4 though, i.e, 1. 
                real_prec kx = _kx * factor;
                real_prec ky = _ky * factor;
                real_prec kx2 = _kx2 * factor2;
                real_prec ky2 = _ky2 * factor2;
		
		kz = factor*q_k*this->deltak_z;
		kz2 = kz * kz;
		
                real_prec kv, kv2;
		if(my_s_box_ave == s_box_linear){          
		  kv2 = kx2 + ky2 + kz2;
		  kv = sqrt(kv2);
		  kmod_g=(int)floor((float)(kv * rDeltaK_data));
		  kmod_r=(int)floor((float)(kv * rDeltaK_window));
		}
		else{
		  if(my_s_box_ave == s_box_log)
		    {
                      real_prec kx_kmin2 = _kx_kmin2 * factor2;
                      real_prec ky_kmin2 = _ky_kmin2 * factor2;
		      
		      kv = kx_kmin2 + ky_kmin2 + kz2 * rkmin * rkmin;
		      if(kv!=0){
			kmod_g=(int)floor((float)( log10(sqrt(kv))/Deltal))+1;
			kmod_r=kmod_g;
		      }
		      else{kmod_g=0;kmod_r=0;}
		    }
		}
		
		
		// Compute the c-ordered index associated to the original fftw array:
		int lp=ijk(o_i,o_j,o_k,this->Nft,this->Nft,this->Nft/2+1);
		
		// Compute correction for MAS:
                real_prec corr=(s_box->use_MAS_correction ? my_correction_MAS(correction_MAS_exp, xx_MAS, yy_MAS, zz_MAS): 1.0);
                real_prec icorr2 = 1.0/(corr*corr);
		
		if(lp!=0)
		  { //Exclude zero-mode
		    
		    
		    // Subtract shot noise and normalize only the power spectrum.
                    real_prec SN_aux;
		    if(true==s_box->use_SN_correction)
		      SN_aux=shot_noise;
		    else
		      SN_aux=0;

		    
                    real_prec F2;
                    real_prec Pk=  (this->data_out_g[lp][0] * this->data_out_g[lp][0]+ this->data_out_g[lp][1]*this->data_out_g[lp][1]) * icorr2;
		    
		    
                    real_prec F2_r=0.5*3.*(  (kx2*this->data_out_g_xx[lp][0])
					  + (ky2*this->data_out_g_yy[lp][0]) 
					  + (kz2*this->data_out_g_zz[lp][0]) 
					  + (2.*kx*ky*this->data_out_g_xy[lp][0]) 
					  + (2.*kx*kz*this->data_out_g_xz[lp][0]) 
					  + (2.*ky*kz*this->data_out_g_yz[lp][0]))/kv2  - 0.5*this->data_out_g[lp][0];
		    
                    real_prec F2_i=0.5*3.*( (kx2*this->data_out_g_xx[lp][1])
					 + (ky2*this->data_out_g_yy[lp][1])
					 + (kz2*this->data_out_g_zz[lp][1])
					 + (2.*kx*ky*this->data_out_g_xy[lp][1]) 
					 + (2.*kx*kz*this->data_out_g_xz[lp][1]) 
					 + (2.*ky*kz*this->data_out_g_yz[lp][1]))/kv2  - 0.5*this->data_out_g[lp][1];
		    
                    real_prec fluc1=  (F2_r*this->data_out_g[lp][0] + F2_i*this->data_out_g[lp][1]) * icorr2;
		    
		    
		    if(statistics=="Pk_ys" || this->statistics=="Pk_ysc")
		      {
			F2= (F2_r * F2_r + F2_i*F2_i) * icorr2; //use this for Scoccimarro's version of Hexadecapole in the the Yamamoto-Blake estimator
		      }	  
		    
		    else if(statistics=="Pk_yb"  ||  statistics=="Pk_ybc")
		      {
                        real_prec F4_r = (35./8.)*((kx2*kx2)*this->data_out_g_xxx[lp][0]
						+(ky2*ky2)*this->data_out_g_yyy[lp][0]  
						+(kz2*kz2)*this->data_out_g_zzz[lp][0]  
						+4.*(kx2*kx)*ky*this->data_out_g_xxy[lp][0]
						+4.*(kx2*kx)*kz*this->data_out_g_xxz[lp][0]
						+4.*(ky2*ky)*kx*this->data_out_g_yyx[lp][0]
						+4.*(ky2*ky)*kz*this->data_out_g_yyz[lp][0]
						+4.*(kz2*kz)*kx*this->data_out_g_zzx[lp][0]
						+4.*(kz2*kz)*ky*this->data_out_g_zzy[lp][0]
						+6.*(kx2)*(ky2)*this->data_out_g_xyy[lp][0]
						+6.*(kx2)*(kz2)*this->data_out_g_xzz[lp][0]
						+6.*(ky2)*(kz2)*this->data_out_g_yzz[lp][0]
						+12.*kx*ky*kz*( kx* this->data_out_g_xyz[lp][0] + ky*this->data_out_g_yxz[lp][0]+kz*this->data_out_g_zxy[lp][0]))/(kv2*kv2)
			  -(15./6.)*F2_r-(7./8.)*this->data_out_g[lp][0];
			
			
                        real_prec F4_i = (35./8.)*((kx2*kx2)*this->data_out_g_xxx[lp][1]
						+(ky2*ky2)*this->data_out_g_yyy[lp][1]  
						+(kz2*kz2)*this->data_out_g_zzz[lp][1]  
						+4.*(kx2*kx)*ky* this->data_out_g_xxy[lp][1]
						+4.*(kx2*kx)*kz* this->data_out_g_xxz[lp][1]
						+4.*(ky2*ky)*kx* this->data_out_g_yyx[lp][1]
						+4.*(ky2*ky)*kz* this->data_out_g_yyz[lp][1]
						+4.*(kz2*kz)*kx* this->data_out_g_zzx[lp][1]
						+4.*(kz2*kz)*ky* this->data_out_g_zzy[lp][1]
						+6.*(kx2)*(ky2)*this->data_out_g_xyy[lp][1]
						+6.*(kx2)*(kz2)* this->data_out_g_xzz[lp][1]
						+6.*(ky2)*(kz2)* this->data_out_g_yzz[lp][1]
						+12.*kx*ky*kz*( kx* this->data_out_g_xyz[lp][1] + ky*this->data_out_g_yxz[lp][1]+kz*this->data_out_g_zxy[lp][1]))/(kv2*kv2)
			  -(15./6.)*F2_i-(7./8.)*this->data_out_g[lp][1];	  
			
			F2=  (F4_r*this->data_out_g[lp][0] + F4_i*this->data_out_g[lp][1])*icorr2; 
		      }
		    
		    
		    // 1d Spherical average:
		    if(kmod_g<nbink)
		      {
#pragma omp atomic              
			powerk_g0[kmod_g]  += Pk;
#pragma omp atomic              
			powerk_g2[kmod_g]  += fluc1;
#pragma omp atomic              
			powerk_g4[kmod_g]  += F2;
#pragma omp atomic              
			mod_g[kmod_g]++;    
		      }
		    
		  }
	      }
	  }
      }
  }    
  
  // Subtract shot noise and normalize only the power spectrum.
  real_prec SN_aux;
  if(true==s_box->use_SN_correction)SN_aux=shot_noise;
  else SN_aux=0;
  
  // For the monopole:
  for(int i=0;i<nbink;++i)
    powerk_g0[i]=(mod_g[i]== 0 ? 0 : powerk_g0[i]/((real_prec)mod_g[i]*this->normal_power)-SN_aux);
  
  // For the quadrupole: no shot noise correction applies
  for(int i=0;i<nbink;++i)
    powerk_g2[i]=(mod_g[i]== 0 ? 0 : 5.*powerk_g2[i]/((real_prec)mod_g[i]*this->normal_power));
  
  if(statistics=="Pk_ys" ||  statistics=="Pk_ysc"){
    // For the hexadecapole: no shot noise correction applies. Scoccimarro recipie.
    for(int i=0;i<nbink;++i)
      powerk_g4[i]=(mod_g[i]== 0 ? 0 : (35./2.)*powerk_g4[i]/((real_prec)mod_g[i]*normal_power)- powerk_g2[i] - (7./2.)*(powerk_g0[i]+SN_aux));
  }
  else if(statistics=="Pk_yb"||  statistics=="Pk_ybc" ){  
    //For the hexadecapole: no shot noise correction applies. Bianchi et al. recipie.
    for(int i=0;i<nbink;++i)
      powerk_g4[i]=(mod_g[i]== 0 ? 0 : 9.0*powerk_g4[i]/((real_prec)mod_g[i]*normal_power));
  }
}

// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************


void FftwFunctions::get_fkp_error_bars(s_parameters_box *s_box, s_data_structure *s_data, vector<real_prec> &kvector, vector<real_prec> &pk, vector<int> &modes, vector<real_prec> &sigma)
{
  //////////////////////////////////////////////////////////  
  // Estimates of the variance for the power spectrum
  // Based on the FKP estimator. Selects between the 
  // exact and the approximate expression for the var(P)
  //////////////////////////////////////////////////////////  


  if(true==this->use_random_catalog)
    {
      if(true==this->FKP_error_bars_exact)
	{//Use the exact FKP formula
	  for(int i=0;i<SN.size();++i)
            this->SN[i]=this->alpha*(1+this->alpha)*SN[i]/this->normal_power;
	  for(int i=0;i<Q.size(); ++i)
            Q[i] =this->alpha*Q[i]/this->normal_power;
          do_fkp_error_bars_exact(s_box, this->Q, this->SN, pk, sigma);
	  for(int i=0;i<sigma.size();++i)
	    sigma[i]=sqrt(2.*sigma[i])/(2.*modes[i]); 
	  // FKP Formula. Multiply the number of modes by a facgor 2 to account also
	  // for those kz<0 that were not explictely counted in the shell average
	}
      else
	{ //If not the exact, use the approximation for large scales, involving the definition of the effective volume
          vector<real_prec> veff(sigma.size(),0);
	  do_fkp_error_bars_veff(s_box,s_data, pk,veff);  
	  for(int i=0;i<veff.size();++i)
	    {
              real_prec deltaVk=4.*M_PI*this->DeltaK_data*pow(kvector[i],2)*(1.+(1./12.)*pow(this->DeltaK_data/kvector[i],2))/(pow(2.*M_PI, 3));
              sigma[i]=sqrt(2./(deltaVk*this->alpha*veff[i]));
	    }
	}      
    }
  else
    {
      for(int i=0;i<sigma.size();++i)
	{
          real_prec deltaVk=4.*M_PI*this->DeltaK_data*pow(kvector[i],2)*(1.+(1./12.)*pow(this->DeltaK_data/kvector[i],2))/(pow(2.*M_PI, 3));
	  sigma[i]=sqrt(2./(this->Lside*this->Lside*this->Lside*deltaVk))*(pk[i]+1./s_data->mean_density);
	  // This has been compared to the variance from N-Body simulations 
	  // giving a good agreement, specially on large scales.
	}
    }
}

// ********************************************************************************************************
// ********************************************************************************************************
// ********************************************************************************************************




void FftwFunctions::do_fkp_error_bars_exact(s_parameters_box *s_box, vector<real_prec>&Q, vector<real_prec> &SN, vector<real_prec> &P, vector<real_prec> &sigma)
{    
  //////////////////////////////////////////////////////////  
  // Estimation of the variance of the power spectrum using 
  // the FKP estimator with the exact expression, equation 
  // 2.4.6  of FKP paper
  //////////////////////////////////////////////////////////  
  
  real_prec temp1,temp2,perp1,perp2;
  real_prec kv, kvp;
  int ik, ikp, lp;

  /*Transform Q and SN*/

  this->data_out_SN=(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));
  this->data_out_Q =(complex_prec *)fftw_malloc(this->NTT*sizeof(real_prec));

  do_fftw_r2c(this->Nft,this->SN, this->data_out_SN);
  do_fftw_r2c(this->Nft,this->Q, this->data_out_Q);

  time_t start;
  time (&start);
  real_prec full=pow((real_prec)this->Nft/2, 2)*pow((real_prec)this->Nft/2, 2)*pow((real_prec)this->Nft/2+1, 2)-1;
  long counter=0;

  /*sum over k' */
  for(int i=0;i<this->Nft/2;++i){   
    for(int j=0;j<this->Nft/2;++j){ 
      for(int k=0;k<this->Nft/2+1;++k){       
	if(s_box->ave=="linear"){
	  kv=sqrt(pow(i*this->deltak_x,2)+pow(j*this->deltak_y,2)+pow(k*this->deltak_z,2));
	  ik=(int)floor((float)(kv/this->DeltaK_data));
	}
	else{
	  if(s_box->ave=="log"){
	    perp1=pow(this->deltak_x*i/(this->kmin),2)+pow(this->deltak_y*j/(this->kmin),2)+pow(this->deltak_z*k/(this->kmin),2);
	    if(perp1!=0){
	      ik= (int)floor( (float) (log10(sqrt(perp1))/this->Deltal))+1;
	    }
	    else ik=0;
	  }
	}
	
	/*sum over k'' */
	for(int l=0;l<this->Nft/2;l++){  
	  for(int m=0;m<this->Nft/2;m++){
	    for(int n=0;n<this->Nft/2;n++){
	      counter++;
	      //	      So.comp_time(start, full, counter);
	      
	      if(s_box->ave=="linear"){
		kvp=sqrt(pow(l*this->deltak_x,2)+pow(m*this->deltak_y,2)+pow(n*this->deltak_z,2));
		ikp=(int)floor((float)(kvp/this->DeltaK_data));
	      }
	      else{
		if(s_box->ave=="log"){
		  kvp=pow(this->deltak_x*l/(this->kmin),2)+pow(this->deltak_y*m/(this->kmin),2)+pow(this->deltak_z*n/(this->kmin),2);
		  if(kvp!=0){
		    ikp= (int)floor( (float) (log10(sqrt(perp1))/this->Deltal))+1;
		  }
		  else ikp=0;
		}
	      }
	      
	      if (ik == ikp){
		int I=fabs(l-i);
		int J=fabs(m-j);
		int K=fabs(n-k);
		lp=ijk(I,J,K,this->Nft,this->Nft,this->Nft/2+1);
		if(lp!=0){
		  temp1 = P[ik]*data_out_Q[lp][0]+this->data_out_SN[lp][0];     
		  temp2 = P[ik]*data_out_Q[lp][1]+data_out_SN[lp][1];     
		  sigma[ik] += (temp1*temp1+temp2*temp2);
		}
		
		if(J > 0 && K>0){
		  lp=ijk(I,this->Nft-J,K,this->Nft,this->Nft,this->Nft/2+1);
		  temp1 = P[ik]*data_out_Q[lp][0]+data_out_SN[lp][0];     
		  temp2 = P[ik]*data_out_Q[lp][1]+data_out_SN[lp][1];     
		  sigma[ik] += (temp1*temp1+temp2*temp2);
		}
		
		if(I > 0 && (J>0 || K>0)){
		  lp=ijk(this->Nft-I,J,K,this->Nft,this->Nft,this->Nft/2+1);
		  temp1 = P[ik]*data_out_Q[lp][0]+data_out_SN[lp][0];     
		  temp2 = P[ik]*data_out_Q[lp][1]+data_out_SN[lp][1];    
		  sigma[ik] += (temp1*temp1+temp2*temp2);
		}
		
		if(I > 0 && J > 0 && K>0){
		  lp=ijk(this->Nft-I,this->Nft-J,K,this->Nft,this->Nft,this->Nft/2+1);
		  temp1 = P[ik]*data_out_Q[lp][0]+data_out_SN[lp][0];     
		  temp2 = P[ik]*data_out_Q[lp][1]+data_out_SN[lp][1];     
		  sigma[ik] += (temp1*temp1+temp2*temp2);
		}
	      } 		
	    }
	  }
	}
      }
    }
  }
  fftw_free(data_out_SN);
  fftw_free(data_out_Q);
  return;
}


// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************


void FftwFunctions::do_fkp_error_bars_veff(s_parameters_box *s_box, s_data_structure *s_data_str, vector<real_prec> &P,vector<real_prec> &Veff)
{    
  
  //////////////////////////////////////////////////////////  
  // Computes the variance in the power spectrum as an 
  // approximation to the original FKP  estimation, introducing 
  // the definition of the effective volume. This is computed 
  // using the random catalogue.
  // Output: effective volume as a function of k:                      
  //////////////////////////////////////////////////////////  
  
  So.message_screen("Computing effetive volume");
  struct s_data_structure * s_data= (struct s_data_structure *)s_data_str;
//  vector< real_prec >prop = (s_data->properties);
//  int n_columns = s_data->n_columns;
  int n_objects = s_data->properties.size() ; // / n_columns;
  omp_set_num_threads(omp_get_max_threads());
#pragma omp parallel
  {
    real_prec *table_private;
    table_private = new real_prec[Veff.size()];
    memset(table_private, 0, sizeof(real_prec) * Veff.size());

#pragma omp for collapse(2)
    for(int i = 0; i < n_objects; ++i){
      for(int k=0;k<Veff.size(); ++k){
        real_prec nbar;
        nbar=s_data->properties[i].mean_density;
        table_private[k]+=nbar*1.0/((1+nbar*P[k])*(1+nbar*P[k]));  /*eff vol*/
      }
    }
#pragma omp critical    
    for(int k=0;k<Veff.size();++k)
      Veff[k]+=table_private[k];
  }
  return;
}


// ****************************************************************************************************



// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


// // Esta tiene los loops grid - galaxias
void FftwFunctions::get_power_moments_fourier_grid_ds_yam(s_parameters_box *s_box,
						   s_data_structure *s_data){
  
  // ***************************************************************************
  // Sampling of the elements used in the Yamamoto et al. estimator            *
  // for the moments of the 3d power spectrum.                                 *
  // We are generating an output with the same structure as that               * 
  // given by the FFTW algorithm, i.e, allocating a complex vector             *
  // with positive and negative frequencies except for the third component     *
  // which has only the positive quadrant. Furthermore, the first two          *
  // components only have the positive Nyquist frequency.                      *
  // The negative z components will be then                                    *
  // recovered using the Hermitian symmetry, as is done for the power spectrum *
  // computed using the FFTW.                                                  *
  // Here we are explicitely computing the multipoles l=0,2,4                  *
  // ***************************************************************************

  cout.precision(12);
 
  // **************************************************************************
  real_prec fac      = M_PI/180.0;
  real_prec area     = s_box->area_survey*pow(fac,2);        /*survey area in rad*/
  real_prec mean_den = s_data->mean_density;
  real_prec DELTAK   = this->DeltaK_data;

  // *****************************************************************************
  // Identify columns in the corresponding catalog                               *
  // *****************************************************************************
  int i_coord1      = (s_data->catalog=="data"? s_box->i_coord1_g : s_box->i_coord1_r);
  int i_coord2      = (s_data->catalog=="data"? s_box->i_coord2_g : s_box->i_coord2_r);
  int i_coord3      = (s_data->catalog=="data"? s_box->i_coord3_g : s_box->i_coord3_r);
  
  int i_weight1      = (s_data->catalog=="data"? s_box->i_weight1_g : s_box->i_weight1_r);
  int i_weight2      = (s_data->catalog=="data"? s_box->i_weight2_g : s_box->i_weight2_r);
  int i_weight3      = (s_data->catalog=="data"? s_box->i_weight3_g : s_box->i_weight3_r);
  int i_weight4      = (s_data->catalog=="data"? s_box->i_weight4_g : s_box->i_weight4_r);

  bool use_weight1   = (s_data->catalog=="data"? s_box->use_weight1_g : s_box->use_weight1_r);
  bool use_weight2   = (s_data->catalog=="data"? s_box->use_weight2_g : s_box->use_weight2_r);
  bool use_weight3   = (s_data->catalog=="data"? s_box->use_weight3_g : s_box->use_weight3_r);
  bool use_weight4   = (s_data->catalog=="data"? s_box->use_weight4_g : s_box->use_weight4_r);

  // *****************************************************************************
  // Arrays used in the interpolation of mean number density                     *
  // *****************************************************************************
  vector< vector<gsl_real> > dndz_m=(s_data->dndz_matrix);
  vector<gsl_real> zzv =  (s_data->zz_v);
  vector<gsl_real> dndz = (s_data->dndz_v);
  int i_mean_density= (s_data->catalog=="data"? s_box->i_mean_density_g : s_box->i_mean_density_r);
  string angles_units = (s_data->catalog=="data" ? s_box->angles_units_g : s_box->angles_units_r);

  // *****************************************************************************
  // Identify catalogue                                                          *
  // *****************************************************************************
//  vector < real_prec >prop = (s_data->properties);
  int nlines=s_data->properties.size();// / s_data->n_columns;
  int n_columns= s_data-> n_columns;
  vector<gsl_real> zz=s_box->zz;
  vector<gsl_real> rc=s_box->rc;

  // *****************************************************************************
  long full=((real_prec)this->Nft/2)*((real_prec)this->Nft/2)*((real_prec)this->Nft/2+1);
  long counter=0;
  
  
  // *****************************************************************************
  // Start loop over the Fourier grid and catalogues                             *
  // *****************************************************************************
  time_t start;
  time (&start);
  
  

  real_prec total_weight=1.0;
  int n_selected=0;
  real_prec normal_power=0;
  real_prec W_r=0;
  real_prec S_r_0=0;
  real_prec a11=0;
  real_prec a12=0;
  real_prec a13=0;
  real_prec a14=0;
  real_prec a15=0;
  real_prec a16=0;
  real_prec a21=0;
  real_prec a22=0;
  real_prec a23=0;
  real_prec a24=0;
  real_prec a25=0;
  real_prec a26=0;
  real_prec a31=0;
  real_prec a32=0;
  real_prec a33=0;
  real_prec a34=0;
  real_prec a35=0;
  real_prec a36=0;
  real_prec a41=0;
  real_prec a42=0;
  real_prec a43=0;
  real_prec a44=0;
  real_prec a45=0;
  real_prec a46=0;
  
  omp_set_num_threads(omp_get_max_threads());
  So.message_screen("Using", omp_get_max_threads(), "threads");

#pragma omp for
  // Loop over the grid

  for(int i=0;i<this->Nft/2;i++){
    int i2=i+this->Nft/2;
    int q_i2= (i==0? this->Nft/2: i-this->Nft/2);
    
    real_prec kx1=((real_prec)i)*deltak_x;
    real_prec kx2=q_i2*deltak_x;
    real_prec kx3=kx2;
    real_prec kx4=kx1;
    
    
    for(int j=0;j<this->Nft/2;j++){
      int j2=j+this->Nft/2;
      int q_j2=(j==0? this->Nft/2 : j-this->Nft/2);      
      
      real_prec ky1=((real_prec)j)*deltak_y;
      real_prec ky2=ky1;
      real_prec ky3=q_j2*deltak_y;
      real_prec ky4=ky3;
      
      for(int k=0;k<=this->Nft/2;k++){ 
	// counter++;
	// So.comp_time(start, full, counter);
        real_prec kz  = ((real_prec)k)*deltak_z;
	
        real_prec kk1=sqrt(kx1*kx1+ky1*ky1+kz*kz);   //Determine magnitude of k-vector
        real_prec kk2=sqrt(kx2*kx2+ky1*ky1+kz*kz);   //Determine magnitude of k-vector
        real_prec kk3=sqrt(kx2*kx2+ky3*ky3+kz*kz);   //Determine magnitude of k-vector
        real_prec kk4=sqrt(kx1*kx1+ky3*ky3+kz*kz);   //Determine magnitude of k-vector
	
	
	// Reseteando cada vex que inicia un loop sobre galaxias
	total_weight=1.0;
	n_selected=0;
	normal_power=0;
	W_r=0;   
	S_r_0=0; 
	a11=0;
	a12=0;
	a13=0;
	a14=0;
	a15=0;
	a16=0;
	a21=0;
	a22=0;
	a23=0;
	a24=0;
	a25=0;
	a26=0;
	a31=0;
	a32=0;
	a33=0;
	a34=0;
	a35=0;
	a36=0;
	a41=0;
	a42=0;
	a43=0;
	a44=0;
	a45=0;
	a46=0;
	
	//omp_set_num_threads(omp_get_max_threads());
	// std::cout<<RED<<"Warning: using one thread in Pk_y_ds"<<RESET<<std::endl;
	// omp_set_num_threads(1);
        real_prec *pprop;
        real_prec ow[4];
	
	for(int Ig=0;Ig<nlines;++Ig){
	  
//	  pprop = &prop[Ig*n_columns];
          real_prec we;
          real_prec x=s_data->properties[i].coord1;
          real_prec y=s_data->properties[i].coord2;
          real_prec z=s_data->properties[i].coord3;
          real_prec nbar=s_data->properties[i].mean_density;
	  
          // real_prec x=prop[I*n_columns+i_coord1];
          // real_prec y=prop[I*n_columns+i_coord2];
          // real_prec z=prop[I*n_columns+i_coord3];
          // real_prec nbar=prop[I*n_columns+i_mean_density];
	  
          real_prec rr= sqrt(x*x+y*y+z*z);
	  
	  // ***************************************************************
	  // Compute weights for each galaxy
	  // Compute the weights for each galaxy
	  
          ow[0] = (use_weight1 && (i_weight1<n_columns))? s_data->properties[i].weight1 : 1.0;
          ow[1] = (use_weight2 && (i_weight2<n_columns))? s_data->properties[i].weight2 : 1.0;
          ow[2] = (use_weight3 && (i_weight3<n_columns))? s_data->properties[i].weight3 : 1.0;
          ow[3] = (use_weight4 && (i_weight4<n_columns))? s_data->properties[i].weight4 : 1.0;
	  
	  total_weight = ow[0] * ow[1] * ow[2] * ow[3];
	  if(true==s_box->FKP_weight){we=1.0/(1+s_box->Pest*nbar);} else we=1.0;
	  total_weight*=we;      
          real_prec ptotal_weight2=total_weight*total_weight;
	  
	  
	  if(kk1!=0){
	    int lp1=ijk(i ,j ,k,this->Nft,this->Nft,this->Nft/2+1);//Find index as in FFTW
	    int lp2=ijk(i2,j ,k,this->Nft,this->Nft,this->Nft/2+1);//Find index as in FFTW
	    int lp3=ijk(i2,j2,k,this->Nft,this->Nft,this->Nft/2+1);//Find index as in FFTW
	    int lp4=ijk(i ,j2,k,this->Nft,this->Nft,this->Nft/2+1);//Find index as in FFTW
	    
	    // *******************************************************************
	    // Now compute r, k and mu and assign to the grid
	    
	    // first octant:
	    if(kk1<=this->kmax_y_ds)
	      {
                real_prec mu    = (x*kx1+y*ky1+z*kz)/(kk1*rr);   //cos of the angle between r and k
                real_prec cc    = total_weight*cos(kk1*rr*mu);
                real_prec ss    =-total_weight*sin(kk1*rr*mu);
                real_prec Leg_2 = 0.5*(3.*mu*mu-1.);
                real_prec Leg_4 = (1./8.)*(35*pow(mu,4)-30.*mu*mu+3);
		
		a11+=cc;
		a12+=ss;
		a13+=cc*Leg_2;
		a14+=ss*Leg_2;
		a15+=cc*Leg_4;
		a16+=ss*Leg_4;
		
		if(s_data->catalog=="data"){
		  this->data_g_out_y0[lp1].real(a11);
		  this->data_g_out_y0[lp1].imag(a12);
		  this->data_g_out_y2[lp1].real(a13);
		  this->data_g_out_y2[lp1].imag(a14);
		  this->SN_g_out_y2[lp1]        += pow(total_weight,2)*Leg_2;     //Shot-noise for multipole
		  this->data_g_out_y4[lp1].real(a15);
		  this->data_g_out_y4[lp1].imag(a16);
		  this->SN_g_out_y4[lp1]        += pow(total_weight,2)*Leg_4;     //Shot-noise for multipoles
		}
		else{
		  this->data_r_out_y0[lp1].real(a11);
		  this->data_r_out_y0[lp1].imag(a12);
		  this->data_r_out_y2[lp1].real(a13);
		  this->data_r_out_y2[lp1].imag(a14);
		  this->SN_r_out_y2[lp1]        += pow(total_weight,2)*Leg_2;     //Shot-noise for multipole
		  this->data_r_out_y4[lp1].real(a15);
		  this->data_r_out_y4[lp1].imag(a16);
		  this->SN_r_out_y4[lp1]        += pow(total_weight,2)*Leg_4;     //Shot-noise for multipoles
		}
		
		
		
	      }
	    
	    
	    // // second octant:
	    if(kk2<=this->kmax_y_ds)
	      {
                real_prec mu    = (x*kx2+y*ky1+z*kz)/(kk2*rr);
                real_prec cc    = total_weight*cos(kk2*rr*mu);
                real_prec ss    =-total_weight*sin(kk2*rr*mu);
                real_prec Leg_2 = 0.5*(3.*mu*mu-1.);
                real_prec Leg_4 = (1./8.)*(35*pow(mu,4)-30.*mu*mu+3);
		a21+=cc;
		a22+=ss;
		a23+=cc*Leg_2;
		a24+=ss*Leg_2;
		a25+=cc*Leg_4;
		a26+=ss*Leg_4;
		
		if(s_data->catalog=="data"){
		  this->data_g_out_y0[lp2].real(a21);
		  this->data_g_out_y0[lp2].imag(a22);
		  this->data_g_out_y2[lp2].real(a23);
		  this->data_g_out_y2[lp2].imag(a24);
		  this->SN_g_out_y2[lp2]        += pow(total_weight,2)*Leg_2;     //Shot-noise for multipole
		  this->data_g_out_y4[lp2].real(a25);
		  this->data_g_out_y4[lp2].imag(a26);
		  this->SN_g_out_y4[lp2]        += pow(total_weight,2)*Leg_4;     //Shot-noise for multipoles
		}
		else{
		  this->data_r_out_y0[lp2].real(a21);
		  this->data_r_out_y0[lp2].imag(a22);
		  this->data_r_out_y2[lp2].real(a23);
		  this->data_r_out_y2[lp2].imag(a24);
		  this->SN_r_out_y2[lp2]        += pow(total_weight,2)*Leg_2;     //Shot-noise for multipole
		  this->data_r_out_y4[lp2].real(a25);
		  this->data_r_out_y4[lp2].imag(a26);
		  this->SN_r_out_y4[lp2]        += pow(total_weight,2)*Leg_4;     //Shot-noise for multipoles
		}
		
	      }
	    
	    
	    // Third octant
	    if(kk3<=this->kmax_y_ds)
	      {
                real_prec mu    = (x*kx2+y*ky3+z*kz)/(kk3*rr);
                real_prec cc    = total_weight*cos(kk3*rr*mu);
                real_prec ss    =-total_weight*sin(kk3*rr*mu);
                real_prec Leg_2 = 0.5*(3.*mu*mu-1.);
                real_prec Leg_4 = (1./8.)*(35*pow(mu,4)-30.*mu*mu+3);
		a31+=cc;
		a32+=ss;
		a33+=cc*Leg_2;
		a34+=ss*Leg_2;
		a35+=cc*Leg_4;
		a36+=ss*Leg_4;
		
		if(s_data->catalog=="data"){
		  this->data_g_out_y0[lp3].real(a31);
		  this->data_g_out_y0[lp3].imag(a32);
		  this->data_g_out_y2[lp3].real(a33);
		  this->data_g_out_y2[lp3].imag(a34);
		  this->SN_g_out_y2[lp3]        += pow(total_weight,2)*Leg_2;     //Shot-noise for multipole
		  this->data_g_out_y4[lp3].real(a35);
		  this->data_g_out_y4[lp3].imag(a36);
		  this->SN_g_out_y4[lp3]        += pow(total_weight,2)*Leg_4;     //Shot-noise for multipoles
		}
		else{
		  this->data_r_out_y0[lp3].real(a31);
		  this->data_r_out_y0[lp3].imag(a32);
		  this->data_r_out_y2[lp3].real(a33);
		  this->data_r_out_y2[lp3].imag(a34);
		  this->SN_r_out_y2[lp3]        += pow(total_weight,2)*Leg_2;     //Shot-noise for multipole
		  this->data_r_out_y4[lp3].real(a35);
		  this->data_r_out_y4[lp3].imag(a36);
		  this->SN_r_out_y4[lp3]        += pow(total_weight,2)*Leg_4;     //Shot-noise for multipoles
		}
		
		
	      }
	    
	    if(kk4<=this->kmax_y_ds)
	      {
		// fourth octant:
                real_prec mu    = (x*kx1+y*ky3+z*kz)/(kk4*rr);   //cos of the angle between r and k
                real_prec cc    = total_weight*cos(kk4*rr*mu);
                real_prec ss    =-total_weight*sin(kk4*rr*mu);
                real_prec Leg_2 = 0.5*(3.*mu*mu-1.);
                real_prec Leg_4 = (1./8.)*(35*pow(mu,4)-30.*mu*mu+3);
		a41+=cc;
		a42+=ss;
		a43+=cc*Leg_2;
		a44+=ss*Leg_2;
		a45+=cc*Leg_4;
		a46+=ss*Leg_4;
		
		if(s_data->catalog=="data"){
		  this->data_g_out_y0[lp4].real(a41);
		  this->data_g_out_y0[lp4].imag(a42);
		  this->data_g_out_y2[lp4].real(a43);
		  this->data_g_out_y2[lp4].imag(a44);
		  this->SN_g_out_y2[lp4]        += pow(total_weight,2)*Leg_2;     //Shot-noise for multipole
		  this->data_g_out_y4[lp4].real(a45);
		  this->data_g_out_y4[lp4].imag(a46);
		  this->SN_g_out_y4[lp4]        += pow(total_weight,2)*Leg_4;     //Shot-noise for multipoles
		}
		else{
		  this->data_r_out_y0[lp4].real(a41);
		  this->data_r_out_y0[lp4].imag(a42);
		  this->data_r_out_y2[lp4].real(a43);
		  this->data_r_out_y2[lp4].imag(a44);
		  this->SN_r_out_y2[lp4]        += pow(total_weight,2)*Leg_2;     //Shot-noise for multipole
		  this->data_r_out_y4[lp4].real(a45);
		  this->data_r_out_y4[lp4].imag(a46);
		  this->SN_r_out_y4[lp4]        += pow(total_weight,2)*Leg_4;     //Shot-noise for multipoles
		}
		
		
		
	      }
	    // ***********************************************************
	    // Parameters for the Pk
	    n_selected++;
	    W_r          +=total_weight;    //weighted number of selected objects
	    S_r_0        +=ptotal_weight2;                  //Shot-noise for multipoles
	    normal_power +=nbar*ptotal_weight2;            //normalization
	    // **********************************************************************
	    
	  }
	}
      }
    }
  }
  
  
  if(s_data->catalog=="data"){
    this->n_gal=n_selected;
    this->w_g=W_r;  
    this->s_g=S_r_0;
  }
  else{
    if(s_data->catalog=="random"){ 
      this->n_ran=n_selected;
      this->normal_p=normal_power;
      this->s_r=S_r_0;
      this->w_r=W_r;
    }
  }
  
  
  time_t end; time (&end);
  real_prec diff = difftime(end,start);
  cout<<BOLDMAGENTA;
  if (diff<60) cout <<"Lapse: "<<diff<<" secs"<<endl;
  else if (diff<3600) cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
  else cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
  cout<<RESET;
  
    
  return;
}






// ****************************************************************************************************************************
// ******************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************

void FftwFunctions::power_yam_1d_ds(s_parameters_box *s_box, vector<real_prec>&powerk0, vector<real_prec>&powerk2, vector<real_prec>&powerk4, vector<int> &mod){
  // ******************************************************************* 
  // Shell average in Fourier space.                                   *
  // Done in the same fashion as the one used in the FFTW. Indeed this *
  // is the same method as that member of the class FFTW_FUNCTION      *
  // *******************************************************************
  
  for(int i=0;i<mod.size();i++)mod[i]=0;
  for(int i=0;i<powerk0.size();i++)powerk0[i]=0;
  for(int i=0;i<powerk2.size();i++)powerk2[i]=0;
  for(int i=0;i<powerk4.size();i++)powerk4[i]=0;
  real_prec DELTAK=DeltaK_data;
  real_prec k_bin_index=(s_box->k_bin_step == 1.? 0.5:0.0);

  for(int i=0;i<=this->Nft;i++){
    int q_i= (i<=this->Nft/2? i: i-(this->Nft+1)); //set coordinate in units of the fundamental mode
    for(int j=0;j<=this->Nft;j++){
      int q_j= (j<=this->Nft/2? j: j-(this->Nft+1));  //set coordinate in units of the fundamental mode
      for(int k=0;k<=this->Nft;k++){
	int q_k= (k<=this->Nft/2? k: k-(this->Nft+1)); //set coordinate in units of the fundamental mode
	
        real_prec kv; int kmod;
	if(s_box->ave=="linear"){
	  kv=sqrt(pow(q_i*this->deltak_x,2)+pow(q_j*this->deltak_y,2)+pow(q_k*this->deltak_z,2));
	  kmod=(int)floor((float)(kv/DELTAK+k_bin_index));
	}
	else{
	  if(s_box->ave=="log"){
	    kv=pow(this->deltak_x*q_i/(this->kmin),2)+pow(deltak_y*q_j/(this->kmin),2)+pow(deltak_z*q_k/(this->kmin),2);
	    if(kv!=0){
	      kmod=(int)floor((float)( log10(sqrt(kv))/this->Deltal))+1;
	    }
	    else{kmod=0;}
	  }
	}
	
	if(kv<=this->kmax_y_ds){
	  int o_i, o_j, o_k;
          real_prec factor;
	  remap(this->Nft, this->Nft, this->Nft, i, j, k, &o_i, &o_j, &o_k, &factor);
	  int lp=ijk(o_i,o_j,o_k,this->Nft,this->Nft,this->Nft/2+1);
	  if(lp!=0){ 
	    if(kmod<powerk0.size()){
	      powerk0[kmod]+=this->data_g_y0[lp];
	      powerk2[kmod]+=this->data_g_y2[lp];
	      powerk4[kmod]+=this->data_g_y4[lp];
	      mod[kmod]++ ;    
	    }
	  }
	}
      }
    }
  }


  //  for(int i=0;i<mod.size();i++)cout<<this->statistics<<"  "<<powerk0[i]<<"  "<<mod[i]<<endl;

  for(int i=0;i<powerk0.size();i++)powerk0[i]=(mod[i]==0 ? 0 : powerk0[i]/mod[i]);
  for(int i=0;i<powerk2.size();i++)powerk2[i]=(mod[i]==0 ? 0 : powerk2[i]/mod[i]);
  for(int i=0;i<powerk4.size();i++)powerk4[i]=(mod[i]==0 ? 0 : powerk4[i]/mod[i]);


} 




// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// BISPECTRUM. TWO APPROACHES ARE IMPLEMENTED. BRUTE FORCE APPROACH AND ONE DESIGNED BY
// JENNIFER POLLAC, EMILIANO SEFUSATTI, WHICH IS FASTER (WAY FASTER)
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************

#define bis_idx(i, j, k) ((i)*this->Nft*this->Nft + (j)*this->Nft + k)

void FftwFunctions::bispectrum_fkp(char d_w, s_parameters_box *s_box, vector<real_prec> &bispect, vector<real_prec> &sn_bispect, vector<real_prec> &p1, vector<int> &mod)
{
  //////////////////////////////////////////////////////////  
  // Count of triangles in Fourier space
  // There is something wrong with the paralelization
  // in this function, Therefore it is commented.
  //////////////////////////////////////////////////////////  
  
  
  real_prec Normal_delta= sqrt(normal_power);
  vector<int> modp(this->Nft,0);
  int nnp= (d_w == 'd' ? this->Nnp_data : this->Nnp_window);
  
  memset(&mod.front(), 0, sizeof(int) * this->Nft * this->Nft * this->Nft);
  memset(&bispect.front(), 0, sizeof(real_prec) * this->Nft * this->Nft * this->Nft);

  real_prec DELTAK=(d_w == 'd' ? this->DeltaK_data : this->DeltaK_window);
  int ik, n_kmod, n_kmod2, n_kmod3;

  real_prec rDeltaK_data;
  real_prec rDeltaK_window;
  real_prec rkmin;
    
  int my_s_box_ave;
  real_prec k_bin_index=(s_box->k_bin_step == 1.? 0.5:0.0);

#define s_box_linear 1
#define s_box_log 2
#define s_box_other 0
  
  my_s_box_ave = s_box_other;
  if(s_box->ave=="linear")
    my_s_box_ave = s_box_linear;
  else if(s_box->ave=="log")
    my_s_box_ave = s_box_log;
    
  if(my_s_box_ave == s_box_linear) {
    rDeltaK_data = 1.0 / this->DeltaK_data;    
    rDeltaK_window = 1.0 / this->DeltaK_window;
  }
  else if(my_s_box_ave == s_box_log)
    rkmin = 1.0/this->kmin;
  
  time_t start;
  time (&start);
  real_prec full=pow((real_prec)this->Nft+1, 2)*pow((real_prec)this->Nft+1, 2)*pow((real_prec)this->Nft+1, 2)-1;
  long counter=0;
  
  
 // #pragma omp parallel
 //   {
  
  real_prec i_deltak_x_2, j_deltak_y_2, k_deltak_z_2;
  real_prec xx_MAS, yy_MAS, zz_MAS;
  
  
  //  First loop k1
  //#pragma omp for nowait
  
  for(int i=0;i<=this->Nft;++i){
    int q_i= (i<=this->Nft/2? i: i-(this->Nft+1)); //set coordinate in units of the fundamental mode
    if(q_i > 0) {
      xx_MAS = q_i * M_PI / this->Nft;
      xx_MAS = sin(xx_MAS)/xx_MAS;
    } else xx_MAS = 1.0;
    i_deltak_x_2 = q_i * q_i * this->deltak_x * this->deltak_x ;
    if(my_s_box_ave == s_box_log)
      i_deltak_x_2 *= rkmin * rkmin;
    
    
    for(int j=0;j<=this->Nft;++j){
      int q_j= (j<=this->Nft/2? j: j-(this->Nft+1)); //set coordinate in units of the fundamental mode
      if(q_j > 0) {
	yy_MAS = q_j * M_PI / this->Nft;
	yy_MAS = sin(yy_MAS)/yy_MAS;
      } else yy_MAS = 1.0;
      j_deltak_y_2 = q_j * q_j * this->deltak_y * this->deltak_y ;
      if(my_s_box_ave == s_box_log)
	j_deltak_y_2 *= rkmin * rkmin;
      
      for(int k=0;k<=this->Nft;++k){ 
	int q_k= (k<=this->Nft/2? k: k-(this->Nft+1)); //set coordinate in units of the fundamental mode
	if(q_k > 0) { 
	  zz_MAS = q_k * M_PI / this->Nft;
	  zz_MAS = sin(zz_MAS)/zz_MAS;
	} else zz_MAS = 1.0;
	k_deltak_z_2 = q_k * q_k * this->deltak_z * this->deltak_z ;
	if(my_s_box_ave == s_box_log)
	  k_deltak_z_2 *= rkmin * rkmin;
        
        real_prec kv, kvv, kvvv;
	int kmod, kmod2, kmod3;
	if(s_box->ave=="linear"){
	  kv=sqrt(i_deltak_x_2 + j_deltak_y_2 + k_deltak_z_2);
	  kmod=(int)floor((float)(kv/DELTAK+k_bin_index)); 
	}
	else{
	  if(s_box->ave=="log"){
	    kv=(i_deltak_x_2 + j_deltak_y_2 + k_deltak_z_2);
	    kmod=(kv==0? 0 : (int)floor((float)( log10(sqrt(kv))/Deltal))+1);
	  }
	}
	
        real_prec corr1 = (MAS_correction ? my_correction_MAS(correction_MAS_exp, xx_MAS, yy_MAS, zz_MAS): 1.00);
	int o_i, o_j, o_k;
        real_prec factor1;
	remap(this->Nft, this->Nft, this->Nft, i, j, k, &o_i, &o_j, &o_k, &factor1);
	
	int lp=ijk(o_i,o_j,o_k,this->Nft,this->Nft,this->Nft/2+1);
        real_prec d1r   = this->data_out_g[lp][0] / corr1;
	
        real_prec d1i   = factor1*this->data_out_g[lp][1] / corr1;
	
        real_prec p1s = pow(d1i,2)+pow(d1r,2);
	
	if(lp!=0){
	  //	  #pragma omp atomic            
	  p1[kmod] += d1r*d1r + d1i*d1i;  //Compute power spectrum
	  
	  //#pragma omp atomic            
	  modp[kmod]++;
	}	    
	
	// *********************************************************************************	    
	// Second loop k2 
	
        real_prec ii_deltak_x_2, jj_deltak_y_2, kk_deltak_z_2;
        real_prec xxx_MAS, yyy_MAS, zzz_MAS;
	
        
	for(int ii=0;ii<=this->Nft;++ii){
	  int q_ii=(ii<=this->Nft/2? ii: ii-(this->Nft+1));
	  if(q_ii > 0) {
	    xxx_MAS = q_ii * M_PI / this->Nft;
	    xxx_MAS = sin(xxx_MAS)/xxx_MAS;
	  } else xxx_MAS = 1.0;
	  ii_deltak_x_2 = q_ii * q_ii * this->deltak_x * this->deltak_x ;
	  if(my_s_box_ave == s_box_log)
	    ii_deltak_x_2 *= rkmin * rkmin;
          
	  
	  for(int jj=0;jj<=this->Nft;++jj){
	    int q_jj=(jj<=this->Nft/2? jj: jj-(this->Nft+1)); 
	    if(q_jj > 0) {
	      yyy_MAS = q_jj * M_PI / this->Nft;
	      yyy_MAS = sin(yyy_MAS)/yyy_MAS;
	    } else yyy_MAS = 1.0;
	    jj_deltak_y_2 = q_jj * q_jj * this->deltak_y * this->deltak_y ;
	    if(my_s_box_ave == s_box_log)
	      jj_deltak_y_2 *= rkmin * rkmin;
	    
            
	    for(int kk=0;kk<=this->Nft;++kk){
	      int q_kk=(kk<=this->Nft/2? kk: kk-(this->Nft+1)); 
	      if(q_kk > 0) {
		zzz_MAS = q_kk * M_PI / this->Nft;
		zzz_MAS = sin(zzz_MAS)/zzz_MAS;
	      } else zzz_MAS = 1.0;
	      kk_deltak_z_2 = q_kk * q_kk * this->deltak_z * this->deltak_z ;
	      if(my_s_box_ave == s_box_log)
		kk_deltak_z_2 *= rkmin * rkmin;
	      
	      
	      counter ++;
	      //	      So.comp_time(start, full, counter);
	      
	      if(my_s_box_ave==s_box_linear){
		kvv=sqrt(ii_deltak_x_2 + jj_deltak_y_2 + kk_deltak_z_2);
		kmod2=(int)floor((float)(kvv/DELTAK+k_bin_index));
	      }
	      else{
		if(my_s_box_ave==s_box_log){
		  kvv=(ii_deltak_x_2 + jj_deltak_y_2 + kk_deltak_z_2);
		  kmod2=(kvv==0? 0 : (int)floor((float)( log10(sqrt(kvv))/Deltal))+1);
		}
	      }	
	      int o_ii, o_jj, o_kk;	      
              real_prec factor2;
	      remap(this->Nft, this->Nft, this->Nft, ii, jj, kk, &o_ii, &o_jj, &o_kk, &factor2);
	      int lp2=ijk(o_ii,o_jj,o_kk,this->Nft,this->Nft,this->Nft/2+1);  
              real_prec corr2 = (MAS_correction ? my_correction_MAS(correction_MAS_exp, xxx_MAS, yyy_MAS, zzz_MAS): 1.000);
              real_prec d2r   = this->data_out_g[lp2][0]/corr2;
              real_prec d2i   = factor2*this->data_out_g[lp2][1]/corr2;

              real_prec p2s = pow(d2i,2)+pow(d2r,2);
	      
	      // Set the components of the vector k3 in units of DELTAk from the components
	      // of the vectors k1 and k2 to fullfil k3=-k1-k2
	      int q_iii=-q_i-q_ii;
	      int q_jjj=-q_j-q_jj;
	      int q_kkk=-q_k-q_kk;
	      
	      // Now proceed if the triangle is in the grid:
	      // the components of the third vector are required to be in the box, 
	      // that is, 0<= |q_i|<=n/2		      
	      if(fabs(q_iii)<=this->Nft/2 && fabs(q_jjj)<=this->Nft/2 && fabs(q_kkk)<=this->Nft/2){
		
		if(my_s_box_ave== s_box_linear){
		  kvvv=sqrt(q_iii*q_iii*deltak_x*deltak_x +
			    q_jjj*q_jjj*deltak_y*deltak_y +
			    q_kkk*q_kkk*deltak_z*deltak_z);
		  //kmod3=(int)floor((float)(kvvv/DELTAK)); //OFFICIAL BINNING, RESTORE IT WHEN COMPARISON SUCCEDS
		  kmod3=(int)floor((float)(kvvv/DELTAK+k_bin_index)); //JENN AND CRIS BINNING, 
		}
		else{
		  if(my_s_box_ave== s_box_log){
		    kvvv = (q_iii*q_iii*deltak_x*deltak_x +
			    q_jjj*q_jjj*deltak_y*deltak_y +
			    q_kkk*q_kkk*deltak_z*deltak_z) * rkmin*rkmin;
		    kmod3=(kvvv==0? 0 : (int)floor((float)( log10(sqrt(kvv))/Deltal))+1);
		  }
		}
		
		// Find the position in the grid from the components of the k3 vector
		int iii = (q_iii>=0 ? q_iii: q_iii+this->Nft+1);
		int jjj = (q_jjj>=0 ? q_jjj: q_jjj+this->Nft+1);
		int kkk = (q_kkk>=0 ? q_kkk: q_kkk+this->Nft+1);
		
		int o_iii, o_jjj, o_kkk;
                real_prec factor3;
		remap(this->Nft, this->Nft, this->Nft, iii, jjj, kkk, &o_iii, &o_jjj, &o_kkk, &factor3);
		int lp3=ijk(o_iii,o_jjj,o_kkk,this->Nft,this->Nft,this->Nft/2+1);
		
		if(lp!=0 &&  lp2!=0 && lp3!=0){			  // Exclude the zero-mode 
                  real_prec corr3=(MAS_correction ? correction_MAS(q_iii,q_jjj,q_kkk): 1.000);
                  real_prec d3r = this->data_out_g[lp3][0]/corr3;
                  real_prec d3i = factor3*this->data_out_g[lp3][1]/corr3;

                  real_prec p3s = pow(d3i,2)+pow(d3r,2);

		  sort_index(kmod, kmod2, kmod3, &n_kmod, &n_kmod2, &n_kmod3);	//Sort
		  
		  //#pragma omp atomic
		  ULONG ind=index_3d(n_kmod, n_kmod2, n_kmod3,this->Nft, this->Nft);
		  //		  mod[bis_idx(n_kmod, n_kmod2, n_kmod3)]++ ;
		  mod[ind]++ ;    		      
		  bispect[ind] +=(d1r*d2r*d3r- d1r*d2i*d3i- d1i*d2r*d3i- d1i*d2i*d3r);
		  sn_bispect[ind] += p1s+p2s+p3s;
		  // ************ Uncomment this to check the imaginary part ********************************//
		  //bispect[n_kmod][n_kmod2][n_kmod3] += (d1r*d2r*d3i + d1r*d2i*d3r + d1i*d2r*d3r - d1i*d2i*d3i); 
		  
		}
	      }
	    }
	  }
	}
      }
    }
    //  } // END of outern parallel for
  } // END of parallel region
  
  //Compute average power in the spherical shells
  for(int i=0;i<this->Nft;++i)p1[i]= (modp[i]==0? 0. : p1[i]/modp[i]); 
  return;
}


// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************

// THIS FUNCTION IS NO LONGER CALLED
void FftwFunctions::get_power_spectrum_for_bispectrum(s_parameters_box *s_box, vector<real_prec> &power_g0)
{
  //////////////////////////////////////////////////////////  
  // Computes the outputs for the P(k)
  // In order to generate the estimates of shot noise for 
  // bispectrum
  //////////////////////////////////////////////////////////  
  
  So.message_screen("Power spectrum for Bispectrum");

  time_t start;
  time (&start);
  

  do_fftw_r2c(this->Nft,this->data_g, this->data_out_g);

  std::cout<<RED;
  time_t end;
  time(&end);
  real_prec diff=difftime(end,start);
  if (diff<60.) std::cout <<"Lapse:  "<<diff<<" secs"<<endl;
  else if (diff<3600) std::cout <<"Lapse: "<<diff/60.0<<" minutes"<<endl;
  else std::cout<<"Lapse: "<<diff/3600.0<<" hours"<<endl;
  std::cout<<RESET;
  
 
  time (&start); 
  std::cout<<CYAN<<"Shell average..."<<RESET<<endl;

  power_spectrum_fkp_for_bispectrum(s_box,power_g0); 
  std::cout<<RED;
  time(&end);
  diff=difftime(end,start);
  if (diff<60.) std::cout <<"Lapse:  "<<diff<<" secs"<<endl;
  else if (diff<3600.) std::cout <<"Lapse: "<<diff/60.<<" minutes"<<endl;
  else std::cout<<"Lapse: "<<diff/3600.<<" hours"<<endl;
  std::cout<<RESET;



}


// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************

void FftwFunctions::get_bispectrum_fkp(char d_w, s_parameters_box *s_box, vector<real_prec> &bispect, vector<real_prec> &sn_bispect, vector< int > &mod)
{

  // *******************************************************
  // This function generates the estimates of Bispectrum   *
  // based on the counts of triangles                      *
  // the normalization and the shot-noise corrections      *
  // *******************************************************
 
  do_fftw_r2c(this->Nft,this->data_g, this->data_out_g);
  
  vector<real_prec> power(this->Nft,0);
  
  memset(&mod.front(), 0, sizeof(int) * this->Nft * this->Nft * this->Nft);
  memset(&bispect.front(), 0, sizeof(real_prec) * this->Nft * this->Nft * this->Nft);
  memset(&sn_bispect.front(), 0, sizeof(real_prec) * this->Nft * this->Nft * this->Nft);

  std::cout << CYAN << "Measuring bispectrum using brute-force approach ..." << RESET << endl;
  time_t start, stop;
  time (&start);
  bispectrum_fkp(d_w, s_box,bispect,sn_bispect, power, mod);
  time (&stop);
  std::cout << "done in " << difftime(stop,start) << " seconds" << endl;
  
  
  real_prec SN_aux;
  if(s_box->use_SN_correction)
    SN_aux=1.0;
  else
    SN_aux=0;
  
  //Estimates of the power spectrum using FKP, as done in the Pk main code.
  for(int i=0;i<power.size();++i)power[i]=power[i]/this->normal_power-shot_noise*SN_aux;
  
#pragma omp parallel
  {
    ULONG I;
#pragma omp for collapse(2) nowait
    for(int i=0;i<this->Nft;++i)
      for(int j=0;j<this->Nft;++j)
	for(int k=0;k<this->Nft;++k)
	  if(i<=j && j<=k) {  
	    I = index_3d(i,j,k, this->Nft, this->Nft);
	    //	    cout<<i<<"  "<<j<<"  "<<k<<"   "<<bispect[I]<<"   "<<mod[I]<<endl;
            //	    bispect[I]=(mod[I]== 0 ? 0 : bispect[I]/((real_prec)mod[I]*this->normal_b) - SN_aux*( (power[i]+power[j]+power[k])*shot_noise_b1+shot_noise_b2));
            bispect[I]=(mod[I]== 0 ? 0 : bispect[I]/((real_prec)mod[I]*this->normal_b) - SN_aux*( (sn_bispect[I]/mod[I])*shot_noise_b1+shot_noise_b2));
	  }
}
  
  fftw_free(data_out_g); 
  
  return ;
}
 
 
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// BISPECTRUM: TRANSLATION FROM THE CODE IN FORTRAN WRITTEN BY JENNIFER POLLACK
// ****************************************************************************************************************************
// ****************************************************************************************************************************

void FftwFunctions::define_kshells(s_parameters_box *s_box){
  
  
  std::cout<<BLUE<<"Defining shells:"<<RESET<<endl;
  time_t start;
  time (&start);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Getting ready for MAS correction
  real_prec xx_MAS, yy_MAS, zz_MAS; //for MAS correction
  real_prec *yy_MAS_array, *zz_MAS_array;
  real_prec *_ky_array, *_ky2_array;
  int q_i, q_j, q_k;

  real_prec k_bin_index=(s_box->k_bin_step == 1.? 0.5:0.0);

  yy_MAS_array = new real_prec[this->Nft+1];
  zz_MAS_array = new real_prec[this->Nft+1];
  real_prec M_PI_rNft = M_PI * this->rNft;


  yy_MAS_array[0] = 1.0;
  for(int j = 1; j<= this->Nft; j++) {
    q_j= (j<=this->Nft/2? j: j-(this->Nft+1)); 
    yy_MAS = q_j * M_PI_rNft;
    yy_MAS_array[j] = (yy_MAS = sin(yy_MAS)/yy_MAS);
  }
  
  zz_MAS_array[0] = 1.0;
  for(int k = 1; k <= this->Nft; k++) {
    q_k= (k<=this->Nft/2? k: k-(this->Nft+1)); 
    zz_MAS = q_k * M_PI_rNft;
    zz_MAS_array[k] = (zz_MAS = sin(zz_MAS)/zz_MAS);
  }
  if(this->MASS=="NGP") {
    MAS_ptr = &FftwFunctions::MAS_NGP;
    correction_MAS_exp = 0;
  }
  else if(this->MASS=="CIC") {
    MAS_ptr = &FftwFunctions::MAS_CIC;
    correction_MAS_exp = 1;
  }
  else if(this->MASS=="TSC") {
    MAS_ptr = &FftwFunctions::MAS_TSC;
    correction_MAS_exp = 2;
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  // Given input parameter of kmin, kmax and the number of shells, 
  // here the k-shells are defined

  int id=0;
  
  // This is a loop over the Fourier grid
  // that goes to the maximum value of k given as an
  // input, kmax_bk, in units of the fundamental mode. 
  // From this value, the sgrid is computed in the set_pars method of this class.
  
  // Viewed as an quarter, the loop if done over the modes below the diagonal.
  // The other modes can be mapped using some Symmetries, designed by Jennifer.

  

  for(int i=0;i<this->sgrid;++i){
    if(i > 0) {
      xx_MAS = i * M_PI_rNft;
      xx_MAS = sin(xx_MAS)/xx_MAS;
    } else xx_MAS = 1.0;



    for(int j=i;j<this->sgrid;++j){
      yy_MAS = yy_MAS_array[j];

      for(int k=j;k<this->sgrid;++k){
	zz_MAS = zz_MAS_array[k];

	// Compute the squared of the magnitude, in units of the fundamental mode.
        real_prec kk= i*i+j*j+k*k;

	// Exclude the zero*mode. Be consistent with line 75 (set_pars()) 
	// in order to get the right dimension of the arrays
	if(kk>0){
	  
	  // %%%%%%%%%%%%%%%%%%%%%%
	  // Assign the MAS correction for this vector
	  this->Array_corr[id]=(s_box->use_MAS_correction ? my_correction_MAS(correction_MAS_exp, xx_MAS, yy_MAS, zz_MAS): 1.000000);

	  // %%%%%%%%%%%%%%%%%%%%%%
	  // Assign to the mode id its x-component in the following array
	  this->Arraykx[id]=i;	
	  

	  // %%%%%%%%%%%%%%%%%%%%%%
	  // Assign to the mode id its y-component in the following array:
	  this->Arrayky[id]=j;
	  
	  // %%%%%%%%%%%%%%%%%%%%%%
	  // Assign to the mode id its z-component in the following array:
	  this->Arraykz[id]=k;
	  
	  
	  // %%%%%%%%%%%%%%%%%%%%%%
	  // Assign to the mode id its own original ID
	  this->ArrayID[id]=id;
	  
	  
	  // %%%%%%%%%%%%%%%%%%%%%%
	  // Assign to the mode id its magnitude squared in the following array:
	  this->Arraykk[id]=kk;
	  
	  // %%%%%%%%%%%%%%%%%%%%%%
	  // Symmetries. These numbers are sacred.
	  {
	    if(i!=0 && i!=j && j!=k)this->VecArray[id]=48;
	    else if(i!=0 && i!=j && j!=0 && j==k)this->VecArray[id]=25;
	    else if(i!=0 && i==j && j!=k)this->VecArray[id]=24;
	    else if(i==0 && i!=j && j!=k)this->VecArray[id]=23;
	    else if(i==0 && j!=0 && j==k)this->VecArray[id]=12;
	    else if(i!=0 && i==j && j==k)this->VecArray[id]=8;
	    else if(i==0 && j==0 && k!=0)this->VecArray[id]=6;
	    else if(i==0 && j==0 && k==0)this->VecArray[id]=0;
	  }
	  
	  // %%%%%%%%%%%%%%%%%%%%%%
	  
	  
	  int kbin = (kk==0? 0 : (int)floor((this->deltak_x*sqrt(static_cast<real_prec>(kk))/this->DeltaK_Bis)+k_bin_index));
	  
	  
	  // %%%%%%%%%%%%%%%%%%%%%%
	  //Assign to a mode identified with the id, the bin where it belongs to
	  this->Kbin[id] = kbin;
	  
	  // %%%%%%%%%%%%%%%%%%%%%%
	  // Count the number of modes in the bin. 
	  if(kbin<this->Nshells_bk)
	    this->Bmodes[kbin] ++;

	  // Count the number of modes 
	  id++;	  
	  
	}
      }  
    }
  }
  

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Define k bins
  for(int i=0;i<this->Nshells_bk;++i)
    this->kbins_bk[i]=(i+k_bin_index)*this->DeltaK_Bis;
  
  int nn_gg=0;
  
  for(int i=0;i<this->Nshells_bk;++i)
    {
      nn_gg=3*(int)((this->kbins_bk[i])/this->DeltaK_data);
      if(nn_gg%2!=0)nn_gg++;
      this->Ngrids_bk[i]=nn_gg;
    }  
  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Sort array with the squared of the magnitude of the wavevector
  // and the other vectors, shuffled accordingly:
  sort_vectors(this->Arraykk, this->Arraykx,this->Arrayky,this->Arraykz,this->Kbin, this->VecArray, this->ArrayID, this->Array_corr);
  
  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  // lopp over sorted indices according to the k squarde value
  // Find, for each shell, the ID associated to the
  // vector with the smallest value of k**2 in that particular shell
  // This can be done easily because the vectors are already sorted.
  int idd=0;
  ArrayID[0]=0;
  
  //#pragma omp parallel for
  for(int i=1;i<this->new_sgrid;++i)
    {  //start from 1 to evaluate i-1 without segfault. That's why we lieave the id++ here instead of writing it below
      idd++;
      this->ArrayID[i]=i; //resort, ask Jennifer
      if(Kbin[i]!=Kbin[i-1]){
	int bin = Kbin[i];
	if(bin<this->Nshells_bk){
	  this->kkminID[bin]=idd;
	}      
      }
    }  
  
  
  
  std::cout<<RED;
  time_t end; time (&end);
  real_prec diff = difftime(end,start);
  if (diff<60) std::cout <<"Lapse: "<<diff<<" secs"<<endl;
  else if (diff<3600) std::cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
  else std::cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
  std::cout<<RESET;
  
}  
 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
 
#define bis_idx_new(i, j, k) ((i)*this->Nshells_bk*this->Nshells_bk + (j)*this->Nshells_bk + k)

void FftwFunctions::get_ift_shells_bispectrum(s_parameters_box *s_box)
{
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  std::cout<<CYAN<<"Construct Inverse Fourier transforms in shells"<<RESET<<std::endl;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Construc the arrays with the inverse Fourier transforms in shells.
  // For each shell we allocate an array with dimension N*N*N, where N is the maximum
  // size of the Nft used, i.e, that of the very last shell.
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  time_t start;
  time (&start);
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  real_prec k_bin_index=(s_box->k_bin_step == 1.? 0.5:0.0);
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  int Nshell_min = (floor)( (this->kmin_bk/this->DeltaK_Bis) +k_bin_index); 
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  this->iFT_output_delta_full.resize(this->Nshells_bk);
  this->iFT_output_triangles_full.resize(this->Nshells_bk);
  this->iFT_shot_noise_p1_cyc_sum_full.resize(this->Nshells_bk);
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  int Ngrid_max=Ngrids_bk[this->Nshells_bk-1];
  ULONG Nt_max=(Ngrid_max*Ngrid_max)*Ngrid_max;

  // for(int io=0;io<3;++io)
  //   this->n2[io]=Ngrid_max;

  vector<real_prec> iFT_output_delta(Nt_max,0);
  vector<real_prec> iFT_output_triangles(Nt_max,0);
  vector<real_prec> iFT_output_p1_cyc_sum(Nt_max,0);
  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // LOOP OVER SHELLS:
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for(int i=0;i<this->Nshells_bk;++i){
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Remap arrays and generate the iFT
    construct_shells(Ngrid_max, this->kkminID[i],  this->kkminID[i] +this->Bmodes[i]-1, iFT_output_delta, iFT_output_triangles, iFT_output_p1_cyc_sum);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // For each shell allocate smae amount of memory
    // to the maximum, though loops will go as far 
    // as the inverse FT needs.

    // For deltas
    this->iFT_output_delta_full[i].resize(Nt_max,0);

    // For number of triangles
    this->iFT_output_triangles_full[i].resize(Nt_max,0);

    // For the shot-noise
    this->iFT_shot_noise_p1_cyc_sum_full[i].resize(Nt_max,0);

    // Reassign arrays to matrices
    for(int ik=0;ik<Nt_max;++ik)
      this->iFT_output_delta_full[i][ik] = iFT_output_delta[ik]/Ngrid_max;
    for(int ik=0;ik<Nt_max;++ik)
      this->iFT_output_triangles_full[i][ik]= iFT_output_triangles[ik]/Ngrid_max;
    for(int ik=0;ik<Nt_max;++ik)
      this->iFT_shot_noise_p1_cyc_sum_full[i][ik]= iFT_output_p1_cyc_sum[ik]/Ngrid_max;
    
  }
  
  iFT_output_delta.shrink_to_fit();
  iFT_output_triangles.shrink_to_fit();
  iFT_output_p1_cyc_sum.shrink_to_fit();
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  std::cout<<RED;
  time_t end; time (&end);
  real_prec diff = difftime(end,start);
  if (diff<60) std::cout <<"Lapse: "<<diff<<" secs"<<endl;
  else if (diff<3600) std::cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
  else std::cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
  std::cout<<RESET;
}
 

 

// **********************************************************************************************************************
// **********************************************************************************************************************
// New version: revision of the version written with Jennifer with an improvement in the time performance/.
// **********************************************************************************************************************
void FftwFunctions::loop_shells_bispectrum(s_parameters_box *s_box, vector<real_prec> &pk, vector<real_prec> &bispect, vector< int > &mod, string fname){
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Next step: loop over shells
  std::cout<<CYAN<<"Loop over shells:"<<RESET<<endl;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ofstream out; 
  out.open(fname.c_str() , ios::out); 
  out.precision(10); 
  out.setf(ios::showpoint); 
  out.setf(ios::scientific); 
  cout.precision(12);
  cout.setf(ios::scientific); 
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  real_prec k_bin_index=(s_box->k_bin_step == 1.? 0.5:0.0);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  time_t start;
  time (&start);
  //  ofstream oi; oi.open("test.dat");
  
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  int Nshell_min = (floor)( (this->kmin_bk/this->DeltaK_Bis) +k_bin_index); 
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  int Ngrid_max=Ngrids_bk[this->Nshells_bk-1];
  long Nt_max=(Ngrid_max*Ngrid_max)*Ngrid_max;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  real_prec SN_aux;
  if(s_box->use_SN_correction)SN_aux=1.0;
  else SN_aux=0;


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for(int i=Nshell_min;i<this->Nshells_bk;++i)
    {  
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      int nstep2 = floor(i/2);  //start at this position for the second loop
      int jstart=(nstep2==0 ? i: nstep2);  
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      for(int j=jstart;j<=i;++j)
	{    
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  // Prepare for loop over k:
	  int nstep3=i-j-1;
	  int kstart=(nstep3>j ? j: nstep3);
	  if(nstep3<=0)
	    kstart=Nshell_min;
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  for(int k=kstart;k<=j;++k)
	    {
	      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	      // Classify triangle bin as type [1] equilateral, [2]=isoceles, [3]Scalene
	      int Triangle_type;
	      if(i==j && j==k)
		Triangle_type=1;
	      else if(i==j && j!=k)
		Triangle_type=2;
	      else if(i!=j && j==k)
		Triangle_type=2;
	      else if(i!=j && j!=k && i!=k)
		Triangle_type=3;
	      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              real_prec ntri=0;
              real_prec c_bispectrum=0;
              real_prec c_shot_noise=0;
	      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	      switch(Triangle_type)
		{	  
		  
		case(1): // Equilateral
		  for(int ik=0;ik<Nt_max;++ik)
		    ntri+=(this->iFT_output_triangles_full[i][ik]*this->iFT_output_triangles_full[i][ik])*this->iFT_output_triangles_full[i][ik];
		  if(ntri>0)
		    for(int ik=0;ik<Nt_max;++ik)
		      c_bispectrum+=(this->iFT_output_delta_full[i][ik]*this->iFT_output_delta_full[i][ik])*this->iFT_output_delta_full[i][ik];
                  c_bispectrum/=(this->normal_b*(real_prec)ntri);
		  
		  for(int ik=0;ik<Nt_max;++ik)
		    c_shot_noise+=(this->iFT_shot_noise_p1_cyc_sum_full[i][ik]*this->iFT_output_triangles_full[i][ik])*this->iFT_output_triangles_full[i][ik];
                  c_shot_noise*=3.0/((real_prec)ntri);
		  
		  break;

		case(2): // Isoceles
		  if(j==k)
		    {
		      for(int ik=0;ik<Nt_max;++ik)
			ntri+=this->iFT_output_triangles_full[i][ik]*this->iFT_output_triangles_full[j][ik]*this->iFT_output_triangles_full[j][ik];
		      if(ntri>0)
			for(int ik=0;ik<Nt_max;++ik)
			  c_bispectrum+=this->iFT_output_delta_full[i][ik]*this->iFT_output_delta_full[j][ik]*this->iFT_output_delta_full[j][ik];
		      
		      
		      for(int ik=0;ik<Nt_max;++ik)
			c_shot_noise+=(this->iFT_shot_noise_p1_cyc_sum_full[i][ik]*this->iFT_output_triangles_full[j][ik])*this->iFT_output_triangles_full[k][ik];
		      
		      for(int ik=0;ik<Nt_max;++ik)
			c_shot_noise+=2.0*(this->iFT_shot_noise_p1_cyc_sum_full[j][ik]*this->iFT_output_triangles_full[i][ik])*this->iFT_output_triangles_full[k][ik];
		      
		      
		    }
		  else if(k!=j)
		    { //i==j
		      for(int ik=0;ik<Nt_max;++ik)
			ntri+=pow(this->iFT_output_triangles_full[i][ik],2)*this->iFT_output_triangles_full[k][ik];
		      if(ntri>0)
			for(int ik=0;ik<Nt_max;++ik)
			  c_bispectrum+=pow(this->iFT_output_delta_full[i][ik],2)*this->iFT_output_delta_full[k][ik];
		      
		      for(int ik=0;ik<Nt_max;++ik)
			c_shot_noise+=2.0*(this->iFT_shot_noise_p1_cyc_sum_full[i][ik]*this->iFT_output_triangles_full[i][ik])*this->iFT_output_triangles_full[k][ik];
		      
		      for(int ik=0;ik<Nt_max;++ik)
			c_shot_noise+=(this->iFT_shot_noise_p1_cyc_sum_full[k][ik]*this->iFT_output_triangles_full[i][ik])*this->iFT_output_triangles_full[i][ik];
		      
		      
		    }
		  c_bispectrum/=(this->normal_b*static_cast<real_prec>(ntri));	
		  
                  c_shot_noise/=(real_prec)ntri;
		  
		  ntri*=3; // Just to compare using the right number of modes in the full Fourier space!
		  
		  
		  break;
		  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
		case(3): // Scaleno
		  for(int ik=0;ik<Nt_max;++ik)
		    ntri+=(this->iFT_output_triangles_full[i][ik]*this->iFT_output_triangles_full[j][ik])*this->iFT_output_triangles_full[k][ik];
		  if(ntri>0)
		    for(int ik=0;ik<Nt_max;++ik)
		      c_bispectrum+=(this->iFT_output_delta_full[i][ik]*this->iFT_output_delta_full[j][ik])*this->iFT_output_delta_full[k][ik];
                  c_bispectrum/=(this->normal_b*(real_prec)ntri);
		  
		  
		  for(int ik=0;ik<Nt_max;++ik)
		    c_shot_noise+=(this->iFT_shot_noise_p1_cyc_sum_full[i][ik]*this->iFT_output_triangles_full[j][ik])*this->iFT_output_triangles_full[k][ik];
		  
		  for(int ik=0;ik<Nt_max;++ik)
		    c_shot_noise+=(this->iFT_shot_noise_p1_cyc_sum_full[j][ik]*this->iFT_output_triangles_full[i][ik])*this->iFT_output_triangles_full[k][ik];
		  
		  for(int ik=0;ik<Nt_max;++ik)
		    c_shot_noise+=(this->iFT_shot_noise_p1_cyc_sum_full[k][ik]*this->iFT_output_triangles_full[i][ik])*this->iFT_output_triangles_full[j][ik];
		  
		  
                  c_shot_noise/=(real_prec)ntri;
		  
		  ntri*=6;// Just to compare usingthe right number of modes in the full Fourier space!
		  
		  break;
		  
		  
		}
	      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	      
	      //for(int ik=0;ik<Nt_n;++ik)cout<<ntri<<"\t"<<ik<<"\t"<<i<<"\t"<<this->iFT_output_triangles_full[i][ik]<<"\t"<<j<<" \t"<<this->iFT_output_triangles_full[j][ik]<<"\t"<<k<<"\t"<<this->iFT_output_triangles_full[k][ik]<<endl;
	      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	      // Normalization of the Bispectrum: divide by a factor (Ngrid**3)*(Number of triangles)*(Normalization of the estimator)
	      // Shot-noise correction
	      // 	c_bispectrum-=(SN_aux*((pk[i]+pk[j]+pk[k])*this->shot_noise_b1+this->shot_noise_b2));
	      c_bispectrum-= SN_aux*(c_shot_noise*this->shot_noise_b1+this->shot_noise_b2);
	      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	      ULONG new_ind=index_3d(i,j,k,this->Nshells_bk,this->Nshells_bk);
	      bispect[new_ind]=c_bispectrum;
	      mod[new_ind]=ntri;
	      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	      out<<i<<"\t"<<j<<"\t"<<k<<"\t"<<c_bispectrum<<"\t"<<c_shot_noise<<"\t"<<ntri<<endl;
	      //cout<<i<<"\t"<<j<<"\t"<<k<<"\t"<<c_bispectrum<<"\t"<<c_shot_noise<<"\t"<<ntri<<endl;
	    }
	}
    }




  
  //  oi.close();
  out.close(); 
  std::cout<<CYAN<<"Output file "<<fname<<RESET<<std::endl;
  
  std::cout<<RED;
  time_t end; 
  time (&end);
  real_prec diff = difftime(end,start);
  if (diff<60) std::cout <<"Lapse: "<<diff<<" secs"<<endl;
  else if (diff<3600) std::cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
  else std::cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
  std::cout<<RESET;
  
  
 }







// // // ****************************************************************************************************************************
// // // ****************************************************************************************************************************
// // // old version, that originally coded with Jennifer. Does the same that the new version, but slower.
// // // ****************************************************************************************************************************

// void FftwFunctions::loop_shells_bispectrum(s_parameters_box *s_box, vector<real_prec> &pk, vector<real_prec> &bispect, vector< int > &mod, string fname){


//   std::cout<<BLUE<<"Loop over shells:"<<RESET<<endl;
//   time_t start;
//   time (&start);
  
//   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   ofstream out; 
//   out.open(fname.c_str() , ios::out); 
//   out.precision(10); 
//   out.setf(ios::showpoint); 
//   out.setf(ios::scientific); 
//   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   ofstream oi; oi.open("test.dat");
//   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   real_prec SN_aux;
//   if(s_box->use_SN_correction)SN_aux=1.0;
//   else SN_aux=0;
//   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   int ID_b;
//   int Nshell_min = (floor)( (this->kmin_bk/this->DeltaK_Bis) +0.5);  //Jennifer's binning
//   //int Nshell_min = (floor)( this->kmin_bk/this->DeltaK_Bis); //Official binning
  
//   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   int Ngrid_max=Ngrids_bk[this->Nshells_bk-1];
//   long Nt_max=(Ngrid_max*Ngrid_max)*Ngrid_max;
  
//   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   for(int i=Nshell_min;i<this->Nshells_bk;++i){  
    
    
//     //set size of grid in Fourier space for this particular shell
//     // The size of the grid is a factor 3 (arbitraty, to be defined)
//     // times the k value of the center of the shell:
//     long Ngrid = Ngrid_max;//this->Ngrids_bk[i];
//     long Nt_n=(Ngrid*Ngrid)*Ngrid;
//     //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     int KminID = this->kkminID[i];  //ID of the vector with kkmin in this shell
//     int KmaxID = KminID+this->Bmodes[i]-1; //ID of the vector with maximum kk in this shell
//     //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     // Allocate vector for IDFT with the new ngrid size
//     this->n2[0]=Ngrid;
//     this->n2[1]=Ngrid;
//     this->n2[2]=Ngrid;

//     //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     // Allocate arrays for IFT
//     // Ouput of the Inverse Fourier transformi

//     vector<real_prec> iFT_output_delta1(Nt_n,0);
//     vector<real_prec> iFT_output_triangles1(Nt_n,0);
//     //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     construct_shells(Ngrid,KminID, KmaxID, iFT_output_delta1, iFT_output_triangles1);

//     //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     int nstep2 = floor(i/2);  //start at this position for the second loop
//     int jstart=(nstep2==0 ? i: nstep2);  
//     //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//     for(int j=jstart;j<=i;++j){    
      
//       //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//       vector<real_prec> iFT_output_delta2(Nt_n,0);
//       vector<real_prec> iFT_output_triangles2(Nt_n,0);
//       //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//       int KminID2 = this->kkminID[j];  //ID of the vector with kkmin in this shell
//       int KmaxID2 = KminID2+this->Bmodes[j]-1; //ID of the vector with maximum kk in this shell
//       if(j!=i){
// 	construct_shells(Ngrid,KminID2, KmaxID2, iFT_output_delta2, iFT_output_triangles2);
//       }
//       else{
// 	if(i==j){
// 	  iFT_output_triangles2=iFT_output_triangles1;
// 	  iFT_output_delta2=iFT_output_delta1;
// 	}
//       }

//     //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     

//       // *************************************************************
//       // real_prec all2=0;//accumulate(iFT_output_triangles2.begin(),iFT_output_triangles2.end(),0);
//       // for(int ik=0;ik<Nt_n;++ik)all2+=(real_prec)iFT_output_triangles2[ik];
//       //   cout<<i<<"  "<<j<<"   "<<all1<<" "<<all2<<endl;
//       // *************************************************************
      

//       // Prepare for loop over k:
//       int nstep3=i-j-1;
      
//       int kstart=(nstep3>j ? j: nstep3);
      
//       if(nstep3<=0)kstart=Nshell_min;
      
//       for(int k=kstart;k<=j;++k){



// 	// Classify triangle bin as type [1] equilateral, [2]=isoceles, [3]Scalene
// 	int Triangle_type;
// 	if(i==j && j==k)Triangle_type=1;
// 	else if(i==j && j!=k)Triangle_type=2;
// 	else if(i!=j && j==k)Triangle_type=2;
// 	else if(i!=j && j!=k && i!=k)Triangle_type=3;
	
// 	real_prec ntri=0;
// 	real_prec c_bispectrum=0;

// 	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 	switch(Triangle_type){
	  
// 	case(1): // Equilateral
// 	  // Get number of triangles
// 	  for(int ik=0;ik<Nt_n;++ik)ntri+=pow(iFT_output_triangles1[ik],3);

// 	  // Get deltas:
// 	  for(int ik=0;ik<iFT_output_delta1.size();++ik)c_bispectrum+=pow(iFT_output_delta1[ik],3);
// 	  c_bispectrum/=(this->normal_b*(real_prec)ntri);

// 	  ntri/=(real_prec)Nt_n; //Just to compare
// 	  break;
	  
// 	case(2): // Isoceles
	  
// 	  if(j==k){
// 	    for(int ik=0;ik<Nt_n;++ik)ntri+=(real_prec)iFT_output_triangles1[ik]*pow((real_prec)iFT_output_triangles2[ik],2);
// 	    if(ntri>0)for(int ik=0;ik<Nt_n;ik++)c_bispectrum+=(real_prec)iFT_output_delta1[ik]*pow((real_prec)iFT_output_delta2[ik],2);
// 	  }
	  
// 	  else if(j!=k){
// 	    vector<real_prec> iFT_output_delta3(Nt_n,0);
// 	    vector<real_prec> iFT_output_triangles3(Nt_n,0);
	    
// 	    int KminID3 = this->kkminID[k];  //ID of the vector with kkmin in this shell
// 	    int KmaxID3 = KminID3+this->Bmodes[k]-1; //ID of the vector with maximum kk in this shell
	    
// 	    construct_shells(Ngrid,KminID3, KmaxID3, iFT_output_delta3, iFT_output_triangles3);
	    
// 	    for(int ik=0;ik<Nt_n;++ik)ntri+=((real_prec)iFT_output_triangles1[ik]*(real_prec)iFT_output_triangles2[ik])*(real_prec)iFT_output_triangles3[ik];
// 	    if(ntri>0)for(int ik=0;ik<Nt_n;++ik)c_bispectrum+=((real_prec)iFT_output_delta1[ik]*(real_prec)iFT_output_delta2[ik])*(real_prec)iFT_output_delta3[ik];
	    
// 	    iFT_output_delta3.clear();
// 	    iFT_output_triangles3.clear();
// 	  }
	  
// 	  c_bispectrum/=(this->normal_b*(real_prec)ntri);

// 	  ntri*=3.0/(real_prec)Nt_n; // Just to compare using the right number of modes in the full Fourier space!
	  
// 	  break;
	  
	  
// 	case(3): // Scaleno
// 	  // Get number of triangles
// 	  vector<real_prec> iFT_output_delta3(Nt_n,0);  //DEFINE ABOVE AND ALLOCATE HERE
// 	  vector<real_prec> iFT_output_triangles3(Nt_n,0);
	  
// 	  int KminID3 = this->kkminID[k];  //ID of the vector with kkmin in this shell
// 	  int KmaxID3 = KminID3+this->Bmodes[k]-1; //ID of the vector with maximum kk in this shell
	  
// 	  // Get deltas
// 	  construct_shells(Ngrid,KminID3, KmaxID3, iFT_output_delta3,iFT_output_triangles3);
// 	  for(int ik=0;ik<Nt_n;++ik)ntri+=((real_prec)iFT_output_triangles1[ik]*(real_prec)iFT_output_triangles2[ik])*iFT_output_triangles3[ik];
// 	  if(ntri>0)for(int ik=0;ik<Nt_n;++ik)c_bispectrum+=((real_prec)iFT_output_delta1[ik]*(real_prec)iFT_output_delta2[ik])*(real_prec)iFT_output_delta3[ik];
	  
// 	  iFT_output_delta3.clear();
// 	  iFT_output_triangles3.clear();
	  
// 	  c_bispectrum/=(this->normal_b*(real_prec)ntri);
	  
// 	  ntri*=6.0/(real_prec)Nt_n;// Just to compare usingthe right number of modes in the full Fourier space!


// 	  break;
// 	}


// 	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 	// &&&&&&&&&&&&&&&&&&&&&&
// 	// Normalization of the Bispectrum: divide by a factor (Ngrid**3)*(Number of triangles)*(Normalization of the estimator)
// 	// Shot-noise correction
// 	c_bispectrum-=(SN_aux*((pk[i]+pk[j]+pk[k])*this->shot_noise_b1+this->shot_noise_b2));
	
	
// 	out<<i<<"\t"<<j<<"\t"<<k<<"\t"<<c_bispectrum<<"\t"<<ntri<<endl;

	
// 	// bispect[ID_b]=c_bispectrum;
// 	// mod[ID_b]=ntri;
// 	// ID_b++;

//       }
      
//       iFT_output_delta2.clear();
//       iFT_output_triangles2.clear();
//     }
    
//     iFT_output_delta1.clear();
//     iFT_output_triangles1.clear();
//   }

//   oi.close();
//   out.close(); 
//   cout<<CYAN;
//   cout<<"Output file "<<fname<<endl;

//   std::cout<<RED;
//   time_t end; time (&end);
//   real_prec diff = difftime(end,start);
//   if (diff<60) std::cout <<"Lapse: "<<diff<<" secs"<<endl;
//   else if (diff<3600) std::cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
//   else std::cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
//   std::cout<<RESET;

  
// }



// ****************************************************************************************************************************
// ****************************************************************************************************************************


void FftwFunctions::construct_shells(int ngrid, int kmnid, int kmxid, vector<real_prec>&iFT_output_delta, vector<real_prec>&iFT_output_triangles, vector<real_prec>&iFT_output_p1_cyc_sum){
  

  long Ntt=ngrid*ngrid*(ngrid/2+1);
  // Inputs of the ift
  // Array used to get the count of triangles
  complex_prec *data_in_ks=(complex_prec *)fftw_malloc(2*Ntt*sizeof(real_prec));
  for(int i=0;i<Ntt;++i)data_in_ks[i][0]=0;
  for(int i=0;i<Ntt;++i)data_in_ks[i][1]=0;
  
  // Inputs of the ift
  // Array used to get the values of delta
  complex_prec *data_in_dk=(complex_prec *)fftw_malloc(2*Ntt*sizeof(real_prec));
  for(int i=0;i<Ntt;++i)data_in_dk[i][0]=0;
  for(int i=0;i<Ntt;++i)data_in_dk[i][1]=0;



  // Inputs of the ift
  // Array used to get shot noise power
  complex_prec *data_pk_shot_noise=(complex_prec *)fftw_malloc(2*Ntt*sizeof(real_prec));
  for(int i=0;i<Ntt;++i)data_pk_shot_noise[i][0]=0;
  for(int i=0;i<Ntt;++i)data_pk_shot_noise[i][1]=0;



  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Create 3d shell mesh to identify all vector components
  // This is a loop over the modes contained in this shell, from the
  // vector with the lowest value of k*k to the vector with the highest value.
  //  for(int i=kmnid;i<=kmxid;++i)cellsym(this->ArrayID[i],ngrid,data_in_ks, data_in_dk);
  for(int i=kmnid;i<=kmxid;++i)
    cellsym(i,ngrid,data_in_ks, data_in_dk, data_pk_shot_noise);
  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Do Inverse Fourier transforms:
  do_fftw_c2r(this->Nft,data_in_ks, iFT_output_triangles);
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // // Do IFT of data_in_dk
  do_fftw_c2r(this->Nft,data_in_dk, iFT_output_delta);
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // // Do IFT of data_pk_sn
  do_fftw_c2r(this->Nft,data_pk_shot_noise, iFT_output_p1_cyc_sum);

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fftw_free(data_in_ks);
  fftw_free(data_in_dk);
  fftw_free(data_pk_shot_noise);
  
}



// ****************************************************************************************************************************


void FftwFunctions::power_spectrum_fkp_for_bispectrum(s_parameters_box *s_box,vector<real_prec> & powerk_g0)
{
  //////////////////////////////////////////////////////////  
  // Shell average in Fourier space. 
  // for the bispectrum
  //////////////////////////////////////////////////////////  

  vector<int>mod_g(powerk_g0.size(),0);
  for(int i=0;i<mod_g.size();++i)mod_g[i]=0;
  for(int i=0;i<mod_g.size();++i)powerk_g0[i]=0;

  real_prec rDeltaK_data;
  real_prec rkmin;
  
  int my_s_box_ave;
  real_prec k_bin_index=(s_box->k_bin_step == 1.? 0.5:0.0);
  
#define s_box_linear 1
#define s_box_log 2
#define s_box_other 0
  
  my_s_box_ave = s_box_other;
  if(s_box->ave=="linear")  my_s_box_ave = s_box_linear;
  else if(s_box->ave=="log")my_s_box_ave = s_box_log;
    
  if(my_s_box_ave == s_box_linear) {
    rDeltaK_data = 1.0 / this->DeltaK_data;    
  }
  else if(my_s_box_ave == s_box_log)
    rkmin = 1.0/this->kmin;
  
  real_prec *yy_MAS_array, *zz_MAS_array;
  yy_MAS_array = new real_prec[this->Nft/2];
  zz_MAS_array = new real_prec[this->Nft/2+1];
  
  
  int NTHREADS = omp_get_max_threads();
  real_prec M_PI_rNft = M_PI * this->rNft;
  
  time_t start;
  time (&start);
  
#pragma omp parallel num_threads(NTHREADS)
  {
    real_prec i_per_fact, j_per_fact;
    real_prec i_deltak_x_2, j_deltak_y_2, k_deltak_z_2;
    real_prec xx_MAS, yy_MAS, zz_MAS;
    real_prec rxx_MAS, ryy_MAS, rzz_MAS;

#pragma omp sections
    {
#pragma omp section
      {
        yy_MAS_array[0] = 1.0;
        for(int j=1;j<this->Nft/2;++j) {
          yy_MAS = j * M_PI_rNft;
          yy_MAS_array[j] = sin(yy_MAS)/yy_MAS;
        }
      }
#pragma omp section
      {
        zz_MAS_array[0] = 1.0;
        for(int k=1;k<=this->Nft/2;++k){  // loop over octant kz>0: this is all FFTW gives.
          zz_MAS = k * M_PI_rNft;
          zz_MAS_array[k] = (zz_MAS = sin(zz_MAS)/zz_MAS);
        }
      }
    }
    
#pragma omp barrier
    
#pragma omp for nowait  //reduction(+:powerk_g0, powerk_g2, powerk_g4, mod_g, power2d_spher, mod_spher, power2d_cart, mod_cart, powerk_r, mod_r)
    for(int i=0;i<this->Nft/2;++i){  // loop over octant kx>0. Half of what FFT gives
      if(i > 0) {
        xx_MAS = i * M_PI_rNft;
        xx_MAS = sin(xx_MAS)/xx_MAS;
      } else xx_MAS = 1.0;
      
      i_deltak_x_2 = i * i * this->deltak_x * this->deltak_x ;
      i_per_fact = i_deltak_x_2;
      if(my_s_box_ave == s_box_log)
        i_deltak_x_2 *= rkmin * rkmin;
      
      for(int j=0;j<this->Nft/2;++j){  // loop over octant ky>0. Half of what FFT gives
        yy_MAS = yy_MAS_array[j];
	
        j_deltak_y_2 = j * j * this->deltak_y * this->deltak_y;
        j_per_fact = i_per_fact + j_deltak_y_2;
        j_per_fact = sqrt(j_per_fact);
        if(my_s_box_ave == s_box_log)
          j_deltak_y_2 *= rkmin * rkmin;
	
        for(int k=0;k<=this->Nft/2;++k){  // loop over octant kz>0: this is all FFTW gives.
          zz_MAS = zz_MAS_array[k];
	  
          k_deltak_z_2 = k * k * this->deltak_z * this->deltak_z;
          if(my_s_box_ave == s_box_log)
            k_deltak_z_2 *= rkmin * rkmin;
	  
          real_prec kv;
          int kmod_g;
        
          // Compute k-shell index 
          if(my_s_box_ave == s_box_linear){
            kv=sqrt(i_deltak_x_2 + j_deltak_y_2 + k_deltak_z_2);
	    kmod_g=(int)floor((float)(kv*rDeltaK_data)+k_bin_index); 
	  }
          else{
            if(my_s_box_ave == s_box_log) {
              kv= i_deltak_x_2 + j_deltak_y_2 + k_deltak_z_2;
              if(kv!=0){
                kmod_g=(int)floor((float)( log10(sqrt(kv))/this->Deltal))+1;
	      }
              else{kmod_g=0;}
            }
          }
	  
	  
          // Compute correction for MAS:
          real_prec corr=(s_box->use_MAS_correction ? my_correction_MAS(correction_MAS_exp, xx_MAS, yy_MAS, zz_MAS): 1.0);
          real_prec icorr2 = 1.0/(corr*corr);
	  
          // Compute index in c-order
          int lp=ijk(i,j,k,this->Nft,this->Nft,this->Nft/2+1);
	  
	  
          // Compute thisgs related to 2d power spectrum and multipole 
          // decomposition for the first octant
          // Define the angle bewtween k and los:       
          // In the FKP we need to specify the LOS direction. Set to z by default.
	  
	  
          //Exclude the zero-frequency onLY when computing the P: for the W keep it:
          if(lp!=0){
            real_prec Pk= (this->data_out_g[lp][0] * this->data_out_g[lp][0] + this->data_out_g[lp][1] * this->data_out_g[lp][1]) * icorr2;

            // 1d Spherical average
            if(kmod_g<powerk_g0.size()){
#pragma omp atomic              
              powerk_g0[kmod_g]  += Pk;
#pragma omp atomic              
              mod_g[kmod_g]++ ;    
            }
	    
	    //  *****************************************
	    //  Now, add modes in other octants:
	    //  Add negative frequencies in y:
	    //  *****************************************
	    if(j>0  && k>0){    
              lp=ijk(i,this->Nft-j,k,this->Nft,this->Nft,this->Nft/2+1);   // Go to the ky<0 octant
              real_prec Pk= (this->data_out_g[lp][0]*this->data_out_g[lp][0] + this->data_out_g[lp][1]*this->data_out_g[lp][1]) * icorr2;
	      
              // 1d Spherical average
              if(kmod_g<powerk_g0.size()){
#pragma omp atomic                              
		powerk_g0[kmod_g]  += Pk;
#pragma omp atomic                              
		mod_g[kmod_g]++ ;    
	      }
	    }	      
	    
	    // ******************************************
	    // Add negative frequencies in x:
	    // ******************************************
	    if(i>0  && (j>0 || k>0)){
	      lp=ijk(this->Nft-i,j,k,this->Nft,this->Nft,this->Nft/2+1);
              real_prec Pk= (this->data_out_g[lp][0]*this->data_out_g[lp][0] + this->data_out_g[lp][1]*this->data_out_g[lp][1]) * icorr2;
	      if(kmod_g<powerk_g0.size()){
#pragma omp atomic                                
		powerk_g0[kmod_g]  += Pk;
#pragma omp atomic                              
		mod_g[kmod_g]++;
	      }    
	    }
	    
	    // ******************************************
	    //  Add negative frequencies in x and y:
	    // ******************************************
	    if(i>0  && j>0  && k>0){
	      lp=ijk(this->Nft-i,this->Nft-j,k,this->Nft,this->Nft,this->Nft/2+1);
              real_prec Pk= (this->data_out_g[lp][0]*this->data_out_g[lp][0] + this->data_out_g[lp][1]*this->data_out_g[lp][1]) * icorr2;
	      if(kmod_g<powerk_g0.size()){
#pragma omp atomic                                  
		powerk_g0[kmod_g] += Pk;
#pragma omp atomic                                  
		mod_g[kmod_g] ++ ;    
	      }
	    }
	  }
	}
      } 
    }
  }
  
  
  // Subtract shot noise and normalize only he 1d power spectrum.
  real_prec SN_aux=0;
  if(s_box->use_SN_correction)
    SN_aux=this->shot_noise;

  
  // For the monopole:
  for(int i=0;i<powerk_g0.size();++i)powerk_g0[i]=(mod_g[i]== 0 ? 0 : powerk_g0[i]/((real_prec)mod_g[i]*this->normal_power)-SN_aux);
  //  for(int i=0;i<powerk_g0.size();++i)cout<<powerk_g0[i]<<"  "<<mod_g[i]<<endl;  
}

// ****************************************************************************************************************************
// ****************************************************************************************************************************

void FftwFunctions::cellsym(int id, int ngrid, complex_prec *data_ks, complex_prec *data_dk,complex_prec *data_pk_shot_noise){


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  int sx=this->Arraykx[id];
  int sy=this->Arrayky[id];
  int sz=this->Arraykz[id];
  int vt=this->VecArray[id];
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // MAS correction for each Fourier mode
  real_prec corr=1.0/this->Array_corr[id];
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  int ngrid_two=ngrid/2+1;
  int n_ip, ip;
  
  switch(vt){
  case(0): //The zero mode
    data_ks[0][0]=1;    
    data_ks[0][1]=0;    
    data_dk[0][0]=this->data_out_g[0][0]*corr;
    data_dk[0][1]=this->data_out_g[0][1]*corr;

    data_pk_shot_noise[0][0]=pow(data_dk[0][0],2)+pow(data_dk[0][1],2);
    data_pk_shot_noise[0][1]=0;

    


    break;
    
  case(6)://Mapping vectors along the poles

    // Positive freqs;
    // 00z
    ip=ijk(0,0,sz,this->Nft,this->Nft,this->Nft/2+1); // Original index. z-axis >0
    n_ip=ijk(0,0,sz,ngrid,ngrid,ngrid_two); //New index. z-axis >0
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    ip=ijk(0,sz,0,this->Nft,this->Nft,this->Nft/2+1); //y-axis
    n_ip=ijk(0,sz,0,ngrid,ngrid,ngrid_two); //y-axis
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;



    ip=ijk(sz,0,0,this->Nft,this->Nft,this->Nft/2+1); //x-axis
    n_ip=ijk(sz,0,0,ngrid,ngrid,ngrid_two); //x-axis
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    
    // Negative freqs:
    // 0-z0
    ip  =ijk(0,this->Nft-sz,0,this->Nft,this->Nft,this->Nft/2+1); //y-axis
    n_ip=ijk(0,ngrid-sz,0,ngrid,ngrid,ngrid_two); //y-axis
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;



    // -z00
    ip=ijk(this->Nft-sz,0,0,this->Nft,this->Nft,this->Nft/2+1); //x-axis
    n_ip=ijk(ngrid-sz,0,0,ngrid,ngrid,ngrid_two); //x-axis
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;  
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;


    break;

  case(8): //x=y=z
    // 1: +++ 
    ip=ijk(sz,sz,sz,this->Nft,this->Nft,this->Nft/2+1); // Original index. z-axis >0
    n_ip=ijk(sz,sz,sz,ngrid,ngrid,ngrid_two); //New index. z-axis >0
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //cout<<"case 1:   "<<sz<<"  "<<sz<<"   "<<sz<<endl;

    //2: +-+
    ip=ijk(sz,this->Nft-sz,sz,this->Nft,this->Nft,this->Nft/2+1); // Original index. z-axis >0
    n_ip=ijk(sz,ngrid-sz,sz,ngrid,ngrid,ngrid_two); //New index. z-axis >0
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //cout<<"case 2:   "<<sz<<"  "<<this->Nft-sz<<"   "<<sz<<endl;

    //3: -++
    ip=ijk(this->Nft-sz,sz,sz,this->Nft,this->Nft,this->Nft/2+1); // Original index. z-axis >0
    n_ip=ijk(ngrid-sz,sz,sz,ngrid,ngrid,ngrid_two); //New index. z-axis >0
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //cout<<"case 3:   "<<this->Nft-sz<<"  "<<sz<<"   "<<sz<<endl;

    // 4: --+
    ip=ijk(this->Nft-sz,this->Nft-sz,sz,this->Nft,this->Nft,this->Nft/2+1); // Original index. z-axis >0
    n_ip=ijk(ngrid-sz,ngrid-sz,sz,ngrid,ngrid,ngrid_two); //New index. z-axis >0
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 4:   "<<this->Nft-sz<<"  "<<this->Nft-sz<<"   "<<sz<<endl;

    break;

  case(12):     // x=0, y=z

    // ++0
    ip=ijk(sz,sz,0,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sz,sz,0,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 1:   "<<sz<<"  "<<sz<<"   "<<0<<endl;

    // +0+
    ip=ijk(sz,0,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sz,0,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 2:   "<<sz<<"  "<<0<<"   "<<sz<<endl;

    // 0++
    ip=ijk(0,sz,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(0,sz,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 3:   "<<0<<"  "<<sz<<"   "<<sz<<endl;

    // 0-+
    ip=ijk(0,this->Nft-sz,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(0,ngrid-sz,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 4:   "<<0<<"  "<<this->Nft-sz<<"   "<<sz<<endl;


    // -0+
    ip=ijk(this->Nft-sz,0,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sz,0,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 5:   "<<this->Nft-sz<<"  "<<0<<"   "<<sz<<endl;



    // +-0
    ip=ijk(sz,this->Nft-sz,0,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sz,ngrid-sz,0,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 6:   "<<sz<<"  "<<this->Nft-sz<<"   "<<0<<endl;

    
    //-+0
    ip=ijk(this->Nft-sz,sz,0,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sz,sz,0,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 7:   "<<this->Nft-sz<<"  "<<sz<<"   "<<0<<endl;


    //--0
    ip=ijk(this->Nft-sz,this->Nft-sz,0,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sz,ngrid-sz,0,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 8:   "<<this->Nft-sz<<"  "<<this->Nft-sz<<"   "<<0<<endl;


    break;

  case(23): //kx=0, kx!=ky, ky!=kz

    //1: 0yz
    ip=ijk(0,sy,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(0,sy,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 1:   "<<0<<"  "<<sy<<"  "<<sz<<endl;



    //2: 0zy
    ip=ijk(0,sz,sy,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(0,sz,sy,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 2:   "<<0<<"  "<<sz<<"  "<<sy<<endl;


    //3: y0z
    ip=ijk(sy,0,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sy,0,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 3:   "<<sy<<"  "<<0<<"  "<<sz<<endl;


    //4: yz0
    ip=ijk(sy,sz,0,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sy,sz,0,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 4:   "<<sy<<"  "<<sz<<"  "<<0<<endl;


    //5: z0y
    ip=ijk(sz,0,sy,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sz,0,sy,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 5:   "<<sz<<"  "<<0<<"  "<<sy<<endl;

    //6: zy0
    ip=ijk(sz,sy,0,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sz,sy,0,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 6:   "<<sz<<"  "<<sy<<"  "<<0<<endl;


    //7: -y-z0
    ip=ijk(this->Nft-sy,this->Nft-sz,0,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sy,ngrid-sz,0,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 7:   "<<this->Nft-sy<<"  "<<this->Nft-sz<<"  "<<0<<endl;


    //8: -z-y 0
    ip=ijk(this->Nft-sz,this->Nft-sy,0,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sz,ngrid-sy,0,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 8:   "<<this->Nft-sz<<"  "<<this->Nft-sy<<"  "<<0<<endl;


    //9: 0-yz
    ip=ijk(0,this->Nft-sy,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(0,ngrid-sy,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 9:   "<<0<<"  "<<this->Nft-sy<<"  "<<sz<<endl;


    //10: z-y0
    ip=ijk(sz,this->Nft-sy,0,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sz,ngrid-sy,0,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 10:   "<<sz<<"  "<<this->Nft-sy<<"  "<<0<<endl;


    //11: -yz0
    ip=ijk(this->Nft-sy,sz,0,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sy,sz,0,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 11:   "<<this->Nft-sy<<"  "<<sz<<"  "<<0<<endl;

    
    //12: -y0z
    ip=ijk(this->Nft-sy,0,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sy,0,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 12:   "<<this->Nft-sy<<"  "<<0<<"  "<<sz<<endl;

    //13: -z0y
    ip=ijk(this->Nft-sz,0,sy,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sz,0,sy,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 13:   "<<this->Nft-sz<<"  "<<0<<"  "<<sy<<endl;


    //14: y-z0 
    ip=ijk(sy,this->Nft-sz,0,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sy,ngrid-sz,0,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 14:   "<<sy<<"  "<<this->Nft-sz<<"  "<<0<<endl;

    //15: 0-zy
    ip=ijk(0,this->Nft-sz,sy,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(0,ngrid-sz,sy,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 15:   "<<0<<"  "<<this->Nft-sz<<"  "<<sy<<endl;


    //16: -zy0
    ip=ijk(this->Nft-sz,sy,0,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sz,sy,0,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 16:   "<<this->Nft-sz<<"  "<<sy<<"  "<<0<<endl;
    break;


  case(24): //x=y!=z not zero

    //1:xxz 
    ip=ijk(sx,sx,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sx,sx,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 1:   "<<sx<<"  "<<sx<<"  "<<sz<<endl;



    //2: xzx
    ip=ijk(sx,sz,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sx,sz,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 2:   "<<sx<<"  "<<sz<<"  "<<sx<<endl;


    //3: zxx
    ip=ijk(sz,sx,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sz,sx,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 3:   "<<sz<<"  "<<sx<<"  "<<sx<<endl;


    //4: z-xx
    ip=ijk(sz,this->Nft-sx,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sz,ngrid-sx,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 4   "<<sz<<"  "<<this->Nft-sx<<"  "<<sx<<endl;


    //5: x-xz 
    ip=ijk(sx,this->Nft-sx,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sx,ngrid-sx,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 5   "<<sx<<"  "<<this->Nft-sx<<"  "<<sz<<endl;


    //6: -xzx 
    ip=ijk(this->Nft-sx,sz,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sx,sz,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 6   "<<this->Nft-sx<<"  "<<sz<<"  "<<sx<<endl;


    //7: -zxx 
    ip=ijk(this->Nft-sz,sx,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sz,sx,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 7   "<<this->Nft-sz<<"  "<<sx<<"  "<<sx<<endl;
    
    //8: x-zx
    ip=ijk(sx,this->Nft-sz,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sx,ngrid-sz,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 8   "<<sx<<"  "<<this->Nft-sz<<"  "<<sx<<endl;


    //9: -xxz
    ip=ijk(this->Nft-sx,sx,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sx,sx, sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 9   "<<this->Nft-sx<<"  "<<sx<<"  "<<sz<<endl;

    //10: -z-xx
    ip=ijk(this->Nft-sz,this->Nft-sx,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sz,ngrid-sx,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;   
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0; 
   //    cout<<"case 10  "<<this->Nft-sz<<"  "<<this->Nft-sx<<"  "<<sx<<endl;



    //11: -x-xz
    ip=ijk(this->Nft-sx,this->Nft-sx,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sx,ngrid-sx,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 11  "<<this->Nft-sx<<"  "<<this->Nft-sx<<"  "<<sz<<endl;

    //12: -x-zx
    ip=ijk(this->Nft-sx,this->Nft-sz,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sx,ngrid-sz,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 12  "<<this->Nft-sx<<"  "<<this->Nft-sz<<"  "<<sx<<endl;
    
    break;
    
  case(25): //x!y y=z

    
    //1 : x z z  
    ip=ijk(sx,sz,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sx,sz,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 1  "<<sx<<"  "<<sz<<"  "<<sz<<endl;


    //2 : z x z  
    ip=ijk(sz,sx,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sz,sx,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 2  "<<sz<<"  "<<sx<<"  "<<sz<<endl;

    //3 : z z x  
    ip=ijk(sz,sz,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sz,sz,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 3  "<<sz<<"  "<<sz<<"  "<<sx<<endl;

    //4 : z -x z  
    ip=ijk(sz,this->Nft-sx,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sz,ngrid-sx,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 4  "<<sz<<"  "<<this->Nft-sx<<"  "<<sz<<endl;


    //5 : -z x z  
    ip=ijk(this->Nft-sz,sx,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sz,sx,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 5  "<<this->Nft-sz<<"  "<<sx<<"  "<<sz<<endl;



    //6 : x -z z  
    ip=ijk(sx,this->Nft-sz,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sx,ngrid-sz,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 6  "<<sx<<"  "<<this->Nft-sz<<"  "<<sz<<endl;



    //7 : -z z x  
    ip=ijk(this->Nft-sz,sz,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sz,sz,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 7  "<<this->Nft-sz<<"  "<<this->Nft-sz<<"  "<<sx<<endl;


    //8 : -x z z  
    ip=ijk(this->Nft-sx,sz,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sx,sz,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 8  "<<this->Nft-sx<<"  "<<sz<<"  "<<sz<<endl;


    //9 : z -z x  
    ip=ijk(sz,this->Nft-sz,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sz,ngrid-sz,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 9  "<<sz<<"  "<<this->Nft-sz<<"  "<<sx<<endl;

    //10 : -z -x z  
    ip=ijk(this->Nft-sz,this->Nft-sx,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sz,ngrid-sx,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 10  "<<this->Nft-sz<<"  "<<this->Nft-sx<<"  "<<sz<<endl;


    //11 : -x -z z  
    ip=ijk(this->Nft-sx,this->Nft-sz,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sx,ngrid-sz,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 11  "<<this->Nft-sx<<"  "<<this->Nft-sz<<"  "<<sz<<endl;

    //12 : -z -zx  
    ip=ijk(this->Nft-sz,this->Nft-sz,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sz,ngrid-sz,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 12  "<<this->Nft-sz<<"  "<<<this->Nft-sz<<"  "<<sx<<endl;
    break;
    
    
  case(48): //x!=y!=z
    
    //1 :xyz 
    ip=ijk(sx,sy,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sx,sy,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 1  "<<sx<<"  "<<sy<<"  "<<sz<<endl;
    

    //2 :xzy 
    ip=ijk(sx,sz,sy,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sx,sz,sy,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 2  "<<sx<<"  "<<sz<<"  "<<sy<<endl;


    //3 :yxz
    ip=ijk(sy,sx,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sy,sx,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 3  "<<sy<<"  "<<sx<<"  "<<sz<<endl;


    //4 :yzx
    ip=ijk(sy,sz,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sy,sz,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 4  "<<sy<<"  "<<sz<<"  "<<sx<<endl;


    //5 :zxy
    ip=ijk(sz,sx,sy,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sz,sx,sy,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 5  "<<sz<<"  "<<sx<<"  "<<sy<<endl;


    //6 :zyx
    ip=ijk(sz,sy,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sz,sy,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 6  "<<sz<<"  "<<sy<<"  "<<sx<<endl;

    //7 :x-yz
    ip=ijk(sx,this->Nft-sy,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sx,ngrid-sy,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 7  "<<sx<<"  "<<this->Nft-sy<<"  "<<sz<<endl;


    //8 :-xyz
    ip=ijk(this->Nft-sx,sy,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sx,sy,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 8  "<<this->Nft-sx<<"  "<<sy<<"  "<<sz<<endl;


    //9 :-x-yz
    ip=ijk(this->Nft-sx,this->Nft-sy,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sx,ngrid-sy,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 9  "<<this->Nft-sx<<"  "<<this->Nft-sy<<"  "<<sz<<endl;


    //10 :x-zy
    ip=ijk(sx,this->Nft-sz,sy,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sx,ngrid-sz,sy,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 10  "<<sx<<"  "<<this->Nft-sz<<"  "<<sy<<endl;



    //11 :-xzy
    ip=ijk(this->Nft-sx,sz,sy,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sx,sz,sy,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 11  "<<this->Nft-sx<<"  "<<sz<<"  "<<sy<<endl;


    //12 :-x-zy
    ip=ijk(this->Nft-sx,this->Nft-sz,sy,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sx,ngrid-sz,sy,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 12  "<<this->Nft-sx<<"  "<<this->Nft-sz<<"  "<<sy<<endl;


    //13 :y-xz
    ip=ijk(sy,this->Nft-sx,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sy,ngrid-sx,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 13  "<<sy<<"  "<<this->Nft-sx<<"  "<<sz<<endl;


    //14 :-yxz
    ip=ijk(this->Nft-sy,sx,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sy,sx,sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 14  "<<this->Nft-sy<<"  "<<sx<<"  "<<sz<<endl;

    //15 :-y-xz
    ip=ijk(this->Nft-sy,this->Nft-sx,sz,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sy,ngrid-sx, sz,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;     
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 15  "<<this->Nft-sy<<"  "<<this->Nft-sx<<"  "<<sz<<endl;   


    //16 :-yzx
    ip=ijk(this->Nft-sy,sz,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sy,sz,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 16  "<<this->Nft-sy<<"  "<<sz<<"  "<<sx<<endl;   


    //17 :y-zx
    ip=ijk(sy,this->Nft-sz,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sy,ngrid-sz,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 17  "<<sy<<"  "<<this->Nft-sz<<"  "<<sx<<endl;   



    //18 :-y-zx
    ip=ijk(this->Nft-sy,this->Nft-sz,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sy,ngrid-sz,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 18  "<<this->Nft-sy<<"  "<<this->Nft-sz<<"  "<<sx<<endl;   

    //19 :z-xy
    ip=ijk(sz,this->Nft-sx,sy,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sz,ngrid-sx,sy,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 19  "<<sz<<"  "<<this->Nft-sx<<"  "<<sy<<endl;   


    //20 :-zxy
    ip=ijk(this->Nft-sz,sx,sy,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sz,sx,sy,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 20  "<<this->Nft-sz<<"  "<<sx<<"  "<<sy<<endl;   


    //21 :-z-xy
    ip=ijk(this->Nft-sz,this->Nft-sx,sy,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sz,ngrid-sx,sy,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
//    cout<<"case 21  "<<this->Nft-sz<<"  "<<this->Nft-sx<<"  "<<sy<<endl;   

    //22 :-zyx
    ip=ijk(this->Nft-sz,sy,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sz,sy,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 22  "<<this->Nft-sz<<"  "<<sy<<"  "<<sx<<endl;   


    //23 :z-yx
    ip=ijk(sz,this->Nft-sy,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(sz,ngrid-sy,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 23  "<<sz<<"  "<<this->Nft-sy<<"  "<<sx<<endl;   


    //24 :-z-yx
    ip=ijk(this->Nft-sz,this->Nft-sy,sx,this->Nft,this->Nft,this->Nft/2+1); 
    n_ip=ijk(ngrid-sz,ngrid-sy,sx,ngrid,ngrid,ngrid_two); 
    data_ks[n_ip][0]=1;    
    data_ks[n_ip][1]=0;    
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;    
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;    
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //cout<<"case 24  "<<this->Nft-sz<<"  "<<this->Nft-sy<<"  "<<sx<<endl;   
      
    break;
  }
  
}
// ****************************************************************************************************************************



void FftwFunctions::get_bispectrum_fkp_fast(s_parameters_box *s_box, vector<real_prec> &pk, vector<real_prec> &bispect, vector<int> &mod, string file)
{
  
  // *******************************************************
  // This function generates the estimates of Bispectrum   *
  // based on Fortan code by Jennifer Pollack
  // *******************************************************
  std::cout<<BLUE<<"Bispectrum a la Jennifer:  "<<std::endl;
  
  cout<<RED<<"The binning in the method that computes the power spectrum"<<endl;
  cout<<"for the shot noise must have the same binning as that used in the bispectrum, (Jenn or official.)"<<endl;
  cout<<"Check that please"<<endl;
  
  // Define the shells:
  define_kshells(s_box);
  
  // Get the inverse FT in the defined shells
  get_ift_shells_bispectrum(s_box);

  // Loop pver k-shells and estimates of bispectrum
  loop_shells_bispectrum(s_box, pk, bispect, mod, file);



}

// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************




// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************

void FftwFunctions::remap(int nn1, int nn2, int nn3, int i, int j, int k, int *ooi, int *ooj, int *ook, real_prec *ofactor)
{
  
  //////////////////////////////////////////////////////////  
  // This auxiliary function returns the indices that are 
  // arguments of the function ijk(), used to retrieve (a C++ 
  // note, ijk is  a method of same class, I do not need to 
  // create another FFTW_FUNCTION object) the amplitudes of 
  // the DTF as given by the FFTW. This takes into account   
  // the fact that we are adding one more frequency, 
  // the -Nyquist to the output and that we use Hermitian 
  // symmetry to retrieve the negative section       
  // of the third component, i.e, 
  // delta(kx, ky, -kz)=delta(-kx, -ky, kz)*      
  // The ofactor is +1 for the real part,
  // -1 for the imaginary part when Hermitian symmetry
  // is explicitely used.
  //////////////////////////////////////////////////////////  

  if(k<=nn3/2){
    if(i<=nn1/2)*ooi=i; else *ooi=i-1;
    if(j<=nn2/2)*ooj=j; else *ooj=j-1;
    *ook=k;
    *ofactor=1.00;
  }
  else{
    if(k>=nn3/2+1){  //negative components of z
      if(i>0 && i<=nn1/2)*ooi=nn1-i;else if(i>=nn1/2+1)*ooi=nn1-i+1;
      if(j>0 && j<=nn2/2)*ooj=nn2-j;else if(j>=nn2/2+1)*ooj=nn2-j+1;
      if(i==0)*ooi=0;
      if(j==0)*ooj=0;
      *ook=nn3-(k-1);
      *ofactor=-1.00; 
    }

  }	  
}

// ********************************************************************************************
// ********************************************************************************************
// ********************************************************************************************


void FftwFunctions::get_total_weight(vector<bool> &uw, vector<real_prec>&ow, real_prec *t_weight)
{
  //////////////////////////////////////////////////////////  
  // Compute the total weight from the available statistical     
  // weights present in the catalogs. The total weight is 
  // taken as the product of those weights                             
  //////////////////////////////////////////////////////////  
  
  bool b1=uw.at(0);
  bool b2=uw.at(1);
  bool b3=uw.at(2);
  bool b4=uw.at(3);
  real_prec w1=ow.at(0);
  real_prec w2=ow.at(1);
  real_prec w3=ow.at(2);
  real_prec w4=ow.at(3);
  real_prec total_weight;
  if(!b1 && !b2 && !b3 && !b4)total_weight=1.000;
  // Use all
  if(b1 && b2 && b3 && b4)total_weight=w1*w2*w3*w4;
  // Use 1 3 4 
  if(b1 && !b2 && b3 && b4)total_weight=w1*w3*w4;
  // Use 1 2 4 
  if(b1 && b2 && !b3 && b4)total_weight=w1*w2*w4;
  // Use 1 2 3 
  if(b1 && b2 && b3 && !b4)total_weight=w1*w2*w3;
  // Use 2 3 4 
  if(!b1 && b2 && b3 && b4)total_weight=w2*w3*w4;
  // Use 3 4 
  if(!b1 && !b2 && b3 && b4)total_weight=w3*w4;
  // Use 2 4 
  if(!b1 && b2 && !b3 && b4)total_weight=w2*w4;
  // Use 2 3 
  if(!b1 && b2 && b3 && !b4)total_weight=w2*w3;
  // Use 1 4 
  if(b1 && !b2 && !b3 && b4)total_weight=w1*w4;
  // Use 1 3 
  if(b1 && !b2 && b3 && !b4)total_weight=w1*w3;
  // Use 1 2
  if(b1 && b2 && !b3 && !b4)total_weight=w1*w2;
  // Use 4
  if(!b1 && !b2 && !b3 && b4)total_weight=w4;
  // Use 3
  if(!b1 && !b2 && b3 && !b4)total_weight=w3;
  // Use 2
  if(!b1 && b2 && !b3 && !b4)total_weight=w2;
  // Use 1
  if(b1 && !b2 && !b3 && !b4)total_weight=w1;
  
  *t_weight=total_weight;

}











// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************

void FftwFunctions::get_mean_density_interpolated(real_prec redshift_min_sample, real_prec redshift_max_sample, real_prec ra, real_prec dec, real_prec redshift, vector< vector<real_prec> > &dndz_m, real_prec *nb)
{
  //*********************************************************************
  // Generate an interpolated value of the mean number density           
  // given the position in the sky of the object, when the depth        *
  // varies with the angular position. This interpolates the matrix     *
  // computed on the method dndz                                        *
  //*********************************************************************
  // COULD WE AVOID THIS?, WE DO NOT USE THE FULL MAP HERE
  long ipix=0;
#ifdef HEALPIX
  Healpix_Map<real_prec>map(log2(this->nside), RING);
  real_prec fac=M_PI/180.0;
  point.phi=ra*fac;  
  point.theta=0.5*M_PI-dec*fac;
  ipix=map.ang2pix(point);
#endif
  real_prec Delta_Z=(redshift_max_sample-redshift_min_sample)/dndz_m.size();
  int iz=(int)floor((float)((redshift-redshift_min_sample)/Delta_Z));
  *nb=dndz_m[iz][ipix];
}


// ****************************************************************************************************************************
// ****************************************************************************************************************************



void FftwFunctions::write_fftw_parameters()
{
  
  So.message_screen("Input values and parameters of the FFT:");
  So.message_screen("Dimension of the grid =",this->Nft,"³");
  So.message_screen("Lenght =",this->Lside," Mpc/h");
  if(this->binning=="linear"){
     So.message_screen("Using linearly-spaced spherical shells");
    So.message_screen("Bin size for P(k) = ",this->DeltaK_data," h/Mpc");
    So.message_screen("Bin size for W(k) = ",this->DeltaK_window," h/Mpc");
    So.message_screen("Nyquist Frequency = ",0.5*this->Nft*(this->deltak_0)," Mpc/h");
  }
  else{
    So.message_screen("Using log-spaced spherical shells ");
    So.message_screen("Bin size in log =", this->Deltal);
  } 
  std::cout<<RESET;

}





// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************

void FftwFunctions::write_fftw_parameters(void *p, string fname)
{
  
  struct s_parameters_box * s_cp= (struct s_parameters_box *)p;
  ofstream out; 
  out.open(fname.c_str() , ios::out); 
  out.precision(12); 
  out.setf(ios::showpoint); 
  out.width(12); 

  out<<"Selected options :"<<endl;
  out<<"Statistics: "<<this->statistics<<endl;
  out<<"MAS: "<<(s_cp->mas)<<endl;
  if(s_cp->use_MAS_correction) out<<"MAS correction: enabled"<<endl;
  else out<<"MAS correction: disabled"<<endl;
  if(s_cp->FKP_weight)out<<"Using FKP weights with Pest = "<<s_cp->Pest<<endl;
  else out<<"Using weights = 1"<<endl;
  if(s_cp->FKP_error_bars)out<<"Computing FKP error bars"<<endl;
  else out<<"Estimate without error bars."<<endl;
  if(s_cp->use_SN_correction)out<<"Shot-noise correction: enabled"<<endl;
  else out<<"Shot-noise correction: disabled.  "<<endl;

  out<<"******************************************************************************"<<endl;
  out<<"Sample information"<<endl;
  out<<"Number of objects = "<<this->n_gal<<endl;
  if(s_cp->use_random_catalog)out<<"Number of random objects = "<<this->n_ran<<endl;

  out<<"******************************************************************************"<<endl;
  out<<"Fourier space information :"<<endl;
  out<<"Lbox = "<<this->Lside<<"  Mpc /h"<<endl;
  out<<"Nft = "<<this->Nft<<"³"<<endl; 

  if(statistics=="Pk_y_ds")out<<"Number of shells to kmax = "<<this->sgrid<<std::endl;

  if(this->binning=="linear")
     {
      out<<"Using linearly-spaced spherical shells "<<endl;
      out<<"Bin size for P(k) = "<<this->DeltaK_data<<" h/Mpc"<<endl;
      out<<"Bin size for W(k) = "<<this->DeltaK_window<<" h/Mpc"<<endl;
      out<<"Nyquist Frequency = "<<0.5*this->Nft*(this->deltak_0)<<" h/Mpc"<<endl;
  }
  else
      {
    out<<"Using log-spaced spherical shells"<<endl;
    out<<"Bin size = "<<this->Deltal<<endl;
  } 

  out<<"****************************************************************************"<<endl;
  out<<"Parameters of the estimator: "<<endl;
  out<<"Normalization = "<<this->normal_power<<endl;
  out<<"Shot_noise (power) = "<<this->shot_noise<<"  (Mpc h^-1)^3"<<endl;
  out<<"Shot_noise (window) = "<<this->shot_noise_window<<endl;
  if(s_cp->use_random_catalog)out<<"Mean number density (from weights) = "<<(this->n_ran-this->w_r)/(s_cp->Pest*this->w_r)<<" (Mpc h^-1)^(-3)"<<endl; 
  out<<"Weighted number of objects = "<<this->w_g<<endl;
  if(s_cp->use_random_catalog)out<<"Weighted number of random objects = "<<this->w_r<<endl;
  out<<"alpha = "<<this->alpha<<endl;
  out<<"Sum nw² galaxies = "<<this->s_g<<endl;
  if(s_cp->use_random_catalog)out<<"Sum nw² random = "<<this->s_r<<endl;

  out<<"******************************************************************************"<<endl;
  time_t rawtime;
  time ( &rawtime );
  out<<"Date: "<<ctime (&rawtime)<<endl;


  out.close(); 
  cout<<YELLOW<<"Log-file written in "<<fname<<RESET<<endl;

}


// *******************************************************************
// *******************************************************************


