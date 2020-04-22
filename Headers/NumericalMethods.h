#ifndef _NMETHODS_
#define _NMETHODS_

# include "def.h"
# include "bstream.h"
# include <ctime>
# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include <cassert>
# include <sstream>

#ifdef _USE_OMP_
#include <omp.h>
#endif

# include <gsl/gsl_sf_gamma.h>
# include <gsl/gsl_bspline.h>
# include <gsl/gsl_multifit.h>
# include <gsl/gsl_statistics.h>
# include <gsl/gsl_histogram.h>
# include <gsl/gsl_histogram2d.h>
# include <gsl/gsl_sf_legendre.h>
# include <gsl/gsl_sf_expint.h>
# include <gsl/gsl_integration.h>
# include <gsl/gsl_errno.h>
# include <gsl/gsl_eigen.h>
# include <gsl/gsl_spline.h>
# include <gsl/gsl_sf_bessel.h>
# include <gsl/gsl_sf_gamma.h>
# include <gsl/gsl_sf_coupling.h>
# include <gsl/gsl_bspline.h>
# include <gsl/gsl_multifit.h>
# include <gsl/gsl_statistics.h>
# include <gsl/gsl_statistics_float.h>
# include <gsl/gsl_linalg.h>
# include <gsl/gsl_sf_erf.h>
# include <gsl/gsl_sf.h>
# include <gsl/gsl_sort_vector.h>
# include <gsl/gsl_sort_vector_float.h>
# include <gsl/gsl_rng.h>
# include <gsl/gsl_randist.h>
# include <gsl/gsl_roots.h>
# include <gsl/gsl_math.h>
# include <gsl/gsl_deriv.h>
# include <gsl/gsl_roots.h>
# include <gsl/gsl_poly.h>
# include <vector>
# include <numeric>
# include <algorithm>
# include <alm.h>
# include <alm_fitsio.h>
# include <healpix_map.h>
# include <healpix_map_fitsio.h>
# include <unistd.h>
# include <fftw3.h>
# include <gsl/gsl_dht.h>
#include "fftw_array.h"
#include "Type_structures_def.h"
using namespace std;



//##################################################################################
//##################################################################################
/**
 *@brief
 */
void randomize_vector(vector<ULONG>&data);

void  get_high_res_id_from_low_res_id(int Nft, int Nft_low,vector<ULONG>&id_info);



//void  get_low_res_id_from_high_res_id(int Nft, int Nft_low,vector<ULONG>&ID_low,vector<s_cell_info>&cell_info_low);
void  get_low_res_id_from_high_res_id(int Nft, int Nft_low,vector<s_cell_info>&cell_info_low);
//##################################################################################
//##################################################################################
/**
 *@brief 
 */
void get_scalar(string FNAME,vector<real_prec>&OUT,ULONG N1,ULONG N2,ULONG N3);
//##################################################################################
//##################################################################################
/**
 *@brief 
 */
void dump_scalar(const vector<real_prec>&A_rm,ULONG N1,ULONG N2,ULONG N3,int sample_number,string fname);
//##################################################################################
//##################################################################################
/**
 *@brief 
 */
//string to_string (real_prec);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
real_prec gsl_integration(gsl_real (*function)(gsl_real, void *) ,void *, gsl_real, gsl_real);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
real_prec gsl_integration3(int N, gsl_real (*function)(gsl_real, void *) ,void *,gsl_real ,gsl_real);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
real_prec gsl_integration2(gsl_real (*function)(gsl_real, void *) ,void *,vector<gsl_real>, vector<gsl_real>);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
real_prec gsl_integration2(gsl_real (*function)(gsl_real, void *) ,void *,vector<gsl_real>, vector<gsl_real>, gsl_real);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void gsl_get_GL_weights(gsl_real,gsl_real,gsl_integration_glfixed_table *, vector<gsl_real>&, vector<gsl_real>&);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
real_prec gsl_inter_pointers(real_prec *, real_prec *, int, real_prec );
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
real_prec gsl_inter(gsl_real *, gsl_real *, int, gsl_real);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
real_prec gsl_inter_new(const vector<gsl_real> &, const vector<gsl_real> &, gsl_real);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief 
 * @returns 
 */
real_prec gsl_inter_new2(vector<real_prec>&, vector<vector<real_prec>>&, int, real_prec);

#ifdef SINGLE_PREC
real_prec gsl_inter_new2(vector<gsl_real>&, vector<vector<gsl_real>>&, int, real_prec);
#endif
//##################################################################################
//##################################################################################
real_prec MAS_CIC_public(real_prec x);
void grid_assignment_cic(real_prec deltax, ULONG N1, real_prec Lside,real_prec x,real_prec y,real_prec z,real_prec weight, vector<real_prec>& field)
;
//##################################################################################
//##################################################################################
/**
 * @brief
 * @returns
 */
real_prec gsl_integration_sin_kernel(gsl_real (*function)(gsl_real, void *), void *p, gsl_real aux_var, gsl_real LowLimit,gsl_real UpLimit);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief
 * @returns
 */
real_prec gsl_integration_sin_kernel_lowlimit_inf(gsl_real (*function)(gsl_real, void *), void *,  gsl_real, gsl_real);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief
 * @returns
 */
real_prec gsl_integration_sin_kernel_lowlimit_inf2(gsl_real (*function)(gsl_real, void *), void *, gsl_real, gsl_real, gsl_integration_workspace *, gsl_integration_workspace *);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief
 * @returns
 */
void gsl_bspline(vector<gsl_real> &, vector<gsl_real>&, vector<gsl_real> &, vector<gsl_real> &);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief
 * @returns
 */
void indexx(vector<int>&arr, vector<long>&indx);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief
 * @returns
 */
void sort_index(int , int , int , int *, int *, int *);
//##################################################################################
//##################################################################################
//##################################################################################//////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
void matrix_inversion(vector< vector<real_prec> > &, vector< vector<real_prec> > &);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief
 * @returns
 */
void matrix_det(vector< vector<real_prec> >&, real_prec &);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief
 * @returns Container Bam::Kernel
 */
void get_eigen(vector<vector<real_prec> >&icov_masses, vector<real_prec>&masses);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void get_det_matrix(vector<vector<real_prec> >&matriz, real_prec &determinant);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void sort_vectors(vector<int>& v1,vector<int>&v2, vector<int>&v3, vector<int>& v4, vector<int>&v5, vector<int>& v6,vector<int>
                  & v7, vector<real_prec>& v8);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
real_prec real_sh(int l, int m, real_prec theta, real_prec phi);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
real_prec imag_sh(int l, int m, real_prec theta, real_prec phi);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
real_prec bessel(real_prec , int);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
real_prec dbessel(real_prec , int );
real_prec ddbessel(real_prec , int );
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
size_t m_countLines(istream& inStream);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
size_t m_getBuffer(istream& inStream);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
string dto_string (real_prec Number);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
real_prec get_mean(const vector<real_prec> &ini);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
real_prec get_var(const vector<real_prec> &ini);
real_prec get_var(real_prec mean, const vector<real_prec> &ini);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void log_smooth(vector<real_prec>&xin, vector<real_prec>&vin);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void gradfindif(ULONG N1,ULONG N2,ULONG N3,real_prec L1,real_prec L2,real_prec L3,vector<real_prec>in,vector<real_prec>&out,int dim);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function computes the pdf of the input field in. builds the kernel from the ratio
 * @brief If type is log, the pdf is measured in log(numinlog+in). The inputs min and max must be either in log(numinlog+"in") or in "in"
 * @returns pdf
 */
void calc_pdf(string type, ULONG N, ULONG Nk, real_prec maxk, real_prec mink, const vector<real_prec>&in, vector<real_prec>&pdfin);
//##################################################################################
//##################################################################################
//##################################################################################
void convert_ngp_to_cic(vector<real_prec>&in, vector<real_prec>&out);
void convert_cic_to_ngp(vector<real_prec>&in, vector<real_prec>&out);

//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void rankorder(int seed, vector<real_prec>dens, ULONG N, ULONG Nk, real_prec maxk, real_prec mink, vector<real_prec>&in, vector<real_prec>&pdfin, vector<real_prec>&pdfout);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief Computes Eigenvalues of the velocity shear
 */
void EigenValuesVweb(ULONG Nft, real_prec L1, vector<real_prec>&inx,vector<real_prec>&iny,vector<real_prec>&inz, vector<real_prec>&diver, vector<real_prec> &out1, vector<real_prec> &out2, vector<real_prec> &out3);

//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief Computes Eigenvalues ofthe tidal tensor of density filed. From CLASSLIN (FS Kitaura)
 */
void EigenValuesTweb(ULONG Nft, real_prec L1, const vector<real_prec> &delta, const vector<real_prec> &phi, vector<real_prec> &out1, vector<real_prec> &out2, vector<real_prec> &out3, vector<real_prec> &S2, vector<real_prec> &S3, vector<real_prec> &N2D);
void EigenValuesTweb_bias(ULONG Nft, real_prec L1, const vector<real_prec> &delta, const vector<real_prec> &phi, vector<real_prec> &S2, vector<real_prec> &S3, vector<real_prec> &N2D);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void PoissonSolver(real_prec Lbox, ULONG Nft,const vector<real_prec>&in, vector<real_prec>&out);


//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief
 * @returns Container Bam::Kernel
 */

ULONG index_12d(int i, int j, int k, int l, int m, int n, int o, int p, int q, int r, int s, int v, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np, int Nq, int Nr, int Ns, int Nv);

//##################################################################################
//##################################################################################
/**
 * @brief
 * @returns Container Bam::Kernel
 */
ULONG index_11d(int i, int j, int k, int l, int m, int n, int o, int p, int q, int r, int s, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np, int Nq, int Nr, int Ns);

//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief
 * @returns Container Bam::Kernel
 */

ULONG index_10d(int i, int j, int k, int l, int m, int n, int o, int p, int q, int r, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np, int Nq, int Nr);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief
 * @returns Container Bam::Kernel
 */
ULONG index_9d(int i, int j, int k, int l, int m, int n, int o, int p, int q, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np, int Nq);

//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief
 * @returns Container Bam::Kernel
 */
ULONG index_8d(int i, int j, int k, int l, int m, int n, int o, int p, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np);

//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief
 * @returns
 */
ULONG index_7d(int i, int j, int k, int l, int m, int n, int o, int Nj, int Nk, int Nl, int Nm, int Nn, int No);
//##################################################################################
//##################################################################################
//##################################################################################

/**
 * @brief
 * @returns Container Bam::Kernel
 */
ULONG index_6d(int i, int j, int k, int l, int m, int n, int Nj, int Nk, int Nl, int Nm, int Nn);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief 
 * @returns 
 */
ULONG index_5d(int i, int j, int k, int l, int m, int Nj, int Nk, int Nl, int Nm);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief 
 * @returns 
 */
ULONG index_4d(int i, int j, int k, int l, int Nj, int Nk, int Nl);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
ULONG index_3d(int i, int j, int k, int Nj, int Nk);
//##################################################################################
//##################################################################################
/**
 * @brief Given cartasian coordiantes and box info, this returns the ID of the cell where the object is located
 */
ULONG grid_ID(s_params_box_mas *params, const real_prec &x, const real_prec &y, const real_prec &z);
//ULONG grid_ID(s_params_box_mas *params, real_prec x, real_prec y, real_prec z);


//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
ULONG index_2d(int i, int j, int Nj);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief Real to Complex FFT
 * @returns 
 */
void do_fftw_r2c(int Nft, vector<real_prec>in, complex_prec *out);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief 
 * @returns 
 */
void do_fftw_c2r(int Nft, complex_prec *in, vector<real_prec>&out);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void do_fftw_3d(ULONG Nft, bool direction, complex_prec *in, complex_prec *out);


//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void get_cumulative(const vector<real_prec>&, vector<real_prec> &, unsigned long &);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
real_prec get_nobjects(const vector<real_prec>&);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
ULONG get_nobjects(const vector<ULONG>&);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void get_overdens(const vector<real_prec>&, vector<real_prec>&);

void get_overdens(const vector<real_prec>&, const vector<real_prec>&, vector<real_prec>&, bool);

//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void get_overdens(const vector<real_prec>&, real_prec mean, vector<real_prec>&);

//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief Given a mesh of N (perdim), this function returns the coordinates of a cell
 * @params N = Number of cells per dimention
 * @params index = C-ordered index of a cell (0,N*N*N-1)
 * @returns Min coordinates (lower-bottom values)
 */
void index2coords(ULONG N, ULONG index, real_prec *XG, real_prec *YG, real_prec *ZG );
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief
 */
void exchange_xy(int, const vector<real_prec>&,vector<real_prec>&);
//##################################################################################
//##################################################################################
/**
 * @brief
 */
void exchange_xz(int, const vector<real_prec>&,vector<real_prec>&);


//##################################################################################
//#################################################################################
real_prec tidal_anisotropy(real_prec lambdba1, real_prec lambda2, real_prec lambda3);

//##################################################################################
//##################################################################################
real_prec invariant_field_I(real_prec lambdba1, real_prec lambda2, real_prec lambda3);
//##################################################################################
//##################################################################################
real_prec invariant_field_II(real_prec lambdba1, real_prec lambda2, real_prec lambda3);
//##################################################################################
//##################################################################################
real_prec invariant_field_III(real_prec lambdba1, real_prec lambda2, real_prec lambda3);
//##################################################################################
//##################################################################################
real_prec ellipticity(real_prec lambdba1, real_prec lambda2, real_prec lambda3);
//##################################################################################
//##################################################################################
real_prec prolaticity(real_prec lambdba1, real_prec lambda2, real_prec lambda3);
//##################################################################################
//##################################################################################
/**
 * @brief Function exported from Patchy (FSK)
 */
#ifdef OMPPARRANGAR
real_prec GR_NUM(gsl_rng ** SEED, real_prec sigma ,int GR_METHOD, int jthread);
#else
real_prec GR_NUM(gsl_rng * SEED, real_prec sigma ,int GR_METHOD);
#endif

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
//##################################################################################
//##################################################################################
/**
 * @brief Function exported from Patchy (FSK)
 */

#ifdef OMPPARRANGAR
void create_GARFIELDR(ULONG N1,ULONG N2,ULONG N3,vector<real_prec> &delta, const vector<real_prec> &Power,gsl_rng ** seed);
#else
void create_GARFIELDR(ULONG N1,ULONG N2,ULONG N3,vector<real_prec> &delta, const vector<real_prec> &Power,gsl_rng * seed);
#endif
//##################################################################################
//##################################################################################
void create_GARFIELDR_from_WHITENOISE(string,ULONG N1,ULONG N2,ULONG N3, vector<real_prec>&in);
//##################################################################################
//##################################################################################
void create_GARFIELD_FIXED_AMP(ULONG N1,ULONG N2,ULONG N3,  vector<real_prec> &delta, const vector<real_prec> &Power, gsl_rng * seed);
//##################################################################################
//##################################################################################
void get_neighbour_cells(int Nft, int N_cells_bf, vector<s_nearest_cells>&);
//##################################################################################
//##################################################################################
real_prec k_squared(ULONG i,ULONG j,ULONG k,real_prec L1,real_prec L2,real_prec L3,ULONG N1,ULONG N2,ULONG N3);
//##################################################################################
//##################################################################################
void kernelcomp(real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3,ULONG N1, ULONG N2, ULONG N3, real_prec smol, int filtertype, string output_dir);
//##################################################################################
//##################################################################################
//##################################################################################
real_prec calc_kx(ULONG j,real_prec L2,ULONG N2);
//##################################################################################
//##################################################################################
real_prec calc_ky(ULONG j,real_prec L2,ULONG N2);
//##################################################################################
//##################################################################################
real_prec calc_kz(ULONG j,real_prec L2,ULONG N2);
//##################################################################################
//##################################################################################
void calc_twolptterm(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, const vector<real_prec>&phiv, vector<real_prec> &m2v);
//##################################################################################
//##################################################################################
void calc_LapPhiv(ULONG N1,ULONG N2,ULONG N3,real_prec L1,real_prec L2,real_prec L3, complex_prec *philv,vector<real_prec>&LapPhiv,int index1,int index2);
//##################################################################################
//##################################################################################
void calc_curlcomp(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, const vector<real_prec>&phiv, vector<real_prec> phiv2, vector<real_prec> &m2v, int comp);
//##################################################################################
//##################################################################################
void calc_mu2term(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, const vector<real_prec> &phiv, vector<real_prec> phiv2, vector<real_prec>&m2v);
//##################################################################################
//##################################################################################
void calc_Det(ULONG N1, ULONG N2, ULONG N3, real_prec L1, real_prec L2, real_prec L3, const vector<real_prec> &in, vector<real_prec> &out);
//##################################################################################
//##################################################################################
void convcomp(real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3,ULONG N1, ULONG N2, ULONG N3, const vector<real_prec> &in, vector<real_prec> &out, int filtertype,real_prec smol,string file_kernel);
//##################################################################################
//##################################################################################
real_prec linearvel3d(int index, real_prec kx, real_prec ky, real_prec kz, real_prec phi);
//##################################################################################
//##################################################################################
void convolvek(ULONG N1, vector<real_prec>&in, vector<real_prec> &kernel, vector<real_prec> &out);
//##################################################################################
//##################################################################################
int get_bin(real_prec x, real_prec xmin,int nbins, real_prec delta, bool);
//##################################################################################
//##################################################################################
int my_gsl_rng_uniform_(gsl_rng *r, int Nmax);
//##################################################################################
//##################################################################################
void sort_2vectors(vector<vector<int> >& v1,vector< vector<int> >&v2);
//##################################################################################
//##################################################################################
void sort_1dvectors(vector<int>& v1,vector<int> &v2);
//##################################################################################

template<class Type> Type get_max(const vector<Type> &in)
{

  if(in.size()==0)
    cerr<<"Error. Empty vector"<<endl;

  Type lka=-static_cast<Type>(LARGE_NUMBER);
  Type lkb;
  for(ULONG i=0;i<in.size();++i)
    {
      lkb=max(static_cast<Type>(in[i]), lka);
      lka=lkb;
    }
  return lkb;
}
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
template<class Type> Type get_min(const vector<Type>&in)
{
  if(in.size()==0)
    cerr<<"Error. Empty vector"<<endl;
  Type lka=static_cast<Type>(LARGE_NUMBER);
  Type lkb;
  for(int i=0;i<in.size();++i)
    {
      lkb=min(static_cast<Type>(in[i]), static_cast<Type>(lka));
      lka=lkb;
    }
  return lka;
}

template<class Type> Type get_min(const vector<Type>&in, bool exclude_zero)
{
  if(in.size()==0)
    cerr<<"Error. Empty vector"<<endl;
  Type lka=static_cast<Type>(LARGE_NUMBER);
  Type lkb;
  for(int i=0;i<in.size();++i)
    {
      if(in[i]!=0)
        {
          lkb=min(static_cast<Type>(in[i]), static_cast<Type>(lka));
      lka=lkb;
            }
          }
  return lka;
}


//##################################################################################
//##################################################################################
template<class Type> Type get_sum(const  vector<Type>&in)
{
  if(in.size()==0)
    cerr<<"Error. Empty vector"<<endl;

  Type lkb=0;
  //#pragma omp for reduction(+:lkb)
  for(int i=0;i<in.size();++i)
    lkb+=static_cast<real_prec>(in[i]);

  return lkb;
}

//##################################################################################
//##################################################################################


#endif

