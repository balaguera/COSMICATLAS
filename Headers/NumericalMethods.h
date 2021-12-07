#ifndef __NMETHODS__
#define __NMETHODS__




# include <boost/tuple/tuple.hpp>
# define GNUPLOT_ENABLE_PTY
# include "../gnuplot-iostream/gnuplot-iostream.h"
//http://stahlke.org/dan/gnuplot-iostream/

#ifdef _USE_PYTHON_
#include <Python.h>
//https://www.codeproject.com/Articles/820116/Embedding-Python-program-in-a-C-Cplusplus-code
#endif

# include "def.h"  // Do not define it before gnuplot or python
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

#ifdef _USE_HEALPIX_
# include <alm.h>
# include <alm_fitsio.h>
# include <healpix_map.h>
# include <healpix_map_fitsio.h>
# include <healpix_data_io.h>
# include <healpix_base.h>
# include <alm_powspec_tools.h>
# include <alm_healpix_tools.h>
# include <xcomplex.h>
# include <sort_utils.h>
# include <sharp_cxx.h>
#endif
# include <unistd.h>
# include <gsl/gsl_dht.h>
# include "fftw_array.h"
# include "Type_structures_def.h"
using namespace std;

//##################################################################################
//##################################################################################
/**
 *@brief
 */
void randomize_vector(vector<ULONG>&data);
//##################################################################################
/**
 *@brief
 */

void  get_high_res_id_from_low_res_id(int Nft, int Nft_low,vector<ULONG>&id_info);
//##################################################################################
/**
 *@brief
 */
ULONG factorial(int);
//##################################################################################
/**
 *@brief
 */
//void  get_low_res_id_from_high_res_id(int Nft, int Nft_low,vector<ULONG>&ID_low,vector<s_cell_info>&cell_info_low);
void  get_low_res_id_from_high_res_id(int Nft, int Nft_low,vector<s_cell_info_reduced>&cell_info_low);
//##################################################################################
/**
 *@brief 
 */
void get_scalar(string FNAME,vector<real_prec>&OUT,ULONG N1,ULONG N2,ULONG N3);
//##################################################################################
/**
 *@brief 
 */
void dump_scalar(const vector<real_prec>&A_rm,ULONG N1,ULONG N2,ULONG N3,int sample_number,string fname);
//##################################################################################
/**
 *@brief 
 */
//string to_string (real_prec);
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
real_prec gsl_integration(gsl_real (*function)(gsl_real, void *) ,void *, gsl_real, gsl_real);
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
real_prec gsl_integration3(int N, gsl_real (*function)(gsl_real, void *) ,void *,gsl_real ,gsl_real);
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
real_prec gsl_integration2(gsl_real (*function)(gsl_real, void *) ,void *,vector<gsl_real>, vector<gsl_real>);
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
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
real_prec gsl_inter(gsl_real *, gsl_real *, int, gsl_real);
//##################################################################################
/**
 * @brief 
 */
real_prec gsl_inter_new(const vector<gsl_real> &, const vector<gsl_real> &, gsl_real);
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
real_prec MAS_CIC_public(real_prec x);
void grid_assignment_cic(real_prec deltax, ULONG N1, real_prec Lside,real_prec x,real_prec y,real_prec z,real_prec weight, vector<real_prec>& field)
;
//##################################################################################
/**
 * @brief
 * @returns
 */
real_prec gsl_integration_sin_kernel(gsl_real (*function)(gsl_real, void *), void *p, gsl_real aux_var, gsl_real LowLimit,gsl_real UpLimit);
//##################################################################################
/**
 * @brief
 * @returns
 */
real_prec gsl_integration_sin_kernel_lowlimit_inf(gsl_real (*function)(gsl_real, void *), void *,  gsl_real, gsl_real);
//##################################################################################
/**
 * @brief
 * @returns
 */
real_prec gsl_integration_sin_kernel_lowlimit_inf2(gsl_real (*function)(gsl_real, void *), void *, gsl_real, gsl_real, gsl_integration_workspace *, gsl_integration_workspace *);
//##################################################################################
/**
 * @brief
 * @returns
 */
void gsl_bspline(vector<gsl_real> &, vector<gsl_real>&, vector<gsl_real> &, vector<gsl_real> &);
//##################################################################################
/**
 * @brief
 * @returns
 */
void indexx(vector<int>&arr, vector<int>&indx);
void indexx_ulong(vector<ULONG>&arr, vector<ULONG>&indx);
//##################################################################################
/**
 * @brief
 * @returns
 */
void sort_index(int , int , int , int *, int *, int *);
//##################################################################################
//##################################################################################
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
void sort_vectors(vector<ULONG>& v1,vector<ULONG>&v2, vector<ULONG>&v3, vector<ULONG>& v4, vector<ULONG>&v5, vector<ULONG>& v6,vector<ULONG>
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
void lin_smooth(vector<real_prec>&xin, vector<real_prec>&vin, int n);
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
void rankorder(int seed, vector<real_prec>dens, ULONG Nk, real_prec maxk, real_prec mink, vector<real_prec>&in, vector<real_prec>&pdfin, vector<real_prec>&pdfout);
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
//##################################################################################
//##################################################################################
void EigenValuesTweb(ULONG Nft, real_prec L1, const vector<real_prec> &delta, const vector<real_prec> &phi, vector<real_prec> &out1, vector<real_prec> &out2, vector<real_prec> &out3);
  /**
   * @brief  
   */
void EigenValuesTweb_bias(ULONG Nft, real_prec L1, const vector<real_prec> &delta, const vector<real_prec> &phi, vector<real_prec> &S2, vector<real_prec> &S3, vector<real_prec> &N2D);
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void PoissonSolver(real_prec Lbox, ULONG Nft, vector<real_prec>&in, vector<real_prec>&out);
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

/**
 * @brief
 * @returns Container Bam::Kernel
 */
ULONG index_6d(int i, int j, int k, int l, int m, int n, int Nj, int Nk, int Nl, int Nm, int Nn);
//##################################################################################
//##################################################################################
/**
 * @brief 
 * @returns ULONG index
 */
ULONG index_5d(int i, int j, int k, int l, int m, int Nj, int Nk, int Nl, int Nm);
//##################################################################################
//##################################################################################
/**
 * @brief 
 * @returns ULONG index
 */
ULONG index_4d(int i, int j, int k, int l, int Nj, int Nk, int Nl);
//##################################################################################
//##################################################################################
/**
 * @brief Computes the 3D index from the label coordinates i, j, k
 * @returns ULONG index
 */
ULONG index_3d(ULONG i, ULONG j, ULONG k, int Nj, int Nk);
//##################################################################################
//##################################################################################
/**
 * @brief Computes the 2D index from the label coordinates i, j
 * @returns ULONG index
 */
ULONG index_2d(ULONG i, ULONG j, ULONG Nj);

//##################################################################################
/**
 * @brief Given a mesh of N (per-dimension), this function returns the coordinates (indices) of the cell
 * @params index = C-ordered (or Fortran) index of a cell (0,N*N*N-1)
 * @params N = Number of cells per dimention
 * @returns indices XG, YG, ZG of the cell if i is in c/row-major order or
 * @returns indices ZG, YG, XG of the cell if i is in Fortran/column-major order
 */
void index2coords(ULONG index, ULONG N, ULONG &XG, ULONG &YG, ULONG &ZG );
//##################################################################################
//##################################################################################

/**
 * @brief Given cartasian coordiantes and box info, this returns the ID of the cell where the object is located
 */
ULONG grid_ID(s_params_box_mas *params, const real_prec &x, const real_prec &y, const real_prec &z);
//ULONG grid_ID(s_params_box_mas *params, real_prec x, real_prec y, real_prec z);
//##################################################################################
//##################################################################################
/**
 * @brief Real to Complex FFT
 * @returns 
 */
void do_fftw_r2c(int Nft, vector<real_prec>&in, complex_prec *out);
//##################################################################################
//##################################################################################
/**
 * @brief Complex to real FFTW
 * @returns 
 */
void do_fftw_c2r(int Nft, complex_prec *in, vector<real_prec>&out);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief  3D FFTW
 */
void do_fftw_3d(ULONG Nft, bool direction, complex_prec *in, complex_prec *out);



//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function applies a low-pass filter (top-hat) to a high resolution field on a mesh Nft_HR**3
 * @returns container with density field with resolution Nft_LRr**3
 */
void downsampling(ULONG Nft_HR, ULONG Nft_LR, int imas, vector<real_prec>&HR_field, vector<real_prec>&LR_field, real_prec Lbox);



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

/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
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
/**
 * @brief
 * @details
 * @param Eigenvalues of the tidal field
 */
real_prec tidal_anisotropy(real_prec lambda1, real_prec lambda2, real_prec lambda3);

//##################################################################################
//##################################################################################
/**
 * @brief Tidal field invariant
 * @details Get Invariant of the tidal field tensir
 * @param Eigenvalues of the tidal field lambda1, lambda2, lambda3 = \f$ \lambda_{1},\lambda_{2},\lambda_{3}\f$
*  @returns \f$ I_{1}= \lambda_{1}+\lambda_{2}+\lambda_{3} \f$
*/
real_prec invariant_field_I(real_prec lambda1, real_prec lambda2, real_prec lambda3);
//##################################################################################
//##################################################################################
/**
 * @brief  Tidal field Invariant
 * @param Eigenvalue of the tidal field lambda1  lambda2, lambda3 = \f$ \lambda_{1},\lambda_{2},\lambda_{3}\f$
 * @param Eigenvalues of the tidal field lambda1, lambda2, lambda3 = \f$ \lambda_{1},\lambda_{2},\lambda_{3}\f$
 * @param Eigenvalues of the tidal field lambda1, lambda2, lambda3 = \f$ \lambda_{1},\lambda_{2},\lambda_{3}\f$
 * @returns  \f$ I_{2}= \lambda_{1}\lambda_{2}+\lambda_{2}\lambda_{3}+\lambda_{1}\lambda_{3} \f$
 */
real_prec invariant_field_II(real_prec lambda1, real_prec lambda2, real_prec lambda3);
//##################################################################################
//##################################################################################
/**
 * @brief  Tidal field invariant
 * @param Eigenvalues of the tidal field \f$ \lambda_{1},\lambda_{2},\lambda_{3}\f$
 * @returns \f$ I_{3}= \lambda_{1}\lambda_{2}\lambda_{3} \f$
 */
real_prec invariant_field_III(real_prec lambda1, real_prec lambda2, real_prec lambda3);

//##################################################################################
/**
 * @brief Tidal field invariant
 * @param Eigenvalues of the tidal field \f$ \lambda_{1},\lambda_{2},\lambda_{3}\f$
 * @returns \f$ I_{3}= \lambda_{1}^{2}+\lambda_{2}^{2}+\lambda_{3}^{2} \f$
 */
real_prec invariant_field_IV(real_prec lambda1, real_prec lambda2, real_prec lambda3);

//##################################################################################
//##################################################################################
/**
* @brief  Ellipticity
* @param Eigenvalues of the tidal field \f$ \lambda_{1},\lambda_{2},\lambda_{3}\f$
* @details \f$ I_{3}= \lambda_{1}^{2}+\lambda_{2}^{2}+\lambda_{3}^{2} \f$
*/
real_prec ellipticity(real_prec lambda1, real_prec lambda2, real_prec lambda3);
//##################################################################################
//##################################################################################
/**
* @brief  Invariant IV
* @param Eigenvalues of the tidal field \f$ \lambda_{1},\lambda_{2},\lambda_{3}\f$
* @details \f$ I_{3}= \lambda_{1}^{2}+\lambda_{2}^{2}+\lambda_{3}^{2} \f$
*/
real_prec prolatness(real_prec lambda1, real_prec lambda2, real_prec lambda3);
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
/**
 * @brief 
 */
void create_GARFIELDR_from_WHITENOISE(string,ULONG N1,ULONG N2,ULONG N3, vector<real_prec>&in);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
void create_GARFIELD_FIXED_AMP(ULONG N1,ULONG N2,ULONG N3,  vector<real_prec> &delta, const vector<real_prec> &Power, gsl_rng * seed);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
void get_neighbour_cells(int Nft, int N_cells_bf, vector<s_nearest_cells>&);
void get_neighbour_cells_cat_analyze(int Nft, int N_cells_bf, vector<s_nearest_cells>&);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
real_prec k_squared(ULONG i,ULONG j,ULONG k,real_prec L1,real_prec L2,real_prec L3,ULONG N1,ULONG N2,ULONG N3);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
void kernelcomp(real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3,ULONG N1, ULONG N2, ULONG N3, real_prec smol, int filtertype, string output_dir);
//##################################################################################
/**
 * @brief 
 */
void kernelcomp_for_boost(real_prec L, real_prec d,ULONG N, vector<real_prec>&kernel, string out_dir);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
real_prec calc_kx(ULONG j,real_prec L2,ULONG N2);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
real_prec calc_ky(ULONG j,real_prec L2,ULONG N2);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
real_prec calc_kz(ULONG j,real_prec L2,ULONG N2);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
void calc_twolptterm(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, vector<real_prec>&phiv, vector<real_prec> &m2v);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
void calc_LapPhiv(ULONG N1,ULONG N2,ULONG N3,real_prec L1,real_prec L2,real_prec L3, complex_prec *philv,vector<real_prec>&LapPhiv,int index1,int index2);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
void calc_curlcomp(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, const vector<real_prec>&phiv, vector<real_prec> phiv2, vector<real_prec> &m2v, int comp);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
void calc_mu2term(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, const vector<real_prec> &phiv, vector<real_prec> phiv2, vector<real_prec>&m2v);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
void calc_Det(ULONG N1, ULONG N2, ULONG N3, real_prec L1, real_prec L2, real_prec L3, const vector<real_prec> &in, vector<real_prec> &out);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
void convcomp(real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3,ULONG N1, ULONG N2, ULONG N3, vector<real_prec> &in, vector<real_prec> &out, int filtertype,real_prec smol,string file_kernel);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
real_prec linearvel3d(int index, real_prec kx, real_prec ky, real_prec kz, real_prec phi);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
void convolvek(ULONG N1, vector<real_prec>&in, vector<real_prec> &kernel, vector<real_prec> &out);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
void give_power(real_prec Lbox, ULONG N1, vector<real_prec>&in,  vector<real_prec> &out);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
int get_bin(real_prec x, real_prec xmin,int nbins, real_prec delta, bool);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
int my_gsl_rng_uniform_(gsl_rng *r, int Nmax);
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
void sort_2vectors(vector<vector<ULONG> >& v1,vector< vector<ULONG> >&v2);
//##################################################################################
//##################################################################################
/**
 * @brief sorts vectors v1 and also returns v2 sorted according to v1. For closest neightbopugh applications, this
 * returns the element vc=v2[0] where v2 is sorted
 */
void sort_1dvectors(vector<ULONG>& v1,vector<ULONG> &v2);

//##################################################################################
/**
 * @brief sorts vectors v1 and also returns v2 sorted according to v1. For closest neightbopugh applications.
 * @return sorted vector v1, v2 and the element vc=v2[0] where v2 is sorted. Returning this value saves one loopkin the funciton.
 */
void sort_1dvectors_v2(vector<ULONG>& v1, vector<ULONG> &v2, ULONG &v1c, ULONG &v2c); // 2nd version of sort_1dvectors
/**
 * @brief 
 */
void sort_1dvectors_v3(vector<ULONG>& v1, vector<ULONG> &v2, ULONG &v2c); // 2nd version of sort_1dvectors
/**
 * @brief 
 */
void sort_1dvectors_iv2(vector<int>& v1, vector<int> &v2, int &v1c, int &v2c); // 2nd version of sort_1dvectors
//##################################################################################
//##################################################################################
/**
 * @brief 
 */
void swap_amp_fourier(ULONG, vector<real_prec>& v1, vector<real_prec> &v2); // 2nd version of sort_1dvectors
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################

//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//#####################TEMPLATES ###################################################

template <class Type>void SWAP_L(Type &a, Type &b)
{
  Type dum=a;
  a=b;
  b=dum;

}

//##################################################################################

template <class Type>void indexx(vector<Type>&arr, vector<Type>&indx)
{
  const int M=7,NSTACK=50;

  ULONG n = indx.size();
  long ir;
  long i,j,k;
  ULONG indxt;
  long jstack=-1;
  long l=0;
  Type a;
  int istack[NSTACK];
  
  ir=n-1;
  for (j=0;j<n;j++) indx[j]=j;
  for (;;) {
    if (ir-l < M)
     {
       for (j=l+1;j<=ir;j++) 
       {
         indxt=indx[j];
	       a=arr[indxt];
	       for (i=j-1;i>=l;i--) 
          {
	          if (arr[indx[i]] <= a) break;
       	  indx[i+1]=indx[i];
          }
	      indx[i+1]=indxt;
       }
      if (jstack < 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
      } 
      else 
      {
       k=(l+ir) >> 1;
       SWAP_L<Type>(indx[k],indx[l+1]);
        if (arr[indx[l]] > arr[indx[ir]]) {
        	SWAP_L<Type>(indx[l],indx[ir]);
        }
      if (arr[indx[l+1]] > arr[indx[ir]]) {
	      SWAP_L<Type>(indx[l+1],indx[ir]);
      }
      if (arr[indx[l]] > arr[indx[l+1]]) {
       	SWAP_L<Type>(indx[l],indx[l+1]);
      }
      i=l+1;
      j=ir;
      indxt=indx[l+1];
      a=arr[indxt];
      for (;;) {
    	 do i++; while (arr[indx[i]] < a);
	     do j--; while (arr[indx[j]] > a);
	     if (j < i) break;
	       SWAP_L<Type>(indx[i],indx[j]);
      }
      indx[l+1]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack >= NSTACK) {
	       cerr << "NSTACK too small in indexx." << endl;
      	exit(1);
      }
      if (ir-i+1 >= j-l) {
	      istack[jstack]=ir;
	      istack[jstack-1]=i;
      	ir=j-1;
      } 
      else {
	     istack[jstack]=j-1;
       istack[jstack-1]=l;
    	l=i;
     }
    }
  }
}


//##################################################################################
/**
* @brief This template functions calculates the rank of the vector v1
* @detAIL The rank is the position of the lowest value in the array v1
*/
template <typename Type> void sort_1d_vectors(vector<Type>&v1, ULONG &rank)
{
   ULONG n=v1.size();
   vector<Type>iwksp(n,0);
   indexx<Type>(v1,iwksp); //av1 must be int
   rank=iwksp[0];
}

//##################################################################################
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
template<typename Type> Type get_max(const vector<Type> &in)
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
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
template<typename Type> Type get_min(const vector<Type>&in)
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
//##################################################################################
//##################################################################################
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */

template<typename Type> Type get_min(const vector<Type>&in, bool exclude_zero)
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
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */

template<typename Type> Type get_mean_occupation_number(real_prec max, const vector<Type>&in)
{
  if(in.size()==0)
    cerr<<"Error. Empty vector"<<endl;
  ULONG Ntot=0;


  int Nd=static_cast<int>(max);
  vector<int>dist(Nd+1,0);
  real_prec delta=1;
  for(int i=0;i<in.size();++i)
      dist[get_bin(in[i],0,Nd,delta,true)]++;

 int mean=0;
 Ntot=0;

 for(int i=0;i<=Nd;++i)
    {
     mean+=i*dist[i];
     Ntot+=dist[i];
  }
 mean/=static_cast<real_prec>(Ntot);

  return static_cast<int>(mean);
}

//##################################################################################
//##################################################################################
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
template<typename Type> Type get_sum(const  vector<Type>&in)
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

