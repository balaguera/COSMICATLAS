#ifndef _NMETHODS_
#define _NMETHODS_

# include "def.h"
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

#include "FileOutput.h"

# include <gsl/gsl_sf_gamma.h>
# include <gsl/gsl_bspline.h>
# include <gsl/gsl_multifit.h>
# include <gsl/gsl_statistics.h>
# include <gsl/gsl_statistics_double.h>
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
# include <gsl/gsl_linalg.h>
# include <gsl/gsl_sf_erf.h>
# include <gsl/gsl_sf.h>
# include <gsl/gsl_sort_vector.h>
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

using namespace std;
//##################################################################################
//##################################################################################
//##################################################################################
/**
 *@brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
string to_string (double);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
double gsl_integration(double (*function)(double, void *) ,void *,double ,double);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
double gsl_integration3(int N, double (*function)(double, void *) ,void *,double ,double);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
double gsl_integration2(double (*function)(double, void *) ,void *,vector<double>, vector<double>);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
double gsl_integration2(double (*function)(double, void *) ,void *,vector<double>, vector<double>, double);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void gsl_get_GL_weights(double,double,gsl_integration_glfixed_table *, vector<double>&, vector<double>&);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
double gsl_inter_pointers(double *, double *, int, double );
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
double gsl_inter(double *, double *, int, double);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
double gsl_inter_new(vector<double> &, vector<double> &, double);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
double gsl_inter_new2(vector<double> &, vector<vector<double> >&, int, double);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
double gsl_integration_sin_kernel(double (*function)(double, void *), void *p, double aux_var, double LowLimit,double UpLimit);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
double gsl_integration_sin_kernel_lowlimit_inf(double (*function)(double, void *), void *, double, double);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
double gsl_integration_sin_kernel_lowlimit_inf2(double (*function)(double, void *), void *, double, double,gsl_integration_workspace *, gsl_integration_workspace *);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void gsl_bspline(vector<double> &, vector<double>&, vector<double> &, vector<double> &);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void indexx(vector<int>&arr, vector<long>&indx);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void sort_index(int , int , int , int *, int *, int *);
//##################################################################################
//##################################################################################
//##################################################################################//////////////////////////////////////////////////////////
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void matrix_inversion(vector< vector<double> > &, vector< vector<double> > &);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void matrix_det(vector< vector<double> >&, double &);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void get_eigen(vector<vector<double> >&icov_masses, vector<double>&masses);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void get_det_matrix(vector<vector<double> >&matriz, double &determinant);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void sort_vectors(vector<int>& v1,vector<int>&v2, vector<int>&v3, vector<int>& v4, vector<int>&v5, vector<int>& v6,vector<int>
		  & v7, vector<double>& v8);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
double real_sh(int l, int m, double theta, double phi);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
double imag_sh(int l, int m, double theta, double phi);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
double bessel(double , int);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
double dbessel(double , int );
double ddbessel(double , int );
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
string dto_string (double Number);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
double get_mean(vector<double> const ini);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
double get_var(vector<double> const ini);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void smooth(vector<double>&xin, vector<double>&vin);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void gradfindif(ULONG N1,ULONG N2,ULONG N3,double L1,double L2,double L3,vector<double>in,vector<double>&out,int dim);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void calc_pdf(ULONG N, ULONG Nk, double maxk, double mink, vector<double>in, vector<double>&pdfin);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void rankorder(int seed, vector<double>dens, ULONG N, ULONG Nk, double maxk, double mink, vector<double>&in, vector<double>&pdfin, vector<double>&pdfout);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void EigenValuesShear(ULONG Nft, double L1, vector<double>inx,vector<double>iny,vector<double>inz, vector<double>&diver, vector<double> &out1, vector<double> &out2, vector<double> &out3);

//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief Computes Eigenvalues ofthe tidal tensor of density filed. From CLASSLIN (FS Kitaura)
 */
void EigenValuesTweb(ULONG Nft, double L1, vector<double>in, vector<double> &out1, vector<double> &out2, vector<double> &out3);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void PoissonSolver(double Lbox, ULONG Nft, vector<double>in, vector<double>&out);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
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
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void do_fftw_r2c(ULONG Nft, vector<double>in, fftw_complex *out);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void do_fftw_c2r(ULONG Nft, fftw_complex *in, vector<double>&out);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void get_cumulative(vector<double>, vector<double> &, unsigned long &);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
double get_nobjects(vector<double>);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
ULONG get_nobjects(vector<ULONG>);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void get_overdens(vector<double>&, vector<double>&);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void get_overdens(vector<double>&, double mean, vector<double>&);

//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief
 * @returns
 */
void index2coords(ULONG N, ULONG index, double *XG, double *YG, double *ZG );
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief
 */
void exchange_xy(int, vector<double>,vector<double>&);


/* //################################################################################## */
/* //################################################################################## */
/* //################################################################################## */
/* /\** */
/*  * @brief Function exported from Patchy (FSK) */
/*  *\/ */
/* #ifdef OMPPARRANGAR */
/* real_prec GR_NUM(gsl_rng ** SEED, real_prec sigma ,int GR_METHOD, int jthread); */
/* #else */
/* real_prec GR_NUM(gsl_rng * SEED, real_prec sigma ,int GR_METHOD); */
/* #endif */

/* //################################################################################## */
/* //################################################################################## */
/* //################################################################################## */
/* /\** */
/*  * @brief Function exported from Patchy (FSK) */
/*  *\/ */

/* #ifdef OMPPARRANGAR */
/* void create_GARFIELDR(ULONG N1,ULONG N2,ULONG N3,real_prec *delta,real_prec * Power,gsl_rng ** seed); */
/* #else */
/* void create_GARFIELDR(ULONG N1,ULONG N2,ULONG N3,real_prec *delta,real_prec * Power,gsl_rng * seed); */
/* #endif */


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
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################

template<class Type> Type get_max(vector<Type>const in)
{
  Type lka=static_cast<Type>(-LARGE_NUMBER);
  Type lkb;
  for(int i=0;i<in.size();++i)
    {
      lkb=max(static_cast<double>(in[i]), static_cast<double>(lka));
      lka=lkb;
    }
  return lka;
}
//##################################################################################
template<class Type> Type get_min(vector<Type>const in)
{
  if(in.size()==0)
    cerr<<"Error. Empty vector"<<endl;
  Type lka=static_cast<Type>(LARGE_NUMBER);
  Type lkb;
  for(int i=0;i<in.size();++i)
    {
      lkb=min(static_cast<double>(in[i]), static_cast<double>(lka));
      lka=lkb;
    }
  return lka;
}


//##################################################################################


#endif
