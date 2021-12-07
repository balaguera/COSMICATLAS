// *****************************************************************************************************
// *****************************************************************************************************
// *****************************************************************************************************
/**
 * @class <FftwFuntions>
 * @brief    Header file for the class FftwFunctions::
 * @file     FftwFunctions.h
 * @title    Manipulation of function in Fourier space
 * @author   Andres Balaguera Antolinez
 * @author   Federico Marulli & Jennifer Pollack
 * @author   Optimization and parallelization by Luca Tornatore
 * @version  1.1a
 * @date     2013-2020
 * @details  NOTES RELATED TO FFTW
 The outputs of the FFTW has dimensions N1*N2*(N3/2+1),
 where only half of frequencies in the third component (kz)
 are stored. The other half of the third component can be found using Hermitian Symmetry (see below).
 The original output of the FFTW has the ordering

     index:    o_i      =  0    1     2     3    ...   N/2     N/2+1       N/2+2     ...   N-1
     freq:     k_i      =  0   k_f   2k_f  3k_f  ...   k_N   -k_n+k_f   -k_n+2k_f   ...   -k_f
     freq:    q_i/k_f   =  0    1     2     3    ...   N/2    -N/2+1      -N/2+2     ...   -1        in units of k_f                                                                                       *
 where k_f is the fundamental mode, N is the number of grid cells per dimension in the DFT and

     k_n= (N/2)k_f

 is the Nyquist frequency. Hence, the slot for the frequency -k_n is not written in the output. We need to explicitely use the negative and positive values
 to account for all posible configurations properly,
 specialy when computing the Bispectrum. For the Power spectrum it is not necessary
 for we do not expect to exploit information at the scales of the Nyquist freq.

 We therfore expand the loops over the wavenumbers by introducing one more
 slot, such that the new ordering (i) and coordinates (q_i/k_f) read as

    index:    i       =   0   1     2      j  ...  N/2,  N/2+1    N/2+2     N/2+3    ...        j            N-2     N-1    N
    freq:     k_i     =   0   k_f  2k_f  jk_f ...  k_N,  -k_N   -k_n+k_f  -k_n+2k_f  ...   -k_n+(j-1)k_f   -3k_f    -2k_f  -k_f
    freq:     q_i/k_f =   0   1     2      j  ...  N/2   -N/2   -N/2+1     -N/2+2    ...     -N/2+(j-1)      -3      -2     -1    in units of k_f                                                          *

  Hence, the new coordinates of the modes are

     q_i/k_f  = i      for i<=N/2,

  and

     q_i/k_f  =  i-(N+1)    for i>=N/2+1

  Similarly, to map the new index to the original output, we have

     o_i  = i       for i<=N/2,

  and

     o_i = i-1    for i>=N/2+1,

 such that for i=N/2, o_i=i and for i= N/2+1, o_i is also o_i=i.

 For the third (kz) component of the DFT we only have N3/2+1 elements displayed originally as

     index     o_k       =  0   1   2   3   i ...  N3/2
     freq     q_k/f_f    =  0   1   2   3   i ...  N3/2

  We reorder the components as we did above for the x-y components. However, in this case, the amplitude
  at a coordinate -k_i is obtained by Hermitian symmetry

        delta(kx,ky,-kz)=delta(-kx,-ky,kz)*.

 That is, when we use the complex conjugate to obtain the negative z-plane, we need to reverse the sign of the kx and ky components.
 Then, the reordered index is

       i3        =  0   1  2   j ...  N3/2   N3/2+1   N3/2+2    N3/2+3  ...     j      ... N3-2  N3-1   N3
    q_i3/k_f     =  0   1  2   j ...  N3/2    N3/2    N3/2-1    N3/2-2  ...  N3-(j-1) ...    3    2     1

 that is,

   q_i3/k_f = i     for i<=N/2,

 and

   q_i3/k_f = N3-(i-1)   for i>=N/2+1.

 The index need in the function ijk() is fabs( q_i / k_f ).
*/

// *****************************************************************************************************
// *****************************************************************************************************
// *****************************************************************************************************
// *****************************************************************************************************

#ifndef __FFTW_FUNCTIONS__
#define __FFTW_FUNCTIONS__


# include <ctime>
# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include <cassert>
# include <vector>
# include <algorithm>
# include <alm.h>
# include <alm_fitsio.h>
# include <healpix_map.h>
# include <healpix_map_fitsio.h>
# include <healpix_data_io.h>
# include <healpix_base.h>
# include <healpix_data_io.h>
# include <alm_powspec_tools.h>
# include <alm_healpix_tools.h>
# include <omp.h>
# include <complex>
# include "fftw_array.h"
# include "NumericalMethods.h"
# include "CosmologicalFunctions.h"
# include "FileOutput.h"
# include "CoordinateSystem.h"
# include "ScreenOutput.h"
# include "Params.h"
using namespace std;



class FftwFunctions{

 private:
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */

  FileOutput File;

  //////////////////////////////////////////////////////////
  /**
   *  @brief Obejct of type GalaxyOperations
   */
  CoordinateSystem Go;
  //////////////////////////////////////////////////////////
  /**
   * @brief Obejct of type ScreenOutput
   */
  ScreenOutput So;
  //////////////////////////////////////////////////////////
  /**
   * @brief Obejct of type ScreenOutput
   */
  Cosmology CosmoF;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Minimum X-coordinate of the sample
   */
  real_prec Xmin;

  //////////////////////////////////////////////////////////
  /**
   *  @brief Maximim X-coordinate of the sample
   */
  real_prec Xmax;

  //////////////////////////////////////////////////////////
  /**
   *  @brief Minimum Y-coordinate of the sample
  */
  real_prec Ymin;

  //////////////////////////////////////////////////////////
  /**
   *  @brief Maximum Y-coordinate of the sample
  */
  real_prec Ymax;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Minimum Z-coordinate of the sample
   */
  real_prec Zmin;
  /////////////////////////////////////////////////////////
  /**
   *  @brief Maximum Z-coordinate of the sample
   */
  real_prec Zmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  real_prec Lbox_data;  
 //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the estimator of the FKP variance of pwoer spectrum
   */
  vector<real_prec> SN;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the estimator of the FKP variance of pwoer spectrum
   */
  vector<real_prec> Q;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of random objects
   */
  ULONG n_ran;
 //////////////////////////////////////////////////////////
  /**
   * @brief Weighted number of galaxies
   */
  real_prec w_g;
  //////////////////////////////////////////////////////////
  /**
   * @brief Weighted number of random objects
     */
  real_prec w_r;
//////////////////////////////////////////////////////////
  /**
   * @brief Parameter alpha, the ratio between the weighted number of galaxies and the weighted number of randoms
   */
  real_prec alpha;

  //////////////////////////////////////////////////////////
  /**
   * @brief Parameter used to compute shot nnoise in power spectrum
   */
  real_prec s_g;
  //////////////////////////////////////////////////////////
  /**
   * @brief Parameter used to compute shot nnoise in power spectrum
   */
  real_prec s_r;
  //////////////////////////////////////////////////////////
  /**
   * @brief Parameter used to compute shot noise in bispectrum
   */
  real_prec sr1;
  //////////////////////////////////////////////////////////
  /**
      @brief Parameter used to compute shot noise in bispectrum
  */
  real_prec sr2;
//////////////////////////////////////////////////////////
  /**
      @brief Parameter used to compute shot noise in bispectrum
  */
  real_prec sg1;
//////////////////////////////////////////////////////////
  /**
      @brief Parameter used to compute shot noise in bispectrum
  */
  real_prec sg2;
  //////////////////////////////////////////////////////////
  /**
   * @brief Used in the normalization in bispectrum
   */
  real_prec normal_p;
//////////////////////////////////////////////////////////
  /**
   *    @brief  Normalization of window function
   */
  real_prec normal_window;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Used in the normalization in bispectrum
   */
  real_prec normal_b;
//////////////////////////////////////////////////////////
  /**
   *   @brief  Normalization in bispectrum
   */
  real_prec normal_bispectrum;
//////////////////////////////////////////////////////////
  /**
   *  @brief  Normalization in power spectrum
   */
  real_prec normal;
//////////////////////////////////////////////////////////
  /**
   * @brief  Poisson Shot noise for window function
   */
  real_prec shot_noise_window;
//////////////////////////////////////////////////////////
  /**
   * @brief Poisson Shot noise for bispectrum
   */
  real_prec shot_noise_b1;
//////////////////////////////////////////////////////////
  /**
   * @brief Poisson Shot noise for bispectrum
   */
  real_prec shot_noise_b2;
//////////////////////////////////////////////////////////
  /** used to avoid if blocks in correction_MAS */
  real_prec correction_MAS_exp;

//////////////////////////////////////////////////////////
  /**
  * \brief pointer to different Mass assignmetn scheme.
  * \details Three different MAS
  * are available for interpolation of the density
  * field. Insted of having "if" blocks inside a massively called function, we preset
  * three different functions; this pointer shall point to the correct one at runtime.
  */
  real_prec (FftwFunctions::*MAS_ptr)(real_prec);

//////////////////////////////////////////////////////////
  /**
  * \brief NGP Mass assignment scheme.
  */
  real_prec MAS_NGP(real_prec);

//////////////////////////////////////////////////////////
  /**
  * \brief CIC Mass assignment scheme.
  */
  real_prec MAS_CIC(real_prec);

//////////////////////////////////////////////////////////
  /**
  * \brief TSC Mass assignment scheme.
  */
  real_prec MAS_TSC(real_prec);

//////////////////////////////////////////////////////////
    /**
    * @brief PSC Mass assignment scheme.
    */
  real_prec MAS_PCS(real_prec);


//////////////////////////////////////////////////////////
  /**
   * @brief  Number of grid cells in each direction for Bispectrum fast
   */
  int sgrid;

//////////////////////////////////////////////////////////
  /**
   * @brief  Total number of grid cells used in arrays defined for the Bispectrum fast
   */
  int new_sgrid;


  //////////////////////////////////////////////////////////
   /**
     * @brief Offset in X direction
     */
  real_prec Xoffset;
  //////////////////////////////////////////////////////////
    /**
     * @brief Offset in Y direction
     */
  real_prec Yoffset;

  //////////////////////////////////////////////////////////
    /**
     * @brief Offset in Z direction
     */
  real_prec Zoffset;

//////////////////////////////////////////////////////////
  /**
   * @brief  inverse of the number of grid cells in each direction for the DFT
   */
  real_prec rNft;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue
   */
  vector <real_prec> cell_x;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue
   */
  vector <real_prec> cell_y;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue
   */
  vector <real_prec> cell_z;

//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_g_xx;

//////////////////////////////////////////////////////////
  /**
   * @brief Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_g_yy;

//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_g_zz;

//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_g_xy;

//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_g_xz;

//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_g_yz;

//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_xxx;

//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_yyy;

//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_zzz;

//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_xxy;

//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_xxz;

//////////////////////////////////////////////////////////
    /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_yyx;

//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_yyz;

//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_zzx;

//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_zzy;

//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_xyy;

//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_xzz;

//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_yzz;

//////////////////////////////////////////////////////////
    /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_xyz;

//////////////////////////////////////////////////////////
    /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_yxz;

//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_zxy;

//////////////////////////////////////////////////////////
    /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_r;

//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_r_xx;

//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_r_yy;

//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_r_zz;

//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_r_xy;

//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_r_xz;

//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_r_yz;

//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_r_xxx;

//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_r_yyy;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_r_zzz;

  //////////////////////////////////////////////////////////

  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_r_xxy;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_r_xxz;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_r_yyx;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_r_yyz;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_r_zzx;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_r_zzy;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_r_xyy;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_r_xzz;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_r_yzz;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_r_xyz;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_r_yxz;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_r_zxy;
  //////////////////////////////////////////////////////////

  /**
   * @brief  Vectors used in the Bispectrum
   */
  vector<ULONG> Arraykx;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum. Contains the MAS correction for each mode
   */
  vector<real_prec> Array_corr;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   */
  vector<ULONG> Arrayky;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   */
  vector<ULONG> Arraykz;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   */
  vector<ULONG> Arraykk;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   */
  vector<ULONG> VecArray;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   * Number of modes per shell for the Bispectrum
   */

  vector<int> Bmodes;


  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   */

  vector<int> kkminID;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   */

  vector<real_prec> kbins_bk;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   */
  vector<int> Ngrids_bk;

  //////////////////////////////////////////////////////////
  /**
   * @brief inverse of deltax/y/z, used in grid_assignment
   */
  real_prec rdeltax;
  //////////////////////////////////////////////////////////
  /**
   * @brief inverse of deltax/y/z, used in grid_assignment
   */

  real_prec rdeltay;
  //////////////////////////////////////////////////////////
  /**
   * @brief inverse of deltax/y/z, used in grid_assignment
   */
  real_prec  rdeltaz;

  //////////////////////////////////////////////////////////
  /**
   * @brief Number of galaxies
   */
  ULONG n_gal;
  //////////////////////////////////////////////////////////
  /**
   * @brief Product of NFT^2
   */
  int Nft2;
  //////////////////////////////////////////////////////////
  /**
      @brief Parameter used to compute the normalization in power spectrum
  */
  real_prec normal_power;
  //////////////////////////////////////////////////////////
  /**
      @brief Parameter used to compute the normalization in power spectrum
  */
  real_prec normal_power_two;

    //////////////////////////////////////////////////////////
  /**
   * @brief  Vector used in Bispectrum
   */

  vector<ULONG> ArrayID;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vector used in Bispectrum
   */
  vector<vector<real_prec>> iFT_output_delta_full;


  //////////////////////////////////////////////////////////
  /**
   * @brief  Vector used in Bispectrum
   */
  vector<vector<real_prec>> iFT_output_triangles_full;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vector used in Bispectrum.
   */
  vector<vector<real_prec>> iFT_shot_noise_p1_cyc_sum_full;

  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  int imcut;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   * Assigns to each vector in K-space the shell in which it's found
   */
  vector<ULONG> Kbin;
  //////////////////////////////////////////////////////////
  /**
   * @brief Complex vector for the output of the FFTW for the fluctuation
   */
  complex_prec *data_out_g;
  //////////////////////////////////////////////////////////
  /**
   * @brief Complex vector for the output of the FFTW for the fluctuation
   */
  complex_prec *data_out_g_rss;

  //////////////////////////////////////////////////////////
  /**
   * @brief Complex vector for the output of the FFTW for the fluctuation, used for the crossed power
   */
  complex_prec *data_out_gp;


  //////////////////////////////////////////////////////////
  /**
      @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_xx;

  //////////////////////////////////////////////////////////
  /**
   *    @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_yy;

  //////////////////////////////////////////////////////////
 /**
      @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_zz;

  //////////////////////////////////////////////////////////
 /**
  *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_xy;

  //////////////////////////////////////////////////////////
 /**
  *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_xz;

  //////////////////////////////////////////////////////////
 /**
  *    @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_yz;

  //////////////////////////////////////////////////////////
 /**
  *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_xxx;

  //////////////////////////////////////////////////////////
 /**
  *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_yyy;

  //////////////////////////////////////////////////////////
 /**
  *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_zzz;

  //////////////////////////////////////////////////////////
 /**
  *  @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_xxy;

  //////////////////////////////////////////////////////////
  /**
   *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
   */
  complex_prec *data_out_g_xxz;

  //////////////////////////////////////////////////////////
  /**
   *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_yyx;

  //////////////////////////////////////////////////////////
 /**
  *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_yyz;

  //////////////////////////////////////////////////////////
 /**
  *    @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_zzx;

  //////////////////////////////////////////////////////////
 /**
  *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_zzy;

  //////////////////////////////////////////////////////////
 /**
  *    @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_xyy;

  //////////////////////////////////////////////////////////
 /**
  *    @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_xzz;

  //////////////////////////////////////////////////////////
 /**
  *    @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_yzz;

  //////////////////////////////////////////////////////////
 /**
  *    @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_xyz;

  //////////////////////////////////////////////////////////
 /**
  *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_yxz;

  //////////////////////////////////////////////////////////
 /**
  *    @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_zxy;

  //////////////////////////////////////////////////////////
  /**
   *   @brief Complex vector for the output of the FFTW
   */
  complex_prec *data_out_r;

  //////////////////////////////////////////////////////////
  /**
   *  @brief Complex vector for the output of the FFTW used in the FKP estimation of the variance
   */
  complex_prec *data_out_SN;

  //////////////////////////////////////////////////////////
    /**
   *  @brief Complex vector for the output of the FFTW used in the FKP estimation of the variance
   */
  complex_prec *data_out_Q;

  //////////////////////////////////////////////////////////
  /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector< complex<real_prec> > data_g_out_y0;

  //////////////////////////////////////////////////////////
  /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector< complex<real_prec> > data_g_out_y2;

  //////////////////////////////////////////////////////////
  /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector< complex<real_prec> > data_g_out_y4;

  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector< complex<real_prec> > data_r_out_y0;

  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector< complex<real_prec> > data_r_out_y2;

  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector< complex<real_prec> > data_r_out_y4;

  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector<real_prec> SN_g_out_y2;

  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector<real_prec> SN_g_out_y4;

  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector<real_prec> SN_r_out_y2;

  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector<real_prec> SN_r_out_y4;

  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector<real_prec> data_g_y0;

  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector<real_prec> data_g_y2;

  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector<real_prec> data_g_y4;

  //////////////////////////////////////////////////////////
    /**
   * @brief Healpix object used when nbar is not tabulated
   */

#ifdef HEALPIX
  pointing point;
#endif
  //////////////////////////////////////////////////////////
  /**
  * @brief Mass assignment scheme.
  * @details Three different MAS
  * are available for interpolation of the density
  * field.
  * @params x Cartesian coordinate
  * @returns Value of the kernel function according to the MAS selected
  */
  real_prec MAS(real_prec x);
  //////////////////////////////////////////////////////////
  /**
   * @brief Estimation of the variance of the power spectrum
   * @details using  the FKP estimator with the exact expression, equation
   * 2.4.6  of FKP paper
   */
  void do_fkp_error_bars_exact(vector<real_prec> &, vector<real_prec> &, vector<real_prec> &,vector<real_prec> &);
  //////////////////////////////////////////////////////////
  /**
   * @brief Computes the variance in the power spectrum
   * @details as an approximation to the original FKP  estimation, introducing
   * the definition of the effective volume. This is computed
   * using the random catalogue.
   * @return Effective volume as a function of k:
   */
  void do_fkp_error_bars_veff(s_data_structure *, vector<real_prec> &,vector<real_prec> &);
  //////////////////////////////////////////////////////////
  /**
   * @brief Object of type Params
   */

  Params params;

  //////////////////////////////////////////////////////////

 public:
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  // THREE TYPES OF CONSTRUCTORS, FOR DIFFERET TYPES OF DEFINITIONS

  /**
   * @brief Constructor
  *  @param Inizialization of private variables
   */
  FftwFunctions(){}

// *****************************************************************************************************************************************************************
// *****************************************************************************************************************************************************************
// *****************************************************************************************************************************************************************
// *****************************************************************************************************************************************************************
// *****************************************************************************************************************************************************************
// *****************************************************************************************************************************************************************
  FftwFunctions(Params _params):params(_params), n_ran(0),w_g(0),w_r(0),alpha(0),s_g(0),s_r(0),sr1(0), sr2(0),sg1(0),sg2(0),normal_p(0),normal_window(0.),normal_b(0),normal(0),shot_noise(0),shot_noise2(0),shot_noise_b1(0),shot_noise_b2(0),correction_MAS_exp(0), sgrid(0),new_sgrid(0),rNft(128),Lside_data(0),DeltaK_Bis(0), Zmin(0),Zmax(0), Ymin(0),Ymax(0),Xmin(0),Xmax(0),nside(1), Nshells_bk(1), npixels(1),area_pixel(0),imcut(0),rdeltax(0),rdeltay(0),rdeltaz(0),n_gal(0),Xoffset(0), Yoffset(0), Zoffset(0)
  {
    this->rdeltax = static_cast<real_prec> (1.0/this->params._d_delta_x());
    this->rdeltay = this->rdeltax;
    this->rdeltaz = this->rdeltax;
    
    this->DeltaK_Bis    = this->params._d_DeltaK_data();// (this->params._kmax_bk()-this->params._kmin_bk())/((real_prec)this->Nshells_bk)
    this->Nshells_bk    = static_cast<int>(this->params._kmax_bk()/this->DeltaK_Bis); 

  ULONG nff;
  if(this->params._statistics()!="Pk_y_ds")
    nff=this->params._Nft();
  // Get the number of grid-cells for the direct sum DFT as a function of kmax
  else if(this->params._statistics()=="Pk_y_ds")
    {
      nff=(int)(this->params._kmax_y_ds()*this->params._Lbox()/M_PI); 
      if(nff%2!=0)nff++;
      this->params.set_Nft(nff);
      this->sgrid=nff;
    }
    
  // Get number of grid cells per dimension used for the bispectrum_fast given kmax
  int n_sgrid= static_cast<int>(this->params._kmax_bk()*this->params._Lbox()/M_PI);
  if(this->params._statistics()=="Bk_fkp_fast")
    {
      if(n_sgrid>this->params._Nft())
        {
         So.message_warning("Warning: kmax greater than Nyquist frequency");
         So.message_screen("Setting kmax to the Nyquist");
         n_sgrid=this->params._Nft();
         this->params.set_kmax_bk(this->params._Nft()*M_PI/this->params._Lbox() );
       }
     if(n_sgrid%2!=0)n_sgrid++;
     this->sgrid=n_sgrid;
     if(this->params._kmin_bk()<2.*M_PI/this->params._Lbox() )
      {
        So.message_warning("Warning: kmin smaller than fundamental mode");
        So.message_screen("Setting kmin as fundamental mode");
        this->params.set_kmin_bk(2.*M_PI/this->params._Lbox()) ;     
      }
    }
  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(this->params._statistics()=="Bk_fkp_fast")
    {
      // Get effective number of grid cells 
      // resulting from looping over half of the positive quadrant.
      // The result will be used to allocate memory for the arrays
      // used in the Bispectrum as performed by Jennifer.
      int new_sd=0;
      for(int i=0;i<this->sgrid;++i)
        for(int j=i;j<this->sgrid;++j)
          for(int k=j;k<this->sgrid;++k)
            if(i*i+j*j+k*k>0)
              new_sd++;    
      this->new_sgrid=new_sd;
    }
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Initialize other private variables for the Pk
  this->Nft2=(int)(nff*nff);
  this->rNft =1.0 / this->params._Nft();

  this->params.set_NGRID(nff*nff*nff);
  this->params.set_NGRID_h(nff*nff*(nff+2));



    time_t time_bam;
    time(&time_bam);
    this->So.initial_time=time_bam;
  }
// *****************************************************************************************************************************************************************
// *****************************************************************************************************************************************************************
// *****************************************************************************************************************************************************************
// *****************************************************************************************************************************************************************
// *****************************************************************************************************************************************************************
// *****************************************************************************************************************************************************************
  /**
   * @brief Destructor
   * @details Free fftw vector
   */
  ~FftwFunctions(){
    free_fftw_vectors();
  }

// *****************************************************************************************************************************************************************
// *****************************************************************************************************************************************************************
// *****************************************************************************************************************************************************************
// *****************************************************************************************************************************************************************
// *****************************************************************************************************************************************************************
// *****************************************************************************************************************************************************************

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue
   */
  vector <real_prec> data_g;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue
   */

  vector <real_prec> data_g_rss;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue
   */
  vector <real_prec> data_g_mw;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue used in case a cross power spectrum is to be measured
   */
  vector <real_prec> data_gp;



  //////////////////////////////////////////////////////////
  /**
   * @brief Lengh (in Mpc/h) of each side the box in configuration space
   computed from the catalog
  */
  real_prec Lside_data;

  //////////////////////////////////////////////////////////
  /**
   * @brief Width of the spherical shell (in h/Mpc) for the window function
   */
  real_prec DeltaK_window;

  //////////////////////////////////////////////////////////
    /**
     * @brief  Poisson Shot noise
     */
  real_prec shot_noise;
  //////////////////////////////////////////////////////////
   /**
     * @brief  Poisson Shot noise
     */
  real_prec shot_noise2;
  //////////////////////////////////////////////////////////
  /**
   * @brief Nside for Healpix
   */
  int nside;

  //////////////////////////////////////////////////////////
  /**
   * @brief Number of shells for the estimation of Bk as done by Jennifer
   */
  int Nshells_bk;

  //////////////////////////////////////////////////////////
  /**
   * @brief Number of pixels according to the Nside value
   */
  long npixels;

  //////////////////////////////////////////////////////////
  /**
   * @brief Area of one Healpix pixel
   */
  real_prec area_pixel;

  //////////////////////////////////////////////////////////
  /**
   * @brief Width of the spherical shell (in h/Mpc) used int he estimates of Bispectrum
   */
  real_prec DeltaK_Bis;


  ///////////////////////////////////////////////////////
  /**
   * @brief Transform the coordinates of the input catalogues to cartessian coordinates
   * @details returning the same input vector in which the three first columns correspond to
   * to the X,Y,Z coordinates and in the i_nbar (see input parameter file)
   * the mean number density is written.
   *
   * @param s_b: structure containing parametrs of the Fourier box
   * @param s_d structure containing information related to the catalogue
   * @param cat: 2d container with the catalog
   * @result Catalogue with position in cartesian coordinates and mean number density
   * tabulated in the corresponding column as stated in the parameter file
   */
  void cart_coordinates(s_parameters_box *s_b, s_data_structure *s_d);
  ///////////////////////////////////////////////////////
  /**
   * @brief Assign the Healpix resolution to the private variable of this class
   */

  void set_healpix_pars(int Nres);

  //////////////////////////////////////////////////////////
  /**
   *@brief Compute the total statistical weight of a given galaxy
   * @details Compute the total weight from the available statistical
   * weights present in the catalogs.
   * @param uw: bool vector indicating whether a weight is used or not.
   * @param ow: real_prec vector with the value of the weight
   * @param t_weight: total weight computed as the product of weights chosen to be used
   * @result Total statistical weight for a particular galaxy
   * @note The code accepts four (4) different weights (per galaxy)
   * plus the FKP weights. See input parameter file.
   */
  void get_total_weight(vector<bool>&uw, vector<real_prec> &ow, real_prec*t_weight);

  //////////////////////////////////////////////////////////
  /**
   * @brief Resizes and initializes the input and
   * output vectors for the FFTW
   */
  void resize_fftw_vectors();
  //////////////////////////////////////////////////////////
  /**
   * @brief Interpolation of the object density field into a grid.
   * @details  Interpolation of the galaxy overdensity field.
   * Periodic bounday conditions are applie to remap objects
   * with coords. outside the range [0,Lx] within the box
   * @param x  x-coordinate of galaxy
   * @param y y-coordinate of galaxy
   * @param z z-coordinate of galaxy
   * @param weight weight
   * @param dat vector containing the interpolated galaxy distribution
   */
#ifdef _USE_VECTORIZED_GRID_ASSIGNMENT_
    void grid_assignment_old(real_prec x, real_prec y, real_prec z, real_prec *weight, pic_storage **data, int, int);
#endif

 #ifdef _MASS_WEIGHT_POWER_
 void grid_assignment(real_prec x, real_prec y, real_prec z, real_prec weight, real_prec mass, vector<real_prec>&field, vector<real_prec>&field_mw);
#else
    void grid_assignment(real_prec x, real_prec y, real_prec z, real_prec weight, vector<real_prec>&field);
#endif

    void grid_assignment_NGP(real_prec x, real_prec y, real_prec z, real_prec weight, vector<real_prec>&field);
    void grid_assignment_CIC(real_prec x, real_prec y, real_prec z, real_prec weight, vector<real_prec>&field);
    void grid_assignment_TSC(real_prec x, real_prec y, real_prec z, real_prec weight, vector<real_prec>&field);
    void grid_assignment_PCS(real_prec x, real_prec y, real_prec z, real_prec weight, vector<real_prec>&field);
    //////////////////////////////////////////////////////////
  /**
   * @brief Sampling of galaxy catalogue
   * @details according to selected MAS
   * @param s_b structure of type s_parameters_box
   * @param s_d structure of type s_data_structure
   * @result Generates private class members: vectors ready to be Fourier transformed
   */
  void get_interpolated_density_field(s_data_structure *s_d);
  void get_interpolated_density_field_real_space(s_data_structure *s_d);
#ifdef _USE_VECTORIZED_GRID_ASSIGNMENT_
  void get_interpolated_density_field_old(s_data_structure *s_d);
#endif
  //////////////////////////////////////////////////////////
  /**
   * @brief Constant density grid assignment
   * @details Function used when no random catalogue is implemented.
   * A vector is filled with the mean number density of
   * the simulation. Such mean number density is computed
   * from the information if the size of the box
   * as given in the parameter file, together with the
   * number of objects.
   * @param vol Volumen of the sample if known
   * @result Private class member containing the galaxy fluctuation interpolated in the mesh of size NFT.
   */
  void raw_sampling(real_prec vol);
  //////////////////////////////////////////////////////////
  /**
   * @brief Build FKP fluctuation
   * @details Build the galaxy fluctuation by subtracting data
   * and random catalogue with factor alpha
   * @params vol Volume of the box
   */
  void get_fluctuation();
  //////////////////////////////////////////////////////////
  /**
   * @brief Computes the parameters associated to the
   * FKP estimator
   * @details e.g., normalization, shot noise, etc.
   * and assign them to public/or private variables of this class.
   * @ parameter ur Use random catalogue (true/false)
   */
  void get_parameters_estimator(bool verbose);

  //////////////////////////////////////////////////////////
  /**
   * @brief Correction for the mass assignment scheme
   * @details Cormode by mode.
   * @param  i kx-coordinate in Fourier space (in units of the fundamental mode)
   * @param  j ky-coordinate in Fourier space (in units of the fundamental mode)
   * @param  k kz-coordinate in Fourier space (in units of the fundamental mode)
   * @return Squared of the Fourier transform of the mass assignment scheme at the position (i,j,k) in Fourier space
   */
  real_prec correction_MAS(int i, int j, int k);

  //////////////////////////////////////////////////////////
  /**
   * @brief Shell average in Fourier space.
   * @details This function returns the spherical average estimate
   * of the monopole, quadrupole, hexadecapole, and the
   * window function of based on the FKP estiamtor.
   * For the spherical averages, we only use a quarter of the full FOURIER box,
   * using Hermitian symmetry F(kx, ky, -kz)=F(-kx, -ky, kz)* to
   * recover the  information in the negative frequencies explicitely.
   * When counting modes and power, we weight by a factor 2
   * order to account for the negative z-quadrant. Although
   * the factor 2 cancels out when computing the average power in each
   * shell this allows comparisons with codes using the full box.
   * This subroutine is ideal also for the multipole decomposition
   * in the modes l=0, 2 and 4. For other moments, the full excusrion
   * through Fourier space has to be done.
   * BINNING: the floor function ensures that we are using intervals
   * of the form [). The spherical shells are such that the first bin
   * starts at the zero frequency (althought that mode is exlcuded
   * for the power spectrum, not for the window function),
   * Note therefore that in the case ndel=1, the fist bin
   * will contain ONLY one Fourier mode, is the zero frequency.
   * Therefore, it is convinient to start with ndel=2
   * This function is called by get_power_spectrum_fkp().
   * @param s_b structure of type s_parameters_box
   * @return p_g0 monopole power spectrum
   * @return p_g2 quadrupole power spectrum
   * @return p_g4 hexadecapole power spectrum
   * @return p_r power spectrum of the window function
   * @return p_2d 2d power spectrum in cartesian coordinates
   * @return p_2s 2d power spectrum in polar coordinates
   * @return nm Number of modes in spherical shells
   */
  void power_spectrum_fkp(s_parameters_box *s_b, vector<real_prec>&p_g0, vector<real_prec>&p_g2, vector<real_prec>&p_g4,vector<real_prec>&p_r, vector< vector<real_prec> >&p_2c, vector< vector<real_prec> >&p_2s, vector<int>&nm);

  void power_spectrum_fkp(s_parameters_box *s_b, vector<real_prec>&p_g0,vector<int>&nm);


  //////////////////////////////////////////////////////////
  /**
   * @brief Shell average in Fourier space for Bispectrum
   * @param s_b structure of type s_parameters_box
   * @return p_g0 monopole power spectrum

   */
  void power_spectrum_fkp_for_bispectrum(s_parameters_box *s_b, vector<real_prec>&p_g0);




  //////////////////////////////////////////////////////////
  /**
   * @brief Shell average in Fourier space.
   * @details This function returns the spherical average estimate
   * of the monopole, quadrupole, hexadecapole, based on the Yamamoto estimator.
   * For the spherical averages, we only use a quarter of the full FOURIER box,
   * using Hermitian symmetry F(kx, ky, -kz)=F(-kx, -ky, kz)* to
   *  recover the  information in the negative frequencies explicitely.
   *  When counting modes and power, we weight by a factor 2
   *  order to account for the negative z-quadrant. Although
   *  the factor 2 cancels out when computing the average power in each
   *  shell this allows comparisons with codes using the full box.
   *  This subroutine is ideal also for the multipole decomposition
   *  in the modes l=0, 2 and 4. For other moments, the full excurion
   *  through Fourier space has to be done.
   *  This function is called by get_power_spectrum_yamammoto().
   * @param s_b structure of type s_parameters_box
   * @return p_g0 monopole power spectrum
   * @return p_g2 quadrupole power spectrum
   * @return p_g4 hexadecapole power spectrum
   * @return mod Number of modes in spherical shells
   */
  void power_spectrum_yamamoto(s_parameters_box *s_b, vector<real_prec>&p0, vector<real_prec>&p2, vector<real_prec>&p4,vector<int>&mod);
  void _power_spectrum_yamamoto(s_parameters_box *s_b, vector<real_prec>&p0, vector<real_prec>&p2, vector<real_prec>&p4,vector<int>&mod);



  //////////////////////////////////////////////////////////
  /**
   *@brief  Compute the multipole decomposition using FKP estimator
   * @param s_b structure of type s_params_box
   * @return p0 monopole power spectrum
   * @return p2 quadrupole power spectrum
   * @return p4 hexadecapole power spectrum
   * @return pr power spectrum of the window function
   * @return p2d TwoD power spectrum in cartesian coordinates
   * @return p2s TwoD power spectrum in polar coordinates
   * @return nm Number of modes in spherical shells
  */
  void get_power_spectrum_for_bispectrum(s_parameters_box *s_b,vector<real_prec>&p0);

  //////////////////////////////////////////////////////////
  /**
   *@brief  Compute the multipole decomposition using FKP estimator
   * @param s_b structure of type s_params_box
   * @return p0 monopole power spectrum
   * @return p2 quadrupole power spectrum
   * @return p4 hexadecapole power spectrum
   * @return pr power spectrum of the window function
   * @return p2d TwoD power spectrum in cartesian coordinates
   * @return p2s TwoD power spectrum in polar coordinates
   * @return nm Number of modes in spherical shells
  */
 void get_power_spectrum_fkp(s_parameters_box *s_b,vector<real_prec>&p0,vector<real_prec>&p2,vector<real_prec>&p4, vector<real_prec>&pr, vector<vector<real_prec> >&p2d,  vector<vector<real_prec> >&p2s,vector<int>&nm);

 void get_power_spectrum_fkp(s_parameters_box *s_b,vector<real_prec>&p0,vector<int>&nm);

 void cross_power_spectrum_fkp(s_parameters_box *s_b,vector<real_prec>&p0,vector<int>&nm);

 //////////////////////////////////////////////////////////
 /**
  *@brief  Compute the multipole decomposition using Yamamoto estimator
  *@details with the FFTW based scheme.
  * @param s_b structure of type s_params_box
  * @result p0 monopole power spectrum
  * @result p2 quadrupole power spectrum
  * @result p4 hexadecapole power spectrum
  * @result p_r power spectrum of the window function
  * @result nm Number of modes in spherical shells
  */
 void get_power_spectrum_yamamoto(s_parameters_box * s_b,vector<real_prec>&p0,vector<real_prec>&p2,vector<real_prec>&p4,vector<int>&nm);

 //////////////////////////////////////////////////////////
  /**
   * @brief Estimates of the variance for the power spectrum
   * @details based on the FKP estimator. Selects between the
   * exact and the approximate expression for the var(P)
   * @param s_b structure of type s_params_box
   * @param s_d structure of type s_data_structure
   * @param kv k-spherical shells
   * @param pk power spectrum in k-spherical shells
   * @param nm number of modes in k-spherical shells
   * @result sig FKP variance of the power spectrum
   */
  void get_fkp_error_bars(s_data_structure *s_d, vector<real_prec> &kv, vector<real_prec> &pk, vector<int>&nm, vector<real_prec> &sig);

  //////////////////////////////////////////////////////////
  /**
   * @brief Generate an interpolated value of the mean number density
   * @details given the position in the sky of the object, when the depth
   * varies with the angular position. This interpolates the matrix
   * computed on the method dndz
   * @param s_b structure of type s_parameters_box
   * @param zmin Minimum redshift of the sample
   * @param zmax Maximum redshift of the sample
   * @param ra Right ascencion of the galaxy
   * @param ra Declination of the galaxy
   * @param zg Redshift of the galaxy
   * @param dndz_m container with the vaues of the mean number density in redshift bins and Healpix pixels tabulated.
   * @result nbar Mean number denisiy tabulated at the position of the galaxy
   */
  void get_mean_density_interpolated(real_prec zmin, real_prec zmax, real_prec ra, real_prec dec, real_prec zg, vector< vector<real_prec> >&dndz_m, real_prec *nbar);



  //////////////////////////////////////////////////////////
  /**
   * @brief Count of triangles in Fourier space
   * @details
   */
  void bispectrum_fkp(char, s_parameters_box *s_p, vector<real_prec> &, vector<real_prec> &, vector<real_prec> &, vector<int> &);

  //////////////////////////////////////////////////////////
  /**
   * @brief Map indices between different DFT output schemes
   * @details This auxiliary function returns the indices that are
   * arguments of the function ijk(), used to retrieve the
   * amplitudes of the DTF as given by the FFTW. This takes into account
   * the fact that we are adding one more frequency,
   * the -Nyquist to the outpu and that we use Hermitian
   * symmetry to retrieve the negative section
   * of the third component, i.e,
   * delta(kx, ky, -kz)=delta(-kx, -ky, kz)*
   * @param s_b structure of type s_parameters_box
   */
  void remap(int, int, int, int , int , int , int *, int *, int *, real_prec *);

  //////////////////////////////////////////////////////////
  /**
   * @brief Estimates of Bispectrum
   * @details  based on the counts of triangles
   * the normalization and the shot-noise corrections.
   * @param s_b structure of type s_parameters_box
   */
  void get_bispectrum_fkp(char, s_parameters_box *s_p, vector<real_prec> &, vector<real_prec> &, vector< int > &);


  //////////////////////////////////////////////////////////
  /**
   * @brief Estimates of Bispectrum
   * @details  based on trick by Scoccimarro and Jennifer
   * @param s_b structure of type s_parameters_box
   */
  void get_bispectrum_fkp_fast(s_parameters_box *s_p, vector<real_prec> &, vector<real_prec> &, vector< int > &, string file);


  //////////////////////////////////////////////////////////
  /**
   * @brief Define k-shells for the estimation of bispectrum
   * @param s_b structure of type s_parameters_box
   */
  void define_kshells(s_parameters_box *s_box);


  //////////////////////////////////////////////////////////
  /**
   * @brief Get the inverse Fourier transform in each Fourier shell
   * @param Pass the definition of binning
   */
  void get_ift_shells_bispectrum(s_parameters_box *s_box);


  //////////////////////////////////////////////////////////
  /**
   * @brief Define k-shells for the estimation of bispectrum
   * @param s_b structure of type s_parameters_box
   */
  void loop_shells_bispectrum(s_parameters_box *s_box, vector<real_prec> &pk, vector<real_prec> &bispect, vector< int > &mod, string file);



  //////////////////////////////////////////////////////////
  /**
   * @brief Define k-shells for the estimation of bispectrum
   * @param s_b structure of type s_parameters_box
   */
  void construct_shells(int ngrid, int kmnid, int kmxid, vector<real_prec> &iFT_output_delta, vector<real_prec> &iFT_output_triangles,vector<real_prec> &iFT_output_p1_cyc_sum);


  //////////////////////////////////////////////////////////
  /**
   * @brief Mapping of vectos from one octant to the octants in the kz>0 sub-volume of Fourier space
   * @param s_b structure of type s_parameters_box
   */

  void cellsym(int id, int ngrid,complex_prec *data_ks, complex_prec *data_dk,complex_prec *data_pk_sn);


  //////////////////////////////////////////////////////////
  /**
   * @brief Evaluates the DFT required for the estimates of
   * multipole decomposition
   * @details of the power spectrum using the Yamamoto-Blake estimator
   * obtained from a direct-sum approach
   * @param s_b structure of type s_parameters_box
   * @param s_d structure of type s_data_structure
   */
  void get_power_moments_fourier_grid_ds_yam(s_data_structure *s_d);

  //////////////////////////////////////////////////////////
  /**
   * @brief Compute shell-averaged multipole decomposition of the
   * power spectrum obtained from the  Yamamoto-Blake estimator
   * implemening a a direct-sum approach.
   * @param s_b structure of type s_parameters_box
   * @result p0 monopole
   * @result p2 quadrupole
   * @result p4 hexadecapole
   * @result nm number of modes in spherical shell
   */


  void power_yam_1d_ds(vector<real_prec>&p0 ,vector<real_prec>&p2 ,vector<real_prec>&p4, vector<int> &nm);


  //////////////////////////////////////////////////////////
  /**
   * @brief Write parameters on screen
   */
  void write_fftw_parameters();


  //////////////////////////////////////////////////////////
  /**
   * @brief Write parameters on .log file
   * @param p void pointer
   * @param log_f log file
   */
  void write_fftw_parameters(void *p, string log_f);



  //////////////////////////////////////////////////////////
  /**
   * @brief Write parameters on .log file
   * @param p void pointer
   * @param log_f log file
   */
  void free_fftw_vectors();


  vector<real_prec>mass_cuts;




  //////////////////////////////////////////////////////////
  /**
   * @brief get the value of private member Number of FftwFunctions::n_ran
   * @return FftwFunctions::n_ran
   */
  ULONG _n_ran(){return n_ran; }
//////////////////////////////////////////////////////////
  /**
   * @brief get the value of private member Number of FftwFunctions::w_g
   * @return FftwFunctions::w_g
   */
  real_prec _w_g(){return w_g; }
//////////////////////////////////////////////////////////
  /**
   * @brief get the value of private member Number of FftwFunctions::w_r
   * @return FftwFunctions::w_r
   */
  real_prec _w_r(){return w_r; }
 //////////////////////////////////////////////////////////
    /**
   * @brief get the value of private member Number of FftwFunctions::alpha
   * @return FftwFunctions::alpha
   */
  /** Return alpha*/
  real_prec _alpha(){return alpha; }
  //////////////////////////////////////////////////////////
  /**
   * @brief get the value of private member Number of FftwFunctions::normal_power
   * @return FftwFunctions::normal_power
   */
  real_prec _normal_power(){return normal_power; }
  void set_normal_power(real_prec new_normal_power){this->normal_power=new_normal_power; }
  //////////////////////////////////////////////////////////
  /**
   * @brief get the value of private member Number of FftwFunctions::normal_power
   * @return FftwFunctions::normal_power
 */
  real_prec _normal_power_two(){return normal_power_two; }
  void set_normal_power_two(real_prec new_normal_power_two){this->normal_power_two=new_normal_power_two; }
//////////////////////////////////////////////////////////
  /**
   *    @brief  Return normalization of window function
   */
  real_prec _normal_window(){return normal_window; }

//////////////////////////////////////////////////////////
  /**
   * @brief Return Poisson Shot noise
   */
  real_prec _shot_noise(){return shot_noise; }

//////////////////////////////////////////////////////////
  /**
   * @brief Return Poisson Shot noise
   */
  void set_imcut(int new_imcut){this->imcut=new_imcut; }
//////////////////////////////////////////////////////////
  /**
   * @brief Return Poisson Shot noise
   */
  real_prec _Lbox_data(){return this->Lbox_data; }
//////////////////////////////////////////////////////////
  /**
   * @brief Return Poisson Shot noise
   */
  void set_Lbox_data(real_prec new_Lbox_data){this->Lbox_data=new_Lbox_data; }

//////////////////////////////////////////////////////////
  /**
   * @brief get the value of private member Number of  FftwFunctions::n_gal
   * @return FftwFunctions::n_gal
  */
  ULONG _n_gal(){return n_gal; }
//////////////////////////////////////////////////////////
  /**
   * @brief get the value of private member Number of  FftwFunctions::n_gal
   * @return FftwFunctions::n_gal
  */
  void set_n_gal(ULONG new_ngal){this->n_gal=new_ngal; }
//////////////////////////////////////////////////////////



};





#endif
