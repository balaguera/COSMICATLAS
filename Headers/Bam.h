/**
 * @class<Bam>
 * @brief Header file for the class Bam::
 * @file Bam.h
 * @title Bias Assignment method for mock catalogs
 * @author Andres Balaguera-Antolínez, Francisco-Shu Kitaura
 * @version   1.0
 * @date      2020

 *@details: This is an example of a main function to call bam. A file called cosmicatlas.cpp
 *@code

# include "../Headers/Bam.h"
# include "../Headers/FileOutput.h"
# include "../Headers/NumericalMethods.h"
# include "Tasks.h"

using namespace std;

  int main(int argc, char *argv[]){

  time_t start_all;
  time(&start_all);
  char temp;
  string par_file_bam;

  ScreenOutput So(start_all);

  if(argc==1){
    So.usage(argv[0]);
    exit(1);
  }
  while((temp =  getopt(argc, argv, "hiam:n:c:t:p:x:r:")) != -1)
    {
      if(temp=='h')
        So.usage(argv[0]);

      else if (temp=='a')  // Show authors
        So.author();

      else if(temp=='c') // Run cosmicatlas
        {

          So.message(start_all);
          par_file_bam = argv[2];
          Params Par(par_file_bam);
          Bam bam(Par);
          bam.So=So;
          bam.bamrunner();
          So.message_time(start_all);
        }
      else if(temp=='i') // displays input parameters
        {
          par_file_bam = argv[2];
          Params params(par_file_bam);
          Bam bam(params);
          bam.show_params();
        }

      else if(temp=='m')   // to Mesure power spectrum
        {
          So.message(start_all);
          par_file_bam = argv[2];
          Params params(par_file_bam);
          PowerSpectrumF cPSF(params);
          cPSF.compute_power_spectrum(true,true);
          So.message_time(start_all);
        }
      else if(temp=='t')   // to read input tracer catalg and analyze it. The operations here should be included in BAM
        {
          So.message(start_all);
          par_file_bam = argv[2];
          Params params(par_file_bam);
          PowerSpectrumF cPSF(params);
        }
     }
  }
 * @endcode


*/
 

// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************

#ifndef _Bam_
#define _Bam_

# include <iostream>
# include <fstream>
# include <vector>
# include <algorithm>
# include <math.h>
# include <sstream>
# include <iomanip>
# include <string>
# include "NumericalMethods.h"
# include "ScreenOutput.h"
# include "Params.h"
# include "PowerSpectrumF.h"
# include "PowerSpectrumTH.h"
# include "McmcFunctions.h"
# include "Patchy.h"
# include "Cwclass.h"
# include "Catalog.h"
using namespace std;


class Bam{
  
 private:

  //////////////////////////////////////////////////////////
  /**
   *  @brief ScreenOutput obejct
   */
  ScreenOutput So;
  //////////////////////////////////////////////////////////
  /**
   * @brief Patchy object to be used in Bam
   */
  PATCHY patchy;
  //////////////////////////////////////////////////////////
  /**
   * @brief Object of class Cwclass, for Cosmic web analysis. 
   */
  Cwclass cwclass;
  //////////////////////////////////////////////////////////
  /**
   * @brief Object of class Cwclass, for Cosmic web analysis. 
   */
  Cwclass cwclass_ref;
  //////////////////////////////////////////////////////////
  /**
   * @brief Object of class Params
   */
  Params params;
  //////////////////////////////////////////////////////////
  /**
   * @brief Object of class FileOutput
   */
  FileOutput File;
  //////////////////////////////////////////////////////////
  /**
   * @brief Object of class Cosmology
   */
  Cosmology Cosmo;
  //////////////////////////////////////////////////////////
  /**
   * @brief Used to compute different statistical quantities ina  likelihood context. Used in _BIAS_MODE_
   */
  McmcFunctions mcmc;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Object of type Catalog used to allocate info of BAM mocks
   */
  Catalog tracer;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Object of type Catalog used to allocate info of reeference catalog
   */
  Catalog tracer_ref;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  s_CosmoInfo s_cosmo_info;
  //////////////////////////////////////////////////////////
  /**
   * @brief Boolean type variable to specify whether we are in bam mode (mock or bias) or something else (e.g, running bam with the -c option)
   * false
   */
  bool bam_mode;
    //////////////////////////////////////////////////////////
  /**
   * @brief Use to solve the behavior of under-biased tracers. 
   * false
   */
  bool aux_kernel;
  //////////////////////////////////////////////////////////
  /**
   * @brief Use to solve the behavior of under-biased tracers
   */
  bool used_once;
  //////////////////////////////////////////////////////////
  /**
   * @brief outoput object
   */
  ofstream sal;
  //////////////////////////////////////////////////////////
  /**
   * @brief outoput object
   */
  ofstream output_res;
  //////////////////////////////////////////////////////////
  /**
   * @brief Mean gas density
   */

  real_prec Mean_rho_gas;
  //////////////////////////////////////////////////////////

  /**
   * @brief Mean density of DM
   */
  real_prec Mean_density_X;
  //////////////////////////////////////////////////////////

  /**
   * @brief Mean density of DM tracer
   */
  real_prec Mean_density_Y;
  //////////////////////////////////////////////////////////
  /**
   * @brief Mean mass-weighted density of dark mattr tracer
   */
  real_prec Mean_density_Y_MASS;
  //////////////////////////////////////////////////////////
  /**
   * @brief Mean f-weighted trcer density. f  is the fraction of satelites
   */
  real_prec Mean_density_Y_SAT_FRAC;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Auxiliary variable
   */
  real_prec mean_aux;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Poisson shot noise of the halo field generated at each iteration
   */
  real_prec shot_noise_new;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Poisson shot noise of the reference halo field
   */
  real_prec shot_noise_ref;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of properties used to characterize the halo bias
   */
  int N_bias_properties;
  //////////////////////////////////////////////////////////
 /**
  * @brief
  */
  real_prec minimum_multiscale_property;
  //////////////////////////////////////////////////////////
   /**
    * @brief function of class BAM. Load the class BAM with parameters
    */
   void warnings();
  //////////////////////////////////////////////////////////
  /**
   * @brief Load cosmological information at the input redshift
   */
  void get_cosmo();
   //////////////////////////////////////////////////////////
  /**
   * @brief Computes the DM density field at each iteration
   * by convolviong the initial DM field with the Kernel
   */
  void get_BAM_DM();
  //////////////////////////////////////////////////////////
 /**
  * @brief Computes the DM density field at eacfh iteration
  * by convolviong the initial DM field with the Kernel
  */
  void get_new_DM_field();
  //////////////////////////////////////////////////////////
 /**
  * @brief Computes the DM density field at eacfh iteration
  * by convolviong the initial DM field with the Kernel
  */
  void get_minimum_multiscale_property();

  //////////////////////////////////////////////////////////
#ifdef MOCK_MODE
  /**
   * @brief Compute the bias as the probability of one cell having a given number of
   * halos conditional to the properties of the DM density field.
   * @details This is accomplished by constructing an histogram in the different variables.
   * string counts gets the bias B(N|{theta}_dm), used to assign number counts in cells,
   * while string mass computs
   * the bias B(M|{theta}_dm), used to assign halo mass to a cell
   * @returns This function loads the class member container Bam::CWT_X_Y_hist

   */
    void get_BIAS(string);
#endif
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
#ifdef MOCK_MODE

  /**
   * @brief Compute the abundance with respect to a X property from the refernece catalog, conditional to the properties of teh DM density field.
   * @details This is accomplished by constructing an histogram in the different variables.
   * This is an alternative to the mass-assignmetn approach given by the the calculation of B(M|{theta}), where M 
   * is the total mass of tracers in cell. Here instead we compute n(Mh|{\theta}) to assing halo masses
   * (i., e, on an object by object basis). This approach  works better. DOES NOT USE THE INFORMATION FROM V-CLASS
   * @implements Dark matter prperties and Vmax
   * @returns This function loads the class member container Bam::CWT_X_Y_hist
   */
    void get_X_function();
  //////////////////////////////////////////////////////////
  /**
   * @brief Same goal as get_X_function, used when two or more references are used to learn the distribution of properties from  
  */
    void get_X_function_two(int );


    //////////////////////////////////////////////////////////
  /**
   * @brief Compute the abundance with respect to a X property from the reference catalog, conditional to the properties of the DM density field.
   * @details This is accomplished by constructing an histogram in the different variables.
   * This is an alternative to the mass-assignment approach given by the the calculation of B(M|{theta}), where M
   * is the total mass of tracers in cell. Here instead we compute n(Mh|{\theta}) to assing halo masses
   * (i., e, on an object by object basis). This approach  works better. DOES NOT USE THE INFORMATION FROM V-CLASS
   * @implements Dark matter prperties and alread inplace assignation of Vmax
   * @returns This function loads the class member container Bam::CWT_X_Y_hist
   */

  void get_X_function_complement(string);

#endif
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
#ifdef MOCK_MODE
  
  /**
   * @brief Generate the tracer mock density field using BAM
   * @details The TR number density field is generating by sampling the conditional 
   * probability distribution obtained by BAM onto the Target DM density field.
   * @param new_dm_field; bool true/false to ask whether a new target DM fiueld is used to creat mock
   * @details In general, this TDMF is another realization of the same IC used to calibrate the kernel
   * @return The output is a TR number density field with its power spectrum
   * used in the iterative procedure of BAM.
   
   */
  void get_mock_grid(string property);


#ifdef _USE_TWO_REFS_MOCKS_
  void get_mock_grid_two(string property);
#endif


#endif


#ifdef _GET_BAM_CAT_
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   * @brief This function is used if _USE_NEW_MASS_ASSIGNMENT_ is defined
  */
  void assign_tracer_property_new_new(bool initial_assignment, string h_prop);
  //////////////////////////////////////////////////////////
#endif
  //////////////////////////////////////////////////////////
  /**
   * @brief This function is used if _USE_NEW_MASS_ASSIGNMENT_ is undefined
  */
  void assign_tracer_mass_new();

  //////////////////////////////////////////////////////////
  /**
   * @brief Generate random positions of tracers
   * @details The function reads the nhumber of tracers per cell and assign random coordinates
   * @return Writes to a file three columns with x, y, z
   */
  void sample_mock();
  //////////////////////////////////////////////////////////
  /**
   * @brief Read density fields using as inputs
   * @param file_dm: path to the DM density field, (file_dm)
   * @param files_tr: path to the DM  tracer (files_tr) 
   * @param files_dm_ref_pdf: path to reference dm density field. Binary files
   * After the call of this functions, the class members
   * delta_Y, delta_X and delta_X_ini contain the information of the TR 
   * and DM DENSITY field
   */

#ifdef _USE_VELOCITIES_
  void read_bam_files(string file_dm, string files_tr, string files_tr_mw, string file_y_mass, string file_y_sf, string files_dm_ref_pdf, string file_Vx, string file_Vy, string file_Vz);
#else
  void read_bam_files(string file_dm, string files_tr, string files_tr_mw, string file_Y_mass, string file_Y_sat_frac, string files_dm_ref_pdf);
#endif
  //////////////////////////////////////////////////////////
  /**
   * @brief Analyze the density fields read by read_bam_files()
   * @details performing some basic statistics and power spectrum. If defiened as pre-processor directive (RUN_TEST), 
   * a TEST is performed here.
   * @details This function redefines the values of delta_X_min etc in case it is requeitsd in parameter file. 
   * Else, input values in parameter file are used.
   * @details The fields Bam::delta_X, Bam::delta_X_ini and Bam::delta_Y are, if requiested, converted to OVERDENSITIES.
   */
  void analyze_input();
  //////////////////////////////////////////////////////////
  /**
   * @brief Compute, if desired, or define from input file, 
   * the min and max values of variables X and Y
   */
  void get_min_max_X_Y(); // to be derpecated
  //////////////////////////////////////////////////////////
  /**
   * @brief Compute, if desired, or define from input file,
   * the min and max values of variables X and Y
   */
  void get_new_min_max_properties();

  //////////////////////////////////////////////////////////
  /**
   * @brief Structure containing the minimum values of the different DM properties used in the measurement of BIAS
   */
  s_minimums s_mins;
  //////////////////////////////////////////////////////////
  /**
   * @brief Structure containing the maximum values of the different DM properties used in the measurement of BIAS
   */
  s_maximums s_maxs;
  //////////////////////////////////////////////////////////
  /**
   * @brief Structure containing the bin sizes values of the different DM properties used in the measurement of BIAS
   */
  s_Deltas s_deltas;
  //////////////////////////////////////////////////////////
  /**
   * @brief Construction of a halo number counts field from a bias model
   * @details Given the DM overdensity, and follopwing the model pre-defined in bias() (def.h), this computs
   * the halo number counts in a mesh (defined by private variables)
   * @return aux: Halo number counts in a mesh
   */
  void dark_matter_to_halos_analytical(const vector<real_prec>& in, vector<real_prec>&aux);
  //////////////////////////////////////////////////////////
  /**
   * @brief Size of bin in Y (often tracer) used to construct the BIAS
   */
  real_prec DELTAY;

  //////////////////////////////////////////////////////////
  /**
   * @brief Size of bin in Y (often tracer) used to construct the BIAS
   */
  real_prec DELTAY_MW;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Size of bin in Y (often DM) used to construct the BIAS
   */
  real_prec DELTAX;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum halo-masses in bins of Mh
   */
  vector<real_prec> MBmin;


  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  bool read_bam_kernel;

  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  bool dm_already_done;

  //////////////////////////////////////////////////////////
  /**
   * @brief Vector of structures. Contains the values of the masses that belongs to a ginven bin in the DM properties
   */

  vector<s_mass_members> dm_properties_bins;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector of structures. Contains the values of the masses that belongs to a ginven bin in the DM properties
   */

  vector<s_mass_members> dm_properties_for_randoms_bins;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector of structures. Contains the values of the masses that belongs to a ginven bin in the DM properties
   */

  vector<s_mass_members> dm_properties_bins_mock;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the maximum halo-masses in bins of Mh
   */
  vector<real_prec> MBmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Container for the X (generally DM) density field
   * During running time it is converted to overdensity
   */
  vector<real_prec>delta_X;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector of structures. Contains the values of the masses that belongs to a ginven bin in the DM properties
   */

  vector<real_prec>Displacement_inicial;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector of structures. Contains the values of the masses that belongs to a ginven bin in the DM properties
   */
  vector<real_prec>delta_dm_aux; // this is used for mass assignmet
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector of structures. Contains the values of the masses that belongs to a ginven bin in the DM properties
   */
  vector<real_prec>delta_dm_aux_mem; // this is used for mass assignmet

  //////////////////////////////////////////////////////////
  /**
   * @brief  Container for the x-component of the velocity field of X (generally DM) density field
   */
  vector<real_prec>Velx_X;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Container for the y-component of the velocity field of X (generally DM) density field
   */
  vector<real_prec>Vely_X;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Container for the z-component of the velocity field of X (generally DM) density field
   */
  vector<real_prec>Velz_X;
  //////////////////////////////////////////////////////////
  /*
     * @brief  Container for the divergence of the DM velocity field
  */
  vector<real_prec>Divergence_VelField;
  //////////////////////////////////////////////////////////
 /**
   * @brief
   */
  void get_shear_velfield();

//////////////////////////////////////////////////////////
  
  /**
   * @brief  Container for the X (generally DM) density field of reference
   * During running time it is converted to overdensity
   */
  vector<real_prec>delta_X_REF_PDF;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector of structures. Contains the values of the masses that belongs to a ginven bin in the DM properties
   */
//  real_prec kmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the X (generally DM) density field. This is used for the calibration procedure, 
   * in which the DM density field is convolved with a kernel to match the TR P(k). In that procedure, the 
   * container Bam::delta_X is modified, so we always use the original Bam::delta_X_ini to do the convolution.
   */
  vector<real_prec>delta_X_ini;
  
  //////////////////////////////////////////////////////////
  /**
   * @brief  Container for the Y (generally TR) density field
   * @brief  During running time it is converted to overdensity, if requested in parameter file
   */
  
  vector<real_prec>delta_Y;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Container for the Y (generally TR) MASS weighted density field
   */
  vector<real_prec>delta_Y_HR;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Container for the Y (generally TR) MASS weighted density field
   */

  vector<real_prec>delta_Y_MASS;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<real_prec>delta_Y_SAT_FRAC;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the number counts or other prop used inside the get_mock function
   */
  vector<real_prec>delta_Y_new;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the REFERENCE TR power spectrum
   */
  vector<real_prec>Power_REF;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<real_prec>weight_mh;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the REFERENCE mass weighted power spectrum
   */
  vector<real_prec>Power_REF_MW;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the mock power spectrum
   */
  vector<real_prec>Power_NEW;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the mock mass-weighted power spectrum
   */
  vector<real_prec>Power_NEW_MW;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the REFERENCE DM power spectrum
   */
  vector<real_prec>Power_DM_REF;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the BAM DM power spectrum
   */
  vector<real_prec>Power_DM_NEW;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the wave vectors
   */
  vector<real_prec>kvec;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of grid cells 
   * @details computed inside Bam::Bam(params _par) as (Bam::Nft)*(Bam::Nft)*(Bam::Nft)
   */
  ULONG NGRID_low_mass_assignment;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of Grid cells N*N*(N/2+1). Computed inside BAM (constructor)
   */
  ULONG NT;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  ULONG NTT;
//////////////////////////////////////////////////////////
  /**
   * @brief Minimum value of x=, delta_x or log10(num_in_log+delta_x). Computed inside BAM
   */
  real_prec Xmin;
//////////////////////////////////////////////////////////
  /**
   * @brief Maximum value of x=, delta_x or log10(num_in_log+delta_x). Computed inside BAM
   */
  real_prec Xmax;
//////////////////////////////////////////////////////////
  /**
   * @brief Minimum value of Y=, delta_Y or log10(num_in_log+delta_Y). Computed inside BAM
   */
  real_prec Ymin;
//////////////////////////////////////////////////////////
  /**
   * @brief Maximum value of y=, delta_t or log10(num_in_log+delta_y). Computed inside BAM
   */
  real_prec Ymax;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of bins in X for BIAS. Computed inside BAM
   */
  ULONG new_nbins_x;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of bins in X for BIAS. Computed inside BAM
   */
  ULONG new_nbins_y;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of bins in X for BIAS. Computed inside BAM
   */
  ULONG new_nbins_y_MW;

  //////////////////////////////////////////////////////////
  /**
   * @brief Aux, refers to index of Mass bin, to be deprecated
   */
  int im;
  //////////////////////////////////////////////////////////
  /**
   * @brief Index specifying the CW type. 
   * @details Set to zero when mock catalog is created
   */
  int tstruct;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the BIAS in the form o Number counts, mass, and satellite fraction 
   */
  vector<ULONG>  BIAS_NCOUNTS;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the BIAS in the form o Number counts, mass, and satellite fraction 
   */
  vector<ULONG>  BIAS_SAT_FRACTION;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the BIAS in the form o Number counts, mass, and satellite fraction 
   */
  vector<ULONG>  ABUNDANCE;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the BIAS in the form o Number counts, mass, and satellite fraction 
   */
  vector<ULONG>  NCELLSperDMBIN;

  //////////////////////////////////////////////////////////
  /**
   * @brief Joint distribution normalized within each bin of DM density
   */
  vector<real_prec>  BIAS_NCOUNTS_normalized;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the BIAS in the form o Number counts, mass, and satellite fraction 
   */
  vector<real_prec>  BIAS_SAT_FRACTION_normalized;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the BIAS in the form o Number counts, mass, and satellite fraction 
   */
  vector<real_prec>  ABUNDANCE_normalized;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the mean values of Y in a bin of DM density
   */
  vector<real_prec> mean_Y;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of tracers Y obtained from the input density field
   * @details Computed in Bam::analyze_input()
   */
  real_prec N_objects_Y;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Number of DM particles obtained from the input density field
   * @details Computed in Bam::analyze_input()
   */
  real_prec N_objects_X;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Number of DM (reference) particles obtained from the input density field
   */
  real_prec N_objects_X_REF;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the PDF of tracers
   */
  vector<ULONG> PDF_NC_Y;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the PDF of DM particles
   */
  vector<ULONG> PDF_NC_X;
 //////////////////////////////////////////////////////////
  /**
   * @brief Number of objects in mock 
   */ 
  ULONG Nobjects;

  //////////////////////////////////////////////////////////
  /**
   * @brief Maximum number of X particles in one cell
   */
  int nmax_x_onecell;
  //////////////////////////////////////////////////////////
  /**
   * @brief Maximum number of Y particles in one cell
   */
  int mean_number_y_onecell;
  //////////////////////////////////////////////////////////
  /**
   * @brief Maximum number of Y particles in one cell
   */
  int nmax_y_onecell;
  //////////////////////////////////////////////////////////
  /**
   * @brief Maximum number of Y particles in one cell
   */
  int nmax_y_sat_onecell;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the bins in X
   */
  vector<real_prec> X_bins;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Container for the bins in X
   */
  vector<real_prec> Y_bins;
  //////////////////////////////////////////////////////////
  /**
   * @brief This vector will keep the maximum number of particles in one cell for the different cwt asked in par file
   */
  vector<int>nmax_y_onecell_cwt;
  //////////////////////////////////////////////////////////

  /**
   * @brief Generate the PDF of Halos given the DM density field
   * @details According to the type of CW requested in parameter file, this function generates, for each CWT
   * the joint and the conditional probablity distribution. This function is not used when creating a mock catalog.
   * This function is active as long as the pre-processor directive define BIAS is found
   * @details Can only be used after having called read_bam_files() and analyze_input()
   * @return The output if this function are files containing
   * @return The mean Y , var Y for bins in X, for a given bin of MK and a given CWT
   * @return All values of Y in a bin of X for a given Mk and CWT
   */

  // For different realizations
#ifdef _SEVERAL_REAL_
  void get_pdf(int);
#elif !_SEVERAL_REAL_
  // Same but for the single files from the input par
  void get_pdf();
#endif
  //////////////////////////////////////////////////////////
  /**
    * @brief Container loaded by the function get_fof_info()
    * containing the bin in Mk in which a given cell is found, according to the results of the FoF.
    */
  vector<ULONG> SKNOT_M_info;
  //////////////////////////////////////////////////////////
  /**
   * @brief Index of the BAM Iteration
   * @details Used to identify output files
   */
  int step;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */

  int step_mass_assignemt;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the BAM isotropic Kernel in Fourier space, 
   * @details Computed as the ratio between the reference power and the mock power, with size Bam::Nft/Bam::ndel
   */
  vector<real_prec>Kernel;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<real_prec>Kernel_mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for power spectrum used in the construction of the BAM kernel
   * @details with size Bam::Nft/Bam::ndel
   */
  vector<real_prec>power_ratio_unsmoothed;
  //////////////////////////////////////////////////////////
  /**
   * @brief This function builds the kernel from the ratio
   * @brief of the ref_tr power and the mock power spectrum
   * @returns Container Bam::Kernel
   */
  void GetKernel(bool, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the pdf of the reference DM field, used when 
   * rank ordering is applied to the DMF used to calibrate the Kernel
   * based on the PDF of a high resolution DM field.
   * This should be depracated, since we should start from one realization of the low res
   * and apply the kernel to ahnother realziation of the low res
   */
  vector<real_prec> pdf_ref;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the pdf of the original DM field
   * This is modified iteration after iteration
   */
  vector<real_prec> pdf_ite;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the pdf of the original DM field
   * Copied from the temporary container pdf_in filled in the iteration 0
   */
  vector<real_prec> pdf_ini;
  //////////////////////////////////////////////////////////
  /**
   * @brief Convolution of the initial DM field with the BAM Kernel
   * @params in: original DM density field
   * @params out: convolved with kernel
   * @returns convolved DM with Kernel, to be used in the determination of the bias inside the BAM iterative process
   */
  void Konvolve(vector<real_prec> &in, vector<real_prec>&out);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void v_konvolve(vector<real_prec> &);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void set_Fourier_vectors();


  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void get_power_spectrum(string type);

  //////////////////////////////////////////////////////////
  /**
   * @brief
   * @brief Read from parameter file
   
   */
  string new_Name_Property_X;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  string new_Name_Property_Y;
  //////////////////////////////////////////////////////////
  /**
   * @brief Read from parameter file
   */
  string Type_of_File_X;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  string Type_of_File_Y;

  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  s_CosmologicalParameters s_cosmo_pars;

  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  bool bin_accumulate_borders;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  bool use_iteration_ini;
  //////////////////////////////////////////////////////////

#ifdef _USE_GNUPLOT_
  /**
   * @brief Objects of Gnuplot type, yused to plot @running time.
   */
  Gnuplot gp_kernel;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  Gnuplot gp_power;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  Gnuplot gp_ratio;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  Gnuplot gp_pdf;

  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  Gnuplot gp_abundance;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  Gnuplot gp_abundance_v;

#endif

// *****************************************************************************************************************************************
// *****************************************************************************************************************************************
// *****************************************************************************************************************************************
// *****************************************************************************************************************************************
  

 public:
  //////////////////////////////////////////////////////////
  /**
   *  @brief default Class Constructor
   *  @return object of class Bam
   */
  Bam(){
      if(true==WARNING_MODELS){
      So.message_warning("More than one model for properties of DM has been allowed. Please check def.h. Cosmicatlass stops here");
      exit(1);
    }

  }
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief Use this constructor to pass an object of type Params
   *  @brief This constructor initializes the value of BAM parameters
   *  @return object of class Bam
   */
  
  Bam(Params _par, bool _bam_mode):params (_par), bam_mode (_bam_mode)
  {
    sal.precision(8);
    sal.setf(ios::showpoint); 
    sal.setf(ios::scientific); 
    this->s_cosmo_pars=this->params.s_cosmo_pars;
    if(true==WARNING_MODELS){
      So.message_warning("More than one model for properties of DM has been allowed. Please check def.h. Cosmicatlass stops here");
      exit(1);
    }

#ifdef _USE_MULTISCALE_LEVEL_4_
    this->NGRID_low_mass_assignment=(params._Nft_low_l4())*(params._Nft_low_l4())*(params._Nft_low_l4());
#endif
    this->NTT=this->params._NGRID_h();

#ifdef _GET_BAM_REALIZATIONS_
#ifdef _MULTISCALE_
    this->get_minimum_multiscale_property();
#endif
#endif

#ifdef USE_GNUPLOT_
    this->gp_power<<"set border linewidth 1.5\n";
    this->gp_kernel<<"set border linewidth 1.5\n";
    this->gp_pdf<<"set border linewidth 1.5\n";
#endif



  this->So.enter(__PRETTY_FUNCTION__);
#ifdef _FULL_VERBOSE_
  So.message_screen("Loading parameters for BAM. THIS HAS TO BE UPDATED> USE THEM DIRECTLY FROM THE PARAMS CLASS");
#endif

  this->used_once = true;


  // --------------------------------------------------------------------------------------------
  if(true==this->bam_mode)
  {
#ifdef _GET_BAM_REALIZATIONS_
    ifstream nxo;
    nxo.open(file_one_cell);
#ifdef _FULL_VERBOSE_
    if(nxo.is_open())
      So.message_screen("Reading maximum number of reference tracers in cells from file ",file_one_cell);
    else
     {
        So.message_screen("File with max number of tracers in one cell does not exist",file_one_cell);
        So.message_screen("BAM stops here");
        exit(0);
      }
#endif
  int iaux; int iaux2;
  nxo>>iaux>>iaux2;
  this->nmax_y_onecell=iaux;
  this->mean_number_y_onecell=iaux2;
  nxo.close();

#ifdef _FULL_VERBOSE_
  if(iaux==0)
    So.message_warning("Maximum number of cells = 0. End");
  else
    So.message_screen("maximum number of reference tracers in cells =",iaux);
  So.message_screen("mean number of reference tracers in cells =",iaux2);
#endif
#endif
  }

  if(true==this->bam_mode)

  this->dm_already_done=false;
  // Feed the structure for cosmological parameters
  // This is also done in patchy, so verify that these lines are also in its init par member
   


#ifdef _BIN_ACCUMULATE_
  this->bin_accumulate_borders = true;
#else
  this->bin_accumulate_borders = false;
#endif
  So.DONE();
  // Determine the number of properties used to characterize the halo bias
  this->N_bias_properties=1;  // This accounts for the dark amtter density.
  vector<string> bias_properties;
  bias_properties.push_back("Local density");
#ifdef _USE_MASS_KNOTS_
  this->N_bias_properties++;
  bias_properties.push_back("Knot-Mass");
#endif
#ifdef _USE_CWC_
  this->N_bias_properties++;
  bias_properties.push_back("T-WEB");
#endif

#ifdef _USE_VEL_KNOTS_V_
  this->N_bias_properties++;
  bias_properties.push_back("V-dispersion in knots");
#endif
#ifdef _USE_CWEB_V_
  this->N_bias_properties++;
  bias_properties.push_back("V-WEB");
#endif
#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_DELTA2_)
  this->N_bias_properties++;
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
  bias_properties.push_back("I-2");
#else
  bias_properties.push_back("ð²");
#endif
#endif
#if defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_DELTA3_)
  this->N_bias_properties++;
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
  bias_properties.push_back("I-3");
#else
  bias_properties.push_back("ð³");
#endif
#endif
#if defined (_USE_INVARIANT_TIDAL_FIELD_IV_) || defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_ELLIPTICITY_)|| defined (_USE_PROLATNESS_)|| defined (_USE_S2_)
  this->N_bias_properties++;
#ifdef _USE_INVARIANT_TIDAL_FIELD_IV_
  bias_properties.push_back("I-4");
#elif defined (_USE_TIDAL_ANISOTROPY_)
  bias_properties.push_back("Tidal anisotropy");
#elif defined (_USE_ELLIPTICITY_)
  bias_properties.push_back("Ellipticity");
#elif defined (_PROLATNESS_)
    bias_properties.push_back("Prolatness");
#else
  bias_properties.push_back("S²");
#endif
#endif
#if defined (_USE_INVARIANT_SHEAR_VFIELD_I_) || defined(_USE_NABLA2DELTA_)
  this->N_bias_properties++;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
  bias_properties.push_back("IV-1");
#else
  bias_properties.push_back("Nabla² ð");
#endif
#endif
#if defined (_USE_INVARIANT_SHEAR_VFIELD_II_) || defined (_USE_S2DELTA_)
  this->N_bias_properties++;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
  bias_properties.push_back("IV-2");
#else
  bias_properties.push_back("s²ð");
#endif
#endif
#if defined (_USE_INVARIANT_SHEAR_VFIELD_III_) || defined (_USE_S3_)
  this->N_bias_properties++;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
  bias_properties.push_back("IV-3");
#else
  bias_properties.push_back("s²");
#endif
#endif
#ifdef _FULL_VERBOSE_
#ifndef _ONLY_PATCHY_
#ifdef _DO_BAM_CALIBRATION_
  So.message_screen("****************************************************************************");
  So.message_screen("BAM is using", this->N_bias_properties,"properties to characterise the bias:");
  for(int i=0;i<bias_properties.size();++i)std::cout<<YELLOW<<bias_properties[i]<<RESET<<endl;
  So.message_screen("****************************************************************************");
  std::cout<<endl;
#endif
#endif
#endif
  this->step_mass_assignemt=0; //initialize it
  }
// *****************************************************************************************************************************************
// *****************************************************************************************************************************************
// *****************************************************************************************************************************************

  //////////////////////////////////////////////////////////
  /**
   *  @brief Default constructor
   *  @return object of class Bam
   */
  ~Bam(){}

  //////////////////////////////////////////////////////////
  /**
   * @brief Runs BAMS. Called from main function
   */
  void bamrunner();
  //////////////////////////////////////////////////////////
  /**
   * @brief Name of the file with mock number density field
   */
  string fnameMOCK;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  string file_residuals;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  vector<real_prec>residuals_power;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */

  vector<real_prec>residuals_power_unsigned;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  vector<int>it_power;

  //////////////////////////////////////////////////////////
  /**
   * @brief
   * In this function we asign coordinates from DM + random particles based on the mock number counts
   * @details Masses are also assigned and the collapse of randoms towards dm (to correct small scale clusterin) is also performed after that
   * Randoms are decided to be collapsed after having assigned mass, for the information on the mass can be used to model the
   * fraction of distance to aproach teh rans to their closest DM particle.
   * Catalog is then written with positions, velocities and masses.
  */

#ifdef MOCK_MODE
#ifdef _USE_OMP_
 void makecat(string stradd,string fnameMOCK,gsl_rng ** gBaseRand,int ir);
#else
  void makecat(string stradd,string fnameMOCK,gsl_rng * gBaseRand,int ir);
#endif

#endif
  //////////////////////////////////////////////////////////
#ifdef MOCK_MODE
  /**
   * @brief
  */
  void collapse_randoms();
#endif
  //////////////////////////////////////////////////////////
#ifdef MOCK_MODE
  /**
   * @brief
  */
  void collapse_randoms_isolated();
#endif
  //////////////////////////////////////////////////////////
#ifdef MOCK_MODE
  /**
   * @brief
  */
  void assign_tracer_mass();
#endif
  //////////////////////////////////////////////////////////
#ifdef MOCK_MODE
  /**
   * @brief
  */
  void move_particles_lpt();
#endif
  //////////////////////////////////////////////////////////
#ifdef MOCK_MODE
  /**
   * @brief
  */
  void correct_for_exclusion(ULONG lenght_dm_bin);
#endif
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  vector<s_nearest_cells> ncells_info;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  vector<gsl_real>mass_bins;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  vector<gsl_real>mfunc;
   //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec prop_min;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec prop_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec min_halo_separation;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  void set_So(ScreenOutput new_So){this->So=new_So;}
  
};
  
#endif
  
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################

