// **************************************************************************************************
// **************************************************************************************************
/**
* @class<Catalog>
* @brief Header file for the class Catalog::
* @file Catalog.h
* @title Functions related to reading and analysing input catalogs
* @details Bias Assignment method for mock catalogs
* @author Andres Balaguera-Antolínez, Francisco-Shu Kitaura
* @version 1.0
* @date    2020
*/
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************

#ifndef _CATALOG_
#define _CATALOG_

#include <vector>
#include <math.h>
#include <fstream>
#include "bstream.h"
#include <iomanip>
#include <iostream>
#include <string>         
#include <string.h>
#include <cassert>
#include <cfloat>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <omp.h>
#include <sstream>
#include <netinet/in.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>


# include "NumericalMethods.h"
# include "ScreenOutput.h"
# include "FileOutput.h"
# include "Params.h"
//# include "Catalog.h"
# include "massFunctions.h"

using namespace std;



class Catalog{
  
 private:

  //////////////////////////////////////////////////////////
  /**
   * @brief  Mean mass of the sample
   */
  real_prec mean_Mass;

  //////////////////////////////////////////////////////////
  /**
   * @brief  vector containing the input parameters
   */
  vector<real_prec> xgal;
  //////////////////////////////////////////////////////////
  /**
   * @brief  vector containing the input parameters
   */
  vector<real_prec> ygal;
  //////////////////////////////////////////////////////////
  /**
   * @brief  vector containing the input parameters
   */
  vector<real_prec> zgal;
  //////////////////////////////////////////////////////////
  /**
   * @brief  vector containing the input parameters
   */

  vector<real_prec> Pxgal; //P meant for vector properties such as velocity or momentum
  //////////////////////////////////////////////////////////
  /**
   * @brief  vector containing the input parameters
   */
  vector<real_prec> Pygal;
  //////////////////////////////////////////////////////////
  /**
   * @brief  vector containing the input parameters
   */
  vector<real_prec> Pzgal;
  //////////////////////////////////////////////////////////
  /**
   * @brief  vector containing the input parameters
   */
  vector<real_prec> VPgal;
  //////////////////////////////////////////////////////////
  /**
   * @brief  vector containing the input parameters
   */
  vector<real_prec> Mass;  // other scalar property such as the mass
  //////////////////////////////////////////////////////////
  /**
   * @brief  vector containing the input parameters
   */
  vector<real_prec> property;  // other scalar
  //////////////////////////////////////////////////////////
  /**
   * @brief  vector containing the input parameters
   */
  vector<real_prec> weight;
  
  //////////////////////////////////////////////////////////
  /**
   * @brief  File
   */
  FileOutput File;
  //////////////////////////////////////////////////////////
  /**
   * @brief  File
   */
  ScreenOutput So;
  //////////////////////////////////////////////////////////
  /**
   * @brief  File
   */
  string Output_directory;
  /**
   * @brief  Class containing the input parameters
   */
  s_Galaxy galaxy;
  
  
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of mass bins to divide the sample 
   */
  int NMASSbins;


  //////////////////////////////////////////////////////////
  /**
   * @brief Number of mass bins to measure mass function
   */
  int NMASSbins_mf;

  //////////////////////////////////////////////////////////
  /**
   * @brief Number of mass bins to measure mass function
   */
  
  int NMBINS;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of mass bins to measure mass function
   */

  Gnuplot gp_abundance;
  //////////////////////////////////////////////////////////
  /**
   * @brief Identifies dark matter "DM" or tracer "TR"
   */
  string type_of_object;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Class containing the input parameters
   */
  s_params_box_mas box;
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
  s_params_box_mas box_low;

  //////////////////////////////////////////////////////////
  /**
   * @brief Number of objects load by the read file or read_file_bin objects
   */
  ULONG NOBJS;

public:
  //////////////////////////////////////////////////////////
  /**
   * @brief Default constructor
   */

  Catalog():logdeltaRS(0),logdeltaVMAX(0),logdeltaM(0),logdeltaM_low(0),var_prop(0),mean_Mass(0),mean_number_density(0),min_halo_separation(0),NOBJS(0),Ntracers_ran(0),Ntracers_dm(0), aux_flag(true)
  {
      time_t time_bam;
      time(&time_bam);
      this->So.initial_time=time_bam;
  }
  //////////////////////////////////////////////////////////
  /**
   * @brief  Constructor
   */
  
 Catalog(Params _params):params (_params), logdeltaRS(0),logdeltaVMAX(0),logdeltaM(0),logdeltaM_low(0),var_prop(0),mean_Mass(0),mean_number_density(0),min_halo_separation(0),NOBJS(0),Ntracers_ran(0),Ntracers_dm(0),aux_flag(true)
  {
    // This constructor can be used when the type is declared
    // and used directly. Otherwise, if a class has an object of
    // thi type, we define the object calling the default
    // constructor and calling set_params_catalog to assign values to this class, specifçically, to a structure
    
    this->type_of_object=this->params._type_of_object();

    this->Output_directory=this->params._Output_directory();
    this->box.masskernel=this->params._masskernel();
    this->box.Lbox=this->params._Lbox();
    this->box.Nft=this->params._Nft();
    this->box.NGRID=this->params._NGRID();
    this->box.d1=this->params._d1();
    this->box.d2=this->params._d2();
    this->box.d3=this->params._d3();
    this->box.min1=this->params._xllc();
    this->box.min2=this->params._yllc();
    this->box.min3=this->params._zllc();

#ifdef _USE_MULTISCALE_LEVEL_4_
     this->box_low.masskernel=this->params._masskernel();
     this->box_low.Lbox=this->params._Lbox();
     this->box_low.Nft=this->params._Nft_low_l4();
     this->box_low.d1=this->params._d1_low();
     this->box_low.d2=this->params._d2_low();
     this->box_low.d3=this->params._d3_low();
     this->box_low.min1=this->params._xllc();
     this->box_low.min2=this->params._yllc();
     this->box_low.min3=this->params._zllc();
#endif
     time_t time_bam;
     time(&time_bam);
     this->So.initial_time=time_bam;


 }
  //////////////////////////////////////////////////////////
  /**
   * @brief  Default constructor
   */
  
  ~Catalog(){}  
  //////////////////////////////////////////////////////////
  /**
   * @brief Main function used to analyze catalog
   */
   void analyze_cat();
   //////////////////////////////////////////////////////////
   /**
    * @brief Writes catalog with N_PROP columns (see def.h) to a binary file
    */
   void write_catalog_bin(string output_file);

   //////////////////////////////////////////////////////////
   /**
    * @brief Writes catalog with N_PROP columns (see def.h) to an ascii file
    */
   void write_catalog(string output_file);

   //////////////////////////////////////////////////////////
  /**
   * @brief Reads particle catalog in ascii format.
   */
#ifdef _USE_MASS_CUTS_PK_
  void read_catalog(string input_file, real_prec mcut);
#elif defined (_USE_MASS_BINS_PK_)
  void read_catalog(string input_file, real_prec m_min, real_prec m_max);
#else
void  read_catalog(string input_file, real_prec aux);
#endif

//////////////////////////////////////////////////////////
  /**
   * @brief Read input catalog in binary format, with x, y, z (and velocities if requested) in three separate files of lenght n
   */

#ifdef _USE_VELOCITIES_
  void read_catalog_bin(ULONG n, string fx, string fy, string fz,string fvx, string fvy, string fvz); // reads binary from fastpm
#else
  void read_catalog_bin(); // reads binary from fastpm
#endif

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   * @brief Get the interpolated prop-density field and write to output_file
   * @bried prop can be _MASS_, _COUNTS_. TBD: must be generalized
   */

  void get_density_field_grid(string prop, string output_file);
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   * @brief Get the interpolated prop-density and write to vector
   */

  void get_density_field_grid(string prop,vector<real_prec>&);

  //////////////////////////////////////////////////////////
  /**
   * @brief
   */

  void get_property_function(string );
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void define_property_bins();


  void get_distribution_reduced_mass_in_cell();



  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void get_distribution_min_separations(vector<s_nearest_cells> &nearest_cells_info);
  //////////////////////////////////////////////////////////
   /**
    * @brief
    */

  void get_min_separation_in_cell();

  //////////////////////////////////////////////////////////
   /**
    * @brief
    */

  void get_neighbour_tracers(vector<s_nearest_cells> &nearest_cells_info );
  //////////////////////////////////////////////////////////
   /**
    * @brief
    */
  void get_masses_of_pairs_in_min_separation_bin_in_theta_bin(real_prec min_sep, vector<s_mass_members> & dm_properties_bins);


  //////////////////////////////////////////////////////////
   /**
    * @brief This function computs occupation numbers in cells
    * @details This function computs occupation numbers in cells in bins of a given property
    */
  void get_pdf_vmax(string);


  //////////////////////////////////////////////////////////


  void get_sep_distribution(real_prec min_value_prop);



  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  vector<s_Halo> Halo;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  vector<s_Halo> tracer;
   //////////////////////////////////////////////////////////
  /**1
   * @brief  
   */
  real_prec logdeltaM;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of mass bins to measure mass function
   */
  real_prec logdeltaM_low;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of mass bins to measure mass function
   */
  real_prec logdeltaVMAX;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of mass bins to measure mass function
   */
  real_prec logdeltaRS;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of mass bins to measure mass function
   */
  real_prec logdeltaSPIN;
  //////////////////////////////////////////////////////////
  /**
   * @brief Retreaves x-coordinate of i-th galaxy in catalog 
   */
  
  real_prec x(ULONG i){return this->Halo[i].coord1;} 
  //////////////////////////////////////////////////////////
  /**
   * @brief Retreaves y-coordinate of i-th galaxyin catalog 
   */
  
  real_prec y(ULONG i){return this->Halo[i].coord2;} 
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  real_prec z(ULONG i){return this->Halo[i].coord3;} 
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  real_prec vx(ULONG i){return this->Halo[i].vel1;} 
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  real_prec vy(ULONG i){return this->Halo[i].vel2;} 
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  real_prec vz(ULONG i){return this->Halo[i].vel3;} 
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  real_prec mass(ULONG i){return this->Halo[i].mass;} 
	    
  //////////////////////////////////////////////////////////
  /**
   * @brief  sample variance of the tracer property
   */
  real_prec var_prop;
  

 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */ 
  vector<gsl_real> MBmin; // +1 to include the full sample
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
  vector<gsl_real> VMAXBmin; // +1 to include the full sample
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
  vector<gsl_real> RSBmin; // +1 to include the full sample
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
  vector<gsl_real> SPINBmin; // +1 to include the full sample
 //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  vector<gsl_real> MBmax;
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
  vector<gsl_real> VMAXBmax;
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
  vector<gsl_real> RSBmax;
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
  vector<gsl_real> SPINBmax;
  //////////////////////////////////////////////////////////
   /**
    * @brief
    */
  vector<gsl_real> MBin;
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
  vector<gsl_real> VMAXBin;
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
  vector<gsl_real> RSBin;
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
  vector<gsl_real> SPINBin;

  //////////////////////////////////////////////////////////
   /**
    * @brief
    */
   vector<real_prec> Dist_Min_Sepatations;
   //////////////////////////////////////////////////////////
    /**
     * @brief
   */
   vector<int>Number_of_neighbours;
   //////////////////////////////////////////////////////////
    /**
     * @brief
     */

   vector<real_prec>local_clustering;


  //////////////////////////////////////////////////////////
  /**
   * @brief  Class containing the input parameters
   */
  int NCOLS;
  //////////////////////////////////////////////////////////
   /**
    * @brief Numbner density of the catalog
    */
  real_prec mean_number_density;

  //////////////////////////////////////////////////////////
   /**
    * @brief Container for the mass function of the catalog
    */
 vector<real_prec> mass_function;
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */

 vector<real_prec> vmax_function;
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */

 vector<real_prec> rs_function;
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */

 vector<real_prec> s_function;

 //////////////////////////////////////////////////////////
  /**
   * @brief
   */

   vector<real_prec>min_separation_in_cell;
 //////////////////////////////////////////////////////////
  /**
   * @brief
   */
 real_prec min_halo_separation;
 //////////////////////////////////////////////////////////
  /**
   * @brief
   */

 vector<s_dist_in_dmbins> masses_in_cells_min_sep;

 //////////////////////////////////////////////////////////
 /**
  * @brief  Class containing the input parameters
  */
 Params params;

 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
 real_prec Prop_threshold;

 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
 real_prec Prop_threshold_multi_scale_1;
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
 real_prec Prop_threshold_multi_scale_2;
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
 real_prec Prop_threshold_multi_scale_3;
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
 real_prec Prop_threshold_multi_scale_4;
 ///////////////////////////////////////////////////////
  /**
   * @brief  Number of properties in the Level 0 of multi-scaling
   */
 ULONG N_props_0;
 ///////////////////////////////////////////////////////
  /**
   * @brief  Number of properties in the Level 1 of multi-scaling
   */
 ULONG N_props_1;
 ///////////////////////////////////////////////////////
  /**
   * @brief  Number of properties in the Level 2 of multi-scaling
   */
 ULONG N_props_2;
 ///////////////////////////////////////////////////////
  /**
   * @brief  Number of properties in the Level 3 of multi-scaling
   */
 ULONG N_props_3;
 ///////////////////////////////////////////////////////
  /**
   * @brief  Number of properties in the Level 4 of multi-scaling
   */
 ULONG N_props_4;
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
 ULONG Ntracers_ran;
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
 ULONG Ntracers_dm;
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
 real_prec fraction_tracer_from_dm;
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
 real_prec fraction_tracer_from_random;
 ///////////////////////////////////////////////////////
  /**
   * @brief  
   */
 real_prec Prop_threshold_rand_dm;
///////////////////////////////////////////////////////
  /**
   * @brief  
   */

 bool aux_flag;
  //////////////////////////////////////////////////////////
  /**
   * @brief Identifies dark matter "DM" or tracer "TR"
   */
  void set_type_of_object(string new_type_of_object){this->type_of_object=new_type_of_object;}

 //////////////////////////////////////////////////////////
  /**
   * @brief 
   */

  ULONG _NOBJS(){return this->NOBJS;}
  void set_NOBJS(ULONG new_NOBJS) {this->NOBJS=new_NOBJS;}
 //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void set_params(Params new_params){
    this->params=new_params;
    this->Output_directory=params._Output_directory();
    this->box.masskernel=this->params._masskernel();
    this->box.Lbox=this->params._Lbox();
    this->box.Nft=this->params._Nft();
    this->box.NGRID=this->params._NGRID();
    this->box.d1=this->params._d1();
    this->box.d2=this->params._d2();
    this->box.d3=this->params._d3();
    this->box.min1=this->params._xllc();
    this->box.min2=this->params._yllc();
    this->box.min3=this->params._zllc();
    this->type_of_object=this->params._type_of_object();

#ifdef _USE_MULTISCALE_LEVEL_4_
     this->box_low.masskernel=this->params._masskernel();
     this->box_low.Lbox=this->params._Lbox();
     this->box_low.Nft=this->params._Nft_low_l4();
     this->box_low.d1=this->params._d1_low();
     this->box_low.d2=this->params._d2_low();
     this->box_low.d3=this->params._d3_low();
     this->box_low.min1=this->params._xllc();
     this->box_low.min2=this->params._yllc();
     this->box_low.min3=this->params._zllc();
#endif


  }

};


#endif
