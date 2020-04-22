// **************************************************************
/**
 * @class<Patchy>
 * @brief Header file for the class Patchy::
 * @file Patchy.h
 * @title Patchy: Approximated gravity solvers
 * @author Francisco-Shu Kitaura, Andres Balaguera-Antolínez
 * @version   1.0
 * @date      2020
*/
// **************************************************************
// **************************************************************
// **************************************************************


 
#ifndef _PATCHY_
#define _PATCHY_


#include "def.h"
#include "bstream.h"
#include <math.h>
#include <fstream>
#include <iomanip>
#include <iostream>
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

#include "FileOutput.h"
#include "NumericalMethods.h"
#include "../Headers/massFunctions.h"

using namespace std;


//##################################################################################

// *****************************************************************************************************
/*!\class PATCHY
  
 * @details   PATCHY
 * @author    Francisco-Shu Kitaura
 * @author    adapted by Andres Balaguera-Antolinez 

 * @brief CLASS PATCHY
 */


//##################################################################################

class PATCHY{

  private:
  
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /* /\** */
  /*  * @brief outoput object */
  /*  *\/ */
  /* ofstream sal; */
  
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   * @brief Object of class FileOutput
   */
  FileOutput File;

  //////////////////////////////////////////////////////////
  /**
   * @brief Input parameter file adapted to Patchy. Ideally, to bemerged with that of BAM
   */
  string par_file;


  string ic_WC_dir;
  string ic_file;

  bool use_ic_file;


  
  /**
   *  @brief ScreenOutput obejct
   */
  ScreenOutput So;

  

  //////////////////////////////////////////////////////////
  /**
   * @brief Structure containing parameters used in patchy
   * this is to be deprecated once homogeneization is complete
   */

  //  struct DATA *data;

  
  int inputmode;
  int seed;
  int seed_ref;
  bool runsim;
  bool runv;
  bool diffcosmorz;
  string Output_directory;
  string ic_power_file;
  string ic_WN_file;
  string ic_WN_dir;
  string dataFileName;
  string fastpmpath;
  bool lognden;
  string dir;
  bool transf;
  real_prec slength;
  real_prec slengthv;
  real_prec vslength;
  real_prec velbias;
  int Nchunk;
  int masskernel;
  real_prec dkbin;
  int N_bin;
  real_prec biasE;
  real_prec biasepsilon; 
  real_prec biasrhoexp;
  real_prec biasone;
  real_prec biassign;
  real_prec biassign2;
  real_prec biasepsilon2;
  real_prec biasrhoexp2;
  real_prec devpois;
  real_prec deltathH;
  real_prec Nmean;
  real_prec deltath;
  real_prec cs2;
  real_prec cs3;
  real_prec cst;
  real_prec cpsi;
  real_prec cdeltas2;
  real_prec biasL;
  real_prec sfac;
  real_prec ep;
  real_prec xllc;
  real_prec yllc;
  real_prec zllc;
  real_prec xobs;
  real_prec yobs;
  real_prec zobs;
  real_prec d1;
  real_prec d2;
  real_prec d3;
  real_prec L1;
  real_prec L2;
  real_prec L3;
  ULONG N1;
  ULONG N2;
  ULONG N3;
  bool readPS;
  ULONG NGRID;
  bool planepar;
  Params params;
  real_prec Redshift_initial;
  bool Normalize_initial_redshift;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<real_prec>tab;
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  
 public:

  /**
   *  @brief default constructor
   *  @brief object of class PATCHY. I use this to defined a PACHY objects as a BAM class member
   *  @brief
   */
  PATCHY(){
  }
  
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   *  @brief constructor 
   *  @param parameters_file parameter file
   *  @return object of class PATCHY
   */
  
  
  
  
  
  
 PATCHY(Params _params, s_CosmoInfo _s_cosmo_info) :s_cosmo_info(_s_cosmo_info)
  {
   this->set_params_patchy(_params);
   this->NGRID=static_cast<ULONG>(this->N1*this->N2*this->N3);
 }

 //////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////
  /**
     *  @brief Default destructor
     *  @return
   */
  ~PATCHY(){}

 //////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////
 /**
  * @brief Reading and initializing structures with params for patchy
  */

  void set_params_patchy(Params pars);
 //////////////////////////////////////////////////////////
 /**
  * @brief Set names of files. Strings are private class members
  */
  void set_fnames();
  //////////////////////////////////////////////////////////

  /**
   *  @brief Generate DM density field according to some SF model and IC
   */
 
#ifdef OMPPARRAN
 void get_dm_field(gsl_rng ** gBaseRand);
#else
  void get_dm_field(gsl_rng * gBaseRand);
#endif


  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  void  displace_part_comp(vector<real_prec> &posi, vector<real_prec>&psii, bool periodic,int comp);

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  void Lag2Eul_comp(real_prec biasLAG2,real_prec kth,int ftype,bool periodic,vector<real_prec> &out,int comp);

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  void Lag2Eul_compB(real_prec kth,int ftype,bool periodic,vector<real_prec> &out,int comp);

  //////////////////////////////////////////////////////////
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief Patchy routine to generate catalog with positions and velocities
   */

#ifdef OMPPARRANRSD
  void makecat(string stradd,string fnameMOCK, ULONG NLL,ULONG NOBJ,gsl_rng ** gBaseRand, int ir);
#else
  void makecat(string stradd,string fnameMOCK, ULONG NLL,ULONG NOBJ,gsl_rng * gBaseRand, int ir);
#endif




  

#ifdef OMPPARRANRSD
  void makecat_withv(string stradd,string fnameMOCK, ULONG NLL,ULONG NOBJ,gsl_rng ** gBaseRand, int ir);
#else
  void makecat_withv(string stradd,string fnameMOCK, ULONG NLL,ULONG NOBJ,gsl_rng * gBaseRand, int ir);
#endif

  void makecat_withv_new(string stradd,string fnameMOCK, ULONG NLL,ULONG NOBJ,gsl_rng ** gBaseRand, int ir);

  
  //////////////////////////////////////////////////////////
  void read_tabulated_power();
  
  
  real_prec linInterp(real_prec xpos, real_prec ypos, real_prec zpos, const vector<real_prec>&);

  //////////////////////////////////////////////////////////
  /**
   * @brief Object of class Cosmology
   */
  void theta2velcomp(const vector<real_prec> & delta, vector<real_prec> &vei, bool zeropad, bool norm, int comp);

  void comp_velbias(const vector<real_prec> &delta, vector<real_prec>&out, bool zeropad, bool norm);


  void normalize_df_z_ini(vector<real_prec>&, vector<real_prec>&, string type, gsl_rng * gBaseRand);
  //////////////////////////////////////////////////////////
  /**
   * @brief Object of class Cosmology
   */
  Cosmology Cosmo;
  //////////////////////////////////////////////////////////
  /**
   * @brief Structure allocating cosmological parameters
   */
  s_CosmologicalParameters s_cosmo_pars;
  
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   * @brief Structure allocating cosmological information
   */
  s_CosmoInfo s_cosmo_info;

  /**
   * @brief Name of file containing the DM produced by PATCHY
   */
  string fnamePOSX;
  string fnameVXpart;
  /**
   * @brief Name of file containing the DM produced by PATCHY
   */
  string fnamePOSY;
  string fnameVYpart;

  /**
   * @brief Name of file containing the DM produced by PATCHY
   */
  string fnamePOSZ;
  string fnameVZpart;
  
  /**
   * @brief Name of file containing the DM produced by PATCHY
   */
  string fnameVX;
  /**
   * @brief Name of file containing the DM produced by PATCHY
   */
  string fnameVY;
  /**
   * @brief Name of file containing the DM produced by PATCHY
   */
  string fnameVZ;

  //////////////////////////////////////////////////////////
  /**
   * @brief Name of file containing the Fourier grid with the 3D power spectrum therein interpolated. No ".dat"
   */
   string fname3DPOWER;

   //////////////////////////////////////////////////////////  
  /**
   * @brief Name of file the White noise in conf space
   */
   string fnameIC;
   
   //////////////////////////////////////////////////////////  
   /**
    * @brief Name of file the initial density field displaying an initial power spectrum
    */
   string fnameICDELTA;

      //////////////////////////////////////////////////////////
   string fnameDM;
   
   
   //////////////////////////////////////////////////////////
   string fname2LPTTERM;


   //////////////////////////////////////////////////////////
   string fnameTHETA;


   //////////////////////////////////////////////////////////
   string fnameDMNGP;

   string fnameS2TERM;

   string fnameS2TERMEUL;
   string fnameP2TERMEUL;

   string fnameS3TERM;
   string fnameS3TERMEUL;
   string fnameSTTERM;
   string fnameSTTERMEUL;
   string fnamePSITERMEUL;
   string fname2LPTTERMEUL;
   string stradd;

   string fnameTRACERCAT;
   string fname_MOCK_NCOUNTS;
   string fname_MOCK_MASS;
   string fname_MOCK_NCOUNTS_SAT;
   
   ULONG Number_of_Tracers;

   int sfmodel;
};
  

#endif
  
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################

