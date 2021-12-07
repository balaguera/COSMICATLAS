// ***************************************************************************************************************
// ***************************************************************************************************************

/**
 * @class<FileOutput>
 * @brief Header file for the class FileOutput::
 * @file FileOutput.h
 * @title Methods for input/output in BAM
 * @author Andres Balaguera-Antolínez
 * @version   1.0
 * @date      2020
*/

// ***************************************************************************************************************
// ***************************************************************************************************************



#ifndef __FILE_OUTPUT__
#define __FILE_OUTPUT__

# include <ctime>
# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include <iomanip>
# include <sstream>
# include <cassert>
# include <vector>
# include "fftw_array.h" 
# include "ScreenOutput.h"
# include <bstream.h> 


using namespace std;


class FileOutput{
  
 private:
  
  bifstream inStream;

  
  
  
 public:
  
  //////////////////////////////////////////////////////////
  /** 
   * @brief Constructor
   */
  FileOutput(){

#ifdef SINGLE_PREC
    this->out.precision(6);
#else
    this->out.precision(10);
#endif
    this->out.setf(ios::showpoint);
    this->out.setf(ios::scientific);
    this->out.width(12);
    
  }
  //////////////////////////////////////////////////////////
  /**
   * @brief  Destructor
   */
  ~FileOutput(){}
  //////////////////////////////////////////////////////////
  /**
   * @brief  Destructor
   */
  string input_type;
  
  
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  ScreenOutput So;
  
  //////////////////////////////////////////////////////////
  /**
   * @brief Reads input file and allocate it in a vector
   * @param ifile input file (ascii)
   */
  int read_file(string ifile,vector<real_prec > &prop);

  //////////////////////////////////////////////////////////
  /**
   * @brief 
   * @param 
   */
  ifstream in;
  
  //##################################################################################
  //##################################################################################
  //##################################################################################
  /**
   * @brief 
   * @param 
   */
  ofstream out;
  
  //##################################################################################
  //##################################################################################
  //##################################################################################
  /** 
   * @brief Reads input file and allocate it in a matrix
   * @param ifile input file (ascii)
   * @result prop 2d container with the number of lines equal to the number of objects in the
   * catalogue and the number of rows equal to the number of columns in the catalogue
   */
  void read_file(string ifile, vector< vector< real_prec > > &prop);
  //  void read_file(string ifile, vector< double > &prop);
  
  //##################################################################################

  bool exists(string);
  //##################################################################################
  //##################################################################################
  /**
   * @brief  
   */
  int read_file_one(string, vector< real_prec > &); // For files with one column
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
//  int read_file_one(string, vector< float > &); // For files with one column
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  void read_file_one_N(long N, string, vector< real_prec > &); // For files with one column

  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  ULONG read_file(string fname, vector<real_prec> &prop_ob, int NTHREADS);


  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  ULONG read_binary_file(string fname, vector<real_prec> &prop_ob);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void write_to_file(string,  vector<real_prec>&);
  void write_to_file_i(string, vector<ULONG>&);
  void write_to_file_i(string,  vector<int>&);

  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void write_to_file(string, vector<real_prec>& , vector<real_prec>&);

  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<real_prec>&, vector<real_prec>&,vector<real_prec> &);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<real_prec>&, vector<real_prec>&,vector<real_prec> &,vector<real_prec> &);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<real_prec>&, vector<real_prec>&, vector<real_prec> &,vector<real_prec> &,vector<real_prec> &);

  
  //////////////////////////////////////////////////////////
  /**
   *@brief Write output
   */
  void write_to_file(string fname, vector<real_prec>& kve, vector<real_prec> &pk, vector<real_prec> &pkk,vector<real_prec> &pkkk, vector<real_prec> &si, vector<int>& nmod);
  
  //////////////////////////////////////////////////////////
  /**
   *@brief Write output
   */
  void write_to_file(string fname, vector<real_prec>& kve, vector<real_prec> &pk, vector<real_prec> &pkk,vector<real_prec> &pkkk, vector<int>& nmod);
  
  //////////////////////////////////////////////////////////
  /** 
   *@brief  Write output
   */
  void write_to_file(string fname, vector<real_prec>& kve, vector<real_prec> &pk, vector<real_prec> &pkk,vector<real_prec> &pkkk, vector<real_prec> &si, vector<real_prec>& nmod);
  
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<real_prec>&, vector< real_prec >&, int);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<real_prec> &, vector<real_prec>&, vector< vector<real_prec> >&);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<real_prec> &, vector<int>&, vector< vector<real_prec> >&);
  //////////////////////////////////////////////////////////

  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<int> &, vector<int>&, vector< vector<real_prec> >&);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<int> &, vector<real_prec>&, vector<real_prec>&, vector<real_prec> &);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<real_prec>&, vector<real_prec> &, vector<int> &);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file2(string, vector<real_prec>&, vector<real_prec> &, vector<int> &);
  void write_to_file2(string, vector<real_prec>&, vector<real_prec> &, vector<int> &, bool);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string fname, vector<real_prec> &kve, vector< real_prec> &bis, vector< real_prec> &sn_bis, vector<int> &mod);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void read_array(string, float *, ULONG);

  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void read_array(string, vector<real_prec>&);
  void read_array(string, vector<ULONG>&);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_array(string, float *, ULONG);
  //////////////////////////////////////////////////////////
  /**
   *@brief   Write output
   */
  template<class Type> void read_array(string, vector<real_prec>&);

  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_array(string fname, vector<real_prec>&out);
  //////////////////////////////////////////////////////////
  /**
   *@brief   Write output
   */
  void write_array(string fname, vector<ULONG>&out);
  void write_array(string fname, vector<int>&out);
  //////////////////////////////////////////////////////////
  /**
   *@brief   Write output
   */

  ULONG count_lines(string);
  //////////////////////////////////////////////////////////
template<class Type> void read_array_t(string fname, vector<real_prec>&OUT){

      ULONG N=OUT.size();
      fftw_array<Type> dummy(N);
#ifdef _FULL_VERBOSE_
      this->So.message_screen("Reading binary file", fname);
#endif
      this->inStream.open(fname.data(),file_is_natural);
      if(!this->inStream.is_open())
        {
          cout<<RED<<"File not found!."<<RESET<<endl;
          cout<<RED<<"Check input parameter file."<<RESET<<endl;
          cout<<RED<<"Code exits here."<<RESET<<endl;
          cout<<endl;
          exit(0);
        }
      assert(this->inStream.is_open());
      this->inStream.get(dummy.data,N);
      this->inStream.close();
    #pragma omp parallel for
      for(ULONG i=0;i<N;i++)
        OUT[i]=static_cast<real_prec>(dummy[i]);
      So.DONE();
     return ;
  }

//////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
template<class Type> void read_array_t(string fname, vector<Type>&OUT){


#ifdef _USE_OMP_
     int NTHREADS=_NTHREADS_;
 omp_set_num_threads(NTHREADS);
#endif

      ULONG N=OUT.size();
      fftw_array<Type> dummy(N);
#ifdef _FULL_VERBOSE_
      this->So.message_screen("Reading binary file", fname);
#endif
      this->inStream.open(fname.data(),file_is_natural);
      if(!this->inStream.is_open())
        {
          cout<<RED<<"File not found!."<<RESET<<endl;
          cout<<RED<<"Check input parameter file."<<RESET<<endl;
          cout<<RED<<"Code exits here."<<RESET<<endl;
          cout<<endl;
          exit(0);
        }
      assert(this->inStream.is_open());
      this->inStream.get(dummy.data,N);
      this->inStream.close();
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<N;i++)
        OUT[i]=static_cast<real_prec>(dummy[i]);
      So.DONE();
     return ;
  }

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
 
 
#ifdef SINGLE_PREC


 void write_to_file(string, vector<gsl_real>&);
//////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void write_to_file(string, vector<gsl_real>& , vector<gsl_real>&);

  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<gsl_real>&, vector<gsl_real>&,vector<gsl_real> &);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<gsl_real>&, vector<gsl_real>&,vector<gsl_real> &,vector<gsl_real> &);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<gsl_real>&, vector<gsl_real>&, vector<gsl_real> &,vector<gsl_real> &,vector<gsl_real> &);

  
  //////////////////////////////////////////////////////////
  /**
   *@brief Write output
   */
  void write_to_file(string fname, vector<gsl_real>& kve, vector<gsl_real> &pk, vector<gsl_real> &pkk,vector<gsl_real> &pkkk, vector<gsl_real> &si, vector<int>& nmod);
  
  //////////////////////////////////////////////////////////
  /** 
   *@brief  Write output
   */
  void write_to_file(string fname, vector<gsl_real>& kve, vector<gsl_real> &pk, vector<gsl_real> &pkk,vector<gsl_real> &pkkk, vector<gsl_real> &si, vector<gsl_real>& nmod);
  
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<gsl_real>&, vector< gsl_real >&, int);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<gsl_real> &, vector<gsl_real>&, vector< vector<gsl_real> >&);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<gsl_real> &, vector<int>&, vector< vector<gsl_real> >&);
  //////////////////////////////////////////////////////////

  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<int> &, vector<int>&, vector< vector<gsl_real> >&);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<int> &, vector<gsl_real>&, vector<gsl_real>&, vector<gsl_real> &);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<gsl_real>&, vector<gsl_real> &, vector<int> &);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file2(string, vector<gsl_real>&, vector<gsl_real> &, vector<int> &);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string fname, vector<gsl_real> &kve, vector< gsl_real> &bis, vector< gsl_real> &sn_bis, vector<int> &mod);
  //////////////////////////////////////////////////////////

  bifstream inStreamp;

#endif

 
 
};








#endif

