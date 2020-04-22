#ifndef FFTW_ARRAY
#define FFTW_ARRAY

#include "fftw3.h"
//#include<iostream>
/* This template allows to avoid writing fftw_free */

template<typename T> class fftw_array
{
 public:
  T *data;
  
  fftw_array(long size)
    {
      
#ifdef SINGLE_PREC
      data = (T *) fftwf_malloc(size*sizeof(T));
#endif
#ifdef DOUBLE_PREC
      data = (T *) fftw_malloc(size*sizeof(T));
#endif
    }
  
  
  ~fftw_array()
    {
#ifdef SINGLE_PREC
      fftwf_free(data);
#endif
#ifdef DOUBLE_PREC
      fftw_free(data);
#endif
    }
  /* implicit conversion: works for functions, but not for templates. For the latter case write: <array_name>.data */
  operator T*()
  {
    return data;
  }
  
};






#endif
