#ifndef _HOD_
#define _HOD_
# include "NumericalMethods.h"


class HOD{

 private:
 public:


  HOD(){}
  ~HOD(){}
  /**
   * Mass of central galaxy as a function of the halo mass.
   */
real_prec CENTRAL(real_prec, s_CosmologicalParameters *);


  /**
   *  HOD for satellite galaxies as a function of the halo mass
   */
  real_prec SATELLITE(real_prec, s_CosmologicalParameters *);
};

#endif
