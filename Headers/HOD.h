#ifndef _HOD_
#define _HOD_


class HOD{
 public:


  HOD(){}
  ~HOD(){}
  /**
   * Mass of central galaxy as a function of the halo mass.
   */
  double CENTRAL(double, s_CosmologicalParameters *);
  
  
  /**
   *  HOD for satellite galaxies as a function of the halo mass
   */
  double SATELLITE(double, s_CosmologicalParameters *);
};

#endif
