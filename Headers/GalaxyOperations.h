#ifndef __GALAXY_OPERATIONS__
#define __GALAXY_OPERATIONS__

# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include <vector>
# include "CosmologicalFunctions.h"




using namespace std;

class GALAXY_OPERATIONS{
 public:
  GALAXY_OPERATIONS(){}
  ~GALAXY_OPERATIONS(){}
  void compute_box_size(void *, void *, double&, double&, double&, double&, double&, double&, double &RA_los, double &Dec_los);
  void equatorial_to_cartesian(double ra, double dec, double r, double &x, double &y, double &z);
  void equatorial_to_galactic(double alpha, double delta, double &b, double &l);
  void equatorial_to_equatorial(double old_alpha, double old_delta, double m_alpha, double m_delta, double &new_alpha, double &new_delta);

  void new_equatorial_to_cartesian(double, double, double, double, double, double *, double *, double*);


};
#endif
