# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include "../Headers/Constants.h"
# include "../Headers/ScalingRelations.h"
# include "../Headers/Marks.h"

double MARKS::CENTRAL(double M, void *p){
  M/=1e14;
  SCALING_RELATIONS sc; /*las relaciones de escala esperan masas en unidades de 10^14*/
  return 1.0; 
}


double MARKS::SATELLITE(double M, void *p){
  return 1.0;
}

double MARKS::HALO(double M, void *p){
  M/=1e14;
  SCALING_RELATIONS sc;
  //  return exp(sc.MOCKS_M2L(m,p)+0.5*sigma_red*sigma_red ) ;  /*marking with the mean luminosity <L|m>*/
  /*this takes into account that ln<L|M> = <ln L|M> + 0.5 sigma^{2}*/   
  return 1.0;
}

