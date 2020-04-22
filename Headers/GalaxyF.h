/**
 *  @file Galaxy.h
 *
 *  @brief The class Galaxy 
 *
 *  This file defines the interface of the class Galaxy.  Note: Our
 *  code may in the future be applied to all kinds of objects
 *  (e.g. clusters, AGN etc.) and we may want to add extra properties
 *  (e.g.Halo mass, Halpha-flux, etc.)
 */

#ifndef __GALAXY__
#define __GALAXY__

#include "Params.h"
#include "CoordinateSystem.h"
#include "CosmologicalFunctions.h"

/**
 *  @class Galaxy Galaxy.h "Lib/Headers/Galaxy.h"
 *
 *  @brief The class Galaxy
 *
 *  This class is used to handle objects of type <EM> Galaxy </EM>
 */
class Galaxy {

 private :

  /// Cartesian coordinate x
  double x;
  
  /// Cartesian coordinate y
  double y;
  
  /// Cartesian coordinate z
  double z;

  /// distance from the observer
  double r;

  /// angular coordinate RA
  double ra;
  
  /// angular coordinate Dec
  double dec;

  /// redshift
  double redshift;

  /// weight
  double weight;

  
 public :

  /**
   *  @brief default constructor
   *  @return object of class Galaxy
   */
  Galaxy () {};

  /**
   *  @brief constructor 
   *  @param _x Cartesian coordinate x
   *  @param _y Cartesian coordinate y
   *  @param _z Cartesian coordinate z
   *  @return object of class Galaxy
   */
  Galaxy (double _x, double _y, double _z) 
    { x = _x; y = _y; z = _z; weight = 1.; r = sqrt(x*x+y*y+z*z); };

  /**
   *  @brief constructor 
   *  @param _x Cartesian coordinate x
   *  @param _y Cartesian coordinate y
   *  @param _z Cartesian coordinate z
   *  @param _weight weight
   *  @return object of class Galaxy
   */
  Galaxy (double _x, double _y, double _z, double _weight) 
    { x = _x; y = _y; z = _z; weight = _weight; r = sqrt(x*x+y*y+z*z); };

  /**
   *  @brief constructor 
   *  @param coord1 Cartesian coordinate x or RA
   *  @param coord2 Cartesian coordinate y or Dec
   *  @param coord3 Cartesian coordinate z or redshift or r
   *  @param _weight weight
   *  @param syst_coord 0 &rarr; [x, y, z]; 1 &rarr; [ra, dec,
   *  redshift]; 2 &rarr; [ra, dec, r]
   *  @param cosmo_par pointer to a structure s_CosmologicalParameters
   *  @return object of class Galaxy
   */
  Galaxy (double coord1, double coord2, double coord3, double _weight, int syst_coord, s_CosmologicalParameters *cosmo_par)
    { 
      if (syst_coord==0) {
	x = coord1; 
	y = coord2; 
	z = coord3; 
	weight = _weight; 
	r = sqrt(x*x+y*y+z*z);
      } else if (syst_coord==1) {
	ra = coord1;
	dec = coord2; 
	redshift = coord3;
	weight = _weight;
	Cosmology c_Cf;
	r = c_Cf.comoving_distance(redshift, cosmo_par);
	CoordinateSystem cs;
	cs.equatorial_to_cartesian(ra, dec, r, &x, &y, &z);
      } else if (syst_coord==2) {
	ra = coord1;
	dec = coord2; 
	redshift = r = coord3;
	weight = _weight;
	CoordinateSystem cs;
	cs.equatorial_to_cartesian(ra, dec, r, &x, &y, &z);
      }
    };

  /**
   *  @brief constructor 
   *  @param _x Cartesian coordinate x
   *  @param _y Cartesian coordinate y
   *  @param _z Cartesian coordinate z
   *  @param _ra RA
   *  @param _dec Dec
   *  @param _redshift redshift
   *  @param _weight weight
   *  @return object of class Galaxy
   */
  Galaxy (double _x, double _y, double _z, double _ra, double _dec, double _redshift, double _weight) 
    { x = _x; y = _y; z = _z; ra = _ra; dec = _dec; redshift = _redshift; weight = _weight; r = sqrt(x*x+y*y+z*z); };

  /**
   *  @brief default destructor
   *  @return none
   */
  ~Galaxy () {};

  /**
   *  @brief get the value of the private member Galaxy::x
   *  @return Galaxy::x
   */
  double _x () { return x; };

  /**
   *  @brief get the value of the private member Galaxy::y
   *  @return Galaxy::y
   */
  double _y () { return y; };

  /**
   *  @brief get the value of the private member Galaxy::z
   *  @return Galaxy::z
   */
  double _z () { return z; };

  /**
   *  @brief get the value of the private member Galaxy::r
   *  @return Galaxy::r
   */
  double _r () { return r; };

  /**
   *  @brief get the value of the private member Galaxy::ra
   *  @return Galaxy::ra
   */
  double _ra () { return ra; };

  /**
   *  @brief get the value of the private member Galaxy::dec
   *  @return Galaxy::dec
   */
  double _dec () { return dec; };

  /**
   *  @brief get the value of the private member Galaxy::redshift
   *  @return Galaxy::redshift
   */
  double _redshift () { return redshift; };

  /**
   *  @brief get the value of the private member Galaxy::weight
   *  @return Galaxy::weight
   */
  double _weight () { return weight; };

};


#endif
