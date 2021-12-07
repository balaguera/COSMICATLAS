/**
 *  @file CoordinateSystem.h
 *
 *  @brief The class CoordinateSystem 
 *
 *  This file defines the interface of the class CoordinateSystem
 */
#ifndef __COORDINATE_SYSTEM__
#define __COORDINATE_SYSTEM__

# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include <vector>

using namespace std;
# include "CosmologicalFunctions.h"



/*!\class GalaxyOperations

  @details   Estimates of power spectrum
  @author    Andres Balaguera Antolinez 
  @author    Optimization and parallelization by Luca Tornatore
  @version   1.1a
  @date      2013-2015
*/

/**
 *  @struct s_CoordinateSystem
 *
 *  @brief Structure containing the data used to compute coordinate
 *  trasformations
 */
struct s_galaxy_operationsF{
  int sys_coord;
  int i_coord1;
  int i_coord2;
  int i_coord3;
  int n_columns;
  string angles_units;
  vector<real_prec>zz;
  vector<real_prec>rc;
  vector< real_prec > properties;
};



/**
 *  @class CoordinateSystem CoordinateSystem.h "Lib/Headers/CoordinateSystem.h"
 *
 *  @brief The class CoordinateSystem
 *
 *  @author Andres Balaguera Antolinez 
 *
 *  @version 1.1a
 *
 *  @date 2013-2015
 *
 * This class is used to handle objects of type <EM> CoordinateSystem
 * </EM>
 */
class CoordinateSystem{




 public:



  /**
   *  @brief default constructor
   *  @return object of class CoordinateSystem
   */
  CoordinateSystem(){}


  /**
   *  @brief default destructor
   *  @return none
   */
  ~CoordinateSystem(){}


  /** Minumum x-coordinate in the sample */
  real_prec XMIN;

  /** Minumum y-coordinate in the sample */
  real_prec YMIN;

  /** Minumum z-coordinate in the sample */
  real_prec ZMIN;


  ///////////////////////////////////////////////////////
  
  /**
   * @brief Get the size of box in configuration space
   * @details Given the galaxy sample, compute the minimum and maximum values 
   * of the cartesian coordinates
   * in order to generate the size of the box (or cuboid) for the DFT
   */
  void compute_box_size(bool, void *, real_prec*);
  ///////////////////////////////////////////////////////

 /**
   * @brief Convert equatorial coordinates, i.e, right ascension, declination and distance  
   * to cartessian coordinates. 
   * @details Input values of RA and Dec are in degrees,          
   * the output is in the units of r       
   */

  void equatorial_to_cartesian(real_prec, real_prec, real_prec, real_prec &, real_prec &, real_prec &);
  
  ///////////////////////////////////////////////////////
   /**
   * @brief Convert equatorial coordinates, i.e, right ascension, declination and distance  
   * to a new spherical coordinate system 
   * @details with the z-axis points towards          
   * the point identified with m_ra and m_dec, that is,  the new e z-direction          
   * points towards the sample baricenter.
   * Input values of RA and Dec are in degrees,                                      
   * the output has the units of r         
  */
  void new_equatorial_to_cartesian(real_prec, real_prec, real_prec , real_prec, real_prec, real_prec &, real_prec &, real_prec &);
  ///////////////////////////////////////////////////////
  /**
   * \brief Convert equatorial coordinates, i.e, right ascension, declination and distance  
   * to galactic coordinates . 
   * \details Input values of RA and Dec are in degrees,           
   * the output are in radians      
   */
  void equatorial_to_galactic(real_prec, real_prec, real_prec *, real_prec *);
  ///////////////////////////////////////////////////////
  /**
   * @brief Convert equatorial coordinates, i.e, right ascension and declination           
   * to a new equatorial coordinate 
   * @details in which the axis north pole                    
   * has angular coordinates  real_prec m_ra, real_prec m_dec.
   * Input values of RA and Dec in Degrees,                                         
   * output in Degrees          
   */
  void equatorial_to_equatorial(real_prec, real_prec, real_prec, real_prec, real_prec *, real_prec *);
  
};
#endif
