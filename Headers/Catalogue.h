/**
 *  @file Catalogue.h
 *
 *  @brief The class Catalogue 
 *  @author Optimization and parallelization by Luca Tornatore
 *
 *  This file defines the interface of the class Catalogue, which is
 *  essentially a vector of galaxies and methods associated with it
 *  (pair counts, correlation function computations, etc.)
 */


#ifndef __CATALOGUE__
#define __CATALOGUE__

#include "GalaxyF.h"

/**
 *  @class Catalogue Catalogue.h "Lib/Headers/Catalogue.h"
 *
 *  @brief The class Catalogue
 *
 *  This class is used to handle objects of type <EM> Catalogue </EM>
 */
class Catalogue {

  private :
    
  /// vector containing the objects of the catalogue
  vector<Galaxy> object;
  
 
 public :

  /**
   *  @brief default constructor
   *  @return object of class Catalogue
   */
  Catalogue () {};

  /**
   *  @brief constructor 
   *  @param nObj number of objects in the catalogue
   *  @return object of class Catalogue
   */
  Catalogue (int nObj) { 
    object.erase(object.begin(), object.end());
    Galaxy galaxy;
    for (int i=0; i<nObj; i++) object.push_back(galaxy);
  };
  
  /**
   *  @brief default destructor
   *  @return none
   */
  ~Catalogue () {};

  /**
   *  @brief add a galaxy into the catalogue
   *  @param galaxy object of class Galaxy
   *  @return none
   */
  void add_galaxy (Galaxy galaxy) { object.push_back(galaxy); };

  /**
   *  @brief add a galaxy into the catalogue using [] operator
   *  @param galaxy object of class Galaxy
   *  @return none
   */
  void add_galaxy (Galaxy galaxy, int i) {object[i] = galaxy;};

  /**
   *  @brief reserve room in memory for the entire catalogue
   *  @param unsigned integer (the number of objects in the catalogue)
   *  @return none
   */
  void resize_catalogue (unsigned n) {object.resize(n);};
  
  /**
   *  @brief read catalogue from a file
   *  @param params object of class Parameters
   *  @param type 0 &rarr; data catalogue; 1 &rarr; random catalogue
   *  @param icatalogue index of the catalogue (for covariances)
   *  @return none
   */
  void read_catalogue (Params &, bool, int icatalogue=1);
  
  /**
   *  @brief number of objects in the catalogue
   *  @return the number of objects
   */
  int number_of_objects () { return object.size(); };
  
  /**
   *  @brief total weight of objects in the catalogue
      @return total weight
   */
  double total_weight ()  
  { 
    double NW = 0.;
    for (size_t i=0; i<object.size(); i++) NW += object[i]._weight();
    return NW;
  };

  /**
   *  @brief get the x coordinate of the i-th object of the catalogue
   *  @return object[i]._x()
   */
  double _x (int i) { return object[i]._x(); };

  /**
   *  @brief get the y coordinate of the i-th object of the catalogue
   *  @return object[i]._y()
   */
  double _y (int i) { return object[i]._y(); };

  /**
   *  @brief get the z coordinate of the i-th object of the catalogue
   *  @return object[i]._z()
   */
  double _z (int i) { return object[i]._z(); };

  /**
   *  @brief get the comoving distance of the i-th object of the catalogue
   *  @return object[i]._r()
   */
  double _r (int i)  { return object[i]._r(); };

  /**
   *  @brief get the ra coordinate of the i-th object of the catalogue
   *  @return object[i]._ra()
   */
  double _ra (int i) { return object[i]._ra(); };

  /**
   *  @brief get the dec coordinate of the i-th object of the catalogue
   *  @return object[i]._dec()
   */
  double _dec (int i) { return object[i]._dec(); };

  /**
   *  @brief get the redshift of the i-th object of the catalogue
   *  @return object[i]._redshift()
   */
  double _redshift (int i) { return object[i]._redshift(); };

  /**
   *  @brief get the weight of the i-th object of the catalogue
   *  @return object[i]._weight()
   */
  double _weight (int i) { return object[i]._weight(); };

};


#endif
