#ifndef _ASTROPHYSICS_
#define _ASTROPHYSICS_


class ASTROPHYSICS{

 public:
  ASTROPHYSICS(){};
  ~ASTROPHYSICS(){};
  
  double mass2temp(double, double, void*);
  double filtering_mass(double, void*);
  double cooling_function(double, void*);
  double vc(double, double, void*); 
  double cooling_radius(double, double, double, void *);
  double cooling_mass_rate(double, double, void *);
  double infall_mass_rate(double, double, void *);
  double baryon_fraction(double, double, void *); 
  double star_formation_rate(double, double, double, void *);
  double baryon_mass_virial_mass_distribution(double, double, double, void *);

  double total_mass_nb_mass_relation(double, double,void *, void *);
  //  double Mobs_Mnb_distribution(double, double, double, void *);
  double cluster_mass_function(double, double, double, void *);
  double mass_lum_distribution(double, double, double, void *, void *);
  double mass_lum_distribution_errors(double, double, double, void *, void*);
  double G(double, double, double, void *, void *); 
  double G_Lmin_Lmax(double, double, double, void *, void *);  
};
#endif
