#ifndef _SCALING_RELALTIONS_
#define _SCALING_RELATIONS_




class SCALING_RELATIONS{
  double m;
public:
  SCALING_RELATIONS(){};
  ~SCALING_RELATIONS(){};
  double M2T(double, void *);
  double RB_M2L(double, void *);
  double FED_M2L(double, void *);
  double STA_M2L(double, void *);
  double MANTZ_BOL_M2L(double, void *);
  double MOCKS_M2L(double, void *);
  double MANTZ_BAND_M2L(double , void *);
};

#endif
