#include <vector>

void cellbound(ULONG N1,ULONG N2,ULONG N3,real_prec L1,real_prec L2,real_prec L3, vector<real_prec>&v1, vector<real_prec>&v2, vector<real_prec>&v3);

void cellboundcomp(ULONG N1,ULONG N2,ULONG N3,vector<real_prec>vi);

void cellbound_comp(ULONG N1,ULONG N2,ULONG N3,vector<real_prec>&,vector<real_prec>&);


void getDensity_NGPv(s_params_box_mas *params, const vector<s_Halo>&Halo, vector<real_prec>&delta, string weightmass);
  
void getDensity_NGP(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3, const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec>&zp, const vector<real_prec>&Par_mass, vector<real_prec>&delta, bool weightmass);

void getDensity_TETCIC(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec>&zp, const vector<real_prec> &Par_mass, ULONG N_OBJ, vector<real_prec>&delta,bool weightmass);

void getDensity_TETCIC_comp(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec> &yp, const vector<real_prec>&zp, const vector<real_prec> &Par_mass, ULONG N_OBJ, vector<real_prec>&delta,bool weightmass);

void getDensity_CIC(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>& xp, const vector<real_prec>&yp, const vector<real_prec>&zp, const vector<real_prec> &Par_mass, vector<real_prec>&delta, bool weightmass);


void getDensity_CICv(s_params_box_mas *params, const vector<s_Halo>& Halo, vector<real_prec>&delta, string weightmass);


void getDensity_CICSHIFT(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3, const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec>&zp, const vector<real_prec>&Par_mass, ULONG N_OBJ, vector<real_prec>&delta, bool weightmass);

void getDensity_TSC(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3, const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec>& zp, const vector<real_prec> &Par_mass, vector<real_prec>&delta, bool weightmass);


void getDensity_TSCv(s_params_box_mas *params, const vector<s_Halo>& Halo, vector<real_prec>&delta, string weightmass);


void MAS_NEW(s_params_box_mas *params, const vector<s_Halo>&, string weight, vector<real_prec> &delta);
