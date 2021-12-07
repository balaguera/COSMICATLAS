#include <vector>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>         
#include <string.h>
#include <cassert>
#include <cfloat>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <omp.h>
#include <sstream>
#include <netinet/in.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>

using namespace std;

//#include "../Headers/def.h"
#include "../Headers/NumericalMethods.h"
#include "../Headers/Type_structures_def.h"
#include "../Headers/FftwFunctions.h"
#include "../Headers/massFunctions.h"
#include "../Headers/ScreenOutput.h"




void getDensity_NGP(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec>&zp, const vector<real_prec> &Par_mass, vector<real_prec>&delta, bool weightmass)
{
  
  ScreenOutput So;

  real_prec xc, yc, zc;
  real_prec dx, dy, dz, tx, tz, ty;    
  ULONG i, j, k, ii, jj, kk, n;
  
  ULONG N_OBJ=xp.size();


#ifndef _GET_INTERPOLATED_FIELDS_FROM_BIN_FILES_
#pragma omp parallel for
  for (i=0;i<delta.size(); i++)
    delta[i]= 0.;  //-1 if we want to calculate overdensity
#endif
  
  ULONG NLOSS=0;
  ULONG count=0;

  for (n=0; n<N_OBJ; n++){
    
    //check if particle is in selected Domain, else discard it
    if((xp[n]>=min1 && xp[n]<min1+L1) && (yp[n]>=min2 && yp[n]<min2+L2) && (zp[n]>=min3 && zp[n]<min3+L3))
      {
	
	i = static_cast<ULONG>(floor((xp[n]-min1)/d1)); // indices of the cell of the particle
	j = static_cast<ULONG>(floor((yp[n]-min2)/d2));
	k = static_cast<ULONG>(floor((zp[n]-min3)/d3));
		
	i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
	j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
	k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));
		
	real_prec mass=num_1;
	if (weightmass==true)
	  mass=Par_mass[n];
	
	delta[k+N3*(j+N2*i)]    +=mass;
    }
   
    else
      NLOSS++;	  
    
  }


#ifndef  _GET_INTERPOLATED_FIELDS_FROM_BIN_FILES_
    count=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:count)
#endif
    for(ULONG i=0;i<delta.size();++i)
      count+=delta[i];
#ifdef _FULL_VERBOSE_
    So.message_screen("Number of objects assigned to grid =",count);
#endif
#endif

#ifdef _FULL_VERBOSE_
  if(NLOSS!=0)
    if(NLOSS!=0) So.message_screen("mass assignment found",NLOSS," particles outside mesh boundary");
#endif

}


// ===========================================================================================================================
// ===========================================================================================================================
// ===========================================================================================================================

void getDensity_NGPv(s_params_box_mas *params, const vector<s_Halo>&Halo, vector<real_prec>&delta, string weight_prop)
{
  
  ScreenOutput So;

  ULONG N_OBJ=Halo.size();

  real_prec min1=params->min1;
  real_prec min2=params->min2;
  real_prec min3=params->min3;
  real_prec L1=params->Lbox;
  real_prec L2=params->Lbox;
  real_prec L3=params->Lbox;
  real_prec N1=params->Nft;
  real_prec N2=params->Nft;
  real_prec N3=params->Nft;
  real_prec d1=params->d1;
  real_prec d2=params->d2;
  real_prec d3=params->d3;
  
  
#pragma omp parallel for
  for (ULONG i=0;i<delta.size(); i++)
    delta[i]= 0.;  //-1 if we want to calculate overdensity
  
  
  ULONG NLOSS=0;
  for (ULONG n=0; n<N_OBJ; n++)
    {
      
      //check if particle is in selected Domain, else discard it
      if((Halo[n].coord1>=min1 && Halo[n].coord1 <=min1+L1) && (Halo[n].coord2>=min2 && Halo[n].coord2<=min2+L2) && (Halo[n].coord3>=min3 && Halo[n].coord3 <=min3+L3))
	{

          int i=get_bin(Halo[n].coord1,min1,N1,d1,true);
          int j=get_bin(Halo[n].coord2,min2,N2,d2,true);
          int k=get_bin(Halo[n].coord3,min3,N3,d3,true);


	  i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
	  j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
	  k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));
	  


          ULONG index=index_3d(i,j,k,N2,N3);
	  real_prec tracer_weight=num_1;
	  if (weight_prop==_MASS_)
	    tracer_weight=Halo[n].mass;
	  else if (weight_prop==_SAT_FRACTION_)
	    tracer_weight=Halo[n].number_sub_structures;
          else if (weight_prop==_RS_)
            tracer_weight=Halo[n].rs;
          else if (weight_prop==_VMAX_)
            tracer_weight=Halo[n].vmax;
          else if (weight_prop==_VIRIAL_)
            tracer_weight=Halo[n].virial;
          else if (weight_prop==_SPIN_)
            tracer_weight=Halo[n].spin;

          delta[index]    += tracer_weight;
	  
	}
      
      else
	NLOSS++;	  
      
    }
  
  
  So.DONE();
  
  if(_COUNTS_ == weight_prop)
    {
      ULONG count=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:count)
#endif
      for(ULONG i=0;i<delta.size();++i)
	count+=delta[i];
#ifdef _FULL_VERBOSE_
      So.message_screen("Number of objects assigned to grid =",count);
      if(NLOSS!=0) So.message_screen("mass assignment found",NLOSS,"particles outside mesh boundary");
#endif
    }
  
}



// ===========================================================================================================================
// ===========================================================================================================================
// ===========================================================================================================================


void getDensity_CICSHIFT(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,vector<real_prec>xp, vector<real_prec>yp, vector<real_prec>zp, vector<real_prec> Par_mass, ULONG N_OBJ, vector<real_prec>&delta, bool weightmass)
{
    real_prec xc, yc, zc;
    real_prec dx, dy, dz, tx, tz, ty;    
    ULONG i, j, k, ii, jj, kk, n;

#define DELTA(i,j,k) delta[k+N3*(j+N2*i)]


#pragma omp parallel for
  for (i=0;i<delta.size(); i++)
    delta[i]= 0.;  //-1 if we want to calculate overdensity


    ULONG NLOSS=0;
    for (n=0; n<N_OBJ; n++){

	//check if particle is in selected Domain, else discard it
      if((xp[n]>=min1 && xp[n]<min1+L1) && (yp[n]>=min2 && yp[n]<min2+L2) && (zp[n]>=min3 && zp[n]<min3+L3))
     {	
        i = static_cast<ULONG>(floor((xp[n]-min1)/d1)); // indices of the cell of the particle
	j = static_cast<ULONG>(floor((yp[n]-min2)/d2));
	k = static_cast<ULONG>(floor((zp[n]-min3)/d3));

	//printf("%f, %f, %f\n", xp[n], yp[n], zp[n]);		

	i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
	j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
	k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));

	ii = static_cast<ULONG>(fmod(real_prec(i+1),real_prec(N1)));
	jj = static_cast<ULONG>(fmod(real_prec(j+1),real_prec(N2)));
	kk = static_cast<ULONG>(fmod(real_prec(k+1),real_prec(N3)));

	xc = static_cast<real_prec>(i); 
	yc = static_cast<real_prec>(j);
	zc = static_cast<real_prec>(k);

	dx = (xp[n]-min1)/d1 - xc; // distance of particle to center of the cell
	dy = (yp[n]-min2)/d2 - yc;
	dz = (zp[n]-min3)/d3 - zc;
	
	tx = num_1 - dx;
	ty = num_1 - dy;
	tz = num_1 - dz;

// Take care this assumes periodic boundary conditions. This conserves Mass when using FFTs to deconvolve with CIC kernel, as the FFT assumes periodicity

        real_prec mass=num_1;
        if (weightmass==true)
          mass=Par_mass[n];

        delta[k+N3*(j+N2*i)]    += mass*tx*ty*tz;
        delta[k+N3*(j+N2*ii)]   += mass*dx*ty*tz;
        delta[k+N3*(jj+N2*i)]   += mass*tx*dy*tz;
        delta[kk+N3*(j+N2*i)]   += mass*tx*ty*dz;
        delta[k+N3*(jj+N2*ii)]  += mass*dx*dy*tz;
        delta[kk+N3*(j+N2*ii)]  += mass*dx*ty*dz;
        delta[kk+N3*(jj+N2*i)]  += mass*tx*dy*dz;
        delta[kk+N3*(jj+N2*ii)] += mass*dx*dy*dz;
	
     }
      //else NLOSS++;	
    }
    if(NLOSS!=0) cout << " >>> ARGO: Mass assignment found "<<NLOSS<<" particles outside mesh boundary...."<<endl<<endl;
}


// ===========================================================================================================================
// ===========================================================================================================================
// ===========================================================================================================================

inline void get_tetrahedron_centroid( const vector<real_prec> allpos, const vector<int> &vertids, real_prec boxsize, real_prec boxhalf, vector<float>&xc)
{
  const real_prec x0[3] = { allpos[3*vertids[0]+0], allpos[3*vertids[0]+1], allpos[3*vertids[0]+2] };
	
  xc[0] = 0.0f;
  xc[1] = 0.0f;
  xc[2] = 0.0f;
	
  int i;
  
  for( i=1; i<4; ++i )
    {
      real_prec dx[3] = { allpos[3*vertids[i]+0]-x0[0], allpos[3*vertids[i]+1]-x0[1], allpos[3*vertids[i]+2]-x0[2] };
      
		
      if( dx[0] < -boxhalf ) dx[0] += boxsize;
      else if( dx[0] > boxhalf ) dx[0] -= boxsize;
      if( dx[1] < -boxhalf ) dx[1] += boxsize;
      else if( dx[1] > boxhalf ) dx[1] -= boxsize;
      if( dx[2] < -boxhalf ) dx[2] += boxsize;
      else if( dx[2] > boxhalf ) dx[2] -= boxsize;
      
      xc[0] += dx[0];
      xc[1] += dx[1];
      xc[2] += dx[2];
		
    }
	
  xc[0] *= static_cast<real_prec>(0.25);
  xc[1] *= static_cast<real_prec>(0.25);
  xc[2] *= static_cast<real_prec>(0.25);
  
  xc[0] = fmodf( xc[0]+x0[0]+boxsize, boxsize );
  xc[1] = fmodf( xc[1]+x0[1]+boxsize, boxsize );
  xc[2] = fmodf( xc[2]+x0[2]+boxsize, boxsize );
}


inline void get_tetrahedron_centroid_comp( real_prec d1, real_prec d2, real_prec d3, const vector<real_prec>&posx, const vector<real_prec>&posy, const vector<real_prec>&posz, const vector<int>&vertids, real_prec boxsize, real_prec boxhalf, vector<float>&xc )
{

   vector<real_prec> x0(3,0);

   x0[0]= posx[vertids[0]]-num_0_5*d1;
   x0[1]= posy[vertids[0]]-num_0_5*d2;
   x0[2]= posz[vertids[0]]-num_0_5*d3;
	
   xc[0] = 0.0;
   xc[1] = 0.0;
   xc[2] = 0.0;
	
  int i;
  
  for( i=1; i<4; ++i )
    {

       vector<real_prec> dx(3,0);
       dx[0]=posx[vertids[i]]-x0[0]-num_0_5*d1;
       dx[1]=posy[vertids[i]]-x0[1]-num_0_5*d2;
       dx[2]=posz[vertids[i]]-x0[2]-num_0_5*d3;
      
		
      if( dx[0] < -boxhalf ) dx[0] += boxsize;
      else if( dx[0] > boxhalf ) dx[0] -= boxsize;
      if( dx[1] < -boxhalf ) dx[1] += boxsize;
      else if( dx[1] > boxhalf ) dx[1] -= boxsize;
      if( dx[2] < -boxhalf ) dx[2] += boxsize;
      else if( dx[2] > boxhalf ) dx[2] -= boxsize;
      
      xc[0] += dx[0];
      xc[1] += dx[1];
      xc[2] += dx[2];
		
    }
	
  xc[0] *= static_cast<real_prec>(0.25);
  xc[1] *= static_cast<real_prec>(0.25);
  xc[2] *= static_cast<real_prec>(0.25);
  
  xc[0] = fmodf( xc[0]+x0[0]+boxsize, boxsize );
  xc[1] = fmodf( xc[1]+x0[1]+boxsize, boxsize );
  xc[2] = fmodf( xc[2]+x0[2]+boxsize, boxsize );
}


// ===========================================================================================================================
// ===========================================================================================================================
// ===========================================================================================================================


void getDensity_TETCIC(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec>& yp,const vector<real_prec>&zp, const vector<real_prec>&Par_mass, ULONG N_OBJ, vector<real_prec> &delta, bool weightmass)
{
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);

  static real_prec to_slab_fac=static_cast<real_prec>(N1)/L1;

  static int slabstart_x, nslab_x, slabstart_y, nslab_y, smallest_slab;

  real_prec xc, yc, zc, xpos, ypos, zpos;
  real_prec dx, dy, dz, tx, tz, ty;    
  ULONG i, j, k, ii, jj, kk, n;

// connectivity for cubic grid, six tetrahedron decomposition
  real_prec boxsize, boxhalf;
	
  boxsize = L1;
  boxhalf = num_0_5 * L1;
    
  const int vert[8][3] = { {0,0,0}, {0,0,1}, {0,1,0}, {0,1,1}, {1,0,0}, {1,0,1}, {1,1,0}, {1,1,1} };
  const int conn[][4] = { {1,0,2,4}, {3,1,2,4}, {3,5,1,4}, {3,6,5,4}, {3,2,6,4}, {3,7,5,6} };
  const int nbase = static_cast<int>(pow(static_cast<real_prec>(N_OBJ),1.0/3.0)+0.1);


  ULONG ix, iy, iz;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
	
  real_prec pweight = static_cast<real_prec>(1.0/6.0);

  vector<real_prec> allpos(3*N_OBJ,0);


#pragma omp parallel for
  for (i=0;i<delta.size(); i++)
    delta[i]= 0.;  //-1 if we want to calculate overdensity

  
#pragma omp parallel for
  for (i=0;i<N_OBJ; i++)
    {
      allpos[3*i+0]=xp[i]-num_0_5*d1;
      allpos[3*i+1]=yp[i]-num_0_5*d2;
      allpos[3*i+2]=zp[i]-num_0_5*d3;
    }
  
  for( ix=0;ix<nbase; ++ix )
    for( iy=0; iy<nbase; ++iy )
      for( iz=0; iz<nbase; ++iz )
				{
  		   int l;
         static vector<int>cube_vertices(8,0);
	       for( l=0; l<8; ++l )
		       cube_vertices[l] = (((ix+vert[l][0])%nbase)*nbase + (iy+vert[l][1])%nbase)*nbase + (iz+vert[l][2])%nbase;
	  
			  for( l=0; l<6; ++l )
	   		{
         static vector<int> vertids(4,0);
	       int m;
	       for( m=0; m<4; ++m )
		vertids[m] = cube_vertices[conn[l][m]];
	      
              static vector<float> xcc(3,0);
	      get_tetrahedron_centroid( allpos, vertids, boxsize, boxhalf, xcc );      
	      
	      xpos=xcc[0];
	      ypos=xcc[1];
	      zpos=xcc[2];
	      
	      if((xpos>=min1 && xpos<min1+L1) && (ypos>=min2 && ypos<min2+L2) && (zpos>=min3 && zpos<min3+L3))
				{	
				  i = static_cast<ULONG>(floor((xpos-min1)/d1)); // indices of the cell of the particle
	      	j = static_cast<ULONG>(floor((ypos-min2)/d2));
	      	k = static_cast<ULONG>(floor((zpos-min3)/d3));

	//printf("%f, %f, %f\n", xpos, ypos, zpos);		
	      
		      i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
		      j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
		      k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));
	      
		      ii = static_cast<ULONG>(fmod(real_prec(i+1),real_prec(N1)));
	 	  	  jj = static_cast<ULONG>(fmod(real_prec(j+1),real_prec(N2)));
	 		    kk = static_cast<ULONG>(fmod(real_prec(k+1),real_prec(N3)));

	      	xc = static_cast<real_prec>(i); 
	      	yc = static_cast<real_prec>(j);
	      	zc = static_cast<real_prec>(k);
	      
	      	dx = (xpos-min1)/d1 - xc; // distance of particle to center of the cell
	      	dy = (ypos-min2)/d2 - yc;
	      	dz = (zpos-min3)/d3 - zc;
	
		      tx = num_1 - dx;
	 	     ty = num_1 - dy;
	 	     tz = num_1 - dz;

// Take care this assumes periodic boundary conditions. This conserves Mass when using FFTs to deconvolve with CIC kernel, as the FFT assumes periodicity

	       real_prec mass=pweight;

	      if (weightmass==true)
		mass*=Par_mass[iz+N1*(iy+N1*ix)];
  
	      delta[k+N3*(j+N2*i)]    += mass*tx*ty*tz;
	      delta[k+N3*(j+N2*ii)]   += mass*dx*ty*tz;
	      delta[k+N3*(jj+N2*i)]   += mass*tx*dy*tz;
	      delta[kk+N3*(j+N2*i)]   += mass*tx*ty*dz;
	      delta[k+N3*(jj+N2*ii)]  += mass*dx*dy*tz;
	      delta[kk+N3*(j+N2*ii)]  += mass*dx*ty*dz;
	      delta[kk+N3*(jj+N2*i)]  += mass*tx*dy*dz;
	      delta[kk+N3*(jj+N2*ii)] += mass*dx*dy*dz;
	
	}
	      /* do CIC assignment of centroids */
	      
	      /*
	      slab_x = to_slab_fac * xc[0];
	      if(slab_x >= N1)
		slab_x = N1 - 1;
	      dx = to_slab_fac * xc[0] - slab_x;
	      slab_x = (slab_x+N1)%N1;
	      slab_xx = (slab_x + 1 + N1)%N1;
	      
	      slab_y = to_slab_fac * xc[1];
	      dy = to_slab_fac * xc[1] - slab_y;
	      slab_y = (slab_y+N1)%N1;
	      slab_yy = (slab_y + 1 +N1)%N1;
	      
	      slab_z = to_slab_fac * xc[2];
	      dz = to_slab_fac * xc[2] - slab_z;
	      slab_z = (slab_z+N1)%N1;
	      slab_zz = (slab_z + 1 +N1)%N1;
	      
	      if( slab_x >= 0 )
		{					
		  delta[slab_x * N1*N1 + slab_y  * N1 + slab_z]  += pweight * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
		  delta[slab_x * N1*N1 + slab_yy * N1 + slab_z]  += pweight * (1.0 - dx) * dy * (1.0 - dz);
		  delta[slab_x * N1*N1 + slab_y  * N1 + slab_zz] += pweight * (1.0 - dx) * (1.0 - dy) * dz;
		  delta[slab_x * N1*N1 + slab_yy * N1 + slab_zz] += pweight * (1.0 - dx) * dy * dz;
		}
	      
	      if( slab_xx >= 0)
		{
		  delta[slab_xx * N1*N1 + slab_y  * N1 + slab_z]  += pweight * (dx) * (1.0 - dy) * (1.0 - dz);
		  delta[slab_xx * N1*N1 + slab_yy * N1 + slab_z]  += pweight * (dx) * dy * (1.0 - dz);
		  delta[slab_xx * N1*N1 + slab_y  * N1 + slab_zz] += pweight * (dx) * (1.0 - dy) * dz;
		  delta[slab_xx * N1*N1 + slab_yy * N1 + slab_zz] += pweight * (dx) * dy * dz;
		}	      
	      */	
	    }					
	}
  //else NLOSS++;	
}



// ===========================================================================================================================
// ===========================================================================================================================
// ===========================================================================================================================

void getDensity_TETCIC_comp(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec>&zp, const vector<real_prec> &Par_mass, ULONG N_OBJ, vector<real_prec>&delta,bool weightmass)
{

  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);

  static real_prec to_slab_fac=static_cast<real_prec>(N1)/L1;

  static int slabstart_x, nslab_x, slabstart_y, nslab_y, smallest_slab;

  real_prec xc, yc, zc, xpos, ypos, zpos;
  real_prec dx, dy, dz, tx, tz, ty;    
  ULONG i, j, k, ii, jj, kk, n;

// connectivity for cubic grid, six tetrahedron decomposition
  real_prec boxsize, boxhalf;
	
  boxsize = L1;
  boxhalf = num_0_5 * L1;
    
  const int vert[8][3] = { {0,0,0}, {0,0,1}, {0,1,0}, {0,1,1}, {1,0,0}, {1,0,1}, {1,1,0}, {1,1,1} };
  const int conn[6][4] = { {1,0,2,4}, {3,1,2,4}, {3,5,1,4}, {3,6,5,4}, {3,2,6,4}, {3,7,5,6} };
  const int nbase = static_cast<int>(pow(static_cast<real_prec>(N_OBJ),1.0/3.0)+0.1);

  ULONG ix, iy, iz;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
	
  real_prec pweight = static_cast<real_prec>(1.0/6.0);


#pragma omp parallel for
  for (i=0;i<delta.size(); i++)
    delta[i]= 0.;  //-1 if we want to calculate overdensity

    
  for( ix=0;ix<nbase; ++ix )
    for( iy=0; iy<nbase; ++iy )
      for( iz=0; iz<nbase; ++iz )
	{
	  int l;
	  static int cube_vertices[8];
	  for( l=0; l<8; ++l )
	    cube_vertices[l] = (((ix+vert[l][0])%nbase)*nbase + (iy+vert[l][1])%nbase)*nbase + (iz+vert[l][2])%nbase;
	  
	  for( l=0; l<6; ++l )
	    {
              static vector<int> vertids(4,0);
	      int m;
	      for( m=0; m<4; ++m )
		vertids[m] = cube_vertices[conn[l][m]];
	      
              static vector<float> xcc(3,0);
	      get_tetrahedron_centroid_comp(d1, d2, d3, xp, yp, zp, vertids, boxsize, boxhalf, xcc );
	      

	      xpos=xcc[0];
	      ypos=xcc[1];
	      zpos=xcc[2];

      if((xpos>=min1 && xpos<min1+L1) && (ypos>=min2 && ypos<min2+L2) && (zpos>=min3 && zpos<min3+L3))
	{	
	      i = static_cast<ULONG>(floor((xpos-min1)/d1)); // indices of the cell of the particle
	      j = static_cast<ULONG>(floor((ypos-min2)/d2));
	      k = static_cast<ULONG>(floor((zpos-min3)/d3));

	//printf("%f, %f, %f\n", xpos, ypos, zpos);		
	      
	      i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
	      j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
	      k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));
	      
	      ii = static_cast<ULONG>(fmod(real_prec(i+1),real_prec(N1)));
	      jj = static_cast<ULONG>(fmod(real_prec(j+1),real_prec(N2)));
	      kk = static_cast<ULONG>(fmod(real_prec(k+1),real_prec(N3)));

	      xc = static_cast<real_prec>(i); 
	      yc = static_cast<real_prec>(j);
	      zc = static_cast<real_prec>(k);
	      
	      dx = (xpos-min1)/d1 - xc; // distance of particle to center of the cell
	      dy = (ypos-min2)/d2 - yc;
	      dz = (zpos-min3)/d3 - zc;
	
	      tx = num_1 - dx;
	      ty = num_1 - dy;
	      tz = num_1 - dz;

// Take care this assumes periodic boundary conditions. This conserves Mass when using FFTs to deconvolve with CIC kernel, as the FFT assumes periodicity

	      real_prec mass=pweight;

	      if (weightmass==true)
		mass*=Par_mass[iz+N1*(iy+N1*ix)];
  

	      delta[k+N3*(j+N2*i)]    += mass*tx*ty*tz;
	      delta[k+N3*(j+N2*ii)]   += mass*dx*ty*tz;
	      delta[k+N3*(jj+N2*i)]   += mass*tx*dy*tz;
	      delta[kk+N3*(j+N2*i)]   += mass*tx*ty*dz;
	      delta[k+N3*(jj+N2*ii)]  += mass*dx*dy*tz;
	      delta[kk+N3*(j+N2*ii)]  += mass*dx*ty*dz;
	      delta[kk+N3*(jj+N2*i)]  += mass*tx*dy*dz;
	      delta[kk+N3*(jj+N2*ii)] += mass*dx*dy*dz;

	}
	      /* do CIC assignment of centroids */
	      
	      /*
	      slab_x = to_slab_fac * xc[0];
	      if(slab_x >= N1)
		slab_x = N1 - 1;
	      dx = to_slab_fac * xc[0] - slab_x;
	      slab_x = (slab_x+N1)%N1;
	      slab_xx = (slab_x + 1 + N1)%N1;
	      
	      slab_y = to_slab_fac * xc[1];
	      dy = to_slab_fac * xc[1] - slab_y;
	      slab_y = (slab_y+N1)%N1;
	      slab_yy = (slab_y + 1 +N1)%N1;
	      
	      slab_z = to_slab_fac * xc[2];
	      dz = to_slab_fac * xc[2] - slab_z;
	      slab_z = (slab_z+N1)%N1;
	      slab_zz = (slab_z + 1 +N1)%N1;
	      
	      if( slab_x >= 0 )
		{					
		  delta[slab_x * N1*N1 + slab_y  * N1 + slab_z]  += pweight * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
		  delta[slab_x * N1*N1 + slab_yy * N1 + slab_z]  += pweight * (1.0 - dx) * dy * (1.0 - dz);
		  delta[slab_x * N1*N1 + slab_y  * N1 + slab_zz] += pweight * (1.0 - dx) * (1.0 - dy) * dz;
		  delta[slab_x * N1*N1 + slab_yy * N1 + slab_zz] += pweight * (1.0 - dx) * dy * dz;
		}
	      
	      if( slab_xx >= 0)
		{
		  delta[slab_xx * N1*N1 + slab_y  * N1 + slab_z]  += pweight * (dx) * (1.0 - dy) * (1.0 - dz);
		  delta[slab_xx * N1*N1 + slab_yy * N1 + slab_z]  += pweight * (dx) * dy * (1.0 - dz);
		  delta[slab_xx * N1*N1 + slab_y  * N1 + slab_zz] += pweight * (dx) * (1.0 - dy) * dz;
		  delta[slab_xx * N1*N1 + slab_yy * N1 + slab_zz] += pweight * (dx) * dy * dz;
		}	      
	      */	
	    }					
	}
  //else NLOSS++;	
}


// ===========================================================================================================================
// ===========================================================================================================================
// ===========================================================================================================================


void getDensity_CIC(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec> &zp, const vector<real_prec> &Par_mass, vector<real_prec>&delta, bool weightmass)
{
  real_prec xc, yc, zc, xpos, ypos, zpos;
  real_prec dx, dy, dz, tx, tz, ty;    
  ULONG i, j, k, ii, jj, kk, n;
  ULONG N_OBJ=xp.size();

  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
  
  
#ifndef  _GET_INTERPOLATED_FIELDS_FROM_BIN_FILES_
#pragma omp parallel for
  for (i=0;i<delta.size(); i++)
    delta[i]= 0.;
#endif  
  
  ULONG NLOSS=0;
  
  real_prec deltax=L1/static_cast<real_prec>(N1);
  
  for (n=0; n<N_OBJ; n++)
    {
      
      
      xpos=xp[n]-num_0_5*d1;
      ypos=yp[n]-num_0_5*d2;
      zpos=zp[n]-num_0_5*d3;
      
      real_prec xnew=xpos;
      real_prec ynew=ypos;
      real_prec znew=zpos;
      
      if (xnew<0.)
	xnew+=L1;
      if (xnew>=L1)
	xnew-=L1;
      
      if (ynew<0.)
	ynew+=L1;
      if (ynew>=L1)
	ynew-=L1;
      
      if (znew<0.)
	znew+=L1;
      if (znew>=L1)
	znew-=L1;
      
      xpos=xnew;
      ypos=ynew;
      zpos=znew;
      
      
      
      //check if particle is in selected Domain, else discard it
      if((xpos>=min1 && xpos<min1+L1) && (ypos>=min2 && ypos<min2+L2) && (zpos>=min3 && zpos<min3+L3))
	{	
	  i = static_cast<ULONG>(floor((xpos-min1)/d1)); // indices of the cell of the particle
	  j = static_cast<ULONG>(floor((ypos-min2)/d2));
	  k = static_cast<ULONG>(floor((zpos-min3)/d3));
	  
	  //printf("%f, %f, %f\n", xpos, ypos, zpos);		
	  
	  i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
	  j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
	  k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));
	  
	  ii = static_cast<ULONG>(fmod(real_prec(i+1),real_prec(N1)));
	  jj = static_cast<ULONG>(fmod(real_prec(j+1),real_prec(N2)));
	  kk = static_cast<ULONG>(fmod(real_prec(k+1),real_prec(N3)));
	  
	  xc = static_cast<real_prec>(i); 
	  yc = static_cast<real_prec>(j);
	  zc = static_cast<real_prec>(k);
	  
	  dx = (xpos-min1)/d1 - xc; // distance of particle to center of the cell
	  dy = (ypos-min2)/d2 - yc;
	  dz = (zpos-min3)/d3 - zc;
	  
	  tx = num_1 - dx;
	  ty = num_1 - dy;
	  tz = num_1 - dz;
	  
	  // Take care this assumes periodic boundary conditions. This conserves Mass when using FFTs to deconvolve with CIC kernel, as the FFT assumes periodicity
	  
	  real_prec mass=num_1;
	  if (true==weightmass)
	    mass=Par_mass[n];
	  
	  delta[k+N3*(j+N2*i)]    += mass*tx*ty*tz;
	  
	  delta[k+N3*(j+N2*ii)]   += mass*dx*ty*tz;

	  delta[k+N3*(jj+N2*i)]   += mass*tx*dy*tz;

	  delta[kk+N3*(j+N2*i)]   += mass*tx*ty*dz;

	  delta[k+N3*(jj+N2*ii)]  += mass*dx*dy*tz;

	  delta[kk+N3*(j+N2*ii)]  += mass*dx*ty*dz;

	  delta[kk+N3*(jj+N2*i)]  += mass*tx*dy*dz;

	  delta[kk+N3*(jj+N2*ii)] += mass*dx*dy*dz;
	}
      //else NLOSS++;	

}
  
  

  if(NLOSS!=0) cout << " >>> Mass assignment found "<<NLOSS<<" particles outside mesh boundary...."<<endl<<endl;
}


//============================================================================================================================
//============================================================================================================================
//============================================================================================================================

void getDensity_CICv(s_params_box_mas *params,const vector<s_Halo>&Halo, vector<real_prec>&delta, string weight_prop)
{
  real_prec xc, yc, zc, xpos, ypos, zpos;
  real_prec dx, dy, dz, tx, tz, ty;    
  ULONG i, j, k, ii, jj, kk, n;
  ULONG N_OBJ=Halo.size();
  real_prec min1=params->min1;
  real_prec min2=params->min2;
  real_prec min3=params->min3;
  real_prec L1=params->Lbox;
  real_prec L2=params->Lbox;
  real_prec L3=params->Lbox;
  real_prec N1=params->Nft;
  real_prec N2=params->Nft;
  real_prec N3=params->Nft;
  real_prec d1=params->d1;
  real_prec d2=params->d2;
  real_prec d3=params->d3;

  //#define DELTA(i,j,k) delta[k+N3*(j+N2*i)]
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);

#pragma omp parallel for
  for (i=0;i<delta.size(); i++)
    delta[i]= 0.;  //-1 if we want to calculate overdensity
  
  
  ULONG NLOSS=0;
  for (n=0; n<N_OBJ; n++)
    {
      
      xpos=Halo[n].coord1-num_0_5*d1;
      ypos=Halo[n].coord2-num_0_5*d2;
      zpos=Halo[n].coord3-num_0_5*d3;
      
      
      //if (periodic==true)
      {	 
	real_prec xnew=xpos;
	real_prec ynew=ypos;
	real_prec znew=zpos;
	
	if (xnew<0.)
	  xnew+=L1;
	if (xnew>=L1)
	  xnew-=L1;
	
	if (ynew<0.)
	  ynew+=L1;
	if (ynew>=L1)
	  ynew-=L1;
	
	if (znew<0.)
	  znew+=L1;
	if (znew>=L1)
	  znew-=L1;
	
	xpos=xnew;
	ypos=ynew;
	zpos=znew;
      }
      
      
      //check if particle is in selected Domain, else discard it
      if((xpos>=min1 && xpos<min1+L1) && (ypos>=min2 && ypos<min2+L2) && (zpos>=min3 && zpos<min3+L3))
	{	
	  i = static_cast<ULONG>(floor((xpos-min1)/d1)); // indices of the cell of the particle
	  j = static_cast<ULONG>(floor((ypos-min2)/d2));
	  k = static_cast<ULONG>(floor((zpos-min3)/d3));
	  
	  //printf("%f, %f, %f\n", xpos, ypos, zpos);		
	  
	  i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
	  j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
	  k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));
	  
	  ii = static_cast<ULONG>(fmod(real_prec(i+1),real_prec(N1)));
	  jj = static_cast<ULONG>(fmod(real_prec(j+1),real_prec(N2)));
	  kk = static_cast<ULONG>(fmod(real_prec(k+1),real_prec(N3)));
	  
	  xc = static_cast<real_prec>(i); 
	  yc = static_cast<real_prec>(j);
	  zc = static_cast<real_prec>(k);
	  
	  dx = (xpos-min1)/d1 - xc; // distance of particle to center of the cell
	  dy = (ypos-min2)/d2 - yc;
	  dz = (zpos-min3)/d3 - zc;
	  
	  tx = num_1 - dx;
	  ty = num_1 - dy;
	  tz = num_1 - dz;
	  
	  // Take care this assumes periodic boundary conditions. This conserves Mass when using FFTs to deconvolve with CIC kernel, as the FFT assumes periodicity
	  
          real_prec tracer_weight=num_1;
          if (weight_prop==_MASS_)
            tracer_weight=Halo[n].mass;
          else if (weight_prop==_SAT_FRACTION_)
            tracer_weight=Halo[n].number_sub_structures;
          else if (weight_prop==_RS_)
            tracer_weight=Halo[n].rs;
          else if (weight_prop==_VMAX_)
            tracer_weight=Halo[n].vmax;
          else if (weight_prop==_VIRIAL_)
            tracer_weight=Halo[n].virial;
          else if (weight_prop==_SPIN_)
            tracer_weight=Halo[n].spin;

          delta[k+N3*(j+N2*i)]    += tracer_weight*tx*ty*tz;
          delta[k+N3*(j+N2*ii)]   += tracer_weight*dx*ty*tz;
          delta[k+N3*(jj+N2*i)]   += tracer_weight*tx*dy*tz;
          delta[kk+N3*(j+N2*i)]   += tracer_weight*tx*ty*dz;
          delta[k+N3*(jj+N2*ii)]  += tracer_weight*dx*dy*tz;
          delta[kk+N3*(j+N2*ii)]  += tracer_weight*dx*ty*dz;
          delta[kk+N3*(jj+N2*i)]  += tracer_weight*tx*dy*dz;
          delta[kk+N3*(jj+N2*ii)] += tracer_weight*dx*dy*dz;
	}
      else NLOSS++;	
    }
  
  ScreenOutput So;
  
  So.DONE();
  
  if(_COUNTS_ == weight_prop)
    {
      ULONG count=get_nobjects(delta);
      So.message_screen("Number of objects assigned to grid =",count);
      if(NLOSS!=0)
	So.message_screen("mass assignment found",NLOSS,"particles outside mesh boundary");    
    }

}

//============================================================================================================================
//============================================================================================================================
//============================================================================================================================


void getDensity_CICWEIGHT(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3, const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec>&zp, const vector<real_prec> &Par_mass, ULONG N_OBJ, vector<real_prec>&delta, bool weightmass)
{

    int NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);

    real_prec xc, yc, zc;
    real_prec dx, dy, dz, tx, tz, ty;    
    ULONG i, j, k, ii, jj, kk, n;

#pragma omp parallel for
  for (i=0;i<delta.size(); i++)
    delta[i]= 0.;  //-1 if we want to calculate overdensity
  


    ULONG NLOSS=0;
    for (n=0; n<N_OBJ; n++){

	//check if particle is in selected Domain, else discard it
      if((xp[n]>=min1 && xp[n]<min1+L1) && (yp[n]>=min2 && yp[n]<min2+L2) && (zp[n]>=min3 && zp[n]<min3+L3))
     {	
       i = static_cast<ULONG>(floor((xp[n]-min1)/d1)); // indices of the cell of the particle
       j = static_cast<ULONG>(floor((yp[n]-min2)/d2));
       k = static_cast<ULONG>(floor((zp[n]-min3)/d3));

	//printf("%f, %f, %f\n", xp[n], yp[n], zp[n]);		

       xc = static_cast<real_prec>(i)+static_cast<real_prec>(0.5); // centers of the cells
       yc = static_cast<real_prec>(j)+static_cast<real_prec>(0.5);
       zc = static_cast<real_prec>(k)+static_cast<real_prec>(0.5);

	dx = (xp[n]-min1)/d1 - xc; // distance of particle to center of the cell
	dy = (yp[n]-min2)/d2 - yc;
	dz = (zp[n]-min3)/d3 - zc;
	  	
	
	if (dx<0) 
	  {
	    dx=+.75; 
	    //dx*=-1.; 
	    	    
	    i-=1;
	    ii-=1;
	  } 

	if (dy<0) 
	  {
	    dy=+.75; 
	    //dy*=-1.; 
	    
	    j-=1;
	    jj-=1;
	  }  
		  
	if (dz<0) 
	  {
	    dz=+.75; 
	    //dz*=-1.; 
	    
	    k-=1;
	    kk-=1;
	  }

	if (i<0)
	  i+=N1;
	if (j<0)
	  j+=N2;
	if (k<0)
	  k+=N3;
	if (ii<0)
	  ii+=N1;
	if (jj<0)
	  jj+=N2;
	if (kk<0)
	  kk+=N3;
		  
	i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
	j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
	k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));

	ii = static_cast<ULONG>(fmod(real_prec(i+1),real_prec(N1)));
	jj = static_cast<ULONG>(fmod(real_prec(j+1),real_prec(N2)));
	kk = static_cast<ULONG>(fmod(real_prec(k+1),real_prec(N3)));
	  
	tx = static_cast<real_prec>(1.) - dx;
	ty = static_cast<real_prec>(1.) - dy;
	tz = static_cast<real_prec>(1.) - dz;
	
	real_prec mass=num_1;
	if (weightmass==true)
	  mass=Par_mass[n];

 
 	
	delta[k+N3*(j+N2*i)]    += mass*tx*ty*tz;
	delta[k+N3*(j+N2*ii)]   += mass*dx*ty*tz;
	delta[k+N3*(jj+N2*i)]   += mass*tx*dy*tz;
	delta[kk+N3*(j+N2*i)]   += mass*tx*ty*dz;
	delta[k+N3*(jj+N2*ii)]  += mass*dx*dy*tz;
	delta[kk+N3*(j+N2*ii)]  += mass*dx*ty*dz;
	delta[kk+N3*(jj+N2*i)]  += mass*tx*dy*dz;
	delta[kk+N3*(jj+N2*ii)] += mass*dx*dy*dz;
    }
     //else NLOSS++;	
  }
if(NLOSS!=0) cout << " >>> mass assignment found "<<NLOSS<<" particles outside mesh boundary...."<<endl<<endl;
}


// ==================================================================================================================
// ==================================================================================================================

void getDensity_TSC(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3, const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec>&zp, const vector<real_prec> & Par_mass, vector<real_prec>&delta, bool weightmass)
{

    int NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);
    real_prec xc, yc, zc;
    real_prec dx, dy, dz, hx0, hz0, hy0, hxp1,hyp1,hzp1, hxm1,hym1,hzm1;    
    ULONG i, j, k, ii, jj, kk,iii,jjj,kkk;
    ULONG n;
    ULONG N_OBJ=xp.size();
    
#define DELTA(i,j,k) delta[k+N3*(j+N2*i)]


#if !defined(_GET_INTERPOLATED_FIELDS_FROM_BIN_FILES_)
#pragma omp parallel for
    for (i=0;i<N1; i++)
      for (j=0; j<N2; j++)
        for (k=0; k<N3; k++)
           delta[k+N3*(j+N2*i)] = 0.;  //-1 if we want to calculate overdensity
#endif


    ULONG NLOSS=0;
    
    for (n=0; n<N_OBJ; n++){

	//check if particle is in selected Domain, else discard it
      if((xp[n]>=min1 && xp[n]<=min1+L1) && (yp[n]>=min2 && yp[n]<=min2+L2) && (zp[n]>=min3 && zp[n]<=min3+L3))
     {

       i = static_cast<ULONG>(floor((xp[n]-min1)/d1)); // indices of the cell of the particle
       j = static_cast<ULONG>(floor((yp[n]-min2)/d2));
       k = static_cast<ULONG>(floor((zp[n]-min3)/d3));


       i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
       j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
       k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));

       ii = static_cast<ULONG>(fmod(real_prec(i+1),real_prec(N1)));
       jj = static_cast<ULONG>(fmod(real_prec(j+1),real_prec(N2)));
       kk = static_cast<ULONG>(fmod(real_prec(k+1),real_prec(N3)));

       iii= static_cast<ULONG>(fmod(real_prec(i-1+N1),real_prec(N1)));
       jjj= static_cast<ULONG>(fmod(real_prec(j-1+N2),real_prec(N2)));
       kkk= static_cast<ULONG>(fmod(real_prec(k-1+N3),real_prec(N3)));
	

	//printf("%f, %f, %f\n", xp[n], yp[n], zp[n]);

       xc = static_cast<real_prec>(i)+static_cast<real_prec>(0.5); // centers of the cells
       yc = static_cast<real_prec>(j)+static_cast<real_prec>(0.5);
       zc = static_cast<real_prec>(k)+static_cast<real_prec>(0.5);

	dx = (xp[n]-min1)/d1 - xc; // distance of particle to center of the cell
	dy = (yp[n]-min2)/d2 - yc;
	dz = (zp[n]-min3)/d3 - zc;
		


	hx0=static_cast<real_prec>(0.75-dx*dx); // fraction of particle assigned
	hy0=static_cast<real_prec>(0.75-dy*dy); // fraction of particle assigned
	hz0=static_cast<real_prec>(0.75-dz*dz); // fraction of particle assigned

	hxp1=static_cast<real_prec>(0.5*(0.5+dx)*(0.5+dx)); // fraction of particle assigned
	hyp1=static_cast<real_prec>(0.5*(0.5+dy)*(0.5+dy)); // fraction of particle assigned
	hzp1=static_cast<real_prec>(0.5*(0.5+dz)*(0.5+dz)); // fraction of particle assigned

	hxm1=static_cast<real_prec>(0.5*(0.5-dx)*(0.5-dx)); // fraction of particle assigned
	hym1=static_cast<real_prec>(0.5*(0.5-dy)*(0.5-dy)); // fraction of particle assigned
	hzm1=static_cast<real_prec>(0.5*(0.5-dz)*(0.5-dz)); // fraction of particle assigned

	real_prec mass=num_1;
	if (weightmass==true)
	  mass=Par_mass[n];

	DELTA(i,j,k)    	+= mass*hx0*hy0*hz0;
	DELTA(ii,jj,kk) 	+= mass*hxp1*hyp1*hzp1;
	DELTA(iii,jjj,kkk) 	+= mass*hxm1*hym1*hzm1;

	DELTA(ii,j,k)   	+= mass*hxp1*hy0*hz0;
	DELTA(i,jj,k)   	+= mass*hx0*hyp1*hz0;
	DELTA(i,j,kk)   	+= mass*hx0*hy0*hzp1;

	DELTA(ii,jj,k)  	+= mass*hxp1*hyp1*hz0;
	DELTA(ii,j,kk)  	+= mass*hxp1*hy0*hzp1;
	DELTA(i,jj,kk)  	+= mass*hx0*hyp1*hzp1;

	DELTA(iii,jj,kk)   	+= mass*hxm1*hyp1*hzp1;
	DELTA(ii,jjj,kk)  	+= mass*hxp1*hym1*hzp1;
	DELTA(ii,jj,kkk)  	+= mass*hxp1*hyp1*hzm1;

	DELTA(iii,jjj,kk)  	+= mass*hxm1*hym1*hzp1;
	DELTA(iii,jj,kkk)  	+= mass*hxm1*hyp1*hzm1;
	DELTA(ii,jjj,kkk)  	+= mass*hxp1*hym1*hzm1;

	DELTA(i,jjj,kkk)   	+= mass*hx0*hym1*hzm1;
	DELTA(iii,j,kkk)  	+= mass*hxm1*hy0*hzm1;
	DELTA(iii,jjj,k)  	+= mass*hxm1*hym1*hz0;

	DELTA(i,j,kkk)  	+= mass*hx0*hy0*hzm1;
	DELTA(i,jjj,k)  	+= mass*hx0*hym1*hz0;
	DELTA(iii,j,k)  	+= mass*hxm1*hy0*hz0;
	
	DELTA(i,jj,kkk)   	+= mass*hx0*hyp1*hzm1;
	DELTA(ii,j,kkk)  	+= mass*hxp1*hy0*hzm1;

	DELTA(iii,jj,k)   	+= mass*hxm1*hyp1*hz0;
	DELTA(ii,jjj,k)  	+= mass*hxp1*hym1*hz0;

	DELTA(i,jjj,kk)   	+= mass*hx0*hym1*hzp1;
	DELTA(iii,j,kk)  	+= mass*hxm1*hy0*hzp1;
    }	
  }
if(NLOSS!=0) cout << " >>> mass assignment found "<<NLOSS<<" particles outside mesh boundary...."<<endl<<endl;
// cout <<"TOTAL NUMBER OF GALAXIES = "<< N_OBJ <<endl;
}



// ==================================================================================================================
// ==================================================================================================================

void getDensity_TSCv(s_params_box_mas *params, const vector<s_Halo>&Halo, vector<real_prec>&delta, string weight_prop)
{
  
    int NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);

  real_prec xc, yc, zc;
    real_prec dx, dy, dz, hx0, hz0, hy0, hxp1,hyp1,hzp1, hxm1,hym1,hzm1;    
    ULONG i, j, k, ii, jj, kk,iii,jjj,kkk;
    ULONG n;
    ULONG N_OBJ=Halo.size();
    
    real_prec min1=params->min1;
    real_prec min2=params->min2;
    real_prec min3=params->min3;
    real_prec L1=params->Lbox;
    real_prec L2=params->Lbox;
    real_prec L3=params->Lbox;
    real_prec N1=params->Nft;
    real_prec N2=params->Nft;
    real_prec N3=params->Nft;

    real_prec d1=params->d1;
    real_prec d2=params->d2;
    real_prec d3=params->d3;



  
#pragma omp parallel for
  for (i=0;i<delta.size(); i++)
    delta[i]= 0.;  //-1 if we want to calculate overdensity
  
    ULONG NLOSS=0;
    
  for (n=0; n<N_OBJ; n++){
      
//check if particle is in selected Domain, else discard it
    if((Halo[n].coord1>=min1 && Halo[n].coord1<=min1+L1) && (Halo[n].coord2>=min2 &&  Halo[n].coord2<=min2+L2) && (Halo[n].coord3>=min3 && Halo[n].coord3<=min3+L3))
		{
	  
	   i = static_cast<ULONG>(floor((Halo[n].coord1-min1)/d1)); // indices of the cell of the particle
	   j = static_cast<ULONG>(floor((Halo[n].coord2-min2)/d2));
	   k = static_cast<ULONG>(floor((Halo[n].coord3-min3)/d3));
	  
	  
	   i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
	   j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
	   k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));
	  
	   ii = static_cast<ULONG>(fmod(real_prec(i+1),real_prec(N1)));
	   jj = static_cast<ULONG>(fmod(real_prec(j+1),real_prec(N2)));
	   kk = static_cast<ULONG>(fmod(real_prec(k+1),real_prec(N3)));
	  
	   iii= static_cast<ULONG>(fmod(real_prec(i-1+N1),real_prec(N1)));
	   jjj= static_cast<ULONG>(fmod(real_prec(j-1+N2),real_prec(N2)));
	   kkk= static_cast<ULONG>(fmod(real_prec(k-1+N3),real_prec(N3)));
	  
	  
	  //printf("%f, %f, %f\n", xp[n], yp[n], zp[n]);
	  
	   xc = static_cast<real_prec>(i)+static_cast<real_prec>(0.5); // centers of the cells
	   yc = static_cast<real_prec>(j)+static_cast<real_prec>(0.5);
	   zc = static_cast<real_prec>(k)+static_cast<real_prec>(0.5);
	  
	   dx = (Halo[n].coord1-min1)/d1 - xc; // distance of particle to center of the cell
	   dy = (Halo[n].coord2-min2)/d2 - yc;
	   dz = (Halo[n].coord3-min3)/d3 - zc;
	  
	  
	   hx0=static_cast<real_prec>(0.75-dx*dx); // fraction of particle assigned
	   hy0=static_cast<real_prec>(0.75-dy*dy); // fraction of particle assigned
	   hz0=static_cast<real_prec>(0.75-dz*dz); // fraction of particle assigned
	  
	   hxp1=static_cast<real_prec>(0.5*(0.5+dx)*(0.5+dx)); // fraction of particle assigned
	   hyp1=static_cast<real_prec>(0.5*(0.5+dy)*(0.5+dy)); // fraction of particle assigned
	   hzp1=static_cast<real_prec>(0.5*(0.5+dz)*(0.5+dz)); // fraction of particle assigned
	  
	   hxm1=static_cast<real_prec>(0.5*(0.5-dx)*(0.5-dx)); // fraction of particle assigned
	   hym1=static_cast<real_prec>(0.5*(0.5-dy)*(0.5-dy)); // fraction of particle assigned
	   hzm1=static_cast<real_prec>(0.5*(0.5-dz)*(0.5-dz)); // fraction of particle assigned
	  
	   real_prec mass=num_1;

	   if (weight_prop==_MASS_)
	     mass=Halo[n].mass;
	   else if (weight_prop==_SAT_FRACTION_)
	     mass=Halo[n].number_sub_structures;
	  
		  
	  DELTA(i,j,k)    	+= mass*hx0*hy0*hz0;
	  DELTA(ii,jj,kk) 	+= mass*hxp1*hyp1*hzp1;
	  DELTA(iii,jjj,kkk) 	+= mass*hxm1*hym1*hzm1;
	  
	  DELTA(ii,j,k)   	+= mass*hxp1*hy0*hz0;
	  DELTA(i,jj,k)   	+= mass*hx0*hyp1*hz0;
	  DELTA(i,j,kk)   	+= mass*hx0*hy0*hzp1;
	  
	  DELTA(ii,jj,k)  	+= mass*hxp1*hyp1*hz0;
	  DELTA(ii,j,kk)  	+= mass*hxp1*hy0*hzp1;
	  DELTA(i,jj,kk)  	+= mass*hx0*hyp1*hzp1;
	  
	  DELTA(iii,jj,kk)   	+= mass*hxm1*hyp1*hzp1;
	  DELTA(ii,jjj,kk)  	+= mass*hxp1*hym1*hzp1;
	  DELTA(ii,jj,kkk)  	+= mass*hxp1*hyp1*hzm1;
	  
	  DELTA(iii,jjj,kk)  	+= mass*hxm1*hym1*hzp1;
	  DELTA(iii,jj,kkk)  	+= mass*hxm1*hyp1*hzm1;
	  DELTA(ii,jjj,kkk)  	+= mass*hxp1*hym1*hzm1;
	  
	  DELTA(i,jjj,kkk)   	+= mass*hx0*hym1*hzm1;
	  DELTA(iii,j,kkk)  	+= mass*hxm1*hy0*hzm1;
	  DELTA(iii,jjj,k)  	+= mass*hxm1*hym1*hz0;
	  
	  DELTA(i,j,kkk)  	+= mass*hx0*hy0*hzm1;
	  DELTA(i,jjj,k)  	+= mass*hx0*hym1*hz0;
	  DELTA(iii,j,k)  	+= mass*hxm1*hy0*hz0;
	  
	  DELTA(i,jj,kkk)   	+= mass*hx0*hyp1*hzm1;
	  DELTA(ii,j,kkk)  	+= mass*hxp1*hy0*hzm1;
	  
	  DELTA(iii,jj,k)   	+= mass*hxm1*hyp1*hz0;
	  DELTA(ii,jjj,k)  	+= mass*hxp1*hym1*hz0;
	  
	  DELTA(i,jjj,kk)   	+= mass*hx0*hym1*hzp1;
	  DELTA(iii,j,kk)  	+= mass*hxm1*hy0*hzp1;
	}	
    }

    ScreenOutput So;
    So.DONE();

    if(_COUNTS_==weight_prop)
      {
	
	ULONG count=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:count)
#endif
	for(ULONG i=0;i<delta.size();++i)
	  count+=delta[i];
	So.message_screen("Number of objects assigned to grid =",count);
      }
    
    if(NLOSS!=0) So.message_screen("mass assignment found",NLOSS,"particles outside mesh boundary");
    
}


// ============================================================================================================================
// ============================================================================================================================
// ============================================================================================================================


void cellbound(ULONG N1,ULONG N2,ULONG N3,real_prec L1,real_prec L2,real_prec L3, vector<real_prec>&v1, vector<real_prec>&v2, vector<real_prec>&v3)
{

    int NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);

    ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);

  /* start: interpolate from cell center to cell boundaries */
  vector<real_prec> vx(N,0),vy(N,0),vz(N,0);
      
      for (ULONG i=0;i<N1; i++)
	for (ULONG j=0; j<N2; j++)
	  for (ULONG k=0; k<N3; k++)
	    {     
	      ULONG l=k+N3*(j+N2*i);
	      ULONG m=k-1+N3*(j-1+N2*(i-1));
	      
	      if (i>0 && j>0 && k>0)
		{
		  vx[l]=static_cast<real_prec>(0.5)*(v1[m]+v1[l]);
		  vy[l]=static_cast<real_prec>(0.5)*(v2[m]+v2[l]);
		  vz[l]=static_cast<real_prec>(0.5)*(v3[m]+v3[l]);
		}
	      
	      /* periodic boundary conditions */
	      
	      ULONG ii=k-1+N3*(j-1+N2*(N1-1));
	      ULONG jj=k-1+N3*(N2-1+N2*(i-1));
	      ULONG kk=N3-1+N3*(j-1+N2*(i-1));
	      ULONG ij=k-1+N3*(N2-1+N2*(N1-1));
	      ULONG ik=N3-1+N3*(j-1+N2*(N1-1));
	      ULONG jk=N3-1+N3*(N2-1+N2*(i-1));
	      ULONG ijk=N3-1+N3*(N2-1+N2*(N1-1));
	      
	      if (i==0 && j>0 && k>0)
		{
		  vx[l]=static_cast<real_prec>(0.5)*(v1[ii]+v1[l]);
		  vy[l]=static_cast<real_prec>(0.5)*(v2[ii]+v2[l]);
		  vz[l]=static_cast<real_prec>(0.5)*(v3[ii]+v3[l]);
		}
	      if (i==0 && j==0 && k>0)
		{
		  vx[l]=static_cast<real_prec>(0.5)*(v1[ij]+v1[l]);
		  vy[l]=static_cast<real_prec>(0.5)*(v2[ij]+v2[l]);
		  vz[l]=static_cast<real_prec>(0.5)*(v3[ij]+v3[l]);
		}
	      if (i==0 && j>0 && k==0)
		{
		  vx[l]=static_cast<real_prec>(0.5)*(v1[ik]+v1[l]);
		  vy[l]=static_cast<real_prec>(0.5)*(v2[ik]+v2[l]);
		  vz[l]=static_cast<real_prec>(0.5)*(v3[ik]+v3[l]);
		}
	      if (i==0 && j==0 && k==0)
		{
		  vx[l]=static_cast<real_prec>(0.5)*(v1[ijk]+v1[l]);
		  vy[l]=static_cast<real_prec>(0.5)*(v2[ijk]+v2[l]);
		  vz[l]=static_cast<real_prec>(0.5)*(v3[ijk]+v3[l]);
		}
	      if (i>0 && j==0 && k==0)
		{
		  vx[l]=static_cast<real_prec>(0.5)*(v1[jk]+v1[l]);
		  vy[l]=static_cast<real_prec>(0.5)*(v2[jk]+v2[l]);
		  vz[l]=static_cast<real_prec>(0.5)*(v3[jk]+v3[l]);
		}
	      if (i>0 && j==0 && k>0)
		{
		  vx[l]=static_cast<real_prec>(0.5)*(v1[jj]+v1[l]);
		  vy[l]=static_cast<real_prec>(0.5)*(v2[jj]+v2[l]);
		  vz[l]=static_cast<real_prec>(0.5)*(v3[jj]+v3[l]);
		}
	      if (i>0 && j>0 && k==0)
		{
		  vx[l]=static_cast<real_prec>(0.5)*(v1[kk]+v1[l]);
		  vy[l]=static_cast<real_prec>(0.5)*(v2[kk]+v2[l]);
		  vz[l]=static_cast<real_prec>(0.5)*(v3[kk]+v3[l]);
		}
	    }
	  
      for (ULONG l=0;l<N; l++)
	{
	  v1[l]=vx[l];
	  v2[l]=vy[l];
	  v3[l]=vz[l];
	}
      /* end: interpolate from cell center to cell boundaries */   
}

void cellboundcomp(ULONG N1,ULONG N2,ULONG N3,vector<real_prec>&vi)
{
  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);

  /* start: interpolate from cell center to cell boundaries */
  vector<real_prec> viout(N,0);
      
      for (ULONG i=0;i<N1; i++)
	for (ULONG j=0; j<N2; j++)
	  for (ULONG k=0; k<N3; k++)
	    {     
	      ULONG l=k+N3*(j+N2*i);
	      ULONG m=k-1+N3*(j-1+N2*(i-1));
	      
	      if (i>0 && j>0 && k>0)
		{
		  viout[l]=num_0_5*(vi[m]+vi[l]);
		}
	      
	      /* periodic boundary conditions */
	      
	      ULONG ii=k-1+N3*(j-1+N2*(N1-1));
	      ULONG jj=k-1+N3*(N2-1+N2*(i-1));
	      ULONG kk=N3-1+N3*(j-1+N2*(i-1));
	      ULONG ij=k-1+N3*(N2-1+N2*(N1-1));
	      ULONG ik=N3-1+N3*(j-1+N2*(N1-1));
	      ULONG jk=N3-1+N3*(N2-1+N2*(i-1));
	      ULONG ijk=N3-1+N3*(N2-1+N2*(N1-1));
	      
	      if (i==0 && j>0 && k>0)
		{
		  viout[l]=num_0_5*(vi[ii]+vi[l]);
		}
	      if (i==0 && j==0 && k>0)
		{
		  viout[l]=num_0_5*(vi[ij]+vi[l]);
		}
	      if (i==0 && j>0 && k==0)
		{
		  viout[l]=num_0_5*(vi[ik]+vi[l]);
		}
	      if (i==0 && j==0 && k==0)
		{
		  viout[l]=num_0_5*(vi[ijk]+vi[l]);
		}
	      if (i>0 && j==0 && k==0)
		{
		  viout[l]=num_0_5*(vi[jk]+vi[l]);
		}
	      if (i>0 && j==0 && k>0)
		{
		  viout[l]=num_0_5*(vi[jj]+vi[l]);
		}
	      if (i>0 && j>0 && k==0)
		{
		  viout[l]=num_0_5*(vi[kk]+vi[l]);
		}
	    }
	  
      for (ULONG l=0;l<N; l++)
	{
	  vi[l]=viout[l];
	}
      /* end: interpolate from cell center to cell boundaries */   
}


// ============================================================================================================================
// ============================================================================================================================
// ============================================================================================================================



void cellbound_comp(ULONG N1,ULONG N2,ULONG N3,vector<real_prec>&vi, vector<real_prec>&viout)
{
  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);

  /* start: interpolate from cell center to cell boundaries */
      
      for (ULONG i=0;i<N1; i++)
	for (ULONG j=0; j<N2; j++)
	  for (ULONG k=0; k<N3; k++)
	    {     
	      ULONG l=k+N3*(j+N2*i);
	      ULONG m=k-1+N3*(j-1+N2*(i-1));
	      
	      if (i>0 && j>0 && k>0)
		{
		  viout[l]=num_0_5*(vi[m]+vi[l]);
		}
	      
	      /* periodic boundary conditions */
	      
	      ULONG ii=k-1+N3*(j-1+N2*(N1-1));
	      ULONG jj=k-1+N3*(N2-1+N2*(i-1));
	      ULONG kk=N3-1+N3*(j-1+N2*(i-1));
	      ULONG ij=k-1+N3*(N2-1+N2*(N1-1));
	      ULONG ik=N3-1+N3*(j-1+N2*(N1-1));
	      ULONG jk=N3-1+N3*(N2-1+N2*(i-1));
	      ULONG ijk=N3-1+N3*(N2-1+N2*(N1-1));
	      
	      if (i==0 && j>0 && k>0)
		{
		  viout[l]=num_0_5*(vi[ii]+vi[l]);
		}
	      if (i==0 && j==0 && k>0)
		{
		  viout[l]=num_0_5*(vi[ij]+vi[l]);
		}
	      if (i==0 && j>0 && k==0)
		{
		  viout[l]=num_0_5*(vi[ik]+vi[l]);
		}
	      if (i==0 && j==0 && k==0)
		{
		  viout[l]=num_0_5*(vi[ijk]+vi[l]);
		}
	      if (i>0 && j==0 && k==0)
		{
		  viout[l]=num_0_5*(vi[jk]+vi[l]);
		}
	      if (i>0 && j==0 && k>0)
		{
		  viout[l]=num_0_5*(vi[jj]+vi[l]);
		}
	      if (i>0 && j>0 && k==0)
		{
		  viout[l]=num_0_5*(vi[kk]+vi[l]);
		}
	    }
	  
#pragma omp parallel for
     for (ULONG l=0;l<N; l++)
	{
	  vi[l]=viout[l];
	}
      /* end: interpolate from cell center to cell boundaries */   
}

// ============================================================================================================================
// ============================================================================================================================
// ============================================================================================================================


void MAS_NEW(s_params_box_mas *params, const vector<s_Halo>&Halo,  string weight_prop, vector<real_prec>&delta)
{
  
  switch (params->masskernel)
    {
    case 0:
      getDensity_NGPv(params, Halo, delta,weight_prop);
      break;
    case 1:
      getDensity_CICv(params, Halo,delta,weight_prop);
      break;
    case 2:
      getDensity_TSCv(params, Halo,delta,weight_prop);
      break;
    }
}


