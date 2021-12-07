/** @file CoordinateSystem.cpp
 *
 *  @brief Methods of the class CoordinateSystem
 *
 *  This file contains the implementation of the methods of the class
 *  CoordinateSystem FftwFunctions, used to transform coordiante system in the catalogues
 */



#include "../Headers/CoordinateSystem.h"
#include "../Headers/ScreenOutput.h"


// ****************************************************************************************
// ****************************************************************************************
void CoordinateSystem::compute_box_size(bool new_los, void *go, real_prec *Lside)
{
  // Given the galaxy sample, compute the minimum and maximum values 
  // of the cartesian coordinates
  // in order to generate the size of the box (or cuboid) for the DFT
  struct s_galaxy_operationsF * go_p= (struct s_galaxy_operationsF *)go;
  Cosmology Cf;
  CoordinateSystem Go;

  real_prec aXMAX=-1e6,aYMAX=-1e6, aZMAX=-1e6,xmax,ymax,zmax;
  real_prec aXMIN=+1e6,aYMIN=+1e6, aZMIN=+1e6,xmin,ymin,zmin;
  //real_prec mean_ra, mean_dec, x, y, z;

  time_t start;
  time (&start);
  
  
  cout<<BLUE"************************************************************************************"<<RESET<<endl; 
  cout<<CYAN<<"Computing size of smaller CUBE fully embedding the real sample"<<RESET<<endl;    
  // ***********************************************************************************
  // ***********************************************************************************
  vector< real_prec > prop=go_p->properties;
  vector<real_prec>zz=go_p->zz;
  vector<real_prec>rc=go_p->rc;
  int n_columns = go_p->n_columns;
  int n_lines=prop.size() / n_columns;
  string angles_units =go_p->angles_units;
  
  //real_prec fac= (angles_units=="R" ? M_PI/180.0:1.0);
  
  real_prec *pprop;
  //real_prec ra_s, dec_s;

  this->XMIN=0;
  this->YMIN=0;
  this->ZMIN=0;

  
  
  
    // if(new_los){
    //   mean_ra=0;
    //   mean_dec=0;
    //   for(int i=0;i<n_lines;i++){
    // 	if(angles_units=="R"){
    // 	  ra_s=prop[i][go_p->i_coord1]*180.0/acos(-1.0);                                                               
    // 	  dec_s=prop[i][go_p->i_coord2]*180.0/acos(-1.0);
    // 	}  
    // 	else{
    // 	  ra_s=prop[i][go_p->i_coord1];                               
    // 	  dec_s=prop[i][go_p->i_coord2];
    // 	}
    // 	mean_ra+=ra_s;                       //use explicitely the variables related to the galaxy catalog, i.e, _g 
    // 	mean_dec+=sin(dec_s*fac); 
    //   }         
      
    //   mean_ra/=(real_prec)n_lines;
    //   mean_dec=asin(mean_dec/n_lines)/(M_PI/180.0);
    //   *RA_los=mean_ra;
    //   *Dec_los=mean_dec;
      
    //   cout<<"Mean RA  [deg] =  "<<mean_ra<<endl;
    //   cout<<"Mean Dec [deg] =  "<<mean_dec<<endl;
    //   cout<<"Rotating the celestial sphere such that the z-direction points towards "<<endl;
    //   cout<<"the baricenter at (mean_alpha, mean_dec) "<<endl;
    //   cout<<RED<<"HERE THERE IS A BUG. MEAN COORDINATES ARE NOTE WELL DEFINED."<<RESET<<endl;
    //   cout<<" "<<endl;
    // }
  
  for(int i=0;i<n_lines;i++){
    pprop = &prop[i * n_columns];
    //real_prec rr;
    real_prec x=pprop[go_p->i_coord1];
    real_prec y=pprop[go_p->i_coord2];
    real_prec z=pprop[go_p->i_coord3];

    // if(go_p->sys_coord==1)rr=prop[i][go_p->i_coord3];
    // else if(go_p->sys_coord==2)rr=gsl_inter_new(zz, rc, prop[i][go_p->i_coord3]);
    // else if(go_p->sys_coord==3)rr=prop[i][go_p->i_coord3];
    // if(angles_units=="R"){
    //   ra_s=prop[i][go_p->i_coord1]*180.0/M_PI;
    //   dec_s=prop[i][go_p->i_coord2]*180.0/M_PI;
    // }
    // else{
    //   ra_s=prop[i][go_p->i_coord1];
    //   dec_s=prop[i][go_p->i_coord2];
    // }
    
    //Convert to cartesian coordinates:
    //  equatorial_to_cartesian(ra_s,dec_s,rr,&x,&y,&z);
    
    // Find min and max coordinates
    xmax=max(x, aXMAX);
    ymax=max(y, aYMAX);
    zmax=max(z, aZMAX);
    xmin=min(x, aXMIN);
    ymin=min(y, aYMIN);
    zmin=min(z, aZMIN);
    aXMAX=xmax;
    aYMAX=ymax;
    aZMAX=zmax;
    aXMIN=xmin;
    aYMIN=ymin;
    aZMIN=zmin; 
  }
  
  this->XMIN=aXMIN;
  this->YMIN=aYMIN;
  this->ZMIN=aZMIN;


  real_prec llx=fabs(aXMAX-aXMIN);
  real_prec lly=fabs(aYMAX-aYMIN);
  real_prec llz=fabs(aZMAX-aZMIN);

  cout<<aXMAX<<" "<<this->XMIN<<"  "<<aYMAX<<"  "<<this->YMIN<<"  "<<aZMAX<<"  "<<this->ZMIN<<endl;
  cout<<llx<<"  "<<lly<<"  "<<llz<<endl;



  llx=max(llx,lly);
  lly=max(lly,llz);
  llz=max(llx,lly);

  *Lside=llz;


  cout<<RED;
  time_t end; time (&end);
  real_prec diff = difftime(end,start);
  if (diff<60) cout <<"Lapse: "<<diff<<" secs"<<endl;
  else if (diff<3600) cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
  else cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
  cout<<RESET;


}



// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
void CoordinateSystem::equatorial_to_cartesian(real_prec ra, real_prec dec, real_prec r, real_prec &x, real_prec &y, real_prec &z)
{
  // ********************************************************************************
  // Convert equatorial coordinates, i.e, right ascension, declination and distance * 
  // to cartessian coordinates. Input values of RA and Dec are in degrees,          *
  // the output is in the units of r                                                *
  // ********************************************************************************
  real_prec fac=M_PI/180.;
  x  =r*cos(ra*fac)*sin(0.5*M_PI-dec*fac);
  y  =r*sin(ra*fac)*sin(0.5*M_PI-dec*fac);
  z  =r*cos(0.5*M_PI-dec*fac);

  // // change x->z 
  // // and z to -x
  // this to mimic the column order that uses Florian
  // in which x varies faster than z
  // *z  =-r*cos(ra*fac)*sa;
  // *y  =r*sin(ra*fac)*sa;
  // *x  =r*ca;
  
}



// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
void CoordinateSystem::new_equatorial_to_cartesian(real_prec ra, real_prec dec, real_prec r, real_prec m_ra, real_prec m_dec, real_prec &x, real_prec &y, real_prec &z){
  // ********************************************************************************
  // Convert equatorial coordinates, i.e, right ascension, declination and distance * 
  // to a new spherical coordinate system with the z-axis pointing towards          *
  // the point identified with m_ra and m_dec. Conversion done such that 
  // that the z-direction          
  // points towards the new line of sight 
  // Input values of RA and Dec are in degrees,                                      * 
  // the output is in the units of r                                                *
  // ********************************************************************************

  real_prec xx, yy, zz;
  // Convert to cartessian
  real_prec new_dec,new_ra;
  equatorial_to_equatorial(ra, dec, m_ra,m_dec,&new_ra,&new_dec); 
  equatorial_to_cartesian(new_ra,new_dec,r, x, y, z);
}


// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
void CoordinateSystem:: equatorial_to_galactic(real_prec alpha, real_prec delta, real_prec *b, real_prec *l)
{
  // ********************************************************************************
  // Convert equatorial coordinates, i.e, right ascension, declination and distance *
  // to galactic coordinates . Input values of RA and Dec are in degrees,           *
  // the output are in radians                                                      *
  // ********************************************************************************
  real_prec fac = M_PI/180.;                         /* Convierte grados a radianes*/
  real_prec lly,llx,bb;
  real_prec alpha_p=fac*192.859508;                        /*Right ascention of galactic north pole in radians*/
  real_prec delta_p=fac*27.128336;                         /*Declination of galactic north pole in radians*/
  real_prec  l_n   =fac*122.932;                           /*Galactic longitude of the celestial pole, in radians*/
  bb = sin(delta)*sin(delta_p)+cos(delta)*cos(delta_p)*cos(alpha-alpha_p);
  lly= cos(delta)*sin(alpha-alpha_p);
  llx=-cos(delta)*sin(delta_p)*cos(alpha-alpha_p)+sin(delta)*cos(delta_p);
  *l = l_n-atan2(lly,llx);
  *b = asin(bb);
}

// ****************************************************************************************
// ****************************************************************************************



void CoordinateSystem:: equatorial_to_equatorial(real_prec old_ra, real_prec old_dec, real_prec m_ra, real_prec m_dec, real_prec *new_ra, real_prec *new_dec)
{

  // ********************************************************************************
  // Convert equatorial coordinates, i.e, right ascension and declination           *
  // to a new equatorial coordinate in which the axis north pole                    *
  // has angular coordinates  real_prec m_ra, real_prec m_dec.                            *
  // Input values of RA and Dec in Degrees,                                         *
  // output in Degrees                                                              *
  // ********************************************************************************
  
  real_prec fac = M_PI/180.0;                         /* Convierte grados a radianes*/

  real_prec aux1=asin(sin(old_dec*fac)*sin(m_dec*fac)-cos(old_dec*fac)*cos(m_dec*fac)*sin(old_ra*fac-m_ra*fac));
  real_prec aux2=cos(old_dec*fac)*cos(old_ra*fac-m_ra*fac)/cos(aux1);
  real_prec aux3=0;
  aux2=aux3-acos(aux2);
  *new_dec=aux1/fac;
  *new_ra =aux2/fac;
}




