// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************
// This file contains some functions, based on GSL subroutines, used  in the measurements  **
// of the power spectrum                                                                   **
// Developer:                                                                              **
// Andres Balaguera Antolinez                                                              **
// abalant@gmail.com                                                                       **
// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************
# include "../Headers/NumericalMethods.h"
# include "../Headers/FileOutput.h"
using namespace std;


// ********************************************************************
// ********************************************************************
// ********************************************************************
// ********************************************************************
// ********************************************************************
// ********************************************************************

// This function generates an integer in the range (0,Nmax-1)
int my_gsl_rng_uniform_(gsl_rng *r, int Nmax){
    real_prec xr=gsl_rng_uniform(r);
    real_prec delta= (1.0)/static_cast<real_prec>(Nmax);
    int val = get_bin(xr, 0., Nmax,delta, true);
    return val;

}

ULONG factorial(int n){
  if(n<=0 || n==1 ) return 1;
  else return tgammal(n+1);
}

// ********************************************************************
// ********************************************************************
// This function does the same as
//gsl_ran_shuffle(rn,&cells_id_still_to_assign[0],cells_id_still_to_assign.size(),sizeof(ULONG));

void randomize_vector(vector<ULONG>&data)
{
    const gsl_rng_type *Tn = gsl_rng_default;
    gsl_rng *rn = gsl_rng_alloc (Tn);
    gsl_rng_default_seed=155;

    int n_aux=data.size();
    vector<ULONG>cells_id_still_to_assign_aux(n_aux,0);  //container to allocate the ID of the cells with particles still to get masses
    vector<bool>cells_id_still_to_assign_chosen(n_aux,false);  //container to allocate the ID of the cells with particles still to get masses

    for(int ide=0; ide < n_aux ;++ide)
      {
        bool flag=false;
        while(false==flag)
          {
            int jk= gsl_rng_uniform_int(rn,n_aux);
            ULONG new_id=data[jk];
            if(cells_id_still_to_assign_chosen[jk]==false)
              {
                cells_id_still_to_assign_aux[ide]=new_id;  //container to allocate the ID of the cells with particles still to get masses
                cells_id_still_to_assign_chosen[jk]=true;  //container to allocate the ID of the cells with particles still to get masses
                flag=true;
              }
          }
      }

    for(ULONG ide=0; ide < n_aux ;++ide)
      data[ide]=cells_id_still_to_assign_aux[ide];  //container to allocate the ID of the cells with particles still to get masses

    cells_id_still_to_assign_aux.clear();
    cells_id_still_to_assign_aux.shrink_to_fit();

    cells_id_still_to_assign_chosen.clear();
    cells_id_still_to_assign_chosen.shrink_to_fit();


}
// ********************************************************************
// ********************************************************************
// ********************************************************************
void  get_high_res_id_from_low_res_id(int Nft_high, int Nft_low,vector<ULONG>&id_info)
{

    double delta_low=static_cast<double>(Nft_high)/static_cast<double>(Nft_low);
    for(int i=0;i<Nft_high ;++i)
      for(int j=0;j<Nft_high ;++j)
        for(int k=0;k<Nft_high ;++k)
         {
            // now get the coords of the i, j , k in the lwo res:
           ULONG il = static_cast<ULONG>(floor((i+0.5)/delta_low)); // indices of the cell of the particle
           ULONG jl = static_cast<ULONG>(floor((j+0.5)/delta_low));
           ULONG kl = static_cast<ULONG>(floor((k+0.5)/delta_low));
           ULONG id_h   = index_3d(i, j, k, Nft_high, Nft_high);
           id_info[id_h]= index_3d(il,jl,kl,Nft_low,Nft_low);
         }
}


// ********************************************************************
// ********************************************************************
// ********************************************************************


//void  get_low_res_id_from_high_res_id(int Nft, int Nft_low, vector<ULONG>&ID_low, vector<s_cell_info>&cell_info_low)
void  get_low_res_id_from_high_res_id(int Nft, int Nft_low, vector<s_cell_info_reduced>&cell_info_low)
{

    double delta_low=static_cast<double>(Nft)/static_cast<double>(Nft_low);
    for(int i=0;i<Nft ;++i)
      for(int j=0;j<Nft ;++j)
        for(int k=0;k<Nft ;++k)
         {
            // now get the coords of the i, j , k in the lwo res:
           ULONG il = static_cast<ULONG>(floor((i+0.5)/delta_low)); // indices of the cell of the particle
           ULONG jl = static_cast<ULONG>(floor((j+0.5)/delta_low));
           ULONG kl = static_cast<ULONG>(floor((k+0.5)/delta_low));
           ULONG id_l = index_3d(il,jl,kl,Nft_low,Nft_low);
           ULONG id   = index_3d(i, j, k, Nft, Nft);
           cell_info_low[id_l].gal_index.push_back(id);  // allocate for each lowres ID all those high_res id living inside
         }
}


// ********************************************************************
// ********************************************************************
// ********************************************************************
void get_neighbour_cells(int Nft, int N_cells_back_forth,vector<s_nearest_cells>&nearest_cells_to_cell){

#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif

#if defined _USE_NEIGHBOURS_ || defined _GET_DIST_MIN_SEP_REF_ || defined _GET_DIST_MIN_SEP_MOCK_
  int max_neigh_per_dim = 2*N_cells_back_forth+1; // Number of neighbouring cells pr dimension, including the centrak cell;
  int N_Neigh_cells=pow(max_neigh_per_dim+1,2)+ pow(max_neigh_per_dim+1,2)  + (pow(max_neigh_per_dim+1,2)); // total number of cells to explore around a cell, excluding the central
#else
  int max_neigh_per_dim =2*N_cells_back_forth+1; // Number of neighbouring cells pr dimension, including the central cell;
  int N_Neigh_cells=pow(max_neigh_per_dim,2)+ pow(max_neigh_per_dim,2)  + (pow(max_neigh_per_dim,2)-1); // total number of cells to explore around a cell, excluding the central
#endif

  vector<int> index_cells_dime(max_neigh_per_dim,0);

  for(int i=0;i< max_neigh_per_dim ;++i)
     index_cells_dime[i]= i-N_cells_back_forth;  // e.g., 1, 0, -1 if N_cells_back_forth=1 and the central excluded


 #pragma omp parallel for collapse(3)
  for(int i=0;i<Nft ;++i)
    for(int j=0;j<Nft ;++j)
      for(int k=0;k<Nft ;++k)
        {
          ULONG ID=index_3d(i,j,k,Nft,Nft);//get ID of cell

          nearest_cells_to_cell[ID].close_cell.resize(N_Neigh_cells,0);

#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
          nearest_cells_to_cell[ID].bc_x.resize(N_Neigh_cells,0);
          nearest_cells_to_cell[ID].bc_y.resize(N_Neigh_cells,0);
          nearest_cells_to_cell[ID].bc_z.resize(N_Neigh_cells,0);
#endif
          ULONG count=0;
          for(int idx =0; idx < max_neigh_per_dim; ++idx)  // loop in the x-direction
            {
              int new_ni = i - index_cells_dime[idx];
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
                 int aux_bc_x=0;
#endif
                 if(new_ni<0)
                  {
                    new_ni+= Nft;
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
                     aux_bc_x=-1;
#endif
                  }
                 if(new_ni>= Nft)
                  {
                     new_ni-=Nft;
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
                     aux_bc_x=1;
#endif
                  }
               for(int idy =0; idy < max_neigh_per_dim; ++idy) // loop in the y-direction
                 {
                   int new_nj = j - index_cells_dime[idy];
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
                   int aux_bc_y=0;
#endif
                   if(new_nj<0)
                     {
                       new_nj+=Nft;
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
                       aux_bc_y=-1;
#endif
                     }//si la celda esta por detrñas, ponga -1
                    if(new_nj>= Nft)
                      {
                       new_nj-=Nft;
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
                       aux_bc_y=1;
#endif
                      }
                    
                    for(int idz =0; idz < max_neigh_per_dim; ++idz) // loop in the z-direction
                      {
                        int new_nk = k - index_cells_dime[idz];
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
                          int aux_bc_z=0;
#endif
                          if(new_nk<0)
                           {
                             new_nk+=Nft;
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
                             aux_bc_z=-1;
#endif
                           }
                          if(new_nk>=Nft)
                            {
                              new_nk-=Nft;
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
                              aux_bc_z=1;
#endif
                          }
                        
                        ULONG new_id=index_3d(new_ni,new_nj,new_nk,Nft,Nft); // Get the ID of the neighbouring cells;
                        if(new_id!=ID) // This explicitely excludes the ID of the same cell where we are now.
                         {
                           nearest_cells_to_cell[ID].close_cell[count]=new_id;                        
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
                           nearest_cells_to_cell[ID].bc_x[count]=aux_bc_x;
                           nearest_cells_to_cell[ID].bc_y[count]=aux_bc_y;
                           nearest_cells_to_cell[ID].bc_z[count]=aux_bc_z;
#endif
                           count++;
                         }
                      }
                   }
                }
              }  
          } 

// ********************************************************************
// ********************************************************************
// ********************************************************************
// ********************************************************************
// ********************************************************************
// ********************************************************************
// ********************************************************************
void get_neighbour_cells_cat_analyze(int Nft, int N_cells_back_forth,vector<s_nearest_cells>&nearest_cells_to_cell){

#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif

  ULONG max_neigh_per_dim = 2*N_cells_back_forth+1; // Number of neighbouring cells pr dimension, including the central cell;
  vector<int> index_cells_dime(max_neigh_per_dim,0);
  for(int i=0;i< max_neigh_per_dim ;++i)
     index_cells_dime[i]= i-N_cells_back_forth;  // e.g., 1, 0, -1 if N_cells_back_forth=1 and the central excluded
 

  // total number of cells to explore around a cell, including the central only once as it should be
// ULONG N_Neigh_cells=pow(max_neigh_per_dim+1,2)+ pow(max_neigh_per_dim+1,2)  + pow(max_neigh_per_dim+1,2); 
  ULONG N_Neigh_cells=pow(max_neigh_per_dim,3)-2;  // subtract 2 to avoid cunting 3 times the central cell
  
  ULONG NGRID = Nft*Nft*Nft;
 
 #pragma omp parallel for
  for(ULONG ID=0;ID<NGRID ;++ID)
    { 
      nearest_cells_to_cell[ID].close_cell.resize(N_Neigh_cells,0);
      nearest_cells_to_cell[ID].bc_x.resize(N_Neigh_cells,0);
      nearest_cells_to_cell[ID].bc_y.resize(N_Neigh_cells,0);
      nearest_cells_to_cell[ID].bc_z.resize(N_Neigh_cells,0);
    }

 #pragma omp parallel for collapse(3)
  for(int i=0;i<Nft ;++i)
    for(int j=0;j<Nft ;++j)
      for(int k=0;k<Nft ;++k)
        {
          ULONG ID=index_3d(i,j,k,Nft,Nft);//get ID of cell
          ULONG count=0;
          for(int idx =0; idx < max_neigh_per_dim; ++idx)  // loop in the x-direction
            {
              int new_ni = i - index_cells_dime[idx];
                 int aux_bc_x=0;
                 if(new_ni<0)
                  {
                    new_ni+= Nft;
                     aux_bc_x=-1;
                  }
                 if(new_ni>= Nft)
                  {
                     new_ni-=Nft;
                     aux_bc_x=1;
                  }
               for(int idy =0; idy < max_neigh_per_dim; ++idy) // loop in the y-direction
                 {
                   int new_nj = j - index_cells_dime[idy];
                   int aux_bc_y=0;
                   if(new_nj<0)
                     {
                       new_nj+=Nft;
                       aux_bc_y=-1;
                     }//si la celda esta por detrñas, ponga -1
                    if(new_nj>= Nft)
                      {
                       new_nj-=Nft;
                       aux_bc_y=1;
                      }
                    
                    for(int idz =0; idz < max_neigh_per_dim; ++idz) // loop in the z-direction
                      {
                        int new_nk = k - index_cells_dime[idz];
                          int aux_bc_z=0;
                          if(new_nk<0)
                           {
                             new_nk+=Nft;
                             aux_bc_z=-1;
                           }
                          if(new_nk>=Nft)
                            {
                              new_nk-=Nft;
                              aux_bc_z=1;
                          }
                        
                        ULONG new_id=index_3d(new_ni,new_nj,new_nk,Nft,Nft); // Get the ID of the neighbouring cells;
                        if(new_id!=ID) // This explicitely excludes the ID of the same cell where we are now.
                         {
                          if(count>N_Neigh_cells) cout<<count<<"  "<<N_Neigh_cells<<endl;
                           nearest_cells_to_cell[ID].close_cell[count]=new_id;                        
                           nearest_cells_to_cell[ID].bc_x[count]=aux_bc_x;
                           nearest_cells_to_cell[ID].bc_y[count]=aux_bc_y;
                           nearest_cells_to_cell[ID].bc_z[count]=aux_bc_z;
                           count++;
                         }
                      }
                   }
                }
              }  
          } 

// ********************************************************************
// ********************************************************************
// ********************************************************************
// ********************************************************************


// ********************************************************************
// ********************************************************************
void get_scalar(string FNAME,vector<real_prec>&OUT,ULONG N1,ULONG N2,ULONG N3)
{
  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
  fftw_array<real_prec> dummy(N);
  string fname=FNAME+string(".dat");
#ifdef _FULL_VERBOSE_
  cout<<YELLOW<<"Reading binary file "<<CYAN<<fname<<RESET<<endl;
#endif

  bifstream inStream(fname.data(),file_is_natural);
  assert(inStream.is_open());
  inStream.get(dummy.data,N);
  inStream.close();
  
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)
    OUT[i]=dummy[i];
#ifdef _FULL_VERBOSE_
  cout<<BOLDGREEN<<"                                               ["<<BOLDBLUE<<"DONE"<<BOLDGREEN<<"]"<<RESET<<endl;
#endif

}
// ********************************************************************
// ********************************************************************
// This is the same function write_array. Defined here to be used in the functions inherited from FSK codes
void dump_scalar(const vector<real_prec>&A_rm,ULONG N1,ULONG N2,ULONG N3,int sample_number,string fname)
{
  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
   fftw_array<real_prec> dummy(N);

#pragma omp parallel for
  for(ULONG i=0;i<N;i++)
    dummy[i]=A_rm[i];
  
  string FNAME=fname+string(".dat");
#ifdef _FULL_VERBOSE_
  cout<<YELLOW<<"Writting binary file "<<CYAN<<FNAME<<RESET<<endl;
#endif
  bofstream outStream(FNAME.data(),file_is_natural);
  assert(outStream.is_open());
  outStream.put(dummy.data,N);
  outStream.close();
#ifdef _FULL_VERBOSE_
  cout<<BOLDGREEN<<"                                               ["<<BOLDBLUE<<"DONE"<<BOLDGREEN<<"]"<<RESET<<endl;
#endif
}


// ********************************************************************
string dto_string (double Number)
{
  stringstream ss;
  ss << Number;
  return ss.str();
}

// ********************************************************************
/*
string to_string (real_prec Number)
{
  stringstream ss;
  ss << Number;
  return ss.str();
}
*/
// ********************************************************************
real_prec gsl_integration(gsl_real (*function)(gsl_real, void *) ,void *p,gsl_real LowLimit,gsl_real UpLimit)
{
  gsl_integration_glfixed_table *wf =  gsl_integration_glfixed_table_alloc (500);
  gsl_function F;
  F.params   = p;  
  F.function = function;
  real_prec result=static_cast<real_prec>(gsl_integration_glfixed(&F,LowLimit,UpLimit,wf));
  gsl_integration_glfixed_table_free(wf);
  return  static_cast<real_prec>(result);
}    

// ********************************************************************

real_prec gsl_integration3(int N, gsl_real(*function)(gsl_real, void *) ,void *p, gsl_real LowLimit, gsl_real UpLimit)
{
  gsl_integration_glfixed_table *wf =  gsl_integration_glfixed_table_alloc (N);
  gsl_function F;
  F.params   = p;  
  F.function = function;
  real_prec result=static_cast<real_prec>(gsl_integration_glfixed(&F,LowLimit,UpLimit,wf));
  gsl_integration_glfixed_table_free(wf);
  return static_cast<real_prec>(result);
}    

// ********************************************************************

real_prec gsl_integration2(gsl_real (*function)(gsl_real, void *) ,void *p,vector<gsl_real>XX, vector<gsl_real>WW)
{
  real_prec result=0;
  int nn=WW.size();
  for(int i=0;i<nn;++i)
      result+=WW[i]*function(XX[i],p);
  return static_cast<real_prec>(result);
}    
// ********************************************************************

void gsl_get_GL_weights(gsl_real LowLimit,gsl_real UpLimit, gsl_integration_glfixed_table *wf, vector<gsl_real> &XX, vector<gsl_real>&WW)
{
  gsl_real xi, wi;
  int nn=XX.size();
#ifdef _USE_OMP_
  omp_set_num_threads(omp_get_max_threads());
#endif

  for(int i=0;i<nn;++i){
    gsl_integration_glfixed_point(LowLimit,UpLimit, i, &xi, &wi, wf);
    XX[i]=xi;
    WW[i]=wi;
  }
}    

// ********************************************************************

real_prec gsl_integration_sin_kernel(gsl_real (*function)(gsl_real, void *), void *p, gsl_real omega, gsl_real LowLimit, gsl_real UpLimit){
  int NIT=1e3;
  gsl_real result, error;
  gsl_real L=UpLimit-LowLimit;
  gsl_real relerr=0.1;
  gsl_real abserr=0.01;

  gsl_integration_workspace  *w  =  gsl_integration_workspace_alloc  (NIT); 
  gsl_integration_workspace  *wc =  gsl_integration_workspace_alloc  (NIT); 
  gsl_integration_qawo_table *wf =  gsl_integration_qawo_table_alloc (omega,L,GSL_INTEG_SINE,1000);   

  // La subrutina qawo integra una funcion f peseada con un sin(omega x), 
  // en un intervalo definido
  // con esto no necesito pasarle parametros al integrando pues en este caso f=kP(k)

  gsl_function F;
  F.params   = p;  
  F.function = function; 
  result=0;
  for(;;){
    if(gsl_integration_qawo(&F,LowLimit,abserr,relerr,NIT,w,wf,&result,&error)){
      abserr*=10;
    }
    else{
      break;
    }
  }
  gsl_integration_workspace_free (w);
  gsl_integration_workspace_free (wc);
  gsl_integration_qawo_table_free(wf);
  return static_cast<real_prec>(result);
}    

// ********************************************************************
real_prec gsl_integration_sin_kernel_lowlimit_inf2(gsl_real (*function)(gsl_real, void *), void *p, gsl_real omega, gsl_real LowLimit, gsl_integration_workspace *w,gsl_integration_workspace *wc ){
  gsl_real result, error;
  int NIT=1e3;
  gsl_real L=1;
  gsl_real relerr=0.1;
  gsl_real abserr=0.01;
  gsl_integration_qawo_table *wf =  gsl_integration_qawo_table_alloc (omega,L,GSL_INTEG_SINE,1000);
  // La subrutina qawo integra una funcion f peseada con un sin(omega x),
  // en un intervalo definido
  // con esto no necesito pasarle parametros al integrando pues en este caso f=kP(k)

  gsl_function F;
  F.params   = p;
  F.function = function;
  result=0;
  for(;;){
    if(gsl_integration_qawf(&F,LowLimit,relerr,NIT,w,wc,wf,&result,&error)){
      abserr*=10;
    }
    else{
      break;
    }
  }
  return static_cast<real_prec>(result);
}


// ********************************************************************
// ********************************************************************

real_prec gsl_integration_sin_kernel_lowlimit_inf(gsl_real (*function)(gsl_real, void *), void *p, gsl_real omega, gsl_real LowLimit){
  int NIT=1e3;
  gsl_real result, error;
  gsl_real L=1;
  gsl_real relerr=1.;
  gsl_integration_workspace  *w  =  gsl_integration_workspace_alloc  (NIT);
  gsl_integration_workspace  *wc =  gsl_integration_workspace_alloc  (NIT);
  gsl_integration_qawo_table *wf =  gsl_integration_qawo_table_alloc (omega,L,GSL_INTEG_SINE,NIT);

  // La subrutina qawo integra una funcion f peseada con un sin( omega x), 
  // el cual es en nuestro caso el que viene de j0(x)
  // con esto no necesito pasarle parametros al integrando pues en este caso f=kP(k)
  // Integramos en k para poder usar esta subrutina 
  // utiliza w = 1 con eta=kr, eta lo llamamos k arriba en las funcciones
  
  gsl_function F;
  F.params   = p;  
  F.function = function; 
  result=0;

   for(;;){
    if(gsl_integration_qawf(&F,LowLimit,relerr,NIT,w,wc,wf,&result,&error)){
      relerr*=10;
    }
    else{
      break;
    }
  }

  gsl_integration_workspace_free (w);
  gsl_integration_workspace_free (wc);
  gsl_integration_qawo_table_free(wf);
  return static_cast<real_prec>(result);
}    




// ********************************************************************
// ********************************************************************
real_prec gsl_inter_pointers(real_prec *x, real_prec *y, int n, real_prec xx){
  real_prec ans;
  gsl_real xa[n];
  gsl_real ya[n];
  for(int i=0;i<n;i++)xa[i]=*(x+i);
  for(int i=0;i<n;i++)ya[i]=*(y+i);
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, n);
  gsl_spline_init (spline, xa, ya, n);
  ans=(xx<*x? *x : gsl_spline_eval (spline, xx, acc));
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return ans;
}


// ********************************************************************
// ********************************************************************

void gsl_bspline(vector<gsl_real>&xx, vector<gsl_real>&yy, vector<gsl_real>&new_xx, vector<gsl_real> &new_yy)
{
  int n=xx.size();
  int new_n=new_yy.size();
  gsl_real x_ini=static_cast<gsl_real>(xx[0]);
  gsl_real x_fin=static_cast<gsl_real>(xx[n-1]);
  gsl_matrix *X, *cov;
  gsl_vector *c, *w;
  gsl_vector *x, *y;
  int k = 4;   /*cubic spline stands for k=4*/
  int nbreak=new_n+2-k;
  gsl_bspline_workspace *bw;
  gsl_multifit_linear_workspace *mw;
  gsl_real chisq;
  gsl_vector *B;
  X= gsl_matrix_alloc(n,new_n);
  cov= gsl_matrix_alloc(new_n,new_n);
  x = gsl_vector_alloc(n);
  y = gsl_vector_alloc(n);
  X = gsl_matrix_alloc(n, new_n);
  c = gsl_vector_alloc(new_n);
  w = gsl_vector_alloc(n);

  mw = gsl_multifit_linear_alloc(n, new_n);
  bw=gsl_bspline_alloc(k,nbreak);
  B=gsl_vector_alloc(new_n);

  gsl_bspline_knots_uniform(x_ini, x_fin,bw);

  for(int i=0;i<n;i++)gsl_vector_set(x,i,xx[i]);
  for(int i=0;i<n;i++)gsl_vector_set(y,i,yy[i]);

  
  for(int i=0;i<n;i++){
    gsl_vector_set(w,i,1e7); /**/
    gsl_real xi=gsl_vector_get(x,i);
    gsl_bspline_eval(xi,B,bw);
    for(int j=0;j<new_n;j++){
      gsl_real Bj=gsl_vector_get(B,j);
      gsl_matrix_set(X,i,j,Bj);
    }
  }
  gsl_multifit_wlinear(X,w,y,c,cov,&chisq,mw);
  {
    gsl_real xi, yi, yerr;
    for(int i=0;i<new_n; i++){
      xi=x_ini+i*(x_fin-x_ini)/(new_n-1);
      gsl_bspline_eval(xi, B, bw);
      gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
      new_xx[i]=xi;
      new_yy[i]=yi;
    }
  }
//  int dof = n - new_n;
//  gsl_real tss = gsl_stats_wtss(w->data, 1, y->data, 1, y->size);
//  gsl_real Rsq = 1.0 - chisq / tss;
//  cout<<"Rsq = "<<Rsq<<endl;
//  cout<<"chisq/dof = "<<chisq/dof<<endl;

  gsl_bspline_free(bw);
  gsl_vector_free(B);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_matrix_free(X);
  gsl_vector_free(c);
  gsl_vector_free(w);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(mw);
  return; 
}

// ********************************************************************
// ********************************************************************
real_prec gsl_inter(gsl_real *x, gsl_real *y, int n, gsl_real xn){
  real_prec ans;
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();

  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, n);
  gsl_spline_init (spline, x, y, n);
  ans = gsl_spline_eval (spline, xn, acc);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return  ans; 
}
// ********************************************************************
// ********************************************************************

real_prec gsl_inter_new(const vector<gsl_real> &xx, const vector<gsl_real> &yy, gsl_real xn){
  real_prec ans;
  int n=xx.size();
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, n);
  gsl_spline_init (spline, &xx[0], &yy[0], n);
  ans = gsl_spline_eval (spline, xn, acc);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return  ans; 
}




// ********************************************************************
// ********************************************************************

real_prec gsl_inter_new2(vector<real_prec> &xx, vector<vector<real_prec> > &yy, int li, real_prec xn){
  real_prec ans;
  int n=xx.size();
  gsl_real x[n],y[n];
  for(int i=0;i<n;++i)
  {
    x[i]=static_cast<gsl_real>(xx[i]);
    y[i]=static_cast<gsl_real>(yy[li][i]);
  }
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, n);
  gsl_spline_init (spline, x, y, n);
  ans = gsl_spline_eval (spline, xn, acc);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return  ans; 
}


#ifdef SINGLE_PREC
real_prec gsl_inter_new2(vector<gsl_real> &xx, vector<vector<gsl_real> > &yy, int li, real_prec xn){
  gsl_real ans;
  int n=xx.size();
  gsl_real x[n],y[n];
  for(int i=0;i<n;++i){x[i]=static_cast<gsl_real>(xx[i]);y[i]=static_cast<gsl_real>(yy[li][i]);}
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, n);
  gsl_spline_init (spline, x, y, n);
  ans = gsl_spline_eval (spline, xn, acc);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return  static_cast<real_prec>(ans); 
}
#endif

// ****************************************************************************************************************************
// ****************************************************************************************************************************
void sort_index(int i, int j, int k, int *ii, int *jj, int *kk){
  vector<real_prec> maa;
  maa.push_back(i);
  maa.push_back(j);
  maa.push_back(k);
  sort(maa.begin(), maa.end());
  *ii=maa.at(0);
  *jj=maa.at(1);
  *kk=maa.at(2);
  return ;
}

// ****************************************************************************************************************************
// ****************************************************************************************************************************

void matrix_inversion(vector< vector<real_prec> > &G, vector< vector<real_prec> > &y){

  int n=G[0].size();
  gsl_matrix *A  = gsl_matrix_alloc(n,n);
  gsl_matrix *Ai = gsl_matrix_alloc(n,n);
  for (int i=0;i<n;++i)
      for(int j=0;j<n;++j)
	gsl_matrix_set(A, i,j,static_cast<gsl_real>(G[i][j]));
  gsl_permutation * p = gsl_permutation_alloc(n);
  int signum;
  
  gsl_linalg_LU_decomp(A, p, &signum);
  gsl_linalg_LU_invert(A, p, Ai);
  gsl_permutation_free(p);
  for (int i=0;i<n;++i)for(int j=0;j<n;++j)y[i][j]=gsl_matrix_get(Ai, i,j);
  gsl_matrix_free (A);
  gsl_matrix_free (Ai);    
  return ;
}

// ****************************************************************************************************************************
// ****************************************************************************************************************************

void matrix_det(vector< vector<real_prec> > &G, real_prec &det){
  int n=G[0].size();
  gsl_matrix *A  = gsl_matrix_alloc(n,n);
  for (int i=0;i<n;++i)
      for(int j=0;j<n;++j)
          gsl_matrix_set(A, i,j,G[i][j]);
  gsl_permutation * p = gsl_permutation_alloc(n);
  int signum;
  gsl_linalg_LU_decomp(A, p, &signum);
  det=gsl_linalg_LU_det(A, signum);
  gsl_permutation_free(p);
  gsl_matrix_free (A);
}


// ************************************************************************************
void get_eigen(vector<vector<real_prec>>&icov_masses, vector<real_prec>&masses){
  gsl_matrix *icova;
  gsl_vector *eig;
  gsl_eigen_symm_workspace *workspace;
  int n_par=masses.size();
  icova = gsl_matrix_alloc(n_par,n_par);
  eig = gsl_vector_alloc(n_par);
  for(int i=0;i<n_par;i++)
      for(int j=0;j<n_par;j++)
          gsl_matrix_set(icova, i,j, icov_masses[i][j]);
  workspace = gsl_eigen_symm_alloc(n_par);
  gsl_eigen_symm(icova, eig, workspace);
  gsl_sort_vector(eig);
  gsl_eigen_symm_free(workspace);
  for(int ip=0;ip<n_par;ip++)masses[ip]=sqrt(gsl_vector_get(eig,ip));
}


// ************************************************************************************
void get_det_matrix(vector<vector<real_prec>>&matriz, real_prec &determinant){
  vector<real_prec>det(matriz.size(),0);
  get_eigen(matriz, det);
  determinant=1;         ;//Para escalar
  for(int i=0;i<matriz.size();i++)
     determinant*=det[i];
}

//##################################################################################

real_prec real_sh(int l, int m, real_prec theta, real_prec phi){
  return  gsl_sf_legendre_sphPlm(l,m,cos(theta))*cos(m*phi);
}

//##################################################################################

real_prec imag_sh(int l, int m, real_prec theta, real_prec phi){
  return gsl_sf_legendre_sphPlm(l,m,cos(theta))*sin(m*phi);
}

//##################################################################################

real_prec bessel(real_prec x, int n){
  return gsl_sf_bessel_jl(n,x);
}


//##################################################################################

real_prec dbessel(real_prec x, int n){
  return  (n*bessel(x,n-1)-(n+1)*bessel(x,n+1))/(2.*n+1);
}

//##################################################################################


real_prec ddbessel(real_prec x, int n){
  return  -(2/x)*dbessel(x,n)-bessel(x,n)*(1.0-n*(n+1)/pow(x,2));
}



// sort (Quicksort) by reindexing adopted from Numerical Recipes in C++
//##################################################################################

inline void SWAP_L(int &a, int &b) {
  int dum=a;
  a=b;
  b=dum;
}

inline void SWAP_LU(ULONG &a, ULONG &b) {
  ULONG dum=a;
  a=b;
  b=dum;
}


//##################################################################################
void indexx(vector<int>&arr, vector<int>&indx)
{
  const int M=7,NSTACK=50;

  int n = indx.size();
  int i,indxt,ir,j,k,jstack=-1,l=0;
  int a;
  int istack[NSTACK];
  
  ir=n-1;
  for (j=0;j<n;j++) indx[j]=j;
  for (;;) {
    if (ir-l < M)
     {
       for (j=l+1;j<=ir;j++) 
       {
         indxt=indx[j];
	       a=arr[indxt];
	       for (i=j-1;i>=l;i--) 
          {
	          if (arr[indx[i]] <= a) break;
       	  indx[i+1]=indx[i];
          }
	      indx[i+1]=indxt;
       }
      if (jstack < 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
      } 
      else 
      {
       k=(l+ir) >> 1;
       SWAP_L(indx[k],indx[l+1]);
        if (arr[indx[l]] > arr[indx[ir]]) {
        	SWAP_L(indx[l],indx[ir]);
        }
      if (arr[indx[l+1]] > arr[indx[ir]]) {
	      SWAP_L(indx[l+1],indx[ir]);
      }
      if (arr[indx[l]] > arr[indx[l+1]]) {
       	SWAP_L(indx[l],indx[l+1]);
      }
      i=l+1;
      j=ir;
      indxt=indx[l+1];
      a=arr[indxt];
      for (;;) {
    	 do i++; while (arr[indx[i]] < a);
	     do j--; while (arr[indx[j]] > a);
	     if (j < i) break;
	       SWAP_L(indx[i],indx[j]);
      }
      indx[l+1]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack >= NSTACK) {
	       cerr << "NSTACK too small in indexx." << endl;
      	exit(1);
      }
      if (ir-i+1 >= j-l) {
	      istack[jstack]=ir;
	      istack[jstack-1]=i;
      	ir=j-1;
      } 
      else {
	     istack[jstack]=j-1;
       istack[jstack-1]=l;
    	l=i;
     }
    }
  }
}


//##################################################################################
void indexx_ulong(vector<ULONG>&arr, vector<ULONG>&indx)
{
  const int M=7,NSTACK=50;

  ULONG n = indx.size();
  ULONG i,indxt,ir,j,k,jstack=-1,l=0;
  ULONG a;
  int istack[NSTACK];
  
  ir=n-1;
  for (j=0;j<n;j++) indx[j]=j;
  for (;;) {
    if (ir-l < M)
     {
       for (j=l+1;j<=ir;j++) 
       {
         indxt=indx[j];
         a=arr[indxt];
         for (i=j-1;i>=l;i--) 
          {
            if (arr[indx[i]] <= a) break;
          indx[i+1]=indx[i];
          }
        indx[i+1]=indxt;
       }
      if (jstack < 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
      } 
      else 
      {
       k=(l+ir) >> 1;
       SWAP_LU(indx[k],indx[l+1]);
        if (arr[indx[l]] > arr[indx[ir]]) {
          SWAP_LU(indx[l],indx[ir]);
        }
      if (arr[indx[l+1]] > arr[indx[ir]]) {
        SWAP_LU(indx[l+1],indx[ir]);
      }
      if (arr[indx[l]] > arr[indx[l+1]]) {
        SWAP_LU(indx[l],indx[l+1]);
      }
      i=l+1;
      j=ir;
      indxt=indx[l+1];
      a=arr[indxt];
      for (;;) {
       do i++; while (arr[indx[i]] < a);
       do j--; while (arr[indx[j]] > a);
       if (j < i) break;
         SWAP_LU(indx[i],indx[j]);
      }
      indx[l+1]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack >= NSTACK) {
         cerr << "NSTACK too small in indexx." << endl;
        exit(1);
      }
      if (ir-i+1 >= j-l) {
        istack[jstack]=ir;
        istack[jstack-1]=i;
        ir=j-1;
      } 
      else {
       istack[jstack]=j-1;
       istack[jstack-1]=l;
      l=i;
     }
    }
  }
}



//##################################################################################
void sort_vectors(vector<ULONG>& v1,vector<ULONG>&v2, vector<ULONG>&v3, vector<ULONG>& v4, vector<ULONG>&v5, vector<ULONG>& v6, vector<ULONG>& v7, vector<real_prec>& v8 ){
  // v1 is to be sorted. The elements of the other vectors
  // are shuffled accordingly.

  
  unsigned long j;
  long n=v1.size();
  vector<ULONG>iwksp(n,0);
  vector<float>wksp(n,0);
  
  indexx_ulong(v1,iwksp);
  for (j=0;j<n;++j) wksp[j]=v1[j];
  for (j=0;j<n;++j) v1[j]=wksp[iwksp[j]];

  for (j=0;j<n;++j) wksp[j]=v2[j];
  for (j=0;j<n;++j) v2[j]=wksp[iwksp[j]];

  for (j=0;j<n;++j) wksp[j]=v3[j];
  for (j=0;j<n;++j) v3[j]=wksp[iwksp[j]];

  for (j=0;j<n;++j) wksp[j]=v4[j];
  for (j=0;j<n;++j) v4[j]=wksp[iwksp[j]];

  for (j=0;j<n;++j) wksp[j]=v5[j];
  for (j=0;j<n;++j) v5[j]=wksp[iwksp[j]];

  for (j=0;j<n;++j) wksp[j]=v6[j];
  for (j=0;j<n;++j) v6[j]=wksp[iwksp[j]];

  for (j=0;j<n;++j) wksp[j]=v7[j];
  for (j=0;j<n;++j) v7[j]=wksp[iwksp[j]];


}


//##################################################################################

size_t m_getBuffer(istream& inStream) {
  
    // reset buffer
    inStream.clear();
    inStream.seekg(0, inStream.beg);

    string line = "";

    // get line from stream until is not a comment or empty
    while (false == inStream.eof() && ((true == line.empty()) || ('#' == line[0]))) {
      getline(inStream, line);
    }

    // even number of lines to read as 3MB buffer
    size_t numLines = 3*1024*1024/line.size();
    numLines += numLines&0x1;


    // reset buffer
    inStream.clear();
    inStream.seekg(0, inStream.beg);

    return numLines;
  }

//##################################################################################

size_t m_countLines(istream& inStream) {
  
  // reset stream to beginning and clear status
  inStream.clear();
  inStream.seekg(0, inStream.beg);
  
  // counter of lines
  size_t counter = -1;
  
  // line to read
  string line = "";
  
  while (false == inStream.eof()) {
    getline(inStream, line);
    counter++;
  }
  
  // reset buffer
  inStream.clear();
  inStream.seekg(0, inStream.beg);
  
  return counter;
}


//##################################################################################
//##################################################################################

real_prec get_mean(const vector<real_prec>& ini)
{
#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif
  double mean=0.;
#pragma opmp parallel for reduction(+:mean)
  for(ULONG i=0;i<ini.size();++i)
    mean+=static_cast<double>(ini[i]);
  mean/=static_cast<double>(ini.size());
  return static_cast<real_prec>(mean);
}



//##################################################################################
//##################################################################################
//##################################################################################

real_prec get_var(const vector<real_prec> & ini)
{
#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif

  real_prec mean=get_mean(ini);
  real_prec vari=0.;
#ifdef _USE_OMP_
#pragma opmp parallel for reduction(+:vari)
#endif
  for(ULONG i=0;i<ini.size();++i)
    vari+=pow(ini[i]-mean,2);
  vari=vari/(static_cast<double>(ini.size())-1.0);
  return vari;
}



//##################################################################################

real_prec get_var(real_prec mean, const vector<real_prec> & ini)
{
#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif

  real_prec vari=0.;
#ifdef _USE_OMP_
#pragma opmp parallel for reduction(+:vari)
#endif
  for(ULONG i=0;i<ini.size();++i)
    vari+=pow(ini[i]-mean,2);
  vari=vari/(static_cast<double>(ini.size())-1.0);
  return vari;
}

//##################################################################################
//##################################################################################
//##################################################################################

void log_smooth(vector<real_prec>&xin, vector<real_prec>&vin)
{

  // **************************************************************************************
  // steps:
  // i ) interpolate the raw ratio in log bins in k
  // ii) smooth the result
  // iii) interpolate the smoothed version and update power_ratio
  
  // Original modes, converted to log
  vector<gsl_real>lkvec(vin.size(),0);
  vector<gsl_real>vin_new(vin.size(),0);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i = 0;i < vin.size(); ++i)
    lkvec[i]=static_cast<gsl_real>(log10(xin[i]));
  for(ULONG i = 0;i < vin.size(); ++i)
    vin_new[i]=static_cast<gsl_real>(vin[i]);

  int n_smooth=vin.size()*2;

  // Begin smooth ker_aux in log k
  // Log bins, in the same range as the original modes
  
  real_prec deltalk= static_cast<real_prec>(log10(xin[xin.size()-1.0]/xin[0])/static_cast<real_prec>(n_smooth-1.0));
  vector<gsl_real>aux_logk(n_smooth,0);
  
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i = 0; i < n_smooth-1; ++i)
    aux_logk[i]=log10(xin[0])+i*deltalk;
  aux_logk[n_smooth-1]=lkvec[xin.size()-1];
  
  // Interpolate original ratio in log k
  vector<gsl_real>aux_oratio(n_smooth,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i = 0;i < n_smooth; ++i)
    aux_oratio[i]=gsl_inter_new(lkvec,vin_new, aux_logk[i]);
  
  // Bspline smooth ratio, already interpolated in log
  int n_smooth_new=  static_cast<int>(floor(vin.size()/2));    //60;
  vector<gsl_real>aux_kvec(n_smooth_new, 0);
  vector<gsl_real>aux_ratio(n_smooth_new, 0);

  gsl_bspline(aux_logk, aux_oratio,  aux_kvec, aux_ratio);

  // Reassign the smoothed version to Kernel_bins via interpolation
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i = 0;i < vin.size(); ++i)
    vin[i]=gsl_inter_new(aux_kvec, aux_ratio, lkvec[i]);

}

//##################################################################################
//##################################################################################
//##################################################################################

void lin_smooth(vector<real_prec>&xin, vector<real_prec>&vin, int nr)
{

  // **************************************************************************************
  // steps:
  // i ) interpolate the raw ratio in log bins in k
  // ii) smooth the result
  // iii) interpolate the smoothed version and update power_ratio

  // Original modes, converted to log
  vector<gsl_real>lkvec(vin.size(),0);
  vector<gsl_real>vin_new(vin.size(),0);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i = 0;i < xin.size(); ++i)
    lkvec[i]=static_cast<gsl_real>(xin[i]);

  for(ULONG i = 0;i < vin.size(); ++i)
    vin_new[i]=static_cast<gsl_real>(vin[i]);



  /*
  int n_smooth=static_cast<int>(floor(vin.size()/3));
  // Begin smooth ker_aux in log k
  // Log bins, in the same range as the original modes

  real_prec deltalk= (xin[xin.size()-1]-xin[0])/static_cast<real_prec>(n_smooth-1);
  vector<gsl_real>aux_logk(n_smooth,xin[0]);
  aux_logk[n_smooth-1]=lkvec[xin.size()-1];

#pragma omp parallel for
  for(ULONG i = 1; i < n_smooth; ++i)
      aux_logk[i]=xin[0]+static_cast<real_prec>(i)*deltalk;

  // Interpolate original ratio in log k
  vector<gsl_real>aux_oratio(n_smooth,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i = 0;i < n_smooth; ++i)
    aux_oratio[i]=gsl_inter_new(lkvec,vin_new, aux_logk[i]);

  // Bspline smooth ratio, already interpolated in log
*/

  int n_smooth_new= static_cast<int>(floor(vin.size())/nr);
  vector<gsl_real>aux_kvec(n_smooth_new, 0);
  vector<gsl_real>aux_ratio(n_smooth_new, 0);

//  gsl_bspline(aux_logk, aux_oratio,  aux_kvec, aux_ratio);


 gsl_bspline(lkvec, vin_new,  aux_kvec, aux_ratio);

  // Reassign the smoothed version to Kernel_bins via interpolation
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i = 0;i < vin.size(); ++i)
    vin[i]=gsl_inter_new(aux_kvec, aux_ratio, lkvec[i]);

}
//##################################################################################
//##################################################################################
//##################################################################################

void gradfindif(ULONG N1,ULONG N2,ULONG N3,real_prec L1,real_prec L2,real_prec L3,vector<real_prec>in,vector<real_prec>&out,int dim)
{
    // Gradient using finite differences
#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif
    real_prec factor=static_cast<real_prec>(N1)/static_cast<real_prec>(num_2*L1);

  int N1i=static_cast<int>(N1);
  int N2i=static_cast<int>(N2);
  int N3i=static_cast<int>(N3);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int x = 0; x < N1i; x++)
    for(int y = 0; y < N2i; y++)
      for(int z = 0; z < N3i; z++)
        {
          int xr, xrr, xl, xll;
          int yr, yrr, yl, yll;
          int zr, zrr, zl, zll;

          xrr = xll = xr = xl = x;
          yrr = yll = yr = yl = y;
          zrr = zll = zr = zl = z;

          switch (dim)
            {
            case 1:
              {
                xr = x + 1;
                xl = x - 1;
                xrr = x + 2;
                xll = x - 2;
                if(xr >= N1)
                  xr -= N1;
                if(xrr >= N1)
                  xrr -= N1;
                if(xl < 0)
                  xl += N1;
                if(xll < 0)
                  xll += N1;
                break;
              }
            case 2:
              {
                yr = y + 1;
                yl = y - 1;
                yrr = y + 2;
                yll = y - 2;
                if(yr >= N2)
                  yr -= N2;
                if(yrr >= N2)
                  yrr -= N2;
                if(yl < 0)
                  yl += N2;
                if(yll < 0)
                  yll += N2;
                break;
              }
            case 3:
              {
                zr = z + 1;
                zl = z - 1;
                zrr = z + 2;
                zll = z - 2;
                if(zr >= N3)
                  zr -= N3;
                if(zrr >= N3)
                  zrr -= N3;
                if(zl < 0)
                  zl += N3;
                if(zll < 0)
                  zll += N3;
                break;
              }
            }


          ULONG iinl = index_3d(xl,yl,zl,N2,N3);
          ULONG iinr = index_3d(xr,yr,zr,N2,N3);
          ULONG iinll= index_3d(xll,yll,zll,N2,N3);
          ULONG iinrr= index_3d(xrr,yrr,zrr,N2,N3);

          ULONG iout = index_3d(x,y,z,N2,N3);
          out[iout] = -static_cast<real_prec>((factor*((4.0/3.0)*(in[iinl]-in[iinr])-(1.0/6.0)*(in[iinll]-in[iinrr]))));
        }
}


//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################

// This function expects delta.
// for type= log, it uses log10(1+delta)
// for type= lin, it uses delta
// It returns delta again

void calc_pdf(string type, ULONG N, ULONG Nk, real_prec maxk, real_prec mink, const vector<real_prec>&in, vector<real_prec>&pdfin)
{
#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif
    real_prec dk=(maxk-mink)/static_cast<double>(Nk);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<pdfin.size();++i)
    pdfin[i]=0;

  ULONG ntot=0;

  if(type=="log")
    {
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ntot)
#endif
      for(ULONG i=0;i<N;i++)
       {
        real_prec den=log10(NUM_IN_LOG+static_cast<real_prec>(in[i]));
        if (den>=mink && den<=maxk)
          {
            ULONG ik=static_cast<ULONG>(floor((den-mink)/dk));
            if(ik==Nk)ik--;
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
            pdfin[ik]+=1.0;
           ntot++;
          }
       }
    }
  else
    {
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ntot)
#endif
    for(ULONG i=0;i<N;i++)
      {
        real_prec den=in[i];
        if(den>=mink && den<=maxk)
          {
            ULONG ik=static_cast<ULONG>(floor((den-mink)/dk));
            if(ik==Nk)ik--;
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
            pdfin[ik]+=1.0;
            ntot++;
          }
      }
 }

 for(ULONG i=0;i<Nk;i++)
    pdfin[i]/=static_cast<double>(ntot);



}


//##################################################################################
//##################################################################################
//Input must be delta, and tjhe limits must be max and mins in log(1+detla)
// The pdfs must be P(log(1+delta))
void rankorder(int seed, vector<real_prec>dens_bins, ULONG Nk, real_prec maxk, real_prec mink, vector<real_prec>&in, vector<real_prec>&pdfin, vector<real_prec>&pdfout)
{
#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif

  ULONG N=in.size();
  double dk=(maxk-mink)/static_cast<double>(Nk);

  for (ULONG i=0;i<N;i++)
    in[i]=log10(NUM_IN_LOG+static_cast<real_prec>(in[i]));

  const gsl_rng_type *  Ta;
  gsl_rng * ra ;
  gsl_rng_env_setup();
  gsl_rng_default_seed=seed;
  Ta = gsl_rng_ranlux;
  ra = gsl_rng_alloc (Ta);

  for (ULONG i=0;i<N;i++)
    {
/*
      ULONG ik=0;
      if (den>=mink)
        ik=static_cast<ULONG>(floor((den-mink)/dk));
       if(ik==Nk)ik--;
*/
      ULONG ik=get_bin(in[i],mink,Nk,dk,true);

      real_prec pdfdi=0.;
      for (ULONG j=0;j<ik+1;j++)
        pdfdi+=pdfin[j];

      real_prec pdfdo=0.;
      ULONG kk=0;
      while (pdfdo<pdfdi)
        {
          pdfdo+=pdfout[kk];
          kk++;
        }

      ULONG l=kk;
      if (kk>0)
        l-=1; 
      if (kk>Nk)
        l=Nk; 

      /* real_prec dr=static_cast<real_prec>(2.*gsl_rng_uniform(ra)*dk); */
      /* real_prec denB=dens[l]+dr-  2.*dk;  // ?? repartirlo aleatoriamente en este bin. OJO ACA VERIFICARLOS */
      // this depends on the way the vector dens has been designed.

      //* / I defined it to lie a the center of the bin */
      // El factor 10 permite que la pdf final sea más suave
      real_prec xr= static_cast<real_prec>(gsl_rng_uniform(ra));
      real_prec denB=dens_bins[l]+ (xr-0.5)*dk;  //s ?? repartirlo aleatoriamente en este bin. OJO ACA VERIFICARLOS */
      in[i]=pow(10.0, denB)-NUM_IN_LOG; // Convert to the format of the input density field
    }

}



//##################################################################################
//##################################################################################
void EigenValuesTweb(ULONG Nft, real_prec L1, const vector<real_prec> &delta, const vector<real_prec> &phi, vector<real_prec> &out1, vector<real_prec> &out2, vector<real_prec> &out3)
{
    // Function from webclass by F. Kitaura. File webclass.cc
    // The potential vctor is the gravitational potential obtained form the Poisson solver
    real_prec L2=L1;
    real_prec L3=L2;

#ifdef _USE_OMP_
    int NTHREADS = omp_get_max_threads();
    omp_set_num_threads(NTHREADS);
#endif

    ULONG NGRID = static_cast<ULONG>(Nft*Nft*Nft);
    ULONG Nhalf=static_cast<ULONG>(Nft)*static_cast<ULONG>(Nft)*static_cast<ULONG>(Nft/2+1);


    vector<real_prec> LapPhivx(NGRID,0), LapPhivy(NGRID,0), LapPhivz(NGRID,0);
    vector<real_prec> LapPhivxy(NGRID,0), LapPhivxz(NGRID,0), LapPhivyz(NGRID,0);
    vector<real_prec> LapPhivzx(NGRID,0), LapPhivzy(NGRID,0), LapPhivyx(NGRID,0);

   
    vector<real_prec> dummy(NGRID,0);


#ifdef _USE_GFINDIFF_EIGENV_
    // Get gradient in x-direction
    gradfindif(Nft,Nft,Nft,L1,L2,L3,phi,dummy,1);
    // Get second derivatives of grad_x
    gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivx,1);
    gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivxy,2);
    gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivxz,3);

    // Get gradient in y-direction
    gradfindif(Nft,Nft,Nft,L1,L2,L3,phi,dummy,2);
    // Get second derivatives of grad_y
    gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivy,2);
    gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivyx,1);
    gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivyz,3);

    // Get gradient in z-direction
    gradfindif(Nft,Nft,Nft,L1,L2,L3,phi,dummy,3);
    // Get second derivatives of grad_z
    gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivz,3);
    gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivzx,1);
    gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivzy,2);
#elif defined _USE_GFFT_EIGENV_
    complex_prec *philv= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<Nhalf;i++)
      {
        philv[i][REAL]=0.0;
        philv[i][IMAG]=0.0;
      }

    do_fftw_r2c(Nft,phi,philv);  // Get Phi(k)

    // Get gradient in x-direction
    calc_LapPhiv(Nft,Nft,Nft,L1,L2,L3,philv,dummy,1,1); // from Phi(k) get \partial Phi(r)/\partial x
    // Get second derivatives of grad_x

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<Nhalf;i++)
    {
      philv[i][REAL]=0.0;
      philv[i][IMAG]=0.0;
    }
  do_fftw_r2c(Nft,dummy,philv); // Get Phi'(k)

  calc_LapPhiv(Nft,Nft,Nft,L1,L2,L3,philv,LapPhivx,1,1);  // From Phi'(k)_x get \partial²Phi/\partial x²
  calc_LapPhiv(Nft,Nft,Nft,L1,L2,L3,philv,LapPhivxy,1,2); // From Phi'(k)_x get \partial²Phi/\partial y \partial x
  calc_LapPhiv(Nft,Nft,Nft,L1,L2,L3,philv,LapPhivxz,1,3); // From Phi'(k)_x get \partial²Phi/\partial z \partial x

    // Get gradient in y-direction
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<Nhalf;i++)
    {
      philv[i][REAL]=0.0;
      philv[i][IMAG]=0.0;
    }
    do_fftw_r2c(Nft,phi,philv);  // Get Phi(k)

    calc_LapPhiv(Nft,Nft,Nft,L1,L2,L3,philv,dummy,2,2); // from Phi(k) get \partial Phi(r)/\partial y

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<Nhalf;i++)
    {
      philv[i][REAL]=0.0;
      philv[i][IMAG]=0.0;
    }
    do_fftw_r2c(Nft,dummy,philv); //Get Phi'(k)
    // Get second derivatives of grad_y
    calc_LapPhiv(Nft,Nft,Nft,L1,L2,L3,philv,LapPhivy,2,2);  // From Phi'(k)_x get \partial²Phi/\partial y²
    calc_LapPhiv(Nft,Nft,Nft,L1,L2,L3,philv,LapPhivyx,2,1);
    calc_LapPhiv(Nft,Nft,Nft,L1,L2,L3,philv,LapPhivyz,2,3);

    // Get gradient in y-direction
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<Nhalf;i++)
    {
      philv[i][REAL]=0.0;
      philv[i][IMAG]=0.0;
    }
    do_fftw_r2c(Nft,phi,philv);

    calc_LapPhiv(Nft,Nft,Nft,L1,L2,L3,philv,dummy,3,3);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<Nhalf;i++)
    {
      philv[i][REAL]=0.0;
      philv[i][IMAG]=0.0;
    }
    do_fftw_r2c(Nft,dummy,philv);
    // Get second derivatives of grad_z
    calc_LapPhiv(Nft,Nft,Nft,L1,L2,L3,philv,LapPhivz,3,3);
    calc_LapPhiv(Nft,Nft,Nft,L1,L2,L3,philv,LapPhivzx,3,1);
    calc_LapPhiv(Nft,Nft,Nft,L1,L2,L3,philv,LapPhivzy,3,2);

       fftw_free(philv);
#endif


    
    dummy.clear();
    dummy.shrink_to_fit();


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0; i < NGRID;++i)
      {
        // M = [ A B C,
        //       D E F,
        //       G H I ]

        // Eqn = aλ³ + bλ² + cλ + d = 0
        // =-λ³+λ²(A+E+I)+λ(DB+GC+FH-AE-AI-EI)+(AEI-AFH-DBI+DCH+GBF-GCE)
        // a=-1
        // b=tr(M)
        // c=1/2[tr(M*M)-tr(M)*tr(M)]
        // d=det(M)

        real_prec A=LapPhivx[i];
        real_prec B=LapPhivxy[i];
        real_prec C=LapPhivxz[i];

        real_prec D=LapPhivyx[i];
        real_prec E=LapPhivy[i];
        real_prec F=LapPhivyz[i];

        real_prec G=LapPhivzx[i];
        real_prec H=LapPhivzy[i];
        real_prec I=LapPhivz[i];

        real_prec trM=A+E+I;
        real_prec detM=A*E*I-A*F*H-D*B*I+D*C*H+G*B*F-G*C*E;

        real_prec trM2=D*B+G*C+F*H-A*E-A*I-E*I;

        // Eqn = a1λ³ + a2λ² + a3λ + a4 = 0
        real_prec a1= -num_1;
        real_prec a2= trM;
        real_prec a3= trM2;
        real_prec a4= detM;

        // Eqn = λ³ + aλ² + bλ + c = 0
        gsl_real a= static_cast<gsl_real>(a2/a1);
        gsl_real b= static_cast<gsl_real>(a3/a1);
        gsl_real c= static_cast<gsl_real>(a4/a1);


        gsl_real x0, x1, x2;
        gsl_poly_solve_cubic (a,  b,  c, &x0, &x1, &x2);

        real_prec Eig1 = x0;
        real_prec Eig2 = x1;
        real_prec Eig3 = x2;


        if (Eig1 >= Eig2)
          {
            if (Eig1>=Eig3)
              {
                out1[i]=Eig1;
                if (Eig2>=Eig3)
                  {
                    out2[i]=Eig2;
                    out3[i]=Eig3;
                  }
                else
                  {
                    out2[i]=Eig3;
                    out3[i]=Eig2;
                  }
              }
            else
              {
                out1[i]=Eig3;
                if (Eig1>=Eig2)
                  {
                    out2[i]=Eig1;
                    out3[i]=Eig2;
                  }
                else
                  {
                    out2[i]=Eig2;
                    out3[i]=Eig1;
                  }
              }
          }
        else
          {
            if (Eig2>=Eig3)
              {
                out1[i]=Eig2;
                if (Eig1>=Eig3)
                  {
                    out2[i]=Eig1;
                    out3[i]=Eig3;
                  }
                else
                  {
                    out2[i]=Eig3;
                    out3[i]=Eig1;
                  }
              }
            else
              {
                out1[i]=Eig3;
                if (Eig1>=Eig2)
                  {
                    out2[i]=Eig1;
                    out3[i]=Eig2;
                  }
                else
                  {
                    out2[i]=Eig2;
                    out3[i]=Eig1;
                  }
              }
          }
      }


}

//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
void downsampling(ULONG Nft_HR, ULONG Nft_LR, int imas, vector<real_prec>&HR_field, vector<real_prec>&LR_field, real_prec Lbox)
{

#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif


  ULONG NTT_HR=static_cast<ULONG>(Nft_HR*Nft_HR*(Nft_HR/2+1));
  ULONG NTT_LR=static_cast<ULONG>(Nft_LR*Nft_LR*(Nft_LR/2+1));
  ULONG NGRID_LR=LR_field.size();

  if(NTT_HR<NTT_LR)
    {
      cout<<RED<<"Error. High resolution mesh must have larger number of Nft cells than low resolution."<<endl;
      exit(0);
    }

#ifdef _DOUBLE_PREC_
    complex_prec *HR_FOURIER= (complex_prec *)fftw_malloc(2*NTT_HR*sizeof(real_prec));
#else
  complex_prec *HR_FOURIER= (complex_prec *)fftwf_malloc(2*NTT_HR*sizeof(real_prec));
#endif

#pragma omp parallel for
   for(ULONG i=0;i<NTT_HR;++i)HR_FOURIER[i][REAL]=0;
#pragma omp parallel for
   for(ULONG i=0;i<NTT_HR;++i)HR_FOURIER[i][IMAG]=0;


   // DO FFT of the high-res density field
   do_fftw_r2c(Nft_HR,HR_field,HR_FOURIER);



    // Get the MAS correction for the H-Res density field
   vector<real_prec> correction(Nft_HR,1.0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int i = 1 ; i < Nft_HR; ++i )
    {
      int  coords= i<=Nft_HR/2? i: i-Nft_HR;
      real_prec xx=static_cast<real_prec>(coords)*M_PI/static_cast<real_prec>(Nft_HR);
      real_prec kernel= sin(xx)/static_cast<real_prec>(xx);
      correction[i]= pow(kernel,imas+1);
    }

#ifdef _DOUBLE_PREC_
  complex_prec *LR_FOURIER= (complex_prec *)fftw_malloc(2*NTT_LR*sizeof(real_prec));
#else
   complex_prec *LR_FOURIER= (complex_prec *)fftwf_malloc(2*NTT_LR*sizeof(real_prec));
#endif



#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<NTT_LR;++i)LR_FOURIER[i][REAL]=0;
#ifdef _USE_OMP_
#pragma omp parallel for
  for(ULONG i=0;i<NTT_LR;++i)LR_FOURIER[i][IMAG]=0;
#endif

  vector<real_prec> kmodes_lr(Nft_LR,0.0);
  real_prec delta_k= 2*M_PI/Lbox;

#pragma omp parallel for
   for(int i = 0 ; i < Nft_LR; ++i )
   {
      int coords= i < Nft_LR/2+1? i: i-Nft_LR;
      kmodes_lr[i]=static_cast<real_prec>(coords)*delta_k;
   }


   // Get the FT of the Low-res density field
real_prec delta_box_hr= Lbox/static_cast<real_prec>(Nft_HR);
real_prec alpha=static_cast<double>(Nft_HR)/static_cast<double>(Nft_LR);
real_prec delta_shift= 0.5*(alpha-1.0)*delta_box_hr; //Yu

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<Nft_LR;++i)
    {
      ULONG i_hr = i < Nft_LR/2 +1 ? i : i+Nft_HR-Nft_LR;
      for(ULONG j=0;j< Nft_LR;++j)
        {
          ULONG j_hr = j < Nft_LR/2 +1 ? j : j+Nft_HR-Nft_LR;
          for(ULONG k=0;k<Nft_LR/2+1;++k)
            {
              ULONG k_hr=  k < Nft_LR/2 +1 ? k : k+Nft_HR-Nft_LR;
              ULONG index_lr=index_3d(i,   j,   k,   Nft_LR, Nft_LR/2+1);
              ULONG index_hr=index_3d(i_hr,j_hr,k_hr,Nft_HR ,Nft_HR/2+1);
              real_prec corr=correction[i_hr]*correction[j_hr]*correction[k_hr];
              real_prec new_real=HR_FOURIER[index_hr][REAL]/corr;
              real_prec new_imag=HR_FOURIER[index_hr][IMAG]/corr;
              real_prec shift=delta_shift*(kmodes_lr[i] + kmodes_lr[j] + kmodes_lr[k]);
              LR_FOURIER[index_lr][REAL] = (cos(shift)*new_real - sin(shift)*new_imag); //This follows the convention delta_2-> exp(iks)*delta_2
              LR_FOURIER[index_lr][IMAG] = (cos(shift)*new_imag + sin(shift)*new_real);
            }
        }
    }


  // Set Ny freq of the FT of the L-res density field to real:
  real_prec a, b;
  int ii, jj, kk;
  ULONG ind;
  real_prec fac=1.0;

  ii=0; jj=0; kk=0;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;

  ii=0; jj=0; kk=Nft_LR/2;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;

  ii=0; jj=Nft_LR/2; kk=0;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;

  ii=0; jj=Nft_LR/2; kk=Nft_LR/2;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;

  ii=Nft_LR/2; jj=0; kk=0;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;

  ii=Nft_LR/2; jj=0; kk=Nft_LR/2;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;

  ii=Nft_LR/2; jj=Nft_LR/2; kk=0;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;

  ii=Nft_LR/2; jj=Nft_LR/2; kk=Nft_LR/2;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;


 // Transform to L-res density field.
  do_fftw_c2r(Nft_LR ,LR_FOURIER,LR_field);


#ifdef _DOUBLE_PREC_
  fftw_free(LR_FOURIER);
  fftw_free(HR_FOURIER);
#else
  fftwf_free(LR_FOURIER);
  fftwf_free(HR_FOURIER);
#endif


}



//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
/*
delta1 = Aexp(i phi1)
delta2 = Bexp(i phi2)
output delta2=A exp(i phi2), that is, conserves the phases but with a different amplitude
*/


void swap_amp_fourier(ULONG Nft, vector<real_prec>&in_ref, vector<real_prec>&out_field)
{

#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif

  ULONG NTT=static_cast<ULONG>(Nft*Nft*(Nft/2+1));
  ULONG NGRID=in_ref.size();

#ifdef _DOUBLE_PREC_
  complex_prec *in_FOURIER= (complex_prec *)fftw_malloc(2*NTT*sizeof(real_prec));
  complex_prec *out_FOURIER= (complex_prec *)fftw_malloc(2*NTT*sizeof(real_prec));
#else
  complex_prec *in_FOURIER= (complex_prec *)fftwf_malloc(2*NTT*sizeof(real_prec));
  complex_prec *out_FOURIER= (complex_prec *)fftwf_malloc(2*NTT*sizeof(real_prec));
#endif

#pragma omp parallel for
  for(ULONG i=0;i<NTT;++i)in_FOURIER[i][REAL]=0;
#pragma omp parallel for
  for(ULONG i=0;i<NTT;++i)in_FOURIER[i][IMAG]=0;
#pragma omp parallel for
  for(ULONG i=0;i<NTT;++i)out_FOURIER[i][REAL]=0;
#pragma omp parallel for
  for(ULONG i=0;i<NTT;++i)out_FOURIER[i][IMAG]=0;
  do_fftw_r2c(Nft,in_ref,in_FOURIER);
  do_fftw_r2c(Nft,out_field,out_FOURIER);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<NTT;++i)
   {
     real_prec Amp_in = sqrt(pow(in_FOURIER[i][REAL],2)+pow(in_FOURIER[i][IMAG],2));
     real_prec phi_o  = atan2(out_FOURIER[i][IMAG],out_FOURIER[i][REAL]);
     out_FOURIER[i][REAL] = Amp_in*cos(phi_o); 
     out_FOURIER[i][IMAG] = Amp_in*sin(phi_o);
   }
  do_fftw_c2r(Nft,out_FOURIER,out_field);
#ifdef _DOUBLE_PREC_
  fftw_free(in_FOURIER);
  fftw_free(out_FOURIER);
#else
  fftwf_free(in_FOURIER);
  fftwf_free(out_FOURIER);
#endif


}



//##################################################################################
//##################################################################################

void EigenValuesTweb_bias(ULONG Nft, real_prec L1, const vector<real_prec> &delta, const vector<real_prec> &phi, vector<real_prec> &S2, vector<real_prec> &S3, vector<real_prec> &N2D)
{
    // Function from webclass by F. Kitaura. File webclass.cc
    // The potential vctor is the gravitational potential obtained form the Poisson solver
    real_prec L2=L1;
    real_prec L3=L2;

#ifdef _USE_OMP_
    int NTHREADS = omp_get_max_threads();
    omp_set_num_threads(NTHREADS);
#endif

    ULONG NGRID = static_cast<ULONG>(Nft*Nft*Nft);


    vector<real_prec> LapPhivx(NGRID,0), LapPhivy(NGRID,0), LapPhivz(NGRID,0);
    vector<real_prec> LapPhivxy(NGRID,0), LapPhivxz(NGRID,0), LapPhivyz(NGRID,0);
    vector<real_prec> LapPhivzx(NGRID,0), LapPhivzy(NGRID,0), LapPhivyx(NGRID,0);



    vector<real_prec> dummy(NGRID,0);
    // Get gradient in x-direction
    gradfindif(Nft,Nft,Nft,L1,L2,L3,phi,dummy,1);
    // Get second derivatives of grad_x
    gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivx,1);
    gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivxy,2);
    gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivxz,3);

    // Get gradient in y-direction
    gradfindif(Nft,Nft,Nft,L1,L2,L3,phi,dummy,2);
    // Get second derivatives of grad_y
    gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivy,2);
    gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivyx,1);
    gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivyz,3);

    // Get gradient in z-direction
    gradfindif(Nft,Nft,Nft,L1,L2,L3,phi,dummy,3);
    // Get second derivatives of grad_z
    gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivz,3);
    gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivzx,1);
    gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivzy,2);


#if defined (_USE_S3_) || defined (_USE_S3_)
    vector<real_prec> La1(NGRID,0), La2(NGRID,0), La3(NGRID,0),Phi(NGRID,0);
    PoissonSolver(L1, Nft,delta,Phi);
    EigenValuesTweb(Nft,L1, delta, Phi, La1, La2, La3);
    Phi.clear();Phi.shrink_to_fit();
#endif



#ifdef _USE_S2_
    // Get S2:


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<NGRID;++i)
       {
        real_prec S11=LapPhivx[i]-(1./3.)*(La1[i]+La2[i]+La3[i]);
        real_prec S12=LapPhivxy[i];
        real_prec S13=LapPhivxz[i];
        real_prec S21=LapPhivyx[i];
        real_prec S22=LapPhivy[i]-(1./3.)*(La1[i]+La2[i]+La3[i]);
        real_prec S23=LapPhivyz[i];
        real_prec S31=LapPhivzx[i];
        real_prec S32=LapPhivzy[i];
        real_prec S33=LapPhivz[i]-(1./3.)*(La1[i]+La2[i]+La3[i]);
        S2[i]= S11 *S11 +S12*S21 + S13*S31+ S21*S12 +S22*S22 + S23*S32 +S31*S13 +S23*S32 +S33*S33;
       }
#endif



#ifdef _USE_S3_
    // Get S3:
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<NGRID;++i)
        {
         real_prec S11=LapPhivx[i]-(1./3.)*(La1[i]+La2[i]+La3[i]);
         real_prec S12=LapPhivxy[i];
         real_prec S13=LapPhivxz[i];
         real_prec S21=LapPhivyx[i];
         real_prec S22=LapPhivy[i]-(1./3.)*(La1[i]+La2[i]+La3[i]);
         real_prec S23=LapPhivyz[i];
         real_prec S31=LapPhivzx[i];
         real_prec S32=LapPhivzy[i];
         real_prec S33=LapPhivz[i]-(1./3.)*(La1[i]+La2[i]+La2[i]);
         real_prec A1=S11*(S11*S11+S12*S21+S12*S31)+S21*(S11*S12+S12*S21+S12*S31)+S31*(S11*S13+S12*S23+S13*S33);
         real_prec A2=S12*(S21*S11+S22*S21+S23*S31)+S22*(S21*S12+S22*S22+S23*S32)+S32*(S21*S13+S22*S23+S23*S33);
         real_prec A3=S13*(S31*S11+S32*S21+S33*S31)+S23*(S31*S12+S32*S22+S33*S32)+S33*(S31*S31+S32*S23+S33*S33);
         S3[i]=A1+A2+A3;
     }

#endif

#if defined (_USE_S3_) || defined (_USE_S3_)
    La1.clear(); La1.shrink_to_fit();
    La2.clear(); La2.shrink_to_fit();
    La3.clear(); La3.shrink_to_fit();
#endif


#ifdef _USE_NABLA2DELTA_

    vector<real_prec> NLapPhivx(NGRID,0), NLapPhivy(NGRID,0), NLapPhivz(NGRID,0);

    // Get partial delta/partial x
    gradfindif(Nft,Nft,Nft,L1,L2,L3,delta,NLapPhivx,1);
    // Get partial delta/partial y
    gradfindif(Nft,Nft,Nft,L1,L2,L3,delta,NLapPhivy,2);
    // Get partial delta/partial z
    gradfindif(Nft,Nft,Nft,L1,L2,L3,delta,NLapPhivz,3);

    // Get partial( delta/partial x) /partial x
    gradfindif(Nft,Nft,Nft,L1,L2,L3,NLapPhivx,NLapPhivx,1);
    // Get partial( delta/partial y) /partial y
    gradfindif(Nft,Nft,Nft,L1,L2,L3,NLapPhivy,NLapPhivy,2);
    // Get partial( delta/partial z) /partial z
    gradfindif(Nft,Nft,Nft,L1,L2,L3,NLapPhivz,NLapPhivz,3);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<NGRID;++i)
      N2D[i]=NLapPhivx[i]+NLapPhivy[i]+NLapPhivz[i];
    NLapPhivx.clear();
    NLapPhivx.shrink_to_fit();
    NLapPhivy.clear();
    NLapPhivy.shrink_to_fit();
    NLapPhivz.clear();
    NLapPhivz.shrink_to_fit();

#endif


    dummy.clear();
    dummy.shrink_to_fit();

}

//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
void EigenValuesVweb(ULONG Nft, real_prec L1, vector<real_prec>&Vinx, vector<real_prec>&Viny,vector<real_prec>&Vinz, vector<real_prec>&diver, vector<real_prec> &out1, vector<real_prec> &out2, vector<real_prec> &out3)
{

    //Get divergence of vel field and eigenvalues of velocity shear tensor
    // Function from webclass by F. Kitaura. File webclass.cc

#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif


    real_prec L2=L1;
    real_prec L3=L2;

    ULONG NGRID = static_cast<ULONG>(Nft*Nft*Nft);
    vector<real_prec> vxy(NGRID,0),vyx(NGRID,0), vzy(NGRID,0);
    vector<real_prec> vxz(NGRID,0),vyz(NGRID,0), vzx(NGRID,0);

    // Sheer tensor, 0.5(dv_i/ dx_j + dv_j/dx_i), symmetric, only 6 independent components

    vector<real_prec> Svxx(NGRID,0), Svyy(NGRID,0),Svzz(NGRID,0);
    vector<real_prec> Svxy(NGRID,0), Svxz(NGRID,0),Svyz(NGRID,0);
    // -------------------------------------------------------------

    // Get gradient in x-direction
    gradfindif(Nft,Nft,Nft,L1,L2,L3,Vinx,Svxx,1); // dVx/dx
    gradfindif(Nft,Nft,Nft,L1,L2,L3,Vinx,vxy,2); //dVx/dy
    gradfindif(Nft,Nft,Nft,L1,L2,L3,Vinx,vxz,3); // dVx/dz

    // -------------------------------------------------------------
    // Get gradient in y-direction
    gradfindif(Nft,Nft,Nft,L1,L2,L3,Viny,Svyy,2); // dVy/dy
    gradfindif(Nft,Nft,Nft,L1,L2,L3,Viny,vyx,1);// dVy/dx
    gradfindif(Nft,Nft,Nft,L1,L2,L3,Viny,vyz,3);// dVy/dz

    // -------------------------------------------------------------

    // Get gradient in z-direction
    gradfindif(Nft,Nft,Nft,L1,L2,L3,Vinz,Svzz,3); // dVz/dz
    gradfindif(Nft,Nft,Nft,L1,L2,L3,Vinz,vzx,1); // dVz/dx
    gradfindif(Nft,Nft,Nft,L1,L2,L3,Vinz,vzy,2); // dVz/dy

    // -------------------------------------------------------------
    // -------------------------------------------------------------
    // Get divergence
    
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0; i<NGRID;++i)
        diver[i]=Svxx[i]+Svyy[i]+Svzz[i];


    // Get shear
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<NGRID;i++)
      {
        Svxy[i]=-0.5*(vxy[i]+vyx[i]);
        Svxz[i]=-0.5*(vxz[i]+vzx[i]);
        Svyz[i]=-0.5*(vyz[i]+vzy[i]);
      }

    vxy.clear();
    vxy.shrink_to_fit();
    vxz.clear();
    vxz.shrink_to_fit();
    vyx.clear();
    vyx.shrink_to_fit();
    vzx.clear();
    vzx.shrink_to_fit();
    vzy.clear();
    vzy.shrink_to_fit();


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<NGRID;i++)
      {
        // M = [ A B C,
        //       D E F,
        //       G H I ]

        // Eqn = aλ³ + bλ² + cλ + d = 0
        // =-λ³+λ²(A+E+I)+λ(DB+GC+FH-AE-AI-EI)+(AEI-AFH-DBI+DCH+GBF-GCE)
        // a=-1
        // b=tr(M)
        // c=1/2[tr(M*M)-tr(M)*tr(M)]
        // d=det(M)

        real_prec A=Svxx[i];
        real_prec B=Svxy[i];
        real_prec C=Svxz[i];

        real_prec D=Svxy[i]; // component yx=xy

        real_prec E=Svyy[i];
        real_prec F=Svyz[i];

        real_prec G=Svxz[i];//Svzx=Svxz


        real_prec H=Svyz[i];//Svzy=Svyz
        real_prec I=Svzz[i];

        real_prec trM=A+E+I;
        real_prec detM=A*E*I-A*F*H-D*B*I+D*C*H+G*B*F-G*C*E;

        real_prec trM2=D*B+G*C+F*H-A*E-A*I-E*I;

        // Eqn = a1λ³ + a2λ² + a3λ + a4 = 0
        real_prec a1= -num_1;
        real_prec a2= trM;
        real_prec a3= trM2;
        real_prec a4= detM;

        // Eqn = λ³ + aλ² + bλ + c = 0
        gsl_real a= static_cast<gsl_real>(a2/a1);
        gsl_real b= static_cast<gsl_real>(a3/a1);
        gsl_real c= static_cast<gsl_real>(a4/a1);


        gsl_real x0, x1, x2;
        gsl_poly_solve_cubic (a,  b,  c, &x0, &x1, &x2);

        real_prec Eig1 = x0;
        real_prec Eig2 = x1;
        real_prec Eig3 = x2;


        if (Eig1 >= Eig2)
          {
            if (Eig1>=Eig3)
              {
                out1[i]=Eig1;
                if (Eig2>=Eig3)
                  {
                    out2[i]=Eig2;
                    out3[i]=Eig3;
                  }
                else
                  {
                    out2[i]=Eig3;
                    out3[i]=Eig2;
                  }
              }
            else
              {
                out1[i]=Eig3;
                if (Eig1>=Eig2)
                  {
                    out2[i]=Eig1;
                    out3[i]=Eig2;
                  }
                else
                  {
                    out2[i]=Eig2;
                    out3[i]=Eig1;
                  }
              }
          }
        else
          {
            if (Eig2>=Eig3)
              {
                out1[i]=Eig2;
                if (Eig1>=Eig3)
                  {
                    out2[i]=Eig1;
                    out3[i]=Eig3;
                  }
                else
                  {
                    out2[i]=Eig3;
                    out3[i]=Eig1;
                  }
              }
            else
              {
                out1[i]=Eig3;
                if (Eig1>=Eig2)
                  {
                    out2[i]=Eig1;
                    out3[i]=Eig2;
                  }
                else
                  {
                    out2[i]=Eig2;
                    out3[i]=Eig1;
                  }
              }
          }
      }
}

//##################################################################################
//##################################################################################
//##################################################################################

void PoissonSolver(real_prec Lbox, ULONG Nft,  vector<real_prec>&in, vector<real_prec>&out)
{


#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif

#ifdef _USE_ZERO_PADDING_POT_
    int factor_ntt=_EXTRA_NFT_FACTOR_*_EXTRA_NFT_FACTOR_*_EXTRA_NFT_FACTOR_;
    ULONG NGRIDnew=factor_ntt*Nft*Nft*Nft;
    vector<real_prec>new_in(NGRIDnew,0);
#pragma omp parallel for
    for(ULONG i=0;i<_EXTRA_NFT_FACTOR_*Nft;++i)
        for(ULONG j=0;j<_EXTRA_NFT_FACTOR_* Nft;++j)
            for(ULONG k=0;k<_EXTRA_NFT_FACTOR_* Nft;++k)
            {
                ULONG index_padd=index_3d(i,j,k,_EXTRA_NFT_FACTOR_*Nft,_EXTRA_NFT_FACTOR_*Nft);
                ULONG index_or=index_3d(i,j,k,Nft,Nft);
                if(index_or<Nft*Nft*Nft)
                   new_in[index_padd]=in[index_or];

            }
    Nft*=_EXTRA_NFT_FACTOR_;
#endif

  ULONG NTT = static_cast<ULONG>(Nft*Nft*(Nft/2+1));
  real_prec deltak=2.*M_PI/Lbox;
  

#ifdef DOUBLE_PREC
   complex_prec *data_out= (complex_prec *)fftw_malloc(2*NTT*sizeof(real_prec));
#else
   complex_prec *data_out= (complex_prec *)fftwf_malloc(2*NTT*sizeof(real_prec));
#endif



#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<NTT;++i)data_out[i][REAL]=0;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<NTT;++i)data_out[i][IMAG]=0;
  
#ifdef _USE_ZERO_PADDING_POT_
  do_fftw_r2c(Nft,new_in,data_out);
  new_in.clear();new_in.shrink_to_fit();
#else
  do_fftw_r2c(Nft,in,data_out);
#endif

  
  vector<real_prec> coords(Nft,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<Nft ;++i)
    coords[i]=deltak*(i<=Nft/2? static_cast<real_prec>(i): -static_cast<real_prec>(Nft-i));
  
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i< Nft ; i++)
    for(ULONG j=0;j< Nft ; j++)
      for(ULONG k=0;k< Nft/2+1; k++)
        {
	  ULONG ind=index_3d(i,j,k, Nft, Nft/2+1);
          real_prec kmod2=pow(coords[i],2)+pow(coords[j],2)+pow(coords[k],2);  // Get k**2
          real_prec term=0.;
	  if(kmod2>0.)
            term=static_cast<real_prec>(-1./kmod2);
          
          data_out[ind][REAL]*=term;   // Divide by k**2
          data_out[ind][IMAG]*=term;
        }
  coords.clear();
  coords.shrink_to_fit();


#ifdef _USE_ZERO_PADDING_POT_
  vector<real_prec> newout(Nft*Nft*Nft,0);
  do_fftw_c2r(Nft ,data_out,newout);
#pragma omp parallel for
  for(ULONG i=0;i<Nft*Nft*Nft/factor_ntt;++i)out[i]=newout[i];
  newout.clear();newout.shrink_to_fit();
 #else
  do_fftw_c2r(Nft ,data_out,out);

#endif


#ifdef _DOUBLE_PREC_
  fftw_free(data_out);
#else
  fftwf_free(data_out);
#endif

}
//##################################################################################
//##################################################################################
void exchange_xy(int Nft,const vector<real_prec>&in, vector<real_prec>&out)
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
for(ULONG i=0;i<Nft;++i)
  for(ULONG j=0;j<Nft;++j)
    for(ULONG k=0;k<Nft;++k)
      out[index_3d(i,j,k,Nft,Nft)]=in[index_3d(j,i,k,Nft,Nft)];
}


//##################################################################################
//##################################################################################
void exchange_xz(int Nft,const vector<real_prec>&in, vector<real_prec>&out)
{
 // this is alsme meant to be column (F) to major (c) order
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<Nft;++i)
      for(ULONG j=0;j<Nft;++j)
         for(ULONG k=0;k<Nft;++k)
             out[index_3d(i,j,k,Nft,Nft)]=in[index_3d(k,j,i,Nft,Nft)];
}



//##################################################################################
//##################################################################################

ULONG index_10d(int i, int j, int k, int l, int m, int n, int o, int p, int q, int r, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np, int Nq, int Nr)
{
  return  static_cast<ULONG>(r) + static_cast<ULONG>(Nr*q) +  static_cast<ULONG>(Nr*Nq*p)+ static_cast<ULONG>(Nr*Nq*Np*o) + static_cast<ULONG>(Nr*Nq*Np*No*n)+static_cast<ULONG>(Nr*Nq*Np*No)*static_cast<ULONG>(Nn*m)+static_cast<ULONG>(Nr*Nq*Np*No)*static_cast<ULONG>(Nn*Nm*l)+static_cast<ULONG>(Nr*Nq*Np*No)*static_cast<ULONG>(Nn*Nm*Nl*k)+static_cast<ULONG>(Nr*Nq*Np*No)*static_cast<ULONG>(Nn*Nm*Nl*Nk*j)+static_cast<ULONG>(Nr*Nq*Np*No)*static_cast<ULONG>(Nn*Nm*Nl*Nk*Nj*i);
}


ULONG index_11d(int i, int j, int k, int l, int m, int n, int o, int p, int q, int r, int s, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np, int Nq, int Nr, int Ns)
{
 return  static_cast<ULONG>(s) + static_cast<ULONG>(Ns*r)+static_cast<ULONG>(Ns*Nq*q)+ static_cast<ULONG>(Ns*Nr*Nq*p)+ static_cast<ULONG>(Ns*Nr*Nq*Np*o) +static_cast<ULONG>(Ns*Nr*Nq*Np*No*n)
            +static_cast<ULONG>(Ns*Nr*Nq)*static_cast<ULONG>(Np*No*Nn*m)+static_cast<ULONG>(Ns*Nr*Nq)*static_cast<ULONG>(Np*No*Nn*Nm*l)+static_cast<ULONG>(Ns*Nr*Nq*Np*No)*static_cast<ULONG>(Nn*Nm*Nl*k)+static_cast<ULONG>(Ns*Nr*Nq*Np)*static_cast<ULONG>(No*Nn*Nm*Nl*Nk*j)+static_cast<ULONG>(Ns*Nr*Nq*Np)*static_cast<ULONG>(No*Nn*Nm*Nl*Nk*Nj*i);
}



ULONG index_12d(int i, int j, int k, int l, int m, int n, int o, int p, int q, int r, int s, int v, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np, int Nq, int Nr, int Ns, int Nv)
{

 return static_cast<ULONG>(v) + static_cast<ULONG>(Nv*s) + static_cast<ULONG>(Nv*Ns*r)+static_cast<ULONG>(Nv*Ns*Nr*q)+static_cast<ULONG>(Nv*Ns*Nr*Nq*p)+ static_cast<ULONG>(Nv*Ns*Nr*Nq*Np*o)+ static_cast<ULONG>(Nv*Ns*Nr*Nq*Np*No*n) +static_cast<ULONG>(Nv*Ns*Nr*Nq*Np*No*Nn*m)
            +static_cast<ULONG>(Nv*Ns*Nr*Nq)*static_cast<ULONG>(Np*No*Nn*Nm*l)+static_cast<ULONG>(Nv*Ns*Nr*Nq)*static_cast<ULONG>(Np*No*Nn*Nm*Nl*k)+static_cast<ULONG>(Nv*Ns*Nr*Nq*Np*No)*static_cast<ULONG>(Nn*Nm*Nl*Nk*j)+static_cast<ULONG>(Nv*Ns*Nr*Nq*Np)*static_cast<ULONG>(No*Nn*Nm*Nl*Nk*Nj*i);
}


//##################################################################################
//##################################################################################

ULONG index_9d(int i, int j, int k, int l, int m, int n, int o, int p, int q, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np, int Nq)
{
  return  static_cast<ULONG>(q) + static_cast<ULONG>(Nq*p) +  static_cast<ULONG>(Nq*Np*o)+ static_cast<ULONG>(Nq*Np*No*n) + static_cast<ULONG>(Nq*Np*No*Nn*m)+static_cast<ULONG>(Nq*Np*No*Nn)*static_cast<ULONG>(Nm*l)+static_cast<ULONG>(Nq*Np*No*Nn)*static_cast<ULONG>(Nm*Nl*k)+static_cast<ULONG>(Nq*Np*No*Nn)*static_cast<ULONG>(Nm*Nl*Nk*j)+static_cast<ULONG>(Nq*Np*No*Nn)*static_cast<ULONG>(Nm*Nl*Nk*Nj*i);
}


//##################################################################################
//##################################################################################
ULONG index_8d(int i, int j, int k, int l, int m, int n, int o, int p, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np)
{
  return  static_cast<ULONG>(p) + static_cast<ULONG>(Np*o) +  static_cast<ULONG>(Np*No*n)+ static_cast<ULONG>(Np*No*Nn*m) + static_cast<ULONG>(Np*No*Nn*Nm*l)+static_cast<ULONG>(Np*No*Nn*Nm)*static_cast<ULONG>(Nl*k)+static_cast<ULONG>(Np*No*Nn*Nm)*static_cast<ULONG>(Nl*Nk*j)+static_cast<ULONG>(Np*No*Nn*Nm)*static_cast<ULONG>(Nl*Nk*Nj*i);
}

//##################################################################################
//##################################################################################
ULONG index_7d(int i, int j, int k, int l, int m, int n, int o, int Nj, int Nk, int Nl, int Nm, int Nn, int No)
{
  return  static_cast<ULONG>(o) + static_cast<ULONG>(No*n) +  static_cast<ULONG>(No*Nn*m)+ static_cast<ULONG>(No*Nn*Nm*l) + static_cast<ULONG>(No*Nn*Nm*Nl*k)+static_cast<ULONG>(No*Nn*Nm*Nl)*static_cast<ULONG>(Nk*j)+static_cast<ULONG>(No*Nn*Nm*Nl)*static_cast<ULONG>(Nk*Nj*i);
}


//##################################################################################
//##################################################################################

ULONG index_6d(int i, int j, int k, int l, int m, int n, int Nj, int Nk, int Nl, int Nm, int Nn)
{
  return  static_cast<ULONG>(n) + static_cast<ULONG>(Nn*m) +  static_cast<ULONG>(Nn*Nm*l)+ static_cast<ULONG>(Nn*Nm*Nl*k) + static_cast<ULONG>(Nn*Nm*Nl*Nk*j)+static_cast<ULONG>(Nn*Nm*Nl)*static_cast<ULONG>(Nk*Nj*i);
}
//##################################################################################
//##################################################################################
ULONG index_5d(int i, int j, int k, int l, int m, int Nj, int Nk, int Nl, int Nm)
{
  return static_cast<ULONG>(m) + static_cast<ULONG>(Nm*l)+ static_cast<ULONG>(Nm*Nl*k)+static_cast<ULONG>(Nm*Nl*Nk*j) +  static_cast<ULONG>(Nm*Nl*Nk*Nj*i);
}
//##################################################################################
//##################################################################################
ULONG index_4d(int i, int j, int k, int l, int Nj, int Nk, int Nl)
{
  return static_cast<ULONG>(l)+ static_cast<ULONG>(Nl*k)+static_cast<ULONG>(Nl*Nk*j)+ static_cast<ULONG>(Nl*Nk*Nj*i);
}
//##################################################################################
//##################################################################################
ULONG index_3d(ULONG i, ULONG j, ULONG k, int Nj, int Nk)
{ // if k denots the z coordinate, then this is row majort (c) order
  return k+ static_cast<ULONG>(Nk*j)+static_cast<ULONG>(Nk*Nj*i);
}
//##################################################################################
//##################################################################################
ULONG index_2d(ULONG i, ULONG j, ULONG Nj)
{
  return j+Nj*i;
}

//##################################################################################
//##################################################################################

// This funciton can be applied to both C like (row-major) or Fortran like (columnm major) order

void index2coords(ULONG index, ULONG N,  ULONG  &XG, ULONG &YG, ULONG &ZG )
{
//  see https://math.stackexchange.com/questions/3758576/how-to-convert-from-a-flattened-3d-index-to-a-set-of-coordinates

  ZG=index % N;  // Fast varying coordinate. Can be z(c++) or x (Fortran)
  index = static_cast<ULONG>(static_cast<real_prec>(index)/static_cast<real_prec>(N));
  YG=index % N;
  XG = static_cast<ULONG>(static_cast<real_prec>(index)/static_cast<real_prec>(N));  // Slow varying coordinate. Can be x(c++) or z(Fortran)
}


//##################################################################################
//#################################################################################
ULONG grid_ID(s_params_box_mas *params, const real_prec &x, const real_prec &y, const real_prec &z)
{
  ULONG i = static_cast<ULONG>(floor((x-params->min1)/params->d1)); // indices of the cell of the particle
  ULONG j = static_cast<ULONG>(floor((y-params->min2)/params->d2));
  ULONG k = static_cast<ULONG>(floor((z-params->min3)/params->d3));

  if(i==params->Nft)
    i--;
  if(j==params->Nft)
    j--;
  if(k==params->Nft)
    k--;
  
  i = static_cast<ULONG>(fmod(real_prec(i),real_prec(params->Nft)));
  j = static_cast<ULONG>(fmod(real_prec(j),real_prec(params->Nft)));
  k = static_cast<ULONG>(fmod(real_prec(k),real_prec(params->Nft)));
  
  return index_3d(i,j,k,params->Nft,params->Nft);
}
//##################################################################################
//##################################################################################
//##################################################################################
void do_fftw_r2c(int Nft, vector<real_prec>&in, complex_prec *out)
{

#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

#ifdef DOUBLE_PREC
  int *n=(int *)fftw_malloc(ic_rank*sizeof(int));
#else
  int *n=(int *)fftwf_malloc(ic_rank*sizeof(int));
#endif


  for(int i=0;i<ic_rank;++i)n[i]=Nft;
#ifdef SINGLE_PREC
#ifdef _USE_OMP_
  fftwf_init_threads();
  fftwf_plan_with_nthreads(NTHREADS);
#endif
  fftwf_plan plan_r2c=fftwf_plan_dft_r2c(ic_rank,n,&in[0],out,FFTW_ESTIMATE);
  fftwf_execute(plan_r2c);
  fftwf_destroy_plan(plan_r2c);
  fftwf_free(n);
#endif
#ifdef DOUBLE_PREC
#ifdef _USE_OMP_
  fftw_init_threads();
  fftw_plan_with_nthreads(NTHREADS);
#endif
  fftw_plan plan_r2c=fftw_plan_dft_r2c(ic_rank,n,&in[0],out,FFTW_ESTIMATE);
  fftw_execute(plan_r2c);
  fftw_destroy_plan(plan_r2c);
  fftw_free(n);
#endif
}

//##################################################################################
//##################################################################################
void do_fftw_c2r(int Nft, complex_prec *in, vector<real_prec>&out)
{

#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  ULONG NGRID = static_cast<ULONG>(Nft*Nft*Nft);


#ifdef DOUBLE_PREC
  int *n=(int *)fftw_malloc(ic_rank*sizeof(int));
#else
  int *n=(int *)fftwf_malloc(ic_rank*sizeof(int));
#endif


  for(int i=0;i<ic_rank;++i)n[i]=Nft;
#ifdef SINGLE_PREC
#ifdef _USE_OMP_
  fftwf_init_threads();
  fftwf_plan_with_nthreads(omp_get_max_threads());
#endif
  fftwf_plan plan_c2r=fftwf_plan_dft_c2r(ic_rank,n,in,&out[0], FFTW_ESTIMATE);
  fftwf_execute(plan_c2r);
  fftwf_destroy_plan(plan_c2r);
  fftwf_free(n);
#endif
#ifdef DOUBLE_PREC
#ifdef _USE_OMP_
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif
  fftw_plan plan_c2r=fftw_plan_dft_c2r(ic_rank,n,in,&out[0], FFTW_ESTIMATE);
  fftw_execute(plan_c2r);
  fftw_destroy_plan(plan_c2r);
  fftw_free(n);
#endif


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<NGRID;i++)  // normalize the c2r transform
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
      out[i]/=static_cast<double>(NGRID);

}


//##################################################################################
void do_fftw_3d(ULONG Nft, bool direction, complex_prec *in, complex_prec *out)
{

#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

  ULONG factor=Nft*Nft*Nft;
  
  int Ni1=static_cast<int>(Nft);
  int Ni2=static_cast<int>(Nft);
  int Ni3=static_cast<int>(Nft);
  

  if (direction == true)//from conf space to fourier
    {
#ifdef SINGLE_PREC
      fftwf_plan fftp;
      fftp = fftwf_plan_dft_3d(Ni1,Ni2,Ni3,in,out,FORWARD,FFTW_OPTION);	
      fftwf_execute(fftp);
#endif
#ifdef DOUBLE_PREC
      fftw_plan fftp;
      fftp = fftw_plan_dft_3d(Ni1,Ni2,Ni3,in,out,FORWARD,FFTW_OPTION);	
      fftw_execute(fftp);
#endif
      
    }
  else // from fourier to conf space
    {
#ifdef SINGLE_PREC
      fftwf_plan fftp;
      fftp = fftwf_plan_dft_3d(Ni1,Ni2,Ni3,in,out,BACKWARD,FFTW_OPTION);
      fftwf_execute(fftp);
#endif
#ifdef DOUBLE_PREC
      fftw_plan fftp;
      fftp = fftw_plan_dft_3d(Ni1,Ni2,Ni3,in,out,BACKWARD,FFTW_OPTION);	
      fftw_execute(fftp);
#endif
      for(ULONG i=0;i<factor;i++)
	{// normalize the c2r transform
          out[i][REAL]/=static_cast<real_prec>(factor);
          out[i][IMAG]/=static_cast<real_prec>(factor);
	}

#ifdef SINGLE_PREC
      fftwf_destroy_plan(fftp);
#endif
#ifdef DOUBLE_PREC
      fftw_destroy_plan(fftp);
#endif


    }

  
}


//##################################################################################
//##################################################################################

void get_cumulative(const vector<real_prec> &dist, vector<real_prec> &cumu, unsigned long &NTOT_h)
{

#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

  NTOT_h=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:NTOT_h)
#endif
  for(int i=0;i<dist.size();++i)
    NTOT_h+=dist[i];
  
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int i=0;i<dist.size();++i)
    for(int j=0;j<=i;++j)
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
        cumu[i]+=dist[j];

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<dist.size();++i)
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
      cumu[i]/=static_cast<double>(NTOT_h);
}

//##################################################################################
//##################################################################################
real_prec get_nobjects(const vector<real_prec> &density_field)
{

#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

  double ans=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ans)
#endif
  for(ULONG i=0;i< density_field.size();++i)
    ans+=static_cast<double>(density_field[i]);

  return static_cast<real_prec>(ans);
}

//##################################################################################
ULONG get_nobjects(const vector<ULONG> &density_field)
{

#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

  ULONG ans=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ans)
#endif
  for(ULONG i=0;i< density_field.size();++i)
    ans+=static_cast<ULONG>(density_field[i]);

  return ans;
}

//##################################################################################
//##################################################################################
void get_overdens(const vector<real_prec>&in, vector<real_prec>&out)
{

#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

#ifdef _FULL_VERBOSE_
    cout<<YELLOW<<"Converting to Overdensity"<<RESET<<endl;
#endif

  double mean=get_nobjects(in);
  mean/=static_cast<double>(in.size());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i< in.size();++i)
    out[i]=in[i]/mean-num_1;
#ifdef _FULL_VERBOSE_
  cout<<BOLDGREEN<<"                                               ["<<BOLDBLUE<<"DONE"<<BOLDGREEN<<"]"<<RESET<<endl;
#endif
}


//##################################################################################
//##################################################################################
void convert_ngp_to_cic(vector<real_prec>&in, vector<real_prec>&out)
{
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

  cout<<YELLOW<<"Transforming NGP to CIC "<<RESET<<endl;
  ULONG N=in.size();
  ULONG N1=static_cast<ULONG>(pow(N,1./3.))+1;
  ULONG Nhalf=static_cast<ULONG>(N1)*static_cast<ULONG>(N1)*static_cast<ULONG>(N1/2+1);  	

#ifdef DOUBLE_PREC
  complex_prec *AUX= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *AUX= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<Nhalf;++i)
    {
      AUX[i][REAL]=0;
      AUX[i][IMAG]=0;
    }


  
  // Convert input field to Fouroer space
  do_fftw_r2c(N1,in,AUX);
  
  vector<real_prec> coords(N1,0);
  vector<real_prec> correction(N1,0);
  
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i = 0 ; i < N1; ++i )
    {
      coords[i]=static_cast<real_prec>( i<=N1/2? i: i-(N1));
      real_prec xx=coords[i]*M_PI/static_cast<real_prec>(N1);
      correction[i]= i==0 ? 1.0 : sin(xx)/static_cast<real_prec>(xx) ;
    }
  
  real_prec we=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:we)
#endif
  for(ULONG i=0;i< N1; i++)
    for(ULONG j=0;j< N1; j++)
      for(ULONG k=0;k< N1/2+1; k++)
	{
          real_prec cor=correction[i]*correction[j]*correction[k];
	  ULONG ind=index_3d(i,j,k, N1, N1/2+1);
	  AUX[ind][REAL]*=cor;
	  AUX[ind][IMAG]*=cor;
	  we+=cor;
	}

  do_fftw_c2r(N1,AUX,out);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<N;i++)
    out[i]=(out[i]<0 ? 0: out[i]);//*(static_cast<real_prec>(N)/static_cast<real_prec>(2.0*we));
  
#ifdef DOUBLE_PREC
  fftw_free(AUX);
#else
  fftwf_free(AUX);
#endif


#ifdef _FULL_VERBOSE_
  cout<<BOLDGREEN<<"                                               ["<<BOLDBLUE<<"DONE"<<BOLDGREEN<<"]"<<RESET<<endl;  
#endif
}






//##################################################################################
//##################################################################################
void convert_cic_to_ngp(vector<real_prec>&in, vector<real_prec>&out)
{

#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
#ifdef _FULL_VERBOSE_
  cout<<YELLOW<<"Transforming CIC to NGP "<<RESET<<endl;
#endif
  ULONG N=in.size();
  ULONG N1=static_cast<ULONG>(pow(N,1./3.))+1;
  ULONG Nhalf=static_cast<ULONG>(N1)*static_cast<ULONG>(N1)*static_cast<ULONG>(N1/2+1);

#ifdef DOUBLE_PREC
  complex_prec *AUX= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *AUX= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<Nhalf;++i)
    {
      AUX[i][REAL]=0;
      AUX[i][IMAG]=0;
    }
  // Convert input field to Fouroer space
  do_fftw_r2c(N1,in,AUX);

  vector<real_prec> coords(N1,0);
  vector<real_prec> correction(N1,0);

#pragma omp parallel for
  for(int i = 0 ; i < N1; ++i )
    {
      coords[i]=static_cast<real_prec>( i<=N1/2? i: i-(N1));
      real_prec xx=coords[i]*M_PI/static_cast<real_prec>(N1);
      real_prec corr =  i==0 ? 1.0 : static_cast<real_prec>(xx)/sin(xx) ;  //this is the inverse of sinc
      correction[i] = corr == 0. ? 1. : corr;
  }

#pragma omp parallel for
  for(ULONG i=0;i< N1; i++)
    for(ULONG j=0;j< N1; j++)
      for(ULONG k=0;k< N1/2+1; k++)
        {
          real_prec cor=correction[i]*correction[j]*correction[k];
          ULONG ind=index_3d(i,j,k, N1, N1/2+1);
          AUX[ind][REAL]*=cor;
          AUX[ind][IMAG]*=cor;
        }

  do_fftw_c2r(N1,AUX,out);

#pragma omp parallel for
  for(ULONG i=0;i<N;i++)
    out[i]=(out[i]<0 ? 0: out[i]);//*(static_cast<real_prec>(N)/static_cast<real_prec>(2.0*we));

#ifdef DOUBLE_PREC
  fftw_free(AUX);
#else
  fftwf_free(AUX);
#endif


#ifdef _FULL_VERBOSE_
  cout<<BOLDGREEN<<"                                               ["<<BOLDBLUE<<"DONE"<<BOLDGREEN<<"]"<<RESET<<endl;
#endif
}

//##################################################################################
//##################################################################################
void get_overdens(const vector<real_prec>&in, real_prec mean, vector<real_prec>&out)
{
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

#ifdef _FULL_VERBOSE_
  cout<<YELLOW<<"Converting to Overdensity"<<RESET<<endl;
#endif

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i< in.size();++i)
    out[i]=in[i]/static_cast<double>(mean)-1.0;

#ifdef _FULL_VERBOSE_
  cout<<BOLDGREEN<<"                                               ["<<BOLDBLUE<<"DONE"<<BOLDGREEN<<"]"<<RESET<<endl;
#endif
}

//##################################################################################
//##################################################################################

void get_overdens(const vector<real_prec>&in,const vector<real_prec>&weight,vector<real_prec>&out, bool wind_binary)
{

#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

#ifdef _FULL_VERBOSE_
  cout<<YELLOW<<"Converting to Overdensity with weights"<<RESET<<endl;
#endif
  double mean=0;
  double veff=0;

#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:mean, veff)
#endif
  for(ULONG i=0;i< in.size();++i)
   {
      mean+=static_cast<real_prec>(in[i]);
     if(false==wind_binary)
        veff+=static_cast<real_prec>(weight[i]);
     else
      if (weight[i]>0.)
        veff++;

  }
  mean/=static_cast<double>(veff);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i< in.size();++i)
    out[i]=in[i]/mean-weight[i];
#ifdef _FULL_VERBOSE_
  cout<<BOLDGREEN<<"                                               ["<<BOLDBLUE<<"DONE"<<BOLDGREEN<<"]"<<RESET<<endl;
#endif
}


//##################################################################################
//##################################################################################

real_prec MAS_CIC_public(real_prec x)
{
  //////////////////////////////////////////////////////////
  // CIC Mass assignment scheme.
  //////////////////////////////////////////////////////////
  x=fabs(x);
  if(x<1)
    return 1.-x;
  else
    return 0.0;
}



// *************************************************************************************************
// *************************************************************************************************
void grid_assignment_cic(real_prec deltax, ULONG Nft, real_prec Lside, real_prec x,real_prec y,real_prec z,real_prec weight, vector<real_prec>& field)
{
  // ******************************************************
  // Boundary conditions:
  // ******************************************************
  real_prec rdelta=1./deltax;
  while(x< 0){x+=Lside;}
  while(y< 0){y+=Lside;}
  while(z< 0){z+=Lside;}
  while(x>=Lside){x-=Lside;}
  while(y>=Lside){y-=Lside;}
  while(z>=Lside){z-=Lside;}

  // *****************************************************
  // *****************************************************
  // Determining the number of cell in each direction:   *
  // Output in the range [0,N-1]                         *
  // *****************************************************
  int xc = get_bin(x,0,Nft,deltax,true);
  int yc = get_bin(y,0,Nft,deltax,true);
  int zc = get_bin(z,0,Nft,deltax,true);

  // ******************************************************
  // ******************************************************
  // Identify the cell to which this particle belongs     *
  // with its centre.                                     *
  // ******************************************************
  real_prec xx  = deltax*(static_cast<real_prec>(xc)+0.5);
  real_prec yy  = deltax*(static_cast<real_prec>(yc)+0.5);
  real_prec zz  = deltax*(static_cast<real_prec>(zc)+0.5);

  // ******************************************************
  // ******************************************************
  // Center of each cell forwards                         *
  // ******************************************************
  real_prec xxf = deltax*(static_cast<real_prec>(xc)+1.5);
  real_prec yyf = deltax*(static_cast<real_prec>(yc)+1.5);
  real_prec zzf = deltax*(static_cast<real_prec>(zc)+1.5);

  // ******************************************************
  // ******************************************************
  // Center of each cell backwards                        *
  // ******************************************************
  real_prec xxb = deltax*(static_cast<real_prec>(xc)-0.5);
  real_prec yyb = deltax*(static_cast<real_prec>(yc)-0.5);
  real_prec zzb = deltax*(static_cast<real_prec>(zc)-0.5);

  // ******************************************************
  // ******************************************************
  // Particular positions (borders)
  int Xb=(xc==0 ? Nft: xc);
  int Yb=(yc==0 ? Nft: yc);
  int Zb=(zc==0 ? Nft: zc);

  int Xf=(xc==Nft-1 ? -1: xc);
  int Yf=(yc==Nft-1 ? -1: yc);
  int Zf=(zc==Nft-1 ? -1: zc);

  // ******************************************************

  vector<int>i_idx(MAX_MAS_DEG,0);
  i_idx[0]=Xb-1;
  i_idx[1]=xc;
  i_idx[2]=Xf+1;

  vector<int>j_idx(MAX_MAS_DEG,0);
  j_idx[0]=Yb-1;
  j_idx[1]=yc;
  j_idx[2]=Yf+1;

  vector<int>k_idx(MAX_MAS_DEG,0);
  k_idx[0]=Zb-1;
  k_idx[1]=zc;
  k_idx[2]=Zf+1;

  vector<real_prec>MAS_xx(MAX_MAS_DEG,0);
  MAS_xx[0]=MAS_CIC_public((xxb- x)*rdelta);
  MAS_xx[1]=MAS_CIC_public((xx - x)*rdelta);
  MAS_xx[2]=MAS_CIC_public((xxf- x)*rdelta);

  vector<real_prec>MAS_yy(MAX_MAS_DEG,0);
  MAS_yy[0]=MAS_CIC_public((yyb- y)*rdelta);
  MAS_yy[1]=MAS_CIC_public((yy - y)*rdelta);
  MAS_yy[2]=MAS_CIC_public((yyf- y)*rdelta);

  vector<real_prec>MAS_zz(MAX_MAS_DEG,0);
  MAS_zz[0]=MAS_CIC_public((zzb- z)*rdelta);
  MAS_zz[1]=MAS_CIC_public((zz - z)*rdelta);
  MAS_zz[2]=MAS_CIC_public((zzf- z)*rdelta);

  for(int ih=0;ih<MAX_MAS_DEG;++ih)
    for(int jh=0;jh<MAX_MAS_DEG;++jh)
      for(int kh=0;kh<MAX_MAS_DEG;++kh)
          field[index_3d(i_idx[ih], j_idx[jh], k_idx[kh],Nft,Nft)]+=weight*MAS_xx[ih]*MAS_yy[jh]*MAS_zz[kh];

}





//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################

#ifdef OMPPARRANGAR
real_prec GR_NUM(gsl_rng ** SEED, real_prec sigma ,int GR_METHOD, int jthread)
{
#else
real_prec GR_NUM(gsl_rng * SEED, real_prec sigma ,int GR_METHOD)
{
#endif

  real_prec val=0.;
  switch (GR_METHOD)
    {
    case 0:
      {
#ifdef OMPPARRANGAR
        val=static_cast<real_prec>(gsl_ran_gaussian(SEED[jthread], sigma));
#else
        val=static_cast<real_prec>(gsl_ran_gaussian(SEED, sigma));
#endif
      }
      break;

    case 1:
      {
#ifdef OMPPARRANGAR
        val=static_cast<real_prec>(gsl_ran_gaussian_ziggurat(SEED[jthread], sigma));
#else
        val=static_cast<real_prec>(gsl_ran_gaussian_ziggurat(SEED, sigma));
#endif
      }
      break;

    case 2:
      {
#ifdef OMPPARRANGAR
        val=static_cast<real_prec>(gsl_ran_gaussian_ratio_method(SEED[jthread], sigma));
#else
        val=static_cast<real_prec>(gsl_ran_gaussian_ratio_method(SEED, sigma));
#endif
      }
      break;

    }
  return(val);
}

//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################

#ifdef OMPPARRANGAR
 void create_GARFIELDR(ULONG N1,ULONG N2,ULONG N3,  vector<real_prec> &delta, const vector<real_prec> &Power, gsl_rng ** seed)
{
#else
  void create_GARFIELDR(ULONG N1,ULONG N2,ULONG N3,  vector<real_prec> &delta, const vector<real_prec> &Power, gsl_rng * seed)
{
#endif
  
  #ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

  //  FileOutput File;
  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
  //  ULONG Nhalf=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3/2+1);
  ULONG Nhalf = static_cast<ULONG>(N1*N2*(N3/2+1));

  #ifdef _FULL_VERBOSE_
  cout<<YELLOW<<"Computing white noise "<<RESET<<endl;
#endif  
#ifdef DOUBLE_PREC
  complex_prec *AUX= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *AUX= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif
  

  for (ULONG i=0;i<N;i++)
    delta[i] = static_cast<real_prec>(GR_NUM(seed,num_1,0));
  
//   string fname=string("white_noise");
  //File.write_array(fname,delta,N);  // no puedo llamar a la clase FO
  
  do_fftw_r2c(N1,delta,AUX);
  //#endif

#ifdef _FULL_VERBOSE_
    cout<<YELLOW<<"Going to the Fourier grid "<<RESET<<endl;
#endif    

#ifdef _USE_OMP_
#pragma omp parallel for// take care!!!
#endif
    for (ULONG i=0 ; i<N1;i++)
      for (ULONG j=0 ; j<N2;j++)
	for (ULONG k=0 ; k<=N3/2;k++)
	  {
	    ULONG ihalf= index_3d(i,j,k,N2,N3/2+1);
	    ULONG iind = index_3d(i,j,k,N2,N3);
	    real_prec mass=0.0;
#ifdef FOURIER_DEF_1
	    mass=Power[iind];//testing
#endif
#ifdef FOURIER_DEF_2
            mass=Power[iind]/static_cast<real_prec>(N);
#endif


            real_prec sigma = sqrt(mass);
            AUX[ihalf][REAL]*=sigma;
            AUX[ihalf][IMAG]*=sigma;

/*
            // usa esta si queremos hacer gausr af en F space,m en ves de crear el delta y tranformarlo
            real_prec phase =2.*M_PI*gsl_rng_uniform(seed);
            real_prec sigma=gsl_ran_rayleigh(seed, sqrt(0.5*Power[iind]));
            AUX[ihalf][REAL]=sigma*cos(phase);
            AUX[ihalf][IMAG]=-sigma*sin(phase);
*/

            }
    
    do_fftw_c2r(N1,AUX,delta); 
#ifdef DOUBLE_PREC
    fftw_free(AUX);
#else
    fftwf_free(AUX);
#endif

  }
  

//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
void create_GARFIELD_FIXED_AMP(ULONG N1,ULONG N2,ULONG N3,  vector<real_prec> &delta, const vector<real_prec> &Power, gsl_rng * seed)
 {


#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif


   ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
   ULONG Nhalf = static_cast<ULONG>(N1*N2*(N3/2+1));
   cout<<YELLOW<<"Computing white noise "<<N2<<RESET<<endl;

#ifdef DOUBLE_PREC
   complex_prec *AUX= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
   complex_prec *AUX= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));

#endif
#ifdef _FULL_VERBOSE_
   cout<<YELLOW<<"Going to the Fourier grid "<<RESET<<endl;
#endif
 #ifdef _USE_OMP_
 #pragma omp parallel for// take care!!!
 #endif
     for (ULONG i=0 ; i<N1;i++)
       for (ULONG j=0 ; j<N2;j++)
         for (ULONG k=0 ; k<=N3/2;k++)
           {
             ULONG iind = index_3d(i,j,k,N2,N3);
             real_prec sigma = sqrt(Power[iind]); // In order to generate WN set sigma=1
             real_prec phase =2.*M_PI*gsl_rng_uniform(seed);

             ULONG ihalf= index_3d(i,j,k,N2,N3/2+1);
             AUX[ihalf][REAL]=sigma*cos(phase);
             AUX[ihalf][IMAG]=sigma*sin(phase);
           }


     do_fftw_c2r(N1,AUX,delta);
     real_prec meanWN=get_mean(delta);
     cout<<YELLOW<<"Mean WN = "<<CYAN<<meanWN<<RESET<<endl;
     real_prec sigma2D=get_var(meanWN, delta);
     cout<<YELLOW<<"Sigma_corr WN = "<<CYAN<<sqrt(sigma2D)<<RESET<<endl;


#ifdef DOUBLE_PREC
      fftw_free(AUX);
#else
     fftwf_free(AUX);
#endif

   }


//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
  
void create_GARFIELDR_from_WHITENOISE(string power_file,ULONG N1,ULONG N2,ULONG N3, vector<real_prec> &in)
{


#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif


  FileOutput File;

  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);	
  ULONG Nhalf=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3/2+1);	

  vector<real_prec>Power(N,0);
  File.read_array(power_file,Power);
#ifdef DOUBLE_PREC
  complex_prec *AUX= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *AUX= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif

  do_fftw_r2c(N1,in,AUX);
  

  real_prec factor=1;
#ifdef FOURIER_DEF_2
          factor=1./static_cast<real_prec>(N);
#endif
  
  cout<<YELLOW<<"Going into Fourier Loop to apply sqrt(P(k)) FT(WN)"<<RESET<<endl;
  
#pragma omp parallel for
  for (ULONG i=0 ; i<N1;i++)
    for (ULONG j=0 ; j<N2;j++)
      for (ULONG k=0 ; k<N3/2+1;k++)
	{
          ULONG ihalf= index_3d(i,j,k,N2,N3/2+1);
	  ULONG iind = index_3d(i,j,k,N2,N3);
          real_prec mass=factor*Power[iind];
          real_prec sigma = sqrt(mass);
	  AUX[ihalf][REAL]*=sigma;
	  AUX[ihalf][IMAG]*=sigma;
         }
  
  do_fftw_c2r(N1,AUX,in);      


#ifdef DOUBLE_PREC
  fftw_free(AUX);
#else
  fftwf_free(AUX);
#endif
}
  

//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
real_prec k_squared(ULONG i,ULONG j,ULONG k,real_prec L1,real_prec L2,real_prec L3,ULONG N1,ULONG N2,ULONG N3)
  {
    real_prec k2=0.;
    real_prec  kfac=static_cast<real_prec>(2.*M_PI);

    real_prec deltak_x=static_cast<real_prec>(2.*M_PI)/L1;
    real_prec deltak_y=static_cast<real_prec>(2.*M_PI)/L2;
    real_prec deltak_z=static_cast<real_prec>(2.*M_PI)/L3;
    
    real_prec kx= i<=N1/2 ? deltak_x*static_cast<real_prec>(i) : -deltak_x*static_cast<real_prec>(N1-i);
    real_prec ky= j<=N2/2 ? deltak_y*static_cast<real_prec>(j) : -deltak_y*static_cast<real_prec>(N2-j);
    real_prec kz= k<=N3/2 ? deltak_z*static_cast<real_prec>(k) : -deltak_z*static_cast<real_prec>(N3-k);

    return kx*kx+ky*ky+kz*kz;

  }

//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
 void kernelcomp(real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3,ULONG N1, ULONG N2, ULONG N3, real_prec smol, int filtertype, string output_dir)
 {

#ifdef _USE_OMP_
   int NTHREADS=_NTHREADS_;
   omp_set_num_threads(NTHREADS);
#endif

#ifdef _FULL_VERBOSE_
   cout<<YELLOW<<"Computing kernel  "<<__PRETTY_FUNCTION__<<RESET<<endl;
#endif
   string fnameR=output_dir+"kernel"+to_string(N1)+"V"+to_string(L1)+"r"+to_string(smol);
   
   ifstream inStream;
   string ffnn=fnameR+".dat";
   inStream.open(ffnn.data());

   if (inStream.is_open() == false )
    {
      bool gauss=false;
      bool errfunc=false;
      bool tophat=false;
       
      switch (filtertype)
     	 {
 	      case 1:
	      gauss=true;
	      break;
	      case 2:
        tophat=true;
	      break;
	      case 3:
	      errfunc=true;
	      break;
  	   }
    
      ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
      ULONG Nhalf=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3/2+1);
      vector<real_prec> out(Nhalf,0);

#ifdef DOUBLE_PREC
      complex_prec *AUX= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
      complex_prec *AUX= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif
     real_prec asmth=num_1;
     real_prec u; 
     real_prec rS=smol;
     real_prec rS2=rS*rS;
     real_prec kcut=smol;//2.*M_PI/rS;
     real_prec sigma=static_cast<real_prec>(.3);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for (ULONG i=0;i<N1;i++)
	 for (ULONG j=0;j<N2;j++)
	   for (ULONG k=0;k<N3/2+1;k++)
	     {

	       ULONG ii=index_3d(i,j,k,N2,N3/2+1);
	       real_prec k2=k_squared(i,j,k,L1,L2,L3,N1,N2,N3/2+1);
	       
	       if (tophat==true)
		 {
		   u = sqrt(k2);
		   
		   if (u>kcut)
		     AUX[ii][REAL]=0.0;
		   else
		     AUX[ii][REAL]=num_1;
		 }
	       
	       if (errfunc==true)
		 {
		   u = static_cast<real_prec>((sqrt(k2)-kcut)/(sqrt(2.)*sigma));
		   real_prec fac = static_cast<real_prec>(erfc(u));
		   AUX[ii][REAL]=fac;
		 }
	       
	       if (gauss==true)
		 AUX[ii][REAL]=static_cast<real_prec>(exp(-k2*rS2/2.));
	       
	       
	       AUX[ii][IMAG]=0.0;

	     }
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(ULONG i=0;i<Nhalf;i++)
   	   out[i]=AUX[i][REAL];

       vector<real_prec> aux(N,0);
       do_fftw_c2r(N1, AUX, aux);
	
        real_prec wtotD=0.;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:wtotD)
#endif
        for(ULONG i=0;i<N;i++)
          wtotD+=static_cast<real_prec>(aux[i]);

	aux.clear();
	aux.shrink_to_fit();

	// Normalize kernel in Fourier space
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
        for(ULONG i=0;i<Nhalf;i++)
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
          out[i]/=static_cast<real_prec>(wtotD);
	
	
	// este out debe ser en Fourier espace

       dump_scalar(out,N1,N2,N3/2+1,0,fnameR);

#ifdef DOUBLE_PREC
       fftw_free(AUX);
#else
       fftwf_free(AUX);
#endif

       }

 }
 

 //##################################################################################
 //##################################################################################
 //##################################################################################
 //##################################################################################
 void kernelcomp_for_boost(real_prec L, real_prec d,ULONG N, vector<real_prec>&kernel,string output_dir)
 {
#ifdef _USE_OMP_
   int NTHREADS=_NTHREADS_;
   omp_set_num_threads(NTHREADS);
#endif

#ifdef _FULL_VERBOSE_
   cout<<YELLOW<<"Computing kernel from ratio of HF/Linear power spectrum"<<RESET<<endl;
#endif
   string fnameR=output_dir+"kernel_boost"+to_string(N);

   ifstream inStream;
   string ffnn=fnameR+".dat";
   inStream.open(ffnn.data());

   FileOutput File;
   vector<real_prec> prop;
   string file_lin="../Output/linear_matter_power_spectrum_EH_redshift_0.000000.txt";
   int nlines=File.read_file(file_lin,prop,1);

   vector<real_prec> l_power(nlines,0);
#pragma omp parallel for
   for(int i=0;i<l_power.size();++i)
       l_power[i]=prop[1+i*2];
   prop.clear(); prop.shrink_to_fit();

   string file_nlin="../Output/non_linear_matter_power_spectrum_EH_halo_fit_redshift_0.000000.txt";
   nlines=File.read_file(file_nlin,prop,1);
   vector<gsl_real> nl_power(nlines,0);
   vector<gsl_real> kvp(nlines,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int i=0;i<nlines;++i)
      kvp[i]=prop[i*2];
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int i=0;i<nlines;++i)
      nl_power[i]=prop[1+i*2];
   prop.clear(); prop.shrink_to_fit();

    for(int i=0;i<nl_power.size();++i)
      nl_power[i]=(nl_power[i]/l_power[i]);

    l_power.clear(); l_power.shrink_to_fit();

    vector<real_prec> coords(N,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<coords.size() ;++i)
      coords[i]= (i<=N/2? static_cast<real_prec>(i): -static_cast<real_prec>(N-i));

    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, nl_power.size());
    gsl_spline_init (spline, &kvp[0], &nl_power[0], nl_power.size());

#ifdef _FULL_VERBOSE_
    cout<<YELLOW<<"Computing kernel"<<RESET<<endl;
#endif

    real_prec delta=2.0*M_PI/static_cast<real_prec>(L);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for (ULONG i=0;i<N;i++)
         for (ULONG j=0;j<N;j++)
           for (ULONG k=0;k<N/2+1;k++)
             {
              real_prec kv=  sqrt(pow(coords[i],2)+pow(coords[j],2)+pow(coords[k],2));
              int kmod=static_cast<int>(floor(kv));
              if(kmod<N/2 && kmod >0)
                kernel[index_3d(i,j,k,N,N/2+1)]=static_cast<real_prec>(gsl_spline_eval (spline,kv*delta, acc));

          }

#ifdef _FULL_VERBOSE_
   cout<<BOLDGREEN<<"                                               ["<<BOLDBLUE<<"DONE"<<BOLDGREEN<<"]"<<RESET<<endl;
#endif

 }



 //##################################################################################
 //##################################################################################
 //##################################################################################
 //##################################################################################

real_prec calc_kx(ULONG i,real_prec L1,ULONG N1)
{
  real_prec kfac=static_cast<real_prec>(2.*M_PI/L1);
  real_prec k1=0.;
  
  if (i<=N1/2)
    k1 = kfac*static_cast<real_prec>(i);
  else
    k1 = -kfac*static_cast<real_prec>(N1-i);		
  
  return(k1);
}

//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################

 
 real_prec calc_ky(ULONG j,real_prec L2,ULONG N2)
 {
   real_prec  kfac=static_cast<real_prec>(2.*M_PI/L2);
   real_prec k2=0.;
   
   if (j<=N2/2)
    k2 = kfac*static_cast<real_prec>(j);
  else
    k2 = -kfac*static_cast<real_prec>(N2-j);		
  
  return(k2);
}

//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################

real_prec calc_kz(ULONG k,real_prec L3,ULONG N3)
{
  real_prec  kfac=static_cast<real_prec>(2.*M_PI/L3);
  real_prec k3=0.;
  
  if (k<=N3/2)
    k3 = kfac*static_cast<real_prec>(k);
  else
    k3 = -kfac*static_cast<real_prec>(N3-k);		
  
  return(k3);
}

 //##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
 void calc_twolptterm(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, vector<real_prec>&phiv, vector<real_prec> &m2v)
{  

#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif


  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
  ULONG Nhalf=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3/2+1);
  
  string fnamef;


#ifndef COMPATCTWOLPT
  #ifdef DOUBLE_PREC
  complex_prec *philv= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *philv= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif


  //  fftw_array<complex_prec> philv(Nhalf);
  
#pragma omp parallel for
  for(ULONG i=0;i<Nhalf;i++)
    {
      philv[i][REAL]=0.0;
      philv[i][IMAG]=0.0;
    }

  
  do_fftw_r2c(N1,phiv,philv);  
   
  vector<real_prec> LapPhivx(N,0), LapPhivy(N,0), LapPhivz(N,0);
  vector<real_prec> LapPhivxy(N,0), LapPhivxz(N,0), LapPhivyz(N,0);
  
#ifdef  GFFT 
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,LapPhivx,1,1);
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,LapPhivy,2,2);
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,LapPhivz,3,3);
	      
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,LapPhivxy,1,2);
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,LapPhivxz,1,3);
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,LapPhivyz,2,3);
#endif


#if defined(GFINDIFF) || defined (GFFT)
  vector<real_prec> dummy(N,0);
#endif



#ifdef GFINDIFF  
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,1);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,LapPhivx,1);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,LapPhivxy,2);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,LapPhivxz,3);

  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,2);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,LapPhivy,2);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,LapPhivyz,3);

  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,3);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,LapPhivz,3);
#endif


#pragma omp parallel for
  for(ULONG i=0;i<N;i++)
     m2v[i]=LapPhivx[i]*LapPhivy[i]-LapPhivxy[i]*LapPhivxy[i]+LapPhivx[i]*LapPhivz[i]-LapPhivxz[i]*LapPhivxz[i]+LapPhivy[i]*LapPhivz[i]-LapPhivyz[i]*LapPhivyz[i];

#else 

#ifdef GFFT

  //  LapPhivx[i]*LapPhivy[i]
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,m2v,1,1);    
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,dummy,2,2);
  fnamef=string("aux2");
  dump_scalar(dummy,N1,N2,N3,0,fnamef);	
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=dummy[i];
  
  // -LapPhivxy[i]*LapPhivxy[i]
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,dummy,1,2);
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i]*dummy[i];

  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	
  
  // +LapPhivx[i]*LapPhivz[i]
  {
    string fname=string("aux1");
    get_scalar(fname,dummy,N1,N2,N3);
  }
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,m2v,3,3);
  fnamef=string("aux3");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=dummy[i];
  {
    string fname=string("aux");
    get_scalar(fname,dummy,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];
  
  // -LapPhivxz[i]*LapPhivxz[i]
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,dummy,1,3);
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i]*dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	
  
  // +LapPhivy[i]*LapPhivz[i]
  {
    string fname=string("aux2");
    get_scalar(fname,dummy,N1,N2,N3);
  }
  
  {
    string fname=string("aux3");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy[i]*=m2v[i];
  {
    string fname=string("aux");
    get_scalar(fname,m2v,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];
  
  // -LapPhivyz[i]*LapPhivyz[i]  
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,dummy,2,3);
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i]*dummy[i];
#endif  
  
  
#ifdef GFINDIFF  
  //  LapPhivx[i]*LapPhivy[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,1);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,m2v,1);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,2);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,2);
  fnamef=string("aux2");
  dump_scalar(dummy,N1,N2,N3,0,fnamef);	
  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	

  // -LapPhivxy[i]*LapPhivxy[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,1);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,2);
  {
    string fname=string("aux");
    get_scalar(fname,m2v,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i]*dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	
  
  // +LapPhivx[i]*LapPhivz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,3);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,m2v,3);
  fnamef=string("aux3");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	
  {
    string fname=string("aux1");
    get_scalar(fname,dummy,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=dummy[i];
  {
    string fname=string("aux");
    get_scalar(fname,dummy,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	

  // -LapPhivxz[i]*LapPhivxz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,1);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,3);
  {
    string fname=string("aux");
    get_scalar(fname,m2v,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i]*dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	

  // +LapPhivy[i]*LapPhivz[i]
  {
    string fname=string("aux2");
    getscalar(fname,dummy,N1,N2,N3);
  }
  {
    string fname=string("aux3");
    getscalar(fname,m2v,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy[i]*=m2v[i];
  {
    string fname=string("aux");
    getscalar(fname,m2v,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	

  // -LapPhivyz[i]*LapPhivyz[i]  
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,2);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,3);
  {
    string fname=string("aux");
    getscalar(fname,m2v,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i]*dummy[i];
#endif 
#endif



#ifdef DOUBLE_PREC
       fftw_free(philv);
#else
       fftwf_free(philv);
#endif

}

 //##################################################################################
//##################################################################################


 void calc_LapPhiv(ULONG N1,ULONG N2,ULONG N3,real_prec L1,real_prec L2,real_prec L3, complex_prec *philv,vector<real_prec>&LapPhiv,int index1,int index2)
{


#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
  ULONG Nhalf=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3/2+1);	
#ifdef DOUBLE_PREC
  complex_prec *LapPhivl= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *LapPhivl= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif


  real_prec k1=0.;
  real_prec k2=0.;

#ifdef FORPAR
#pragma omp parallel for
#endif  
  for (ULONG i=0;i<N1;i++)
    for (ULONG j=0;j<N2;j++)
      for (ULONG k=0;k<=N3/2;k++)
	{	  
	  ULONG ihalf=k+(N3/2+1)*(j+N2*i);

	  real_prec kx=calc_kx(i,L1,N1);
	  real_prec ky=calc_ky(j,L2,N2);
	  real_prec kz=calc_kz(k,L3,N3);

	  switch (index1)
	    {	      
	    case 1:
	      k1=kx;
	      break;
	    case 2:
	      k1=ky;
	      break;
	    case 3:
	      k1=kz;
	      break;
	    }

	  
	  switch (index2)
	    {	      
	    case 1:
	      k2=kx;
	      break;
	    case 2:
	      k2=ky;
	      break;
	    case 3:
	      k2=kz;
	      break;
	    }

	  LapPhivl[ihalf][REAL]=-k1*k2*(philv[ihalf][REAL]);	
	  LapPhivl[ihalf][IMAG]=-k1*k2*(philv[ihalf][IMAG]);
	}

  do_fftw_c2r(N1,LapPhivl,LapPhiv);

}




 //##################################################################################
//##################################################################################


 void calc_curlcomp(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, const vector<real_prec> &phiv, vector<real_prec> phiv2, vector<real_prec> &m2v, int comp)
{  

#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif


  int pcomp, mcomp;

  switch(comp)
    {
    case 1:
      {
	pcomp=2;
	mcomp=3;
	break;
      }
    case 2:
      {
	pcomp=1;
	mcomp=3;
	break;
      }
    case 3:
      {
	pcomp=1;
	mcomp=2;
	break;
      }
    }

  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
 
  string fnamef;     

  vector<real_prec> dummy(N,0);

#define num_0_2_5 static_cast<real_prec>(0.25)
  
#ifdef GFINDIFF  
  //LapPhiv1l[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,comp);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,L1,L2,L3,0,fnamef);	
  //LapPhiv2pp[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,pcomp);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,pcomp);

  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }

  //1/4*LapPhiv1l[i]*LapPhiv2pp[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_2_5*dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	

  //LapPhivl[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,comp);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv2mm[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,mcomp);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,mcomp);
 
  {
    string fname=string("aux1");
    getscalar(fname,m2v,N1,N2,N3);
  }
  //1/4*LapPhiv1l[i]*LapPhiv2mm[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_2_5*dummy[i];

  {
    string fname=string("aux");
    getscalar(fname,dummy,N1,N2,N3);
  }
 
  //1/4*LapPhiv1l[i]*LapPhiv2pp[i]+1/4*LapPhiv1l[i]*LapPhiv2mm[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];
  fnamef=string("auxs");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
 
  /////

  //LapPhiv2l[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,comp);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv1pp[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,pcomp);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,pcomp);

  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }

  //1/4*LapPhiv2l[i]*LapPhiv1pp[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_2_5*dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	

  //LapPhiv2l[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,comp);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv1mm[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,mcomp);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,mcomp);
 
  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //1/4*LapPhiv2l[i]*LapPhiv1mm[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_2_5*dummy[i];

  {
    string fname=string("auxs");
    get_scalar(fname,dummy,N1,N2,N3);
  }
 
  //1/4*LapPhiv1l[i]*LapPhiv2pp[i]+1/4*LapPhiv1l[i]*LapPhiv2mm[i]
  //-1/4*LapPhiv2l[i]*LapPhiv1pp[i]+1/4*LapPhiv2l[i]*LapPhiv1mm[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i];
  fnamef=string("auxs");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
 
  /////

  //LapPhiv1p[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,pcomp);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv2lp[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,comp);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,pcomp);

  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }

  //1/4*LapPhiv1p[i]*LapPhiv2lp[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_2_5*dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	

  //LapPhiv1m[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,mcomp);
  fnamef=string("aux1");
  // dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv2lm[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,comp);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,mcomp);
 
  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //1/4*LapPhiv1m[i]*LapPhiv2lm[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_2_5*dummy[i];

  {
    string fname=string("auxs");
    get_scalar(fname,dummy,N1,N2,N3);
  }
 
  // 1/4*LapPhiv1l[i]*LapPhiv2pp[i]+1/4*LapPhiv1l[i]*LapPhiv2mm[i]
  //-1/4*LapPhiv2l[i]*LapPhiv1pp[i]+1/4*LapPhiv2l[i]*LapPhiv1mm[i]
  //-1/4*LapPhiv1l[i]*LapPhiv2lp[i]+1/4*LapPhiv1l[i]*LapPhiv2lm[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i];
  fnamef=string("auxs");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
 
  /////

  //LapPhiv2p[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,pcomp);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv1lp[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,comp);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,pcomp);

  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }

  //1/4*LapPhiv2p[i]*LapPhiv1lp[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_2_5*dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	

  //LapPhiv2m[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,mcomp);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv1lm[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,comp);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,mcomp);
 
  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //1/4*LapPhiv2m[i]*LapPhiv1lm[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_2_5*dummy[i];

  {
    string fname=string("auxs");
    get_scalar(fname,dummy,N1,N2,N3);
  }
 
  // 1/4*LapPhiv1l[i]*LapPhiv2pp[i]+1/4*LapPhiv1l[i]*LapPhiv2mm[i]
  //-1/4*LapPhiv2l[i]*LapPhiv1pp[i]+1/4*LapPhiv2l[i]*LapPhiv1mm[i]
  //-1/4*LapPhiv1l[i]*LapPhiv2lp[i]+1/4*LapPhiv1l[i]*LapPhiv2lm[i]
  //+1/4*LapPhiv2l[i]*LapPhiv1lp[i]+1/4*LapPhiv2l[i]*LapPhiv1lm[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];
 
#endif 
}



 //##################################################################################
//##################################################################################


 void calc_mu2term(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, const vector<real_prec> &phiv, vector<real_prec> phiv2, vector<real_prec>&m2v)
{  

#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
 
  string fnamef;     

  vector<real_prec> dummy(N);

#ifdef GFINDIFF  
  //LapPhiv1xx[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,1);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,m2v,1);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv2yy[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,2);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,2);

  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }

  //0.5*LapPhiv1xx[i]*LapPhiv2yy[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_5*dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	

  //LapPhiv2xx[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,2);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,m2v,2);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv1yy[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,1);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,1);
 
  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //0.5*LapPhiv2xx[i]*LapPhiv1yy[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_5*dummy[i];

  {
    string fname=string("aux");
    get_scalar(fname,dummy,N1,N2,N3);
  }
 
  //0.5*LapPhiv1xx[i]*LapPhiv2yy[i]+0.5*LapPhiv2xx[i]*LapPhiv1yy[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
 
  /////

  //LapPhiv1xx[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,1);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,m2v,1);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv2zz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,3);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,3);

  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //0.5*LapPhiv1xx[i]*LapPhiv2zz[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_5*dummy[i];
  fnamef=string("auxb");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	

  //LapPhiv2xx[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,3);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,m2v,3);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv1zz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,1);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,1);

  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //0.5*LapPhiv2xx[i]*LapPhiv1zz[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_5*dummy[i];

 {
    string fname=string("auxb");
    get_scalar(fname,dummy,N1,N2,N3);
  }
 
  //0.5*LapPhiv1xx[i]*LapPhiv2zz[i]+0.5*LapPhiv2xx[i]*LapPhiv1zz[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];  

  {
    string fname=string("aux");
    get_scalar(fname,dummy,N1,N2,N3);
  }
 
  // 0.5*LapPhiv1xx[i]*LapPhiv2yy[i]+0.5*LapPhiv2xx[i]*LapPhiv1yy[i]
  //+0.5*LapPhiv1xx[i]*LapPhiv2zz[i]+0.5*LapPhiv2xx[i]*LapPhiv1zz[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	


  ////

  //LapPhiv1yy[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,2);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,m2v,2);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv2zz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,3);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,3);
  fnamef=string("aux2");
  dump_scalar(dummy,N1,N2,N3,9,fnamef);	

  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //0.5*LapPhiv1yy[i]*LapPhiv2zz[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_5*dummy[i];
  fnamef=string("auxb");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	

  //LapPhiv2yy[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,3);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,m2v,3);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv1zz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,2);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,2);
  fnamef=string("aux2");
  dump_scalar(dummy,N1,N2,N3,9,fnamef);	

  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //0.5*LapPhiv2yy[i]*LapPhiv1zz[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_5*dummy[i];

  {
    string fname=string("auxb");
    getscalar(fname,dummy,N1,N2,N3);
  }
 
  //0.5*LapPhiv1yy[i]*LapPhiv2zz[i]+0.5*LapPhiv2yy[i]*LapPhiv1zz[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];  

  {
    string fname=string("aux");//!!!
    getscalar(fname,dummy,N1,N2,N3);
  }
 
  // 0.5*LapPhiv1xx[i]*LapPhiv2yy[i]+0.5*LapPhiv2xx[i]*LapPhiv1yy[i]
  //+0.5*LapPhiv1xx[i]*LapPhiv2zz[i]+0.5*LapPhiv2xx[i]*LapPhiv1zz[i]
  //+0.5*LapPhiv1yy[i]*LapPhiv2zz[i]+0.5*LapPhiv2yy[i]*LapPhiv1zz[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	


  ////


  // LapPhiv1xy[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,1);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,2);
  fnamef=string("aux1");
  dump_scalar(dummy,N1,N2,N3,9,fnamef);	
  // LapPhiv2xy[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,1);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,2);
 
  {
    string fname=string("aux1");
    getscalar(fname,m2v,N1,N2,N3);
  }

 // LapPhiv1xy[i]*LapPhiv2xy[i] 
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy[i]*=m2v[i];

  {
    string fname=string("aux");
    getscalar(fname,m2v,N1,N2,N3);
  }
 
  // 0.5*LapPhiv1xx[i]*LapPhiv2yy[i]+0.5*LapPhiv2xx[i]*LapPhiv1yy[i]
  //+0.5*LapPhiv1xx[i]*LapPhiv2zz[i]+0.5*LapPhiv2xx[i]*LapPhiv1zz[i]
  //+0.5*LapPhiv1yy[i]*LapPhiv2zz[i]+0.5*LapPhiv2yy[i]*LapPhiv1zz[i]
  //-LapPhiv1xy[i]*LapPhiv2xy[i] 
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  
  ////


  // LapPhiv1xz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,1);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,3);
  fnamef=string("aux1");
  dump_scalar(dummy,N1,N2,N3,9,fnamef);	
  // LapPhiv2xz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,1);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,3);
 
  {
    string fname=string("aux1");
    getscalar(fname,m2v,N1,N2,N3);
  }

 // LapPhiv1xz[i]*LapPhiv2xz[i] 
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy[i]*=m2v[i];

  {
    string fname=string("aux");
    getscalar(fname,m2v,N1,N2,N3);
  }
 
  // 0.5*LapPhiv1xx[i]*LapPhiv2yy[i]+0.5*LapPhiv2xx[i]*LapPhiv1yy[i]
  //+0.5*LapPhiv1xx[i]*LapPhiv2zz[i]+0.5*LapPhiv2xx[i]*LapPhiv1zz[i]
  //+0.5*LapPhiv1yy[i]*LapPhiv2zz[i]+0.5*LapPhiv2yy[i]*LapPhiv1zz[i]
  //-LapPhiv1xy[i]*LapPhiv2xy[i] 
  //-LapPhiv1xz[i]*LapPhiv2xz[i] 
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  
  ////


  // LapPhiv1yz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,2);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,3);
  fnamef=string("aux1");
  dump_scalar(dummy,N1,N2,N3,9,fnamef);	
  // LapPhiv2yz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,2);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,3);
 
  {
    string fname=string("aux1");
    getscalar(fname,m2v,N1,N2,N3);
  }

 // LapPhiv1yz[i]*LapPhiv2yz[i] 
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy[i]*=m2v[i];

  {
    string fname=string("aux");

    getscalar(fname,m2v,N1,N2,N3);
  }
 
  // 0.5*LapPhiv1xx[i]*LapPhiv2yy[i]+0.5*LapPhiv2xx[i]*LapPhiv1yy[i]
  //+0.5*LapPhiv1xx[i]*LapPhiv2zz[i]+0.5*LapPhiv2xx[i]*LapPhiv1zz[i]
  //+0.5*LapPhiv1yy[i]*LapPhiv2zz[i]+0.5*LapPhiv2yy[i]*LapPhiv1zz[i]
  //-LapPhiv1xy[i]*LapPhiv2xy[i] 
  //-LapPhiv1xz[i]*LapPhiv2xz[i] 
  //-LapPhiv1yz[i]*LapPhiv2yz[i] 
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i];
 

#endif 
}



 //##################################################################################
 //##################################################################################
 
 void calc_Det(ULONG N1, ULONG N2, ULONG N3, real_prec L1, real_prec L2, real_prec L3, const vector<real_prec>&in, vector<real_prec> &out)
 {  

#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif


  ULONG N=N1*N2*N3;

  string fname;
  vector<real_prec> dummy(N,0), dummy2(N,0);

  //phi,11
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,1);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,dummy2,1);
  fname=string("aux1");
  //  dump_scalar(dummy2,N1,N2,N3,9,fname);
	
  //phi,22
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,2);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,dummy2,2);
	
  //phi,33
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,3);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,out,3);
 
  //phi,22*phi,33
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=out[i];
 
  fname=string("aux1");
  get_scalar(fname,out,N1,N2,N3);
  //phi,11*phi,22*phi,33
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=out[i];
  
  fname=string("auxs");
  //  dump_scalar(dummy2,N1,N2,N3,9,fname);
  
  //phi,23
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,2);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,dummy2,3);
  
  //phi,23*phi,23
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=dummy2[i];
 
  fname=string("aux1");
  get_scalar(fname,out,N1,N2,N3);
  //phi,11*phi,23*phi,23
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=out[i];
  
  fname=string("auxs");
  get_scalar(fname,out,N1,N2,N3);
  //phi,11*phi,22*phi,33-phi,11*phi,23*phi,23
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    out[i]-=dummy2[i];
  //  dump_scalar(out,N1,N2,N3,9,fname);
  
  ////
  
  //phi,12
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,1);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,dummy2,2);
  fname=string("aux1");
  //  dump_scalar(dummy2,N1,N2,N3,9,fname);
  
  //phi,12*phi,12
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=dummy2[i];
  
  //phi,33
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,3);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,out,3);
  
  //phi,12*phi,12*phi,33
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=out[i];
 
  fname=string("aux2");
  //  dump_scalar(dummy2,N1,N2,N3,9,fname);

  //phi,23
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,2);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,dummy2,3);

  //phi,13
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,1);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,out,3);

  //phi,23*phi,13
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=out[i];

  fname=string("aux1");
  get_scalar(fname,out,N1,N2,N3);

  //phi,12*phi,23*phi,13
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=num_2*out[i];//this term appears twice in the determinant
 

  fname=string("aux2");
  get_scalar(fname,out,N1,N2,N3);
  //-phi,12*phi,12*phi,33+phi,12*phi,23*phi,13//sign must be changed below
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]-=out[i];

  fname=string("auxs");
  get_scalar(fname,dummy2,N1,N2,N3);
  //    phi,11*phi,22*phi,33-phi,11*phi,23*phi,23
  //-(2*phi,12*phi,12*phi,33-phi,12*phi,23*phi,13)
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]+=out[i];//here comes a plus sign see above
  //  dump_scalar(dummy2,N1,N2,N3,9,fname);


  ////

  //phi,13
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,1);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,dummy2,3);
  fname=string("aux1");
  //  dump_scalar(dummy2,N1,N2,N3,9,fname);
	
  //phi,13*phi,13
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=dummy2[i];

  //phi,22
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,2);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,out,2);

  //phi,13*phi,13*phi,22
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=out[i]; 
  
  fname=string("auxs");
  get_scalar(fname,out,N1,N2,N3);
  //  phi,11*phi,22*phi,33-phi,11*phi,23*phi,23
  //-(phi,12*phi,12*phi,33-phi,12*phi,23*phi,13)
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    out[i]-=dummy2[i];
  
}



//##################################################################################
//##################################################################################

 void convcomp(real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3,ULONG N1, ULONG N2, ULONG N3,  vector<real_prec>&in, vector<real_prec> &out, int filtertype,real_prec smol, string file_kernel)
{


#ifdef _FULL_VERBOSE_
  cout<<YELLOW<<"Convolution"<<RESET<<endl;
#endif
  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
  ULONG Nhalf=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3/2+1);  	
#ifdef DOUBLE_PRC
  complex_prec *AUX= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *AUX= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<Nhalf;++i)
    {
      AUX[i][REAL]=0;
      AUX[i][IMAG]=0;
    }


  // Convert input field to Fourier space

  //read the kernel in 3D
  string fname=file_kernel+"kernel"+to_string(N1)+"V"+to_string(L1)+"r"+to_string(smol);

  // read kernel in Fourier space
  vector<real_prec> kern(Nhalf,0);
  get_scalar(fname,kern,N1,N2,N3/2+1);
   

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG i=0;i<Nhalf;i++)
    {
      real_prec cor=kern[i];
      AUX[i][REAL]*=cor;
      AUX[i][IMAG]*=cor;
    }

  do_fftw_r2c(N1,in,AUX);

  kern.clear();
  kern.shrink_to_fit();

  do_fftw_c2r(N1,AUX,out);

#ifdef DOUBLE_PRC
  fftw_free(AUX);
#else
  fftwf_free(AUX);
#endif


 }

 //##################################################################################
 //#################################################################################
 void convolvek(ULONG N1, vector<real_prec>&in, vector<real_prec> &kernel, vector<real_prec> &out)
{

#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

#ifdef _FULL_VERBOSE_
  cout<<YELLOW<<"Convolution"<<RESET<<endl;
#endif
  ULONG N=in.size();
  ULONG Nhalf=kernel.size();
#ifdef DOUBLE
  complex_prec *AUX= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *AUX= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<Nhalf;++i)
    {
      AUX[i][REAL]=0;
      AUX[i][IMAG]=0;
    }

  // Convert input field to Fouroer space
  do_fftw_r2c(N1,in,AUX);
   
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG i=0;i<Nhalf;i++)
    {
      real_prec cor=kernel[i];
      AUX[i][REAL]*=cor;
      AUX[i][IMAG]*=cor;
      }

  do_fftw_c2r(N1,AUX,out);
 
#ifdef DOUBLE
  fftw_free(AUX);
#else
  fftwf_free(AUX);
#endif

#ifdef _FULL_VERBOSE_
   cout<<BOLDGREEN<<"                                               ["<<BOLDBLUE<<"DONE"<<BOLDGREEN<<"]"<<RESET<<endl;
#endif

}




 //##################################################################################
 //#################################################################################
 void give_power(real_prec Lbox, ULONG N1, vector<real_prec>&in, vector<real_prec> &out)
{

#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif


  ULONG N=in.size();
  ULONG Nhalf=static_cast<ULONG>(N1)*static_cast<ULONG>(N1)*static_cast<ULONG>(N1/2+1);
#ifdef DOUBLE_PREC
  complex_prec *AUX= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *AUX= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif

#pragma omp parallel for
  for(ULONG i=0;i<Nhalf;++i)
      AUX[i][REAL]=0;
#pragma omp parallel for
  for(ULONG i=0;i<Nhalf;++i)
      AUX[i][IMAG]=0;

  // Convert input field to Fouroer space
  do_fftw_r2c(N1,in,AUX);

  real_prec deltak=2.*M_PI/Lbox;

  vector<real_prec> coords(N1,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<N1 ;++i)
    coords[i]=deltak*(i<=N1/2? static_cast<real_prec>(i): -static_cast<real_prec>(N1-i));


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i< N1 ; i++)
    for(ULONG j=0;j< N1 ; j++)
      for(ULONG k=0;k< N1/2+1; k++)
        {
          ULONG ind=index_3d(i,j,k, N1, N1/2+1);
          real_prec kmod2=pow(coords[i],2)+pow(coords[j],2)+pow(coords[k],2);  // Get k**2
          real_prec f_k =0.;
          if(kmod2>0.)
//            f_k=exp(0.5*pow(kmod2, 0.2));
             f_k=2.*(1.0+tanh(2.*kmod2));

          AUX[ind][REAL]*=f_k;
          AUX[ind][IMAG]*=f_k;
        }
  coords.clear();
  coords.shrink_to_fit();

  do_fftw_c2r(N1,AUX,out);
#ifdef DOUBLE_PREC
  fftw_free(AUX);
#else
  fftwf_free(AUX);
#endif

 }


//##################################################################################
//########################################################################### 
 real_prec linearvel3d(int index, real_prec kx, real_prec ky, real_prec kz, real_prec phi)
{
  real_prec out, kl;

  switch(index)
    {
    case 1:
      {
	kl=kx;
	break;
      }
    case 2:
      {
	kl=ky;
	break;
      }
    case 3:
      {
	kl=kz;
	break;
      }
    }
  real_prec kmod2=kx*kx+ky*ky+kz*kz;

  real_prec fackern=0.0;
  if (kmod2>eps)
    fackern = kl/kmod2;
  out=fackern*phi;

  return out;
}



 //##################################################################################
 //#################################################################################
 int get_bin(real_prec x, real_prec xmin, int nbins, real_prec delta, bool accumulate)
 {
   int ibin  =  static_cast<int>(floor((static_cast<double>(x) - static_cast<double>(xmin))/static_cast<double>(delta)));
   if(ibin==nbins)ibin--;
   if(true==accumulate)
     {
       if (ibin>nbins)ibin=nbins-1; // All values above the maximum are placed in the last bin
       if (ibin<0)ibin=0; // All values below the minimum are placed in the first bin
     }
   return ibin;
 }


 //##################################################################################
//##################################################################################

real_prec tidal_anisotropy(real_prec lambda1, real_prec lambda2, real_prec lambda3)
{
//  real_prec tidal= (sqrt(0.5*(pow(lambda3-lambda1,2)+pow(lambda3-lambda2,2)+pow(lambda2-lambda1,2)))/(1.+lambda1+lambda2+lambda3));
  real_prec tidal= (sqrt(0.5*(pow(lambda3-lambda1,2)+pow(lambda3-lambda2,2)+pow(lambda2-lambda1,2)))/(1.));
  return tidal;

}

//##################################################################################
//##################################################################################
real_prec invariant_field_II(real_prec lambda1, real_prec lambda2, real_prec lambda3)
{
#ifdef _USE_EIGENVALUES_
    real_prec inv= lambda1;
#else
    real_prec inv= lambda1*lambda2 +lambda2*lambda3 + lambda1*lambda3;
#endif
    return inv;
 }

//##################################################################################
//##################################################################################

real_prec invariant_field_III(real_prec lambda1, real_prec lambda2, real_prec lambda3)
{
#ifdef _USE_EIGENVALUES_
    real_prec inv= lambda2;
#else
  real_prec inv= (lambda1*lambda2)*lambda3;
#endif
  return inv;

}

//##################################################################################
//##################################################################################

real_prec invariant_field_I(real_prec lambda1, real_prec lambda2, real_prec lambda3)
{
  real_prec inv= lambda1+lambda2+lambda3;
  return inv;
}


//##################################################################################
//##################################################################################

real_prec invariant_field_IV(real_prec lambda1, real_prec lambda2, real_prec lambda3)
{
#ifdef _USE_EIGENVALUES_
    real_prec inv= lambda3;
#else
//  real_prec inv=  pow(lambda1,3) + pow(lambda2,3) + pow(lambda3,3);
  real_prec inv=  pow(lambda1,2) + pow(lambda2,2) + pow(lambda3,2);
//  real_prec inv=  (lambda1+lambda2+lambda3)*(lambda1*lambda2 +lambda2*lambda3 + lambda1*lambda3)  ;// pow(lambda1,2) + pow(lambda2,2) + pow(lambda3,2);
#endif
  return inv;
}


//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################

real_prec ellipticity(real_prec lambda1, real_prec lambda2, real_prec lambda3)
{
  real_prec inv= (lambda1-lambda3)/(2.*(lambda1+lambda2+lambda3));
  return inv;

}

//##################################################################################
//##################################################################################

real_prec prolat(real_prec lambda1, real_prec lambda2, real_prec lambda3)
{
  real_prec inv= (lambda1+lambda3-2.*lambda2)/(2.*(lambda1+lambda2+lambda3));
  return inv;

}

//##################################################################################
//##################################################################################

 void sort_2vectors(vector<vector<ULONG> >& v1,vector< vector<ULONG> >&v2){
   // v1 is to be sorted. The elements of the other vectors
   // are shuffled accordingly.
   

#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif


   ULONG i, j;
   ULONG NN=v1.size();
   ULONG n=v1[0].size();

   for (i=0;i<NN;++i)
     {
       vector<ULONG>iwksp(n,0);
       vector<float>wksp(n,0);
       vector<ULONG>av1(n,0);
       for (j=0;j<n;++j) av1[j]=v1[i][j];
       indexx_ulong(av1,iwksp); //av1 must be int
       
       for (j=0;j<n;++j) wksp[j]=av1[j];
       for (j=0;j<n;++j) av1[j]=wksp[iwksp[j]];
       for (j=0;j<n;++j) v1[i][j]=av1[j];
       
       for (j=0;j<n;++j) wksp[j]=v2[i][j];
       for (j=0;j<n;++j) v2[i][j]=wksp[iwksp[j]];
     }
 }
 //##################################################################################
 //##################################################################################
 void sort_1dvectors(vector<ULONG> & v1, vector<ULONG> &v2){
   // v1 is to be sorted. The elements of the other vectors
   // are shuffled accordingly. 
   ULONG i, j;
   ULONG n=v1.size();
   vector<ULONG>iwksp(n,0);
   vector<ULONG>wksp(n,0);
   vector<ULONG>av1(n,0);
   // Prepare working space
   for (j=0;j<n;++j) av1[j]=v1[j];

  // sort v1. iwksp contains the order of

   //v1= [ 10 , 8, 25, 54, 1] -> v1=  returns the smae and
  //        0   1   2   3  4    iwksp= [4,1,0,2,3] is the rank
   indexx_ulong(v1,iwksp); 
    // order av1
   for (j=0;j<n;++j) wksp[j]=av1[j];
   for (j=0;j<n;++j) v1[j]=wksp[iwksp[j]];
    // order v2
   for (j=0;j<n;++j) wksp[j]=v2[j];
    // allocate values of v2 sorted
   for (j=0;j<n;++j) v2[j]=wksp[iwksp[j]];
   
 }


 //##################################################################################
 void sort_1dvectors_v2(vector<ULONG> &v1, vector<ULONG> &v2, ULONG &v1cero, ULONG &v2cero){
   // v1 is to be sorted. The elements of the other vectors
   // are shuffled accordingly.
   ULONG n=v1.size();
   vector<ULONG>iwksp(n,0);
   indexx_ulong(v1,iwksp); //av1 must be int
   v1cero=v1[iwksp[0]];
   v2cero=v2[iwksp[0]];
  }

 //##################################################################################
 void sort_1dvectors_iv2(vector<int> &v1, vector<int> &v2, int &v1cero, int &v2cero){
   // v1 is to be sorted. The elements of the other vectors
   // are shuffled accordingly.
   ULONG n=v1.size();
   vector<int>iwksp(n,0);
   indexx(v1,iwksp); //av1 must be int
   v1cero=v1[iwksp[0]];
   v2cero=v2[iwksp[0]];
  }


 //##################################################################################
 void sort_1dvectors_v3(vector<ULONG> & v1, vector<ULONG> &v2,  ULONG &v2cero){
  // v1 is to be sorted. The elements of the other vectors
   // are shuffled accordingly.
  // No ordered vector is returned, only the zero element of v2 sorted according to v1
   vector<ULONG>iwksp(v1.size(),0);
   indexx_ulong(v1,iwksp); 
   v2cero=v2[iwksp[0]];
 }

