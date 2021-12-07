// ##################################################################################
// #################################################################################
// ##################################################################################
// #################################################################################
// ##################################################################################
// #################################################################################
// ##################################################################################
// #################################################################################
// COSMICATLASS
// COSMologIcal
// CATalogs for
// LArge
// Scale
// Structure
// ##################################################################################
// #################################################################################
// ##################################################################################
// #################################################################################


/** @file main_cosmicatlassu.cpp
 *
 *  @brief Generation of mock catalosgs of DM tracers
 *  based on the BAM method + ALPT.
 *  @authors Andrés Balaguera-Antolinez, Francisco-Shu Kitaura
 */
# include "../Headers/Bam.h"
# include "../Headers/FileOutput.h"
# include "../Headers/NumericalMethods.h"
# include "../Headers/massFunctions.h"
//# include "Tasks.h"

using namespace std;

// ##################################################################################
// ##################################################################################

int main(int argc, char *argv[]){

  time_t start_all;
  time(&start_all);
  char temp;
  string par_file_bam;

  int NTHREADS=omp_get_max_threads();
  omp_set_num_threads(NTHREADS);

  string logfile="logfile.log";
  ScreenOutput So(start_all, logfile);

  if(argc==1){
    So.usage(argv[0]);
    exit(1);
  }
  while((temp =  getopt(argc, argv, "hadi:b:p:c:t:u:f:")) != -1)
    {
      if(temp=='h')
        So.usage(argv[0]);

      else if (temp=='a')  // Show authors
        So.author();

      else if(temp=='i') // displays input parameters
        {
          par_file_bam = argv[2];
          Params params(par_file_bam);
          params.show_params();
        }
      else if(temp=='d') // displays input parameters
        {
          So.show_preproc();
        }
      else if(temp=='b') // Run BAM
        {
#if !defined _POWER__
          So.message(start_all);
          par_file_bam = argv[2];
          Params Par(par_file_bam);
          Bam bam(Par,true);
          bam.set_So(So);
          bam.bamrunner();
#else
          So.message_warning("_POWER_ is defined in def.h.");
         exit(0);
#endif
          }

      else if(temp=='p')   // to mesure Power spectrum
        {
          So.message(start_all);
          par_file_bam = argv[2];
          Params params(par_file_bam);
          PowerSpectrumF cPSF(params);
	        if(params._measure_cross()==true && params._input_type()!="catalog")
	          cPSF.compute_cross_power_spectrum_grid(true, params._delta_grid_file(), params._delta_grid_file2());
	        else if(params._measure_cross()==true && params._input_type()=="catalog")
	          cPSF.compute_power_spectrum(false,false);
	        else if(params._measure_cross()==false && ( params._input_type()=="density_grid" || params._input_type()=="delta_grid" ) )
	          cPSF.compute_power_spectrum_grid();
          else if(params._measure_cross()==false && params._input_type()=="catalog")
            cPSF.compute_power_spectrum(true,false);

#ifdef _USE_REDSHIFT_BINS_
          if(params._sys_of_coord_g()==0)
          {
              So.message_warning("Asking for redshift bins while providing cartessian coordinates. Check input parameter file and def.h");
              exit(0);
          }
#endif

          So.message_time(start_all);
        }
      else if(temp=='c')   // to read input tracer catalg and analyze it. 
        {

          So.message_screen("cosmicatlass running under the option -c");
#ifdef _REDSHIFT_SPACE_
          So.message_warning("Option _REDSHIFT_SPACE_ in def.h is enabled. If velocities are not to be used, please desable it to save space. Code ends here.");
          So.message_screen("After editing the def.h file, please do make clean; make -j bam;");
          exit(0);
#endif
          So.message(start_all);
          par_file_bam = argv[2];
          Params params(par_file_bam);
          Bam bam(params,false);
          Catalog catalog(params); // this is redundant
          catalog.analyze_cat();
          So.message_time(start_all);
        }

      else if(temp=='f')   // convert the xyz binary files from fastpm to a density field. Mianly used for IC
        {
/*          So.message(start_all);
          par_file_bam = argv[2];
          Params params(par_file_bam);
          Catalog tracer(params, "TRACER");
          tracer.read_catalog_bin();*/
          }
      else if(temp=='t')   // Run Patchy
        {

        }

      else if(temp=='m')   // Run Patchy
        {
        }


      else if('u'==temp)  // to make various tasks
        {
              //          tasks(argv[2]);
        }
      else if('?'==temp)
       {
         cout<<endl;
         cout<<"Argument not recognized."<<endl;
         cout<<"Please run ./cosmicatlass.exe -h."<<endl;
         exit(1);
       }

    }
  exit(0);
  return 0;

}
// ##################################################################################
// #################################################################################
// ##################################################################################
// #################################################################################
// ##################################################################################
// #################################################################################
// ##################################################################################
// #################################################################################


