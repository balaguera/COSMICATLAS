
// ##################################################################################
// ################################################################################# 
// COSMICATLASS 
// COSMologIcal
// CATalogs for
// LArge
// Scale
// Structure


/** @file main_cosmicatlassu.cpp
 *
 *  @brief Generation of mock catalosgs of DM tracers
 *  based on the BAM method + ALPT.
 *  @authors Andrés Balaguera-Antolinez, Francisco-Shu Kitaura
 */
# include "../Headers/Bam.h"
# include "../Headers/FileOutput.h"
# include "../Headers/NumericalMethods.h"
# include "Tasks.h"

using namespace std;

// ##################################################################################
// ##################################################################################

int main(int argc, char *argv[]){

  time_t start_all;
  time(&start_all);
  char temp;
  string par_file_bam;
  

  ScreenOutput So(start_all);
  
  if(argc==1){
    So.usage(argv[0]);
    exit(1);
  }
  while((temp =  getopt(argc, argv, "hiam:n:c:t:p:x:r:")) != -1)
    {
      if(temp=='h')
	So.usage(argv[0]);
      
      else if (temp=='a')  // Show authors
        So.author();
      
      else if(temp=='c') // Run cosmicatlas 
	{  

	  So.message(start_all);      
          par_file_bam = argv[2];
	  Params Par(par_file_bam); 
	  Bam bam(Par); 
	  bam.So=So;
          bam.bamrunner();
	  So.message_time(start_all);
	}
      else if(temp=='i') // displays input parameters
	{  
	  par_file_bam = argv[2];
	  Params params(par_file_bam); 
	  Bam bam(params);
	  bam.show_params();
	}
      
      else if(temp=='m')   // to Mesure power spectrum
	{
	  So.message(start_all);
	  par_file_bam = argv[2];
	  Params params(par_file_bam);
	  PowerSpectrumF cPSF(params);
          cPSF.compute_power_spectrum(true,true);
	  So.message_time(start_all);
	}
      else if(temp=='t')   // to read input tracer catalg and analyze it. The operations here should be included in BAM
	{
	  So.message(start_all);
	  par_file_bam = argv[2];
	  Params params(par_file_bam);
          PowerSpectrumF cPSF(params);

          /*
          Catalog tracer(params, "tracer");
	  tracer.read_catalog(params._file_catalogue());
	  // vector<real_prec>dummy(params._NGRID(),0);
	  // tracer.get_number_density_field_grid(bool weight_mass=false, dummy);
	  tracer.get_mass_function();
	  So.message_time(start_all);
*/
	}
      

      else if(temp=='r')   // to read input tracer catalg and analyze it. The operations here should be included in BAM
	{
	  // string par_file = "parameters_cosmolib.ini";
	  // CosmoLib Clib(par_file);
	  // Clib.feed_cosmo_struct();
	  // Clib.get_cosmolib();
	  
	}
      
      
      else if('x'==temp)
          {
          // to make various tasks

          par_file_bam = argv[2];
          Params params(par_file_bam);
          ULONG NTT=params._Nft()*params._Nft()*params._Nft();
          string fileo=params.Output_directory;
          string filen=params.Output_directory;
          for(int j=0;j<6;++j)
          {
            vector<real_prec>dm(NTT,0);
            real_prec expo=1./9.;
            FileOutput File;
            ofstream aou;
            switch(j)
             {
              case(0):
                File.read_array(params.Output_directory+"INVARIANT_TIDAL_I_original.dat", dm);
                fileo=params.Output_directory+"Tinv1.txt";
                filen=params.Output_directory+"Tinv1_new";
               break;
              case(1):
                File.read_array(params.Output_directory+"INVARIANT_TIDAL_II_original.dat", dm);
                fileo=params.Output_directory+"Tinv2.txt";
                filen=params.Output_directory+"Tinv2_new";
              break;
              case(2):
                File.read_array(params.Output_directory+"INVARIANT_TIDAL_III_original.dat", dm);
                fileo=params.Output_directory+"Tinv3.txt";
                filen=params.Output_directory+"Tinv3_new";
             break;
             case(3):
                File.read_array(params.Output_directory+"INVARIANT_SHEAR_I_original.dat", dm);
                fileo=params.Output_directory+"Vinv1.txt";
                filen=params.Output_directory+"Vinv1_new";
              break;
             case(4):
               File.read_array(params.Output_directory+"INVARIANT_SHEAR_II_original.dat", dm);
               fileo=params.Output_directory+"Vinv2.txt";
               filen=params.Output_directory+"Vinv2_new";
             break;
             case(5):
                File.read_array(params.Output_directory+"INVARIANT_SHEAR_III_original.dat", dm);
               fileo=params.Output_directory+"Vinv3.txt";
               filen=params.Output_directory+"Vinv3_new";
             break;
            }
          int Nbins=200;
//          get_overdens(dm,dm);
          vector<real_prec>hist(Nbins,0);
          real_prec delta=2./static_cast<real_prec>(Nbins);
#pragma omp parallel for
          for(ULONG i=0;i<NTT;++i)
            {
              real_prec signx=sign(dm[i]);
              dm[i]=signx*pow(abs(dm[i]), expo);
            }
          real_prec xmin=get_min(dm);
          real_prec xmax=get_max(dm);
          for(ULONG i=0;i<NTT;++i)
            {
              dm[i]=2.*(dm[i]-xmin)/(xmax-xmin)-1.0;
              hist[get_bin(dm[i],-1.0, Nbins, delta,true)]++;
            }
          File.write_array(filen,dm);
          aou.open(fileo.c_str());
          real_prec maxh=get_max(hist);
          for(int i=0;i<Nbins;++i)
             aou<<-1.0+(i+0.5)*delta<<"  "<<hist[i]/static_cast<real_prec>(maxh)<<endl;
          aou.close();

          }

      }

//          tasks(argv[2]);
      
    }
  
  
  return 0;
}
  
  
