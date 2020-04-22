/**
 * @file def.h
 * @brief Pre-processor directives for cosmicatlass
 * @author Andres Balaguera-Antolínez, Francisco-Shu Kitaura
 * @version   1.0
 * @date      2020
 * */


// ********************************************
// ********************************************
// ********************************************
// PREPROCESSOR DIRECTIVES FOR BAM-PATCHY (COSMICATLAS)
// ********************************************
// ********************************************
// ********************************************
// ********************************************
// DEFINE SOME COLORS
/**
 * @brief Reset Color
*/
#define RESET   "\033[0m"
/**
 * @brief Color Black
*/
#define BLACK   "\033[30m"      /* Black */
/**
 * @brief Color Red
*/
#define RED     "\033[31m"      /* Red */
/**
 * @brief Color Green
*/
#define GREEN   "\033[32m"      /* Green */
/**
 * @brief Color Yellow
*/
#define YELLOW  "\033[33m"      /* Yellow */
/**
 * @brief Color Blue
*/
#define BLUE    "\033[34m"      /* Blue */
/**
 * @brief Color Magenta
*/
#define MAGENTA "\033[35m"      /* Magenta */
/**
 * @brief Color Cyan
*/
#define CYAN    "\033[36m"      /* Cyan */
/**
 * @brief Color White
*/
#define WHITE   "\033[37m"      /* White */
/**
 * @brief Color Boldblack
*/
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
/**
 * @brief Color BoldRed
*/
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
/**
 * @brief Color BoldGreen
*/
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
/**
 * @brief Color BoldYellow
*/
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
/**
 * @brief Color BoldBlue
*/
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
/**
 * @brief Color BoldMagenta
*/
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
/**
 * @brief Color BoldCyan
*/
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
/**
 * @brief Color BoldWhite
*/
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ***************************************************PRECISION************************************************************************
// ************************************************************************************************************************************

//
// /**
// * @brief Define double precision for cosmicatlas
// * @details Applies all over the code except for gsl-type defined variables/containers
// */
// #define DOUBLE_PREC


#ifndef DOUBLE_PREC
/**
 * @brief Define single precision for cosmicatlas
 * @details Applies all over the code except for gsl-type defined variables/containers
 * @details Defined when DOUBLE_PREC is undefined
*/
#define SINGLE_PREC
#endif


#define gsl_real double  // Even with single precision, gsl will use double precision for all its calculations

#ifdef SINGLE_PREC
/**
 * @brief Precision at output
*/
#define _PREC_OUTPUT_ 4
/**
 * @brief Precision for FFT operations
*/
#define fftwf_real float


#ifndef real_prec
/**
 * @brief Effective Precision for cosmicatlas
*/
#define real_prec fftwf_real
#endif

/**
 * @brief Effective Precision for complex containers used in FFTW
*/
#define complex_prec fftwf_complex
#endif

#ifdef DOUBLE_PREC
#define _PREC_OUTPUT_ 6
#define fftw_real double
#define real_prec fftw_real
#define complex_prec fftw_complex
#endif


#define BIG_NUMBER 1e7

#define ULONG unsigned long
#define LONG long
#define UULONG unsigned long long
#define ASCII "ascii"

// Precision type of the input par file of the DM field
#define PrecType_X float

// Precision type of the input par file of the Tracer field
#define PrecType_Y float

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************CONSTANTS***************************************************************
// ************************************************************************************************************************************
// Define some constant factors
#define REAL 0
#define IMAG 1
#define ONE static_cast<int>(1)
#define KNOT_MAX_SIZE static_cast<ULONG>(100000000)
#define V_MAX_SIZE 100000000
#define LARGE_NUMBER 1e20
#define MINUS_LARGE_NUMBER -100000000
#define NOCELL -2
#define NOVEL -999

#define I_KNOT 1
#define I_FILAMENT 2
#define I_SHEET 3
#define I_VOID 4

#define NMIN_X_ONECELL 0
#define NMIN_Y_ONECELL 0
#define num_1   static_cast<double>(1.)
#define num_2   static_cast<double>(2.)
#define num_0_1 static_cast<double>(0.1)
#define num_0_5 static_cast<double>(0.5)

// Number of mdoes used to get an estimate of the average bias on large
#define N_MODES 10

// ********************************************
// NUmbers used in the class measuring power spectrum
//#define CHUNK 27
#define ic_rank 3 // For FFTW in 3 dimensions
#define MAX_NUMBER_WEIGHTS 4

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************FUNCTIONS***************************************************************
// ************************************************************************************************************************************

#define myfabs(x) (*((int*)&x + 1) &= 0x7FFFFFFF)

// Bias function used when TEST is defined
#define bias_test(x,alpha, rhoep, ep) (pow(1+x, alpha)*exp(-pow( (1+x)/rhoep, ep)))

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ********************************************************GENERAL STUFF AND OMP*******************************************************
// ************************************************************************************************************************************

// Define this if if's in the loop over the cells for each property is to be avoid
// This allows the binning function to place in the extremes values outside the ranges
#define _BIN_ACCUMULATE_

//#define _SHOW_ISSUES_

// Use OMP parallelization
#define _USE_OMP_

#ifdef _USE_OMP_
#define OMPPARRAN
#define OMPPARRANRSD
#define OMPPARGAR
#endif

// Define VERBOSE to see some values on the screen
//#define _VERBOSE_
// ************************************************************************************************************************************
// ************************************************************************************************************************************
//#define _ADD_POISSON_TO_CIC_
#define MEAN_NEW static_cast<real_prec>(50.0)
// ************************************************************************************************************************************
// ************************************************************************************************************************************




// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ******************************************************MAIN MODE OF  BAM CODE*********************************************
// ************************************************************************************************************************************
// Define MOCK_MODE if calibration or mock production is to be run with cosmicatlas.
// If not defined, then the option BIAS follows, for which the bias statistics is computed


#define MOCK_MODE

// If MOCK_MODE is undefined, then define BIAS
#ifndef MOCK_MODE
#define BIAS_MODE
#endif

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// This is defiend in the case when, in Catalog::read_input_bin, several .dat (binary) files are read as input in positions and velocities,
// and the interpolation of density fields on a grid are perfoemd at while reading these files. This also affects the function, get_density_CIC and NGC in massFunctions
// since we should not initialize the delta arrays there, as they are being filled in parallel with the reading.

//#define _GET_INTERPOLATED_FIELDS_FROM_BIN_FILES_
#ifdef _GET_INTERPOLATED_FIELDS_FROM_BIN_FILES_
//#define _GET_NGP_DENS_FIELD_
#define _GET_CIC_DENS_FIELD_
//#define _GET_TSC_DENS_FIELD_
//#define _GET_VEL_FIELD_
#endif
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// *******************************************************REFERENCE TRACER CATALOG ****************************************************
// ************************************************************************************************************************************

// Use the reference catalog. This is available under _DO_BAM_CALIBRATION_
// and allows to generate the tracer density field from the refenrece.
// The reference catalog is in the Catalog_file variable of the param file 
// If undef, the code shall read TR_DENS_FIELD.txt

//#define _READ_REF_CATALOG_


// If defined: The code selects objects from the input ref catalog using a minimum MASS, even if the observable is other quantity such as the vmax
//#define _SET_CAT_WITH_MASS_CUT_


#ifndef _SET_CAT_WITH_MASS_CUT_
// Define this when power spectrum (-m opiton) is to me measured in cuts of VMAX
// In this case undef  _USE_MASS_TRACERS_
// Recall to undef _SET_GLOBAL_MASS_CUT_ below to use the -m option
// For -c option, undefine, and let mass define the cut
//#define _SET_CAT_WITH_VMAX_CUT_
#endif

#define MASS_SCALE static_cast<double>(1e12)

// This allows to read (if already created from the ref cat) the mass density field of the tracers
// Define also when  _READ_REF_CATALOG_ is defined in order to allow the reading of the mass tracer
#define _USE_MASS_TRACERS_

//#define _USE_SAT_FRACTION_
//#define _USE_VELOCITIES_TRACERS_

// ***********************************************

// Defining something as observables, asks the code
//to compute all possible diagnosis as a function of cuts or bins
// in that particular quantuity
// The user must indicate the position of the VMAX in the -ini file in the slot for i_mass
// and specify mins and max in logMMIN and logMMX

//#define _USE_MASS_AS_OBSERVABLE_
#define _USE_VMAX_AS_OBSERVABLE_

#ifdef _USE_VMAX_AS_OBSERVABLE_
#define _USE_VMAX_TRACERS_
#endif




// ***********************************************
// ***********************************************
// ***********************************************
// ***********************************************
#define test_vmax_mass
//if this is defined, the code *only* explores the vmax-mass relation (i.e, ignoring DM properties such as ð, CWC, Mk)
// and uses it to assign masses, once vmax has been assigned with the usual dm properties.
// In the end we obtain the right vmax-mass relation, the correct vmax function and
// 3% residuals (over all mass range) in the mass function. Doing this allows us to increase the number of vmax bins
// as is explicitely shown below in the definition of N_VMAX_BINS

#ifdef test_vmax_mass
#define _add_dm_density_
#endif

// ***********************************************
// ***********************************************
// ***********************************************
// ***********************************************
// ***********************************************


// ***********************************************
// This asks the code to assign halo masses after having assingned vmax
#ifdef _USE_VMAX_AS_OBSERVABLE_
#define _ASSIGN_MASS_POST_
#ifdef test_vmax_mass
#define N_VMAX_BINS static_cast<int>(500)
#else
#define N_VMAX_BINS static_cast<int>(50)
#endif
#endif
// ***********************************************

#if defined _USE_MASS_TRACERS_ || defined _USE_VMAX_TRACERS_
#define _ASSIGN_PROPERTY_
#endif


#ifdef _USE_MASS_TRACERS_
//#define _USE_MASS_FIELD_
//#define _test_mass_assign_    // UNDEF WHEN MASS ARE DONE!! THIS IS MEANT TO SPEED THE CODE UNDER TESTINT THE COLLAPSE OF RANDOMS
#endif





//#define _GET_DIST_MIN_SEP_REF_

//#define _GET_DIST_MIN_SEP_MOCK_


// ***********************************************
// Use the Log10 of the trcer property under study
//#define _USE_LOG_MASS_
// Use the Xroot of the tracer property under study
#ifdef _USE_MASS_AS_OBSERVABLE_
//#define _USE_X_ROOT_MASS_   // To be deprecated
#endif

#ifdef _USE_X_ROOT_MASS_
// Mass to be transformed as M->M**(exponent_mass_tracer) in order to get better rsults on mass assignment
// At assignation time, masses are reconverted to M->M**(1/exponent_mass_tracer)
// the exponent 1./4. seems to be better than 0.5, as long as high bins in theta are used.
#define exponent_mass_tracer static_cast<real_prec>(1./4.)
#else
#define exponent_mass_tracer static_cast<real_prec>(1.0)
#endif
// ***********************************************


// Maximum of the property (log, or sqrt, or whatever is selected above)
#define MAX_PROPERTY 16
#define _COUNTS_ "COUNTS"
#define _MASS_ "mass"

#define _SAT_FRACTION_ "sat_fraction"

// Please fix this
#ifdef _USE_MASS_AS_OBSERVABLE_
#define _MASS_LOG_   //ueful to measaure mass (or vmax) function
#define MBINS
#elif defined (_USE_VMAX_AS_OBSERVABLE_)
#define MBINS
#define _MASS_LOG_   //ueful to measaure mass (or vmax) function
#endif
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ***********************************************************GENERATION OF MOCK CATALOG **********************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************

#ifdef MOCK_MODE
#define _GET_BAM_CAT_
#endif

// If defined, the code prints the correlation betweeen the halo mass and the properties used to measure the conditional mass functio (or make the list of halo masses in a theta bin)
#define _GET_HALO_MASS_PROPERTIES_CORRELATIONS_
// ************************************************************************************************************************************
// INFORMATION ON PROPERTIES USED TO ASSIGN MASSES:            NEWIGHBOURS
// ************************************************************************************************************************************

//#define _USE_NUMBER_OF_SATELLITES_ //not usedful, ince the mockhas not assigned mass

#define N_CELLS_BACK_FORTH 1

//#define _USE_NEIGHBOURS_

#ifdef _USE_NEIGHBOURS_

//The following three option are to be used disjointly::

#define _USE_NUMBER_OF_NEIGHBOURS_
//Number of close tracers within a range SCALE_MAX_N_NEIGHBOUR

//#define _USE_MIN_DISTANCE_TO_NEIGHBOURS_
//The following three option are to be used disjointly

// Information for local clusering. In the end it si the same as using the number of neighbours
//#define _USE_LOCAL_CLUSTERING_


//#define MIN_LOCAL_CLUSTERING static_cast<real_prec>(0)
//#define MAX_LOCAL_CLUSTERING static_cast<real_prec>(1.2)

// This number is used for _USE_MIN_DISTANCE_TO_NEIGHBOURS_ and _USE_LOCAL_CLUSTERING_
// these numbers have to be tuned, e.g for UNITSIM (N_NEIGH_MAX,SCALE_MAX_N_NEI)=(62,8), (28,4)
// Also we can fix N_NEIGHBOURS_MAX and if a tracer happens to have more neighb, put them in the last bin
// When local clustering is used, SCALE_MAX-- =8Mpc
// N_NEIGHBOURS_MAX acts as the number of bins when the number of is used
// and takes the slot of I_C1

#define N_BINS_MIN_DIST_TO_NEI static_cast<int>(200)

#define N_NEIGHBOURS_MAX static_cast<int>(10)


//#define _USE_LOG_DIST_
#ifdef _USE_LOG_DIST_
#define MAX_OF_MIN_SEPARATION static_cast<real_prec>(log10(20.0))
#define MIN_OF_MIN_SEPARATION static_cast<real_prec>(-2)
#else
#define MAX_OF_MIN_SEPARATION static_cast<real_prec>(15.0)
#define MIN_OF_MIN_SEPARATION static_cast<real_prec>(0)
#endif
#endif

#define SCALE_MAX_N_NEIGHBOURS static_cast<real_prec>(12.0)

// ************************************************************************************************************************************
// *************************************************SEPARATION WITHIN CELLS ***********************************************************

//#define _USE_MIN_SEPARATIONS_IN_CELLS_  //When defined for mocks results are not that good. High Vmax are spoiled.
#define N_BINS_MIN_SEP_IN_CELLS  static_cast<int>(20)  // use 20 for tests, 100 for optimal mass assignment

//#define _USE_LOG_MSIC_
// This information goes into the slot for I_C_BIN3
#ifdef _USE_LOG_MSIC_
#define MIN_SEP_IN_CELLS static_cast<real_prec>(-2.0) // set 0 y linear, -2 if log This is used when mass is assigned to mocks.
#define MAX_SEP_IN_CELLS static_cast<real_prec>(log10(7.0))
#else
#define MIN_SEP_IN_CELLS static_cast<real_prec>(0) // set 0 y linear, -2 if log This is used when mass is assigned to mocks.
#define MAX_SEP_IN_CELLS static_cast<real_prec>(7.0)
#endif


#define DELTA_MIN_SEP  (MAX_SEP_IN_CELLS-MIN_SEP_IN_CELLS)/static_cast<real_prec>(N_BINS_MIN_SEP_IN_CELLS)

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// *************************************************************REDUCED MASS **********************************************************
//#define _GET_DIST_REDUCED_MASS_
#define N_BINS_REDUCED_MASS 200  // use 20 for tests, 100 for optimal mass assignment
#define MAX_REDUCED_MASS 300  // use 20 for tests, 100 for optimal mass assignment
#define MIN_REDUCED_MASS 100  // use 20 for tests, 100 for optimal mass assignment
#define DELTA_REDUCED_MASS    (MAX_REDUCED_MASS-MIN_REDUCED_MASS)/static_cast<real_prec>(N_BINS_REDUCED_MASS)   //200 is roughly th maximum of mu in the sample
#define N_BINS_DIST_EX 10   //NUmber of bis in separation to get dist of reduced masses

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// *****************************************************TRACERS IN CELLS **************************************************************
//#define _USE_TRACERS_IN_CELLS_
// Use the number of tracers in cells as a proxy for mass assigment.
#define N_TRACERS_IN_CELLS_MAX  15   // This is the maximumn number of tracers in a cell + 1
// ************************************************************************************************************************************
// *****************************************************TOTAL MASS IN CELLS **************************************************************
#define _USE_TOTAL_MASS_IN_CELL_ // this causes problems with the mock.
#define N_BINS_TOTAL_MASS_IN_CELL  static_cast<int>(20)
#define MAX_TOTAL_MASS_IN_CELL static_cast<real_prec>(1e15)

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ********************************************************MASS ASSIGNMENT*************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************

// ************************************************************************************************************************************
// ********************************************************TESTING *********************************************************************
// Define this in case we want to reassign masses to the reference catalog, in
// order to see whether the bias, from mass cuts, is still that from the reference
// or if the mass assignment induces the wrong clustering
// This option does not allow to use any propbabilistic

#define _MASS_ASSIGNMENT_TO_REFERENCE_
// ************************************************************************************************************************************
// ***************************************************COLLAPSE RANDOMS ***************************************************************************

#ifndef _MASS_ASSIGNMENT_TO_REFERENCE_
// Move random  particles towads their closest dark matter particle when assigning coordinates to Bam outputs with patchy
//#define _COLLAPSE_RANDOMS_    // This is implemented *after assigning* masses. Mass depencendies can be applied in this case
// However, if some properties below the scale of the cell are to be measured, we reccoment to collapse before mass assignment
//#define _COLLAPSE_RANDOMS_AUX_ // This was used in order to move the randoms at the same time when positions of DM are assigned
#endif

// ************************************************************************************************************************************
// ************************************************************************************************************************************

#define _USE_MASS_ASSIGNMENT_READING_REF_MASSES_
// This allows the code to assign masses ussing the actual values from the reference tracer catalog.
// abd uses the function  Bam::assign_tracer_mass_new_new()
// If not defiend, the code uses the assignement based on mass bins and uses  Bam::assign_tracer_mass_new()
// Note: the fnction  Bam::assign_tracer_mass_new() is marginally sensitive to the bins of mass in which the mass function os measured
// The function  Bam::assign_tracer_mass_new_new() is still less sensitive for the latter keeps the list of masses in a bin of DM properties.
// Note: this is planned to be the stabndard assignment. Whehn other propertioes are to be assigned, we shall use a probabilistic approach
// ************************************************************************************************************************************
// ***************************************************MULTISCALE MASS ASSIGNMENT  *****************************************************
// ************************************************************************************************************************************
// If this order is defined, the masses above a threshold Ms are assigned using a loop over the grids, chosen aribtrarely,
// ion order to prevent that the higherst value of mass are assigned in the same cell and avoid exlcusion.
// If not defiedn, the residuals at MCUT is 48%. When defined, it reduces to 9%
#define  _USE_MULTISCALE_MASS_ASSIGNMENT_
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Mass threshold: for mass es above thes values, masses are assigned in a grid-exclusion based approach. THis aims at accounting for exlcusion
// FOr masses below this values, masses are taken form the allowed lists, in a tracer-based approach. This neglects exclusion

// Mass threshold above which a low-resolution cell is used in order to assign masses
// If this number is larger than the maximum tracer mass, PLEASE undef _USE_MULTISCALE_MASS_ASSIGNMENT_

#define MASS_THRESHOLD static_cast<real_prec>(1e11)  //this is now mass threshold_leve1
// -----------------------------------------------------------------------------------------

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifdef _USE_MULTISCALE_MASS_ASSIGNMENT_

// USE_EXCLUSION_CELL DOES NOT OWRK. BEST RESULTS ARE OBTAINED WITHOUT IT!!!!!!!!
// -----------------------------------------------------------------------------------------
#define _USE_MULTISCALE_LEVEL_1_ // Original mesh resolution. ALWAYS defined

// LEVEL 1 HAS ALSO THE OPTION TO APPLY EXCLUSION BY MASS-BIN
// demanding that in every cell there can be only a number (MAX_NUMBER_OF_MASSES_IN_CELL_SAME_MBIN)
// of objects. The code gets flexible when it still needs to
// assign N_MASSES_ABOVE_M_EX masses in this limit.
// Exclusion can be made stroner either decreasing MAX_NUMBER_OF_MASSES_IN_CELL_SAME_MBIN
// or decreasing N_MASS_BINS_EXC
// The parameter N_MASSES_ABOVE_M_EX hs to be increased if strong exlcusion is in. DOES NOT WORK

#ifdef _USE_MULTISCALE_LEVEL_1_

//#define _USE_EXCLUSION_WITH_MASS_BINS_

#ifdef _USE_EXCLUSION_WITH_MASS_BINS_
#define N_MASS_BINS_EXC 80
#define DELTA_MASS_EXC (MASS_MAX)
#define MAX_NUMBER_OF_MASSES_IN_CELL_SAME_MBIN  2
#define N_MASSES_ABOVE_M_EX 200
#endif

//#define _USE_EXCLUSION_CELL_LEVEL1_
#ifdef _USE_EXCLUSION_CELL_LEVEL1_
#define NUMBER_OF_NEIGHBOUR_CELL_FOR_EXCLUSION_LEVEL1 static_cast<int>(1)
#define N_MASSES_ABOVE_THRESHOLD 50
#endif
// If this mass is below the min of the catalog, all mass will be assigned with
// from a multiscale approach
#ifdef _USE_MASS_AS_OBSERVABLE_
#define MASS_THRESHOLD_MULTI_SCALE_1 static_cast<real_prec>(2e13)
#elif defined (_USE_VMAX_AS_OBSERVABLE_)
#define MASS_THRESHOLD_MULTI_SCALE_1 static_cast<real_prec>(400.0)  // 400 km/s seems to be an optimal value up to which particle-based assignmet can proceed
#endif
#endif // end use_multiscale_level_1
// -----------------------------------------------------------------------------------------


//#define _USE_MULTISCALE_LEVEL_2_   // Lowest resolution NOT YET INCLUDED IN Bam.cpp
// LEVEL 2 HAS THE OPTION TO INCREASE THE GRID SIZE IN ORDER TO APPLY
// EXLCUSION ON SCALES LOWER THAN THE FIDUCIAL RESOLUTION SCALE. DOES NOT HELP
#ifdef _USE_MULTISCALE_LEVEL_2_
//#define _USE_EXCLUSION_CELL_LEVEL2_
#ifdef _USE_EXCLUSION_CELL_LEVEL2_
#define NUMBER_OF_NEIGHBOUR_CELL_FOR_EXCLUSION_LEVEL2 static_cast<int>(1)
#endif
#ifdef _USE_MASS_AS_OBSERVABLE_
#define MASS_THRESHOLD_MULTI_SCALE_2 static_cast<real_prec>(1e13)
#elif defined _USE_VMAX_AS_OBSERVABLE_
#define MASS_THRESHOLD_MULTI_SCALE_2 static_cast<real_prec>(640)
#endif

#define _HIGHEST_RES_LEVEL2_
#ifdef _HIGHEST_RES_LEVEL2_
#define  NFT_LOW_2 static_cast<int>(512)// this mut be even factors of the fiducial value
#else
#define  NFT_LOW_2 static_cast<int>(160)
#endif
#define  NGRID_MS_2 static_cast<ULONG>(NFT_LOW_2*NFT_LOW_2*NFT_LOW_2)
// NGRID_MS_2 corresponds to the original resolution NGRID
#endif //end use_multiscale_level_2
// -----------------------------------------------------------------------------------------

//#define _USE_MULTISCALE_LEVEL_3_
#ifdef _USE_MULTISCALE_LEVEL_3_
//#define _USE_EXCLUSION_CELL_LEVEL3_
#ifdef _USE_EXCLUSION_CELL_LEVEL3_
#define NUMBER_OF_NEIGHBOUR_CELL_FOR_EXCLUSION_LEVEL3 static_cast<int>(1)
#endif
#ifdef _USE_MASS_AS_OBSERVABLE_
#define MASS_THRESHOLD_MULTI_SCALE_3 static_cast<real_prec>(1e13)
#elif defined _USE_VMAX_AS_OBSERVABLE_
#define MASS_THRESHOLD_MULTI_SCALE_3 static_cast<real_prec>(700)
#endif
#define  NFT_LOW_3 static_cast<int>(64)
#define  NGRID_MS_3 static_cast<ULONG>(NFT_LOW_3*NFT_LOW_3*NFT_LOW_3)
#endif  //end use_multiscale_level_3
// -----------------------------------------------------------------------------------------
//#define _USE_MULTISCALE_LEVEL_4_
#ifdef _USE_MULTISCALE_LEVEL_4_
//#define _USE_EXCLUSION_CELL_LEVEL4_  // do not use. When defined, large scale pwoer is spolied
#ifdef _USE_EXCLUSION_CELL_LEVEL4_
#define NUMBER_OF_NEIGHBOUR_CELL_FOR_EXCLUSION_LEVEL4 static_cast<int>(1)
#define N_MASSES_ABOVE_THRESHOLD_LEVEL4 static_cast<int>(30)
#endif
#ifdef _USE_MASS_AS_OBSERVABLE_
#define MASS_THRESHOLD_MULTI_SCALE_4 static_cast<real_prec>(6e13)
#elif defined _USE_VMAX_AS_OBSERVABLE_
#define MASS_THRESHOLD_MULTI_SCALE_4 static_cast<real_prec>(640.0)
#endif
#define NFT_LOW_4 static_cast<int>(256)
#define NGRID_MS_4 static_cast<ULONG>(NFT_LOW_4*NFT_LOW_4*NFT_LOW_4)
#endif //end use_multiscale_level_4

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif  // end ifdef _USE_MULTISCALE_MASS_ASSIGNMENT_
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// This avoids the convolution of the DM field with the kernel, in order to see whether the convolution makes an effeect on the preciscion
//#define _DO_NOT_CONVOLVE_
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// **************************************************BOUNDARIES FOR HISTOGRAMS *******************************************************
//This is important when assignning masses to tracers, in order to use at its most all inforamtion
#define _MODIFY_LIMITS_


// Define this if the glocbal mass funciton is to be used in the assignment of mass
// For a coprrect mass assignment, leave it undef
//#define _USE_GLOBAL_MASS_FUNCTION_
// ************************************************************************************************************************************



// ************************************************************************************************************************************
// ************************************************************************************************************************************
// **************************************************EXCLUSION ************************************************************************
//#define _CORRECT_FOR_EXCLUSION_   // Toooooooo slowwwwwwwww
// This is used in order to alloca the list of masses in pairs between the mini sep of the reference and MAXIMUM_DISTANCE_EXCLUSION
#define MAXIMUM_DISTANCE_EXCLUSION static_cast<real_prec>(7.0) // 1.5 Mpc/h as maximum separation explored??
#define MINIMUM_DISTANCE_EXCLUSION static_cast<real_prec>(0.1) // 1.5 Mpc/h as maximum separation explored??
#define EXCLUSION_SCALE static_cast<real_prec>(2.5) // Correct masses of pairs within this separation

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************

#define _VEL_UNITS_MPC_PER_h_    // define if the velocities in the final cat are to be given in Mpc/h

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ***********************************************************POWER SPECTRUM **********************************************************
// ************************************************************************************************************************************
// Directives for the section of the code measuring Pk,(called with -m at compillation time)
// This definition has to work together with the one in the input parameter file.
// This only works on the method PowerSpectrumF::compute_power_spectrum(book, bool)
// to get power spectrum reading a catalog

//#define _CHECK_PARTICLES_IN_BOX_
//#define _USE_WEIGHTS_IN_POWER_              //This must be defined if weighs, as specified in the parameter file, are to be used in the determination of power spectrum
//#define _GET_BISPECTRUM_NUMBERS_   // This doe snot compute bispctrum, but allows to compute numbers used in the shot noise and normalziation of that statistics. Verify that this is de fined when that option is demanded from parfile

//#define _MASS_WEIGHT_POWER_


// These three disjoint options are to be used with the -m compilling flag.

#define _USE_ALL_PK_
//#define _USE_MASS_CUTS_PK_
//#define _USE_MASS_BINS_PK_


// DEFINE when the option -m is to be used
//#define _POWER_
#ifdef _POWER_
#ifdef _USE_X_ROOT_MASS_
#undef _USE_X_ROOT_MASS_
#endif
//#define _USE_MASS_CUTS_IN_MARKED_PK_
#endif



// This is meant to measure power spectrum when passing the structure tracer.Halo
// This is due to the fact that when one reads the cat and measures P(k), one can pass the minimum cut
#ifdef _USE_MASS_CUTS_PK_
#define _SET_GLOBAL_MASS_CUT_
#ifdef _USE_MASS_AS_OBSERVABLE_
#define MINIMUM_PROP_CUT_original static_cast<real_prec>(2e13) // put the largest of the mcuts used in power, for comparison purposes
#define MINIMUM_PROP_CUT static_cast<real_prec>(pow(MINIMUM_PROP_CUT_original, exponent_mass_tracer))
#elif defined _USE_VMAX_AS_OBSERVABLE_
#define MINIMUM_PROP_CUT static_cast<real_prec>(500.0)
#endif
#endif


//#define _WRITE_2DPOWER_
//#define _WRITE_MULTIPOLES_

// Define this to use the vectorized version of grid assignment, written by L. Tornatore.
// That version has some memmory issues depending on the computer the code is run at.
// If undef, the code uses my own version which is  efficient enough
// IMPORTANT: when using the standar assignment, I demand here to specify whether PSC
// is to be used, in order to speed the code. I will set a warning in the parameter  file
//#define _USE_VECTORIZED_GRID_ASSIGNMENT_
#undef _USE_VECTORIZED_GRID_ASSIGNMENT_
#ifndef _USE_VECTORIZED_GRID_ASSIGNMENT_

//#define _USE_FOURTH_ORDER_
#ifdef _USE_FOURTH_ORDER_
#define MAX_MAS_DEG 5
#else
#define MAX_MAS_DEG 3
#endif
#endif
#define CHUNK MAX_MAS_DEG*MAX_MAS_DEG*MAX_MAS_DEG   // this is 5*5*5



#define KMAX_RESIDUALS_high static_cast<real_prec>(0.5)
#define KMAX_RESIDUALS_low static_cast<real_prec>(0.3)

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// *************************************************PATCHY ****************************************************************************
// ********************************************
// Use patchy. Define if we want to produce DM fields with Patchy.
// If undef, the input DM file must be provided in the parameter file.

//#define _USE_PATCHY_

#ifdef _USE_PATCHY_
// Define if ony a simuatoipn with patchy is to be run. The code stopes before BAM starts doing the
// calibration or the mock production.
//#define _ONLY_PATCHY_

#define _GET_VELOCITIES_PATCHY_  //get the interpolated velocity field for patchy


#endif

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// *************************************************BAM GENERAL ***********************************************************************
// ********************************************


// avoid writing or computing files already existing
// Not all files are under this def, COMPLETE!
//#define _VERIFY_FILES_

// ********************************************
// Change coordiantes x<->y in the input DM  density field
//#define _EXCHANGE_X_Y_DENSITY_
// ********************************************
// Change coordiantes x<->y in the input DM  density field
//#define _EXCHANGE_X_Y_DENSITY_
// ********************************************
// Change coordiantes x<->y in the input DM velocity fields
//#define _EXCHANGE_X_Y_VEL_

// ********************************************
// Compute and write to putput files PDF
//#define WRITE_PDF

// ********************************************
// Define if the TR density field is to be written (bin) during the iteration
//#define _WRITE_ALL_FIELDS_
#ifndef _WRITE_ALL_FIELDS_
#define _WRITE_TR_DENSITY_FIELD_
//#undef _WRITE_TR_DENSITY_FIELD_
#endif
// ********************************************
// Define if the DM density field is to be written (bin)
#define _WRITE_DM_DENSITY_FIELD_

// ********************************************
// Define this if DM spectra duwing iterations is to be produced
#define _GET_POWER_REFS_

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************

// ********************************************
// Define if several realiations are used to compute bias statistics from. Narems in this case are hard-coded
//#undef _SEVERAL_REAL_
// ********************************************
#ifdef BIAS_MODE
// This isused in the context of chaning the reference DM when bias are below one
// situation in which BAM seems to be out of order. 

#define NUM_IN_LOG static_cast<real_prec>(1.0000001)
//#define NUM_IN_LOG static_cast<real_prec>(2)
// Define this if the Joint B is to be computed from several realizations of the DM field (alpt so far) and the reference
//#define _SEVERAL_REAL_BIAS_
#endif

// ********************************************
// Compute smoothed version of the BAM kernel.
// Verified that smoothing is best. 
//#define _SMOOTHED_KERNEL_
// ********************************************
// Define if smoothed is to be done only in the last iteration
//#define _SMOOTHED_KERNEL_LAST_ITERATION_

// ********************************************
// ********************************************
// This is the number ised in log(NuUM_IN_LOG+delta). If
// Undesirable situations with log10 are to be avoidded, I suggest this number to be set to 2.
// specially when mocks are to be built. If plots are to be done after using the -n option, set 1
#ifdef MOCK_MODE
#define NUM_IN_LOG static_cast<real_prec>(2)
// --------------------------------------------
// Define to run the test in which a Halo density field is created from the input DM density
// field using a particular halo bias, given by bias_test
//#define _RUN_TEST_
#endif


//#define _UNDER_BIASED_

// ********************************************
// Dynamical way of populating bins of density, MK, etc with the reference number of cells
#define _DYNAMICAL_SAMPLING_  //optimal, faster


// ********************************************************************************************************************************************
#ifdef MOCK_MODE

// ********************************************
// Generate realizations by applying the BAM kernel to independent realizations of DM density fields. 
// If this is not defined, then one (or seveal from different IC according to _SEVERAL_REAL) calibrations are done
//#define _GET_BAM_REALIZATIONS_

// ********************************************
// If no mocks is to be created, the code assumes that we watn to perform the calibration of kernel.

#ifndef _GET_BAM_REALIZATIONS_
#define _DO_BAM_CALIBRATION_
#endif


// ********************************************
#ifdef _GET_BAM_REALIZATIONS_

#define _ONLY_POST_PROC_     //This is meant to be used when the positions in catalog have been *already* assigned and post-processing, i.e, assign masses, follows
//#define _EXTRAPOLATE_VOLUME_
#endif


// THE OPION use_tracer_hr is allowed to be on only if we are doing the calibration
#ifndef _EXTRAPOLATE_VOLUME_
//This option is not working, The kernel, when checked with UNITsim, goes above 1 and
// no convergece is acquired.
//#define _USE_TRACER_HR_
#endif
// ********************************************

#ifdef _DO_BAM_CALIBRATION_
// This is useful when calibration over different ref catalogs are to be done
// Note taht the paths to input files are hard coded in Bam.cpp
//#define _SEVERAL_REAL_CAL_
#endif


#endif   // end of get_bam_realizations 791


// ********************************************
// Define if the limits in the DM properties used to measure the bias are to be updated after each iteration
// This option is not working, For some reason it makes the method
// to get stuck.
// Update: bug related with DM: when turned off CWC and MK, got WN in the 0th iteration.
// Solved: new_nbins_x was not assigned, then when computing the new value of delta, it was undef.
// It worked before becuase DELTAX was computed before.
#ifdef _DO_BAM_CALIBRATION_
//#define _MODIFY_LIMITS_
#endif

// ********************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ******************************************** COSMIC WEB ****************************************************************************
// ********************************************

// Define this if the CWC is to be used as variable for the bias.
// This does not define though whether CWC is computed or not!!!
// If undef, BAM could still do the classification in order to use the eigenvalues
// of the Tidal field as bias variables. This happens if USE_INVARIANTS is defnied
// RECALL TO SET UP THE VARIABLE v_CWT_used in the parameter file
// in case this variable is defined

//#define _USE_CWC_
// If defined, the CWC is done in every iteration of the loop
// If not, the CWC is done only in the first iteration, and that particular classification
// is used in the iterative process.

#define _USE_CWC_INSIDE_LOOP_

// ***********************************************************************************************************************************
// ***********************************************************************************************************************************

// This is meant in case we have CIC DM fields and we want to transform it to NGP
//#define _CONVERT_CIC_TO_NGP_

// ***********************************************************************************************************************************

// ********************************************INVARIANTS OF TIDAL FIELD
// ***********************************************************************************************************************************
// This is only meant to write the invariant in a binary file.
// The invariant I is the same delta, so it is already there
//#define _USE_INVARIANT_TIDAL_FIELD_I_
#ifdef _USE_INVARIANT_TIDAL_FIELD_I_
#define _USE_EXPONENT_INVARIANT_I_
#define _USE_SIGN_INVARIANT_I_
#define _MAP_TO_INTERVAL_INV_I_
#define NEWMAX_INV_I static_cast<real_prec>(1.0)
#define NEWMIN_INV_I static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_INVARIANT_I_
#define EXPONENT_INVARIANT_I  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_INVARIANT_I  static_cast<real_prec>(1.0)
#endif
#endif
// ************************************************************************************************************************************
#define _USE_INVARIANT_TIDAL_FIELD_II_
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
#define _USE_EXPONENT_INVARIANT_II_
#define _USE_SIGN_INVARIANT_II_
//#define _MAP_TO_INTERVAL_INV_II_
#define NEWMAX_INV_II static_cast<real_prec>(1.0)
#define NEWMIN_INV_II static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_INVARIANT_II_
#define EXPONENT_INVARIANT_II  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_INVARIANT_II  static_cast<real_prec>(1.0)
#endif
#endif
// OR
//#define _USE_DELTA2_   // ð².  This uses the memory space of _USE_INVARIANT_TIDAL_FIELD_I_
//#define _WRITE_DELTA2_
// ************************************************************************************************************************************
#define _USE_INVARIANT_TIDAL_FIELD_III_
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
#define _USE_EXPONENT_INVARIANT_III_
//#define _USE_SIGN_INVARIANT_III_
//#define _MAP_TO_INTERVAL_INV_III_
#define NEWMAX_INV_III static_cast<real_prec>(1.0)
#define NEWMIN_INV_III static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_INVARIANT_III_
#define EXPONENT_INVARIANT_III  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_INVARIANT_III  static_cast<real_prec>(1.0)
#endif
#endif

// OR
//#define _USE_DELTA3_   //  ð³  This uses the memory space of _USE_INVARIANT_TIDAL_FIELD_II_
//#define _WRITE_DELTA3_
// ************************************************************************************************************************************


// ********************************************
//#define _USE_TIDAL_ANISOTROPY_
// OR
//#define _USE_S2_   // s²  This takes the memory space from _USE_TIDAL_ANISOTROPY_
//#define _WRITE_S2_

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ******************************************** COLLAPSING REGIONS ********************************************************************
// ********************************************
// Define if MKnots is to be used. In this case, the CWC classification is done, even if _USE_CWC_ is not defined 
//#define _USE_MASS_KNOTS_   // Mk.

// Min and max for the log of the mass of SuperKnots found by the FoF

#define MKMAX 1e4  // I used 1e4 for UNITSIM and Minerva. After applciations with Minerva, set it to 1e5.
#define MKMIN 0.0

//#define MKMAX 1e-27
//#define MKMIN 1e-30

#ifdef _USE_MASS_KNOTS_
//#define _DM_NEW_UNITS_   //define this if youn have dm with cgs units. If a DM is done with counts, undef
#define _WRITE_MKNOTS_
#endif

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ******************************************** VELOCITIES ****************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ********************************************
// Use information from the DM peculiar velocity field.
// in order to characterize the bias
// In particular, use 2 eigenvalues of the shear of the V field from DM

//#define _USE_VELOCITIES_

#ifdef _USE_VELOCITIES_
#ifndef _USE_PATCHY_
#define _READ_VELOCITIES_  // if we use vels, we can read them from files
#endif
#endif

// ************************************************************************************************************************************
// ******************************************** V--WEB ********************************************************************************
// ************************************************************************************************************************************

#ifdef _USE_VELOCITIES_
//#define _USE_CWC_V_  // Use the V_ classification
//#define _USE_VEL_KNOTS_V_  //Use the velocity dispersion of the cells classified as knots
#endif
#define VKMIN  0.0
#define VKMAX  1e10


// ************************************************************************************************************************************
//#ifndef _GET_BAM_REALIZATIONS_
//#define _TEST_THRESHOLDS_RESIDUALS_   // Define only when the thresholds of the T and V eigenvalues are to be explored sistematicaly in thw first (raw) interation
//#endif

// ************************************************************************************************************************************
// ******************************************** V-INVARIANTS **************************************************************************
// ************************************************************************************************************************************
#ifdef _USE_VELOCITIES_

// ************************************************************************************************************************************
//#define _USE_INVARIANT_SHEAR_VFIELD_I_  // def or undef here
#ifdef _USE_EXPONENT_INVARIANT_VS_I_
//#define _USE_EXPONENT_INVARIANT_VS_I_
//#define _USE_SIGN_INVARIANT_VS_I_
#define _MAP_TO_INTERVAL_INV_SHEAR_I_
#define NEWMAX_INV_SHEAR_I static_cast<real_prec>(1.0)
#define NEWMIN_INV_SHEAR_I static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_INVARIANT_VS_I_
#define EXPONENT_INVARIANT_VS_I  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_INVARIANT_VS_I  static_cast<real_prec>(1.0)
#endif
#endif
// ************************************************************************************************************************************
//#define _USE_INVARIANT_SHEAR_VFIELD_II_  // def or undef here
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
#define _USE_EXPONENT_INVARIANT_VS_II_
//#define _USE_SIGN_INVARIANT_VS_II_
#define _MAP_TO_INTERVAL_INV_SHEAR_II_
#define NEWMAX_INV_SHEAR_II static_cast<real_prec>(1.0)
#define NEWMIN_INV_SHEAR_II static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_INVARIANT_VS_II_
#define EXPONENT_INVARIANT_VS_II  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_INVARIANT_VS_II  static_cast<real_prec>(1.0)
#endif
#endif
// ************************************************************************************************************************************
//#define _USE_INVARIANT_SHEAR_VFIELD_III_  // def or undef here
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
#define _USE_EXPONENT_INVARIANT_VS_III_
//#define _USE_SIGN_INVARIANT_VS_III_
#define _MAP_TO_INTERVAL_INV_SHEAR_III_
#define NEWMAX_INV_SHEAR_III static_cast<real_prec>(1.0)
#define NEWMIN_INV_SHEAR_III static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_INVARIANT_VS_III_
#define EXPONENT_INVARIANT_VS_III  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_INVARIANT_VS_III  static_cast<real_prec>(1.0)
#endif
#endif

// ************************************************************************************************************************************
#endif

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ******************************************** BIAS TERMS ****************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************

//#define _USE_NABLA2DELTA_  //  Nabla ²ð   This takes the memory space from _USE_INVARANT_SHEAR_VFIELD_I
//#define _WRITE_NABLA2DELTA_


// #define _USE_S2DELTA_  //  s²ð   This takes the memory space from _USE_INVARANT_SHEAR_VFIELD_II
//#define _WRITE_S2DELTA_


//#define _USE_S3_        //   s³    This takes the memory space from _USE_INVARANT_SHEAR_VFIELD_III
//#define _WRITE_S3_
// ************************************************************************************************************************************
// ************************************************************************************************************************************


// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ********************************************
// When BAM finishes, generate cats with particle posistions according to PATCHY
// ********************************************
// If catalog is to be produced, def if we want velocities or not
#ifdef _GET_BAM_CAT_
#define _GET_VELOCITIES_

#define CATWITHV
#endif

// ********************************************
#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_DELTA2_)
#ifdef _USE_NEIGHBOURS_
#undef _USE_NEIGHBOURS_
#endif
#define N_C_BIN1 static_cast<int>(100)
#define C1_MIN  static_cast<real_prec>(-1e8)
#define C1_MAX  static_cast<real_prec>(1e9)
#define DELTA_C1 (C1_MAX-C1_MIN)/(static_cast<double>(N_C_BIN1))
#else
#define N_C_BIN1 static_cast<int>(1)
#define C1_MIN -200
#define C1_MAX  200
#define DELTA_C1 static_cast<real_prec>(1.)
#endif

// ********************************************
#if defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_DELTA3_)
#define N_C_BIN2 static_cast<int>(100)
#define C2_MIN  static_cast<real_prec>(-1e8)
#define C2_MAX  static_cast<real_prec>(1e9)
#define DELTA_C2 (C2_MAX-C2_MIN)/(static_cast<double>(N_C_BIN2))
#else
#define N_C_BIN2 static_cast<int>(1)
#define C2_MIN -200
#define C2_MAX  200
#define DELTA_C2 static_cast<real_prec>(1.)

#endif
// ********************************************
// This slot is dedicaqtd for the tital anisotriu OR the S2 term
#if defined (_USE_TIDAL_ANISOTROPY_)  || defined (_USE_S2_)
#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
#undef _USE_MIN_SEPARATIONS_IN_CELLS_
#endif
#define N_C_BIN3 static_cast<int>(10)
#define C3_MIN  static_cast<real_prec>(-9e10)
#define C3_MAX  static_cast<real_prec>(9e10)
#define DELTA_C3 (C3_MAX-C3_MIN)/(static_cast<double>(N_C_BIN3))
#else
#define N_C_BIN3 static_cast<int>(1)
#define C3_MIN -200
#define C3_MAX  200
#define DELTA_C3 static_cast<real_prec>(1.)
#endif

// ********************************************
// This slot is dedicaqtd for the invariant odf the vshear field OR the Nabla² delta term
#if defined (_USE_INVARIANT_SHEAR_VFIELD_I_) || defined (_USE_NABLA2DELTA_)
#define N_CV_BIN1 static_cast<int>(50)
#define CV1_MIN  static_cast<real_prec>(-40)
#define CV1_MAX  static_cast<real_prec>(40)
#define DELTA_CV1 (CV1_MAX-CV1_MIN)/(static_cast<double>(N_CV_BIN1))
#else
#define N_CV_BIN1 static_cast<int>(1)
#define CV1_MIN  static_cast<real_prec>(-40)
#define CV1_MAX  static_cast<real_prec>(40)
#define DELTA_CV1 static_cast<real_prec>(1.)

#endif
// ********************************************
#if defined (_USE_INVARIANT_SHEAR_VFIELD_II_) || defined (_USE_S2DELTA_)
#define N_CV_BIN2 static_cast<int>(20)
#define CV2_MIN  static_cast<real_prec>(-40)
#define CV2_MAX  static_cast<real_prec>(40)
#define DELTA_CV2 (CV2_MAX-CV2_MIN)/(static_cast<double>(N_CV_BIN2))
#else
#define N_CV_BIN2 static_cast<int>(1)
#define CV2_MIN  static_cast<real_prec>(-40)
#define CV2_MAX  static_cast<real_prec>(40)
#define DELTA_CV2 static_cast<real_prec>(1.)
#endif

// ********************************************
#if defined (_USE_INVARIANT_SHEAR_VFIELD_III_) || defined (_USE_S3_)
#define N_CV_BIN3 static_cast<int>(50)
#define CV3_MIN  static_cast<real_prec>(-40)
#define CV3_MAX  static_cast<real_prec>(40)
#define DELTA_CV3 (CV3_MAX-CV3_MIN)/(static_cast<double>(N_CV_BIN3))
#else
#define N_CV_BIN3 static_cast<int>(1)
#define CV3_MIN  static_cast<real_prec>(-40)
#define CV3_MAX  static_cast<real_prec>(40)
#define DELTA_CV3 static_cast<real_prec>(1.)
#endif





// Transform the vels from km/smto Mpc/h  / sec
// If input velñs are already in Mpc/h/sec, undef
#define _VEL_KMS_TO_MPCHS_




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//IMPORTANT DEFINITIONS FOR PATCHY
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//ONLY GALAXY NUMBER COUNTS CALCULATION
#define GONLY
//ONLY DM FIELD CALCULATION

#define DMONLY
//#undef DMONLY

#define TET
//#undef TET

///////////////////////////////////////////////7
//STRUCTURE FORMATION MODEL
////////////////////////////////////////////////

#undef MUSCLE
#undef MUSCLE2LPT
#undef PTCURL
#undef SMOOTHCURL

#undef LPT3A
#undef LPT3B


///////////////////////////////////////////////
//RSDs
////////////////////////////////////////////////

#define FOGS
#undef NORSD



///////////////////////////////////////////////
//BASIC DEFINITIONS
////////////////////////////////////////////////

#define SAMPRAN

// ********************************************
#define SAVEMEM

// ********************************************
#undef PHOTOZ

// ********************************************
#define PARTRAN // distribute halos/galaxies inside cells using the position of the available DM particles

// ********************************************
#undef UNIRAN   // distribute halos/galaxies inside cells using uniform random distribution

// ********************************************
#undef CICRAN   // distribute halos/galaxies inside cells using pseudo-cic

// ********************************************
#define SAVECOMP // safe calculations for RSD at the expense of writing more fields

// ********************************************
#define WCAT

// ********************************************
#undef COMPACTTWOLPT

// ********************************************
#undef CELLBOUND

// ********************************************
//#define SAVEIC
#undef SAVEIC

// ********************************************
#undef GFINDIFF// gradient with finite differences
#define  GFFT    // gradient with FFTs
// ********************************************

#define  _USE_GFINDIFF_EIGENV_    // gradient with FFTs FOR eIGENVALUES
//#define  _USE_GFFT_EIGENV_    // gradient with FFTs  FOR EIGENVALUES

// ********************************************
//#define _USE_ZERO_PADDING_POT_
#define _EXTRA_NFT_FACTOR_   static_cast<int>(2)
// ********************************************
#undef TRANSF
// ********************************************
#undef TRANSFSC
// ********************************************
#undef TRANSFDENS
// ********************************************
#define LOGRO
// ********************************************
#undef SCSMOO

///////////////////////////////////////////////7
///////////////////////////////////////////////7
///////////////////////////////////////////////7
//BASIC DEFINITIONS
////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FOURIER
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//#define FOURIER_DEF_1  // only tested with DEF_2
#define FOURIER_DEF_2  // only tested with DEF_2
#define FFTW_OPTION FFTW_ESTIMATE
#define FORWARD FFTW_FORWARD 
#define BACKWARD FFTW_BACKWARD
//#define fwd true
//#define inv false
#define fourier_space true
#define real_space false
#define to_Fspace true
#define to_Rspace false
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define _CONVERT_MPC_per_H_TO_KM_per_SECOND_

#define cgs_Mpc num_1
#define cgs_sec num_1
#define cgs_km static_cast<real_prec>(0.3240779290e-19 * cgs_Mpc) // 1 kilometer in 1 Mpc
#define cgs_clight static_cast<real_prec>(0.9715611892e-14 * cgs_Mpc/cgs_sec)

#define eps static_cast<real_prec>(1.e-14)
#define epsint static_cast<real_prec>(1.0e-6) /* numerical accuracy for integrations */
#define NEVAL 1000    /* numerical integration a2com */

#define _MESS_ So.message_screen("ABA", __LINE__);
#define _ACHTUNG_ So.message_warning("ACHTUNG in line ", __LINE__);



// ========================================================================================

// Define if Y is in NGP and want to be converted to CIC in order to smooth contours.
//By default
//#define  _KONV_


#ifdef BIAS_MODE
// Use this to convert NGP to CIC when using number counts in the tracers and want to display the bias
//#define  _NGP2CIC_Y_

//Uset this to convolve the DM with a kernel, which is read from the output directory
//#define  _KONV_




//define if a mask containig redshifts is to be used for bias analysis involving reconstrunction
//#define _USE_REDSHIFT_MASK_

//#define _USE_BINARY_MASK_  // if the file binary_mask is 1 or 0, define this. If not, do not use it



#ifdef _USE_REDSHIFT_MASK_
//#define REDSHIFT_MIN static_cast<real_prec>(0.35)
//#define REDSHIFT_MAX static_cast<real_prec>(2.34)
//#define N_REDSHIFT_BINS static_cast<int>(9)
#define REDSHIFT_MIN static_cast<real_prec>(0)
#define REDSHIFT_MAX static_cast<real_prec>(5)
#define N_REDSHIFT_BINS static_cast<int>(1)

#define DELTA_Z static_cast<real_prec>((REDSHIFT_MAX-REDSHIFT_MIN)/static_cast<real_prec>(N_REDSHIFT_BINS))
#endif

#ifndef _USE_REDSHIFT_MASK_
#define N_REDSHIFT_BINS 1
#endif

#else
#define N_REDSHIFT_BINS 1
#endif


