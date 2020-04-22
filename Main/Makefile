##################################################################################
##################################################################################
##################################################################################
# Makefile to compile COSMICATLASS
##################################################################################
##################################################################################
# objetivo: dependencias
#         instrucciones	
##################################################################################
.SUFFIXES:

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


NO_COLOR=\x1b[0m
OK_COLOR=\x1b[32;01m
ERROR_COLOR=\x1b[31;01m
WARN_COLOR=\x1b[33;01m

OK_STRING=$(OK_COLOR)[OK]$(NO_COLOR)
ERROR_STRING=$(ERROR_COLOR)[ERRORS]$(NO_COLOR)
WARN_STRING=$(WARN_COLOR)[WARNINGS]$(NO_COLOR)


##################################################################################
FLAGS    = -std=c++11 -fopenmp -fPIC
OPTIMIZE = -g3 -unroll 
# Use CFLAGS to point out warnings
CFLAGS   = -Wall -Wextra -pedantic -Wno-deprecated -Wno-write-strings  
DEBUGG   =
##################################################################################
# Use /net/vaina if the library is to be used in diode, diosa or dardo
#AUX_DIR = ../Lib/
AUX_DIR = /home/andres/data/Numerics/EUCLID/CODES/GitHub/
##################################################################################
##################################################################################
# math and gsl libraries
LIBS = -lm -lgsl -lgslcblas
#########################################################################
# fftw libraries
dir_FFTW =  /usr/lib64/
LIBS_FFTW = -L$(dir_FFTW) -lfftw3_threads -lfftw3  -lfftw3f -lfftw3f_threads
#########################################################################
# healpix libraries
dir_CFITSIO_LIB=$(AUX_DIR)cfitsio/
dir_HEALPIX_INC=$(AUX_DIR)Healpix_2.20/src/cxx/generic_gcc/include/
dir_HEALPIX_LIB=$(AUX_DIR)Healpix_2.20/src/cxx/generic_gcc/lib/ 
LIBS_HEALPIX = -L$(dir_CFITSIO_LIB)  -L$(dir_HEALPIX_LIB) -lhealpix_cxx -lcxxsupport -lpsht -lc_utils -lfftpack -lcfitsio -lgomp
# This flag is detected by make and uses it to compile the dependencies
CXXFLAGS  = -I$(dir_HEALPIX_INC) -g

#########################################################################
SOURCE_DIR = ../Source/
HEADER_DIR = ../Headers/
#########################################################################
#########################################################################
# These dependencies are compiled by make using the implicit variables and CXXFLAGS
# Here I did not write a rule for compiling every dependency. Make does it when asking these dependencies

#########################################################################
# BAM
OBJS = cosmicatlass.o $(SOURCE_DIR)cxxutils.o $(SOURCE_DIR)massFunctions.o\
$(SOURCE_DIR)Params.o $(SOURCE_DIR)Parameters_CosmoLib.o $(SOURCE_DIR)NumericalMethods.o $(SOURCE_DIR)CosmologicalFunctions.o\
$(SOURCE_DIR)ScreenOutput.o $(SOURCE_DIR)CoordinateSystem.o $(SOURCE_DIR)FftwFunctions.o\
$(SOURCE_DIR)PowerSpectrumF.o $(SOURCE_DIR)PowerSpectrumTH.o  $(SOURCE_DIR)CorrelationFunction.o\
$(SOURCE_DIR)DensityProfiles.o\
$(SOURCE_DIR)Statistics.o $(SOURCE_DIR)HOD.o $(SOURCE_DIR)BiasFunctions.o $(SOURCE_DIR)DnDz.o\
$(SOURCE_DIR)FileOutput.o $(SOURCE_DIR)McmcFunctions.o $(SOURCE_DIR)Catalog.o $(SOURCE_DIR)Bam.o\
$(SOURCE_DIR)Patchy.o $(SOURCE_DIR)Cwclass.o  

HDRS = $(HEADER_DIR)cxxutils.h $(HEADER_DIR)massFunctions.h $(HEADER_DIR)CosmoLib.h\
$(HEADER_DIR)Params.h $(HEADER_DIR)Parameters_CosmoLib.h $(HEADER_DIR)NumericalMethods.h $(HEADER_DIR)CosmologicalFunctions.h\
$(HEADER_DIR)ScreenOutput.h $(HEADER_DIR)CoordinateSystem.h $(HEADER_DIR)FftwFunctions.h\
$(HEADER_DIR)PowerSpectrumF.h $(HEADER_DIR)PowerSpectrumTH.h  $(HEADER_DIR)DensityProfiles.h\
$(HEADER_DIR)Statistics.h $(HEADER_DIR)HOD.h $(HEADER_DIR)BiasFunctions.h $(HEADER_DIR)DnDz.h\
$(HEADER_DIR)FileOutput.h $(HEADER_DIR)McmcFunctions.h $(HEADER_DIR)Catalog.h $(HEADER_DIR)Bam.h\
$(HEADER_DIR)Patchy.h $(HEADER_DIR)Cwclass.h 


SRCS_BAM = $(OBJS:.o=.cpp)
TARGET_BAM = cosmicatlass.exe

#########################################################################
#########################################################################
# COSMIC_LIB
OBJS_CL = cosmolib.o $(SOURCE_DIR)cxxutils.o \
$(SOURCE_DIR)Parameters_CosmoLib.o $(SOURCE_DIR)NumericalMethods.o $(SOURCE_DIR)CosmologicalFunctions.o\
$(SOURCE_DIR)ScreenOutput.o $(SOURCE_DIR)CoordinateSystem.o $(SOURCE_DIR)FftwFunctions.o\
$(SOURCE_DIR)PowerSpectrumTH.o  $(SOURCE_DIR)CorrelationFunction.o $(SOURCE_DIR)CosmoLib.o\
$(SOURCE_DIR)AngularPowerSpectrum.o $(SOURCE_DIR)DensityProfiles.o\
$(SOURCE_DIR)Statistics.o $(SOURCE_DIR)HOD.o $(SOURCE_DIR)BiasFunctions.o $(SOURCE_DIR)FileOutput.o

HDRS_CL = $(HEADER_DIR)cxxutils.h $(HEADER_DIR)CosmoLib.h\
$(HEADER_DIR)Parameters_CosmoLib.h $(HEADER_DIR)NumericalMethods.h $(HEADER_DIR)CosmologicalFunctions.h\
$(HEADER_DIR)ScreenOutput.h $(HEADER_DIR)PowerSpectrumTH.h $(HEADER_DIR)AngularPowerSpectrum.h\
$(HEADER_DIR)DensityProfiles.h $(HEADER_DIR)CorrelationFunction.h\
$(HEADER_DIR)Statistics.h $(HEADER_DIR)HOD.h $(HEADER_DIR)BiasFunctions.h $(HEADER_DIR)FileOutput.h 
SRCS_CL = $(OBJS_CL:.o=.cpp)
TARGET_COSMOLIB = cosmolib.exe
#########################################################################
#H-GAPS
##################################################################################
##################################################################################
##################################################################################	
# make will execute this rule by default
all: $(TARGET_BAM) $(TARGET_COSMOLIB) 
##################################################################################	
bam: $(TARGET_BAM)
##################################################################################	
cosmolib: $(TARGET_COSMOLIB) 

#########################################################################
$(TARGET_BAM): $(OBJS) $(HDRS)
	@echo "*******************************************"
	@echo "*******************************************"	
	@echo "Compiling COSMICATLASS"
	@date
	@echo "*******************************************"
	@echo "*******************************************"	
	$(CXX) $(FLAGS) $(SRCS_BAM) $(LIBS_HEALPIX) $(CXXFLAGS) $(LIBS_FFTW) $(LIBS) $(OPTIMIZE) -o $(TARGET_BAM)
	@echo "Done!"
params:
	./$(TARGET) -i input_deltastat_minerva.ini


#########################################################################
$(TARGET_COSMOLIB): $(OBJS_CL) $(HDRS_CL)
	@echo "*******************************************"
	@echo "*******************************************"	
	@echo "Compiling COSMOLIB"
	@date
	@echo "*******************************************"
	@echo "*******************************************"	
	$(CXX) $(FLAGS) $(SRCS_CL) $(LIBS_HEALPIX) $(CXXFLAGS) $(LIBS_FFTW) $(LIBS) $(OPTIMIZE) -o $(TARGET_COSMOLIB)
	@echo "Done!"


#########################################################################
clean:
	@echo "Cleaning:"
	rm -f  core.* *o *exe *~ ../Source/*o  ../Source/*~  ../Headers/*~
#########################################################################
