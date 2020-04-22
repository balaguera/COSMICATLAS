#ifndef _SCREEN_OUTPUT_
#define _SCREEN_OUTPUT_

# include <ctime>
# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include <cassert>
# include "def.h"
# include "CosmologicalFunctions.h"
using namespace std;


//*****************************************************************************************/
// CLASS 
//*****************************************************************************************/

class ScreenOutput{
  
 private:
  time_t initial_time;
  
 public:
  ScreenOutput(){}
  
 ScreenOutput(time_t _init_t): initial_time(_init_t){}
  ~ScreenOutput(){}
  void message(string);
  void welcome_message();
  void welcome_message_cl();
  void welcome_message_fb();
  void welcome_message_yama();
  void welcome_message_bispectrum();
  void welcome_message_bispectrum_fast();
  void message_interrupt();
  void done(string);
  void enter(string);  
  void error_ncolumns(string);
  void ending_message();
  void write_fftw_parameters(void *);
  void write_parameters_estimator(void *);
  void write_cosmo_parameters(void *, void *);
  void write_cosmo_parameters(void *);
  void write_parameters_b_estimator(void *);
  void write_yam_parameters(void *, bool);
  void comp_time(time_t, unsigned long, unsigned long);
  void message_screen(string ss, double d2, time_t time);
  void message_screen(string ss, double d2, time_t time, time_t time2);
  void message_warning(string ss);
  void message_screen_flush(string ss, int s2);
  void message_screen_flush(string ss, real_prec s2, string sa, real_prec s3);
  void message_warning(string ss, ULONG);
  int message_error(string ss);
  int message_error(string ss, double a);
  void message_error(string ss, double a, string sa);
  void message_screen(string ss);
  void message_screen(string ss, int i, string sa,double s2);
  void message_screen(string ss, int s2);
  void message_screen(string ss, ULONG s2);
  void message_screen(string ss, double s2);
  void message_screen(string ss, string s2);
  void message_screen(string ss, double d2, string s2);
  void message_screen(string ss, string s2, string s3, double d3);
  void message_time2(time_t start_all);
  void message_time_mock(time_t start_all);
  void message_time(time_t start_all);
  void message_screen(string ss, double a, string s2, string f);
  void usage(string );
  void author();
  void message(time_t);
  void message_output_file(string s, ULONG l);
  void message_output_file(string s, int l, int n);
  void DONE();
};

#endif
