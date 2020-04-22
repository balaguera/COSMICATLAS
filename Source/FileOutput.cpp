/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// This file contains the methods of the class FFTW_FUNCTION used in the estimates     
// of clustering for EUCLID                                                             
// Developer:  Andres Balaguera Antolinez                                              
// e-mail:     abalant@gmail.com  
// Afiliation: INAF, OAR. Uni. Roma3
// 2013-2015
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



# include "../Headers/FileOutput.h"
# include <iostream>
# include <fstream>
# include <omp.h>

unsigned long FileOutput::read_file(string fname, vector< real_prec > &prop_ob, int NTHREADS)
{
  
  time_t start;
  time(&start);
  
  
  cout<<CYAN;
  So.message_screen("Reading input file",fname);
  So.message_screen("Nthreads = ",NTHREADS);

  int items_per_line = 0;
  unsigned nElem = 0;
  
  ifstream inputf (fname.c_str(), ifstream::in);
  if (!inputf) {
   std::cerr<<"Error in opening the input file " << fname << "!" << endl; exit(1);
  }
  
  // now counts how many doubles are in the first line
  // assuming that (1) that number is the same for all lines
  // and (2) that there is no header in the file
  
  nElem=m_countLines(inputf);
  //  cout<<YELLOW<<"Number of lines = "<<CYAN<<nElem <<endl;
  So.message_screen("Number of lines =",static_cast<int>(nElem));  
  string firstline;
  real_prec aa;
  if(!inputf.eof()) {
    getline (inputf, firstline);
    stringstream ss(firstline);
    while(ss>>aa)
      items_per_line++;
  }
  inputf.seekg (0, inputf.beg);
  So.message_screen("Number of columns =",items_per_line);  
  
  prop_ob.resize(nElem * items_per_line);
  
  size_t bufSize = min(static_cast<int>(m_getBuffer(inputf)), static_cast<int>(nElem));
  
  vector<string> tmplines(bufSize, "");
  int line_idx = 0;
  int ii;
  
  
  size_t doneElem = 0;
  size_t load = 0; // buffer to load
  volatile size_t nfail = 0; // bad catalog lines
  string line = ""; //line to read
  
  // parallel section to read file
  
#pragma omp parallel num_threads(NTHREADS)
  {
    
    
    while ( (doneElem < nElem) && (false == inputf.eof()) ){
#pragma omp single
      {
	load = min(static_cast<int>(bufSize), static_cast<int>(nElem - doneElem));
	for (size_t i=0; i<load; ++i)
	  getline(inputf, tmplines[i]);
	nfail = 0;
      }
      
      
      // -- process buffer - parallel -
#pragma omp for reduction (+:nfail)
      for (size_t i=0; i<load; ++i) 
	{
	  
	  // skip empty lines (none of them should be empty)
	  if (true == tmplines[i].empty())
	    {
	      nfail ++;
	      continue;
	    }
	  int jj = 0;
	  stringstream ss(tmplines[i]);
	  real_prec bb;
	  while (ss>>bb) prop_ob[(doneElem+i)*items_per_line + jj++] = static_cast<real_prec>(bb);
	}
      
#pragma omp single
      {
	size_t effective = load - nfail;
	doneElem += effective;
      }
    }
  } 
  
  
  inputf.clear();
  inputf.close();
  
  cout<<CYAN;
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (diff<60) cout <<CYAN<<diff<<" secs"<<endl;
  else if (diff<3600) cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
  else cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
  cout<<RESET;
  So.DONE();
  return static_cast<unsigned long>(nElem);
}



// *********************************************************************************************
// **********************************************************************************************
// *********************************************************************************************
// **********************************************************************************************

// Read file with one column
int FileOutput::read_file_one(string fname, vector<real_prec> &prop_ob){

  time_t start;
  time(&start);
  
  cout<<"Input file "<<fname<<endl;
  int number_of_lines = -1;
  
  //  ifstream inputf(fname.c_str());
  in.open(fname.c_str());
  if (!in) { 
   std::cerr<<"Error in opening the input file " << fname << "!" << endl; exit(1); 
  }
  
  real_prec aa;
  while(!in.eof())
    {
      number_of_lines++;
      in>>aa;
      prop_ob.push_back (aa);         //->fill array with all galaxies
    }
  in.close();
  
  cout<<CYAN;
  cout<<"Number of lines = " << number_of_lines <<RESET<<endl;
  
  cout<<RED;
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (diff<60) cout <<"Lapse: "<<diff<<" secs"<<endl;
  else if (diff<3600) cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
  else cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
  cout<<RESET;

  return number_of_lines;
}

// *********************************************************************************************
// **********************************************************************************************
// *********************************************************************************************
/*

// Read file with one column
int FileOutput::read_file_one(string fname, vector<float> &prop_ob){

  time_t start;
  time(&start);

  cout<<"Input file "<<fname<<endl;
  int number_of_lines = -1;

  //  ifstream inputf(fname.c_str());
  this->in.open(fname.c_str());
  if (!this->in) {
    std::cerr<<"Error in opening the input file " << fname << "!" << endl; exit(1);
  }
  
  double aa;
  while(!in.eof())
    {
      number_of_lines++;
      this->in>>aa;
      prop_ob.push_back (aa);         //->fill array with all galaxies
    }
  this->in.close();
  
  cout<<CYAN;
  cout<<"Number of lines = " << number_of_lines <<RESET<<endl;
  
  cout<<RED;
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (diff<60) cout <<"Lapse: "<<diff<<" secs"<<endl;
  else if (diff<3600) cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
  else cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
  cout<<RESET;

  return number_of_lines;
}

*/


// *********************************************************************************************
// **********************************************************************************************



void FileOutput::read_file_one_N(long N, string fname, vector<real_prec> &prop_ob)
{
  time_t start;
  time(&start);

  So.message_screen("Input file",fname," with ",N);

  //  ifstream inputf(fname.c_str());
  in.open(fname.c_str());
  if (!this->in) {
    std::cerr<<"Error in opening the input file " << fname << "!" << endl; exit(1);
  }
  
  prop_ob.resize(N,0);
  for(int i=0; i<N;++ i)
    this->in>>prop_ob[i];
  this->in.close();
  
  
  cout<<RED;
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (diff<60) cout <<"Lapse: "<<diff<<" secs"<<endl;
  else if (diff<3600) cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
  else cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
  cout<<RESET;

}




// *********************************************************************************************
// **********************************************************************************************

// THis
int FileOutput::read_file(string fname, vector<real_prec > &prop_ob){
  ScreenOutput So;

  time_t start;
  time(&start);

  cout<<"Input UUUU file "<<fname<<endl;

  char BUFF[20000];
  int items_per_line = 0;
  unsigned line_count = 0;

  //  ifstream inputf (fname.c_str());
  this->in.open(fname.c_str());

  if (!this->in)
    {
      std::cerr<<"Error in opening the input file " << fname << "!" << endl; exit(1);
    }
  
  cout<<CYAN;
  
  // now counts how many doubles are in the first line
  // assuming that (1) that number is the same for all lines
  // and (2) that there is no header in the file
  
  string firstline;
  real_prec aa;
  if(!this->in.eof()) {
    getline (this->in, firstline);
    stringstream ss(firstline);
    while(ss>>aa)
      items_per_line++;
  }
  this->in.seekg (0, this->in.beg);
  
  while (!this->in.eof())
    {
      this->in.read(BUFF, 20000);
      int count = this->in.gcount();
      for(int jj = 0; jj < count; jj++)
        if ( BUFF[jj] == '\n' )
          line_count++;
    }
  this->in.clear();
  this->in.close();
  
  prop_ob.resize(line_count * items_per_line);
  
  cout<<"Number of lines = " << prop_ob.size() / items_per_line <<endl;
  
  ifstream dina;
  dina.open(fname.c_str(), ios::in);
  assert(dina);



  string line2;
  double bb;



  int NTHREADS = 1;//omp_get_max_threads();

  cout<<CYAN;
  cout<<"Reading using "<<NTHREADS<<" thread(s)"<<endl;


  vector<string> line(1000*NTHREADS);

  // cada thread va a leer 1000 lineas. line entonces va a almacenar las lineas
  // para todos los threads y despues estos vendrana a buscarlas de modo paralelo

  bool sstop = false;
  int line_idx = 0;
  int ii;



#pragma omp parallel num_threads(NTHREADS)
  {
    double bb;

    while (!sstop) {  //mientras no se niege sstop, es decir, mientas sstop sea falso, dele

#pragma omp single
      {
        for(ii = 0; ii < line.size() && line_idx+ii < line_count; ii++)getline(dina,line[ii]);

        if(line_idx+ii == line_count) sstop = true; // si esto pasa, detiene todo porque ya llegamos al numero de lineas del file

      }


#pragma omp barrier



#pragma omp for
      for(int ll = 0; ll < ii; ll++)
        {
          int jj = 0;
          real_prec bb;
	  stringstream ss(line[ll]);
          while (ss>>bb) {prop_ob[(line_idx+ll)*items_per_line + jj++] = bb;}
        }


#pragma omp single
      line_idx += ii;   //line_idx counts the number of lines

    }

  }

  cout<<"Number of columns = "<< items_per_line <<endl;
  if(input_type=="catalog") if(items_per_line < 3) So.error_ncolumns(fname);
  dina.clear();
  dina.close();

  cout<<RED;
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (diff<60) cout <<"Lapse: "<<diff<<" secs"<<endl;
  else if (diff<3600) cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
  else cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
  cout<<RESET;
  return items_per_line;
}





// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************

/* This is the old version, in which a 2sim container is returned. The rest of the code
 * is not prepared for this 2dim, instead use the functiona above
 */

void FileOutput::read_file(string fname, vector< vector<real_prec> > &prop_ob){

  So.message_screen("Input file using FileOutput, OLD version ");

  time_t start;
  time(&start);

  string line;
  this->in.open(fname.c_str(), ios::in);
  assert(this->in);
  real_prec aa;
  cout<<CYAN;
  cout<<"Reading input file  "<<fname<<endl;
  while(getline (this->in, line)){
    stringstream ss(line);
    vector <real_prec> nc;
    while(ss>>aa)nc.push_back (aa); //->vector for each galaxy with properties
    prop_ob.push_back (nc);         //->fill array with all galaxies
  }
  
  cout<<"Number of lines = "<<prop_ob.size()<<endl;
  cout<<"Number of columns = "<<prop_ob[0].size()<<endl;
  if(prop_ob[0].size()<3)So.error_ncolumns(fname);
  this->in.clear();
  this->in.close();
  
  cout<<RED;
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (diff<60) cout <<"Lapse: "<<diff<<" secs"<<endl;
  else if (diff<3600) cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
  else cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
  cout<<RESET;
}


// **********************************************************************************************
// **********************************************************************************************
void FileOutput::write_to_file(string fname, vector<real_prec>& kve)
{
  this->out.open(fname.c_str() , ios::out);
    So.message_output_file(fname, kve.size(),1);

    for(int i=0;i<kve.size();++i)
    this->out<<kve[i]<<endl;
  this->out.close();
  So.DONE();
}

// **********************************************************************************************
// **********************************************************************************************
void FileOutput::write_to_file_i(string fname, vector<ULONG>&kve)
{
  this->out.open(fname.c_str() , ios::out);
    So.message_output_file(fname, kve.size(),2);

    for(int i=0;i<kve.size();++i)
    this->out<<i<<"  "<<kve[i]<<endl;
  this->out.close();
    So.DONE();
}


void FileOutput::write_to_file_i(string fname, vector<int>&kve)
{
  this->out.open(fname.c_str() , ios::out);
    So.message_output_file(fname, kve.size(),2);

    for(int i=0;i<kve.size();++i)
    this->out<<i<<"  "<<kve[i]<<endl;
  this->out.close();
    So.DONE();
}


// **********************************************************************************************
void FileOutput::write_to_file(string fname, vector<real_prec>& kve, vector<real_prec> &pk)
{
    So.message_output_file(fname, kve.size(),2);

  this->out.open(fname.c_str() , ios::out); 
  for(int i=0;i<kve.size();++i)out<<kve[i]<<"\t"<<pk[i]<<endl; 
  this->out.close();
  So.DONE();
}


// **********************************************************************************************
// In case we are in single precision, we might still want to write the gsl_real vectors (double prec)

// **********************************************************************************************
// **********************************************************************************************
void FileOutput::write_to_file(string fname, vector<real_prec>& kve, vector<real_prec> &pk, vector<real_prec> &si)
{
  this->out.open(fname.c_str() , ios::out); 
    So.message_output_file(fname, kve.size(),2);

    for(int i=0;i<kve.size();++i)out<<kve[i]<<"\t"<<pk[i]<<"\t"<<si[i]<<endl;
  this->out.close(); 
    So.DONE();

}
// **********************************************************************************************
// **********************************************************************************************

void FileOutput::write_to_file2(string fname, vector<real_prec>& kve, vector<real_prec> &pk, vector<int> &si)
{
  So.message_output_file(fname, kve.size(),3);
  this->out.open(fname.c_str() , ios::out);
  for(int i=0;i<kve.size();++i)this->out<<kve[i]<<"\t"<<pk[i]<<"\t"<<si[i]<<endl; 
  this->out.close();
  So.DONE();

}

// **********************************************************************************************
// **********************************************************************************************

void FileOutput::write_to_file(string fname, vector<real_prec>& kve, vector<real_prec> &pk, vector<real_prec> &si, vector<real_prec> &sai)
{
  this->out.open(fname.c_str() , ios::out); 
    So.message_output_file(fname, kve.size(),4);
    for(int i=0;i<kve.size();++i)this->out<<kve[i]<<"\t"<<pk[i]<<"\t"<<si[i]<<"\t"<<sai[i]<<endl;
  this->out.close();
  So.DONE();
}

// **********************************************************************************************
// **********************************************************************************************

void FileOutput::write_to_file(string fname, vector<real_prec>& kve, vector<real_prec> &pk, vector<real_prec> &pkk,vector<real_prec> &pkkk, vector<real_prec> &si)
{
  this->out.open(fname.c_str() , ios::out); 
  So.message_output_file(fname, kve.size(),5);
   for(int i=0;i<kve.size();++i)this->out<<kve[i]<<"\t"<<pk[i]<<"\t"<<pkk[i]<<"\t"<<pkkk[i]<<"\t"<<si[i]<<endl;
  this->out.close(); 
  So.DONE();
}


// **********************************************************************************************

void FileOutput::write_to_file(string fname, vector<real_prec>& kve, vector<real_prec> &pk, vector<real_prec> &pkk,vector<real_prec> &pkkk, vector<real_prec> &si, vector<real_prec>& nmod)
{
  this->out.open(fname.c_str() , ios::out); 
  So.message_output_file(fname, kve.size(),7);
  for(int i=0;i<kve.size();++i)
      this->out<<kve[i]<<"\t"<<pk[i]<<"\t"<<pkk[i]<<"\t"<<pkkk[i]<<"\t"<<si[i]<<"\t"<<nmod[i]<<endl;
  So.DONE();
  this->out.close(); 
}

// **********************************************************************************************
void FileOutput::write_to_file(string fname, vector<real_prec>& kve, vector<real_prec> &pk, vector<real_prec> &pkk,vector<real_prec> &pkkk, vector<real_prec> &si, vector<int>& nmod)
{
  
  this->out.open(fname.c_str() , ios::out);
    So.message_output_file(fname, kve.size(), 6);

  for(int i=0;i<kve.size();++i)
    this->out<<kve[i]<<"\t"<<pk[i]<<"\t"<<pkk[i]<<"\t"<<pkkk[i]<<"\t"<<si[i]<<"\t"<<nmod[i]<<endl;
  
  this->out.close();
  So.DONE();
}

// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
void FileOutput::write_to_file(string fname, vector<real_prec> &kvx, vector<real_prec>&kvy,vector< vector<real_prec> > &pk)
{
  So.message_output_file(fname, kvx.size()*kvy.size());
  this->out.open(fname.c_str() , ios::out); 
  for(int i=0;i<kvx.size();++i)
    for(int j=0;j<kvy.size();++j)
      this->out<<kvx[i]<<"\t"<<kvy[j]<<"\t"<<pk[i][j]<<endl; 
  this->out.close(); 
  So.DONE();
}

// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
void FileOutput::write_to_file(string fname, vector<int> &kvx, vector<int>&kvy,vector< vector<real_prec> > &pk)
{
  this->out.open(fname.c_str() , ios::out); 
  So.message_output_file(fname, kvx.size()*kvy.size());
    for(int i=0;i<kvx.size();++i)
    for(int j=0;j<kvy.size();++j)
      this->out<<kvx[i]<<"\t"<<kvy[j]<<"\t"<<pk[i][j]<<endl; 
  this->out.close();
    So.DONE();

}

// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************



void FileOutput::write_to_file(string fname, vector<int> &kvx, vector<real_prec>&kvy, vector<real_prec> &pk, vector<real_prec>&pkd)
{
  this->out.open(fname.c_str() , ios::out); 
  So.message_output_file(fname, kvy.size(), 3);

  for(int i=0;i<kvx.size();++i)
    this->out<<kvx[i]<<"\t"<<kvy[i]<<"\t"<<pk[i]<<"\t"<<pkd[i]<<endl; 
  this->out.close(); 
  So.DONE();
  return;
}


// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************


void FileOutput::write_to_file(string fname, vector<real_prec> &kvx, vector<int>&kvy,vector< vector<real_prec> > &pk)
{
  this->out.open(fname.c_str() , ios::out); 
  So.message_output_file(fname, kvx.size()*kvy.size());

    for(int i=0;i<kvx.size();++i)
    for(int j=0;j<kvy.size();++j)
      this->out<<kvx[i]<<"\t"<<kvy[j]<<"\t"<<pk[i][j]<<endl; 
  this->out.close(); 
  So.DONE();
}


// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************


void FileOutput::write_to_file(string fname, vector<real_prec>& kve, vector<real_prec> &pl, int n_columns)
{
  So.message_output_file(fname, kve.size(), pl.size());
  this->out.open(fname.c_str() , ios::out);
  for(int i=0;i<kve.size();++i){
    this->out<<kve[i]<<"    ";
    for(int l=0;l<pl.size()/n_columns;l++)
      this->out<<pl[l*n_columns+i]<<"   ";
    this->out<<endl;
  }
  this->out.close(); 
  cout<<CYAN; 
  cout<<RESET;
  So.DONE();
  return;
}

// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************


void FileOutput::write_to_file(string fname, vector<real_prec> &kve, vector<real_prec> &bis, vector<int> &mod)
{

  this->out.open(fname.c_str() , ios::out); 
  this->out.precision(12); 
  this->out.setf(ios::showpoint); 
  this->out.setf(ios::scientific); 

#define index(a,b,c) ((a)*i_n*i_n + (b)*i_n + (c))
  So.message_screen("Writting output file",fname);
  int nc=0;
  int i_n=kve.size();
  for(int i=0;i<i_n;++i){
    for(int j=0;j<i_n;++j){
      for(int k=0;k<i_n;++k){
  	if(i<=j && j<=k){
  	  if(mod[index(i,j,k)]!=0){
	    this->out<<i<<"\t"<<j<<"\t"<<k<<"\t"<<bis[index(i,j,k)]<<"\t\t"<<mod[index(i,j,k)]<<endl;
  	  }
  	}
      }
    }
  }

  this->out.close(); 

  cout<<BOLDGREEN;
  cout<<RESET;
  So.DONE();
  return;
}

// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************

void FileOutput::write_to_file(string fname, vector<real_prec> &kve, vector<real_prec> &bis, vector<real_prec> &sn_bis, vector<int> &mod)
{
  
  ScreenOutput So;
  this->out.open(fname.c_str() , ios::out);
  this->out.precision(12); 
  this->out.setf(ios::showpoint); 
  this->out.setf(ios::scientific); 

#define index(a,b,c) ((a)*i_n*i_n + (b)*i_n + (c))
  So.message_screen("Writting output file",fname);

  int nc=0;
  int i_n=kve.size();
  for(int i=0;i<i_n;++i){
    for(int j=0;j<i_n;++j){
      for(int k=0;k<i_n;++k){
	//if(i<=j && j<=k){
	  if(mod[index(i,j,k)]!=0){
	    this->out<<i<<"\t"<<j<<"\t"<<k<<"\t"<<bis[index(i,j,k)]<<"\t"<<sn_bis[index(i,j,k)]/mod[index(i,j,k)]<<"\t"<<mod[index(i,j,k)]<<endl;
	  }
	  //}
      }
    }
  }

  this->out.close(); 
  cout<<BOLDGREEN;
  cout<<RESET;
  So.DONE();
  return; 
}



//##################################################################################
void FileOutput::read_array(string fname,float *OUT,ULONG N)
{
  ScreenOutput So;
  fftw_array<float> dummy(N);
  this->So.message_screen("Reading binary file", fname);
  this->inStream.open(fname.data(),file_is_natural); 
  assert(this->inStream.is_open()); 
  this->inStream.get(dummy.data,N); 
  this->inStream.close(); 
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)
    OUT[i]=dummy[i];
  So.DONE();
  return ;
}



//##################################################################################
void FileOutput::read_array(string fname,vector<real_prec>&OUT)
{
  ULONG N=OUT.size();
  ScreenOutput So;

  fftw_array<float> dummy(N);
  this->So.message_screen("Reading binary file", fname);
  this->inStream.open(fname.data(),file_is_natural); 
  assert(this->inStream.is_open()); 
  this->inStream.get(dummy.data,N); 
  this->inStream.close(); 
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)
    OUT[i]=static_cast<real_prec>(dummy[i]);
  So.DONE();
  return ;
}

//##################################################################################
void FileOutput::read_array(string fname,vector<ULONG>&OUT)
{
  ULONG N=OUT.size();
  fftw_array<float> dummy(N);
  this->So.message_screen("Reading binary file", fname);
  this->inStream.open(fname.data(),file_is_natural); 
  assert(this->inStream.is_open()); 
  this->inStream.get(dummy.data,N); 
  this->inStream.close(); 
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)
    OUT[i]=static_cast<ULONG>(dummy[i]);
  So.DONE();
  return ;
}

//##################################################################################
bool FileOutput::exists(string fname)
{
  fname=fname+".dat";
  this->So.message_screen("Verifying existing file", fname);
  this->inStream.open(fname.data(), file_is_natural); 
  bool open = this->inStream.is_open();
  this->inStream.close(); 
  if(false==open)
    this->So.message_screen("File not found");
  else
    this->So.message_screen("File found");
    
  return open;
}



//##################################################################################
//##################################################################################
//##################################################################################

void FileOutput::write_array(string fname,float *A_rm,ULONG N)
{
  fftw_array<float> dummy(N);
  
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)
    dummy[i]=A_rm[i];
  
  string FNAME=fname+string(".dat");
  this->So.message_screen("Writting in binary file", FNAME);
  bofstream outStream(FNAME.data(),file_is_natural);
  assert(outStream.is_open());
  outStream.put(dummy.data,N);
  outStream.close();
  So.DONE();
  return;
}

//##################################################################################
//##################################################################################


void FileOutput::write_array(string fname, vector<real_prec>&out)
{

  ULONG N=out.size();
  fftw_array<float> dummy(N);


#pragma omp parallel for
  for(ULONG i=0;i<N;i++)
    dummy[i]=static_cast<float>(out[i]);
  
  string FNAME=fname+string(".dat");
  this->So.message_screen("Writting in binary file", FNAME);
  bofstream outStream(FNAME.data(),file_is_natural);
  assert(outStream.is_open());
  outStream.put(dummy.data,N);
  outStream.close();
  So.DONE();
  return;
}


//##################################################################################
//##################################################################################


void FileOutput::write_array(string fname, vector<ULONG>&out)
{
  ULONG N=out.size();
  fftw_array<float> dummy(N);
  ScreenOutput So;

#pragma omp parallel for
  for(ULONG i=0;i<N;i++)
    dummy[i]=static_cast<float>(out[i]);
  
  string FNAME=fname+string(".dat");
  this->So.message_screen("Writting binary file", FNAME);
  bofstream outStream(FNAME.data(),file_is_natural);
  assert(outStream.is_open());
  outStream.put(dummy.data,N);
  outStream.close();
  So.DONE();
  return;
}

//##################################################################################
//##################################################################################
ULONG FileOutput::count_lines(string fname)
{
 ifstream inStream;
 inStream.open(fname.data());
 assert(inStream.is_open());
 string line;
 ULONG linenum = 0;
 while (getline (inStream, line))
   {
     linenum++;     
   }
 inStream.close();
 
 return linenum;
}



// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************


#ifdef SINGLE_PREC
// **********************************************************************************************
// In case we are in single precision, we might still want to write the gsl_real vectors (double prec)
// **********************************************************************************************
// **********************************************************************************************
void FileOutput::write_to_file(string fname,  vector<gsl_real>& kve)
{
  this->out.open(fname.c_str() , ios::out);
    So.message_output_file(fname, kve.size(),1);

    for(int i=0;i<kve.size();++i)
    this->out<<kve[i]<<endl;
  this->out.close();
  So.DONE();
}

// **********************************************************************************************
void FileOutput::write_to_file(string fname, vector<gsl_real>& kve, vector<gsl_real> &pk)
{
    So.message_output_file(fname, kve.size(),2);

  this->out.open(fname.c_str() , ios::out); 
  for(int i=0;i<kve.size();++i)out<<kve[i]<<"\t"<<pk[i]<<endl; 
  this->out.close();
  So.DONE();
}


// **********************************************************************************************
void FileOutput::write_to_file(string fname, vector<gsl_real>& kve, vector<gsl_real> &pk, vector<gsl_real> &si)
{
  this->out.open(fname.c_str() , ios::out); 
    So.message_output_file(fname, kve.size(),2);

    for(int i=0;i<kve.size();++i)out<<kve[i]<<"\t"<<pk[i]<<"\t"<<si[i]<<endl;
  this->out.close(); 
    So.DONE();

}
// **********************************************************************************************
// **********************************************************************************************

void FileOutput::write_to_file2(string fname, vector<gsl_real>& kve, vector<gsl_real> &pk, vector<int> &si)
{
  So.message_output_file(fname, kve.size(),3);
  this->out.open(fname.c_str() , ios::out);
  for(int i=0;i<kve.size();++i)this->out<<kve[i]<<"\t"<<pk[i]<<"\t"<<si[i]<<endl; 
  this->out.close();
  So.DONE();

}

// **********************************************************************************************
// **********************************************************************************************

void FileOutput::write_to_file(string fname, vector<gsl_real>& kve, vector<gsl_real> &pk, vector<gsl_real> &si, vector<gsl_real> &sai)
{
  this->out.open(fname.c_str() , ios::out); 
    So.message_output_file(fname, kve.size(),4);
    for(int i=0;i<kve.size();++i)this->out<<kve[i]<<"\t"<<pk[i]<<"\t"<<si[i]<<"\t"<<sai[i]<<endl;
  this->out.close();
  So.DONE();
}

// **********************************************************************************************
// **********************************************************************************************

void FileOutput::write_to_file(string fname, vector<gsl_real>& kve, vector<gsl_real> &pk, vector<gsl_real> &pkk,vector<gsl_real> &pkkk, vector<gsl_real> &si)
{
  this->out.open(fname.c_str() , ios::out); 
  So.message_output_file(fname, kve.size(),5);
   for(int i=0;i<kve.size();++i)this->out<<kve[i]<<"\t"<<pk[i]<<"\t"<<pkk[i]<<"\t"<<pkkk[i]<<"\t"<<si[i]<<endl;
  this->out.close(); 
  So.DONE();
}


// **********************************************************************************************

void FileOutput::write_to_file(string fname, vector<gsl_real>& kve, vector<gsl_real> &pk, vector<gsl_real> &pkk,vector<gsl_real> &pkkk, vector<gsl_real> &si, vector<gsl_real>& nmod)
{
  this->out.open(fname.c_str() , ios::out); 
  So.message_output_file(fname, kve.size(),7);
  for(int i=0;i<kve.size();++i)
      this->out<<kve[i]<<"\t"<<pk[i]<<"\t"<<pkk[i]<<"\t"<<pkkk[i]<<"\t"<<si[i]<<"\t"<<nmod[i]<<endl;
  So.DONE();
  this->out.close(); 
}

// **********************************************************************************************
void FileOutput::write_to_file(string fname, vector<gsl_real>& kve, vector<gsl_real> &pk, vector<gsl_real> &pkk,vector<gsl_real> &pkkk, vector<gsl_real> &si, vector<int>& nmod)
{
  
  this->out.open(fname.c_str() , ios::out);
    So.message_output_file(fname, kve.size(), 6);

  for(int i=0;i<kve.size();++i)
    this->out<<kve[i]<<"\t"<<pk[i]<<"\t"<<pkk[i]<<"\t"<<pkkk[i]<<"\t"<<si[i]<<"\t"<<nmod[i]<<endl;
  
  this->out.close();
  So.DONE();
}

// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
void FileOutput::write_to_file(string fname, vector<gsl_real> &kvx, vector<gsl_real>&kvy,vector< vector<gsl_real> > &pk)
{
  So.message_output_file(fname, kvx.size()*kvy.size());
  this->out.open(fname.c_str() , ios::out); 
  for(int i=0;i<kvx.size();++i)
    for(int j=0;j<kvy.size();++j)
      this->out<<kvx[i]<<"\t"<<kvy[j]<<"\t"<<pk[i][j]<<endl; 
  this->out.close(); 
  So.DONE();
}

// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
void FileOutput::write_to_file(string fname, vector<int> &kvx, vector<int>&kvy,vector< vector<gsl_real> > &pk)
{
  this->out.open(fname.c_str() , ios::out); 
  So.message_output_file(fname, kvx.size()*kvy.size());
    for(int i=0;i<kvx.size();++i)
    for(int j=0;j<kvy.size();++j)
      this->out<<kvx[i]<<"\t"<<kvy[j]<<"\t"<<pk[i][j]<<endl; 
  this->out.close();
    So.DONE();

}

// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************



void FileOutput::write_to_file(string fname, vector<int> &kvx, vector<gsl_real>&kvy, vector<gsl_real> &pk, vector<gsl_real>&pkd)
{
  this->out.open(fname.c_str() , ios::out); 
  So.message_output_file(fname, kvy.size(), 3);

  for(int i=0;i<kvx.size();++i)
    this->out<<kvx[i]<<"\t"<<kvy[i]<<"\t"<<pk[i]<<"\t"<<pkd[i]<<endl; 
  this->out.close(); 
  So.DONE();
  return;
}


// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************


void FileOutput::write_to_file(string fname, vector<gsl_real> &kvx, vector<int>&kvy,vector< vector<gsl_real> > &pk)
{
  this->out.open(fname.c_str() , ios::out); 
  So.message_output_file(fname, kvx.size()*kvy.size());

    for(int i=0;i<kvx.size();++i)
    for(int j=0;j<kvy.size();++j)
      this->out<<kvx[i]<<"\t"<<kvy[j]<<"\t"<<pk[i][j]<<endl; 
  this->out.close(); 
  So.DONE();
}


// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************


void FileOutput::write_to_file(string fname, vector<gsl_real>& kve, vector<gsl_real> &pl, int n_columns)
{
  So.message_output_file(fname, kve.size(), pl.size());
  this->out.open(fname.c_str() , ios::out);
  for(int i=0;i<kve.size();++i){
    this->out<<kve[i]<<"    ";
    for(int l=0;l<pl.size()/n_columns;l++)
      this->out<<pl[l*n_columns+i]<<"   ";
    this->out<<endl;
  }
  this->out.close(); 
  cout<<CYAN; 
  cout<<RESET;
  So.DONE();
  return;
}

// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************


void FileOutput::write_to_file(string fname, vector<gsl_real> &kve, vector<gsl_real> &bis, vector<int> &mod)
{

  this->out.open(fname.c_str() , ios::out); 
  this->out.precision(12); 
  this->out.setf(ios::showpoint); 
  this->out.setf(ios::scientific); 

#define index(a,b,c) ((a)*i_n*i_n + (b)*i_n + (c))
  So.message_screen("Writting output file",fname);
  int nc=0;
  int i_n=kve.size();
  for(int i=0;i<i_n;++i){
    for(int j=0;j<i_n;++j){
      for(int k=0;k<i_n;++k){
  	if(i<=j && j<=k){
  	  if(mod[index(i,j,k)]!=0){
	    this->out<<i<<"\t"<<j<<"\t"<<k<<"\t"<<bis[index(i,j,k)]<<"\t\t"<<mod[index(i,j,k)]<<endl;
  	  }
  	}
      }
    }
  }

  this->out.close(); 

  cout<<BOLDGREEN;
  cout<<RESET;
  So.DONE();
  return;
}

// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************
// **********************************************************************************************

void FileOutput::write_to_file(string fname, vector<gsl_real> &kve, vector<gsl_real> &bis, vector<gsl_real> &sn_bis, vector<int> &mod)
{
  
  ScreenOutput So;
  this->out.open(fname.c_str() , ios::out);
  this->out.precision(12); 
  this->out.setf(ios::showpoint); 
  this->out.setf(ios::scientific); 

#define index(a,b,c) ((a)*i_n*i_n + (b)*i_n + (c))
  So.message_screen("Writting output file",fname);

  int nc=0;
  int i_n=kve.size();
  for(int i=0;i<i_n;++i){
    for(int j=0;j<i_n;++j){
      for(int k=0;k<i_n;++k){
	//if(i<=j && j<=k){
	  if(mod[index(i,j,k)]!=0){
	    this->out<<i<<"\t"<<j<<"\t"<<k<<"\t"<<bis[index(i,j,k)]<<"\t"<<sn_bis[index(i,j,k)]/mod[index(i,j,k)]<<"\t"<<mod[index(i,j,k)]<<endl;
	  }
	  //}
      }
    }
  }

  this->out.close(); 
  cout<<BOLDGREEN;
  cout<<RESET;
  So.DONE();
  return; 
}






#endif
