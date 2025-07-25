#include <fvar.hpp>
#include "SeapodymCohort.h"

/*void Hyperspace_projection(SeapodymCohort& sc, dvar_vector x);
void Taylor_derivative_test(const char* parfile);
void Hessian_comp(const char* parfile);*/
//void buffers_init(long int &mv, long int &mc, long int &mg, const bool grad_calc);
//void buffers_set(long int &mv, long int &mc, long int &mg);

/*!
\brief The first function to be executed.
\details
This is the main routine that calls upper-level functions such as
1. Parfile reading and initialization of model parameters and optimization variables;
2. Running the application in different regimes: 
   a) (default) running the model in forward (simulation) and backward (gradient computation) mode;
   b) simulation only with offline coupling with forage sub-model;
   c) running coupled simulation for tuna-forage model;
   d) computing Hessian;
   e) sensitivity analysis;
   f) computing 2d projection of likelihood function the pair of parameters (should be specified in parfile).
*/

SeapodymCohort* seapodym_cohort(const char* parfile, int cmp_regime, const bool reset_buffers, int cohort_id, gradient_structure& gs)
{
	time_t time_sec;
	time(&time_sec);
	const time_t time0 = time_sec;

	int out_hessian = 0;
	gradient_structure::set_USE_FOR_HESSIAN(out_hessian);

	cout << "\nstarting time: " << ctime(&time_sec) << endl;

	//read parfile
	SeapodymCohort* scp = new SeapodymCohort(parfile, cohort_id);
	SeapodymCohort& sc = *scp;

	//iniitalize variables of optimization
	const int nvar = sc.nvarcalc();
	independent_variables x(1, nvar);
	adstring_array x_names(1,nvar);

	sc.xinit(x, x_names);
	cout << "Total number of variables: " << nvar << '\n'<<'\n';

	//initialization of simulation
	sc.prerun_model();

	//the function is invoked in the coupled simulation only
	string tempparfile = "tempparfile.xml";
	string newparfile  = "newparfile.xml";

	//after minimization is finished one simulation will 
	//be run with estimated parameters; outputs will be saved
	gradient_structure::set_NO_DERIVATIVES();
	sc.run_cohort((dvar_vector)x, true);
	sc.write(newparfile.c_str());

	remove(tempparfile.c_str());

	//writes new parameters on the screen
	sc.param->outp_param(x_names,nvar);

	time(&time_sec);
	cout << "\nfinished time: " << ctime(&time_sec) << endl;

	time_t time2 = time_sec;
	double total_elapsed_time = (double)(time2-time0) / 60.0;
	cout << "\ntotal time: " << total_elapsed_time << " minutes" << endl;

	return scp;
}


void verify_identifier_string2(char* str1) //ASSUME str1 is not null
{

  // Back up the stream and read the number of bytes written in the
  // ``write function'' corresponding to this ``read function''
  long int num_bytes=strlen(str1);
  char* str = new char[num_bytes+1]; 
  str[num_bytes]='\0';
  gradient_structure::get_fp()->fread(str,num_bytes);
  if(strcmp(str1,str))
  {
    cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: \"" << str << "\" != \"" << str1 << "\"\n";
    ad_exit(1);
  }

  if (str) delete [] str;

}

int save_identifier_string2(char* str)
{
  int length=strlen(str);
  gradient_structure::get_fp()->fwrite(str,length);
  return 0;
}

void save_long_int_value(unsigned long int x)
{
  int num_bytes = sizeof(unsigned long int); 
  void* y = (void*)&x; 
  gradient_structure::get_fp()->fwrite(y,num_bytes);
}

unsigned long int restore_long_int_value(void)
{
  void* tmpout;
  int num_bytes = sizeof(unsigned long int);
  gradient_structure::get_fp()->fread(&tmpout,num_bytes);
  return (unsigned long int)tmpout;
}

