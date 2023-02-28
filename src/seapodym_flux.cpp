#include <fvar.hpp>
#include "SeapodymCoupled.h"
//#include <functional> // to pass function as an argument 

string get_path(const char* full_path);
void Hyperspace_projection(SeapodymCoupled& sc, dvar_vector x);
void Sensitivity_analysis(const char* parfile, const int sftype);
void Hessian_comp(const char* parfile);
void buffers_init(long int &mv, long int &mc, long int &mg, const bool grad_calc);
void buffers_set(long int &mv, long int &mc, long int &mg);

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
int seapodym_flux(const char* parfile, int cmp_regime, int FLAG, const bool reset_buffers)
{
	time_t time_sec;
	time(&time_sec);
	const time_t time0 = time_sec;

	//-----Memory stack sizes for dvariables and derivatives storage------
	gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
	long int gradstack_buffer, cmpdif_buffer, gs_var_buffer;
	bool grad_calc = false;
	if (cmp_regime==-1 || cmp_regime==2) grad_calc = true;
	buffers_init(gs_var_buffer, gradstack_buffer, cmpdif_buffer, grad_calc);
	if (reset_buffers)
		buffers_set(gs_var_buffer, gradstack_buffer, cmpdif_buffer);

	gradient_structure::set_GRADSTACK_BUFFER_SIZE(gradstack_buffer);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(cmpdif_buffer);
	gradient_structure gs(gs_var_buffer);
	//--------------------------------------------------------------------		

	cout << "\nstarting time: " << ctime(&time_sec) << endl;

	//redirect to other run modes
	if (cmp_regime == 2){
		Hessian_comp(parfile);
		return 0;
	} else if (cmp_regime == 3){
		Sensitivity_analysis(parfile,FLAG);
		return 0;
	}

	//Main model class instance
	SeapodymCoupled sc(parfile);
	
	//Initialize variable parametres
	const int nvar = sc.nvarcalc();
	independent_variables x(1, nvar);
	adstring_array x_names(1,nvar);

	sc.xinit(x, x_names);
	cout << "Total number of variables: " << nvar << '\n'<<'\n';

	//function minimizer class
	fmm fmc(nvar);

	//flags for function minimizer
	fmc.iprint = 1;
	fmc.crit = sc.get_crit();//0.1
	fmc.imax = 30;
	fmc.scroll_flag = 1;
	fmc.ifn = 0;
	fmc.maxfn = sc.get_maxfn();//2000;
	if (fmc.maxfn <= 0)
		fmc.ireturn = -1; 

	//if this flag is 0 then the gradient will not be computed
	int compute_gradient = 1;

	//coupled simulation for tuna-predator model, no minimization
	if (cmp_regime == 0 || sc.param->flag_coupling){
		if (sc.param->flag_coupling && compute_gradient!=0){
			cout << "Warning: COUPLING regime is ON, no parameter estimation will be performed!" << endl;
		}
		sc.param->set_gradcalc(false);
		compute_gradient = 0;
	}

	//simulation regime to compute 2d projection of likelihood function
	//over any two variable parameters (should be specified through parfile)
	if (cmp_regime == 1){
		//sc.param->set_gradcalc(false);
		gradient_structure::set_NO_DERIVATIVES();
		Hyperspace_projection(sc,(dvar_vector)x);
		return 0;
	}

	double likelihood = 0;
	double elapsed_time = 0;

	//gradient vector allocation and initialization
	dvector g(1, nvar); g.initialize();

	int idx = 0;
	int itr = 0;
	ios::sync_with_stdio();

	//initialization of simulation parameters
	sc.OnRunFirstStep();
	//the function is invoked in the coupled simulation only
	string tempparfile = "tempparfile.xml";
	string newparfile  = "newparfile.xml";


	//clock_t time1 = clock();
	time(&time_sec);
	time_t time1 = time_sec;
	//'run_flux' runs the model in forward (simulation) mode
	//'gradcalc' runs backward (adjoint) mode
	//'save_statistics' stores current information on function minimization
	if (compute_gradient){
		string dirout = get_path(parfile);
		tempparfile = dirout +"/tempparfile.xml";
		newparfile  = dirout +"/newparfile.xml";
		cout << "\nentering minimization loop" << endl;
		while (fmc.ireturn >= 0) {
			//call function minimizer
			fmc.fmin(likelihood, x, g);

			//update the statistics.out file if the solution has been improved
			int itn = fmc.itn;
			if ((itn==itr+1 && (idx>=1))){
				time(&time_sec);
				time_t time2 = time_sec;
				elapsed_time += (double)(time2-time1)/60.0;
				time1 = time2;
				sc.save_statistics(dirout,x_names,likelihood,g,elapsed_time,idx-1,itr,nvar);
				sc.write(tempparfile.c_str());
				itr = itn;
			}
			//reset control parameters and run the model and its adjoint
			if (fmc.ireturn > 0) {
				likelihood = sc.run_flux((dvar_vector)x);
				gradcalc(nvar,g);
				cout << "function evaluation " << idx++ << endl;
			}
		}
		//test: write the optimisation stats on exit with final values of like and parameters
		//when it exists the 'while' loop right after gradcalc returns exit code = 1
		if (fmc.ireturn==-1){
			time(&time_sec);
			time_t time2 = time_sec;
			elapsed_time += (double)(time2-time1)/60.0;
			sc.save_statistics(dirout,x_names,likelihood,g,elapsed_time,idx-1,itr,nvar);
		}
		
	}

	//after minimization is finished one simulation will 
	//be run with estimated parameters; outputs will be saved
	gradient_structure::set_NO_DERIVATIVES();
	sc.run_flux((dvar_vector)x, true);
	sc.write(newparfile.c_str());

	remove(tempparfile.c_str());

	//writes new parameters on the screen
	sc.param->outp_param(x_names,nvar);

	time(&time_sec);
	cout << "\nfinished time: " << ctime(&time_sec) << endl;

	time_t time2 = time_sec;
	double total_elapsed_time = (double)(time2-time0) / 60.0;
	cout << "\ntotal time: " << total_elapsed_time << " minutes" << endl;

	return 0;
}

double run_model(SeapodymCoupled& sc, dvar_vector x, dvector& g, const int nvar)
{
	double like = 0.0;
	like = sc.run_flux((dvar_vector)x);
	gradcalc(nvar,g); 
	return like;
}

int seapodym_phases(const char* parfile)
{
	time_t time_sec;
	time(&time_sec);
	const time_t time0 = time_sec;

	//initialize stack sizes for derivative information storage
	gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
	//-----ONLY OPTIMIZATION MODE HERE------
	long int gs_var_buffer    =  500000000L;
        long int gradstack_buffer =  110000000L;
        long int cmpdif_buffer    = 4000000000L;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(gradstack_buffer);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(cmpdif_buffer);
	gradient_structure gs(gs_var_buffer);

	cout << "\nstarting time: " << ctime(&time_sec) << endl;

	//read parfile
	SeapodymCoupled sc(parfile);

	//iniitalize variables of optimization
	const int nvar = sc.nvarcalc();
	adstring_array x_names;
	string* pnames;
	pnames=new string[nvar];
	
	//read parameter sensitivities from sa_res_all.txt and
	//initialize a table 4 x nvar of phases flags for parameters
	dvector gmax;
	gmax.allocate(1,nvar);
	gmax.initialize();
	string filename = "sa_res_all.txt";
	ifstream littxt(filename.c_str());
	if (littxt){

		cout << endl << "Reading sensitivity analysis data from the file: " << filename.c_str() << endl;
		//skip one line
		string s;
		getline(littxt,s,'\n');
		int id;
		double func,G, parval;
		string parname;
		string line;
		while (getline(littxt, line)) {
			istringstream ss(line);
			ss >> id >> parname >> parval;
			pnames[id-1] = parname;
			//ss.ignore(1024,'\t'); to skip 
			for (int ivar=1; ivar<=25; ivar++){
				ss >> G;
				if (gmax(ivar)<sqrt(G*G)) gmax(ivar) = sqrt(G*G);
			}
			ss >> func;		
		}
	}
	else {
		cerr << "The file sa_res_all.txt was not found! Will exit now..." << endl;
		exit(1);
	}
	littxt.close();


	ivector phase_par_flags;
	phase_par_flags.allocate(1,nvar);
	phase_par_flags.initialize();
	sc.param->set_all_false(pnames);

	cout << endl << "Starting optimization in phases... "<< endl;
	for (int phase_no=0; phase_no<4; phase_no++){
	
		cout << endl << "Entering phase "<< phase_no <<  endl;

		const int crit = pow(10,4-phase_no);
		std::ostringstream ostr;
		ostr << phase_no;

		for (int ivar=1; ivar<=nvar; ivar++){
			if (gmax[ivar]>=crit){
				//cout<< gmax[ivar] << " " << pnames[ivar-1] << endl;
				//TEMPORALLY will require the same initpafile being used 
				//in the sensitivity analysis. 
				//TODO: use the pnames in order to 
				//find the index of parfile parameter to be chosen 
				//as a control variable in the current phase
				//It DOES NOT work at the moment since parfile_names are 
				//initializes in xinit routine!!!
				/*
				int idx = 1;
       				for (int idx=1; idx<=nvar; idx++){
		                        string pname = sc.param->parfile_names[idx-1];
					//pname.erase(0,1);
					cout << pname << " " << pnames[ivar-1] << endl;
					int res = strcmp(pname.c_str(),pnames[ivar-1].c_str());
        	                	if (!res) break;
				}
				phase_par_flags[idx] = 1;
				*/
				phase_par_flags[ivar] = 1;
			}
		}
		int nvar_phase = sc.param->set_var_parameters(phase_par_flags, pnames);
		
		string dirout = get_path(parfile);
		string initparfile = dirout +"/initparfile_" + ostr.str() + ".xml";
		sc.write(initparfile.c_str());

		independent_variables x(1, nvar_phase);
		x_names.allocate(1,nvar_phase);	
		sc.xinit(x, x_names);

		//gradient vector allocation and initialization
		dvector g(1, nvar); 

		g.initialize();
		cout << "Total number of variables in phase "<< phase_no << " is: " << nvar_phase << '\n'<<'\n';

		//function minimizer class
		fmm fmc(nvar_phase);

		//flags for function minimizer
		fmc.iprint = 1;
		fmc.crit = crit;
		fmc.imax = 30;
		fmc.scroll_flag = 1;
		fmc.ifn = 0;
		fmc.maxfn = sc.get_maxfn();
		if (fmc.maxfn <= 0)
		fmc.ireturn = -1; 

		double likelihood = 0;
		double elapsed_time = 0;

		int idx = 0;
		int itr = 0;
		ios::sync_with_stdio();

		//initialization of simulation parameters
		sc.OnRunFirstStep();

		time(&time_sec);
		time_t time1 = time_sec;
		//'run_flux' runs the model in forward (simulation) mode
		//'gradcalc' runs backward (adjoint) mode
		//'save_statistics' stores current information on function minimization
		string tempparfile = dirout +"/tempparfile_" + ostr.str() + ".xml";
		string newparfile  = dirout +"/newparfile_" + ostr.str() + ".xml";
		cout << "\nentering minimization loop" << endl;
		while (fmc.ireturn >= 0) {
			//call function minimizer
			fmc.fmin(likelihood, x, g);

			//update the statistics.out file if the solution has been improved
			int itn = fmc.itn;
			if ((itn==itr+1 && (idx>=1))){
				time(&time_sec);
				time_t time2 = time_sec;
				elapsed_time += (double)(time2-time1)/60.0;
				time1 = time2;
				sc.save_statistics(dirout,x_names,likelihood,g,elapsed_time,idx-1,itr,nvar_phase);
				sc.write(tempparfile.c_str());
				itr = itn;
			}
			//reset control parameters and run the model and its adjoint
			if (fmc.ireturn > 0) {
				likelihood = sc.run_flux((dvar_vector)x);
				gradcalc(nvar_phase,g); 
				cout << "function evaluation " << idx++ << endl;
			}
		}
		string file1 = "optim.rep"; 
		string file2 = "optim_" + ostr.str() + ".rep"; 
		string CMD = "cp " + file1 + " " + file2; 
		int out = system(CMD.c_str());
		if (out == -1) cerr << "Could not copy the file " << file1 << " to " << file2 << endl; 

		sc.write(newparfile.c_str());
		remove(tempparfile.c_str());
	}

	//after minimization is finished one simulation will 
	//be run with estimated parameters; outputs will be saved
//	gradient_structure::set_NO_DERIVATIVES();
//	sc.run_flux((dvar_vector)x, true);
//	sc.write(newparfile.c_str());

	//writes new parameters on the screen
	sc.param->outp_param(x_names,nvar);

	time(&time_sec);
	cout << "\nfinished time: " << ctime(&time_sec) << endl;

	time_t time2 = time_sec;
	double total_elapsed_time = (double)(time2-time0) / 60.0;
	cout << "\ntotal time: " << total_elapsed_time << " minutes" << endl;
	delete [] pnames;
	return 0;
}


///1. Option for computing likelihood projection in 2D parametric space.
void Hyperspace_projection(SeapodymCoupled& sc, dvar_vector x)
{
//	sc.param->set_gradcalc(false);
	const int Npars = sc.param->nb_varproj-1;
	ivector ix(0,Npars); ix.initialize(); 
	int n1 = sc.param->varproj_nsteps[0];
	int n2 = sc.param->varproj_nsteps[1];
	int nmax = max(n1,n2)-1;  
	dmatrix xvalues(0,Npars,0,nmax); xvalues.initialize();
	dmatrix pars(0,Npars,0,nmax); pars.initialize();
	dmatrix Lproj(0,n1-1,0,n2-1); 
	Lproj.initialize();

	sc.param->get_param_index(ix, xvalues, pars);

	sc.OnRunFirstStep();
	clock_t time1 = clock();
	cout << "\nstarting hyperspace projection computation for ";
	for (int n=0; n<sc.param->nb_varproj; n++) cout << sc.param->varproj[n] << " ";
	cout << endl;

	ofstream ofs;
	const char* filename = "hyperproj.out";
	ofs.open(filename, ios::out);
	for (int n=0; n<sc.param->nb_varproj; n++)
		ofs << sc.param->varproj[n] << " ";
	ofs << "\n" << n1 << " " << n2 << "\n"; 
	for (int n=0; n<=Npars; n++){
		for (int i=0; i<sc.param->varproj_nsteps(n); i++)
			ofs << pars(n,i)<< " ";
		ofs << "\n";
	}
	ofs.close();
		
	for (int i=0; i<n1; i++){
		for (int j=0; j<n2; j++){
			cout << j+i*n2+1<< ": ";
			for (int n=0; n<=Npars; n++){
				int k;
				if (n==0)  k = i; 
				if (n==1)  k = j;
				x[ix(n)] = xvalues(n,k); cout << pars(n,k) << " "; 
			}
			Lproj(i,j) = sc.run_flux(x);
		}
		ofs.open(filename, ios::app);
		for (int j=0; j<n2; j++)
			ofs << Lproj(i,j) << " ";
		ofs << "\n";
		ofs.close();
	}

	time_t time2 = clock();
	double total_elapsed_time = (double)((time2-time1)/CLOCKS_PER_SEC)/60.0;
	cout << "\ntotal time: " << total_elapsed_time << " minutes" << endl;
	//cleanup_temporary_files();
}

///3. Option for parametric sensitivity analysis
void Sensitivity_analysis(const char* parfile, const int sftype)
{

	SeapodymCoupled sc(parfile);
	sc.param->set_scalc(true);
	
	//values to solve for
	const int nvar = sc.nvarcalc();
	independent_variables x(1, nvar);
	adstring_array x_names(1,nvar);
	sc.xinit(x, x_names);
	cout << "Total number of variables: " << nvar << '\n'<<'\n';

	sc.OnRunFirstStep();

	dvector s(1, nvar); s.initialize();

	clock_t time1 = clock();

	if (sftype==0){
		sc.param->like_types[0] = 7;
		//cout << "\nstarting computing sensitivities using model predictions only" << endl;

		cout << "\nstarting computing sensitivities using likelihood and gradient comp" << endl;

		dvector g(1, nvar); g.initialize();

		double func_predict = sc.run_flux((dvar_vector)x);
		gradcalc(nvar,g); 
		s = g/func_predict;

	}
	else if (sftype==1){
		cout << "\nstarting computing sensitivities using likelihood" << endl;

		const double eps = 1e-3;
		gradient_structure::set_NO_DERIVATIVES();
		cout << "At current parameters..." << endl;
		double like_opt = sc.run_flux((dvar_vector)x);
		cout << "At the boundaries..." << endl;
		for (int i=1; i<=nvar; i++){
			double xs = x(i);
			x(i) = sc.param->par_init_lo(i,eps);
			double like = sc.run_flux((dvar_vector)x);
			double s_lo = like/like_opt - 1.0;

			x(i) = sc.param->par_init_up(i,eps);
			like = sc.run_flux((dvar_vector)x);
			double s_up = like/like_opt - 1.0;

			x(i) = xs;
			s(i) = max(abs(s_lo),abs(s_up));
			cout << i << " " << x_names[i] << " " << s(i) << endl;
		}
	}
	else if (sftype==2){//ONE-AT-a-TIME sensitivity analysis
		cout << "\nstarting computing likelihoods for OAT sensitivity analysis" << endl;
		string dirout = get_path(parfile);
		string newparfile  = dirout +"/newparfile.xml";

		gradient_structure::set_NO_DERIVATIVES();
		dvector xr;
		int nxr = 25;//used 50 for first experiments with tags, it's too much
		xr.allocate(0,nxr);
		xr.initialize();
		double like = 1e2*sc.param->get_parval(1);//to be used for the seed
		double fmin = 1e10; //uncomment to always have non-increasing (<=) fmin between iterations

		for (int i=1; i<=nvar; i++){
		//for (int i=nvar; i>=1; i--){
		//for (int i=1; i<=2; i++){
		
			//double fmin = 1e10; //uncomment to allow increase in function value in iterations
			
			double xmin = x(i);//parameter value at start
			int n=(int)like;
			random_number_generator r(n);
			randu(r);
			xr.fill_randu(r);

			for (int k=0; k<nxr; k++){
			//for (int k=0; k<1; k++){
				x(i) = sc.param->par_init_step(i,xr[k]);

				like = sc.run_flux((dvar_vector)x);
				if (fmin>like){
					fmin = like;
					xmin = x(i);
				}
				cout << i << "." << k+1 << " \t" << x_names[i] << " \t" << sc.param->get_parval(i) << " " << like << endl;
			}
			x(i) = xmin; //if fmin not improved, xmin contains value at start 	
		}
		//Note, in case if xmin wasn't updated in the last iteration, 
		//the instruction 'x(i)=xmin' is useless, then need to reset 
		//parameters, i.e. to pass them to the VarParam class:
		sc.param->reset(x);
		//parameters corresponding to the minimal function value:
		sc.param->outp_param(x_names,nvar);
		cout << "Minimal function value: " << fmin << endl; 
		//write parfile with "best" parameters:
		sc.param->total_like = fmin;
		sc.write(newparfile.c_str());
	}
	else if (sftype==3){//just a forward run, usually to be used in ALL-AT-a-TIME sensitivity analysis

		gradient_structure::set_NO_DERIVATIVES();
		cout << "\nComputing likelihood only: " << endl << endl;
		double like = sc.run_flux((dvar_vector)x);
		cout << like << endl;	

	}	
	if (sftype<2){
		cout << endl << "N \t" << "parameter \t" << "\trelative sensitivity" << endl;
		cout << "-------------------------------------------------------" << endl;
	
		for (int i=1; i<=nvar; i++){
			int l = length(x_names[i]);
			string tab = "\t";
			if (l<16) tab += "\t";
			if (l<7) tab += "\t";
			cout << i <<  " \t" << x_names[i] << tab << s[i]<< endl;	
		}
	}
	
	time_t time2 = clock();
	double total_elapsed_time = (double)((time2-time1)/CLOCKS_PER_SEC)/60.0;
	cout << "\ntotal time: " << total_elapsed_time << " minutes" << endl;
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

string get_path(const char* parfile)
{
  size_t pos;
  string full_path = string(parfile);
  pos = full_path.find_last_of("/");
  if (pos==string::npos) 
	return ".";
  string dir = full_path.substr(0,pos);
  return dir;
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

