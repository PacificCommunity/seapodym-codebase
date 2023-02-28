#include <iostream>
#include <cstring>
#include <cstdlib>
  #include <pthread.h>
  #include <fvar.hpp>
  #include <adthread.h>

using std::cout;
void help(char* argv0);
int OptionToCode(char* Option,int &sub_option);
void seapodym_threaded(void *threadid);
void init_outfile();

#define NUM_THREADS 2

int main(int argc, char** argv) {

	int cmp_regime = -1;
	int sub_option = 0;
	int k=1;
	char *cmdLineOption = argv[k];
	cmp_regime = OptionToCode(cmdLineOption, sub_option);
	if (cmp_regime==-1) k=0;
	if (cmp_regime==-3) help(argv[0]);
	if (argc < k+2) {cout << "Too few parameters... \n"; help(argv[0]);}

	init_outfile();

	int nthread = NUM_THREADS;
	int ngroups = 1;
	ivector ng(1,ngroups);
	ng(1) = nthread;
	ad_comm::pthread_manager = new adpthread_manager(ngroups,ng,500);

	new_thread_data* data1 = new new_thread_data[nthread+1];
	for (int i=1;i<=nthread;i++) {
		data1[i].thread_no = i; // only the thread number is important
		data1[i].m=0;           // not used
	}
   
	ad_comm::pthread_manager->attach_code(&seapodym_threaded);
	ad_comm::pthread_manager->create_all(data1);

	for (int kk=1;kk<=nthread;kk++){
		int ipar = kk+50;

		// take control of the constant buffer for sending
		ad_comm::pthread_manager->cwrite_lock_buffer(kk);
		ad_comm::pthread_manager->send_int(ipar,kk);
		// release the constant buffer
		ad_comm::pthread_manager->cwrite_unlock_buffer(kk);
	}
		//delete ad_comm::pthread_manager;
		//ad_comm::pthread_manager=0;
		//pthread_exit(NULL);	
//	}
	pthread_exit(NULL);	
/* does not work for the moment
      for(t=0; t<nthread; t++){
      cout << "In main: creating thread " <<  t << "\n";
      //rc = pthread_create(&threads[t], NULL, seapodym_threaded, (void *)t);
ad_comm::pthread_manager->write_lock_buffer(t);
     // signal the thread that we are not finished
     ad_comm::pthread_manager->send_int(1,t); 
     // send current values of parameters
     ad_comm::pthread_manager->send_dvariable(a,t); 
     ad_comm::pthread_manager->send_dvariable(b,t); 
     // release the variable buffer
     ad_comm::pthread_manager->write_unlock_buffer(t);      
   }
*/
  //dmatrix s;
  //s.allocate(0,nthread,1,nvar); s.initialize();
      
//  dvector g;
    /*
  for (int kk=1;kk<=nthread;kk++)
  {
      // take control of the variable buffer for reading
      ad_comm::pthread_manager->read_lock_buffer(kk);
      //g = ad_comm::pthread_manager->get_dvector(kk);
      //dvector g = ad_comm::pthread_manager->get_dvector(kk);
     // release the variable buffer
      ad_comm::pthread_manager->read_unlock_buffer(kk);
      //cout << "got the gradient from thread "<< kk << ": " << g << endl;
  }   
  */

 //   delete ad_comm::pthread_manager;
 //   ad_comm::pthread_manager=0;
   /* Last thing that main() should do */
//	return 0;
}

void init_outfile(){

	ofstream wtxt;
	char* outfile = "./sa_res.txt";
	wtxt.open(outfile, ios::out);
			
	if (wtxt){
		wtxt << "id" << "\t" << "name" << "\t" << "value" << "\t" << "G" << "\t" << "F" << endl;
	}
	wtxt.close();
}

int OptionToCode(char* op, int &sub_option) {
	const int N = 22;
	char *cmdop[N] = {"-s","-p","-H","-sa","-t","-h","-v","--simulation","--likelihood-projection","--hessian","--sensitivity-analysis","--twin-experiment","--help","--version","-sa=0","-t=0","--sensitivity-analysis=0","--twin-experiment=0","-sa=1","-t=1","--sensitivity-analysis=1","--twin-experiment=1"};
	int cmpCode[N] = {0,1,2,3,4,-3,-2,0,1,2,3,4,-3,-2,3,4,3,4,3,4,3,4};
	for (int i=0; i<N; i++)
		if (strcmp(op,cmdop[i])==0){
			if (i>17) sub_option = 1;
			if (cmpCode[i]==-2) {
				cout << "SEAPODYM with parameter estimation 2.0 \n";
				cout << "Copyright (C) 2008, Patrick Lehodey and Inna Senina\n";
				exit(0);
			}
			return cmpCode[i];
		}

	return -1;
}

void help(char* argv0) {
	cout << "Usage:" << argv0 << " [options] parfile \n";
	cout << "Options: \n";
	cout << "  -h, --help \t\t\t Print this message and exit.\n";
	cout << "  -H, --hessian \t\t Compute Hessian matrix.\n";
	cout << "  -p, --projection \t\t Compute 2D-projection of the likelihood on a grid specified in parfile.\n";
	cout << "  -s, --simulation \t\t Run simulation without optimization.\n";
	cout << "  -sa[=FLAG]  \t\t\t Make sensitivity analysis. By default[=0] sensitivity function takes model predictions only.\n";

	cout << "   --sensitivity-analysis[=FLAG] If FLAG=1 the sensitivity function takes both predictions and observations.\n";
	cout << "  -t[=FLAG]   \t\t\t Perform identical (by default, or FLAG=0) twin experiment.\n";
	cout << "   --twin-experiment=[FLAG] \t If FLAG=1 the noise will be added to the artificial data.\n";
	cout << "  -v, --version \t\t Print version number and exit.\n";
	exit(0);
}



/*
int main(int argc, char** argv) {
	int cmp_regime = 4;
	if (argc > 2) cmp_regime = atoi(argv[2]);

	if ((argc > 3)||(argc < 2)) {
		cerr << "Usage["__FILE__ << ':' << __LINE__ << "]: " << argv[0] << " parfile\n";
		return EXIT_SUCCESS;
	}
	return seapodym_coupled(argv[1],cmp_regime);
	//return seapodym_coupled("newparfile.xml",0);
}
*/
