#include <fvar.hpp>
#include "SeapodymCoupled.h"

double run_model(SeapodymCoupled& sc, dvar_vector x, dvector& g, const int nvar);
double run_sim(SeapodymCoupled& sc, dvar_vector x);

///2. Option for computing Hessian matrix
void Hessian_comp(const char* parfile)
{
	SeapodymCoupled sc(parfile);

	//values to solve for
	const int nvar = sc.nvarcalc();
	independent_variables x(1, nvar);
	adstring_array x_names(1,nvar);
	sc.xinit(x, x_names);
	cout << "Total number of variables: " << nvar << '\n'<<'\n';

	double delta = 1e-6;
	double epsilon = 0.1;
	dvector g1(1, nvar); g1.initialize();
	dvector g2(1, nvar); g2.initialize();
	dvector H1(1, nvar); H1.initialize();
	dvector H2(1, nvar); H2.initialize();
	dmatrix H(1,nvar,1,nvar); H.initialize();

	sc.prerun_model();

	double likelihood = 0.0;

	clock_t time1 = clock();
	cout << "\nstarting computing Hessian" << endl;
	likelihood = run_model(sc,x,g1,nvar);
	
	cout << "Likelihood and Gradient for estimated vector: \n" << likelihood << "; " << g1 << endl;
	
	//two-point finite-difference approximation of the Hessian
	for (int ix=1; ix<=nvar; ix++){
		double xs = x(ix);

		x(ix) = xs + delta; 
		likelihood = run_model(sc,x,g2,nvar);
		
		H1 = (g2-g1)/delta; 
			
		g2.initialize();

		x(ix) = xs + epsilon*delta; //1. step correction
		likelihood = run_model(sc,x,g2,nvar);

		H2 = (g2-g1)/(epsilon*delta); // 1. step correction

		H(ix) = (H2-epsilon*H1)/(1-epsilon); // 1. step correction

		cout << ix << ".\t"<< H(ix) << endl;

		x(ix) = xs;
		g2.initialize();
	}
	//making Hessian symmetric as it's subject to truncation error
	for (int ix=1; ix<nvar; ix++)
		for (int kx=ix+1; kx<=nvar; kx++){
			double av = 0.5*(H(ix,kx)+H(kx,ix));
			H(ix,kx) = av; H(kx,ix) = av;
		}
	
	dmatrix Cov = inv(H);
	double determ = det(H);
	dvector evalues = eigenvalues(H);
	ofstream ofs;
	const char* filename = "Hessian.out";

	ofs.open(filename, ios::out);
	ofs << nvar << "\n"; 

	dvector pars = sc.param->get_parvals();
	ofs << "Parameter\t" << "Est. value\t" << "Gradient" << "\n";   
	for (int i=1; i<=nvar; i++)
		ofs << x_names[i] << "\t" << pars(i) << "\t" << g1(i) << "\n";
	ofs << "\n";

	ofs << determ << "\n"; 
	ofs << evalues << "\n\n"; 

	for (int i=1; i<=nvar; i++){
		for (int j=1; j<=nvar; j++)
			ofs << H(i,j) << " ";
		ofs << "\n";
	}
	ofs << "\n";
	for (int i=1; i<=nvar; i++){
		for (int j=1; j<=nvar; j++)
			ofs << Cov(i,j) << " ";
		ofs << "\n";
	}
	ofs << "\n";
	ofs.close();

	time_t time2 = clock();
	double total_elapsed_time = (double)((time2-time1)/CLOCKS_PER_SEC)/60.0;
	cout << "\ntotal time: " << total_elapsed_time << " minutes" << endl;

}

///4. Option for Taylor derivative test
void Taylor_derivative_test(const char* parfile)
{

	SeapodymCoupled sc(parfile);
	sc.param->set_scalc(true);
	
	//values to solve for
	const int nvar = sc.nvarcalc();
	independent_variables x(1, nvar);
	adstring_array x_names(1,nvar);
	sc.xinit(x, x_names);
	cout << "Total number of variable parameters: " << nvar << '\n'<<'\n';

	sc.prerun_model();	

	//analytical derivatives
	dvector adv(1, nvar); adv.initialize();

	//finite difference derivatives at step
	const int nbs = 15;
	dvector step(0,nbs);
	for (int n=0; n<nbs; n++)
		step(n) = pow(10,0.5*(n-(nbs+1)));//1e-8,..,1e-1

	dmatrix fdv;
	fdv.allocate(1,nvar,0,nbs-1); fdv.initialize();

	//relative errors
	dmatrix err;
	err.allocate(1,nvar,0,nbs-1); err.initialize();

	clock_t time1 = clock();

	cout << "\nEntering Taylor derivative test" << endl;
	cout << "\n1. Computing analytical derivative(s)" << endl;

	dvector g(1, nvar); g.initialize();

	double func = run_model(sc,x,g,nvar);
	//gradcalc(nvar,g); 
	adv = g;
	cout << "Done." << endl;

	cout << "\n2. Computing finite-difference(s) with different stepping" << endl << endl;
	cout << "--- Parameters labelled with '^' - trend warning, '*' - elevated error warning, '!' - failed Taylor test ---" << endl <<endl;
	cout << setw(4) <<  left << " # "
                     << setw(10) << "X value" << " "
                     << setw(10) << "Analytical" << " ";
	for (int n=0; n<nbs; n+=2)
		cout << left << "rel.err"<<setw(7) << step(n) <<" ";
	cout << endl;
	for (int n=0; n<nbs; n+=2)
		cout << "-----------------";
	cout << endl;

	//Write detailed results to the text file
	ofstream ofs;
	const char* filename = "Taylor.out";

	ofs.open(filename, ios::out);
	ofs << "Function value: " << func << "\n"; 
	ofs << "Parameter names: ";
	for (int i=1; i<=nvar; i++)
		ofs << x_names[i] << "; ";
	ofs << "\n";
	ofs << "Step sizes: " << step << "\n";

	ofs << "N;" << "X.Value;" << "Analytical;";
	for (int n=0; n<nbs; n++)
		ofs << "fd"<< n+1 << ";";
	for (int n=0; n<nbs; n++)
		ofs << "rel.err"<< n+1 << ";";
	ofs << "\n";
	ofs.close();

	
	gradient_structure::set_NO_DERIVATIVES();

	string xindex;
	dvector derr;
	derr.allocate(0,nbs-1);

	//Compute central finite difference
	for (int i=1; i<=nvar; i++){

		if (adv(i)==0) 
			cout << "WARNING: analytical derivative is zero, test will return absolute error for this variable..." <<endl;

		std::ostringstream ostr;
		ostr << i;
		xindex = ostr.str();
		bool adv0 = false;
		double xval = x[i];
		for (int n=0; n<nbs; n++){
			x[i] = xval + step(n);
			double func_up = run_sim(sc,x);//sc.run_coupled((dvar_vector)x);

			x[i] = xval - step(n);
			double func_lo = run_sim(sc,x);//sc.run_coupled((dvar_vector)x);

			x[i] = xval;
	
			fdv(i,n) = (func_up-func_lo)/(2.0*step(n));

			if (adv(i)!=0)
				err(i,n) = abs((fdv(i,n) - adv[i])/adv[i]);
			else {  //return absolute errors
				err(i,n) = abs(fdv(i,n) - adv[i]);
				adv0 = true;
			}
		}

		//Catch failed test (by minimal value and shape of relative error):
		double errmin = min(err(i));
		derr = err(i) - errmin;
		int ind_min = nbs+1;
		int n=0;
		while (ind_min>nbs) {
			if (derr(n)==0)
				ind_min = n;
			n++;
		}

		double absminerr = abs(fdv(i,ind_min) - adv[i]);
		if (errmin>=1e-3){

			if (absminerr>=1e-5 && absminerr<1e-3) {
				xindex += "*";
				//cout << "Minimal rel. error " << errmin
				//	<<" > 0.001! Corresponding absolute error = " 
				//	<< absminerr << endl;
			}

			if (absminerr>=1e-3) {
				xindex += "!";
				cout << "Minimal relative, " << errmin 
					<<", and absolute, "<< absminerr 
					<< ", errors are > 0.001!" << endl;
			}
		}

		for (int n=0; n<nbs-1; n++){
			
			/*if (!adv0 && derr(n)<=derr(n+1) && ind_min > n){ 
				xindex += "^";
				cout << "Left-hand relative error not decreasing!" << endl;
			}*/
			if (!adv0 && derr(n)>=derr(n+1) && ind_min < n){ 
				xindex += "^";
				//cout << "Right-hand relative error not increasing!" << endl;
			}
		}

		cout << setw(4)  << left << xindex
                     << setw(10) << xval << " "
                     << setw(10) << adv[i] << " ";
		for (int n=0; n<nbs; n+=2)
			cout << setw(14) << left << err(i,n) << " ";
		cout << endl;

		//Writing to the file
		ofs.open(filename, ios::app);
		ofs << i << ";" << xval << ";" << setprecision(15) << adv[i] << ";";
		for (int n=0; n<nbs; n++)
			ofs << fdv(i,n) << ";";
		for (int n=0; n<nbs; n++)
			ofs << err(i,n) << ";";
		ofs << "\n";
		ofs.close();
	}

	cout << "Taylor test done.\n";

	
	time_t time2 = clock();
	double total_elapsed_time = (double)((time2-time1)/CLOCKS_PER_SEC)/60.0;
	cout << "\ntotal time: " << total_elapsed_time << " minutes" << endl;
}

