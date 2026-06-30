#include <fvar.hpp>
#include "SeapodymCoupled.h"

string get_path(const char* full_path);
double run_model(SeapodymCoupled& sc, dvar_vector x, dvector& g, const int nvar);
double run_sim(SeapodymCoupled& sc, dvar_vector x);

///1. Option for computing Hessian matrix
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
	
	//bound-aware one-sided FD Hessian with Richardson extrapolation.
	//scaled x in [-1,1] (ADMB arcsin transform): step INWARD -- forward by default,
	//backward when a forward step would cross the upper bound (|x|=1).
	for (int ix=1; ix<=nvar; ix++){
		double xs = x(ix);
		double s  = (xs + delta > 1.0) ? -1.0 : 1.0;   //-1 => backward (near upper bound)

		x(ix) = xs + s*delta;
		likelihood = run_model(sc,x,g2,nvar);
		H1 = (g2-g1)/(s*delta);
		g2.initialize();

		x(ix) = xs + s*epsilon*delta; //1. step correction
		likelihood = run_model(sc,x,g2,nvar);
		H2 = (g2-g1)/(s*epsilon*delta); //1. step correction

		H(ix) = (H2-epsilon*H1)/(1-epsilon); //1. step correction

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

	double max_abs_eig = fabs(evalues(1)), min_abs_eig = fabs(evalues(1));
	int    n_negative  = (evalues(1) < 0.0) ? 1 : 0;
	for (int i = 2; i <= nvar; i++){
		double a = fabs(evalues(i));
		if (a > max_abs_eig) max_abs_eig = a;
		if (a < min_abs_eig) min_abs_eig = a;
		if (evalues(i) < 0.0) n_negative++;
	}
	double condition_number = max_abs_eig / min_abs_eig;   // largest / smallest eigenvalue by magnitude (spectral condition number; = lambda_max/lambda_min when PD Hessian)


	dvector pars = sc.param->get_parvals();

	//=== convergence / identifiability diagnostics ========================
	// Convergence is NOT judged by the raw gradient norm max|g_i|: it is not
	// invariant under reparametrisation, so a fixed threshold confounds proximity
	// to the minimum with the scaling of the parameters and of the objective.

	double gmax = 0.0;
	for (int i=1; i<=nvar; i++){ double a = fabs(g1(i)); if (a>gmax) gmax = a; }

	// Newton decrement lambda2 = g'H^{-1}g. Half of it is the objective decrease
	// predicted on stepping to the local quadratic minimum; 0.5*lambda2/L is the
	// scale-invariant relative improvement still available (<<1 => at the minimum).
	dvector Hinv_g = Cov*g1;
	double lambda2 = g1*Hinv_g;
	dvector nstep  = -Hinv_g;                 // Newton step
	double maxrel_step = 0.0;
	for (int i=1; i<=nvar; i++)
		if (pars(i)!=0.0){ double a = fabs(nstep(i)/pars(i)); if (a>maxrel_step) maxrel_step = a; }

	// Eigenvalues of H: all positive => positive definite => genuine local minimum.
	// Condition number = max/min eigenvalue flags ill-conditioning (near-flat dirs).
	double emin = min(evalues), emax = max(evalues);
	int is_pd = (emin > 0.0);

	// Standard errors = sqrt(diag(Cov)); NaN/huge here => non-PD or unidentified.
	dvector SE(1,nvar);
	for (int i=1; i<=nvar; i++) SE(i) = sqrt(Cov(i,i));

	// Correlation matrix (built once; scale-free; reused below).
	dmatrix Corr(1,nvar,1,nvar);
	for (int i=1; i<=nvar; i++)
		for (int j=1; j<=nvar; j++)
			Corr(i,j) = Cov(i,j)/(SE(i)*SE(j));

	// (1) FLATTEST direction: dominant eigenvector of Cov (= smallest-eigenvalue
	//     eigenvector of H), the direction the likelihood curves least in absolute
	//     terms. Flags a SATURATED / individually unidentified parameter sitting in
	//     a flat region of its functional form (small gradient AND small curvature).
	//     var_flat = variance along it = 1/min eigenvalue of H.
	dvector vsat(1,nvar); vsat = 1.0/sqrt((double)nvar);
	for (int it=0; it<500; it++){ vsat = Cov*vsat; double nv = norm(vsat); if (nv>0.0) vsat /= nv; }
	double var_flat = vsat*(Cov*vsat);

	// (2) MOST-COLLINEAR direction: dominant eigenvector of the correlation matrix
	//     (standardized, scale-free). Flags a TRADE-OFF -- standardized parameters
	//     the data constrain only in combination, not individually. vinfl =
	//     variance-inflation factor along it. The two directions are complementary:
	//     (1) finds saturation, (2) finds confounding.
	dvector vcol(1,nvar); vcol = 1.0/sqrt((double)nvar);
	for (int it=0; it<500; it++){ vcol = Corr*vcol; double nv = norm(vcol); if (nv>0.0) vcol /= nv; }
	double vinfl = vcol*(Corr*vcol);

	cout << "\n--- convergence / identifiability diagnostics ---" << endl;
	cout << "L = " << likelihood << " ; Gmax = " << gmax << " (scale-dependent; not used for convergence)" << endl;
	cout << "Newton decrement^2 = " << lambda2
	     << " ; relative remaining = " << 0.5*lambda2/likelihood << " (<<1 => at the minimum)" << endl;
	cout << "max relative Newton step |dx/x| = " << maxrel_step << endl;
	cout << (is_pd ? "PD (local min)" : "NOT PD -> saddle/flat")
	     << " ; condition number = " << condition_number << endl;

	ofstream ofs;
	const char* filename = "Hessian.out";
	ofs.open(filename, ios::out);
	ofs << nvar << "\n\n";

	ofs << "Parameter\tEstimate\tGradient\tStdErr\tCV\n";
	ofs << "# StdErr = sqrt(diag(inverse Hessian)); CV = StdErr/|Estimate| (relative uncertainty)\n";
	for (int i=1; i<=nvar; i++)
		ofs << x_names[i] << "\t" << pars(i) << "\t" << g1(i) << "\t"
		    << SE(i) << "\t" << SE(i)/fabs(pars(i)) << "\n";
	ofs << "\n";

	ofs << "CONVERGENCE DIAGNOSTICS\n";
	ofs << "Note, negative eigenvalues at a converged minimum may indicate FD noise (|min_eig| large, scales ~1/h)\n";
	ofs << "Run Hessian with different FD steps and verify it's a minimum IF: \n";
	ofs << "- Newton_decrement_sq, relative_remaining and variance_along are stable across FD steps\n";
	ofs << "- relative_remaining is small (<<1)\n";
	ofs << "If so, the non-PD is attributable to the ill-conditioning (condition_number > 1e6), not to a real saddle\n";
	and positive, if relative_remaining is small variance_along and condition_number (if >10^6 then ill-conditioning)\n"
	ofs << "likelihood\t"           << likelihood             << "\t# objective (neg. log-likelihood) \n";
	ofs << "Gmax\t"                 << gmax                   << "\t# max|gradient|, scale-dependent - NOT a reliable convergence test\n";
	ofs << "Newton_decrement_sq\t"  << lambda2                << "\t# g'H^-1g; curvature-weighted distance to the minimum\n";
	ofs << "pred_remaining_dL\t"    << 0.5*lambda2            << "\t# objective decrease predicted to reach the quadratic minimum\n";
	ofs << "relative_remaining\t"   << 0.5*lambda2/likelihood << "\t# pred_remaining_dL / L; scale-invariant; <<1 => at the minimum\n";
	ofs << "max_rel_newton_step\t"  << maxrel_step            << "\t# largest |dx/x|; how far parameters still want to move\n";
	ofs << "determinant\t"          << determ                 << "\t# det(H); >0 consistent with positive definite\n";
	ofs << "min_eigenvalue\t"       << emin                   << "\t# smallest curvature (flattest direction)\n";
	ofs << "max_eigenvalue\t"       << emax                   << "\t# largest curvature (stiffest direction)\n";
	ofs << "nb_neg_eigenvalues\t"   << n_negative             << "\t# count of negative eigenvalues; 0 => PD";
	ofs << "positive_definite\t"    << (is_pd ? "yes" : "no") << "\t# yes => (local) minimum; no => saddle / not a minimum\n";
	ofs << "condition_number\t"     << condition_number       << "\t# |max|/|min| eigenvalue magnitude (spectral); valid whether PD or not; high => ill-conditioned\n\n";
	ofs << "variance_along\t" 	<< var_flat               << "\t# = 1 / smallest-magnitude eigenvalue = variance along the flattest direction; stable & positive => real minimum\n";

	ofs << "FLATTEST direction (dominant eigenvector of covariance = smallest-curvature direction of the likelihood)\n";
	ofs << "# Shows a SATURATED / individually unidentified parameter: a flat region of its functional form.\n";
	{
		ivector idx(1,nvar);
		for (int i=1;i<=nvar;i++) idx(i)=i;
		for (int a=1;a<nvar;a++){ int best=a;
			for (int b=a+1;b<=nvar;b++) if (fabs(vsat(idx(b)))>fabs(vsat(idx(best)))) best=b;
			int tmp=idx(a); idx(a)=idx(best); idx(best)=tmp; }
		bool marked=false;
		for (int a=1;a<=nvar;a++){ int i=idx(a);
			ofs << x_names[i] << "\t" << vsat(i);
			if (!marked && fabs(vsat(i))<0.01){ ofs << "\t< 0.01 onward"; marked=true; }
			ofs << "\n"; }
	}
	ofs << "\n";

	ofs << "MOST-COLLINEAR direction (dominant eigenvector of correlation matrix; standardized units)\n";
	ofs << "# Shows a TRADE-OFF: standardized parameters constrained only in combination, not individually.\n";
	ofs << "variance_inflation\t" << vinfl << "\t# joint variance / uncorrelated-direction variance\n";
	{
		ivector idx(1,nvar);
		for (int i=1;i<=nvar;i++) idx(i)=i;
		for (int a=1;a<nvar;a++){ int best=a;
			for (int b=a+1;b<=nvar;b++) if (fabs(vcol(idx(b)))>fabs(vcol(idx(best)))) best=b;
			int tmp=idx(a); idx(a)=idx(best); idx(best)=tmp; }
		bool marked=false;
		for (int a=1;a<=nvar;a++){ int i=idx(a);
			ofs << x_names[i] << "\t" << vcol(i);
			if (!marked && fabs(vcol(i))<0.01){ ofs << "\t< 0.01 onward"; marked=true; }
			ofs << "\n"; }
	}
	ofs << "\n";

	// --- strongly cross-correlated parameter pairs (|rho| > 0.8) ---
	ofs << "Cross-correlated pairs (|correlation| > 0.8)\n";
	ofs << "# |r|>0.8 (rho^2>0.64, >64% shared variance); * = |r|>0.9, ** = |r|>0.95 (effectively non-separable)\n";
	for (int i=1; i<=nvar; i++)
		for (int j=i+1; j<=nvar; j++)
			if (fabs(Corr(i,j)) > 0.8){
				const char* mark = (fabs(Corr(i,j)) > 0.95) ? "**" : (fabs(Corr(i,j)) > 0.9) ? "*" : "";
				char rbuf[16];
				snprintf(rbuf, sizeof rbuf, "%.2f", Corr(i,j));
				ofs << x_names[i] << "\t" << x_names[j] << "\t" << rbuf << mark << "\n";
			}
	ofs << "\n";

	ofs << "Eigenvalues:\n" << evalues << "\n\n";

	ofs << "Hessian:\n";
	for (int i=1; i<=nvar; i++){
		for (int j=1; j<=nvar; j++)
			ofs << H(i,j) << " ";
		ofs << "\n";
	}
	ofs << "\n";

	ofs << "Covariance (inverse Hessian):\n";
	for (int i=1; i<=nvar; i++){
		for (int j=1; j<=nvar; j++)
			ofs << Cov(i,j) << " ";
		ofs << "\n";
	}
	ofs << "\n";

	ofs << "Correlation matrix:\n";
	for (int i=1; i<=nvar; i++){
		for (int j=1; j<=nvar; j++)
			ofs << Corr(i,j) << " ";
		ofs << "\n";
	}
	ofs << "\n";
	ofs.close();

	time_t time2 = clock();
	double total_elapsed_time = (double)((time2-time1)/CLOCKS_PER_SEC)/60.0;
	cout << "\ntotal time: " << total_elapsed_time << " minutes" << endl;

}

///2. Option for computing likelihood projection in 2D parametric space.
void Hyperspace_projection(SeapodymCoupled& sc, dvar_vector x)
{
	if (sc.param->get_doc_empty_status("/hyperspace_projection")){
		cout << "\nDeclare two variables for computing hyperspace projection:" <<endl;
		cout << "<hyperspace_projection>" << endl;
		cout << "    <variables nb=\"2\"/>" << endl;
		cout << "    <var1 name=\"selected_parameter_name\" nsteps=\"number_of_steps\"/>" << endl;
		cout << "    <var2 name=\"selected_parameter_name\" nsteps=\"number_of_steps\"/>" << endl;
		cout << "</hyperspace_projection>" << endl;
		cout << "Exit now..." <<endl;
		exit(1);
	}

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
			Lproj(i,j) = run_sim(sc,x);
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

///3. Options for parametric sensitivity analyses: FLAG 0 - local sensitivity, 1 - edge sensitivity, 2 - OAT of global sensitivity analysis, 3 - AAT of global sensitivity analysis.
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

	sc.prerun_model();	

	dvector s(1, nvar); s.initialize();

	clock_t time1 = clock();

	if (sftype == 0){

		cout << "\nComputing local sensitivities using likelihood gradient" << endl;
		dvector g(1, nvar); g.initialize();

		
		double func_predict = run_model(sc,x,g,nvar);
		cout << "\nLikelihood at current parameters: " << func_predict << endl;

		dvector parderivative = sc.param->dpar_dx(x,nvar);
			
		s = elem_prod(parderivative, g) / func_predict;

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
	else if (sftype==1){//Edge sensitivity metric
		
		cout << "\nComputing likelihood change at parameters boundaries, L_at_boundary - L_cur" << endl;

		const double eps = 1e-3;
		gradient_structure::set_NO_DERIVATIVES();

		double like_cur = run_sim(sc,x);
		cout << "Likelihood at current parameters: " << like_cur << endl;
		cout << endl;

		int wname = 24;
		int wnum  = 14;
		cout << setw(2) << left << "N" << " " 
			<< setw(wname) << "parameter "
			<< right
			<< setw(wnum) << "lower dL" << " "
			<< setw(wnum) << "upper dL" << " " 
			<< setw(wnum) << "max(|dL|)/Lcur" << endl;
		int wline = wname + 3*wnum + 5;
		cout << string(wline, '-') << '\n';

		for (int i=1; i<=nvar; i++){
			double xs = x(i);
			x(i) = sc.param->par_init_lo(i,eps);
			double like = run_sim(sc,x);
			double s_lo = like - like_cur;

			x(i) = sc.param->par_init_up(i,eps);
			like = run_sim(sc,x);
			double s_up = like - like_cur;

			x(i) = xs;

			//get the relative sensitivity
			s(i) = max(abs(s_lo),abs(s_up))/like_cur;

			int l = length(x_names[i]);
			string tab = " ";
			for (int k=1; k<wname-l; k++) tab += " ";
			cout << setw(2) << left << i << " "
				<< x_names[i] << tab 
				<< right << setprecision(6)
				<< setw(wnum) << s_lo << " "
				<< setw(wnum) << s_up << " "
				<< setprecision(4)
				<< setw(wnum) << s(i) << endl;
		}	
		cout << string(wline, '-') << '\n';
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

				like = run_sim(sc,x);
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
		double like = run_sim(sc,x);//sc.run_coupled((dvar_vector)x);
		//cout << like << endl;	
		
		// Get likelihood components
		double clike = sc.get_clike();
		double lflike = sc.get_lflike();
		double stocklike = sc.get_stocklike();
		double taglike = sc.get_taglike();
		double larvaelike = sc.get_larvaelike();

		// Likelihood breakdown
		cout << "Total:	 " << like << endl;
		cout << "Catch:	 " << clike << endl;
		cout << "LF:	 " << lflike << endl;
		cout << "Stock:	 " << stocklike << endl;
		cout << "Tags:	 " << taglike << endl;
		cout << "Larvae: " << larvaelike << endl;

	}	
	
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

	//finite difference derivatives at steps 1e-10,..,1e-0
	int nbs = 11;
	do {
		cout << "Will compute derivatives with 11 steps. Press ENTER if you want to continue, otherwise enter number 1, 3, or 6 and press ENTER." << endl;
		if (cin.peek() != '\n')
    			cin >> nbs;
	} while (cin.get() != '\n');

	if (nbs != 1 && nbs != 3 && nbs != 6 && nbs != 11){
		cerr << "Need to enter correct number of steps. Exit now!" << endl;
		exit(1);
	}

	int k = 1;
	dvector step(0,nbs-1);
	if (nbs == 3) k = 4;
	if (nbs == 6) k = 2;
	for (int n=0; n<nbs; n++)
		step(n) = pow(10,n*k - 10);
	
	if (nbs == 1) step(0) = 1e-6;
	int err_skip = 0;
	if (nbs == 11) err_skip = 1;

	dmatrix fdv;
	fdv.allocate(1,nvar,0,nbs-1); fdv.initialize();

	//relative errors
	dmatrix err;
	err.allocate(1,nvar,0,nbs-1); err.initialize();

	clock_t time1 = clock();

	cout << "\nEntering Taylor derivative test" << endl;
	cout << "\n1. Computing analytical derivative(s) and finite differences with " << nbs << " step(s)" << endl;

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
	for (int n=0; n<nbs; n+=(err_skip+1))
		cout << left << "rel.err"<<setw(7) << step(n) <<" ";
	cout << endl;
	for (int n=0; n<nbs; n+=(err_skip+1))
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
		for (int n=0; n<nbs; n+=(err_skip+1))
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

