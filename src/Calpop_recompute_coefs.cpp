#include "calpop.h"

void f_accessibility(dvector& l_access, dvector& lf_access, const dvector forage, const dvector O2, const dvector T, double twosigsq, double temp_age, double oxy_teta, double oxy_cr, const int nl, const int nb_forage, const ivector day_layer, const ivector night_layer, const double DL);
void f_accessibility(dvector& l_access, dvector& lf_access, const dvector forage, const dvector O2, const dvector T, double temp_mean, 
		double temp_max, double delta1, double delta2, double delta3, double oxy_teta, double oxy_cr, const int nl, 
		const int nb_forage, const ivector day_layer, const ivector night_layer, const double DL);


//constants section
const double Vmax_diff = 1.25;//1.9;
//const double rmax = 0.0;
const double rmax = 0.99;
const double rc = 0.0005;//0.002;//0.003;
const double rho = 0.99;


void CCalpop::RecompDiagCoef_juv(const PMap& map, CMatrices& mat, const int t_count, const dmatrix mortality, dmatrix& aa, dmatrix& bbm, dmatrix& cc, dmatrix& dd, dmatrix& ee, dmatrix& ff)
{
	dvector lat_correction(map.jmin,map.jmax);
	lat_correction.initialize();
	dmatrix u(map.imin, map.imax, map.jinf, map.jsup);
	dmatrix v(map.imin, map.imax, map.jinf, map.jsup);
	dmatrix diffusion_x, diffusion_y, advection_x, advection_y;
	advection_x.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	advection_y.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	diffusion_x.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	diffusion_y.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	diffusion_x.initialize();
	diffusion_y.initialize();
	advection_x.initialize();
	advection_y.initialize();

	u = mat.un[t_count][0];
	v = mat.vn[t_count][0];
	lat_correction = mat.lat_correction;

	for (int i = map.imin; i <= map.imax; i++){	
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			if (map.carte(i,j)){
				diffusion_x[i][j] = sigma_fcte * lat_correction[j];
				diffusion_y[i][j] = sigma_fcte ;
		
				advection_x[i][j] = u(i,j);
				advection_y[i][j] = v(i,j);
			}
		}
	}


	CBord bord;	
	double dx = deltax;
	double dxdx    = dx*dx;
	double twodxdx = 2*dxdx;

	for (int j=map.jmin; j <= map.jmax; j++){
		for (int i=map.iinf[j] ; i<=map.isup[j] ; i++){                            
			bord.b		= map.bord_cell[i][j];
			char pos	= bord.cotex();
			double sigma	= diffusion_x[i][j];
			double sigmam	= diffusion_x[i-1][j];
			double sigmap	= diffusion_x[i+1][j];
			double uinf	= advection_x[i-1][j];
			double ucur	= advection_x[i][j];
			double usup	= advection_x[i+1][j];

			aa[j][i] = d1(pos, sigmam, sigma, uinf, twodxdx, dxdx, dx);
			bbm[j][i] = d2(pos, sigmam, sigma, sigmap, ucur, twodxdx, dxdx, dx, 0);
			cc[j][i] = d3(pos, sigma, sigmap, usup, twodxdx, dxdx, dx);
		}
	}
	double dy = deltay;
	double dydy	= dy*dy;
	double twodydy	= 2*dydy;

	for (int i=map.imin; i <= map.imax; i++){
		for (int j=map.jinf[i] ; j<=map.jsup[i] ; j++){
			bord.b		= map.bord_cell[i][j];
			char pos	= bord.cotey();
			double sigma	= diffusion_y[i][j];
			double sigmam	= diffusion_y[i][j-1];
			double sigmap	= diffusion_y[i][j+1];
			double vinf	= advection_y[i][j-1];
			double vcur 	= advection_y[i][j];
			double vsup	= advection_y[i][j+1];

			dd[i][j] = d1(pos, sigmam, sigma, vinf, twodydy, dydy, dy);
			ee[i][j] = d2(pos, sigmam, sigma, sigmap, vcur, twodydy, dydy, dy, 0);
			ff[i][j] = d3(pos, sigma, sigmap, vsup, twodydy, dydy, dy);
		}
	}

	for (int j=map.jmin; j <= map.jmax; j++){
		for (int i=map.iinf[j] ; i<=map.isup[j] ; i++){  
			bbm[j][i] += mortality[i][j];
		}
	}
}

void CCalpop::Recomp_abc_coef(const PMap& map, CMatrices& mat, const int t_count, const dmatrix& mortality, dmatrix& aa, dmatrix& bbm, dmatrix& cc)
{
	DVECTOR lat_correction(map.jmin,map.jmax);
	lat_correction.initialize();
	DMATRIX u(map.imin, map.imax, map.jinf, map.jsup);
	dmatrix diffusion_x, advection_x;
	advection_x.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	diffusion_x.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	diffusion_x.initialize();
	advection_x.initialize();

	u = mat.un[t_count][0];
	lat_correction = mat.lat_correction;

	//Precaldia part
	for (int i = map.imin; i <= map.imax; i++){	
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			if (map.carte(i,j)){
				diffusion_x[i][j] = sigma_fcte * lat_correction[j];
				advection_x[i][j] = u[i][j];
			}
		}
	}

	//Caldia part
	CBord bord;	
	double dx      = deltax;
	double dxdx    = dx*dx;
	double twodxdx = 2*dxdx;

	for (int j=map.jmin; j <= map.jmax; j++){
		for (int i=map.iinf[j] ; i<=map.isup[j] ; i++){                            
			bord.b		= map.bord_cell[i][j];
			char pos	= bord.cotex();
			double sigma	= diffusion_x[i][j];
			double sigmam	= diffusion_x[i-1][j];
			double sigmap	= diffusion_x[i+1][j];
			double uinf	= advection_x[i-1][j];
			double ucur	= advection_x[i][j];
			double usup	= advection_x[i+1][j];

			aa[j][i]  = d1(pos, sigmam, sigma, uinf, twodxdx, dxdx, dx);
			bbm[j][i] = d2(pos, sigmam, sigma, sigmap, ucur, twodxdx, dxdx, dx, 0);
			cc[j][i]  = d3(pos, sigma, sigmap, usup, twodxdx, dxdx, dx);
		}
	}
	//Precalrec part
	for (int j=map.jmin; j <= map.jmax; j++){
		for (int i=map.iinf[j] ; i<=map.isup[j] ; i++){  
			bbm[j][i] += mortality[i][j];
		}
	}
}

void CCalpop::Recomp_DEF_coef(const PMap& map, CParam& param, CMatrices& mat, const int t_count, const int jday, const dmatrix& habitat, dmatrix& dd, dmatrix& ee, dmatrix& ff, dmatrix& advection_x, dmatrix& advection_y, const int sp, const int age, const double MSS, const double c_diff_fish, const double sigma_species)
{
	dmatrix diffusion_y;
	diffusion_y.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	diffusion_y.initialize();

	//PRECALDIA SECTION
	const double MSS_size_slope = param.MSS_size_slope[sp];
	const double length  = param.length[sp][age]*0.01;
	const double lmax    = param.length[sp][param.sp_nb_cohorts[sp]-1]*0.01;
	const double dx	     = param.deltaX;
	const double dy	     = param.deltaY;
	const double dt	     = param.deltaT;
	const double CHI_x   = MSS*pow(length,MSS_size_slope)*(3600*24.0*dt/1852)*dx;
	const double CHI_y   = MSS*pow(length,MSS_size_slope)*(3600*24.0*dt/1852)*dy;
	const double Dspeed = Vmax_diff-0.25*length/lmax;
	const double Dinf   = pow(Dspeed*length*3600*24.0*dt/1852,2)/(4.0*dt);


	const double Dmax    = sigma_species*Dinf;

	const double nb_layer = param.nb_layer;

	CBord bord;	
	for (int i = map.imin; i <= map.imax; i++){	
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			const int nl = map.carte(i,j);
			if (nl>0){
				bord.b	   = map.bord_cell[i][j];
				//bord.b   = map.nbl_bord_cell[i][j];
				char pos_x = bord.cotex();
				char pos_y = bord.cotey();
				double dHdx = 0;
				double dHdy = 0;


				if (nl >=nb_layer){

				if (pos_x == SANS)
                    			dHdx = (habitat[i+1][j] - habitat[i-1][j])/(2*dx);
				else if (pos_x == D_FERME)
					dHdx = (habitat[i][j] - habitat[i-1][j])/(dx);
				else if (pos_x == G_FERME)
					dHdx = (habitat[i+1][j] - habitat[i][j])/(dx);
				
				if (pos_y == SANS)
					dHdy = (habitat[i][j+1] - habitat[i][j-1])/(2*dy);
				else if (pos_y == D_FERME)
					dHdy = (habitat[i][j] - habitat[i][j-1])/(dy);
				else if (pos_y == G_FERME)
					dHdy = (habitat[i][j+1] - habitat[i][j])/(dy);

				}

				//double diff_habitat = 1 - habitat(i,j)/(c_diff_fish + habitat(i,j));
				double diff_habitat = 1.0 - c_diff_fish*pow(habitat(i,j),3);
				double D = Dmax * diff_habitat;

				double sfunc = mat.season_switch(sp,jday,j);
				D = (0.9*D*sfunc + D*(1.0-sfunc));


				double rho_y = 0.0;
				rho_y = 1.0- rho * sqrt(dHdy*dHdy) * dy;

				double U = mat.un[t_count][0][i][j]; double V = mat.vn[t_count][0][i][j];

				const double r = rmax/(1+rc*sqrt(U*U+V*V));
				const double c = 1.0-r*length/lmax;

				double v_x = CHI_x * dHdx;
				double v_y = CHI_y * dHdy;

				//limit maximal velocity my Vinf
				if (v_x<0) v_x = -Vinf*param.func_limit_one(-v_x/Vinf);
				if (v_x>=0) v_x = Vinf*param.func_limit_one(v_x/Vinf);
				if (v_y<0) v_y = -Vinf*param.func_limit_one(-v_y/Vinf);
				if (v_y>=0) v_y = Vinf*param.func_limit_one(v_y/Vinf);


				advection_x(i,j) = c*U + v_x*mat.lat_correction[j];// for computing derivatives in dv_caldia
				advection_y(i,j) = c*V + v_y;
				diffusion_y(i,j) = rho_y * D;
			}
		}
	}
	//CALDIA SECTION
	double dydy	= dy*dy;
	double twodydy	= 2*dydy;

	for (int i=map.imin; i <= map.imax; i++){
		for (int j=map.jinf[i] ; j<=map.jsup[i] ; j++){
			bord.b		= map.bord_cell[i][j];
			char pos	= bord.cotey();
			double sigma	= diffusion_y[i][j];
			double sigmam	= diffusion_y[i][j-1];
			double sigmap	= diffusion_y[i][j+1];
			double vinf	= advection_y[i][j-1];
			double vcur 	= advection_y[i][j];
			double vsup	= advection_y[i][j+1];

			dd[i][j] = d1(pos, sigmam, sigma, vinf, twodydy, dydy, dy);
			ee[i][j] = d2(pos, sigmam, sigma, sigmap, vcur, twodydy, dydy, dy, 0);
			ff[i][j] = d3(pos, sigma, sigmap, vsup, twodydy, dydy, dy);
		}
	}
}

void CCalpop::Recomp_DEF_UV_coef(const PMap& map, CParam& param, CMatrices& mat, dmatrix& u, dmatrix& v, const dmatrix& habitat, dmatrix& dd, dmatrix& ee, dmatrix& ff, dmatrix& advection_x, dmatrix& advection_y, const int sp, const int age, const double MSS, const double c_diff_fish, const double sigma_species, const int jday)
{
	DVECTOR lat_correction(map.jmin,map.jmax);
	lat_correction.initialize();
        for (int j=map.jmin;j<=map.jmax;j++){
                double lastlat = param.lastlat(j);
                lat_correction[j] = param.correction_lat(lastlat);
        }

	dmatrix diffusion_y;
	diffusion_y.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	diffusion_y.initialize();


	//PRECALDIA SECTION
	const double MSS_size_slope = param.MSS_size_slope[sp];
	const double length  = param.length[sp][age]*0.01;
	const double lmax    = param.length[sp][param.sp_nb_cohorts[sp]-1]*0.01;
	const double dx	     = param.deltaX;
	const double dy	     = param.deltaY;
	const double dt	     = param.deltaT;
	const double CHI_x   = MSS*pow(length,MSS_size_slope)*(3600*24.0*dt/1852)*dx;
	const double CHI_y   = MSS*pow(length,MSS_size_slope)*(3600*24.0*dt/1852)*dy;
	const double Dspeed  = Vmax_diff-0.25*length/lmax;
	const double Dinf    = pow(Dspeed*length*3600*24.0*dt/1852,2)/(4.0*dt);

	const double Dmax  = sigma_species*Dinf;

	const double nb_layer = param.nb_layer;

	CBord bord;	
	for (int i = map.imin; i <= map.imax; i++){	
		
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			const int nl = map.carte(i,j);
			if (nl>0){
				bord.b	   = map.bord_cell[i][j];
				//bord.b   = map.nbl_bord_cell[i][j];
				char pos_x = bord.cotex();
				char pos_y = bord.cotey();
				double dHdx = 0;
				double dHdy = 0;


				if (nl >= nb_layer){
					if (pos_x == SANS)
        	            			dHdx = (habitat[i+1][j] - habitat[i-1][j])/(2*dx);
					else if (pos_x == D_FERME)
						dHdx = (habitat[i][j] - habitat[i-1][j])/(dx);
					else if (pos_x == G_FERME)
						dHdx = (habitat[i+1][j] - habitat[i][j])/(dx);
				
					if (pos_y == SANS)
						dHdy = (habitat[i][j+1] - habitat[i][j-1])/(2*dy);
					else if (pos_y == D_FERME)
						dHdy = (habitat[i][j] - habitat[i][j-1])/(dy);
					else if (pos_y == G_FERME)
						dHdy = (habitat[i][j+1] - habitat[i][j])/(dy);

				}

				//double diff_habitat = 1 - habitat(i,j)/(c_diff_fish + habitat(i,j));
				double diff_habitat = 1.0 - c_diff_fish*pow(habitat(i,j),3);
				double D = Dmax * diff_habitat;

				double sfunc = mat.season_switch(sp,jday,j);
				D = (0.9*D*sfunc + D*(1.0-sfunc));

				double rho_y = 0.0;
				rho_y = 1.0- rho * sqrt(dHdy*dHdy) * dy;

				double U = u(i,j); double V = v(i,j);
				double C = sqrt(U*U+V*V);
				const double r = rmax/(1+rc*C);
				const double c = 1.0-r*length/lmax;

				double v_x = CHI_x * dHdx;
				double v_y = CHI_y * dHdy;

				//limit maximal velocity my Vinf
				if (v_x<0) v_x = -Vinf*param.func_limit_one(-v_x/Vinf);
				if (v_x>=0) v_x = Vinf*param.func_limit_one(v_x/Vinf);
				if (v_y<0) v_y = -Vinf*param.func_limit_one(-v_y/Vinf);
				if (v_y>=0) v_y = Vinf*param.func_limit_one(v_y/Vinf);


				advection_x(i,j) = c*U + v_x*lat_correction[j];// for computing derivatives in dv_caldia
				advection_y(i,j) = c*V + v_y;

				double fV = 1.0 - C/(500.0*dt/30.0+C);
				rho_y = rho_y * fV;

				diffusion_y(i,j) = rho_y * D;
			}
		}
	}
	//CALDIA SECTION
	double dydy	= dy*dy;
	double twodydy	= 2*dydy;

	for (int i=map.imin; i <= map.imax; i++){
		for (int j=map.jinf[i] ; j<=map.jsup[i] ; j++){
			bord.b		= map.bord_cell[i][j];
			char pos	= bord.cotey();
			double sigma	= diffusion_y[i][j];
			double sigmam	= diffusion_y[i][j-1];
			double sigmap	= diffusion_y[i][j+1];
			double vinf	= advection_y[i][j-1];
			double vcur 	= advection_y[i][j];
			double vsup	= advection_y[i][j+1];

			dd[i][j] = d1(pos, sigmam, sigma, vinf, twodydy, dydy, dy);
			ee[i][j] = d2(pos, sigmam, sigma, sigmap, vcur, twodydy, dydy, dy, 0);
			ff[i][j] = d3(pos, sigma, sigmap, vsup, twodydy, dydy, dy);
		}
	}
}


void CCalpop::RecompDiagCoef_adult(const PMap& map, CParam& param, CMatrices& mat, const int t_count, const int jday, const dmatrix& mortality, const dmatrix& habitat, dmatrix& aa, dmatrix& bbm, dmatrix& cc, dmatrix& dd, dmatrix& ee, dmatrix& ff, const int sp, const int age, const double MSS, const double c_diff_fish, const double sigma_species)
{
	dmatrix diffusion_x, diffusion_y, advection_x, advection_y;
	advection_x.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	advection_y.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	diffusion_x.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	diffusion_y.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	diffusion_x.initialize();
	diffusion_y.initialize();
	advection_x.initialize();
	advection_y.initialize();

	//PRECALDIA SECTION
	const double MSS_size_slope = param.MSS_size_slope[sp];
	const double length  = param.length[sp][age]*0.01;
	const double lmax    = param.length[sp][param.sp_nb_cohorts[sp]-1]*0.01;
	const double dx	     = param.deltaX;
	const double dy	     = param.deltaY;
	const double dt	     = param.deltaT;
	const double CHI_x   = MSS*pow(length,MSS_size_slope)*(3600*24.0*dt/1852)*dx;
	const double CHI_y   = MSS*pow(length,MSS_size_slope)*(3600*24.0*dt/1852)*dy;
	const double Dspeed  = Vmax_diff-0.25*length/lmax;
	const double Dinf    = pow(Dspeed*length*3600*24.0*dt/1852,2)/(4.0*dt);
	const double Dmax    = sigma_species*Dinf;

	const double nb_layer = param.nb_layer;

	CBord bord;	

	for (int i = map.imin; i <= map.imax; i++){	
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			const int nl = map.carte(i,j);
			if (nl>0){
				bord.b	   = map.bord_cell[i][j];
				//bord.b   = map.nbl_bord_cell[i][j];
				char pos_x = bord.cotex();
				char pos_y = bord.cotey();
				double dHdx = 0;
				double dHdy = 0;

				if (nl >= nb_layer){

					if (pos_x == SANS)
        	            			dHdx = (habitat[i+1][j] - habitat[i-1][j])/(2*dx);
					else if (pos_x == D_FERME)
						dHdx = (habitat[i][j] - habitat[i-1][j])/(dx);
					else if (pos_x == G_FERME)
						dHdx = (habitat[i+1][j] - habitat[i][j])/(dx);
				
					if (pos_y == SANS)
						dHdy = (habitat[i][j+1] - habitat[i][j-1])/(2*dy);
					else if (pos_y == D_FERME)
						dHdy = (habitat[i][j] - habitat[i][j-1])/(dy);
					else if (pos_y == G_FERME)
						dHdy = (habitat[i][j+1] - habitat[i][j])/(dy);
					
				}

				//double diff_habitat = 1 - habitat(i,j)/(c_diff_fish + habitat(i,j));
				double diff_habitat = 1.0 - c_diff_fish*pow(habitat(i,j),3);
				double D = Dmax * diff_habitat;

				double sfunc = mat.season_switch(sp,jday,j);
				D = (0.9*D*sfunc + D*(1.0-sfunc));

				double rho_x = 0.0;
				double rho_y = 0.0;

				rho_x = 1.0- rho * sqrt(dHdx*dHdx) * dx;
				rho_y = 1.0- rho * sqrt(dHdy*dHdy) * dy;

				double U = mat.un[t_count][0][i][j]; double V = mat.vn[t_count][0][i][j];

				double C = sqrt(U*U+V*V);
				const double r = rmax/(1+rc*C);
				const double c = 1.0-r*length/lmax;


				double v_x = CHI_x * dHdx;
				double v_y = CHI_y * dHdy;

				//limit maximal velocity by Vinf
				if (v_x<0) v_x = -Vinf*param.func_limit_one(-v_x/Vinf);
				if (v_x>=0) v_x = Vinf*param.func_limit_one(v_x/Vinf);
				if (v_y<0) v_y = -Vinf*param.func_limit_one(-v_y/Vinf);
				if (v_y>=0) v_y = Vinf*param.func_limit_one(v_y/Vinf);

				advection_x(i,j) = c*U + v_x * mat.lat_correction[j];
				advection_y(i,j) = c*V + v_y;
				diffusion_x(i,j) = rho_x * D * mat.lat_correction[j];
				diffusion_y(i,j) = rho_y * D;
			}
		}
	}
	//CALDIA SECTION
	double dxdx    = dx*dx;
	double twodxdx = 2*dxdx;

	for (int j=map.jmin; j <= map.jmax; j++){
		for (int i=map.iinf[j] ; i<=map.isup[j] ; i++){                            
			bord.b		= map.bord_cell[i][j];
			char pos	= bord.cotex();
			double sigma	= diffusion_x[i][j];
			double sigmam	= diffusion_x[i-1][j];
			double sigmap	= diffusion_x[i+1][j];
			double uinf	= advection_x[i-1][j];
			double ucur	= advection_x[i][j];
			double usup	= advection_x[i+1][j];

			aa[j][i]  = d1(pos, sigmam, sigma, uinf, twodxdx, dxdx, dx);
			bbm[j][i] = d2(pos, sigmam, sigma, sigmap, ucur, twodxdx, dxdx, dx, 0);
			cc[j][i]  = d3(pos, sigma, sigmap, usup, twodxdx, dxdx, dx);
		}
	}
	double dydy	= dy*dy;
	double twodydy	= 2*dydy;

	for (int i=map.imin; i <= map.imax; i++){
		for (int j=map.jinf[i] ; j<=map.jsup[i] ; j++){
			bord.b		= map.bord_cell[i][j];
			char pos	= bord.cotey();
			double sigma	= diffusion_y[i][j];
			double sigmam	= diffusion_y[i][j-1];
			double sigmap	= diffusion_y[i][j+1];
			double vinf	= advection_y[i][j-1];
			double vcur 	= advection_y[i][j];
			double vsup	= advection_y[i][j+1];

			dd[i][j] = d1(pos, sigmam, sigma, vinf, twodydy, dydy, dy);
			ee[i][j] = d2(pos, sigmam, sigma, sigmap, vcur, twodydy, dydy, dy, 0);
			ff[i][j] = d3(pos, sigma, sigmap, vsup, twodydy, dydy, dy);
		}
	}

	//PRECALREC SECTION
	for (int j = bbm.rowmin(); j <= bbm.rowmax(); j++){
		for (int i = bbm[j].indexmin(); i <= bbm[j].indexmax(); i++){
			bbm[j][i] += mortality[i][j];
		}
	}
}

void CCalpop::RecompDiagCoef_UV_adult(const PMap& map, CParam& param, CMatrices& mat, const int t_count, const int jday, const dmatrix& mortality, const dmatrix& habitat, dmatrix& aa, dmatrix& bbm, dmatrix& cc, dmatrix& dd, dmatrix& ee, dmatrix& ff, const int sp, const int age, const double MSS, const double c_diff_fish, const double sigma_species)
{
	dmatrix diffusion_x, diffusion_y, advection_x, advection_y;
	advection_x.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	advection_y.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	diffusion_x.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	diffusion_y.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	diffusion_x.initialize();
	diffusion_y.initialize();
	advection_x.initialize();
	advection_y.initialize();

	//Parameters to recompute accessibility to layers 
	const double oxy_teta = param.a_oxy_habitat[sp];
	const double oxy_cr   = param.b_oxy_habitat[sp];
	const double sigma_ha = param.sigma_ha[sp][age];
	const double temp_age = param.temp_age[sp][age];
	const double twosigsq = 2.0*sigma_ha*sigma_ha;
	const double temp_max = param.b_sst_spawning(sp);
	const double delta1   = param.thermal_func_delta[0][sp];
	const double delta2   = param.thermal_func_delta[1][sp];
	const double delta3   = param.thermal_func_delta[2][sp];
	
	const int Tfunc_Gaussian = param.gaussian_thermal_function[sp];

	const int nb_forage = param.get_nbforage();
	ivector day_layer(0,nb_forage-1); day_layer = param.day_layer;
	ivector night_layer(0,nb_forage-1); night_layer = param.night_layer;
	
	const int nb_layer = param.nb_layer;
	dvector lf_access(0,nb_layer-1);
	dvector l_access(0,nb_layer-1);
	dvector F(0,nb_forage-1);
	dvector O2(0,nb_forage-1); 
	dvector T(0,nb_forage-1);  
	//end of accessiblity parameters section

	//PRECALDIA SECTION
	const double MSS_size_slope= param.MSS_size_slope[sp];
	const double length  = param.length[sp][age]*0.01;
	const double lmax    = param.length[sp][param.sp_nb_cohorts[sp]-1]*0.01;
	const double dx	     = param.deltaX;
	const double dy	     = param.deltaY;
	const double dt	     = param.deltaT;
	const double CHI_x   = MSS*pow(length,MSS_size_slope)*(3600*24.0*dt/1852)*dx;
	const double CHI_y   = MSS*pow(length,MSS_size_slope)*(3600*24.0*dt/1852)*dy;

	const double Dspeed  = Vmax_diff-0.25*length/lmax;
	const double Dinf    = pow(Dspeed*length*3600*24.0*dt/1852,2)/(4.0*dt);
	const double Dmax    = sigma_species*Dinf;

	CBord bord;	

	for (int i = map.imin; i <= map.imax; i++){	
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			const int nl = map.carte(i,j);
			if (nl>0){
				bord.b	   = map.bord_cell[i][j];
				//bord.b   = map.nbl_bord_cell[i][j];
				char pos_x = bord.cotex();
				char pos_y = bord.cotey();
				double dHdx = 0;
				double dHdy = 0;


				if (nl >= nb_layer){

					if (pos_x == SANS)
        	            			dHdx = (habitat[i+1][j] - habitat[i-1][j])/(2*dx);
					else if (pos_x == D_FERME)
						dHdx = (habitat[i][j] - habitat[i-1][j])/(dx);
					else if (pos_x == G_FERME)
						dHdx = (habitat[i+1][j] - habitat[i][j])/(dx);
				
					if (pos_y == SANS)
						dHdy = (habitat[i][j+1] - habitat[i][j-1])/(2*dy);
					else if (pos_y == D_FERME)
						dHdy = (habitat[i][j] - habitat[i][j-1])/(dy);
					else if (pos_y == G_FERME)
						dHdy = (habitat[i][j+1] - habitat[i][j])/(dy);

				}

				//double diff_habitat = 1 - habitat(i,j)/(c_diff_fish + habitat(i,j));
				double diff_habitat = 1.0 - c_diff_fish*pow(habitat(i,j),3);
				double D = Dmax * diff_habitat;

				double sfunc = mat.season_switch(sp,jday,j);
				D = (0.9*D*sfunc + D*(1.0-sfunc));

				double rho_x = 0.0;
				double rho_y = 0.0;

				rho_x = 1.0- rho * sqrt(dHdx*dHdx) * dx;
				rho_y = 1.0- rho * sqrt(dHdy*dHdy) * dy;

				//recalculate U and V
				l_access.initialize();
				lf_access.initialize();
				O2.initialize();
				T.initialize();
				for (int n=0; n<nb_forage; n++)
					F(n) = mat.forage(t_count,n,i,j);
				//for (int l=0; l<nl; l++){
				for (int l=0; l<nb_layer; l++){
					O2(l) = mat.oxygen(t_count,l,i,j);
					T(l) = mat.tempn(t_count,l,i,j);
				}		
				const double DL = mat.daylength(jday,j)/24.0;

				//need to recompute average currents, attn, accessibility is computed for 
				//ages param.age_compute_habitat[sp][age] only (passed here as age)
				if (Tfunc_Gaussian){
					f_accessibility(l_access,lf_access,F,O2,T,twosigsq,temp_age,oxy_teta,oxy_cr,
							nb_layer,nb_forage,day_layer,night_layer,DL);
				} else {
					f_accessibility(l_access,lf_access,F,O2,T,temp_age,temp_max,delta1,delta2,delta3,oxy_teta,oxy_cr,
							nb_layer,nb_forage,day_layer,night_layer,DL);
				}
						
				double U = 0.0; double V = 0.0;
				//for (int l=0; l<nl; l++){
				if (nl<=nb_layer){
					for (int l=0; l<nb_layer; l++){
						U += mat.un(t_count,l,i,j) * lf_access(l);
						V += mat.vn(t_count,l,i,j) * lf_access(l);
					}	
				}
				//end of recalculation section

				double C = sqrt(U*U+V*V);
				const double r = rmax/(1+rc*C);
				const double c = 1.0-r*length/lmax;

				double v_x = CHI_x * dHdx;
				double v_y = CHI_y * dHdy;

				//limit maximal velocity by Vinf
				if (v_x<0) v_x = -Vinf*param.func_limit_one(-v_x/Vinf);
				if (v_x>=0) v_x = Vinf*param.func_limit_one(v_x/Vinf);
				if (v_y<0) v_y = -Vinf*param.func_limit_one(-v_y/Vinf);
				if (v_y>=0) v_y = Vinf*param.func_limit_one(v_y/Vinf);

				advection_x(i,j) = c*U + v_x * mat.lat_correction[j];
				advection_y(i,j) = c*V + v_y;

				double fV = 1.0 - C/(500.0*dt/30.0+C);
				rho_x = rho_x * fV;
				rho_y = rho_y * fV;

				diffusion_x(i,j) = rho_x * D * mat.lat_correction[j];
				diffusion_y(i,j) = rho_y * D;

			}
		}
	}
//cout << "1: " << t_count << " " << age << " "<<  l_access << " " << lf_access << " " << temp_age << " " << temp_max << " " << delta1 << " " << delta2 << " " << delta3 << " " << oxy_teta << " " << oxy_cr << " " << nb_layer << endl;			
//cout << norm(advection_x) << endl;
	//CALDIA SECTION
	double dxdx    = dx*dx;
	double twodxdx = 2*dxdx;

	for (int j=map.jmin; j <= map.jmax; j++){
		for (int i=map.iinf[j] ; i<=map.isup[j] ; i++){                            
			bord.b		= map.bord_cell[i][j];
			char pos	= bord.cotex();
			double sigma	= diffusion_x[i][j];
			double sigmam	= diffusion_x[i-1][j];
			double sigmap	= diffusion_x[i+1][j];
			double uinf	= advection_x[i-1][j];
			double ucur	= advection_x[i][j];
			double usup	= advection_x[i+1][j];

			aa[j][i]  = d1(pos, sigmam, sigma, uinf, twodxdx, dxdx, dx);
			bbm[j][i] = d2(pos, sigmam, sigma, sigmap, ucur, twodxdx, dxdx, dx, 0);
			cc[j][i]  = d3(pos, sigma, sigmap, usup, twodxdx, dxdx, dx);
		}
	}
	double dydy	= dy*dy;
	double twodydy	= 2*dydy;

	for (int i=map.imin; i <= map.imax; i++){
		for (int j=map.jinf[i] ; j<=map.jsup[i] ; j++){
			bord.b		= map.bord_cell[i][j];
			char pos	= bord.cotey();
			double sigma	= diffusion_y[i][j];
			double sigmam	= diffusion_y[i][j-1];
			double sigmap	= diffusion_y[i][j+1];
			double vinf	= advection_y[i][j-1];
			double vcur 	= advection_y[i][j];
			double vsup	= advection_y[i][j+1];

			dd[i][j] = d1(pos, sigmam, sigma, vinf, twodydy, dydy, dy);
			ee[i][j] = d2(pos, sigmam, sigma, sigmap, vcur, twodydy, dydy, dy, 0);
			ff[i][j] = d3(pos, sigma, sigmap, vsup, twodydy, dydy, dy);
		}
	}

	//PRECALREC SECTION
	for (int j = bbm.rowmin(); j <= bbm.rowmax(); j++){
		for (int i = bbm[j].indexmin(); i <= bbm[j].indexmax(); i++){
			bbm[j][i] += mortality[i][j];
		}
	}
}

void CCalpop::RecompM_sp(const PMap& map, const CParam& param, dmatrix& M, const dmatrix& H, const double age, const int sp)
{
	double Mp_max   = param.Mp_mean_max[sp];
	double Mp_exp   = param.Mp_mean_exp[sp];
	double Ms_slope = param.Ms_mean_slope[sp];
	double Ms_max   = param.Ms_mean_max[sp];
	double range    = param.M_mean_range[sp]; 
	
	//const double Mp = (Ms_max+Mp_max) * exp(- Mp_exp * age);
	const double Mp = Mp_max * exp(- Mp_exp * age);
	const double Ms = Ms_max * pow(age,Ms_slope);

	double	Rage = 1.0/(age+2.5) + range; 
 	//this function gives +- 33%, 25% and 18%

	double Hval = 0.5;

	M = Mp + Ms;
	for (int i = map.imin; i <= map.imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			if (map.carte(i,j)){
				M(i,j) *= pow(1.0+Rage,1-H(i,j)/Hval); 
			}
		}
	}
}



