#include "calpop.h"

///Forward functions for: 
///precaldia and caldia functions, which perform discrete approximation
///of advection-diffusion equation and precomputes diagonal coefficients.
	
const double Vmax_diff  = 1.25;//1.9;
//const double rmax = 0.0;// no reduction of currents due to vertical movement here!
const double rmax = 0.99;//0.99;
const double rc = 0.0005;//0.002;//0.003;
const double rho = 0.99;

void CCalpop::precaldia_comp(const PMap& map, CParam& param, CMatrices& mat, const dmatrix& habitat, const dmatrix& total_pop, double MSS, double MSS_size_slope, double sigma_species, double c_diff_fish, const int sp, const int age, const int jday)
{
	mat.diffusion_x.initialize();
	mat.diffusion_y.initialize();
	mat.advection_x.initialize();
	mat.advection_y.initialize();

	CBord bord;
	const double length  = param.length[sp][age]*0.01; //convert to meters;
	const int    agemax  = param.sp_nb_cohorts[sp]-1;
	const double lmax    = param.length[sp][agemax]*0.01;
	const double dx	     = deltax;
	const double dy	     = deltay;
	const double dt	     = param.deltaT;
	//Inna, Nov 2014:
	//Theoretical V_mss = MSS_BL* pow(L,MSS_BL_slope) 
	//power constant suggested by Weihs in 1977 and 1981 is MSS_BL_slope = 0.43
	//Ex. parameters MSS=2.0 and b=0.30103 give (10,3.25,2,1.23) BL for 
	//BL=(0.1,0.5,1.0,2.0) m. correspondingly which results in speeds (1,1.6,2,2.5) m/sec. 
	const double CHI_x   = MSS*pow(length,MSS_size_slope)*(3600*24.0*dt/1852)*dx;
	const double CHI_y   = MSS*pow(length,MSS_size_slope)*(3600*24.0*dt/1852)*dy;
	const double Dspeed  = Vmax_diff-0.25*length/lmax;//fixed, given in 'body length' units
	const double Dinf    = pow(Dspeed*length*3600*24.0*dt/1852,2)/(4.0*dt);
	const double Dmax    = sigma_species*Dinf;
cout << age << " " << Dmax << endl;
	double rho_x = 0.0;
	double rho_y = 0.0;
	double v_x = 0.0;
	double v_y = 0.0;
	const double nb_layer = param.nb_layer;
	int nti = param.get_nbi();

	for (int i = map.imin; i <= map.imax; i++){	
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		dvector& diffusion_x = mat.diffusion_x[i];
		dvector& diffusion_y = mat.diffusion_y[i];
		dvector& advection_x = mat.advection_x[i];
		dvector& advection_y = mat.advection_y[i];
		dvector& speed	     = mat.speed[i];
		for (int j = jmin; j <= jmax; j++){
			const int nl = map.carte[i][j];
			if (nl>0) {
				bord.b   = map.bord_cell[i][j];
				int pos_x = bord.cotex();
				int pos_y = bord.cotey();
				double dHdx = 0;
				double dHdy = 0;
if (habitat(i,j)==0){cout << "Zero habitat at " << age << " " << i << " "<< j << " " << endl; exit(1);}
				if (nl>=nb_layer){// only if all layers exist! Make sure the mask has value 3 for all layers!

					if (pos_x == SANS){
                    				dHdx = (habitat[i+1][j] - habitat[i-1][j])/(2*dx);
					}
					else if (pos_x == D_FERME)
						dHdx = (habitat[i][j] - habitat[i-1][j])/(dx);
					else if (pos_x == G_FERME)
						dHdx = (habitat[i+1][j] - habitat[i][j])/(dx);
					
					if (map.global){
						if ((i == 1)&&(map.carte[nti-2][j])){
							if (pos_x == SANS)
								dHdx = (habitat[i+1][j] - habitat[nti-2][j])/(2*dx);
							else if (pos_x == D_FERME)
								dHdx = (habitat[i][j] - habitat[nti-2][j])/(dx);
						}
						if ((i == nti-2)&&(map.carte[1][j])){
							if (pos_x == SANS)
								dHdx = (habitat[1][j] - habitat[i-1][j])/(2*dx);
							else if (pos_x == G_FERME)
								dHdx = (habitat[1][j] - habitat[i][j])/(dx);
						}
					}
				
					if (pos_y == SANS){
						dHdy = (habitat[i][j+1] - habitat[i][j-1])/(2*dy);
					}
					else if (pos_y == D_FERME)
						dHdy = (habitat[i][j] - habitat[i][j-1])/(dy);
					else if (pos_y == G_FERME)
						dHdy = (habitat[i][j+1] - habitat[i][j])/(dy);

				}
				double diff_habitat = 1.0 - c_diff_fish*pow(habitat(i,j),3);

				double D = Dmax * diff_habitat;
				rho_x = 1.0- rho * sqrt(dHdx*dHdx) * dx;
				rho_y = 1.0- rho * sqrt(dHdy*dHdy) * dy;

				//Reduction of diffusion during seasonal migrations, i.e.
				//when spawning migrations take place
				double sfunc = mat.season_switch(sp,jday,j);
				D = (0.9*D*sfunc + D*(1.0-sfunc));

				double U = mat.u(i,j); 
				double V = mat.v(i,j);
				double C = sqrt(U*U+V*V);

				//reduction of upper layer current speed due to vertical behaviour 
				//of tuna within layers(used for test)
				//(1-r) the coefficient of reduction for ~0 currents and largest fish.
				const double r = rmax/(1.0+(rc*30.0/dt)*C);
				const double c = 1.0-r*length/lmax;

				//species velocity is proportional to the gradient of habitat index:
				v_x = CHI_x * dHdx;
				v_y = CHI_y * dHdy;
				//limit maximal velocity by Vinf to avoid approximation errors with 
				//finite differences in case of strong gradients
				if (v_x<0) v_x = -Vinf*param.func_limit_one(-v_x/Vinf);
				if (v_x>=0) v_x = Vinf*param.func_limit_one(v_x/Vinf);
				if (v_y<0) v_y = -Vinf*param.func_limit_one(-v_y/Vinf);
				if (v_y>=0) v_y = Vinf*param.func_limit_one(v_y/Vinf);
				advection_x[j] = c*U + v_x*mat.lat_correction[j];
				advection_y[j] = c*V + v_y;

				if (param.vert_movement[sp]){
					//correction of rho by passive advection:
					double fV = 1.0 - C/(500.0*dt/30.0+C);
					rho_x = rho_x * fV;
					rho_y = rho_y * fV;
				}

				diffusion_x[j] = rho_x * D * mat.lat_correction[j];
				diffusion_y[j] = rho_y * D;

				//to save predicted speed for some age group (currently the last one)
				if (age == agemax)
					speed[j] = sqrt(v_x*v_x+v_y*v_y);

			}
		}
	}

}
