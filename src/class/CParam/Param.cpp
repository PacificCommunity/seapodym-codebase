// Param.cpp: implementation of the CParam class.
//
//////////////////////////////////////////////////////////////////////
//#include "StdAfx.h"
#include "Param.h"
#include "Utilities.h"
#include "ReadWrite.h"
//#include "sepodym.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

const double pi = 3.14159265358979;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
CParam::CParam()
{
	
}
CParam::~CParam()
{
	delete_param(true);
}
void  CParam::delete_param(bool flag)
{
	if ( (flag) &&(nb_species != 0)&& (sp_nb_cohort_ad[0]!=0) && (nb_fishery))
	{
		name_sp_by_file.clear();
		Utilities::delete1d(list_fishery_name);
		Utilities::delete1d(parfile_names);
	}
	if ((flag) && (nb_region))
	{
		for (int i=0 ; i<nb_region ; i++)
			delete [] area[i];
		delete [] area;
	}
	if (nb_EEZ) Utilities::delete1d(EEZ_name);

	if (flag)
	{
		frg_name.clear();
		sp_name.clear();
		strfile_u.clear();
		strfile_v.clear();
		strfile_t.clear();
		strfile_F.clear();
		strfile_umc.clear();
		strfile_vmc.clear();
		strfile_tmc.clear();
	}
}
//////////////////////////////////////////////////////////////
void  CParam::init_param()
{
	deltaX=1;		// dx (constant)
	deltaY=1;		// dy (constant)
	deltaT=1;		// dt (constant)
	iterationNumber=1;	// nbre d'iteration dans calcul derivees
	longitudeMin=0;		// longitude (0-359) minimale de la grille
	longitudeMax=0;		// longitude (0-359) maximale de la grille
	latitudeMin=0;		// latitude minimale (S negatif)  de la grille
	latitudeMax=0;		// latitude maximale  de la grille
	E = 0.0f;
	//max_NPP = 0.0f;
	save_first_yr=1950;
	inv_lambda_max = 0;
	inv_lambda_curv = 0.0f;
	Tr_max=0;		//Time (in days) before recruitment in the forage population
	Tr_exp=0;
	lambda=0;		//Forage mortality through time transfer
	//tstep_forage = 1;	//time step for the forage (days)
	sigma_fcte = 0.0f;
	//name_index = 0;
	nb_fishery =0;
	nb_region=0;
	nb_EEZ=0;
	m_f = 0.f;		//multiplicateur d'effort
	
}

void CParam::init_param_dym()
{
	//Generate the name of the forage files 
	for (int n=0; n< nb_forage; n++)	{
		//time series of prod		
		string dymFileS=str_dir_forage;
		dymFileS+="Fprod_";
		dymFileS+=frg_name[n];
		dymFileS+=".dym";
		strfile_S.push_back(dymFileS);

		//climatological series of prod
		string dymFileSmc=str_dir_forage;
		dymFileSmc+="mclFprod_";
		dymFileSmc+=frg_name[n];
		dymFileSmc+=".dym";
		strfile_Smc.push_back(dymFileSmc);

		//time series of forage biomass			
		string dymFileF=str_dir_forage;
		dymFileF+="Fbiom_";
		dymFileF+=frg_name[n];
		dymFileF+=".dym";
		strfile_F.push_back(dymFileF);

		//climatological series of forage biomass
		string dymFileFmc=str_dir_forage;
		dymFileFmc+="mclFbiom_";
		dymFileFmc+=frg_name[n];
		dymFileFmc+=".dym";
		strfile_Fmc.push_back(dymFileFmc);
	}
}

/////////////////////////////////////////////////////////
// calcul du facteur correctif fonction de la latitude
// la longueur d'un degre de longitude selon la latitude est proportionnelle au cosinus de la latitude
double CParam::correction_lat(double lat)
{
//	const double pi = 3.14159265358979;

	double g = (lat * pi) / 180.0;// transform lat(deg) in lat(radian)
	
	// For a given latitude, return the ratio to the distance of 1 deg at equator
	double function = 1.0 / cos(g);
	
	return function;
}

double CParam::cell_surface_area(int j)
{
//	const double pi = 3.14159265358979;
	double lat = (double) (latitudeMax - (j*deltaY/60.0) + deltaY/120.0);
	double R = 6378.1;//6371.0; //Radius of the Earth in km;
	double Phi1 = lat*pi/180.0;
	double Phi2 = (lat+deltaY/60.0)*pi/180.0;
	double dx_radian = (deltaX/60.0)*pi/180;

	double S = R*R*dx_radian*(sin(Phi2)-sin(Phi1));

	return S;
}

/////////////////////////////////////////////////////////
// calcul de la latitude a l'indice j
double CParam::lastlat(int j)
{
	double lastlat = (double) (latitudeMax - (j*deltaY/60.0) + deltaY/60.0);
	return lastlat;
}

//IMPORTANT: in order to get correct indices, (latitudeMax,longitudeMin) in parfile
//should give the NORTH-WEST corner of the grid cell and not the center, where
//all variables are calculated
double CParam::jtolat(int j)
{//returns latitudinal coordinate for Seapodym j index
	//latmax is the center: double lat = (double) (latitudeMax - (deltaY/60.0) * (j - 1));
	//latmax is the corner
	double lat = (double) (latitudeMax - (deltaY/60.0) * (j-0.5));
	return lat;
}

double CParam::itolon(int i)
{//returns longitudinal coordinate for Seapodym i index
	//lonmin is the center: double lon = (double) (longitudeMin + (deltaX/60.0) * (i - 1));
	//lonmin is the corner:
	double lon = (double) (longitudeMin + (deltaX/60.0) * (i-0.5));
	return lon;
}

int CParam::lattoj(double lat)
{//returns Seapodym index for the given latitude
	//if latmax is the center of the cell, i.e j = (latmax+.5*dy-lat)/dy + 1; 
	//if latmax is the NORTH corner of the first cell
	int j = (int) ((latitudeMax - lat)*60.0/deltaY + 1.0);
	return j;
}

int CParam::lontoi(double lon)
{//returns Seapodym index for the given longitude
	//lonmin is the center of the cell, i.e.i=(lon-(lonmin-.5*dx))/dx + 1;
	//if (lon<0) lon = lon+360;
	//if  lonmin is the WEST corner of the first cell
	int i = (int) ((lon - longitudeMin)*60.0/deltaX + 1.0);
//if (i==92) {cout << lon << endl; exit(1);}
	return i;
}

//lontoi<-function(lon) return(as.integer((lon - lonmin)*1/dx + 1.0))
//lattoj<-function(lat) return(as.integer((latmax - lat)*1/dy + 1.0))
//itolon<-function(i) return(lonmin + dx * (i-0.5))
//jtolat<-function(j) return(latmax - dy * (j-0.5))


/////////////////////////////////////////////////
//This function is the continuous version of
//if (x<=1) y = x; if (x>1) y = 1;
//uses: 1) catch equation (see dv_predicted_catch) 
//and 2) velocities (see dv_caldia)
//It uses south-opening rotated hyperbola with center 
//at (1,1) and 135 degree angle between asymptotes 

double CParam::func_limit_one(const double x)
{
	//parameters of hyperbola
	double phi = 22.5*pi/180.0;
	double a = 0.07;
	double e = 1.0/cos(phi);
	double b = a*sqrt(e*e-1.0);

	//coordinate center
	//shift is to have all y>=0
	double x0 = 1.0-0.00101482322788;
	double y0 = 1.0;

	//equation for hyperbola
	double sinsq = sin(phi)*sin(phi);
	double cossq = 1.0-sinsq;
	double rasq  = 1.0/(a*a);
	double rbsq  = 1.0/(b*b);
	double A = sinsq*rasq - cossq*rbsq;
	double B = -2.0*(x-x0)*cos(phi)*sin(phi)*(rasq+rbsq);
	double C = 1.0-(x-x0)*(x-x0)*(sinsq*rbsq-cossq*rasq);

	return(y0+(B+sqrt(B*B-4.0*A*C))/(2*A));
}

double CParam::dffunc_limit_one(const double x, const double dfy)
{
	//parameters of hyperbola
	double phi = 22.5*pi/180.0;
	double a = 0.07;
	double e = 1.0/cos(phi);
	double b = a*sqrt(e*e-1.0);

	//coordinate center
	//shift is to have all y>=0
	double x0 = 1.0-0.00101482322788;
	//double y0 = 1.0;

	//precompute 
	double sinsq = sin(phi)*sin(phi);
	double cossq = 1.0-sinsq;
	double rasq  = 1.0/(a*a);
	double rbsq  = 1.0/(b*b);
	double A = sinsq*rasq - cossq*rbsq;
	double B = -2.0*(x-x0)*cos(phi)*sin(phi)*(rasq+rbsq);
	double C = 1.0-(x-x0)*(x-x0)*(sinsq*rbsq-cossq*rasq);
	double D = sqrt(B*B-4.0*A*C);

	//derivatives
	double dfx = 0.0;
	//double y = y0+(B+sqrt(B*B-4.0*A*C))/(2*A);
	double dfC = -dfy/D;
	double dfB = dfy*(1.0+B/D)/(2*A);

	//double C = 1.0-(x-x0)*(x-x0)*(sinsq*rbsq-cossq*rasq);
	dfx -= 2.0*(x-x0)*(sinsq*rbsq-cossq*rasq)*dfC;

	//double B = -2.0*(x-x0)*cos(phi)*sin(phi)*(rasq+rbsq);
	dfx -= 2.0*cos(phi)*sin(phi)*(rasq+rbsq)*dfB;

	return(dfx);
}

/*
double CParam::func_limit_one(const double m)
{
        const double r = 0.1;
        double phi = 135.0*pi/180.0;
        double alpha = pi-phi;
        double eps = r*tan(0.5*alpha);
        double xinf = 1-eps/sqrt(2);
        double xsup = 1+eps;
        if (m<xinf)
                return(m);
        if (m>=xinf && m<=xsup)
                return(1-r+sqrt(r*r-(m-xsup)*(m-xsup)));
        if (m>xsup)
                return(1);
}

double CParam::dffunc_limit_one(const double m)
{
        const double r = 0.1;
        double phi = 135.0*pi/180.0;
        double alpha = pi-phi;
        double eps = r*tan(0.5*alpha);
        double xinf = 1-eps/sqrt(2);
        double xsup = 1+eps;
	if (m>xsup)
		return(0.0);
	if (m>=xinf && m<=xsup) 
		return(-(m-xsup)/sqrt(r*r-(m-xsup)*(m-xsup)));
	if (m<xinf)
		return(1.0);	
}
*/

void CParam::afcoef(const double lon, const double lat, dmatrix& a, int& ki, int& kj, const int f)
{// returns matrix of coefficients for the biomass aggregation over fishing reso*reso area
 // used in predicted_catch computation:
 //  |a(ki,kj) ...       | 
 //  | ...                   |
 //  | ...  a(ki+i,kj+j) ... |
 //  | ...               ... |
 //  | ...       a(ki+m,kj+n)|, 

 // calculated as the fishing area inside a given grid cell
 
	a.initialize(); 
	//coordinates of fishing area
	const float reso = fishery_reso(f);
	const float xw = lon-.5*reso;
	const float yn = lat+.5*reso;
	const float xe = lon+.5*reso;
	const float ys = lat-.5*reso;
	const int iw = lontoi(xw); const int ie = lontoi(xe);
	const int jn = lattoj(yn); const int js = lattoj(ys);

	ki = iw; kj = jn;
	double y0 = yn;
	double dx=0; double dy=0;
	for (int jj=jn; jj<=js; jj++){
		double yg = jtolat(jj)-.5*deltaY/60.0;

		dy = y0-yg;
		if (jj==js) dy = y0-ys;
		double x0 = xw;
		for (int ii=iw; ii<=ie; ii++){

			double xg = itolon(ii)+.5*deltaX/60.0;
			dx = xg-x0;
			if (ii==ie) dx = xe-x0;
			a(ii-iw,jj-jn) = dx*dy;		
			x0 = xg;
		}
		y0 = yg;
	}
}


double CParam::selectivity_comp(const int sp, const int age, const int f, const int k/*, const int step_count*/)
{
	double selectivity = 0;
	double sigma_sq;
	double selected_length = s_length_sp_fishery(sp,k);


    	switch (s_func_type[f]) {
                case 1://logistic
                        selectivity = length(sp,age)/(s_slope_sp_fishery(sp,k) + length(sp,age));
                        break;

		case 2://sigmoid
			selectivity = 1/(1+exp(s_slope_sp_fishery(sp,k)*(selected_length-length(sp,age))));
			break;
		
		case 3://asymmetric Gaussian
			sigma_sq = pow(s_slope_sp_fishery(sp,k),2);
			if (length(sp,age)<=selected_length)
				selectivity = exp(-pow(length(sp,age)-selected_length,2)/sigma_sq);
			else
				selectivity = (1-s_asympt_sp_fishery(sp,k))*exp(-pow(length(sp,age)-selected_length,2)/sigma_sq)+s_asympt_sp_fishery(sp,k);

			break;

		case 4://two-modal Gaussian
			const double s_length1 = selected_length;
			const double s_length2 = s_asympt_sp_fishery(sp,k);//temporal
			sigma_sq = pow(s_slope_sp_fishery(sp,k),2);//the same sigma - temporal
			const double g1 = 0.15*exp(-pow(length(sp,age)-s_length1,2)/sigma_sq);
			const double g2 = exp(-pow(length(sp,age)-s_length2,2)/sigma_sq);
			selectivity = (g1+g2)/(1+exp(-pow(s_length2-s_length1,2)/sigma_sq)+exp(-pow(s_length1-s_length2,2)/sigma_sq));
			break;

	}
	return selectivity;
}

void CParam::rbin_input2d(string file_in, const imatrix& carte, DMATRIX& mat2d, int nbi, int nbj, int nbytetoskip)
{
	//CReadWrite::rdym_input2d(file_in, carte, mat2d, nbi, nbj, nbytetoskip);
	ifstream litbin(file_in.c_str(), ios::binary | ios::in);
	
	if (!litbin)
	{
		cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: Unable to read file \"" << file_in << "\"\n";
		exit(1);
	}

	litbin.seekg(nbytetoskip, ios::cur);

 	const int sizeofDymInputType = sizeof(float);
	float buf;
	for (int j=1;j<nbj-1;j++){
		for (int i=1;i<nbi-1;i++){
			litbin.read(( char *)&buf,sizeofDymInputType);
			if (carte(i,j)){
				mat2d[i][j]= buf;
				//if (mat2d[i][j]<=-99)//indicates invalid value
				if (mat2d[i][j]<=-9)//tmp: ecco
					mat2d[i][j] = 0;
			}
		}
	}
	litbin.close();

}

void CParam::rbin_input2d(string file_in, DMATRIX& mat2d, int nbi, int nbj, int nbytetoskip)
{
	ifstream litbin(file_in.c_str(), ios::binary | ios::in);
	
	if (!litbin)
	{
		cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: Unable to read file \"" << file_in << "\"\n";
		exit(1);
	}

	litbin.seekg(nbytetoskip, ios::cur);

 	const int sizeofDymInputType = sizeof(float);
	float buf;
	for (int j=1;j<nbj-1;j++){
		for (int i=1;i<nbi-1;i++){
			litbin.read(( char *)&buf,sizeofDymInputType);
			mat2d[i][j]= buf;
		}
	}

	litbin.close();
}


void CParam::rbin_mat2d(string file_in, const imatrix& carte, DMATRIX& mat2d, int nlat, int nlong, int nbytetoskip)
{
	ifstream litbin(file_in.c_str(), ios::binary | ios::in);
	
	if (!litbin)
	{
		cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: Unable to read file \"" << file_in << "\"\n";
		exit(1);
	}

	litbin.seekg(nbytetoskip, ios::cur);

 	const int sizeofTmpDataType = sizeof(double);
	double buf;

	for (int i=0;i<nlat;i++)
	{
		for (int j=0;j<nlong;j++)
		{
			litbin.read(( char *)&buf,sizeofTmpDataType);
			if (carte(i,j))
				mat2d[i][j]= buf;
		}
	}

	litbin.close();
}

void CParam::rbin_mat2d(string file_in, DMATRIX& mat2d, int nlat, int nlong, int nbytetoskip)
{
	ifstream litbin(file_in.c_str(), ios::binary | ios::in);
	
	if (!litbin)
	{
		cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: Unable to read file \"" << file_in << "\"\n";
		exit(1);
	}

	litbin.seekg(nbytetoskip, ios::cur);

 	const int sizeofTmpDataType = sizeof(double);
	double buf;

	for (int i=0;i<nlat;i++)
	{
		for (int j=0;j<nlong;j++)
		{
			litbin.read(( char *)&buf,sizeofTmpDataType);
			mat2d[i][j]= buf;
		}
	}

	litbin.close();
}


