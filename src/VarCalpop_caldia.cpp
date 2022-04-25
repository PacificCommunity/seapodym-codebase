#include "calpop.h"

void CCalpop::caldia(const PMap& map, const VarParamCoupled& param, dvar_matrix& diffusion_x,  dvar_matrix& advection_x, dvar_matrix& diffusion_y, dvar_matrix& advection_y)
{
	dvarsA.initialize();	
	dvarsB.initialize();	
	dvarsBM.initialize();	
	dvarsC.initialize();	
	dvarsD.initialize();	
	dvarsE.initialize();	
	dvarsF.initialize();
	Ybet.initialize();

	CBord bord;
	const double dx = param.deltaX;
	const double dxdx = dx*dx;
	const double twodxdx = 2*dxdx;
	for (int j=map.jmin; j <= map.jmax; j++) {
		for (int i=map.iinf[j] ; i<=map.isup[j] ; i++) {                            
			bord.b			= map.bord_cell[i][j];
			const char pos		= bord.cotex();
			dvariable sigma		= diffusion_x[i][j];
			dvariable sigmam	= diffusion_x[i-1][j];
			dvariable sigmap	= diffusion_x[i+1][j];
			dvariable uinf		= advection_x[i-1][j];
			dvariable u		= advection_x[i][j];
			dvariable usup		= advection_x[i+1][j];

			dvarsA[j][i] = d1(pos, sigmam, sigma, uinf, twodxdx, dxdx, dx);
			dvarsB[j][i] = d2(pos, sigmam, sigma, sigmap, u, twodxdx, dxdx, dx);
			dvarsC[j][i] = d3(pos, sigma, sigmap, usup, twodxdx, dxdx, dx);
		}
	}

	const double dy = param.deltaY;
	const double dydy = dy*dy;
	const double twodydy = 2*dydy;
	for (int i=map.imin; i <= map.imax; i++) {
		for (int j=map.jinf[i]; j<=map.jsup[i]; j++) {
			bord.b			= map.bord_cell[i][j];
			const char pos		= bord.cotey();
			dvariable sigma		= diffusion_y[i][j];
			dvariable sigmam	= diffusion_y[i][j-1];
			dvariable sigmap	= diffusion_y[i][j+1];
			dvariable vinf		= advection_y[i][j-1];
			dvariable v		= advection_y[i][j];
			dvariable vsup		= advection_y[i][j+1];

			dvarsD[i][j] = d1(pos, sigmam, sigma, vinf, twodydy, dydy, dy);
			dvarsE[i][j] = d2(pos, sigmam, sigma, sigmap, v, twodydy, dydy, dy);
			dvarsF[i][j] = d3(pos, sigma, sigmap, vsup, twodydy, dydy, dy);
		}
	}
	Ybet_comp(map,2*iterationNumber);
}


dvariable CCalpop::d1(const char pos, dvariable& sigmam, dvariable& sigma, dvariable& uinf, const double twodd, const double dd, const double d)
{
	dvariable val = 0;
	switch (pos){
        case SANS:			// si frontieres ouvertes
		case D_FERME:			// ou fermee a droite seulement
        	if (uinf > 0)
                    	val =	-((sigmam+sigma)/twodd)
				-(uinf/d);
		else
			val =	-((sigmam+sigma)/twodd);
		break;
	}
	return val;
}
/*
dvariable CCalpop::d1(const char pos, dvariable& sigmam, dvariable& sigma, dvariable& uinf, const double twodd, const double dd, const double d)
{
	dvariable val = 0;
	switch (pos){
        case SANS:		// si frontieres ouvertes
		case D_FERME:	// ou fermee a droite seulement
                    	val = -((sigmam+sigma)/twodd)-0.5*(uinf+sfabs(uinf))/d;
		break;
	}
	return val;
}
*/


dvariable CCalpop::d2(const char pos, dvariable& sigmam, dvariable& sigma, dvariable& sigmap, dvariable& u, const double twodd, const double dd, const double d)
{
	dvariable val = 0;
	switch (pos){
	        case SANS:
			if (u > 0)
				val = ((sigmam+2*sigma+sigmap)/twodd)
					+(u/d);
			else
				val = ((sigmam+2*sigma+sigmap)/twodd)
					-(u/d);
			break;

		case G_FERME:
			if (u > 0)
				val = ((sigma+sigmap)/twodd)
					+(u/d);
			else
				val = ((sigma+sigmap)/twodd);
			break;

		case D_FERME:
			if (u > 0)
				val = ((sigmam+sigma)/twodd);
			else
				val = ((sigmam+sigma)/twodd)
					-(u/d);
			break;

		case TERRE:
			val = 0;
			break;
	}
	return val;
}
/*
dvariable CCalpop::d2(const char pos, dvariable& sigmam, dvariable& sigma, dvariable& sigmap, dvariable& u, const double twodd, const double dd, const double d)
{
	dvariable val = 0;
	switch (pos){
	        case SANS:
			val = ((sigmam+2*sigma+sigmap)/twodd)+0.5*(u+sfabs(u))/d-0.5*(u-sfabs(u))/d;
		break;

		case G_FERME:
			val = ((sigma+sigmap)/twodd)+0.5*(u+sfabs(u))/d;
		break;

		case D_FERME:
			val = ((sigmam+sigma)/twodd)-0.5*(u-sfabs(u))/d;
		break;

		case TERRE:
			val = 0;
		break;
	}
	return val;
}
*/



dvariable CCalpop::d3(const char pos, dvariable& sigma, dvariable& sigmap, dvariable& usup,const double twodd, const double dd, const double d)
{
	dvariable val = 0;
    	switch (pos) {
        	case SANS:
			case G_FERME:
                		if (usup > 0)
					val = -((sigma+sigmap)/twodd);
				else
					val = -((sigma+sigmap)/twodd)
							+(usup/d);
			break;
	}
	return val;
}
/*
dvariable CCalpop::d3(const char pos, dvariable& sigma, dvariable& sigmap, dvariable& usup,const double twodd, const double dd, const double d)
{
	dvariable val = 0;
    	switch (pos) {
        	case SANS:
			case G_FERME:
				val = -((sigma+sigmap)/twodd)+0.5*(usup-sfabs(usup))/d;
			break;
	}
	return val;
}
*/



