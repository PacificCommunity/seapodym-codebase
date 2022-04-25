#include "calpop.h"

void CCalpop::precaldia(const CParam& param, const PMap& map, CMatrices& mat)
{
	const double sigma_fcte = param.sigma_fcte;

	for (int i = map.imin; i <= map.imax; i++){	
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		dvector& diffusion_x = mat.diffusion_x[i];
		dvector& diffusion_y = mat.diffusion_y[i];
		dvector& advection_x = mat.advection_x[i];
		dvector& advection_y = mat.advection_y[i];
		for (int j = jmin; j <= jmax; j++){
			if (map.carte[i][j]){
				diffusion_x[j] = sigma_fcte * mat.lat_correction[j];
				diffusion_y[j] = sigma_fcte ;

				// sur l'axe des x : u (deja corrige par rapport a la latitude a la lecture)
				advection_x[j] = mat.u[i][j];
				advection_y[j] = mat.v[i][j];
			}
		}
	}
}
