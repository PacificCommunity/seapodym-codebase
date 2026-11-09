#include "calpop.h"

///Forward function for:
///predicting catch by fishery based on fishing effort. 

const double sq_degree   = 12347.65; // sq.km in sq.degree (without lat_correction)

///Fishing areas of the effort records of fishery f on the model grid: the coefficients of
///CParam::afcoef on the ocean mask and their sum afr. They do not depend on age, so they are
///computed once per fishery and month and used for all age classes, in the forward and the
///adjoint codes (same expressions and summation order as when computed for each age).
void CCalpop::fishing_area_comp(const PMap& map, CParam& param, const dmatrix& effort, const dmatrix& efflon, const dmatrix& efflat, const int f, std::vector<fishing_area>& fa)
{
	fa.clear();
	const double cell_area_deg = param.deltaX*param.deltaY/(60*60);
	const float reso = param.fishery_reso(f);
	int m = reso*60.0/param.deltaX+2;
	int n = reso*60.0/param.deltaY+2;
	dmatrix af(0,m-1,0,n-1);
	for (int i = map.imin; i <= map.imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			if (map.carte(i,j) && effort(i,j)){
				fishing_area e;
				e.i = i; e.j = j; e.effort = effort(i,j);
				int ki = 0; int kj = 0;
				//returns fishing area in the cell to be used in biomass computation
				param.afcoef(efflon(i,j),efflat(i,j), af, ki, kj, f);

				//the factor 1/(rx*rw) = sum(af)/cell_area_deg
				e.afr = 0.0;
				for (int ii=0; ii<m; ii++)
					for (int jj=0; jj<n; jj++)
						if (af(ii,jj) && map.carte(ii+ki,jj+kj)){
							e.afr += af(ii,jj)/cell_area_deg;
							e.ii.push_back(ii+ki); e.jj.push_back(jj+kj); e.af.push_back(af(ii,jj));
						}
				fa.push_back(e);
			}
		}
	}
}

void CCalpop::predicted_catch_fishery_comp(const PMap& map, CParam& param, VarMatrices& mat, /* dmatrix total_mort, */const int f, const int k, const int sp, const int age, const dmatrix& uu, const int step_count)
{
	double catchability = param.q_sp_fishery[sp][k] + param.q_slope_fishery[k] * step_count;
	if (catchability<0) catchability = 0; // in case if negative trend has been chosen
	//const double selectivity = param.selectivity_comp(sp,age,f,k);
	//const double selectivity = param.selectivity_comp(sp,age,f,k,step_count);
	const double selectivity = Selectivity(sp,f,age);

	//catch_units=0 if catch is in nb; =1 if in metric tons
	const int catch_units = param.fishery_catch_units[f];
	//Note, variable weight is either =real weight in mt or =1 if predicted catch is in numbers
	const double weight      = pow(1e-3*param.weight[sp][age],catch_units);

	//compute catch in tonnes per cell Nb/km^2 * mt = W(mt)*dx*dy*(1.852)^2/cell
	//to translate length frequencies into number of fish caugth in cell: Nb/km^2 = Nb*(1.852)^2*dx*dy/cell

	const int nb_region  = param.nb_region;
	const int nb_region_sp   = param.nb_region_sp_B[sp];
	const double s_c     = selectivity*catchability;
	const double w_area  = weight*sq_degree;
	//Catch equation: we will compute and fit predicted catch at the resolution of the fishing data
	//assume we have E(m,n) at coarser resolution than ours: rx>=deltaX, ry>=deltaY
	//C(m,n) = q*sum_ij(E(i,j)*B(i,j)) = [if E(i,j)=const=E(m,n)/(rx*ry)] = q*E/(rx*ry)*sum_ij(B(i,j))
	//where rx*ry is the size of fisheries data cell in deltaX and deltaY correspondingly
	//(e.g. if deltaX=deltaY=2deg and fishery_reso=5deg => rx*ry = (5/2)*(5/2) = 6.25)

	//the fishing areas of the effort cells (fa) are computed by fishing_area_comp
	//before the loop over age classes
	for (size_t nc = 0; nc < fa.size(); nc++){
		const fishing_area& e = fa[nc];
		const int i = e.i; const int j = e.j;
		const double afr = e.afr;

		double fish  = 0;
		//Nb. of fish in a given (lon,lat)+/-freso/2 = cell_area*sum_ij(af*N(i,j)):
		for (size_t p = 0; p < e.af.size(); p++)
			fish += e.af[p]*uu(e.ii[p],e.jj[p])/(afr*mat.lat_correction[e.jj[p]]);

		const double F = afr*param.smin1_c(e.effort*s_c/afr);

		// total catch in weight (tonnes)
		mat.dvarCatch_est(sp,k).elem_value(i,j) += F * fish * w_area ;
		mat.catch_est[sp][k][i][j] += F * fish * w_area;

		//Predicted Catch-at-age by species, fleet and region (the regions containing the cell)
		if (nb_region){
			const std::vector<int>& regs = map.regions_of_cell(i,j);
			for (size_t n=0; n<regs.size(); n++){
				const int a = regs[n];
				if (a >= nb_region_sp) continue;
				//attention: area_sp_B should not contain zero IDs!
				int r = param.area_sp_B[sp][a]-1;
				mat.C_N_sp_age_fishery[sp][age][k][r] += F * fish * sq_degree;
			}
		} else if (nb_region_sp)
			mat.C_N_sp_age_fishery[sp][age][k][0] += F * fish * sq_degree;
	}
}
