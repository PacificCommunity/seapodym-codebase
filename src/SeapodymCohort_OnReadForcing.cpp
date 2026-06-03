#include "SeapodymCoupled.h"
#include "SeapodymCohort.h"

void SeapodymCohort::ReadTimeSeriesData(int t, int t_series, DataProvider dataProvider)	
{
	int nbytetoskip = (9 +(3* nlat * nlon) + param->nlevel + ((nlat *nlon)* (t_series-1))) * 4;
	rw.rbin_input2d(param->strfile_pp, map, mat.np1[t], nbi, nbj, nbytetoskip);

	for (int k=0;k<nb_layer;k++){
		rw.rbin_input2d(param->strfile_u[k], map, mat.un[t][k],   nbi, nbj, nbytetoskip);
		rw.rbin_input2d(param->strfile_v[k], map, mat.vn[t][k],   nbi, nbj, nbytetoskip);
		rw.rbin_input2d(param->strfile_t[k], map, mat.tempn[t][k],nbi, nbj, nbytetoskip);
		if (!param->type_oxy)
			rw.rbin_input2d(param->strfile_oxy[k], map, mat.oxygen[t][k],nbi, nbj, nbytetoskip);

	}
	UnitConversions(t);
	if (param->use_sst)
		rw.rbin_input2d(param->strfile_sst, map, mat.sst[t], nbi, nbj, nbytetoskip);
	else 
		mat.sst[t] = mat.tempn[t][0];

	if (param->use_vld){
		rw.rbin_input2d(param->strfile_vld, map, mat.vld[t], nbi, nbj, nbytetoskip);
		mat.vld[t] = mat.vld[t]/1000.0; //in km (mtl is in mt/km2). DO NOT change the units here without changing the use of vld in Ha computation
	} else 
		mat.vld[t] = 1.0; //wont be used

	if (!param->flag_coupling){
		for (int n=0; n<nb_forage; n++){
			rw.rbin_input2d(param->strfile_F[n], map, mat.forage[t][n], nbi, nbj, nbytetoskip);
			
			mat.forage(t,n) = mat.forage(t,n)+0.000001; //to avoid zero habitat index, expecially in fine resolution simulations F can be zero!!!
		}
	} else {
		for (int n=0; n<nb_forage; n++)
			rw.rbin_input2d(param->strfile_S[n], map, mat.mats[n], nbi, nbj, nbytetoskip);	
	}
}


void SeapodymCohort::ReadAll(int tstart, int tend, int offset, DataProvider dataProvider)
{
	//cout << "Reading all forcing variables at once... " << endl;
	int jday = 0;
	int t = tstart;

	//need to read oxygen in case if month==past_month
	//(otherwise we may not have it for the first time steps)
	if (param->type_oxy==1 && month==past_month)
		ReadClimatologyOxy(1, month);
	//need to read oxygen in case if qtr==past_qtr 
	if (param->type_oxy==2 && qtr==past_qtr)
		ReadClimatologyOxy(1, qtr);

	for (; t<=tend; t++){
		getDate(jday,t);
		//----------------------------------------------//
		//	DATA READING SECTION: U,V,T,O2,PP	//
		//----------------------------------------------//
		if (t > nbt_building) {
			//TIME SERIES 
			t_series = t - nbt_building + nbt_start_series + offset;
			ReadTimeSeriesData(t, t_series, dataProvider);	
		}
		else if ((t <= nbt_building) && (month != past_month)) {
			//AVERAGED CLIMATOLOGY DATA
			ReadClimatologyData(t, month);
		}
		if (param->type_oxy==1 && month != past_month) {
			//MONTHLY O2
			ReadClimatologyOxy(t, month);
		}
		if (param->type_oxy==2 && qtr != past_qtr) {
			//QUARTERLY O2
			ReadClimatologyOxy(t, qtr);
		}
	}
}