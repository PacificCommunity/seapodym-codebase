#include "SeapodymCoupled.h"
#ifdef SEAPODYM_WITH_DATAPROVIDER
#  include <mpi.h>
#  include <cstdio>
#endif

void SeapodymCoupled::ReadTimeSeriesData(int t, int t_series)
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

void SeapodymCoupled::ReadClimatologyData(int t, int month)
{
	int nbytetoskip = (9 +(3* nlat * nlon) + 12 + ((nlat *nlon)* (month-1))) * 4;
	rw.rbin_input2d(param->strfile_ppmc, map, mat.np1[t], nbi, nbj, nbytetoskip);
	for (int k=0;k<nb_layer;k++){
		rw.rbin_input2d(param->strfile_umc[k], map, mat.un[t][k],   nbi, nbj, nbytetoskip);
		rw.rbin_input2d(param->strfile_vmc[k], map, mat.vn[t][k],   nbi, nbj, nbytetoskip);
		rw.rbin_input2d(param->strfile_tmc[k], map, mat.tempn[t][k],nbi, nbj, nbytetoskip);
		if (!param->type_oxy)
			rw.rbin_input2d(param->strfile_oxymc[k], map, mat.oxygen[t][k],nbi, nbj, nbytetoskip);
	}
	UnitConversions(t);
	if (param->use_sst)
		rw.rbin_input2d(param->strfile_sstmc, map, mat.sst[t], nbi, nbj, nbytetoskip);
	else
		mat.sst[t] = mat.tempn[t][0];
	if (param->use_vld){
		rw.rbin_input2d(param->strfile_vldmc, map, mat.vld[t], nbi, nbj, nbytetoskip);
		mat.vld[t] = mat.vld[t]/1000.0; //in km
	} else
		mat.vld[t] = 1.0;

	if (!param->flag_coupling){
		for (int n=0; n<nb_forage; n++)
			rw.rbin_input2d(param->strfile_Fmc[n], map, mat.forage[t][n], nbi, nbj, nbytetoskip);
	} else {//only simulation mode (no optimization yet)
		for (int n=0; n<nb_forage; n++)
			rw.rbin_input2d(param->strfile_Smc[n], map, mat.mats[n], nbi, nbj, nbytetoskip);
	}
}

void SeapodymCoupled::UnitConversions(int t)
{//Convert units of U and V from m/s to nmi/dt
	mat.un(t) = mat.un(t) * 3600 * 24 * deltaT / 1852;
	mat.vn(t) = mat.vn(t) * (-1) * 3600 * 24 * deltaT / 1852;
	const int imin = map.imin;
	const int imax = map.imax;
	for (int i = imin; i <= imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin ; j <= jmax; j++){
			if (map.carte[i][j]){
				int nlayer = map.carte(i,j);
					for (int k=0 ; k<nb_layer ; k++){
					if (k < nlayer)
						mat.un[t][k][i][j] = mat.un[t][k][i][j] * mat.lat_correction[j];
					else {
						mat.un[t][k][i][j] = 0.0;
						mat.vn[t][k][i][j] = 0.0;
					}
				}
			}
		}
	}
}

void SeapodymCoupled::ReadClimatologyOxy(int t, int t_clm)
{
	int nlevel_oxy = 12;
	if (param->type_oxy==2) nlevel_oxy = 4;
	int nbytetoskip = (9 +(3* nlat * nlon) + nlevel_oxy + ((nlat *nlon)* (t_clm-1))) * 4;
	for (int k=0;k<nb_layer;k++)
		rw.rbin_input2d(param->strfile_oxy[k], map, mat.oxygen[t][k], nbi, nbj, nbytetoskip);
}

#ifdef SEAPODYM_WITH_DATAPROVIDER
// -----------------------------------------------------------------------
// Step 2 helpers: serialise / deserialise all forcing arrays between the
// ADMB matrices and the DataProvider flat shared-memory buffers.
//
// Flat buffer layout (consistent between the two directions):
//   3D field arr[t][i][j]: buf[(t-t0)*mapCells + cell_idx]
//   4D field arr[t][k][i][j]: buf[((t-t0)*ndim + k)*mapCells + cell_idx]
// where cell_idx advances as i: imin..imax, j: jinf[i]..jsup[i].
//
// Note: coupling mode (flag_coupling) reads into mat.mats (a non-time-
// dependent D3_ARRAY) rather than mat.forage. mat.mats is not currently
// wired through DataProvider; a future step can add a "mats" entry.
// -----------------------------------------------------------------------

void SeapodymCoupled::copyForcingToDP(DataProvider* dp, int t0, int nbt)
{
	const int imin = map.imin;
	const int imax = map.imax;

	// --- 3D field helper: arr[t][i][j] ---
	auto ser3 = [&](const char* name, D3_ARRAY& arr) {
		double* buf = dp->getDataPtr(name);
		if (!buf) { printf("[copyForcingToDP] WARNING: no buffer for '%s'\n", name); return; }
		std::size_t idx = 0;
		for (int t = t0; t <= nbt; ++t)
			for (int i = imin; i <= imax; ++i)
				for (int j = map.jinf[i]; j <= map.jsup[i]; ++j)
					buf[idx++] = arr[t][i][j];
	};

	// --- 4D field helper: arr[t][k][i][j] ---
	auto ser4 = [&](const char* name, D4_ARRAY& arr, int ndim) {
		double* buf = dp->getDataPtr(name);
		if (!buf) { printf("[copyForcingToDP] WARNING: no buffer for '%s'\n", name); return; }
		std::size_t idx = 0;
		for (int t = t0; t <= nbt; ++t)
			for (int k = 0; k < ndim; ++k)
				for (int i = imin; i <= imax; ++i)
					for (int j = map.jinf[i]; j <= map.jsup[i]; ++j)
						buf[idx++] = arr[t][k][i][j];
	};

	ser3("np1",  mat.np1);
	ser3("sst",  mat.sst);
	ser3("vld",  mat.vld);
	ser4("un",     mat.un,     nb_layer);
	ser4("vn",     mat.vn,     nb_layer);
	ser4("tempn",  mat.tempn,  nb_layer);
	ser4("oxygen", mat.oxygen, nb_layer);
	if (!param->flag_coupling)
		ser4("forage", mat.forage, nb_forage);
}

void SeapodymCoupled::copyForcingFromDP(DataProvider* dp, int t0, int nbt)
{
	const int imin = map.imin;
	const int imax = map.imax;

	// --- 3D field helper: arr[t][i][j] ---
	auto des3 = [&](const char* name, D3_ARRAY& arr) {
		const double* buf = dp->getDataPtr(name);
		if (!buf) { printf("[copyForcingFromDP] WARNING: no buffer for '%s'\n", name); return; }
		std::size_t idx = 0;
		for (int t = t0; t <= nbt; ++t)
			for (int i = imin; i <= imax; ++i)
				for (int j = map.jinf[i]; j <= map.jsup[i]; ++j)
					arr[t][i][j] = buf[idx++];
	};

	// --- 4D field helper: arr[t][k][i][j] ---
	auto des4 = [&](const char* name, D4_ARRAY& arr, int ndim) {
		const double* buf = dp->getDataPtr(name);
		if (!buf) { printf("[copyForcingFromDP] WARNING: no buffer for '%s'\n", name); return; }
		std::size_t idx = 0;
		for (int t = t0; t <= nbt; ++t)
			for (int k = 0; k < ndim; ++k)
				for (int i = imin; i <= imax; ++i)
					for (int j = map.jinf[i]; j <= map.jsup[i]; ++j)
						arr[t][k][i][j] = buf[idx++];
	};

	des3("np1",  mat.np1);
	des3("sst",  mat.sst);
	des3("vld",  mat.vld);
	des4("un",     mat.un,     nb_layer);
	des4("vn",     mat.vn,     nb_layer);
	des4("tempn",  mat.tempn,  nb_layer);
	des4("oxygen", mat.oxygen, nb_layer);
	if (!param->flag_coupling)
		des4("forage", mat.forage, nb_forage);
}
#endif // SEAPODYM_WITH_DATAPROVIDER

// -----------------------------------------------------------------------

void SeapodymCoupled::ReadAll(int tstart, int tend, int offset, void* dpVoid)
{
	// Step 2: only the shmRoot on each node reads from disk.
	// dpVoid is non-null only when called from SeapodymCohort::OnRunFirstStep
	// with the MPI cohort target.
#ifdef SEAPODYM_WITH_DATAPROVIDER
	DataProvider* dp = static_cast<DataProvider*>(dpVoid);
	const bool doRead = (!dp || dp->isShmRoot());
#else
	const bool doRead = true;
#endif

	if (doRead) {
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
				ReadTimeSeriesData(t, t_series);
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
	} // end doRead

	// Step 2: broadcast forcing data within each shared-memory node.
#ifdef SEAPODYM_WITH_DATAPROVIDER
	if (dp) {
		if (dp->isShmRoot())
			copyForcingToDP(dp, tstart, tend);
		MPI_Barrier(dp->getShmComm());
		if (!dp->isShmRoot())
			copyForcingFromDP(dp, tstart, tend);
	}
#endif
}

void SeapodymCoupled::RestoreDistributions(ivector& nb_age_built)
{
	//to read initial distributions
	t_count  += nbt_spinup_tuna;
	for (int sp=0; sp<nb_species; sp++){
		const int nb_ages = param->sp_nb_cohorts[sp];
		string fileCohorts = param->str_dir_init + param->sp_name[sp] + "_cohorts.dym";
		int nlevel = 0;
		int nlat = param->nlat;
		int nlon = param->nlong;
		rw.rbin_headpar(fileCohorts, param->nlong, param->nlat, nlevel);
		/*if (nlevel != nb_ages)
			cout << "WARNING: in file " << fileCohorts << " number of cohorts is " <<
				nlevel << " != " << nb_ages << " in the current parfile!" << endl;*/
		for (int a=0; a<nb_ages; a++){
			int nbytetoskip = (9 +(3* nlat * nlon) + nlevel + ((nlat *nlon)* a)) * 4;
			rw.rbin_input2d(fileCohorts, map, mat.init_density_species[sp][a], nbi, nbj, nbytetoskip);
		}
		if (param->nlat!=nlat) {
			cerr << "WARNING: different nlat in file "<< fileCohorts << " " << param->nlat <<  endl;
			param->nlat = nlat;
		}
		if (param->nlong!=nlon)  {
			cerr << "WARNING: different nlon in file "<< fileCohorts << " " << param->nlong <<  endl;
			param->nlong = nlon;
		}
		nb_age_built[sp] = nb_ages-1;
	}
}
