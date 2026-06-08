#include "SeapodymCoupled.h"
#ifdef SEAPODYM_WITH_DATAPROVIDER
#  include <mpi.h>
#  include <cstdio>
#endif

// ---------------------------------------------------------------------------
// Helper: read a binary 2D grid into a FlatField3D time slice.
// rbin_input2d writes into a ragged dmatrix; we use the field's pre-
// allocated stage_ dmatrix and pack into the flat buffer afterwards.
// ---------------------------------------------------------------------------
// do-while(0) ensures the macro is a single statement safe inside if/else and loops.
#define READ3(field, file, t, skip) \
    do { \
        rw.rbin_input2d(file, map, mat.field.stage(), nbi, nbj, skip); \
        mat.field.commit(t); \
    } while(0)

// Same for a (t, k) slice of a FlatField4D.
#define READ4(field, file, t, k, skip) \
    do { \
        rw.rbin_input2d(file, map, mat.field.stage(), nbi, nbj, skip); \
        mat.field.commit(t, k); \
    } while(0)


void SeapodymCoupled::ReadTimeSeriesData(int t, int t_series)
{
	int nbytetoskip = (9 +(3* nlat * nlon) + param->nlevel + ((nlat *nlon)* (t_series-1))) * 4;
	READ3(np1, param->strfile_pp, t, nbytetoskip);

	for (int k=0;k<nb_layer;k++){
		READ4(un,     param->strfile_u[k], t, k, nbytetoskip);
		READ4(vn,     param->strfile_v[k], t, k, nbytetoskip);
		READ4(tempn,  param->strfile_t[k], t, k, nbytetoskip);
		if (!param->type_oxy)
			READ4(oxygen, param->strfile_oxy[k], t, k, nbytetoskip);
	}
	UnitConversions(t);

	if (param->use_sst)
		READ3(sst, param->strfile_sst, t, nbytetoskip);
	else
		mat.sst[t] = mat.tempn[t][0];  // RawSlice2D memcpy

	if (param->use_vld){
		READ3(vld, param->strfile_vld, t, nbytetoskip);
		mat.vld[t].scaleInPlace(1.0/1000.0); //in km (mtl is in mt/km2)
	} else
		mat.vld[t] = 1.0; //won't be used

	if (!param->flag_coupling){
		for (int n=0; n<nb_forage; n++){
			READ4(forage, param->strfile_F[n], t, n, nbytetoskip);
			mat.forage[t][n].addInPlace(0.000001); //avoid zero habitat index
		}
	} else {
		for (int n=0; n<nb_forage; n++)
			rw.rbin_input2d(param->strfile_S[n], map, mat.mats[n], nbi, nbj, nbytetoskip);
	}
}

void SeapodymCoupled::ReadClimatologyData(int t, int month)
{
	int nbytetoskip = (9 +(3* nlat * nlon) + 12 + ((nlat *nlon)* (month-1))) * 4;
	READ3(np1, param->strfile_ppmc, t, nbytetoskip);

	for (int k=0;k<nb_layer;k++){
		READ4(un,     param->strfile_umc[k], t, k, nbytetoskip);
		READ4(vn,     param->strfile_vmc[k], t, k, nbytetoskip);
		READ4(tempn,  param->strfile_tmc[k], t, k, nbytetoskip);
		if (!param->type_oxy)
			READ4(oxygen, param->strfile_oxymc[k], t, k, nbytetoskip);
	}
	UnitConversions(t);

	if (param->use_sst)
		READ3(sst, param->strfile_sstmc, t, nbytetoskip);
	else
		mat.sst[t] = mat.tempn[t][0];

	if (param->use_vld){
		READ3(vld, param->strfile_vldmc, t, nbytetoskip);
		mat.vld[t].scaleInPlace(1.0/1000.0);
	} else
		mat.vld[t] = 1.0;

	if (!param->flag_coupling){
		for (int n=0; n<nb_forage; n++)
			READ4(forage, param->strfile_Fmc[n], t, n, nbytetoskip);
	} else {
		for (int n=0; n<nb_forage; n++)
			rw.rbin_input2d(param->strfile_Smc[n], map, mat.mats[n], nbi, nbj, nbytetoskip);
	}
}

void SeapodymCoupled::UnitConversions(int t)
{//Convert units of U and V from m/s to nmi/dt
	mat.un.scaleTimeStep(t,  3600.0 * 24 * deltaT / 1852);
	mat.vn.scaleTimeStep(t, -3600.0 * 24 * deltaT / 1852);
	// Apply latitude correction and zero sub-surface cells beyond nlayer.
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
		READ4(oxygen, param->strfile_oxy[k], t, k, nbytetoskip);
}

// -----------------------------------------------------------------------

void SeapodymCoupled::ReadAll(int tstart, int tend, int offset, void* dpVoid)
{
	// Step 3: each FlatField already aliases the DataProvider shared-memory
	// window directly, so only the shmRoot needs to read from disk.
	// After reading, a single MPI_Barrier lets non-shmRoot workers proceed.
#ifdef SEAPODYM_WITH_DATAPROVIDER
	DataProvider* dp = static_cast<DataProvider*>(dpVoid);
	const bool doRead = (!dp || dp->isShmRoot());
#else
	const bool doRead = true;
#endif

	if (doRead) {
		int jday = 0;
		int t = tstart;

		if (param->type_oxy==1 && month==past_month)
			ReadClimatologyOxy(1, month);
		if (param->type_oxy==2 && qtr==past_qtr)
			ReadClimatologyOxy(1, qtr);

		for (; t<=tend; t++){
			getDate(jday,t);
			if (t > nbt_building) {
				t_series = t - nbt_building + nbt_start_series + offset;
				ReadTimeSeriesData(t, t_series);
			}
			else if ((t <= nbt_building) && (month != past_month)) {
				ReadClimatologyData(t, month);
			}
			if (param->type_oxy==1 && month != past_month) {
				ReadClimatologyOxy(t, month);
			}
			if (param->type_oxy==2 && qtr != past_qtr) {
				ReadClimatologyOxy(t, qtr);
			}
		}
	} // end doRead

	// Ensure shmRoot has finished writing into the shared-memory window
	// before non-shmRoot workers read from it.
#ifdef SEAPODYM_WITH_DATAPROVIDER
	if (dp) MPI_Barrier(dp->getShmComm());
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
