//#include "StdAfx.h"
#include "ReadWrite.h"
#include "SaveTimeArea.h"
#include "Numfunc.h"
#include "Date.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#ifndef DYM_OUTPUT_TYPE
	#define DYM_OUTPUT_TYPE float
#endif
#ifndef DYM_OUTPUT_MY_TYPE
	#define DYM_OUTPUT_MY_TYPE double
#endif

#ifndef DYM_INPUT_TYPE
	#define DYM_INPUT_TYPE float
#endif

///////////////////////////***************//////////////////////
////////////////////////// LECTURE BINAIRE//////////////////////
///////////////////////////***************//////////////////////
// lecture des parametres du fichier permettant de creer leas tableaux
void CReadWrite::rbin_headpar(string file_in, int &nlong, int &nlat, int &nlevel)
{
		
//cerr << "Force Exit[" << __FILE__ << ':' << __LINE__ << "]: function hasn't been fixed yet\n";
	ifstream litbin(file_in.c_str(), ios::binary|ios::in);
	if (!litbin)
	{
		cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: Unable to read file \"" << file_in << "\"\n";
		exit(1);
	}

 	const int nbytes= 4;

	// for new dym format, skip the first 4 parameters (Idformat, Idfunc, Minval, Maxval)
	const int nbytetoskip = 4 * nbytes;  // Patrick 21Oct2004 
	litbin.seekg(nbytetoskip, ios::beg); // Patrick 21Oct04 

	int bufin;
	litbin.read((char *)&bufin,nbytes);
	nlong = bufin;
 	litbin.read((char *)&bufin,nbytes);
	nlat = bufin;
	litbin.read((char *)&bufin,nbytes);
	nlevel = bufin;

	litbin.close();
}
/*
void CReadWrite::set_par_domain(CParam& param)
{//initialization of domain parameters

	latitudeMax  = param.latitudeMax;
	longitudeMin = param.longitudeMin;
	deltaX = (param.deltaX/60);
	deltaY = (param.deltaY/60);
}
*/

// reading the min/maxval and update them if necessary
void CReadWrite::rwbin_minmax(string file_io, double minvalstep, double maxvalstep)//Gael Nov04
{
	fstream rwbin(file_io.c_str(), ios::binary|ios::in|ios::out);

	if (!rwbin) {
		cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: Unable to read-write file \"" << file_io << "\"\n";
		exit(1);
	}

 	const int nbytes= 4;
	//skip the first 2 parameters (Idformat, Idfunc)
	const int nbytetoskip = 2 * nbytes; 
	rwbin.seekg(nbytetoskip, ios::beg); 

 	const int sizeofDymOutputType = sizeof(DYM_OUTPUT_TYPE);
	DYM_OUTPUT_TYPE buffl;

	rwbin.read(( char *)&buffl, sizeofDymOutputType);	
	const double minval = buffl;

	if (minvalstep<minval) {

		rwbin.seekp(-sizeofDymOutputType, ios::cur);
		buffl = minvalstep;
		rwbin.write(( char *)&buffl, sizeofDymOutputType);
	}

	rwbin.read(( char *)&buffl, sizeofDymOutputType);	
	const double maxval = buffl;

	if (maxvalstep>maxval) {

		rwbin.seekp(-sizeofDymOutputType, ios::cur);
		buffl = maxvalstep;
		rwbin.write(( char *)&buffl, sizeofDymOutputType);
	}

	rwbin.close();
}

//---------------------------------------
// Lecture de l'entete
// nlong, nlat, nlevel, startdate, enddate
// val des longitudes, latitudes, dates
// mask
//---------------------------------------
void CReadWrite::rbin_header(string file_in, string &idformat, int &idfunc, double &minval, double &maxval,  
					int nlong, int nlat, int nlevel, double &startdate, double &enddate,
					DMATRIX& xlon, DMATRIX& ylat, DVECTOR& zlevel, IMATRIX& msksp)
{
//cerr << "Force Exit[" << __FILE__ << ':' << __LINE__ << "]: function hasn't been fixed yet\n";
	ifstream litbin(file_in.c_str(), ios::binary|ios::in);
	if (!litbin)
	{
		cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: Unable to read file \"" << file_in << "\"\n";
		exit(1);
	}
	const int nbytes= 4;
	
	int bufin;

 	const int sizeofDymInputType = sizeof(DYM_INPUT_TYPE);
	DYM_INPUT_TYPE buffl;

        char *bufch;
        int bufchlen = nbytes+1;
        bufch = new char [bufchlen];
        bufch[bufchlen-1] = '\0';

	  	
	/////////////////////////////////////////////
	// reading new dym format	: Patrick 21Oct04 
	litbin.read(( char *)bufch,nbytes);
	idformat = bufch;
	litbin.read(( char *)&bufin,nbytes);	
	idfunc = bufin;
	litbin.read(( char *)&buffl,sizeofDymInputType);	
	minval = buffl;
	litbin.read(( char *)&buffl,sizeofDymInputType);	
	maxval = buffl;
	/////////////////////////////: Patrick 21Oct04	
	litbin.read(( char *)&bufin,nbytes);
	nlong = bufin;
	litbin.read(( char *)&bufin,nbytes);
	nlat = bufin;
	litbin.read(( char *)&bufin,nbytes);
	nlevel = bufin;
	litbin.read(( char *)&buffl,sizeofDymInputType);
	startdate = buffl;
	litbin.read(( char *)&buffl,sizeofDymInputType);
	enddate = buffl;

	///////////////////////////////////////////
	// reading new dym format:  Patrick 21Oct04
	//---------------------------------------
	//xlon	valeurs des longitudes
	//---------------------------------------
	for (int i=0;i<nlat;i++)
	{
		for (int j=0;j<nlong;j++)
		{
			litbin.read(( char *)&buffl,sizeofDymInputType);
			xlon[i][j] = buffl;
		}
	}

	//---------------------------------------
	//ylat	valeurs des latitudes
	//---------------------------------------
	for (int i=0;i<nlat;i++)
	{
		for (int j=0;j<nlong;j++)
		{
			litbin.read(( char *)&buffl,sizeofDymInputType);
			ylat[i][j] = buffl;
		}
	}

	///////////////////////////: Patrick 21Oct04
	//---------------------------------------
	//zlevel: contient les valeurs des dates
	//---------------------------------------
	for (int k=0;k<nlevel;k++)
	{
		litbin.read(( char *)&buffl,sizeofDymInputType);
		zlevel[k] = buffl;
	}
	enddate = zlevel[nlevel-1];

	//---------------------------------------
	//msksp
	//---------------------------------------
	for (int i=0;i<nlat;i++)
	{
		for (int j=0;j<nlong;j++)
		{
			litbin.read(( char *)&bufin,nbytes);
			msksp[i][j] = bufin;
		}
	}
	delete [] bufch;

	litbin.close();
}
//////////////////////////////////////////////////////////////////////
void CReadWrite::rbin_mat2d(string file_in, PMap& map, DMATRIX& mat2d, int nlat, int nlong, int nbytetoskip)
{
	ifstream litbin(file_in.c_str(), ios::binary | ios::in);
	
	if (!litbin)
	{
		cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: Unable to read file \"" << file_in << "\"\n";
		exit(1);
	}

	//---------------------------------------
	// Reading the 2d matrix, written by seapodym_coupled
	//---------------------------------------
	litbin.seekg(nbytetoskip, ios::cur);

 	const int sizeofDymInputType = sizeof(DYM_OUTPUT_MY_TYPE);
	DYM_OUTPUT_MY_TYPE buf;

	mat2d.initialize();
	for (int i=0;i<nlat;i++)
	{
		for (int j=0;j<nlong;j++)
		{
			litbin.read(( char *)&buf,sizeofDymInputType);
			if (map.carte(i,j))
				mat2d[i][j]= buf;
		}
	}

	litbin.close();
}

//////////////////////////////////////////////////////////////////////
void CReadWrite::rbin_input2d(string file_in, PMap& map, DMATRIX& mat2d, int nbi, int nbj, int nbytetoskip)
{
//cerr << "Force Exit[" << __FILE__ << ':' << __LINE__ << "]: function hasn't been fixed yet\n";
////////////// lecture fichier binaire
	ifstream litbin(file_in.c_str(), ios::binary | ios::in);
	
	if (!litbin)
	{
		cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: Unable to read file \"" << file_in << "\"\n";
		exit(1);
	}

	//---------------------------------------
	// Reading the 2d matrix
	//---------------------------------------
	litbin.seekg(nbytetoskip, ios::cur);


 	const int sizeofDymInputType = sizeof(DYM_INPUT_TYPE);
	DYM_INPUT_TYPE buf;

	mat2d.initialize();
	for (int j=1;j<nbj-1;j++)
	{
		for (int i=1;i<nbi-1;i++)
		{
			litbin.read(( char *)&buf,sizeofDymInputType);
			if (map.carte(i,j)){
				mat2d[i][j]= buf;
				//if (mat2d[i][j]<=-99)//indicates invalid value
				if (mat2d[i][j]<=-9)//tmp: ecco
					mat2d[i][j] = 0;
			}
		}
	}

	litbin.close();
}
/*
void CReadWrite::rbin_input2d(string file_in, const imatrix& carte, DMATRIX& mat2d, int nbi, int nbj, int nbytetoskip)
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
			if (carte(i,j)){
				mat2d[i][j]= buf;
				if (mat2d[i][j]<=-99)//indicates invalid value
					mat2d[i][j] = 0;
			}
		}
	}
	litbin.close();
}
*/
///////////////////////////***************//////////////////////
////////////////////////// ECRITURE BINAIRE//////////////////////
///////////////////////////***************//////////////////////
void CReadWrite::wbin_header(string file_out, string &idformat, int &idfunc, double &minval, double &maxval, 
					int nlong, int nlat, int nlevel, double &startdate, double &enddate,
					const DMATRIX& xlon, const DMATRIX& ylat, const DVECTOR& zlevel, const IMATRIX& msksp)
{
	ofstream ecritbin(file_out.c_str(), ios::binary|ios::out);

	if (!ecritbin)
	{
		cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: Unable to write file \"" << file_out << "\"\n";
		exit(1);
	}

	const int nbytes= 4;
	int bufin;
 	const int sizeofDymOutputType = sizeof(DYM_OUTPUT_TYPE);
	DYM_OUTPUT_TYPE buf;

	//---------------------------------------
	// Write idformat, idfunc, minval, maxval	//Patrick 21Oct04
	//---------------------------------------
	//ecritbin.write(idformat.c_str(),nbytes);
	ecritbin << idformat.c_str();
	bufin = idfunc;
	ecritbin.write(( char *)&bufin,nbytes);
	buf = minval;
	ecritbin.write(( char *)&buf, sizeofDymOutputType);
	buf = maxval;
	ecritbin.write(( char *)&buf, sizeofDymOutputType);

	//---------------------------------------
	// Write nlong, nlat, nlevel, stardate, enddate
	//---------------------------------------
	bufin = nlong;
	ecritbin.write(( char *)&bufin,nbytes);
	bufin = nlat;
	ecritbin.write(( char *)&bufin,nbytes);
	bufin = nlevel;
	ecritbin.write(( char *)&bufin,nbytes);
	buf = startdate;
	ecritbin.write(( char *)&buf, sizeofDymOutputType);
	buf = enddate;
	ecritbin.write(( char *)&buf, sizeofDymOutputType);
	//---------------------------------------
	//xlon	//Patrick 21Oct04
	//---------------------------------------
	for (int i=0;i<nlat;i++)
	{ 
		for (int j=0;j<nlong;j++)
		{
			buf = xlon[i][j];
			ecritbin.write(( char *)&buf, sizeofDymOutputType);
		}
	}
	//---------------------------------------
	//ylat //Patrick 21Oct04
	//---------------------------------------
	for (int i=0;i<nlat;i++)
	{
		for (int j=0;j<nlong;j++)
		{
			buf = ylat[i][j];
			ecritbin.write(( char *)&buf, sizeofDymOutputType);	
		}
	}
	//---------------------------------------
	//zlevel
	//---------------------------------------
	for (int k=0;k<nlevel;k++)
	{
		buf = zlevel[k];
		ecritbin.write(( char *)&buf, sizeofDymOutputType);
	}

	//---------------------------------------
	//msksp
	//---------------------------------------
	for (int i=0;i<nlat;i++)
	{
		for (int j=0;j<nlong;j++)
		{
			bufin = msksp[i][j];
			ecritbin.write(( char *)&bufin,nbytes);
		}
	}
	//---------------------------------------

	ecritbin.close();
}

////////////// ecriture matrice
void CReadWrite::wbin_mat2d(string file_out, const DMATRIX& mat2d, int nlat, int nlong, bool FILEMODE)
{
	ofstream ecritbin;
	if (FILEMODE)
		ecritbin.open(file_out.c_str(), ios::binary|ios::app);
	else
		ecritbin.open(file_out.c_str(), ios::binary|ios::out);

        if (!ecritbin)
        {
                cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: Unable to write file \"" << file_out << "\"\n";
                exit(1);
        }
                                                                                                                                                                                                     
 	const int sizeofDymOutputType = sizeof(DYM_OUTPUT_MY_TYPE);
	DYM_OUTPUT_MY_TYPE buf;
		
	//---------------------------------------
	// Write 2d matrix
	//---------------------------------------
	for (int i=0;i<nlat;i++)
	{
		for (int j=0;j<nlong;j++)
		{
			buf = mat2d[i][j];
			ecritbin.write(( char *)&buf, sizeofDymOutputType);
		}
	}

	ecritbin.close();
}

//////////////////////////////////////////////////////////
//Ecriture de matrice transposee
/////////////////////////////////////////////////////////////
void CReadWrite::wbin_transpomat2d(string file_out, const DMATRIX& mat2d, int nlong, int nlat, bool FILEMODE)
{
//cerr << "Force Exit[" << __FILE__ << ':' << __LINE__ << "]: function hasn't been fixed yet\n";
////////////// ecriture matrice
	ofstream ecritbin;
	if (FILEMODE)
		ecritbin.open(file_out.c_str(), ios::binary|ios::app);
	else
		ecritbin.open(file_out.c_str(), ios::binary|ios::out);

        if (!ecritbin)
        {
                cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: Unable to write file \"" << file_out << "\"\n";
                exit(1);
        }
                                                                                                                                                                                                     
 	const int sizeofDymOutputType = sizeof(DYM_OUTPUT_TYPE);
	DYM_OUTPUT_TYPE buf;
		
	//---------------------------------------
	// Write 2d matrix
	//---------------------------------------
	for (int j=0;j<nlat;j++)
	{
		for (int i=0;i<nlong;i++)
		{
			buf = (float)mat2d[i][j];
			ecritbin.write(( char *)&buf, sizeofDymOutputType);
		}
	}

	ecritbin.close();
}

////////////////////////////////////
void CReadWrite::InitSepodymFileDym(CParam& param, CMatrices& mat, int nb_mo, DVECTOR& zlevel, const IMATRIX& msksp)
{
//const int nlon = param.nlong;
//const int nlat = param.nlat;
//EFR.allocate(0, nlon+1, 0, nlat+1);

	double minval = (double) 3.4e38;
	double maxval = (double) -3.4e38;
	const int nb_species = param.get_nbspecies();
	const int nb_fishery = param.get_nbfishery();

	// Recording 2 variables by forage component:
	// 0- Fprod predicted spatial distribution of forage production
	// 1- Fbiom predicted spatial distribution of forage biomass
	//variable name
/*	int nb_var_F= 1;
	vector<string> var_F_name;
	//var_F_name.push_back("Fprod");
	var_F_name.push_back("Fbiom");

	for (int n=0; n<nb_forage;n++){
		for (int var=0; var<nb_var_F;var++){
			//Generate the name of the files
			dymFileFPred.push_back(param.strdir_output+var_F_name[var]+"_"+param.frg_name[n]+".dym");
			// Call the function to write the header of each file
			wbin_header(dymFileFPred[var+(nb_var_F*n)], param.idformat, param.idfunc, minval, maxval, 
				param.nlong, param.nlat, nb_mo, zlevel[0], zlevel[nb_mo -1],
				mat.xlon, mat.ylat, zlevel, msksp);
		}
	}
*/
	vector<string> var_name;

	if (nb_species){
		for (int sp=0; sp<nb_species;sp++){
			int nb_var = 12;
			var_name.push_back("Ha_first_maturity");
			var_name.push_back("Ha_oldest");
			var_name.push_back("Vtot_x");
			var_name.push_back("Vtot_y");
			var_name.push_back("speed");
			var_name.push_back("diffusion");
			var_name.push_back("larve");
			var_name.push_back("juvnl");
			var_name.push_back("young");
			var_name.push_back("recru");
			var_name.push_back("adult");
			var_name.push_back("totbm");
			for (int var=0; var<nb_var; var++){
				dymFileSpPred.push_back(param.strdir_output+param.sp_name[sp]+"_"+var_name[var]+".dym");
				wbin_header(dymFileSpPred[var], param.idformat, param.idfunc, minval, maxval, 
						param.nlong, param.nlat, nb_mo, zlevel[0], zlevel[nb_mo -1],
						mat.xlon, mat.ylat, zlevel, msksp);
			}
		}
	}
	int nfiles = 0;
	if (nb_fishery && !param.flag_no_fishing){
		//Generate the name of the files
		for (int sp =0; sp<nb_species;sp++){
			dymFileSpC.push_back(param.strdir_output+param.sp_name[sp]+"_Cobs.dym");
			dymFileSpC.push_back(param.strdir_output+param.sp_name[sp]+"_Cpred.dym");
			nfiles += 2;
			if ( param.write_all_fisheries_dym ){
				int k = 0;
				for (int f=0;f<nb_fishery; f++){
					if (param.mask_fishery_sp[sp][f]){
						//Inna 09.09: need to write prediction and observation file by fishery to be used for validation
						dymFileSpC.push_back(param.strdir_output + param.sp_name[sp]+"_Cobs" + param.list_fishery_name[f] +".dym");
						dymFileSpC.push_back(param.strdir_output + param.sp_name[sp]+"_Cpred" + param.list_fishery_name[f] +".dym");
						nfiles += 2;
						//dymFileSpC.push_back(param.strdir_output + param.sp_name[sp]+"_CPUEobs" + param.list_fishery_name[f] +".dym");
						//dymFileSpC.push_back(param.strdir_output + param.sp_name[sp]+"_CPUEpred" + param.list_fishery_name[f] +".dym");
						//nfiles += 4;
						k++;
					} 
				}
			}
		}
	}
	for (int nf=0; nf<nfiles;nf++){

		wbin_header(dymFileSpC[nf], param.idformat, param.idfunc, minval, maxval, 
			param.nlong, param.nlat, nb_mo,	zlevel[0], zlevel[nb_mo -1], 
			mat.xlon, mat.ylat, zlevel, msksp);
	}
}
///////////////////////////////////////////////////////////////////////////////
void CReadWrite::SaveSepodymFileDym(CParam& param, PMap& map, CMatrices& mat)
{	
	const int nlon = param.nlong;
	const int nlat = param.nlat;
	const int nb_species = param.get_nbspecies();
	const int nb_fishery = param.get_nbfishery();

	int ixv = 0;

	// SAVE SPECIES DISTRIBUTIONS
	for (int sp=0;sp<nb_species;sp++)	{
		SaveDymFile(map, mat, dymFileSpPred[ixv], mat.Hs, nlon, nlat); ixv++;
		SaveDymFile(map, mat, dymFileSpPred[ixv], mat.Ha, nlon, nlat); ixv++;
		SaveDymFile(map, mat, dymFileSpPred[ixv], mat.advection_x, nlon, nlat); ixv++;
		SaveDymFile(map, mat, dymFileSpPred[ixv], mat.advection_y, nlon, nlat); ixv++;
		SaveDymFile(map, mat, dymFileSpPred[ixv], mat.speed, nlon, nlat); ixv++;
		SaveDymFile(map, mat, dymFileSpPred[ixv], mat.diffusion_y, nlon, nlat); ixv++;
		//SaveDymFile(map, mat, dymFileSpPred[ixv], mask_catch, nlon, nlat); ixv++;
	//	SaveDymFile(map, mat, dymFileSpPred[ixv], EFR, nlon, nlat); ixv++;
		
		//larvae
		SaveDymFile(map, mat, dymFileSpPred[ixv], mat.larvae[sp], nlon, nlat); ixv++;

		//juvenile
		SaveDymFile(map, mat, dymFileSpPred[ixv], mat.juvenile[sp], nlon, nlat);ixv++;

		//young
		SaveDymFile(map, mat, dymFileSpPred[ixv],  mat.young[sp], nlon, nlat);ixv++;

		//recruit
		SaveDymFile(map, mat, dymFileSpPred[ixv],  mat.recruit[sp], nlon, nlat); ixv++;
//		SaveDymFile(map, mat, dymFileSpPred[ixv],  mat.Hj, nlon, nlat); ixv++;

		//adults
		SaveDymFile(map, mat, dymFileSpPred[ixv],  mat.adult[sp], nlon, nlat);ixv++;

		//total population
		SaveDymFile(map, mat, dymFileSpPred[ixv],  mat.total_pop[sp], nlon, nlat);ixv++;
	

		// PREDICTED AND OBSERVED CATCH and CPUE
		if (nb_fishery && !param.flag_no_fishing){
			int ixv = 0;
			SaveDymFile(map, mat, dymFileSpC[ixv],  mat.total_obs_catch[sp],  nlon, nlat); ixv++;
			SaveDymFile(map, mat, dymFileSpC[ixv],  mat.total_pred_catch[sp], nlon, nlat); ixv++;
			if ( param.write_all_fisheries_dym ){
				int k = 0;
				dmatrix cpue_obs,cpue_est;
				cpue_obs.allocate(0, nlon+1, 0, nlat+1); 
				cpue_est.allocate(0, nlon+1, 0, nlat+1);
				for (int f=0;f<nb_fishery; f++){
					if (param.mask_fishery_sp[sp][f]){
						//Catch by fishery
						SaveDymFile(map, mat, dymFileSpC[ixv],  mat.catch_obs[sp][k],  nlon, nlat); ixv++;
						SaveDymFile(map, mat, dymFileSpC[ixv],  mat.catch_est[sp][k],  nlon, nlat); ixv++;	
						//SaveDymFile(map, mat, dymFileSpC[ixv],  mat.effort[f],  nlon, nlat); ixv++;	
/*
						//CPUE by fishery
						cpue_obs.initialize(); cpue_obs = -1;
						cpue_est.initialize();
						for (int i=map.imin; i <= map.imax; i++){
							for (int j=map.jinf[i] ; j<=map.jsup[i] ; j++){
								if (map.carte[i][j] && mat.effort[f][i][j]){
	double pred = mat.catch_est(sp,k,i,j);
	double obs = mat.catch_obs(sp,k,i,j);
	double like = 0.0;
	if (pred>0){
		like = pred - obs*log(pred) + gammln(obs+1.0);
	}
//	dvector like_lf = LFlike(mat.LF_qtr_obs(sp),mat.C_N_sp_age_fishery[0],param.sp_a0_adult[sp],
//			param.sp_cohorts(sp),param.nb_region_sp_B[0],k,param.frq_like_param(0));
        

								//cpue_obs(i,j) = like;//mat.effort(f,i,j);
								cpue_obs(i,j) = mat.catch_obs(sp,k,i,j)/mat.effort(f,i,j);
								cpue_est(i,j) = mat.catch_est(sp,k,i,j)/mat.effort(f,i,j);
							}
							}
						}						
					SaveDymFile(map, mat, dymFileSpC[ixv],  cpue_obs,  nlon, nlat); ixv++;
					SaveDymFile(map, mat, dymFileSpC[ixv],  cpue_est,  nlon, nlat); ixv++;	
*/
					k++;
					} 
				}
			}
		}
	}
}

void CReadWrite::SaveDymFile(PMap& map, CMatrices& mat, string file, const dmatrix& data, const int nlon, const int nlat)
{
	bool FileMode = true; 
	double minval = 3.4e38;
	double maxval = -3.4e38;

	mat2d_NoBorder.shallow_copy(mat.mat2d_NoBorder);
	mat2d_NoBorder.initialize();

	for (int i=map.imin; i <= map.imax; i++){
		for (int j=map.jinf[i] ; j<=map.jsup[i] ; j++){
			if (map.carte[i][j]){
				mat2d_NoBorder(i-1,j-1) = data(i,j);
				if (minval>data(i,j)) minval = (DYM_OUTPUT_TYPE)data(i,j);
				if (maxval<data(i,j)) maxval = (DYM_OUTPUT_TYPE)data(i,j);
			}
		}
	}
//if (file=="run4804x30dx2deg/output/simout/Young_Habitat.dym") cout << minval << " " << maxval << endl;
	wbin_transpomat2d(file, mat2d_NoBorder,	nlon, nlat, FileMode);
	rwbin_minmax(file, minval, maxval);
}


