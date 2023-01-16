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


double Etot = 0.0;
double Ctot = 0.0;
int nrec_oceanmask_original = 0;

//Function reads MONTHLY fishing data and distributes it to the model time step upon reading
//Note, we don't duplicate the data at every model time step, which is <= 30days, but only
//scale the effort and catch with the factor deltaT/30. Then the use of this data in accounting
//for fishing mortality or the prediction of catch is straighforward
void CReadWrite::rtxt_fishery_data(CParam& param, const PMap& map, const int nbt, const int jday_spinup)
{


	//IMPORTANT: data for each fishery should be sorted by date!
	const int nb_fishery = param.get_nbfishery();
	const int nb_species = param.get_nbspecies();
	const int nbi = param.get_nbi();
	const int nbj = param.get_nbj();

	//IMPORTANT:
	//1. use of monthly data only
	//2. final date doesn't take into account forecast period
	const int yearini  = param.ndatini/10000;
	const int monthini = ( param.ndatini - (yearini * 10000) ) / 100 ;
	const int yearfin  = param.ndatfin/10000;
	const int monthfin = ( param.ndatfin - (yearfin * 10000) ) / 100 ;
	const int start_fishing = yearini*10000+monthini*100;
	const int end_fishing   = yearfin*10000+monthfin*100;
	const int nby = yearfin-yearini+1;

	double deltaX = (param.deltaX/60);
	double deltaY = (param.deltaY/60);

	nbsp_file.allocate(0,nb_fishery-1);
	nbsp_file.initialize();
	position.allocate(0,nb_fishery-1); 
	numrec.allocate(0,nb_fishery-1);

	//currently we will save all catch data in one file
	string filename = param.file_catch_data[0];

	all_rec = 0;
	int nb_fishery_file = 0;

	//this variable will be used to count number of records in the domain
	nrec_oceanmask = 0;
	fishery_reso.allocate(0,nb_fishery-1);
	fishery_reso.initialize(); 
	param.fishery_reso = fishery_reso;

	double SUM = 0.0;

	ifstream littxt(filename.c_str());
	if (littxt){

		cout << endl << "Reading fishing data from the file: " << endl << filename.c_str() << endl;
		littxt >> nb_fishery_file;
		ivector nb_recs_fishery(0,nb_fishery_file-1);
		nb_recs_fishery.initialize();

		for (int f=0; f<nb_fishery_file; f++)
			littxt >> nb_recs_fishery(f);

		all_rec = sum(nb_recs_fishery);

		int nobs_tmp = all_rec;
		frec_tmp = new fishery_record[nobs_tmp];

		int fishery, year, month, day;
		float lon, lat, freso, eff;
		string gear;

		for (int f=0; f<nb_fishery; f++){
			if (param.mask_fishery_sp[0][f]){
				position(f).allocate(0,nby-1,1,12);
				numrec(f).allocate(0,nby-1, 1,12);
				position(f).initialize();
				numrec(f).initialize();
				nbsp_file(f) = nb_species;
			}
		}

		string s;
		//skip one line
		getline(littxt,s,'\n');
		getline(littxt,s,'\n');

		dvector mean_fishery_cpue;
		ivector nrec_fishery;
		mean_fishery_cpue.allocate(0,nb_fishery-1);
		nrec_fishery.allocate(0,nb_fishery-1);
		mean_fishery_cpue.initialize();
		nrec_fishery.initialize();

		int rec = 0;
		bool DataExist = false;
		while (rec<all_rec) {

			dvector harv(0,nb_species-1);
			harv.initialize();

			//read the data in following order:
			//f yr mm dd gr lat lon res E C
			littxt >> fishery >> year >> month >> day >> gear >> lat >> lon >> freso >> eff >> harv;
			if (lon<0) lon += 360.0;
			

			std::ostringstream ostr;
			ostr << fishery;
			string fishery_name = gear+ostr.str();

			int f = 0;
			for (;f<nb_fishery; f++)
			    if (param.mask_fishery_sp[0][f]){
				if (fishery_name == param.list_fishery_name[f]){
					fishery_reso(f) = freso;

					int date = year*10000 + month*100;
					int yr = year-yearini; 

					if (date >= start_fishing && date <= end_fishing){
						DataExist = true;

						//convert LL E units according to this parameter
						eff *= param.eff_units_converter[f];

						//scale eff and harv from monthly (30 days) to model temporal resolution
						//cout << "Scale fishing effort by " << (double)param.deltaT / 30 << endl;
						eff  *= (double)param.deltaT / 30;
						harv *= (double)param.deltaT / 30;

						//const 
						int i = param.lontoi(lon);
						//const 
						int j = param.lattoj(lat);	

						if (i==0 || j==0) {
							cout << "Attn: fishing data index 0 for (lon,lat) = (" 
								<< lon << "," << lat << "), (i,j) = (" << i << "," 
								<< j << ")" <<endl;
						}
						const int hr_x = (int)(freso/(2.0*deltaX));
						const int hr_y = (int)(freso/(2.0*deltaY));
						if ((i>hr_x && i+hr_x<nbi-1) && (j>hr_y && j+hr_y<nbj-1)){
							if(map.carte(i,j) && eff){
								SUM += harv[0];
	
								mean_fishery_cpue(f) += harv(0)/eff;
								nrec_fishery(f)++;
								frec_tmp[nrec_oceanmask].set_record(lon,lat,i,j,eff,harv[0]);
								
								if (numrec(f,yr,month)==0)
									position(f,yr,month) = nrec_oceanmask;
								nrec_oceanmask++;
								numrec(f,yr,month)++;
							}
						}
					}
				}
			}
			rec++;
		}		   
		if (! DataExist && sum(param.mask_fishery_sp(0))!=0){
			param.mask_fishery_sp.initialize();
			cout << "WARNING: no fishery data found: fishery mask forced to zero" << endl;
		}

		littxt.close();

		for (int f=0; f<nb_fishery; f++)
			mean_fishery_cpue(f) /= nrec_fishery(f);

		for (int f=0; f<nb_fishery; f++){
			double fact = 1.0;
			if (param.sp_name[0].find("bet")==0) fact = 0.4;
			if (param.sp_name[0].find("alb")==0) fact = 25.0;
			param.cpue_mult(f) = fact*max(mean_fishery_cpue)/mean_fishery_cpue(f);
		}
		cout << "Mean cpue by fishery: "<< mean_fishery_cpue << endl;
		cout << "weights for fisheries: "<< param.cpue_mult << endl;
	} else {
		cout << endl << "WARNING : Cannot read file " << filename.c_str() << ". Will exit now!" << endl;
		exit(1);
	}
	cout << "Number of cells with fishing data: " << nrec_oceanmask << ", sumt(Ct): " << SUM << endl;
	nrec_oceanmask_original = nrec_oceanmask;

	//reallocate the memory for this class
	frec = new fishery_record[nrec_oceanmask];
	for (int rec=0; rec<nrec_oceanmask; rec++)
		frec[rec] = frec_tmp[rec];

	delete [] frec_tmp;

	param.fishery_reso = fishery_reso;

	cout << param.fishery_reso << endl;

	//Redistribute fishing data on model resolution if 
	//1. There are 'no effort' fisheries. In this case
	//only these fisheries data will be redistributed
	//and stored on model resolution
	if (sum(param.mask_fishery_sp_no_effort))
		set_frec_rm_no_effort_fisheries(param,map,nbt,jday_spinup);		

	//2.
	if (param.mpa_simulation || param.nb_EEZ){
		set_frec_rm(param,map,nbt,jday_spinup);		
		param.fdata_rm = 1;
		ofstream wtxt;
		string file_out = param.strdir_output+"mpa_stats.txt";
		wtxt.open(file_out.c_str(), ios::out);
		
		if (wtxt){
		        wtxt << "year"<<'\t' << "month" << '\t' << "C.tot.wcpfc" << '\t'<< 
				"C.fselected.wcpfc" << '\t' << "C.fselected.mpa" << endl;
		}
		wtxt.close();	
	}		
}

void CReadWrite::degrade_fishery_reso(CParam& param, PMap& map, const int nbt, const int jday_spinup)
{
	const int nb_fishery = param.get_nbfishery();
	double deltaX = (param.deltaX/60);
	double deltaY = (param.deltaY/60);
	const int yearini  = param.ndatini/10000;

	//precompute the vector of lon,lat coords for degraded resolution
	float creso = param.catch_reso;
	int nflon = (int) ((param.longitudeMax - param.longitudeMin + deltaX)/creso+.5);
	int nflat = (int) ((param.latitudeMax - param.latitudeMin + deltaY)/creso+.5);
	flon.allocate(0,nflon-1); flon.initialize();
	flat.allocate(0,nflat-1); flat.initialize();
	float cwest = param.longitudeMin+0.5*deltaX+0.5;
	float cnorth= param.latitudeMax-0.5*deltaY+0.5;
	for (int i=0; i<nflon; i++) flon[i] = cwest + creso*(i+0.5);
	for (int j=0; j<nflat; j++) flat[j] = cnorth - creso*(j+0.5); 

	//Interested to see only in optimisation mode, which does not support currently high resolutions
	if (deltaX>=1) cout << "Centers of degraded C-cells: " << flon << endl << flat << endl; //exit(1);

	int year, month, day, jday, xx;
	for (int f=0; f<nb_fishery; f++){
	   int count = 0;
           if (!param.mask_fishery_sp_no_effort[0][f]){
		if (param.mask_fishery_sp[0][f]){
			if (fishery_reso[f] < creso) { // we want only to degrade resolution!
			//note: if length(data)=0, i.e. fishery_reso[f]==0 code below will be executed but 'p' loop
				fishery_reso[f] = creso;
				int m_old = 0;
				for (int tcount=1; tcount<=nbt; tcount++){
					//we need year and month only
					Date::update_time_variables(tcount, param.deltaT, param.date_mode, jday_spinup, jday, day, month, year, xx);	
					int y = year - yearini;
					int m = month;
					if (m != m_old){
						int nstart = position(f,y,m);
						int nend = nstart+numrec(f,y,m);
						m_old = m;
						for (int p=nstart; p<nend; p++){
							double lon = frec[p].get_efflon();
							double lat = frec[p].get_efflat();
				
							int fi = (int)((lon-(cwest+0.5*creso))/creso);
							int fj = (int)((cnorth-0.5*creso-lat)/creso);
							int i = param.lontoi(flon[fi]);
							int j = param.lattoj(flat[fj]);

							//NOTICE: some records may be lost if the C-cell center is on land
							if (map.carte(i,j)==0) { i = -1;count++;}
							frec[p].change_coord(flon[fi],flat[fj],i,j);
						}
					}
				}
			}
		}
            }
	}
	param.fishery_reso = fishery_reso;
	cout << "After degradation the fishing data resolutions (by fishery) are: " << 
		endl  << param.fishery_reso << endl;
}

void CReadWrite::delete_fisheries_rec() 
{ 
	if (fisheries_data_read)
		delete [] frec; 
}

// the fishing effort redistributed to the resolution of the model is used to compute fishing mortality and in MPA simulations
void CReadWrite::set_effort_rm(CParam& param, PMap& map, const int nbt, const int jday_spinup)
{
	const int nb_fishery = param.get_nbfishery();
	const int nb_species = param.get_nbspecies();
	const int nbi = param.get_nbi();
	const int nbj = param.get_nbj();
	const int yearini  = param.ndatini/10000;
	const int yearfin  = param.ndatfin/10000;
	const int nby = yearfin-yearini+1;

	//this mask will be used in comparison with MFCL estimations
        mask_catch.allocate(0,nbi-1,0,nbj-1);
	mask_catch.initialize(); 

	double deltaX = (param.deltaX/60);
	double deltaY = (param.deltaY/60);

	position_rm.allocate(0,nb_fishery-1); 
	numrec_rm.allocate(0,nb_fishery-1);
	for (int f=0; f<nb_fishery; f++){
		if (param.mask_fishery_sp[0][f]){
			position_rm(f).allocate(0,nby-1,1,12);
			numrec_rm(f).allocate(0,nby-1, 1,12);
			position_rm(f).initialize();
			numrec_rm(f).initialize();
		}
	}

	int af_m = max(fishery_reso)*60.0/param.deltaX; 
	int af_n = max(fishery_reso)*60.0/param.deltaY;  
	if (param.deltaX>60.0) af_m += 2;
	if (param.deltaY>60.0) af_n += 2;

	long int maxn = nrec_oceanmask;//corrected on 20210719, double check that it works in all cases!
	if (!param.fdata_rm)
        	maxn = nrec_oceanmask_original*af_m*af_n;
	if (maxn<0){ cerr << "Wrong size of memory to allocate for fishing effort at model resolution, nrec_oceanmask is too big:" << nrec_oceanmask_original << "; Will exit now!" << endl; exit(1);}


	//temporal containers as we don't know here the size of redistributed effort
	//fishing_effort *efr_tmp;
	efr_tmp = new fishing_effort[maxn];
	dvector C;
	C.allocate(0,nb_species-1);
	dmatrix af;
	nrec_oceanmask = 0;
	int year, month, day, jday, xx;

	for (int f=0; f<nb_fishery; f++)
	    if (param.mask_fishery_sp[0][f]){
		af_m = fishery_reso(f)*60.0/param.deltaX+2; 
		af_n = fishery_reso(f)*60.0/param.deltaY+2;  
		af.allocate(0,af_m-1,0,af_n-1);
		af.initialize();

		//same as in set_frec_rm: rearrange fishing data 
		//classes once a month, which is the time reso of fishing data
		int m_old = 0;
		for (int tcount=1; tcount<=nbt; tcount++){
			//we need year and month only
			Date::update_time_variables(tcount, param.deltaT, param.date_mode, jday_spinup, jday, day, month, year, xx);	
			int y = year - yearini;
			int m = month;
			if (m != m_old){
				int nstart = position(f,y,m);
				int nend = nstart+numrec(f,y,m);
				double Elon, Elat, E;
				m_old = m;

				for (int p=nstart; p<nend; p++){
					Elon = frec[p].get_efflon();
					Elat = frec[p].get_efflat();
					E = frec[p].get_effort();

					C = frec[p].get_catch();
					double sumE = 0;
					// if original resolution > model resolution
					if (deltaX<fishery_reso(f) || deltaY<fishery_reso(f)){
						int ki = 0; int kj = 0;
						//this function returns area of fishing within one cell
						param.afcoef(Elon,Elat, af, ki, kj, f);
			
						//move the effort to the ocean if some cells within its area are on land
						double afw = 0.0;
						for (int ii=0; ii<af_m; ii++)
							for (int jj=0; jj<af_n; jj++)
								if (af(ii,jj) && map.carte(ii+ki,jj+kj)) afw += af(ii,jj);
						afw = 1.0/afw;
						for (int ii=0; ii<af_m; ii++)
							for (int jj=0; jj<af_n; jj++){
								if (af(ii,jj) && map.carte(ii+ki,jj+kj)){
									double eff = afw*af(ii,jj)*E;	
									if (C[0]>0)
										mask_catch(ii+ki,jj+kj) += afw*af(ii,jj);//*C[0];
									efr_tmp[nrec_oceanmask].set_effort(ii+ki,jj+kj,eff);

									if (numrec_rm(f,y,m)==0) //position to start reading
										position_rm(f,y,m) = nrec_oceanmask;
									nrec_oceanmask++;
									numrec_rm(f,y,m)++;

									sumE += eff; 
								}
							}
					} else { // if original resolution < model resolution
						const int i = param.lontoi(Elon);
						const int j = param.lattoj(Elat);
						mask_catch(i,j) ++;//= C[0];
						efr_tmp[nrec_oceanmask].set_effort(i,j,E);
						if (numrec_rm(f,y,m)==0) 
							position_rm(f,y,m) = nrec_oceanmask;
						nrec_oceanmask++;
						numrec_rm(f,y,m)++;
						sumE+= E;
					}
					if (E-sumE>1e-2) {
						cout << "E was not conserved! " << Elon << " " << Elat << 
						"E before redistribution = " << E << "!= " << sumE << "after!" <<  endl; 
					}
				}//end of reading E for a given y and m
			}
		}
	}
	if (!param.use_mask_catch)
		mask_catch = 1e3;
	//reallocate the memory for this class knowing
	//the true number of records at model resolution
	efr = new fishing_effort[nrec_oceanmask];
	for (int rec=0; rec<nrec_oceanmask; rec++)
		efr[rec] = efr_tmp[rec];
	//deleting temporal container
	delete [] efr_tmp;
}


//the catch redistribution (optional) is necessary for MPA simulations || EEZ/MFCL regional extraction only
//i.e. if (deltaX<fishery_reso(f) | deltaY<fishery_reso(f))
void CReadWrite::set_frec_rm(CParam& param, const PMap& map, const int nbt, const int jday_spinup)
{
	cout << "Notification: EEZ or MPA option activated..." << endl;
	cout << "Begin redistribution of the fishing data, resolutions (by fishery) are:" << endl 
		<< param.fishery_reso << endl;

	const int nb_fishery = param.get_nbfishery();
	const int nb_species = param.get_nbspecies();
	const int nbi = param.get_nbi();
	const int nbj = param.get_nbj();
	const int yearini  = param.ndatini/10000;
	
	double deltaX = (param.deltaX/60);
	double deltaY = (param.deltaY/60);


	int af_m = max(fishery_reso)*60.0/param.deltaX+2; 
	int af_n = max(fishery_reso)*60.0/param.deltaY+2;  
	int maxn = nrec_oceanmask*af_m*af_n;

	//to account for fisheries redistributed to model reso
	ivector reso_model(0,nb_fishery-1);
	reso_model.initialize();

	frec_tmp = new fishery_record[maxn];

	nrec_oceanmask = 0;

	dmatrix Elon, Elat, effort;
	d3_array Cobs;

	Elon.allocate(0, nbi-1, 0, nbj-1); Elon.initialize();
	Elat.allocate(0, nbi-1, 0, nbj-1); Elat.initialize();
	effort.allocate(0, nbi-1, 0, nbj-1); 
	Cobs.allocate(0,nb_species-1);
	for (int sp=0; sp<nb_species; sp++)
		Cobs(sp).allocate(0, nbi-1, 0, nbj-1);

	for (int f=0; f<nb_fishery; f++)
	    if (param.mask_fishery_sp[0][f]){
		af_m = fishery_reso(f)*60.0/param.deltaX+2; 
		af_n = fishery_reso(f)*60.0/param.deltaY+2;  
		dmatrix af(0,af_m-1,0,af_n-1);
		dvector harv(0,nb_species-1);

		int year, month, day, jday, xx;
		int m_old = 0;
		for (int tcount=1; tcount<=nbt; tcount++){
			//we need year and month only
			Date::update_time_variables(tcount, param.deltaT, param.date_mode, jday_spinup, jday, day, month, year, xx);	
			int y = year - yearini;
			int m = month;
			//rearrange frec classes only once a month - reso of catch data
			if (m != m_old){
				int nstart = position(f,y,m);
				int nend = nstart+numrec(f,y,m);
				effort.initialize();
				Cobs.initialize();
				m_old = m;
				for (int p=nstart; p<nend; p++){
					int i = frec[p].get_i();
					int j = frec[p].get_j();
					Elon(i,j) = frec[p].get_efflon();
					Elat(i,j) = frec[p].get_efflat();
					effort(i,j) += frec[p].get_effort();
					//to be implemented later: 
					//int nbsp = nbsp_file[f];	
					harv = frec[p].get_catch(); 
					for (int sp=0; sp<nb_species; sp++)
						Cobs(sp,i,j) += harv[sp];
				}	
				double sumE = 0;
				double sumC = 0;
				position(f,y,m) = 0;
				numrec(f,y,m) = 0;

				for (int i=0; i<nbi-1; i++){
				for (int j=0; j<nbj-1; j++){
					if (effort(i,j)){
						// if original resolution > model resolution
						if (deltaX<fishery_reso(f) || deltaY<fishery_reso(f)){
							int ki = 0; int kj = 0;
							//this function returns area of fishing within one cell
							param.afcoef(Elon(i,j),Elat(i,j), af, ki, kj, f);
			
							//move the effort to the ocean if some cells within its area are on land
							double afw = 0.0;
							for (int ii=0; ii<af_m; ii++)
								for (int jj=0; jj<af_n; jj++)
									if (af(ii,jj) && map.carte(ii+ki,jj+kj)) afw += af(ii,jj);
							
							afw = 1.0/afw;
	
							for (int ii=0; ii<af_m; ii++){
								for (int jj=0; jj<af_n; jj++){
									if (af(ii,jj) && map.carte(ii+ki,jj+kj)){
										double rmf = afw*af(ii,jj);
										double eff = rmf*effort(i,j);	
										double harv = rmf*Cobs(0,i,j);
										//dvector harv(0,nb_species-1);
										//for (int sp=0; sp<nb_species; sp++)
										//	harv[sp] = rmf*Cobs(sp,i,j);	
									
										const int iilon	= param.itolon(ii+ki);			
										const int jjlat	= param.jtolat(jj+kj);	
										frec_tmp[nrec_oceanmask].set_record(iilon,jjlat,ii+ki,jj+kj,eff,harv);
		
										if (numrec(f,y,m)==0) //position to start reading
											position(f,y,m) = nrec_oceanmask;
										nrec_oceanmask++;
										numrec(f,y,m)++;
										sumE += eff; 
										sumC += harv; 
										//sumC += harv(0); 
									}
								}
							}
							reso_model(f) = 1;
						} else { // if original resolution < model resolution
							//dvector harv(0,nb_species-1);
							//for (int sp=0; sp<nb_species; sp++)
							//	harv[sp] = Cobs(sp,i,j);
							double harv = Cobs(0,i,j);
							frec_tmp[nrec_oceanmask].set_record(Elon(i,j),Elat(i,j),i,j,effort(i,j),harv);
	
							if (numrec(f,y,m)==0) 
								position(f,y,m) = nrec_oceanmask;
							nrec_oceanmask++;
							numrec(f,y,m)++;
							sumE+= effort(i,j);
							sumC+= harv;
							//sumC+= harv(0);
						}
					}
				}}
				if (sum(effort)-sumE>1e-2 ) {
					cout << "WARNING: for the fishery " << f << " "<< 
					" E was not conserved due to redistribution! " << 
					"total E before redistribution = " << sum(effort) << 
					"!= " << sumE << "after!" <<  endl; exit(1);
				}
				if (sum(Cobs[0])-sumC>1e-2 ) {
					cout << "WARNING: for the fishery " << f << " "<< 
					" C was not conserved due to redistribution! " << 
					"total C before redistribution = " << sum(Cobs[0]) << 
					"!= " << sumC << "after!" <<  endl; exit(1);
				}
			}
		}
	    }

	for (int f=0;f<nb_fishery; f++){
		if (reso_model(f)) 
			fishery_reso(f) = max(param.deltaX/60.0,param.deltaY/60.0);
	}
	param.fishery_reso = fishery_reso;

	cout << "Finished redistribution of the fishing data to the resolutions:" << endl 
		<< param.fishery_reso << endl;

	//clean-up the old class, reallocate and rewrite
	delete [] frec;
	frec = new fishery_record[nrec_oceanmask];
	for (int rec=0; rec<nrec_oceanmask; rec++)
		frec[rec] = frec_tmp[rec];
	
	delete [] frec_tmp;	
}

void CReadWrite::set_frec_rm_no_effort_fisheries(CParam& param, const PMap& map, const int nbt, const int jday_spinup)
{
	cout << "Notification: The redistribution of fisheries 'without effort' to the model resolution..." << endl;

	const int nb_fishery = param.get_nbfishery();
	const int nb_species = param.get_nbspecies();
	const int nbi = param.get_nbi();
	const int nbj = param.get_nbj();
	double deltaX = (param.deltaX/60);
	double deltaY = (param.deltaY/60);
	const int yearini  = param.ndatini/10000;


	int af_m = max(fishery_reso)*60.0/param.deltaX+2; 
	int af_n = max(fishery_reso)*60.0/param.deltaY+2;  
	int maxn = nrec_oceanmask*af_m*af_n;

	//to account for fisheries redistributed to model reso
	ivector reso_model(0,nb_fishery-1);
	reso_model.initialize();

	frec_tmp = new fishery_record[maxn];

	nrec_oceanmask = 0;

	dmatrix Elon, Elat, effort;
	d3_array Cobs;

	Elon.allocate(0, nbi-1, 0, nbj-1); Elon.initialize();
	Elat.allocate(0, nbi-1, 0, nbj-1); Elat.initialize();
	effort.allocate(0, nbi-1, 0, nbj-1); 
	Cobs.allocate(0,nb_species-1);
	for (int sp=0; sp<nb_species; sp++)
		Cobs(sp).allocate(0, nbi-1, 0, nbj-1);

	int year, month, day, jday, xx;

	for (int f=0; f<nb_fishery; f++)
	    if (param.mask_fishery_sp[0][f]){
		af_m = fishery_reso(f)*60.0/param.deltaX+2; 
		af_n = fishery_reso(f)*60.0/param.deltaY+2;  
		dmatrix af(0,af_m-1,0,af_n-1);
		int m_old = 0;
		for (int tcount=1; tcount<=nbt; tcount++){
			Date::update_time_variables(tcount, param.deltaT, param.date_mode, jday_spinup, jday, day, month, year, xx);	

			int y = year - yearini;
			int m = month;
			if (m != m_old){
				int nstart = position(f,y,m);
				int nend = nstart+numrec(f,y,m);
				m_old = m;
				effort.initialize();
				Cobs.initialize();
				for (int p=nstart; p<nend; p++){
					int i = frec[p].get_i();
					int j = frec[p].get_j();
					Elon(i,j) = frec[p].get_efflon();
					Elat(i,j) = frec[p].get_efflat();
					effort(i,j) += frec[p].get_effort();
					Cobs(0,i,j) += frec[p].get_catch();
					//int nbsp = nbsp_file[f];	
					//dvector harv(0,nbsp-1);
					//harv = frec[p].get_catch(); 
					//for (int sp=0; sp<nb_species; sp++)
					//	Cobs(sp,i,j) += harv[sp];
				}	
				
				double sumE = 0;
				double sumC = 0;
				position(f,y,m) = 0;
				numrec(f,y,m) = 0;
				for (int i=0; i<nbi-1; i++){
				for (int j=0; j<nbj-1; j++){
					if (effort(i,j)){
						// if original resolution > model resolution
						if ((deltaX<fishery_reso(f) || deltaY<fishery_reso(f)) 
								&& param.mask_fishery_sp_no_effort[0][f]){//Attn: '0' for species!
							int ki = 0; int kj = 0;
							//this function returns area of fishing within one cell
							param.afcoef(Elon(i,j),Elat(i,j), af, ki, kj, f);
				
							//move the effort to the ocean if some cells within its area are on land
							double afw = 0.0;
							for (int ii=0; ii<af_m; ii++)
								for (int jj=0; jj<af_n; jj++)
									if (af(ii,jj) && map.carte(ii+ki,jj+kj)) afw += af(ii,jj);
							afw = 1.0/afw;
		
							for (int ii=0; ii<af_m; ii++){
							for (int jj=0; jj<af_n; jj++){
								if (af(ii,jj) && map.carte(ii+ki,jj+kj)){
									double rmf = afw*af(ii,jj);
									double eff = rmf*effort(i,j);	
									//dvector harv(0,nb_species-1);
									//for (int sp=0; sp<nb_species; sp++)
									//	harv[sp] = rmf*Cobs(sp,i,j);	
									double harv = rmf*Cobs(0,i,j);	
									
									const int iilon	= param.itolon(ii+ki);			
									const int jjlat	= param.jtolat(jj+kj);		
									//frec_tmp[nrec_oceanmask].create_catch(nb_species);	
									frec_tmp[nrec_oceanmask].set_record(iilon,jjlat,ii+ki,jj+kj,eff,harv);

									if (numrec(f,y,m)==0) //position to start reading
										position(f,y,m) = nrec_oceanmask;
									nrec_oceanmask++;
									numrec(f,y,m)++;
									sumE += eff; 
									sumC += harv; 
									//sumC += harv(0); 
								}
							}}
							reso_model(f) = 1;
						} else { // if original resolution < model resolution or 'normal' fishery
							//dvector harv(0,nb_species-1);
							//for (int sp=0; sp<nb_species; sp++)
							//	harv[sp] = Cobs(sp,i,j);
							double harv = Cobs(0,i,j);
							//frec_tmp[nrec_oceanmask].create_catch(nb_species);	
							frec_tmp[nrec_oceanmask].set_record(Elon(i,j),Elat(i,j),i,j,effort(i,j),harv);
	
							if (numrec(f,y,m)==0) 
								position(f,y,m) = nrec_oceanmask;
							nrec_oceanmask++;
							numrec(f,y,m)++;
							sumE+= effort(i,j);
							sumC+= harv;
							//sumC+= harv(0);
						}
					}
				}}
				if (sum(effort)-sumE>1e-2 ) {
					cout << "WARNING: for the fishery " << f << " "<< 
						" E was not conserved due to redistribution! " << 
						"total E before redistribution = " << sum(effort) << 
						"!= " << sumE << "after!" <<  endl; exit(1);
				}
				if (sum(Cobs[0])-sumC>1e-2 ) {
					cout << "WARNING: for the fishery " << f << " "<< 
						" C was not conserved due to redistribution! " << 
						"total C before redistribution = " << sum(Cobs[0]) << 
						"!= " << sumC << "after!" <<  endl; exit(1);
				}
			}
		}
    	    }

	for (int f=0;f<nb_fishery; f++){
		if (reso_model(f)) 
			fishery_reso(f) = max(param.deltaX/60.0,param.deltaY/60.0);
	}
	param.fishery_reso = fishery_reso;

	cout << "Finished redistribution of the fishing data to the resolutions:" << endl 
		<< param.fishery_reso << endl;

	//clean-up the old class, reallocate and rewrite
	delete [] frec;
	frec = new fishery_record[nrec_oceanmask];
	for (int rec=0; rec<nrec_oceanmask; rec++)
		frec[rec] = frec_tmp[rec];
	
	delete [] frec_tmp;	
}


void CReadWrite::get_fishery_data(CParam& param, D3_ARRAY& effort, D4_ARRAY& catch_obs, D3_ARRAY& efflon, D3_ARRAY& efflat, int y, const int m)
{	
	y -= (int)param.save_first_yr;
	int k = 0;
	const int nb_fishery = param.get_nbfishery();
	for (int f=0; f<nb_fishery; f++){
		effort(f).initialize();
		efflon(f).initialize();
		efflat(f).initialize();
		if (param.mask_fishery_sp[0][f]){
			catch_obs(0,k).initialize();
			int nstart = position(f,y,m);
			int nend = nstart+numrec(f,y,m); 
			for (int p=nstart; p<nend; p++){
				int i = frec[p].get_i(); 
				if (i<0) continue;
				int j = frec[p].get_j(); 
				effort(f,i,j) += frec[p].get_effort();
				efflon(f,i,j) = frec[p].get_efflon();
				efflat(f,i,j) = frec[p].get_efflat();
				catch_obs(0,k,i,j) += frec[p].get_catch();
				//int nb_species = nbsp_file[f];	
				//dvector harv(0,nb_species-1);
				//harv = frec[p].get_catch(); 
				//for (int sp=0; sp<nb_species; sp++){
				//	catch_obs(sp,k,i,j) += harv[sp];
				//}	
			} k++;
		}
	}
}

void CReadWrite::mpa_areas_comp(PMap& map, CParam& param)
{
	double area_mpa = 0.0;
	double area_els = 0.0;

	for (int i = map.imin; i <= map.imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			if (map.carte(i,j)){
				if (map.maskMPA(i,j)<0) area_mpa += param.cell_surface_area(j);
	                        if (map.maskMPA(i,j)>0) area_els += param.cell_surface_area(j);

	                        //compute areas of MPAs and WCPFC Convention Area:
                                //Northern part
                                if (param.itolon(i)>=115 && param.itolon(i)<=210 && 
                                    param.jtolat(j)>=-5 && param.jtolat(j)<=45){
                                        if (map.maskMPA(i,j)<0) area_mpa += param.cell_surface_area(j);
                                        if (map.maskMPA(i,j)>0) area_els += param.cell_surface_area(j);

                                }
                                //Southern part
                                if (param.itolon(i)>=115 && param.itolon(i)<=230 && 
                                    param.jtolat(j)>=-45 && param.jtolat(j)<-5){
                                        if (map.maskMPA(i,j)<0) area_mpa += param.cell_surface_area(j);
                                        if (map.maskMPA(i,j)>0) area_els += param.cell_surface_area(j);
                                }
			}
		}
	}
	cout << endl<<"Areas: mpa = " << area_mpa << 
		", total = " << area_mpa+area_els << 
		", percentage of mpa = " << 100*area_mpa/(area_els+area_mpa) << endl;
}

void CReadWrite::get_fishery_data_mpa(PMap& map, CParam& param, D3_ARRAY& effort, D4_ARRAY& catch_obs, D3_ARRAY& efflon, D3_ARRAY& efflat, int y, const int m)
{
        //WARNING: mpa scenarios work for 'sp=0' only
	//mpa mask should be flagged with negative values for MPA and positive for nonMPA zones
        y -= (int)param.save_first_yr;
        int nb_mpa = param.nb_mpa;
        int k = 0;
        const int nb_fishery = param.get_nbfishery();

	//Some statistics
	double Ctot = 0.0;
	double Ctot_mpa_fishery = 0.0;
	double Ctot_mpa = 0.0;


	dvector sum_mpa(0,nb_mpa-1);
	dvector nb_cells_mpa(0,nb_mpa-1);
	dvector sumE1(0,nb_mpa-1);
	double sum_nonmpa, sum_cpue_nonmpa;

	nb_cells_mpa.initialize();
	if (!param.actual_eff){
	  for (int i = map.imin; i <= map.imax; i++){
            const int jmin = map.jinf[i];
	    const int jmax = map.jsup[i];
	    for (int j = jmin; j <= jmax; j++){
	      if (map.carte(i,j)){
		for (int n=0; n<nb_mpa; n++)
		  if (map.maskMPA(i,j)==param.mpa_ID[n])
			nb_cells_mpa[n]++;
	      }
	    }
	  }
	}
	
        for (int f=0; f<nb_fishery; f++){

                effort(f).initialize();

                if (param.mask_fishery_sp[0][f]){

                        catch_obs(0,k).initialize();

                        int nstart = position(f,y,m);
                        int nend = nstart+numrec(f,y,m);

			sum_mpa.initialize();
			sum_nonmpa = 0.0;
			sum_cpue_nonmpa = 0.0;
                        int count = 0;

			sumE1.initialize();
			double sumE2 = 0;

                        for (int p=nstart; p<nend; p++){
                                int i = frec[p].get_i();
                                int j = frec[p].get_j();
                                if (!map.carte(i,j)) continue;

                                effort(f,i,j) += frec[p].get_effort();
				double harv = frec[p].get_catch();
				catch_obs(0,k,i,j) += harv;
                                //dvector harv(0,nb_species-1);
                                //harv = frec[p].get_catch();
                                //for (int sp=0; sp<nb_species; sp++){
                                //        catch_obs(sp,k,i,j) += harv[sp];
                                //}
                                
				//compute Ctot in WCPFC
                                if ((param.itolon(i)>=115 && param.itolon(i)<=210 &&
                                    param.jtolat(j)>=-5 && param.jtolat(j)<=45)
				    || (param.itolon(i)>=115 && param.itolon(i)<=230 &&
                                    param.jtolat(j)>=-45 && param.jtolat(j)<-5)){
                                	Ctot += harv;
                                	//Ctot += harv[0];
	                                if (param.mpa_fishery[f]){
						Ctot_mpa_fishery += harv;
						//Ctot_mpa_fishery += harv[0];
						for (int n=0; n<nb_mpa; n++)
	        	                       		if (map.maskMPA(i,j)==param.mpa_ID[n])
								Ctot_mpa += harv;
								//Ctot_mpa += harv[0];
					}
				}

                                //sum_mpa - total effort inside MPAs
                                //sum_nonmpa - total effort outside of MPAs; 
                                //sum_cpue_nonmpa - total CPUE outside of MPAs
                                //count- nb of cells visited by fishery outside of MPAs
				//nb_cells_mpa - number of cells within MPA
				if (param.mpa_fishery[f]){
					for (int n=0; n<nb_mpa; n++)
		                                if (map.maskMPA(i,j)==param.mpa_ID[n])
        		                               	sum_mpa[n] += frec[p].get_effort();

					if (map.maskMPA(i,j)>0){
						count ++;
						sum_nonmpa += frec[p].get_effort();
					}
				}
                        }

                        //need to compute sum(CPUE) after reading all effort and catch data:
                        //after redistribution their might be more than one record for ij
                        for (int i = map.imin; i <= map.imax; i++){
                          const int jmin = map.jinf[i];
                          const int jmax = map.jsup[i];
                          for (int j = jmin; j <= jmax; j++){
                            if (map.carte(i,j)){
			      if (param.mpa_fishery[f]){
	                        if (map.maskMPA(i,j)>0){
                              		double E = effort(f,i,j);
	                           	if (E>0)
	                                  sum_cpue_nonmpa += catch_obs(0,k,i,j)/E;
				}
			      }
                            }
			  }
                        }

                        double sumtot = sum(sum_mpa)+sum_nonmpa; //total effort at the time step
                        //nothing to do if
                        //1. sum = 0
                        //2. fishery is not considered in MPA scenarios
	                if (sumtot==0 || !param.mpa_fishery[f]){
	                	k++;
	                        continue;
	                }
			double CHECKSUM = 0.0;
			if (!param.actual_eff){
			//will redistribute the effort homogeneously
                          for (int i = map.imin; i <= map.imax; i++){
                            const int jmin = map.jinf[i];
                            const int jmax = map.jsup[i];
                            for (int j = jmin; j <= jmax; j++){
                              if (map.carte(i,j)){
			        for (int n=0; n<nb_mpa; n++){
				  if (map.maskMPA(i,j)==param.mpa_ID[n]){
					effort(f,i,j) = sum_mpa[n]/nb_cells_mpa[n];
					efflon(f,i,j) = param.itolon(i);
					efflat(f,i,j) = param.jtolat(j);	
					CHECKSUM += effort(f,i,j);
				  }
			        }
			      }
			    }
			  }
			  if (sqrt(pow(sum(sum_mpa)-CHECKSUM,2))>1e-3) cout << sum_mpa << " " << sum(sum_mpa) << " != "<< CHECKSUM << endl;
			}

			double X = 1.0;

                        for (int i = map.imin; i <= map.imax; i++){
                         const int jmin = map.jinf[i];
                         const int jmax = map.jsup[i];
                          for (int j = jmin; j <= jmax; j++){
                           if (map.carte(i,j)){
                            if (effort(f,i,j)){
			     for (int n=0; n<nb_mpa; n++){
	
				int S = param.mpa_scenario[n];

                                //scenarios with MPA closure:
                                if (map.maskMPA(i,j)==param.mpa_ID[n] && (S==2 || S==3 || S==4))
                                        effort(f,i,j) = 0.0;
                                switch (S) {
                                        case 1: //+X% in MPA    
						X = 1+param.mpa_S1_X[n]/100.0;
                                                if (map.maskMPA(i,j)==param.mpa_ID[n])
                                                        effort(f,i,j) *= X;
                                                break;
                                        case 2: //close MPA for fishing
                                                //see above, E is set =0 in 'negative' cells 
                                                break;
                                        case 3: //close MPA and redistribute effort equally outside of MPAs
                                                if (map.maskMPA(i,j)>0)
                                                        effort(f,i,j) += sum_mpa[n]/count;
                                                break;
                                        case 4: //close MPA and redistribute effort ~ to:
                                                if (map.maskMPA(i,j)>0){
                                                        //CPUE outside of MPAs
                                                        double CPUE = catch_obs(0,k,i,j)/effort(f,i,j);
                                                        if (sum_cpue_nonmpa>0) {//if total CPUE!=0
                                                                effort(f,i,j) += CPUE*sum_mpa[n]/sum_cpue_nonmpa;
							}
                                                        else if (sum_cpue_nonmpa==0){
                                                                //if sum(CPUE)=0, redistr ~ E
                                                                effort(f,i,j) *= (1+sum_mpa[n]/sum_nonmpa);                    
							}
                                                }
                                                break;
                                        case 5: //effort proportional reduction equivalent to MPA closure
                                                effort(f,i,j) *= 1-sum_mpa[n]/sumtot;
						break;
                                }//end of switch
                                if (map.maskMPA(i,j)==param.mpa_ID[n]) sumE1[n] += effort(f,i,j);
			     }//end of 'n' loop
                             if (map.maskMPA(i,j)>0) sumE2 += effort(f,i,j);
                        }}}}
                        //Control:
                        double beforsum = sumtot;
			for (int n=0; n<nb_mpa; n++){
                            double aftersum = sumE1[n]+sumE2;
			    int S = param.mpa_scenario[n];
                            switch(S) {
                                case 1:
					X = param.mpa_S1_X[n]/100.0;
                                        if (sqrt(pow(sum_nonmpa+(1+X)*sum_mpa[n]-aftersum,2))>1e-3) 
						cerr << "WARNING: An error in the scenario " << S << " for the fishery "<< f <<", "
                                                << sum_nonmpa+(1+X)*sum_mpa[n] << " != " << aftersum << endl;
                                	break;
                                case 2:
                                        if (sqrt(pow(aftersum-sum_nonmpa,2))>1e-3) 
						cerr << "WARNING: An error in the scenario " << S << " for the fishery "<< f <<", "
                                                << sum_nonmpa << " != " << aftersum << endl;
                                	break;
                                case 3:
                                        if (sqrt(pow(aftersum-beforsum,2))>1e-3) 
						cerr << "WARNING: total effort of the fishery "<< f <<" is not preserved in the scenario "
                                                << S << ", " << beforsum << " != " << aftersum << endl;
                                	break;
                                case 4:
                                        if (sqrt(pow(aftersum-beforsum,2))>1e-3) 
						cerr << "WARNING: total effort of the fishery "<< f <<" is not preserved in the scenario "
                                                << S << ", " << beforsum << " != " << aftersum << endl;
                                	break;
                                case 5:
                                        if (sqrt(pow(aftersum-sum_nonmpa,2))>1e-3) 
						cerr << "WARNING: total effort of the fishery "<< f <<" is not reduced by Empa in the scenario "
                                                << S << ", " << sum_nonmpa << " != " << aftersum << endl;
                               		break;

                            }
			}
                        k++;
                }
        }//end of 'f' loop
//to write some MPA statistics in the ASCII file:
ofstream wtxt;
string file_out = param.strdir_output+"mpa_stats.txt";
wtxt.open(file_out.c_str(), ios::app);
if (wtxt){
        wtxt << y+(int)param.save_first_yr << '\t' << m << '\t' << Ctot << '\t'<< Ctot_mpa_fishery << '\t' << Ctot_mpa << endl;
}
wtxt.close();
}

void CReadWrite::inc_obs_catch_mpa(PMap& map, CParam& param, dmatrix& catch_obs, const int sp)
{
        //Will multiply the observed catch in MPA 
	//for the fisheries declared without effort
	//temporarily using the effort multipliers
	int nb_mpa = param.nb_mpa;
        for (int i = map.imin; i <= map.imax; i++){
        	const int jmin = map.jinf[i];
                const int jmax = map.jsup[i];
                for (int j = jmin; j <= jmax; j++){
                	if (map.carte(i,j)){
				for (int n=0; n<nb_mpa; n++){
					if (map.maskMPA(i,j)==param.mpa_ID[n]){
						double X = 1+param.mpa_S1_X[n]/100.0;
						catch_obs(i,j) *= X;
					}
				}
			}
		}
	}
}


void CReadWrite::get_average_effort(CParam& param, D3_ARRAY& effort, D3_ARRAY& efflon, D3_ARRAY& efflat, const int nby, const int m)
{	
	int y1 = max(param.save_first_yr,param.save_last_yr-nby)-param.save_first_yr;
	int y2   = (int)param.save_last_yr-(int)param.save_first_yr; 

	int k = 0;
	const int nb_fishery = param.get_nbfishery();


	for (int f=0; f<nb_fishery; f++){
		effort(f).initialize();
		efflon(f).initialize();
		efflat(f).initialize();
		if (param.mask_fishery_sp[0][f]){
			//catch_obs(0,k).initialize();
			for (int y=y1; y<y2; y++){
				int nstart = position(f,y,m);
				int nend = nstart+numrec(f,y,m); 
				for (int p=nstart; p<nend; p++){
					int i = frec[p].get_i(); 
					if (i<0) continue;
					int j = frec[p].get_j(); 
					effort(f,i,j) += frec[p].get_effort();
					efflon(f,i,j) = frec[p].get_efflon();
					efflat(f,i,j) = frec[p].get_efflat();				
				} k++;
			}
			effort(f) /= (double)nby;	 
			//effort(f) = 2.0*effort(f);
		}
	}
}

void CReadWrite::get_average_effort_rm(CParam& param, dmatrix& effort, const int f, const int nby, const int m)
{	
	int y1 = max(param.save_first_yr,param.save_last_yr-nby+1)-param.save_first_yr;
	int y2   = (int)param.save_last_yr-(int)param.save_first_yr; 

	int k = 0; 
	effort.initialize();
	for (int y=y1; y<=y2; y++){
		int nstart = position_rm(f,y,m);
		int nend = nstart+numrec_rm(f,y,m); 
		for (int p=nstart; p<nend; p++){
			int i = efr[p].get_i(); 
			if (i<0) continue;
			int j = efr[p].get_j(); 
			effort(i,j) += efr[p].get_effort();
		} 
		k++;
	}
	effort /= (double)min(nby,k);	 
	//effort *= 2.0;
}

//this effort is given on original data resolution, is used to predict catch
void CReadWrite::get_catch(CParam& param, dmatrix& catch_obs, const int f, int y, const int m, const int sp)
{	
	catch_obs.initialize();
	y -= (int)param.save_first_yr;
	int nstart = position(f,y,m);
	int nend = nstart+numrec(f,y,m); 
	//dvector harv(0,nb_species-1);
	//harv.initialize();
	for (int p=nstart; p<nend; p++){
		int i = frec[p].get_i();
		int j = frec[p].get_j();
		catch_obs(i,j) += frec[p].get_catch();
		//harv = frec[p].get_catch(); 
		//catch_obs(i,j) += harv[sp];
	}
}


//this effort is given on original data resolution, is used to predict catch
void CReadWrite::get_effort(CParam& param, dmatrix& effort, const int f, int y, const int m)
{	
	effort.initialize();
	y -= (int)param.save_first_yr;
	int nstart = position(f,y,m);
	int nend = nstart+numrec(f,y,m); 
	for (int p=nstart; p<nend; p++){
		int i = frec[p].get_i();
		int j = frec[p].get_j();
		effort(i,j) += frec[p].get_effort();
	}
}

//this effort is given on original data resolution, is used to predict catch
void CReadWrite::get_effort_lonlat(CParam& param, dmatrix& effort, dmatrix& efflon, dmatrix& efflat, const int f, int y, const int m)
{	
	effort.initialize();
	efflon.initialize();
	efflat.initialize();
	y -= (int)param.save_first_yr;
	int nstart = position(f,y,m);
	int nend = nstart+numrec(f,y,m); 
	for (int p=nstart; p<nend; p++){
		int i = frec[p].get_i();
		if (i<0) continue;
		int j = frec[p].get_j();
		effort(i,j) += frec[p].get_effort();
		efflon(i,j) = frec[p].get_efflon();
		efflat(i,j) = frec[p].get_efflat();
	}
}

//redistiributed to the resolution of the model fishing effort is used to compute fishing mortality
void CReadWrite::get_effort_rm(CParam& param, dmatrix& effort, const int f, int y, const int m)
{
	effort.initialize();
	string fname = param.list_fishery_name[f];
	if (fname.find("Ac") == 0) return;
	if (param.nb_yr_forecast) 
		if (y>param.save_last_yr) {//temporal, no month
			get_average_effort_rm(param,effort,f,10,m); //temporal, 10 must be global
			return;
		}
	y -= (int)param.save_first_yr;
	int nstart = position_rm(f,y,m);
	int nend = nstart+numrec_rm(f,y,m); 
	for (int p=nstart; p<nend; p++){
		int i = efr[p].get_i();
		int j = efr[p].get_j();
		effort(i,j) += efr[p].get_effort();
	}
}

void CReadWrite::get_average_selectivity(PMap& map, CParam& param, dvector& swa, const ivector fisheries, const int nbf, const int nbt, const int nb_ages, const int sp, const int step_count)
{
        swa.initialize();
        const int nb_fishery = param.get_nbfishery();
        dvector Ctot;
        Ctot.allocate(0,nb_fishery-1);
        Ctot.initialize();
        int firstyear = (int)param.save_first_yr; int year = firstyear;
        int firstmonth = param.get_month(param.save_first_yr); int month = firstmonth;
        for (int tcount=0; tcount<nbt; tcount++){
                int y = year - firstyear;
                int m = month;
                for (int f=0; f<nb_fishery; f++){
                        if (fisheries(f)){
                                int nstart = position(f,y,m);
                                int nend = nstart+numrec(f,y,m);
                                for (int p=nstart; p<nend; p++){
                                        int i = frec[p].get_i();
                                        int j = frec[p].get_j();
					if (map.maskEEZ[i][j]<0)
	                                        Ctot(f) += frec[p].get_catch();
                                        //dvector harv(0,nb_species-1);
                                        //harv = frec[p].get_catch();
                                        //for (int sp=0; sp<nb_species; sp++)
                                        //        Ctot(f) += harv[sp];
                                }
                        }
                }
        }
        double Ctotf = sum(Ctot);
        int k = 0;
        dmatrix s(0, nb_fishery-1, 0, nb_ages-1);
        s.initialize();
        for (int f=0; f<nb_fishery; f++){
                if (param.mask_fishery_sp[0][f]){
                        if (fisheries(f)){
                                for (int a=0; a<nb_ages; a++)
                                        s(f,a) = param.selectivity_comp(sp,a,f,k);

                        }k++;
                }
//                s(f) /= max(s(f)); //scale so maximum=1
        }
        for (int a=0; a<nb_ages; a++){
                for (int f=0; f<nb_fishery; f++)
                        if (fisheries(f))
                                swa(a) += s(f,a)*Ctot(f)/Ctotf;
        }
        swa /= max(swa);
	cout << endl << "Total catches for fisheries used to compute exploited biomass = " << Ctot << endl;
	cout << "average selectivity = " << swa << endl;
        //To activate F=0 scenario
        //param.mask_fishery_sp.initialize();
}


void CReadWrite::read_lf_WCPO(CParam& param, string filename, const float startdate, const float enddate, const int sp)
{ //data used by Multifan-CL, file changed only in header

	const int nb_fishery = param.get_nbfishery();
	const int a0	     = param.sp_a0_adult[sp];
	const int nb_ages    = param.sp_nb_cohorts[sp];
	int nbsd = 0; //count the number of size distributions to be used

	ifstream littxt(filename.c_str());
	if (littxt){
		int nb_regions, nb_fleets, nb_records, nb_intervals, l1, dl;
		//float w1, wl;
		string* fcodes;

		littxt >> nb_regions >> nb_fleets >> nb_records;

		ivector region(0,nb_fleets-1);
		region.initialize();
		fcodes = Utilities::create1d(fcodes, nb_fleets);

		for (int f=0; f<nb_fleets; f++)
			littxt >> fcodes[f];
		for (int f=0; f<nb_fleets; f++)
			littxt >> region[f];
		//read the line 
		//flexible to different file contents, namely to the presence of weight data
		string s;
		getline(littxt,s,'\n');
		getline(littxt,s,'\n');
		std::istringstream iss(s);
		iss >> nb_intervals >> l1 >> dl;

		dvector lf(0,nb_intervals-1);
		lf.initialize();
		ivector Length(0,nb_intervals-1);
		for (int l=0; l< nb_intervals; l++)
			Length(l) = (int)(l1+dl*(l+0.5)); 
			//Length(l) = (int)(l1+dl*l); 

		//skip the line
		getline(littxt,s,'\n');
		getline(littxt,s,'\n');
		int year, month, week, f_num;
		float E, C;

		for (int n=0; n<nb_records; n++){
			littxt >> year >> month >> week >> f_num >> C >> E;
			float frequency;
			lf.initialize();
			//read the whole line (some files contain weight frequencies data)
			getline(littxt,s,'\n');
			std::istringstream iss(s);
			for (int l=0; l<nb_intervals; l++){
				iss >> frequency;
				if (frequency == -1) break;
				lf[l] = frequency;	
			}
			bool fishery_to_throw = true;
			int f;
			for (f=0;f<nb_fishery;f++){
				if (param.list_fishery_name[f] == fcodes[f_num-1]){
					fishery_to_throw = false; 
					break;
				}
				
			}
			if (fishery_to_throw) continue;

			float date = year+month/12.0;
			if ((date >= startdate) && (date <=enddate)){

				nbsd ++;

				int qtr = month/3;
				///double L_pr = param.juv_length(sp,param.sp_nb_age_class_jv[sp]-1);
				double L_pr = param.length(sp,a0-1);
				int y=year-(int)startdate;
				double right;
				for (int a=a0; a< nb_ages-1; a++){

					double left  = 0.5*(L_pr+param.length(sp,a));
					right = 0.5*(param.length(sp,a)+param.length(sp,a+1));
					L_pr = param.length(sp,a);

					for (int l=0; l<nb_intervals; l++){
						if ((Length[l]>left) && (Length[l]<=right))
							frq(f,region(f_num-1)-1,y,qtr,a) += lf[l];
					}
				}
				for (int l=0; l<nb_intervals; l++)
					if ((Length[l]>right)&& (Length[l]<param.length[sp][nb_ages-1]+10))
						frq(f,region(f_num-1)-1,y,qtr,nb_ages-1) += lf[l];
			}
		}
		littxt.close();
	} else {
		cout << endl << "WARNING : Cannot read file " << filename.c_str() << endl;
	}
	cout << "Number of size distributions used from LF data: " << nbsd << endl;
	cout << endl;
}

void CReadWrite::read_lf_EPO(CParam& param, string filename, const float startdate, const float enddate, const int sp)
{//data provided by IATTC, file changed only in header

	const int nb_fishery = param.get_nbfishery();
	///const int nb_ages    = param.sp_nb_age_class_ad[sp];
	const int a0	     = param.sp_a0_adult[sp];
	const int nb_ages    = param.sp_nb_cohorts[sp];
	int nbsd = 0; //count the number of size distributions to be used

	ifstream littxt(filename.c_str());
	if (littxt){
		int region, nb_fleets, nb_records, nb_intervals, l1, dl;
		string* fcodes;

		littxt >> region >> nb_fleets >> nb_records;
		fcodes = Utilities::create1d(fcodes, nb_fleets);

		for (int f=0; f<nb_fleets; f++)
			littxt >> fcodes[f];

		littxt >> nb_intervals >> l1 >> dl;
		dvector lf(0,nb_intervals-1);
		lf.initialize();
		dvector Length(0,nb_intervals-1);
		for (int l=0; l< nb_intervals; l++)
			Length(l) = l1+(float)dl*(l+0.5); 
		string s;
		getline(littxt,s,'\n');
		getline(littxt,s,'\n');

		int year, qtr, f_num;
		for (int n=0; n<nb_records; n++){
			littxt >> year >> qtr >> f_num >> s;

			int frequency;
			lf.initialize();

                        for (int l=0; l<nb_intervals; l++){
                        	littxt >> frequency;
                                lf[l] = frequency;
                        }

			bool fishery_to_throw = true;
			int f;
			for (f=0;f<nb_fishery;f++){
				if (param.list_fishery_name[f] == fcodes[f_num-1]){
					fishery_to_throw = false; 
					break;
				}
				
			}
			if (fishery_to_throw) continue;

			float date = year+(qtr*3.0-1)/12;
			if ((date >= startdate) && (date <=enddate)){
				nbsd++;

				double L_pr = param.length(sp,a0-1);
				int y=year-(int)startdate;
				double right = 0;
				for (int a=a0; a< nb_ages-1; a++){
					double left  = 0.5*(L_pr+param.length(sp,a));
					right = 0.5*(param.length(sp,a)+param.length(sp,a+1));
					L_pr = param.length(sp,a);

					for (int l=0; l<nb_intervals; l++){
						if ((Length[l]>left) && (Length[l]<=right))
							frq(f,region-1,y,qtr-1,a) += lf[l];
					}
				}
				for (int l=0; l<nb_intervals; l++)
					if ((Length[l]>right)&& (Length[l]<param.length[sp][nb_ages-1]+10))
						frq(f,region-1,y,qtr-1,nb_ages-1) += lf[l];
			}
		}
		littxt.close();
	} else {
		cout << endl << "WARNING : Cannot read file " << filename.c_str() << endl;
	}
	cout << "Number of size distributions used from LF data: " << nbsd << endl;
	cout << endl;
}

void CReadWrite::read_lf_fine(CParam& param, string filename, const float startdate, const float enddate, const int sp)
{ 
//////////////////////////////////////////////////////////
// data used by Multifan-CL, file header has been changed:
// yr qtr mo fishery region lf1..lfK
	const int nb_fishery = param.get_nbfishery();
	const int a0	     = param.sp_a0_adult[sp];
	const int nb_ages    = param.sp_nb_cohorts[sp];
	int nbsd = 0; //count the number of size distributions to be used

	ifstream littxt(filename.c_str());
	if (littxt){
		cout << endl << "Reading LF data from the file: " << endl << filename.c_str() << endl;
		int nb_regions, nb_fleets, nb_records, nb_intervals, l1, dl;

		littxt >> nb_regions >> nb_fleets >> nb_records;
		ivector region(0,nb_regions-1);
		region.initialize();
		//temporal
		int num;
		float c1,c2,c3,c4; // not used for the moment
		for (int r=0; r<nb_regions; r++)
			littxt >> num >> c1 >> c2 >> c3 >> c4;	    
    
		//read the line 
		//flexible to different file contents, namely to the presence of weight data
		//read number of size groups, first size and step
		littxt >> nb_intervals >> l1 >> dl;

		dvector lf(0,nb_intervals-1);
		lf.initialize();
		dvector Length(0,nb_intervals-1);
		dvector bin_bound(0,nb_intervals);
		Length.initialize();
		bin_bound.initialize();
		for (int l=0; l< nb_intervals; l++)
			Length(l) = l1+dl*(l+0.5); 
		for (int l=0; l<= nb_intervals; l++)
			bin_bound(l) = l1+dl*l; 

		// skip the header
		string s;
		getline(littxt,s,'\n');
		getline(littxt,s,'\n');

		int year, qtr, month, reg;
		string fishery;
		for (int n=0; n<nb_records; n++){
			littxt >> year >> qtr >> month >> fishery >> reg ;
			if (month==0) month = qtr*3-1;

			float frequency;
			lf.initialize();
			//read the whole line (some files contain weight frequencies data)
			getline(littxt,s,'\n');
			std::istringstream iss(s);
			for (int l=0; l<nb_intervals; l++){
				iss >> frequency;
				if (frequency == -1) break;
				lf[l] = frequency;	
			}

			bool fishery_to_throw = true;
			int f;
			for (f=0;f<nb_fishery;f++){
				if (param.list_fishery_name[f] == fishery){
					fishery_to_throw = false; 
					break;
				}		
			}
			if (fishery_to_throw) continue;
			float date = param.fdate(year, month);			
			if ((date >= startdate) && (date <=enddate) && qtr != 0){
				nbsd++;

				//double L_pr = param.length(sp,a0-1);
				int y=year-(int)startdate;
				frq(f,reg-1,y,qtr-1).initialize();

				//interpolate observed LF on modeled length intervals
				//taking into account the bin sizes
				for (int l=0;l<nb_intervals; l++){
				   if (bin_bound[l]>=param.length_bins(sp,a0) && bin_bound[l+1]<=param.length_bins(sp,nb_ages)){
					if (lf[l]>0){
						for (int a=a0; a<nb_ages; a++){
							bool left_in = false;
							bool right_in = false;
							bool both_in = false;
							bool part_in = false;
							if (bin_bound[l]>=param.length_bins(sp,a) && 
							    bin_bound[l]<=param.length_bins(sp,a+1))
								left_in = true;
							if (bin_bound[l+1]>=param.length_bins(sp,a) && 
							     bin_bound[l+1]<=param.length_bins(sp,a+1))
								right_in = true;
							if (left_in && right_in) {
								both_in = true;
								left_in = !both_in;
								right_in= !both_in;
							}
							if (bin_bound[l]<param.length_bins(sp,a) && 
							     bin_bound[l+1]>param.length_bins(sp,a+1))
								part_in = true;
				
							if (!left_in && !right_in && !both_in && !part_in) continue;
	
							double deltal = 0;
							if (left_in)
								deltal = param.length_bins(sp,a+1)-bin_bound[l];
							if (right_in)
								deltal = bin_bound[l+1]-param.length_bins(sp,a);
							if (both_in) 
								deltal = dl;
							if (part_in) 
								deltal = param.length_bins(sp,a+1)-param.length_bins(sp,a);
							frq(f,reg-1,y,qtr-1,a) += lf[l] * deltal/dl;
						}
					}
				    }
				    else if (bin_bound[l]<param.length_bins(sp,a0)) 
					frq(f,reg-1,y,qtr-1,a0) += lf[l];
				    else if (bin_bound[l+1]>param.length_bins(sp,nb_ages)) 
					frq(f,reg-1,y,qtr-1,nb_ages-1) += lf[l];

				}

				//control check						
				if (sqrt(pow(sum(lf) - sum(frq(f,reg-1,y,qtr-1)),2))>0.001) {
					cout << f << " "<< reg<< " "<<y << " " << qtr << " " << 
					sum(lf)<< " " << sum(frq(f,reg-1,y,qtr-1)) << endl;
					cout << lf << " " <<  endl << frq(f,reg-1,y,qtr-1) << endl << endl;
					exit(1);
				}

				//!!!TEST here: constrain the sample size by 1000 if S>1000 (see robust LF likelihood)
				//must be taken out to pre-processing
				double smax = 1000.0; //do not allow samples greater 1000
				double sumQ = sum(frq(f,reg-1,y,qtr-1));
				if (sumQ>smax){
					//put warning message HERE
					for (int a=a0; a<nb_ages; a++)
						frq(f,reg-1,y,qtr-1,a) *= smax/sumQ;
				}
			}
		}
		littxt.close();
	} else {
		cout << endl << "WARNING : Cannot read file " << filename.c_str() << endl;
	}
	cout << "Number of fine-resolution size distributions used from LF data: " << nbsd << endl;
	cout << endl;
}


void CReadWrite::read_frq_data(CParam& param, PMap& map, const float startdate, const float enddate, const int sp)
{
	const int nby = (int)(enddate) -(int)(startdate)+1;
	const int nbq = 4; 
	int nb_fishery = param.get_nbfishery();
	const int a0	     = param.sp_a0_adult[sp];
	const int nb_ages    = param.sp_nb_cohorts[sp];
	int nb_region  = param.nb_region;
	frq.allocate(0,nb_fishery-1,0,nb_region-1,0,nby-1,0,nbq-1,a0,nb_ages-1);
	frq.initialize();

	if (param.frq_like[sp]){
		if (param.flag_twin) {
			read_pred_frq_data(param,param.file_frq_data[0],startdate,enddate,sp);
			return;
		}
	
		if (!param.use_lf_regstruc){
			read_lf_WCPO(param,param.file_frq_data[0],startdate,enddate,sp);
			if (param.nb_frq_files==2)
				read_lf_WCPO(param,param.file_frq_data[1],startdate,enddate,sp);
				//read_lf_EPO(param,param.file_frq_data[1],startdate,enddate,sp); //for skj only!
		} else {
			for (int f=0; f<param.nb_frq_files; f++)
				read_lf_fine(param,param.file_frq_data[f],startdate,enddate,sp);
			}
	}
	//after reading the data do the check: 
	//if nbbins for which LF>0 is < 0.1*nb_ages, delete the record
	for (int f=0; f<nb_fishery; f++){
		for (int r=0; r<nb_region; r++)
			for (int y=0; y<nby; y++)
				for (int q=0; q<nbq; q++){
					int nbinsfil = 0;
					for (int a=a0; a<nb_ages; a++)
						if (frq(f,r,y,q,a)>0) nbinsfil++;
					if (nbinsfil>0 && nbinsfil<=2){
						frq(f,r,y,q).initialize();
					}
				}
	}
	//write the file with observations
	d4_array sum_frq(0,nb_fishery-1,0,nb_region-1,0,nbq,a0,nb_ages-1);
	d3_array sum_frq_fishery(0,nb_fishery-1,0,nbq,a0,nb_ages-1);
	d3_array sum_frq_region(0,nb_region-1,0,nbq,a0,nb_ages-1);
	sum_frq.initialize();
	sum_frq_fishery.initialize();
	sum_frq_region.initialize();
	//we shall aggregate only those LF data which will be taken into account in the likelihood, i.e. on the mask...
        for (int i=map.imin; i <= map.imax; i++){
        for (int j=map.jinf[i] ; j<=map.jsup[i] ; j++){
        if (map.carte[i][j]){
	for (int f=0; f<nb_fishery; f++){
		for (int reg=0; reg<nb_region; reg++){
		        int r = param.area_sp_B[sp][reg]-1;
			if ((i>=map.regimin[reg]) && (i<map.regimax[reg]) &&(j>=map.regjmin[reg]) && (j<map.regjmax[reg])){
			for (int q=0; q<nbq; q++){
				for (int y=0; y<nby; y++){
					if (norm(frq(f,r,y,q))){
						for (int a=a0; a< nb_ages; a++){
							sum_frq(f,r,q,a) += frq(f,r,y,q,a);
							sum_frq_fishery(f,q,a) += frq(f,r,y,q,a);
							sum_frq_region(r,q,a) += frq(f,r,y,q,a);
						}
					}
				}
			}
		}
}
	}}}}
        for (int f=0; f<nb_fishery; f++){
                for (int reg=0; reg<nb_region; reg++){
                        int r = param.area_sp_B[sp][reg]-1;
			for (int q=0; q<nbq; q++){
				for (int a=a0; a< nb_ages; a++){
					sum_frq(f,r,nbq,a) += sum_frq(f,r,q,a);  
					sum_frq_fishery(f,nbq,a) += sum_frq_fishery(f,q,a);
					sum_frq_region(r,nbq,a) += sum_frq_region(r,q,a);
				}
			}
		}
	}

	//write the file with read observations aggregated over Seapodym's age classes during simulation period
	string filename	= param.strdir_output + param.sp_name[0] + "_LF_obs.txt";
	ofstream ecritLF(filename.c_str());

	if (ecritLF) {
		ecritLF << "length" <<'\t' ;
		for (int reg=0; reg< param.nb_region_sp_B[sp]; reg++)
			for (int f=0; f< nb_fishery; f++)
				if (param.mask_fishery_sp[sp][f])
					ecritLF<<param.list_fishery_name[f]<<"_"<<param.sp_name[sp]<<"_region_"<<param.area_sp_B[sp][reg] <<'\t';
			
		for (int f=0; f< nb_fishery; f++)
			if (param.mask_fishery_sp[sp][f])
				ecritLF<<"sum_" << param.list_fishery_name[f]<<"_"<<param.sp_name[sp]<<'\t';
					

		for (int reg=0; reg< param.nb_region_sp_B[sp]; reg++)	
			ecritLF<<"sum_" << param.sp_name[sp]<<"_region_"<<param.area_sp_B[sp][reg] <<'\t';		

		ecritLF << endl;
		
		for (int q=0; q<=nbq; q++){
			ecritLF <<"Quarter "<< q+1 << endl;
			for (int a=a0; a< nb_ages; a++){
				ecritLF << param.length(sp,a) << '\t';
				for (int r=0; r< param.nb_region_sp_B[sp]; r++)
					for (int f=0; f< nb_fishery; f++)
						if (param.mask_fishery_sp[sp][f]){
							int reg = param.area_sp_B[sp][r]-1;
							ecritLF << sum_frq(f,reg,q,a) << '\t'; 
						}
				for (int f=0; f< nb_fishery; f++)
					if (param.mask_fishery_sp[sp][f])				
						ecritLF << sum_frq_fishery(f,q,a) << '\t';

				for (int r=0; r< param.nb_region_sp_B[sp]; r++){
					int reg = param.area_sp_B[sp][r]-1;
					ecritLF << sum_frq_region(reg,q,a) << '\t';
				}

				ecritLF << endl;
					
			}
		}
			
		ecritLF  << endl;
		ecritLF.close();
	}
	else
	{
		cout<<endl<<"WARNING : Cannot write file "<< filename.c_str()<<endl;
	}
}

void CReadWrite::write_frq_data(CParam& param, int sp, int year, int qtr, d3_array frq, bool FILEMODE)
{
	///int nb_ages = param.sp_nb_age_class_ad[sp];
	const int a0	     = param.sp_a0_adult[sp];
	const int nb_ages    = param.sp_nb_cohorts[sp];
	int nb_fishery = param.nb_fishery_by_sp(sp);
	int nb_region  = param.nb_region;


	string filename	= param.str_dir_fisheries + "skj_LF_pred.txt";
	ofstream ecritLF;
	if (!FILEMODE){
		ecritLF.open(filename.c_str(), ios::out);
		if (ecritLF) {
			ecritLF << "year" <<'\t' << "quarter" <<'\t'<< "fisheries"<<'\t'<<"region" <<'\t';
			for (int a=0; a<nb_ages; a++)
				ecritLF << a <<'\t';
			ecritLF << endl;
		}
		else cout<<endl<<"WARNING : Cannot write file "<< filename.c_str()<<endl;
	}
	else {
		ecritLF.open(filename.c_str(), ios::app);
		if (ecritLF) {
			for (int f=0; f<nb_fishery; f++){
				for (int r=0; r<nb_region; r++){
					int reg = param.area_sp_B[sp][r]-1;
					ecritLF << year << '\t' << qtr <<'\t'<< f <<'\t'<< reg <<'\t';
					for (int a=a0; a<nb_ages; a++)
						ecritLF << frq(a,f,reg) <<'\t';
					ecritLF << endl;
				}
			}
		}
		else cout<<endl<<"WARNING : Cannot write file "<< filename.c_str()<<endl;
	}
}

void CReadWrite::read_pred_frq_data(CParam& param, string filename, const float startdate, const float enddate, const int sp)
{
	///int nb_ages = param.sp_nb_age_class_ad[sp];
	const int a0	     = param.sp_a0_adult[sp];
	const int nb_ages    = param.sp_nb_cohorts[sp];
	int nb_fishery = param.nb_fishery_by_sp(sp);
	int nb_region  = param.nb_region;

	ifstream littxt(filename.c_str());
	if (littxt){
		string s;
		getline(littxt,s,'\n');
		
		int year, qtr, f, reg;
		int nb_records = (int)((enddate-startdate)*4)*nb_fishery*nb_region;
		for (int n=0; n<nb_records; n++){
			littxt >> year >> qtr >> f >> reg;
			int y = year-(int)startdate;
			for (int a=a0; a<nb_ages; a++)
				littxt >> frq(f,reg,y,qtr-1,a);
		}	
	}
	else cout<<endl<<"WARNING : Cannot read file "<< filename.c_str()<<endl;
}


void CReadWrite::get_LF_qtr_data(CParam& param, d4_array LF_qtr_obs, int y, const int q)
{
	LF_qtr_obs.initialize();
	int nb_species = param.get_nbspecies();
	int nb_fishery = param.get_nbfishery();
	y -= (int)param.save_first_yr;

	for (int sp=0; sp<nb_species; sp++){
		if (param.frq_like[sp]){
			const int a0	     = param.sp_a0_adult[sp];
			const int nb_ages    = param.sp_nb_cohorts[sp];
			int nb_region = param.nb_region_sp_B[sp];
			int k=0;
			for (int f=0; f<nb_fishery; f++)
				if (param.mask_fishery_sp[sp][f]){
					for (int age=a0; age<nb_ages; age++)
						for (int r=0; r<nb_region; r++){
							int reg = param.area_sp_B[sp][r]-1;
							LF_qtr_obs(sp,age,k,r) = frq(f,reg,y,q,age); 
						}
					k++;
				}
		}
	}	
}

