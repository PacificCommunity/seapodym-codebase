#include "SeapodymCoupled.h"

//string get_date_recaptures(const char* rec_file_name);
int get_date_recaptures(string full_name);


tag_release ***rel;

/*!
\brief Reading tagging data function.
*/
const double pi = 3.14159265358979;


int get_date_recaptures(string full_name)
{
  size_t pos_from, pos_to;
  
  pos_from = full_name.find_last_of("_")+1;
  pos_to = full_name.find_last_of(".")-1;
  string rec_date_str = full_name.substr(pos_from,pos_to);
  int rec_date = atoi(rec_date_str.c_str());
  return rec_date;
}


void SeapodymCoupled::create_tag_recaptures(){

	float deltax, deltay;
	//DX, DY - model resolution in degrees
	float DX = (param->deltaX/60.0);
	float DY = (param->deltaY/60.0);
	//deltax, deltay - resolution of tagging data grid
	deltax = param->dx_tags; deltay = param->dy_tags;
	xr_tags = deltax/DX;
	yr_tags = deltay/DY;
	if (xr_tags < deltax/DX){
		deltax = xr_tags*DX;
		cout << "WARNING: tagging data resolution should be divisible of model resolution, set dx = " << deltax << " degree(s) "<< endl;
	}
	if (yr_tags < deltay/DY){
		deltay = yr_tags*DY;
		cout << "WARNING: tagging data resolution should be divisible of model resolution, set dy = " << deltay << " degree(s) " << endl;
	}

	//Correcting lonmin_tags and latmin_tags so the tag cells
	//would include the model grid cells entirely:
	float dx_deg = param->deltaX/60.0;
	float dy_deg = param->deltaY/60.0;
	int lon_i = (param->lonmin_tags- param->longitudeMin)/dx_deg + .5;
	int lat_j = (param->latitudeMax-param->latmax_tags)/dy_deg +.5;
	if (param->longitudeMin+lon_i*dx_deg != param->lonmin_tags ||
	    param->latitudeMax-lat_j*dy_deg != param->latmax_tags){
		param->lonmin_tags = param->longitudeMin+lon_i*dx_deg;
		param->latmax_tags = param->latitudeMax -lat_j*dy_deg;
		cout << "\n Reset tagging data grid corner to " << param->lonmin_tags << ", " << param->latmax_tags<< endl;
		cout << " so that the larger tagging data cells would include model cells entirely." << endl;
	}

	nx_obs = (param->lonmax_tags-param->lonmin_tags)/deltax;
	ny_obs = (param->latmax_tags-param->latmin_tags)/deltay;

	xlon.allocate(0,nx_obs);
	ylat.allocate(0,ny_obs);
	for (int ii=0; ii<=nx_obs; ii++)
		xlon(ii) = param->lonmin_tags + deltax*ii;
	for (int jj=0; jj<=ny_obs; jj++)
		ylat(jj) = param->latmax_tags - deltay*jj;
        if (!param->gcalc()){
                cout << endl << "GRID being used to aggregate tag recaptures:" << endl;
                cout << "nx=" << nx_obs << ", ny="<< ny_obs << endl;
                cout << "xlon: " << xlon << endl; cout << "ylat: " << ylat << endl;
        }

	rec_obs.allocate(0,nb_tagpops-1);
	tlib_obs.allocate(0,nb_tagpops-1);
	rec_pred.allocate(0,nb_tagpops-1);
	for (int n=0; n<nb_tagpops; n++){	
		rec_obs(n).allocate(0,nx_obs-1,0,ny_obs-1);
		tlib_obs(n).allocate(0,nx_obs-1,0,ny_obs-1);
		rec_pred(n).allocate(0,nx_obs-1,0,ny_obs-1);
		rec_obs(n).initialize();
		tlib_obs(n).initialize();
		rec_pred(n).initialize();
	}
	//aggregated observations and predictions being used 
	//in the likelihood calculations only
	//aggregation time step is given by the parameter
	rec_obs_like.allocate(0,nx_obs-1,0,ny_obs-1);
	rec_pred_like.allocate(0,nx_obs-1,0,ny_obs-1);
	rec_obs_like.initialize();
	rec_pred_like.initialize();

}

void SeapodymCoupled::ReadTaggingData(imatrix& nb_rel, ivector& t_count_rec)
{
	const int yyini  = param->ndatini/10000;
	const int mmini  = (param->ndatini - (yyini * 10000) ) / 100 ;
	int ddini  =  param->ndatini - ( yyini * 10000 ) - ( mmini *100 );
	//month starts at first day, but the time step is set in the middle 
	if (param->date_mode==3) ddini = ddini -14; 

	dvector Xlon, Ylat;
	Xlon.allocate(0,param->nlong-1);
	Ylat.allocate(0,param->nlat-1);
	Xlon = mat.xlon[0];
	Ylat = trans(mat.ylat)[0];


	int gauss_kernel = param->tag_gauss_kernel_on;
	dmatrix gkernel(map.imin, map.imax, map.jinf, map.jsup); 
	gkernel.initialize();

	//in kernel functions will need the length of the cohort
	//attn, violation of multi-species
	const int    agemax  = param->sp_nb_cohorts[0]-1;
	//const double lmax    = param->length[0][agemax]*0.01;

	t_count_rec.allocate(0,nb_tagpops-1);
	t_count_rec.initialize();

	tag_release ***rel_rtxt;

	int ntags_max = 2000; //in one file
	rel_rtxt = new tag_release** [nb_tagpops];
	for (int p=0; p<nb_tagpops; p++){
		rel_rtxt[p] = new tag_release* [nbt_total];
		for (int n=0; n<nbt_total; n++){
			rel_rtxt[p][n] = new tag_release[ntags_max];//Attn!
		}
	}
	ifstream rtxt;

	int drec, mrec, yrec;
	int jday_ini;
	//with default date_mode=3 use 360-days calendar:
        jday_ini = Date::clmjulday(ddini,mmini,yyini);
        if (param->date_mode<3) 
		jday_ini = Date::julday(ddini,mmini,yyini);

	int nbtot_tags_files = 0;
	for (int p=0; p<nb_tagpops; p++){

		string file_in = param->file_tag_data[p];
		//date of all recaptures in the cohort will be read from the name of the file
		int date_rec = get_date_recaptures(file_in);	
		Date::idatymd(date_rec,yrec,mrec,drec);
		int jdate_rec;
		//with default date_mode=3 use 360-days calendar:
	        jdate_rec = Date::clmjulday(drec,mrec,yrec);
	        if (param->date_mode<3) 
			jdate_rec = Date::julday(drec,mrec,yrec);
		int tt_rec = jdate_rec-jday_ini;//in days
		if (tt_rec>0){
			t_count_rec(p) = (int)(tt_rec/deltaT)+1; //in time steps
		} else {cerr << "Something wrong with reading the date from file with recaptures " << file_in.c_str() << ": " << date_rec << endl; exit(1);}
		
		rtxt.open(file_in.c_str(), ios::in);
	
		int nbtot_tags = 0;
		rtxt >> nbtot_tags;
		nbtot_tags_files += nbtot_tags;
		string s;
		//skip one line
		getline(rtxt,s,'\n');
		getline(rtxt,s,'\n');
		int yy,mm,dd,yy_rec,mm_rec,dd_rec;
		float lon,lat,lon_rec,lat_rec;//,len_rel,len_rec;
		string id,tag_no,len_rel,len_rec;//,dd_rec;
		//string id,tag_no;

		int i_mod = 0;
		int j_mod = 0;
		int age_mod = 0;
		if (rtxt){
			int n = 0;
			for (; n<nbtot_tags; n++){
				PMap* map_ptr=nullptr;
				if (param->use_tag_masks){
					if (n < static_cast<int>(tagmaps.size()))
						map_ptr = &tagmaps[n];
				}else{
					map_ptr = &map;
				}
				rtxt >> id >> tag_no >> yy >> mm >> dd >> lat >> lon >> yy_rec >> mm_rec >> dd_rec >> lat_rec >> lon_rec >> len_rel >> len_rec;				
				int jday_rel, jday_rec;
				//with default date_mode=3 use 360-days calendar:
				jday_rel = Date::clmjulday(dd,mm,yy);
				jday_rec = Date::clmjulday(dd_rec,mm_rec,yy_rec);
				if (param->date_mode<3) {
					jday_rel = Date::julday(dd,mm,yy);
					jday_rec = Date::julday(dd_rec,mm_rec,yy_rec);
				}

				//1. tag release time in (t_count-1)!
			 	// = the number of full deltaT since first date
				int tt = (int)(jday_rel-jday_ini)/deltaT;
				//time at liberty
				int days_liberty = jday_rec-jday_rel;
				float dtlib = (float)days_liberty/deltaT;

				//2. tag length to model's age index
				age_mod = 0;
				//for (int a=0; a<param->sp_nb_cohorts[0]; a++){
				for (int a=a0_adult(0); a<aN_adult(0); a++){
					float len= strtof(len_rel.c_str(),NULL);
					if (len >= param->length_bins[0][a] &&
					    len < param->length_bins[0][a+1]) 
						age_mod = a ;
				}

				int rec_added = 0;
				i_mod = param->lontoi(lon_rec);
				j_mod = param->lattoj(lat_rec);	
				if (lon_rec>xlon[0] && lon_rec<=xlon[nx_obs] 
				    && lat_rec<=ylat[0] && lat_rec>ylat[ny_obs]){
					if (map_ptr->carte(i_mod,j_mod))
						rec_added = 1;
				}

				//4. rel positions in ij of model grid
				int rel_added = 0;
				i_mod = param->lontoi(lon);
				j_mod = param->lattoj(lat);
				if (i_mod>=map.imin && i_mod<=map.imax && j_mod>=map.jmin && j_mod<=map.jmax){
					if (map_ptr->carte(i_mod,j_mod))
						rel_added = 1;
					//attn: generic way is to search over 8 surrounding cells
					//if (!map.carte(i_mod,j_mod) && map.carte(i_mod,j_mod-1)) { 
					//	j_mod = j_mod-1; rel_added = 1;
					//} else if (!map.carte(i_mod,j_mod) && map.carte(i_mod+1,j_mod)) { 
					//	i_mod = i_mod+1; rel_added = 1;
					//}
				}
				//TESTS of data validity
				if (tt<0) { 
					if (!param->gcalc())
						cout<< "WARNING: tag releases out of model time range: " << 
							p << " " << yy << " " << mm << " " 
							<< yy_rec << " "<< mm_rec<< endl;
					continue;
				}

				if (!rec_added){
					int i = param->lontoi(lon_rec);
					int j = param->lattoj(lat_rec);	
					if (!param->gcalc()){
				        	if (lon_rec>xlon[0] && lon_rec<=xlon[nx_obs] && 
					    	lat_rec<=ylat[0] && lat_rec>ylat[ny_obs]){	
							cout << "WARNING: tag recapture is on land: " << 
								p << " " << id << " " << yy_rec << " "<< mm_rec 
								<< " " << lat_rec << " " << lon_rec << "; mask: " 
								<< map_ptr->carte(i,j)<< endl;
						} else 
							cout<< "WARNING: tag recapture out of observational space: " << 
								p << " " << id << " " << yy_rec << " "<< mm_rec 
								<< " " << lat_rec << " " << lon_rec << endl;
					}
					continue;
				}

				if (!rel_added){
					if (!param->gcalc())
						cout<< "WARNING: tag release out of model domain: " << 
							p << " " << id << " " << yy << " " << mm << " "
							<< lat << " " << lon << "; mask: " 
							<< map_ptr->carte(i_mod,j_mod) << endl;
					continue;
				}
				//One more check might be necessary: discrepancy between length(age) at recapture
				//which depends on the implemented growth model
				if (age_mod < a0_adult(0) || age_mod >aN_adult(0)) {
					if (!param->gcalc())
						cout << "WARNING: age of tag release out of modelled lifespan " << 
							p << " " << id << " " << yy_rec << " "<< mm_rec 
							<< " "<< age_mod << endl;
					continue;
				}
				if (days_liberty<0){
					if (!param->gcalc())
						cout<< "WARNING: tag releases time is post time recapture: " << 
							p << " " << yy << " " << mm << " " 
							<< yy_rec << " "<< mm_rec<< endl;
					continue;
				}
				//end of TESTS
				
				//Now remember releases, which will be used in optimization
				
				int rel_index = nb_rel[p][tt];
				rel_rtxt[p][tt][rel_index].set_release(i_mod, j_mod, age_mod);
				nb_rel[p][tt]++;

				//observed recaptures
				if (!gauss_kernel){
					for (int ii=0; ii<nx_obs; ii++){
						for (int jj=0; jj<ny_obs; jj++){
							//3. rec position in ij of observational grid
							if (lon_rec>xlon[ii] && lon_rec<=xlon[ii+1] 
							&& lat_rec<=ylat[jj] && lat_rec>ylat[jj+1]){
								rec_obs(p,ii,jj) ++;
								tlib_obs(p,ii,jj) += dtlib;
							}
						}
					}
				} else {//Gaussian Kernel
					double sigx = sqrt(2*pow(lon_distance(lon,lon_rec,lat,lat_rec),2)/(2*(dtlib+1)));
					double sigy = sqrt(2*pow(lat_rec-lat,2)/(2*(dtlib+1)));
					//minimal error (100 nmi)
					//test BET 20200525: 2 degrees, i.e. 120 nmi
					if (sigx<120/(60*1.852)) sigx = 120/(60*1.852);
					if (sigy<120/(60*1.852)) sigy = 120/(60*1.852);
					float length_rel= strtof(len_rel.c_str(),NULL);
					float length_rec= strtof(len_rec.c_str(),NULL);
					//bound length at age at recapture
					if (length_rec==0 || length_rec<length_rel){
						int age_rec = age_mod+dtlib;
						if (age_rec>agemax) age_rec = agemax;
						length_rec = param->length(0,age_rec);
					}
					//finally, do not allow latitudinal error to be larger than longitudinal
					if (sigy>sigx) sigy = sigx;
					sigx /= cos(lat_rec*pi/180.0);
					gaussian_kernel(gkernel,Xlon,Ylat,lon_rec,lat_rec,sigx,sigy);
					
					for (int ii=0; ii<nx_obs; ii++){
						int ifrom = param->lontoi(xlon[ii]);
						int ito = ifrom + xr_tags;
						for (int jj=0; jj<ny_obs; jj++){
							int jfrom = param->lattoj(ylat[jj]);
							int jto   = jfrom + yr_tags;
							for (int i=ifrom; i<ito; i++){
								for (int j=jfrom; j<jto; j++){
									if (map_ptr->carte[i][j]){
										rec_obs(p,ii,jj) += gkernel(i,j);
										tlib_obs(p,ii,jj) += gkernel(i,j)*dtlib;
									}
								}
							}
						}
					}
				}
			}
			rtxt.close();
		} else 	cout<<endl<<"WARNING : Cannot read file "<<file_in<<endl;	
	}

	int nbtot_tags_used = 0;
	rel = new tag_release** [nb_tagpops];
	for (int p=0; p<nb_tagpops; p++){
		rel[p] = new tag_release* [nbt_total];
		for (int n=0; n<nbt_total; n++){
			int nbr = nb_rel[p][n];
			nbtot_tags_used += nbr;
			rel[p][n] = new tag_release[nbr];
		}
	}
        if (!param->gcalc()){
		cout << "Number of tag releases being read in files: "<< nbtot_tags_files << endl;
		cout << "Number of tag releases being used in the model: "<< nbtot_tags_used << endl;
		if (gauss_kernel!=0)
			cout << "TOTAL number of tag recaptures after Gaussian kernel transform: " << sum(rec_obs) << endl;
	}

	for (int p=0; p<nb_tagpops; p++){
		for (int n=0; n<nbt_total; n++){
			for (int nr=0; nr<nb_rel[p][n]; nr++){
				rel[p][n][nr] = rel_rtxt[p][n][nr];
			}
		}
	}
	//delete temporary class containers rel_rtxt
	for (int p=0; p<nb_tagpops; p++){
		for (int n=0; n<nbt_total; n++)
			delete [] rel_rtxt[p][n];
		delete [] rel_rtxt[p];
	}
	delete [] rel_rtxt;
}

double SeapodymCoupled::lon_distance(const double lon_rel, const double lon_rec, const double lat_rel, const double lat_rec){
//computed the longitudinal distance by grid method integration with 0.1deg latitudinal step
//returns longitudinal distance in 'degrees' corrected by latitude.
    
	if (lon_rel==lon_rec) return 0;
	if (lat_rel==lat_rec) return (lon_rec-lon_rel)*cos(lat_rel*pi/180.0);
	    
	double lat = lat_rel;
	int nlat = (int) abs(lat_rec-lat_rel)/0.5+1;
	double lat_step = (lat_rec-lat_rel)/nlat;
	double lon_step = (lon_rec-lon_rel)/((lat_rec-lat_rel)/lat_step);
	double lon_dist = 0;
	for (int i=0; i<nlat; i++){		
		lon_dist += lon_step*cos(lat*pi/180.0);
		lat += lat_step;
    	}
	return lon_dist;
}


void SeapodymCoupled::get_tag_releases(dvar4_array& density, i3_array& tagpop_age_solve, imatrix& nb_rel){
	for (int p=0; p<nb_tagpops; p++){
		//Read tagging data and initialize populations
		for (int nr = 0; nr < nb_rel[p][t_count-1]; nr++){
			int irel = rel[p][t_count-1][nr].get_i(); 
			int jrel = rel[p][t_count-1][nr].get_j(); 
			int age_cohort = rel[p][t_count-1][nr].get_age();
			double lat_corrected_area = cell_area / mat.lat_correction(jrel);
			tagpop_age_solve(p,t_count-1,age_cohort) = 1;
			density(p+1,age_cohort,irel,jrel) += 1.0/lat_corrected_area;
			//density(p+1,age_cohort,irel,jrel) += 1.0;
		}
	}
}

void SeapodymCoupled::delete_tag_releases(){
	for (int p=0; p<nb_tagpops; p++){
		for (int n=0; n<nbt_total; n++)
			delete [] rel[p][n];
		delete [] rel[p];
	}
	delete [] rel;
}

void SeapodymCoupled::gaussian_kernel(dmatrix& gauss_kernel, dvector x, dvector y, double lon, double lat, double rx, double ry){
	gauss_kernel.initialize();
	const int imin = map.imin; 
	const int imax = map.imax; 
	double outside = 0.0;
	double A  = 1/(2.0*pi*rx*ry);
	for (int i = imin; i <= imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin ; j <= jmax; j++){
			double G = A*exp((-(x(i-1)-lon)*(x(i-1)-lon)/(rx*rx)-(y(j-1)-lat)*(y(j-1)-lat)/(ry*ry))/2.0);
			if (map.carte[i][j]){
				gauss_kernel[i][j] = G;
			} else 
				outside += G;
		}
	}

	//scale sum to 1
	double gsum = sum(gauss_kernel);
//	cout << "for (lon,lat) = ("<< lon <<  "," << lat << ") and (sigx,sigy) = (" << rx<<","<< ry <<"), A = "<< A <<"; at sea: " << gsum << " on land: " << outside << endl;
	
	gauss_kernel = gauss_kernel/gsum;
}

