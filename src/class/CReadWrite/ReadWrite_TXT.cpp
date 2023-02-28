//#include "StdAfx.h"
#include "ReadWrite.h"
#include "SaveTimeArea.h"
#include "Numfunc.h"

//////////////////////////////////////////////////////////
//WRITE HEADERS OF FR_SP_AREA.TXT
//ecriture des noms de variable du fichier species_FR_area.txt
//////////////////////////////////////////////////////////

/*
void CReadWrite::InitFRspeciesFileTxt(CParam& param)
{
	const int nb_species = param.get_nbspecies();
	const int nb_forage = param.get_nbforage();

        //Generate the name of the files
        for (int sp=0; sp<nb_species;sp++)
                FileSpFR.push_back(param.strdir_output+param.sp_name[sp]+"_FR_area.txt");

        int nb_reg =0;
        ofstream ecritFR;
        for (int sp=0; sp<nb_species;sp++)
        {
                ecritFR.open(FileSpFR[sp].c_str(), ios::out);
                if (ecritFR)
                {
                        ecritFR << param.sp_name[sp]<<'\t'<<'\t'<<'\t'<<endl;
                        ecritFR << "young: from age "<< param.age_autonomous[sp]<<" to age "<<param.age_mature[sp]<<'\t'<<'\t'<<'\t';
                        ecritFR<<endl;
                        ecritFR << "Adult: from age "<<param.age_mature[sp]<<" to age " <<param.sp_nb_cohorts[sp] <<'\t'<<'\t'<<'\t';
                        ecritFR<<endl<<endl;

                        ecritFR << "year" <<'\t'<< "month"<< '\t';

			int nb_reg = param.nb_region_sp_B[sp];
			if (nb_reg==0) nb_reg = 1;

                        for (int a=0; a<nb_reg; a++)
                          for (int n=0; n< nb_forage; n++)
                                  ecritFR << param.sp_name[sp] << " young: area " << param.area_sp_B[sp][a] << " FR " << param.frg_name[n] << '\t';

                        for (int a=0; a<nb_reg; a++)
                          for (int n=0; n< nb_forage; n++)
                                  ecritFR << param.sp_name[sp] << " adult: area " << param.area_sp_B[sp][a] << " FR " << param.frg_name[n] << '\t';      

                        ecritFR << endl;
                        ecritFR.close();
                }
                else
                        cout<<endl<<"WARNING : Failed to initialize file "<<FileSpFR[sp]<<endl;
        }

}

/// Recording FR by species by area in files Species_FR.txt
void CReadWrite::SaveFRspeciesFileTxt(CParam& param, CMatrices& mat, PMap& map, int mois2, int yr2)
{

        CSaveTimeArea save;
	const int nb_species = param.get_nbspecies();
	const int nb_forage =  param.get_nbforage();

//	dmatrix FR_area_average;

        ofstream ecritFR;

        for (int sp=0; sp< nb_species;sp++){

		int nb_reg = param.nb_region_sp_B[sp];
		if (nb_reg==0) nb_reg = 1;


//		FR_area_average.allocate(0,nb_forage-1,0,nb_reg-1);
//		FR_area_average.initialize();


                ecritFR.open(FileSpFR[sp].c_str(), ios::app);

                if (ecritFR)
                {
                        ecritFR << yr2 <<'\t' << mois2 <<'\t';
//                        for (int n=0;n<nb_forage;n++){
//                                save.AverageByArea(map,FR_species,FR_area_average[n],cell_area,nb_reg);
                        for (int reg=0; reg< param.nb_region_sp_B[sp]; reg++)
				for (int n=0;n<nb_forage;n++)
                        		ecritFR  << mat.FR_young(sp,n,reg) <<'\t';

      //                  for (int n=0;n<nb_forage;n++)
        //                        save.AverageByArea(map,FR_species,FR_area_average[n],cell_area,nb_reg);

			for (int reg=0; reg<nb_reg; reg++)
				for (int n=0;n<nb_forage;n++)
                        		ecritFR  << mat.FR_adult(sp,n,reg) <<'\t';
//			

                        ecritFR << endl;
                        ecritFR.close();
                }
                else
                        cout<<endl<<"WARNING : Cannot write file "<<FileSpFR[sp]<<endl;
        }
}
*/

///Put here all variables that need to be initialized for writing
void CReadWrite::init_writing(CParam& param)
{
	const int nbi = param.get_nbi();
	const int nbj = param.get_nbj();
	
	//this mask will be used in comparison with MFCL estimations
        mask_catch.allocate(0,nbi-1,0,nbj-1);
	mask_catch.initialize(); 

	//Note, if no fishing, the use_mask_catch=false also
	if (!param.use_mask_catch)
		mask_catch = 1e3;
}

//////////////////////////////////////////////////////////
// ECRITURE EN-TETE FICHIERS DE SAUVEGARDE SEPODYM 'TEXT'
//////////////////////////////////////////////////////////
void CReadWrite::InitSepodymFileTxt(CParam& param)
{
//cerr << "Force Exit[" << __FILE__ << ':' << __LINE__ << "]: function hasn't been fixed yet\n";
	//////////////////////////////////////
	//WRITE HEADERS OF SUMDYM.TXT
	//ecriture des noms de variable du fichier SumDym.txt
	// fichier sommes mensuelles
	string fileSumDym	= param.strdir_output + "SumDym.txt";
	ofstream ecritSumDym;
	ecritSumDym.open(fileSumDym.c_str(), ios::out);

	const int nb_species = param.get_nbspecies();
	const int nb_fishery = param.get_nbfishery();
	const int nb_forage = param.get_nbforage();

	if (ecritSumDym)
	{
		//ecritSumDym <<"date"<< '\t'<<"tstep"<< '\t' << "P" <<'\t';
		ecritSumDym <<"date"<< '\t'<<"tstep"<< '\t' << "P in 10N-45N" << '\t' << "P in 10S-10N" <<'\t'<< "P in 35S-10S" <<'\t'<< "P total" <<'\t';
		for(int n = 0; n< nb_forage;n++)
		{
			//ecritSumDym <<"Fprime_"<<param.frg_name[n]<<'\t'
			ecritSumDym <<"F_"<<param.frg_name[n]<<'\t';
				//<<"F_"<<param.frg_name[n]<<" in pred. area"<<'\t'
				//<<"F_"<<param.frg_name[n]<<" required by pred."<<'\t'
				//<<"Mean Omega_"<<param.frg_name[n]<<'\t' ;
		}
		for (int sp=0; sp<nb_species;sp++)
		{
			ecritSumDym <<"B Larvae "<<param.sp_name[sp] <<'\t'
				<<"B Juv "<<param.sp_name[sp] <<'\t'
				<<"B Rec "<<param.sp_name[sp] <<'\t'
				<<"B Young "<<param.sp_name[sp] <<'\t'
				<<"B Adult "<<param.sp_name[sp] <<'\t'
				<<"B Total"<<param.sp_name[sp] <<'\t';
		}

		if (nb_fishery && !param.flag_no_fishing){
			for (int f=0;f<nb_fishery; f++){
				for (int sp=0;sp<nb_species;sp++){
					if (param.mask_fishery_sp[sp][f]){

						ecritSumDym <<"effort " << param.list_fishery_name[f] << '\t' ;
						ecritSumDym <<"obs C_"  << param.sp_name[sp]<<"_"<<param.list_fishery_name[f] <<'\t' ;
						ecritSumDym <<"pred C_" << param.sp_name[sp]<<"_"<<param.list_fishery_name[f] <<'\t' ;
						ecritSumDym <<"obs CPUE_" << param.sp_name[sp]<<"_"<<param.list_fishery_name[f] <<'\t' ;
						ecritSumDym <<"pred CPUE_" << param.sp_name[sp]<<"_"<<param.list_fishery_name[f] <<'\t' ;
					}
				}
			}
		}
		ecritSumDym << endl;
	}
	else
	{
		cout<<endl<<"WARNING : Cannot write file "<<fileSumDym<<endl;
	}

	//////////////////////////////////////
	//WRITE HEADERS OF SUMQAREA.TXT
	//ecriture des noms de variable du fichier SumQArea.txt
	// fichier sum by quarter and area (MFCL)
	string fileSumQA   = param.strdir_output + "SumQArea.txt";
	ofstream ecritSumQA;
	ecritSumQA.open(fileSumQA.c_str(), ios::out);
	
	//NOTE: cohorts larvae-recruits will be written in numbers, young-adult in tonnes
	if (ecritSumQA)
	{
		for (int sp=0; sp<nb_species;sp++)
			ecritSumQA << param.sp_name[sp]<< " regions: "<< '\t'<<'\t'<<'\t';
		ecritSumQA <<endl;
		for (int r=0; r<param.nb_region;r++)
			ecritSumQA << param.area[r]->lgmin <<" "<< param.area[r]->lgmax<< " " << param.area[r]->ltmin << " " << param.area[r]->ltmax << endl;
		
		for (int sp=0; sp<nb_species;sp++)
			ecritSumQA << "juvenile: from age 0 to " <<param.age_autonomous[sp]<<'\t'<<'\t'<<'\t';
		ecritSumQA <<endl;
		for (int sp=0; sp<nb_species;sp++)
			ecritSumQA << "young: from age "<< param.age_autonomous[sp]<<" to age" << param.age_mature[sp] << " (see maturity at age)"<<'\t'<<'\t'<<'\t';
		ecritSumQA <<endl;
		for (int sp=0; sp<nb_species;sp++)
			ecritSumQA << "Recruit: at age "<< param.age_recruit[sp] <<'\t'<<'\t'<<'\t';
		ecritSumQA <<endl;
		for (int sp=0; sp<nb_species;sp++)
			ecritSumQA << "Adult: from age "<<param.age_mature[sp]<<" to age " <<param.sp_a0_adult[sp] <<'\t'<<'\t'<<'\t';
		ecritSumQA <<endl;
		for (int sp=0; sp<nb_species;sp++)
			ecritSumQA << "Total biomass: from age "<<param.age_autonomous[sp]<< " to " <<param.sp_nb_cohorts[sp] <<'\t'<<'\t'<<'\t';
		ecritSumQA <<endl<<endl;
		ecritSumQA << "year" <<'\t'<<"month" <<'\t' << "day" << '\t';

		if (param.nb_region)
		{
			for (int sp=0; sp<nb_species;sp++)
			{
				int nb_reg = param.nb_region_sp_B[sp];
				if (nb_reg>0){

					for (int a=0; a<nb_reg; a++)
						ecritSumQA<<param.sp_name[sp] << " N Larvae region " << param.area_sp_B[sp][a] << '\t';
					ecritSumQA<< "Total N Larvae"<<'\t';
					for (int a=0; a<nb_reg; a++)
						ecritSumQA<<param.sp_name[sp] << " N Juv. region " << param.area_sp_B[sp][a] << '\t';
					ecritSumQA<< "Total N Juv."<<'\t';
					for (int a=0; a<nb_reg; a++)
						ecritSumQA<<param.sp_name[sp]<< " N Rec. region " << param.area_sp_B[sp][a] << '\t';
					ecritSumQA<< "Total N. Rec."<<'\t';
					for (int a=0; a<nb_reg; a++)
						ecritSumQA<<param.sp_name[sp]<< " B Young region " << param.area_sp_B[sp][a] << '\t';
					ecritSumQA<< "Total B. Young"<<'\t';
					for (int a=0; a<nb_reg; a++)
						ecritSumQA<<param.sp_name[sp]<< " B Adult region " << param.area_sp_B[sp][a] << '\t';
					ecritSumQA<< "Total B Adult"<<'\t';
					for (int a=0; a<nb_reg; a++)
						ecritSumQA<<param.sp_name[sp]<< " B tot. region " << param.area_sp_B[sp][a] << '\t';
					ecritSumQA<< "Total B"<<'\t';
				}
			}
		}
		ecritSumQA<<endl;
		ecritSumQA.close();
	}
	else
	{
		cout<<endl<<"WARNING : Cannot write file "<<fileSumQA<<endl;
	}

       //////////////////////////////////////
        // WRITE HEADERS OF SumEEZ.txt files (Inna 06/2011)
        // one file by EEZ:
        // date, N.Larvae, N.Juveniles, N.Recruites, B.Young, B.Adult, B.tot, E[f], Cobs[f], Cpred[f] 
        if (param.nb_EEZ!=0){
                for (int a=0; a<param.nb_EEZ; a++){
                        string fileSumEEZ   = param.strdir_output + "SumEEZ" + "_" + param.EEZ_name[a] + ".txt";
                        ofstream ecritSumEEZ;
                        ecritSumEEZ.open(fileSumEEZ.c_str(), ios::out);

                        if (ecritSumEEZ){

                                ecritSumEEZ << "year" <<'\t'<<"month" <<'\t' << "day" << '\t';
                                for (int sp=0; sp<nb_species;sp++){
                                //NOTE: cohorts larvae-recruits will be written in numbers, young-adult in tonnes

                                        ecritSumEEZ<<param.sp_name[sp]<< " N Larvae " << param.EEZ_name[a] << '\t';
                                        ecritSumEEZ<<param.sp_name[sp]<< " N Juv. " << param.EEZ_name[a] << '\t';
                                        ecritSumEEZ<<param.sp_name[sp]<< " N Rec. " << param.EEZ_name[a] << '\t';
                                        ecritSumEEZ<<param.sp_name[sp]<< " B Young " << param.EEZ_name[a] << '\t';
                                        ecritSumEEZ<<param.sp_name[sp]<< " B Adult " << param.EEZ_name[a] << '\t';
                                        ecritSumEEZ<<param.sp_name[sp]<< " B tot. " << param.EEZ_name[a] << '\t';
                                        if (nb_fishery && !param.flag_no_fishing){
                                                for (int f=0;f<nb_fishery; f++){
                                                        if (param.mask_fishery_sp[sp][f]){

                                                                ecritSumEEZ <<"Nb.obs " << param.list_fishery_name[f] << '\t' ;
                                                                ecritSumEEZ <<"effort " << param.list_fishery_name[f] << '\t' ;
                                                                ecritSumEEZ <<"obs C_"  << param.sp_name[sp]<<"_"<<param.list_fishery_name[f] <<'\t' ;
                                                                ecritSumEEZ <<"pred C_" << param.sp_name[sp]<<"_"<<param.list_fishery_name[f] <<'\t' ;
                                                                ecritSumEEZ <<"obs CPUE_"  << param.sp_name[sp]<<"_"<<param.list_fishery_name[f] <<'\t' ;
                                                                ecritSumEEZ <<"std obs CPUE_"  << param.sp_name[sp]<<"_"<<param.list_fishery_name[f] <<'\t' ;
                                                                ecritSumEEZ <<"pred CPUE_" << param.sp_name[sp]<<"_"<<param.list_fishery_name[f] <<'\t' ;
                                                        }
                                                }
                                        }
                                }
                                ecritSumEEZ<<endl;
                                ecritSumEEZ.close();
                        }
                        else
                                cout<<endl<<"WARNING : Cannot write file "<<fileSumEEZ<<endl;
                }
        }

        ////////////////////////////////////
        //WRITE HEADER OF MEANVAR.TXT
        for (int sp=0; sp<nb_species;sp++){
        string fileMeanVar  = param.strdir_output + param.sp_name[sp] + "_MeanVar.txt";
        ofstream ecritMeanVar;
        ecritMeanVar.open(fileMeanVar.c_str(), ios::out);

        if (ecritMeanVar){

                ecritMeanVar << "year" <<'\t'<<"month" <<'\t' << "day" << '\t';
                for (int sp=0; sp<nb_species;sp++){
                        int nb_ages = param.sp_nb_cohorts(sp);
                        for (int a=0; a< nb_ages; a++)
                                ecritMeanVar<<param.sp_name[sp]<< " mortality-at-age " << a << '\t';
                        for (int a=0; a< nb_ages; a++)
                                ecritMeanVar<<param.sp_name[sp]<< " speed-at-age " << a << '\t';
                        for (int a=0; a< nb_ages; a++)
                                ecritMeanVar<<param.sp_name[sp]<< " diffusion-at-age " << a << '\t';
                        for (int a=0; a< nb_ages; a++)
                                ecritMeanVar<<param.sp_name[sp]<< " temperature-at-age " << a << '\t';

                }
                ecritMeanVar<<endl;
        } else
                cout<<endl<<"WARNING : Cannot write file "<<fileMeanVar<<endl;
        }

	//////////////////////////////////////
	// WRITE HEADERS OF SPATIALCORR.TXT
	// ecriture des noms de variable du fichier SpatialCorr.txt
	// un fichier par espece
	// fichier correlations captures predites-observees
if (!param.flag_no_fishing){
	if (nb_fishery)	{
		//Generate the name of the files
		for (int sp=0; sp<nb_species;sp++)
		{
			dymFileSpCorr.push_back(param.strdir_output+param.sp_name[sp]+"_Spatial_Corr.txt");
		}
	}
	
	ofstream ecritcor;

	if (nb_fishery)	{ 
		for (int sp=0; sp<nb_species;sp++) {  

			ecritcor.open(dymFileSpCorr[sp].c_str(), ios::out);

			if (ecritcor) {
				ecritcor <<"date"<< '\t' ;
				for (int f=0;f<nb_fishery; f++){
					if (param.mask_fishery_sp[sp][f])	
						ecritcor <<"n" << '\t'<<"r " << param.list_fishery_name[f]<<"_"<<param.sp_name[sp]<< '\t'<<"prob" << '\t'
								<< "cpue_r " << param.list_fishery_name[f]<<"_"<<param.sp_name[sp]<< '\t'<<"prob" << '\t';
					
				}
				ecritcor<< endl;
				ecritcor.close();
			}
			else
			{
				cout<<endl<<"WARNING : Cannot write file "<<dymFileSpCorr[sp]<<endl;
			}
		}
	}
}
	//////////////////////////////////////
	//WRITE HEADERS OF SP_LF_Q_FISHERY.TXT
	//ecriture des noms de variable du fichier LF_QFishery.txt
    //fichier Frequences de longueurs par trimestre
if (!param.flag_no_fishing){	
	if (nb_fishery>0)
	{
		//Generate the name of the files
		for (int sp=0; sp<nb_species;sp++)
		{
			dymFileSpLF.push_back(param.strdir_output+param.sp_name[sp]+"_LF_Q_fishery.txt");
		}
	}

	ofstream ecritQtr_N;

	if (nb_fishery > 0){  
		for (int sp=0; sp<nb_species;sp++){
			const int a0	     = param.sp_a0_adult[sp];
			const int nb_ages    = param.sp_nb_cohorts[sp];
  
			ecritQtr_N.open(dymFileSpLF[sp].c_str(), ios::out);
			if (ecritQtr_N) {
				int nb_reg = param.nb_region_sp_B[sp];
				ecritQtr_N << nb_reg << " " << nb_fishery << endl;
				for (int a=a0;a<nb_ages;a++)
					ecritQtr_N  <<param.length[sp][a]<<'\t';
				ecritQtr_N << endl;
				ecritQtr_N << "Year "<< '\t' << "Quarter "<< '\t'<< "Month " << '\t'<< "Fishery " << '\t' << 
					      "Region " << '\t' << "LF[1] .. " << "LF[" << nb_ages << "]" << endl;
				

/*
				ecritQtr_N <<"Quarter "<< '\t'<< "Year "<< endl;
				ecritQtr_N << "length" <<'\t' ;
				for (int reg=0; reg< param.nb_region_sp_B[sp]; reg++){
					for (int f=0; f<nb_fishery; f++){
						if (param.mask_fishery_sp[sp][f])
							ecritQtr_N<<param.list_fishery_name[f]<<"_"<<param.sp_name[sp]<<"_region_"<<param.area_sp_B[sp][reg] <<'\t';
					}
				}
				ecritQtr_N  << endl;
*/
				ecritQtr_N.close();
			}
			else
			{
				cout<<endl<<"WARNING : Cannot write file "<<dymFileSpLF[sp]<<endl;
			}
		}
	}
}
}

void CReadWrite::InitFluxesCohortsFileTxt(CParam& param)
{
	int nb_species = param.get_nbspecies();
	int nb_reg = param.nb_EEZ;
	if (!nb_reg){ 
		nb_reg = param.nb_region_sp_B[0];//attn sp =0!
	}

	if (param.nb_region || param.nb_EEZ){
	    for (int sp=0; sp<nb_species;sp++){
		double mean_age_cohort = (double).5*param.sp_unit_cohort[sp][0];
		for (int a=0; a<param.sp_nb_cohorts[sp]; a++){
			std::ostringstream ostr;
			ostr << a;
			string fileJreg   = param.strdir_output + param.sp_name[sp] + "_FluxesRegion_age" + ostr.str() + ".txt";
			ofstream ecritJreg;
			ecritJreg.open(fileJreg.c_str(), ios::out);
			if (ecritJreg){

				ecritJreg << "Regional fluxes predicted for "<< param.sp_name[sp] << endl; 
				ecritJreg << "Mean age(months) " << mean_age_cohort << endl;
				ecritJreg << "Mean length(cm) " << param.length(sp,a) << endl;
				ecritJreg << "Mean weight(kg) " << param.weight(sp,a) << endl;
				ecritJreg << "Nb.regions " << nb_reg << endl;
				ecritJreg<<endl;
				ecritJreg.close();
			}
			else cout<<endl<<"WARNING : Cannot write file "<<fileJreg<<endl;
			if (a<param.sp_nb_cohorts[sp]-1)
	                	mean_age_cohort = mean_age_cohort + .5*param.sp_unit_cohort[sp][a]+.5*param.sp_unit_cohort[sp][a+1];
		}
	    }
	}
}

void CReadWrite::SaveFluxesCohortsFileTxt(CParam& param, CMatrices& mat, PMap& map, int day, int month, int year)
{
	int nb_species = param.get_nbspecies();
	int fluxes_between_polygons = 1;
	int nb_reg = param.nb_EEZ;
	if (!nb_reg){ 
		nb_reg = param.nb_region_sp_B[0];//attn sp =0!
		fluxes_between_polygons = 0;
	}
	
	if (param.nb_region || param.nb_EEZ){
	    for (int sp=0; sp<nb_species;sp++){
		for (int a=0; a<param.sp_nb_cohorts[sp]; a++){
			std::ostringstream ostr;
			ostr << a;
			string fileJreg   = param.strdir_output + param.sp_name[sp] + "_FluxesRegion_age" + ostr.str() + ".txt";
			ofstream ecritJreg;
			ecritJreg.open(fileJreg.c_str(), ios::app);
			if (ecritJreg){

				ecritJreg << "# "<< year << '-' << month << '-' << day << endl;
				if (!fluxes_between_polygons){
					int nbreg = param.nb_region_sp_B[sp];
					if (nbreg>0){
						for (int ri=0; ri<nb_reg; ri++){
							int reg_i = param.area_sp_B[sp][ri]-1;
							for (int rj=0; rj<nb_reg; rj++){
								int reg_j = param.area_sp_B[sp][rj]-1;
								ecritJreg << mat.fluxes_region(a,reg_i,reg_j) << '\t';
							}
							ecritJreg << endl;
						}
					}
				} else {
					for (int ri=0; ri<nb_reg; ri++){
						for (int rj=0; rj<nb_reg; rj++){
							ecritJreg << mat.fluxes_region(a,ri,rj) << '\t';
						}
						ecritJreg << endl;
					}

				}
				ecritJreg.close();
			}
			else cout<<endl<<"WARNING : Cannot append to the file "<<fileJreg<<endl;
		}
	    }
	}
}

///////////////////////////////////////////////////////////////////////////////
void CReadWrite::SaveSepodymFileTxt(CParam& param, CMatrices& mat, PMap& map,
				dvector sumP, DVECTOR& sumF, DVECTOR& sumFprime,
				DVECTOR& sumF_area_pred, DVECTOR& sumF_required_by_sp,
				DVECTOR& mean_omega_sp,int day, int mois2, int yr2, int t_total,
				/*int ntq,*/ int qtr1, int qtr2, int nbi, int nbj)
{
//cerr << "Force Exit[" << __FILE__ << ':' << __LINE__ << "]: function hasn't been fixed yet\n";
	CSaveTimeArea save;
	CNumfunc fonction;

	const int nb_species = param.get_nbspecies();
	const int nb_fishery = param.get_nbfishery();
	const int nb_forage = param.get_nbforage();

	
	//////////////////////////////////////
	//UPDATE SUMDYM.TXT
	double total=0.f;
	//opening the existing file in mode append
	string fileSumDym	= param.strdir_output +"SumDym.txt";
	ofstream ecritSumDym;
	ecritSumDym.open(fileSumDym.c_str(), ios::app);
	
	if (ecritSumDym)
	{
		// writing the date // Patrick 23Nov04 //Modified Inna 27Jan11
		ecritSumDym  << yr2 << "-" << mois2 << "-" << day <<'\t';

		//writing the data
		ecritSumDym << t_total << '\t'<< sumP[0] << '\t' << sumP[1] << '\t' << sumP[2] << '\t' << sumP[3] << '\t';
		for (int n=0;n<nb_forage;n++)
		{
			//ecritSumDym << sumFprime[n] <<'\t'
			ecritSumDym <<sumF[n] <<'\t';
				//<< sumF_area_pred[n]<<'\t'
				//<<sumF_required_by_sp[n] <<'\t'
				//<< mean_omega_sp[n] <<'\t';
		}
		for (int sp=0; sp<nb_species;sp++)
		{
			ecritSumDym << mat.sum_B_larvae[sp]  <<'\t'
				<< mat.sum_B_juv[sp]<<'\t'
				<< mat.sum_B_recruit[sp]<<'\t'
				<< mat.sum_B_young[sp]  <<'\t'
				<< mat.sum_B_adult[sp]  <<'\t'
				<< mat.sum_total_pop[sp]<<'\t';
		}
		if (nb_fishery>0 && !param.flag_no_fishing){ 
			int k = 0;
			for (int f=0; f<nb_fishery; f++){
				for (int sp=0;sp<nb_species;sp++){
					if (param.mask_fishery_sp[sp][f]){
	  
						// observed fishing effort
						// calcule et enregistre l'effort total de la pecherie
						total=0.f;
						for (int i=1; i< nbi-1; i++){
							for (int j=1; j < nbj-1 ; j++){
								if (map.carte(i,j))
									total += mat.effort[f][i][j];
							}
						}
						ecritSumDym << total << '\t';
//cout <<"total effort by fishery "<<f<< " " << total << endl;

				//for (int sp=0; sp< nb_species;sp++){ 
				//	if (param.mask_fishery_sp[sp][f]){
						// observed catch
						total=0.f;
						for (int i=1; i< nbi-1; i++){
							for (int j=1; j < nbj-1 ; j++){
								if (map.carte(i,j))
									total += mat.catch_obs(sp,k,i,j);
							}
						}
						ecritSumDym << total << '\t';
	
						// predicted catch
						total=0.f;
						for (int i=1; i< nbi-1; i++){
							for (int j=1; j < nbj-1 ; j++){
								if (map.carte(i,j))
									total += mat.catch_est(sp,k,i,j);
							}
						}
						ecritSumDym << total << '\t';	

//start test
						// observed CPUE
						total=0.f;
						int N = 0;
						double cpue = 0;
						for (int i=1; i< nbi-1; i++){
							for (int j=1; j < nbj-1 ; j++){
								if (map.carte(i,j)){
									double efr = mat.effort[f][i][j];
									double Cobs = mat.catch_obs(sp,k,i,j);
									//double Cest = mat.catch_est(sp,k,i,j);
									if (Cobs){
									  total += Cobs/efr;
									N++;
									}
								}
							}
						}
						if (N>0) cpue = total/N;
//cout << f << " OBSERVED: "<< N << " CPUE: " << total << " MEAN: " << total/(N+1) << " " << cpue;
						ecritSumDym << cpue << '\t';
						//ecritSumDym << total << '\t';
	
						// predicted CPUE
						total=0.f;
						cpue = 0;
						for (int i=1; i< nbi-1; i++){
							for (int j=1; j < nbj-1 ; j++){
								if (map.carte(i,j)){
									double efr = mat.effort[f][i][j];
									double Cobs = mat.catch_obs(sp,k,i,j);
									double Cest = mat.catch_est(sp,k,i,j);
									if (Cobs){
									//if (efr && mat.catch_obs(sp,k,i,j)){
									  total += Cest/efr;
									}
								}
							}
						}
						if (N>0) cpue = total/N;
//cout << "\t PREDICTED CPUE: " << total << " MEAN: " << total/(N+1) << " " << cpue << endl;
						ecritSumDym << cpue << '\t';	
						//ecritSumDym << total << '\t';	
//end test

			
						k++;
					}// else {
					//	ecritSumDym << 0 << '\t' << 0 << '\t';				
					//}
				}
			}
		}
		ecritSumDym << endl;
		ecritSumDym.close();
	}
	else
	{
		cout<<endl<<"WARNING : Cannot write file "<<fileSumDym<<endl;
	}

	//will be necessary in both SumQArea and SumByEEZ files
	dvector cell_area;
	dvector no_area;
	cell_area.allocate(map.jmin, map.jmax);
	no_area.allocate(map.jmin, map.jmax);
	cell_area.initialize();
	no_area.initialize();
	for (int j=map.jmin;j<=map.jmax;j++)
		cell_area[j] = param.cell_surface_area(j);
	no_area = 1.0;

	//////////////////////////////////////
	//UPDATE SUMQAREA.TXT
	//opening the existing file in mode append
	string fileSumQA   = param.strdir_output + "SumQArea.txt";
	ofstream ecritSumQA;
	ecritSumQA.open(fileSumQA.c_str(), ios::app);

		
	if (ecritSumQA)
	{
		// writing the date
		ecritSumQA << yr2 << '\t' << mois2 << '\t' << day << '\t';
/*
		dvector cell_area;//_corr;
		cell_area.allocate(map.jmin, map.jmax);
		cell_area.initialize();
		const double cell_area = 1.852*param.deltaX*1.852*param.deltaY;
		for (int j=map.jmin;j<=map.jmax;j++)
			cell_area[j] = cell_area / mat.lat_correction[j];
*/
		//int nbt_yr = (int) (365.25/param.deltaT);
		const int nbt  = param.get_nbt();//(int) (nbt_yr*(param.save_last_yr - param.save_first_yr + 1.0/nbt_yr) );
		
		//sum by region (rectangle)
		if (param.nb_region)
		{
			for (int sp=0; sp< nb_species;sp++)
			{
				int nb_reg = param.nb_region_sp_B[sp];
				if (nb_reg>0) {

					DVECTOR sum_area(0, nb_reg);

					save.SumByArea(map, mask_catch, mat.larvae[sp], sum_area, cell_area, nb_reg, nbt);
					for (int z =0; z< nb_reg+1;z++)
						ecritSumQA << sum_area[z] <<'\t';

					save.SumByArea(map, mask_catch, mat.juvenile[sp], sum_area, cell_area, nb_reg, nbt);
					for (int z =0; z< nb_reg+1;z++)
						ecritSumQA << sum_area[z] <<'\t';

					save.SumByArea(map, mask_catch, mat.recruit[sp], sum_area, cell_area, nb_reg, nbt);
					for (int z =0; z< nb_reg+1;z++)
						ecritSumQA << sum_area[z] <<'\t';

					save.SumByArea(map, mask_catch, mat.young[sp], sum_area, cell_area, nb_reg, nbt);
					for (int z =0; z< nb_reg+1;z++)
						ecritSumQA << sum_area[z] <<'\t';

					save.SumByArea(map, mask_catch, mat.adult[sp], sum_area, cell_area, nb_reg, nbt);
					for (int z =0; z< nb_reg+1;z++)
						ecritSumQA << sum_area[z] <<'\t';

					save.SumByArea(map, mask_catch, mat.total_pop[sp], sum_area, cell_area, nb_reg, nbt);
					for (int z =0; z< nb_reg+1;z++)
						ecritSumQA << sum_area[z] <<'\t';
				}
			}
		}
		ecritSumQA<<endl;
		ecritSumQA.close();
	}
	else
	{
		cout<<endl<<"WARNING : Cannot write file "<<fileSumQA<<endl;
	}

        //////////////////////////////////////
        // Update SumEEZ.txt files (Inna 06/2011)
        if (param.nb_EEZ!=0){

                int nlon = param.nlong;
                int nlat = param.nlat;
                for (int a=0; a<param.nb_EEZ; a++){

                        int EEZ_ID = param.EEZ_ID[a];
                        string fileSumEEZ   = param.strdir_output + "SumEEZ" + "_" + param.EEZ_name[a] + ".txt";
                        ofstream ecritSumEEZ;
                        ecritSumEEZ.open(fileSumEEZ.c_str(), ios::app);

                        if (ecritSumEEZ){
                                double sumEEZ = 0.0;
                                ecritSumEEZ << yr2 << '\t' << mois2 << '\t' << day << '\t';
                                for (int sp=0; sp<nb_species;sp++){
                                //NOTE: cohorts larvae-recruits will be written in numbers, young-adult in tonnes
                                        sumEEZ = save.SumByEEZ(map,EEZ_ID,mat.larvae[sp],cell_area,nlon,nlat);
                                        ecritSumEEZ << sumEEZ <<'\t';
                                        sumEEZ = save.SumByEEZ(map,EEZ_ID,mat.juvenile[sp],cell_area,nlon,nlat);
                                        ecritSumEEZ << sumEEZ <<'\t';
                                        sumEEZ = save.SumByEEZ(map,EEZ_ID,mat.recruit[sp],cell_area,nlon,nlat);
                                        ecritSumEEZ << sumEEZ <<'\t';
                                        sumEEZ = save.SumByEEZ(map,EEZ_ID,mat.young[sp],cell_area,nlon,nlat);
                                        ecritSumEEZ << sumEEZ <<'\t';
                                        sumEEZ = save.SumByEEZ(map,EEZ_ID,mat.adult[sp],cell_area,nlon,nlat);
                                        ecritSumEEZ << sumEEZ <<'\t';
                                        sumEEZ = save.SumByEEZ(map,EEZ_ID,mat.total_pop[sp],cell_area,nlon,nlat);
                                        ecritSumEEZ << sumEEZ <<'\t';

                                        if (nb_fishery && !param.flag_no_fishing){
                                                int k = 0; int Nb_obs = 0;
                                                for (int f=0;f<nb_fishery; f++){
                                                        if (param.mask_fishery_sp[sp][f]){
                                                                Nb_obs = save.NobsByEEZ(map,EEZ_ID,mat.effort[f],nlon,nlat);
                                                                ecritSumEEZ << Nb_obs <<'\t';
                                                                sumEEZ = save.SumByEEZ(map,EEZ_ID,mat.effort[f],no_area,nlon,nlat);
                                                                ecritSumEEZ << sumEEZ <<'\t';
                                                                sumEEZ = save.SumByEEZ(map,EEZ_ID,mat.catch_obs(sp,k),no_area,nlon,nlat);
                                                                ecritSumEEZ << sumEEZ <<'\t';
                                                                sumEEZ = save.SumByEEZ(map,EEZ_ID,mat.catch_est(sp,k),no_area,nlon,nlat);
                                                                ecritSumEEZ << sumEEZ <<'\t';
                                                                sumEEZ = save.SumByEEZ(map,EEZ_ID,mat.catch_obs(sp,k),mat.effort(f),nlon,nlat);
                                                                ecritSumEEZ << sumEEZ <<'\t';

								double mean_CPUE = 0;
								if (Nb_obs>0) mean_CPUE = sumEEZ/Nb_obs;
                                                                sumEEZ = save.StdCPUEByEEZ(map,EEZ_ID,mat.catch_obs(sp,k),mat.effort(f),mean_CPUE,Nb_obs,nlon,nlat);
                                                                ecritSumEEZ << sumEEZ <<'\t';

                                                                sumEEZ = save.SumByEEZ(map,EEZ_ID,mat.catch_est(sp,k),mat.effort(f),nlon,nlat);
                                                                ecritSumEEZ << sumEEZ <<'\t';

                                                                k++;
                                                        }
                                                }
                                        }
                                }
                                ecritSumEEZ<<endl;
                                ecritSumEEZ.close();
                        }
                        else
                                cout<<endl<<"WARNING : Cannot write file "<<fileSumEEZ<<endl;
                }
        }



        ////////////////////////////////////
        //UPDATE MEANVAR.TXT
        for (int sp=0; sp<nb_species;sp++){
        string fileMeanVar  = param.strdir_output + param.sp_name[sp] + "_MeanVar.txt";
        ofstream ecritMeanVar;
        ecritMeanVar.open(fileMeanVar.c_str(), ios::app);

        if (ecritMeanVar){

                ecritMeanVar <<  yr2 << '\t' << mois2 << '\t' << day << '\t';
                for (int sp=0; sp<nb_species;sp++){
                        int nb_ages = param.sp_nb_cohorts(sp);
                        for (int a=0; a< nb_ages; a++)
                                ecritMeanVar << mat.mean_mortality(sp,a) << '\t';
                        for (int a=0; a< nb_ages; a++)
                                ecritMeanVar << mat.mean_speed(sp,a) << '\t';
                        for (int a=0; a< nb_ages; a++)
                                ecritMeanVar << mat.mean_diffusion(sp,a) << '\t';
                        for (int a=0; a< nb_ages; a++)
                                ecritMeanVar << mat.mean_temperature(sp,a) << '\t';

                }
                ecritMeanVar<<endl;
        }
        else
                cout<<endl<<"WARNING : Cannot write file "<<fileMeanVar<<endl;
        }


	//////////////////////////////////////
	//UPDATE SPATIALCORR.TXT
	ofstream ecritcor;
	if (nb_fishery>0 && !param.flag_no_fishing)
	{	
		for (int sp=0;sp<nb_species;sp++)
		{	
			ecritcor.open(dymFileSpCorr[sp].c_str(), ios::app);
			if (ecritcor)
			{
				double cor_catch =0.f;
				double cor_cpue =0.f;
				double z_catch =0.f;
				double z_cpue =0.f;
				double prob_catch= 0.f;
				double prob_cpue= 0.f;
				int nn=0;
				ecritcor << mois2 <<"-" << yr2 <<'\t';
				int k=0;
				for (int f=0; f<nb_fishery; f++)
					if (param.mask_fishery_sp[sp][f]){

					fonction.corcatch(mat.catch_est[sp][k], mat.catch_obs[sp][k], map.imin, map.imax, map.jinf, map.jsup, nn, cor_catch, z_catch, prob_catch, map.carte, 0) ;
					fonction.corcpue(mat.catch_est[sp][k], mat.catch_obs[sp][k], mat.effort[f], map.imin, map.imax, map.jinf, map.jsup, nn, cor_cpue, z_cpue, prob_cpue, map.carte, 0) ;
					ecritcor << nn << '\t'<< cor_catch << '\t' << prob_catch << '\t' << cor_cpue << '\t' << prob_cpue << '\t';
					k++;
				}
				ecritcor << endl;
				ecritcor.close();
			}
			else
			{
				cout<<endl<<"WARNING : Cannot write file "<<dymFileSpCorr[sp]<<endl;
			}
		}
	}
	
	//////////////////////////////////////
	//UPDATE SP_LF_Q_FISHERY.TXT
	if (!param.flag_no_fishing){
	// quarterly catch in number by age(size)by species and fishery
	ofstream ecritQtr_N;

//exit(1);
	// si nouveau trimestre ENREGISTRE a condition qu'il y ait eu un passage (ntq) au minimum
	//if ( (nb_fishery>0) && (qtr2  != qtr1) && (ntq>0) )
//cout << param.nb_region_sp_B[0] << endl; exit(1);
	if (sum(param.mask_fishery_sp[0])!=0) {year0_fishing = (int)(param.ndatini/10000);}
	int year = yr2-year0_fishing;//(int)param.save_first_yr;
	if ( (nb_fishery>0) && (mois2==3 || mois2==6 || mois2==9 || mois2==12)) {
		for (int sp=0; sp< nb_species;sp++)
		{
			const int a0	     = param.sp_a0_adult[sp];
			const int nb_ages    = param.sp_nb_cohorts[sp];

			ecritQtr_N.open(dymFileSpLF[sp].c_str(), ios::app);
			if (ecritQtr_N)	{
				int k = 0;
				for (int f=0; f<nb_fishery; f++){
					if (param.mask_fishery_sp[sp][f]){
						for (int reg=0; reg< param.nb_region_sp_B[sp]; reg++){
							double sumpred = 0.0; double sumobs = 0.0;
							for (int a=a0;a<nb_ages;a++){
								sumpred += mat.C_N_sp_age_fishery_qtr[sp][a][reg][k];
								if (!param.nb_yr_forecast)
									sumobs  += frq(f,reg,year,mois2/3-1,a);
							}
							//if (sumpred>0 && sumobs>0){
							if (sumpred>0){
								ecritQtr_N << yr2 << '\t'<< qtr1 <<'\t'<< mois2 << '\t'<< param.list_fishery_name[f] << '\t' << reg+1 << '\t';
								for (int a=a0;a<nb_ages;a++)
									ecritQtr_N  << mat.C_N_sp_age_fishery_qtr[sp][a][reg][k] <<'\t';
								ecritQtr_N << endl;
							}
						} k++;
					}
				}
						
/*
				//ecritQtr_N <<"Quarter "<< qtr1 <<'\t'<< "Year "<< yr2 <<endl;
				for (int a=0;a<param.sp_nb_age_class_ad[sp];a++)
				{
					ecritQtr_N  <<param.length[sp][a]<<'\t';
					for (int reg=0; reg< param.nb_region_sp_B[sp]; reg++)
					{
						for (int k=0;k<param.nb_fishery_by_sp[sp]; k++)
						{
							ecritQtr_N  << mat.C_N_sp_age_fishery_qtr[sp][a][reg][k] <<'\t';
						}
					}
					ecritQtr_N << endl;
				}
*/
				ecritQtr_N.close();
			}
			else
			{
				cout<<endl<<"WARNING : Cannot write file "<<dymFileSpLF[sp]<<endl;
			}
		}
	}
	}
	//////////////////////////////////////
	// Write/Rewrite the file (s)
	//////////////////////////////////////
	//UPDATE SP_LF_Q_SUM.TXT
if (!param.flag_no_fishing){
	//sum of catch in number by age(size) by species, by region and fishery
	//sum for each of four quarter over all the series + sum of the four quarters

	//if ( (nb_fishery>0) && (qtr2  != qtr1) && (ntq>0) )
	if ( (nb_fishery>0) && (mois2==3 || mois2==6 || mois2==9 || mois2==12)) {
  
		for (int sp=0 ; sp<nb_species ; sp++) {

			const int a0	     = param.sp_a0_adult[sp];
			const int nb_ages    = param.sp_nb_cohorts[sp];

			for (int a=a0 ; a<nb_ages; a++) {
				
				const int nb_region_sp_B = param.nb_region_sp_B[sp];
				for (int r=0; r< nb_region_sp_B; r++) {
					int reg = param.area_sp_B[sp][r]-1;
					const int nb_fishery_by_sp = param.nb_fishery_by_sp[sp];
					for (int k=0;k<nb_fishery_by_sp; k++) { 

						// somme pour la climatologie trimestrielle par region et pecherie
						mat.Sum_C_N_sp_age_fishery_area[sp][a][qtr1-1][reg][k]
						+= mat.C_N_sp_age_fishery_qtr[sp][a][reg][k];
						// somme pour la climatologie trimestrielle par pecherie
						mat.Sum_C_N_sp_age_fishery_area[sp][a][qtr1-1][nb_region_sp_B][k]
						+= mat.C_N_sp_age_fishery_qtr[sp][a][reg][k];
						// somme pour la climatologie trimestrielle par region
						mat.Sum_C_N_sp_age_fishery_area[sp][a][qtr1-1][reg][nb_fishery_by_sp]
						+= mat.C_N_sp_age_fishery_qtr[sp][a][reg][k];

						// somme totale par pecherie et par region
						mat.Sum_C_N_sp_age_fishery_area[sp][a][5-1][reg][k]
						+= mat.C_N_sp_age_fishery_qtr[sp][a][reg][k];
						// somme totale par pecherie
						mat.Sum_C_N_sp_age_fishery_area[sp][a][5-1][nb_region_sp_B][k]
						+= mat.C_N_sp_age_fishery_qtr[sp][a][reg][k];
						// somme totale par region
						mat.Sum_C_N_sp_age_fishery_area[sp][a][5-1][reg][nb_fishery_by_sp]
						+= mat.C_N_sp_age_fishery_qtr[sp][a][reg][k];

						// somme grand total par trimestre
						mat.Sum_C_N_sp_age_fishery_area[sp][a][qtr1-1][nb_region_sp_B][nb_fishery_by_sp]
						+= mat.C_N_sp_age_fishery_qtr[sp][a][reg][k];

					}
				}
			}
		}
	}
}
	//////////////////////////////////////
	//WRITE HEADERS OF SP_LF_Q_SUM.TXT
	//ecriture des noms de variable du fichier LF_QFishery.txt
    //fichier Frequences de longueurs par trimestre
if (!param.flag_no_fishing){
	//if ( (nb_fishery>0) && (qtr2  != qtr1) && (ntq>0) )
	if ( (nb_fishery>0) && (mois2==3 || mois2==6 || mois2==9 || mois2==12)) {
		//Generate the name of the files
		for (int sp=0; sp<nb_species;sp++) {
	
			dymFileSumSpLF.push_back(param.strdir_output+param.sp_name[sp]+"_LF_Q_sum.txt");
	  	}
	}

	ofstream ecritQtr_S;

	//if ( (nb_fishery>0) && (qtr2  != qtr1) && (ntq>0) ){  
	if ( (nb_fishery>0) && (mois2==3 || mois2==6 || mois2==9 || mois2==12)) {  
		for (int sp=0; sp<nb_species;sp++){  

			ecritQtr_S.open(dymFileSumSpLF[sp].c_str(), ios::out);
			if (ecritQtr_S){

				ecritQtr_S << "length" <<'\t' ;
				for (int reg=0; reg< param.nb_region_sp_B[sp]; reg++){
					for (int f=0; f< nb_fishery; f++){
						if (param.mask_fishery_sp[sp][f])
							ecritQtr_S << param.list_fishery_name[f]<<"_" << param.sp_name[sp]<<"_region_"<<param.area_sp_B[sp][reg]<<'\t';
					}
				}
				for (int f=0; f< nb_fishery; f++){
					if (param.mask_fishery_sp[sp][f])
						ecritQtr_S  <<"sum_"<< param.list_fishery_name[f]<<"_" << param.sp_name[sp]<<'\t' ;
				}

				for (int reg=0; reg< param.nb_region_sp_B[sp]; reg++){

					ecritQtr_S  <<"sum_"<< param.sp_name[sp]<<"_region_" <<param.area_sp_B[sp][reg]<<'\t'  ;
				}
				ecritQtr_S  << endl;
				ecritQtr_S.close();
			}
			else
			{
				cout<<endl<<"WARNING : Cannot write file "<<dymFileSumSpLF[sp]<<endl;
			}
		}

	}
	//////////////////////////////
	// ECRITURE des sommes des captures par age, region pecherie pour les 4 trimestres et le total des 4
	double totalregion; 
	double totalpecherie;
	//if ( (nb_fishery>0) && (qtr2  != qtr1) && (ntq>0) )
	if ( (nb_fishery>0) && (mois2==3 || mois2==6 || mois2==9 || mois2==12)) {	
		for (int sp=0; sp< nb_species;sp++) {	

			const int a0	     = param.sp_a0_adult[sp];
			const int nb_ages    = param.sp_nb_cohorts[sp];

			ecritQtr_S.open(dymFileSumSpLF[sp].c_str(), ios::app);

			if (ecritQtr_S) {
				for (int q=0; q<5 ;q++) {	

					ecritQtr_S <<"Quarter "<< q+1 <<'\t'<<endl;

					for (int a=a0;a<nb_ages;a++) {	

						ecritQtr_S  <<param.length[sp][a]<<'\t';
						totalpecherie=totalregion=0.0;

						const int nb_fishery_by_sp = param.nb_fishery_by_sp[sp];
						const int nb_region_sp_B = param.nb_region_sp_B[sp];
						for (int r=0; r< nb_region_sp_B; r++) {
							int reg = param.area_sp_B[sp][r]-1;
							for (int k=0; k< nb_fishery_by_sp; k++) { 
								// somme pour la climatologie trimestrielle par region et pecherie
								ecritQtr_S  << mat.Sum_C_N_sp_age_fishery_area[sp][a][q][reg][k] <<'\t'  ;
							}
						}
						for (int k=0; k< nb_fishery_by_sp; k++) {
							// somme pour la climatologie trimestrielle par pecherie
							ecritQtr_S  << mat.Sum_C_N_sp_age_fishery_area[sp][a][q][nb_region_sp_B][k] <<'\t' ;
							totalpecherie+= mat.Sum_C_N_sp_age_fishery_area[sp][a][q][nb_region_sp_B][k];
						}
						for (int r=0; r< nb_region_sp_B; r++) {
							int reg = param.area_sp_B[sp][r]-1;
							// somme pour la climatologie trimestrielle par region
							ecritQtr_S  <<  mat.Sum_C_N_sp_age_fishery_area[sp][a][q][reg][nb_fishery_by_sp] <<'\t'  ;
							totalregion+= mat.Sum_C_N_sp_age_fishery_area[sp][a][q][reg][nb_fishery_by_sp];
						}
						ecritQtr_S <<totalpecherie <<'\t'  ;
						ecritQtr_S <<totalregion <<'\t'  ;
						ecritQtr_S << endl;
					}
				}	//trimestre suivant
				ecritQtr_S.close();
			} else {
				cout<<endl<<"WARNING : Cannot write file "<<dymFileSumSpLF[sp]<<endl;
			}
		}//espece et fichier suivant
	}
}
}
///////////////////////////

