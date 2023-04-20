/// ReadWrite class containes routines to read/write input/output files
///////////////////////////////////////////////////////////////////////

#ifndef __ReadWrite_h__
#define __ReadWrite_h__

#include <string>
#include <vector>

#include "Matrices.h"
#include "Utilities.h"

using std::string;
using std::vector;

/*!
\brief Class that reads and stores redistributed fishing effort data
\details Fishing effort redistributed to the model resolution will be read and
used in Calpop class for each i,j to compute fishing mortality rates.
*/

class fishing_effort {
   public:
    fishing_effort(){/*DoNothing*/};
    virtual ~fishing_effort(){/*DoNothing*/};

   private:
    int i, j;
    double efrt;  /// fishing effort

   public:
    int get_i() { return i; }
    int get_j() { return j; }
    double get_effort(void) { return efrt; };

    void set_effort(int ii, int jj, double ee) {
        i = ii;
        j = jj;
        efrt = ee;
    }
};

/*!
\brief Class that reads and stores all fishing data
\details Fishing data to be stored in SEAPODYM: i,j indices, lon/lat coordinates
(center of the fishing area), effort and catch
*/
// TO DO: make class fishery_record inheritant of class fishing_effort
class fishery_record {
   public:
    fishery_record(){/*DoNothing*/};
    virtual ~fishery_record(){/*DoNothing*/};

   private:
    int i, j;
    double lon, lat;  /// lon-lat coordinates
    double efrt;      /// fishing effort
    double ctch;      /// 2015: catch for a single species, can be modified to
                  /// multispecies if needed

   public:
    int get_i() { return i; }
    int get_j() { return j; }
    double get_lon() { return lon; }
    double get_lat() { return lat; }
    double get_effort(void) { return efrt; };
    double get_efflon(void) { return lon; };
    double get_efflat(void) { return lat; };
    double get_catch(void) { return ctch; };

    // void create_catch(int nbsp){ctch.allocate(0,nbsp-1);};
    void set_record(
        double longitude, double latitude, int ii, int jj, double ee,
        double cc) {
        lon = longitude;
        lat = latitude;
        i = ii;
        j = jj;
        efrt = ee;
        ctch = cc;
    }
    // void add_record(double ee, dvector cc){efrt+=ee; ctch+=cc;}
    void change_coord(double longitude, double latitude, int ii, int jj) {
        lon = longitude;
        lat = latitude;
        i = ii;
        j = jj;
    }
};

/*!
\brief IO class
\details Class functions are accessible through all computational classes. All
types of input data are read here, any new output writing routines must be
placed here as well.
*/
class CReadWrite {
   public:
    CReadWrite() { fisheries_data_read = false; /*DoNothing*/ };
    virtual ~CReadWrite() { delete_fisheries_rec(); };

    friend class fishery_record;
    friend class fishing_effort;

   private:
    DMATRIX mat2d_NoBorder;

    // fishery file variables
    int year0_fishing;
    IMATRIX position_file, nrecords;
    IVECTOR nbsp_file, rec_count;
    DVECTOR startdate_fishery, enddate_fishery;

    dmatrix mask_catch;

    i3_array position, numrec;
    i3_array position_rm, numrec_rm;

    int all_rec;
    int nrec_oceanmask;
    d5_array frq;
    bool fisheries_data_read;
    dvector fishery_reso;  // spatial resolution of fishing data
    dvector flon, flat;    // the grid for degraded fishing data (to be used in
                         // likelihood and catch_pred)

    fishery_record* frec;
    fishery_record* frec_tmp;
    fishing_effort* efr;
    fishing_effort* efr_tmp;

   public:
    //	vector<string> dymFileFPred ;
    vector<string> dymFileSpPred;
    vector<string> dymFileSpC;
    vector<string> dymFileSpLF;
    vector<string> dymFileSumSpLF;
    vector<string> dymFileSpCorr;
    vector<string> FileSpFR;
    //	string fileS ;
    //	string fileF ;

   public:
    void init_writing(CParam& param);

    void rbin_headpar(string file_in, int& nlong, int& nlat, int& nlevel);

    void rtxt_headpar(string file_in, int& nlong, int& nlat, int& nlevel);

    void rwbin_minmax(
        string file_io, double minvalstep, double maxvalstep);  // Gael Nov04

    void rtxt_mat2d(string file_in, DMATRIX& mat2d, int& nlong, int& nlat);

    void rbin_header(
        string file_in, string& idformat, int& idfunc, double& minval,
        double& maxval, int nlong, int nlat, int nlevel, double& startdate,
        double& enddate, DMATRIX& xlon, DMATRIX& ylat, DVECTOR& zlevel,
        IMATRIX& msksp);

    void wbin_header(
        string file_out, string& idformat, int& idfunc, double& minval,
        double& maxval, int nlong, int nlat, int nlevel, double& startdate,
        double& enddate, const DMATRIX& xlon, const DMATRIX& ylat,
        const DVECTOR& zlevel, const IMATRIX& msksp);

    void rtxt_header(
        string file_in, int nlong, int nlat, int nlevel, double& startdate,
        double& enddate, DVECTOR& xlon, DVECTOR& ylat, DVECTOR& zlevel,
        IMATRIX& msksp);

    // string wtxt_headpar(string varname, const int nb_reg, ivector area_sp_B,
    // string sp_name);//, string& header_part);

    void rbin_mat2d(
        string file_out, PMap& map, DMATRIX& mat2d, int nlat, int nlong,
        int nbytetoskip);

    void rbin_input2d(
        string file_in, PMap& map, DMATRIX& mat2d, int nbi, int nbj,
        int nbytetoskip);

    //	void rbin_input2d(string file_in, const imatrix& carte, DMATRIX& mat2d,
    //int nbi, int nbj, int nbytetoskip) ;

    void wbin_mat2d(
        string file_out, const DMATRIX& mat2d, int nlat, int nlong,
        bool FILEMODE);

    void wbin_transpomat2d(
        string file_out, const DMATRIX& mat2d, int nlong, int nlat,
        bool FILEMODE);

    void wtxt_header(
        string file_out, int nlong, int nlat, int nlevel, double& startdate,
        double& enddate, const DVECTOR& xlon, const DVECTOR& ylat,
        const DVECTOR& zlevel, const IMATRIX& msksp);

    void wtxt_mat2d(
        string file_out, const DMATRIX& mat2d, int nlat, int nlong,
        bool FILEMODE);

    void rtxt_col_lonlat(
        string file_in, DMATRIX& mat2d, int nlong, int nlat, DVECTOR& xlon,
        DVECTOR& ylat, int nbvar, int var);

    void wbin_fishery(string file_in, string file_out, int nbvar);

    void rbin_fishery(
        string file_in, DMATRIX& mat2d, CParam& param, int nbvar, int nvar,
        int yyyy, int mm);

    void InitSepodymFileDym(
        CParam& param, CMatrices& mat, int nb_mo, DVECTOR& zlevel,
        const IMATRIX& msksp);

    void SaveSepodymFileDym(CParam& param, PMap& map, CMatrices& mat);

    void SaveDymFile(
        PMap& map, CMatrices& mat, string file, const dmatrix& data,
        const int nlon, const int nlat);

    void InitFluxesCohortsFileTxt(CParam& param);

    void SaveFluxesCohortsFileTxt(
        CParam& param, CMatrices& mat, PMap& map, int day, int month, int year);

    void InitSepodymFileTxt(CParam& param);

    void SaveSepodymFileTxt(
        CParam& param, CMatrices& mat, PMap& map, dvector sumP, DVECTOR& sumF,
        DVECTOR& sumFprime, DVECTOR& sumF_area_pred,
        DVECTOR& sumF_required_by_sp, DVECTOR& mean_omega_sp, int day,
        int mois2, int yr2, int t_total, int qtr1, int qtr2, int nbi, int nbj);
    // Inna 05Jul05: reading GMB fisheries files format and destroying variables
    // after all
    void rbin_fishery_header(CParam& param);
    //	void rbin_fishery(CParam& param, D3_ARRAY& effort, D4_ARRAY& catch_obs,
    //float date, bool& fishing); 	void rbin_fishery(CParam& param, D3_ARRAY&
    //effort, float date, const int sp); 	void rbin_fishery(CParam& param, PMap&
    //map, dmatrix& effort, float date, const int f);
    //-----------------------------------------------------------------------------------
    // void read_fishery_data(CParam& param, const PMap& map);
    // void rbin_fishery_effort_rm(CParam& param, PMap& map);
    // void write_fishery_data(CParam& param, const PMap& map, CMatrices& mat,
    // const int sp, const int year, const int month, bool FILEMODE);
    void rtxt_fishery_data(
        CParam& param, const PMap& map, const int nbt, const int jday_spinup);
    void set_effort_rm(
        CParam& param, PMap& map, const int nbt, const int jday_spinup);
    void degrade_fishery_reso(
        CParam& param, PMap& map, const int nbt, const int jday_spinup);
    // void frec_rm(CParam& param, const PMap& map, fishery_record *frec, const
    // int f, const float lon, const float lat, const int y, const int m, const
    // float E, const dvector C);
    void set_frec_rm(
        CParam& param, const PMap& map, const int nbt, const int jday_spinup);
    void set_frec_rm_no_effort_fisheries(
        CParam& param, const PMap& map, const int nbt, const int jday_spinup);
    void delete_fisheries_rec(void);
    void get_catch(
        CParam& param, dmatrix& catch_obs, const int f, int y, const int m,
        const int sp);
    void get_effort(
        CParam& param, dmatrix& effort, const int f, int y, const int m);
    void get_effort_lonlat(
        CParam& param, dmatrix& effort, dmatrix& efflon, dmatrix& efflat,
        const int f, int y, const int m);
    void get_effort_rm(
        CParam& param, dmatrix& effort, const int f, int y, const int m);
    // void get_effort_rm(CParam& param, dmatrix& effort, const imatrix carte,
    // const int f, int y, const int m);
    void get_fishery_data(
        CParam& param, D3_ARRAY& effort, D4_ARRAY& catch_obs, int y,
        const int m);
    void get_fishery_data(
        CParam& param, D3_ARRAY& effort, D4_ARRAY& catch_obs, D3_ARRAY& efflon,
        D3_ARRAY& efflat, int y, const int m);
    void get_average_effort(
        CParam& param, D3_ARRAY& effort, D3_ARRAY& efflon, D3_ARRAY& efflat,
        const int nby, const int m);
    void get_average_effort_rm(
        CParam& param, dmatrix& effort, const int f, const int nby,
        const int m);
    void get_average_selectivity(
        PMap& map, CParam& param, dvector& swa, const ivector fisheries,
        const int nbf, const int nbt, const int nb_ages, const int sp,
        const int step_count);
    // void get_fishery_data_mpa(PMap& map, CParam& param, D3_ARRAY& effort,
    // D4_ARRAY& catch_obs, int y, const int m);
    void get_fishery_data_mpa(
        PMap&, CParam&, d3_array&, d4_array&, d3_array&, d3_array&, int, int);
    void mpa_areas_comp(PMap&, CParam&);
    void inc_obs_catch_mpa(
        PMap& map, CParam& param, dmatrix& catch_obs, const int sp);
    int get_numrec(const int f, const int y, const int m) {
        return numrec(f, y, m);
    }

    void read_lf_WCPO(
        CParam& param, string filename, const float startdate,
        const float enddate, const int sp);
    void read_lf_EPO(
        CParam& param, string filename, const float startdate,
        const float enddate, const int sp);
    void read_lf_fine(
        CParam& param, string filename, const float startdate,
        const float enddate, const int sp);
    void read_frq_data(
        CParam& param, PMap& map, const float startdate, const float enddate,
        const int sp);
    void get_LF_qtr_data(
        CParam& param, d4_array LF_qtr_obs, int y, const int q);

    void write_frq_data(
        CParam& param, int sp, int year, int qtr, d3_array frq, bool FILEMODE);
    void read_pred_frq_data(
        CParam& param, string filename, const float startdate,
        const float enddate, const int sp);

    //	dvector get_reso(){return fishery_reso;}
};

#endif
