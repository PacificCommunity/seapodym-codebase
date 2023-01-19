#ifndef __Date_h__
#define __Date_h__

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream> 
#include <cmath>
#include <cstdlib>
#include <Param.h> 


using  namespace std;

const unsigned char daysInMonth[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

/*!
\brief Class written by J.Jouanno to handle date format.
\details The model supports three date formats depending on the calendar: 360-day year (most frequently used), 365-day year and standard calendar with leap years. 
*/
class Date 
{
public:
        Date() {/*DoNothing*/};
        virtual ~Date() {/*DoNothing*/};

public:
        static void init_time_variables(CParam& param, int& Tr_step, int& nbt_spinup_tuna,
					int& jday_run, int& jday_spinup, int& nbstot, 
				       const int info, const int flagsimu) ;
        static void update_time_variables(const int t_count, const int deltaT, 
				const int date_mode, const int jday_spinup, int& jday, 
				int& day, int& month, int& year, int& newyear);
	static int get_nbstot(const int ndat000, const int ndatfin, const int jdays_run, const int deltaT, const int date_mode, ivector& rundates);
	static int get_nbt_before_first_recruitment(const int first_recruitment_date, const int ndatini, const int deltaT, const int date_mode);
        static int dym_startdate_run(CParam& param, const dvector zlevel_dym, const int nbstot);
        //static void zlevel_run(CParam& param, const dvector zlevel_dym, const int nbstot, dvector& zlevel, const int nbt_start_series);
        static int leapYear(int year) ;
	static int dayWithinMonth(int day, int month, int year);
	static unsigned long julday(int day, int month, int year);
    	static unsigned long clmjulday(int day, int month, int year);
	static unsigned long nlyjulday(int day, int month, int year);
	static unsigned long juldayy(int day, int month, int year);
    	static unsigned long clmjuldayy(int day, int month, int year);
	static unsigned long nlyjuldayy(int day, int month, int year);
 	static void dmy(unsigned long julnum, int& d, int& m, int& y);
	static void clmdmy(unsigned long julnum, int& d, int& m, int& y);
	static void nlydmy(unsigned long julnum, int& d, int& m, int& y);
	static void idatymd(const int ndat, int& year, int& month, int& day);
	static string MakeDate(int yr, int mo, int jr);
	static string MakeDate(int yr, int mo);
	static string MonthName(int mo);
	static int    Update_now_time(int yr, int month, int day);
	static string Update_now_time_str(int yr, int month, int day);
        static string Update_now_time_str_spinup(int month);
};
#endif
