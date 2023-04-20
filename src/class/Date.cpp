#include "Date.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#define new DEBUG_NEW
#endif

void Date::init_time_variables(
    CParam& param, int& Tr_step, int& nbt_spinup_tuna, int& jday_run,
    int& jday_spinup, int& nbstot, const int info, const int flagsimu) {
    //-------------------------------------------------------------------
    //  Initialisation of some variables for the time stepping. Should work
    //  for several mode of time stepping :
    //  date_mode=1 : input and output data each "deltaT" days, leap years
    //                are taken into account
    //  date_mode=2 : input and output data each "deltaT" days, leap years
    //                are NOT taken into account
    //  date_mode=3 : monthly input and output data
    //  data_mode=4: year of 360 days, month of 30 days
    //
    //  Initialized variables are :
    //
    //  Tr_sptep          :
    //  jday_run          : julian date at which runs start
    //  jday_spinup       : julian date at which spin-up start
    //
    //   MODIFICATIONS:
    //   --------------
    //      Julien (May 2008)

    //  (Inna, March 2011):
    //  NEW SEAPODYM DATE STANDARD
    //  (employed by Beatriz for 2002-2011 dataset with 7days step)
    //  The date in given for the first day of the time step
    //-------------------------------------------------------------------
    //        int nbstot=param.nbstot;
    int ndat000 = param.ndatini;
    int ndatfin = param.ndatfin;

    int deltaT = param.deltaT;

    int year, month, day;
    idatymd(ndat000, year, month, day);

    // Spinup to be REMOVED
    nbt_spinup_tuna = 0;

    int nbt_spinup_total = nbt_spinup_tuna;
    //	int nbt_no_forecast = t_count + nbt_spinup_tuna + datafit_period - 1;

    int dayspinup, monthspinup, yearspinup;
    switch (param.date_mode) {
        case 1: {
            jday_run = julday(day, month, year);
            jday_spinup = jday_run - (deltaT * (nbt_spinup_total));
            cout << "-------------Standard Time Stepping-----------------"
                 << endl;
            dmy(jday_spinup, dayspinup, monthspinup, yearspinup);
            break;
        }

        case 2: {
            jday_run = nlyjulday(day, month, year);
            jday_spinup = jday_run - (deltaT * (nbt_spinup_total));
            cout << "------------Noleap year Time Stepping---------------"
                 << endl;
            nlydmy(jday_spinup, dayspinup, monthspinup, yearspinup);
            break;
        }

        case 3: {
            jday_run = clmjulday(day, month, year);
            jday_spinup = jday_run - (deltaT * (nbt_spinup_total));
            cout << "-----Standard Time Stepping - 360 days calendar-----"
                 << endl;
            clmdmy(jday_spinup, dayspinup, monthspinup, yearspinup);
            break;
        }
    }
    int ndatspinup = yearspinup * 10000 + monthspinup * 100 + dayspinup;

    nbstot = get_nbstot(
        ndat000, ndatfin, jday_run, deltaT, param.date_mode, param.rundates);

    cout << "Number of time steps in simulation : " << nbstot << endl;
    cout << "Number of time steps in spin up    : " << nbt_spinup_total << endl;
    cout << "Delta T                            : " << deltaT << endl;
    if (nbt_spinup_total)
        cout << "Jday_spinup                        : " << jday_spinup << endl;
    cout << "Jday_run                           : " << jday_run << endl;
    if (nbt_spinup_total)
        cout << "Start Date SpinUp                  : " << ndatspinup << endl;
    cout << "Start Date Run                     : " << ndat000 << endl;
    cout << "End Date Run                       : " << ndatfin << endl;
    cout << "----------------------------------------------------" << endl;
}

void Date::idatymd(const int ndat, int& year, int& month, int& day) {
    // date = year*10000+month*100+day;
    year = ndat / 10000;
    month = (ndat - (year * 10000)) / 100;
    day = ndat - (year * 10000) - (month * 100);
}

int Date::get_nbstot(
    const int ndat000, const int ndatfin, const int jday_start,
    const int deltaT, const int date_mode, ivector& rundates) {
    int day = 1;
    int month = 1;
    int year = 1900;
    int jday = 1;
    int newyear = 1900;
    int t_count = 0;
    int ndatcur = ndat000;
    rundates.initialize();

    // CHECK when spinup is fixed that rundates (used for fishing data)
    // do not include spinup!
    // we will need one more date for fishing data!
    while (ndatcur < ndatfin + deltaT) {
        t_count++;
        // idatymd(ndatcur, year, month, day);
        update_time_variables(
            t_count, deltaT, date_mode, jday_start, jday, day, month, year,
            newyear);
        ndatcur = year * 10000 + month * 100 + day;
        rundates[t_count - 1] = ndatcur;
    }
    return t_count - 1;
}

int Date::get_nbt_before_first_recruitment(
    const int first_recruitment_date, const int ndatini, const int deltaT,
    const int date_mode) {
    int day, month, year, jday, lastday;
    int t_count = 0;
    int ndatcur = ndatini;
    int jday_run = 1;

    idatymd(ndatini, year, month, day);

    switch (date_mode) {
        case 1: {
            jday_run = julday(day, month, year);
            break;
        }

        case 2: {
            jday_run = nlyjulday(day, month, year);
            break;
        }

        case 3: {
            jday_run = clmjulday(day, month, year);
            break;
        }
    }

    while (ndatcur <= first_recruitment_date) {
        t_count++;
        update_time_variables(
            t_count + 1, deltaT, date_mode, jday_run, jday, day, month, year,
            lastday);
        ndatcur = year * 10000 + month * 100 + day;
    }
    return t_count - 1;
}

void Date::update_time_variables(
    const int t_count, const int deltaT, const int date_mode,
    const int jday_spinup, int& jday, int& day, int& month, int& year,
    int& newyear) {
    //-------------------------------------------------------------------
    //  Update variables related to NOW time step :
    //   - jday_step
    //   - jday = number of the day in the current year
    //   - day
    //   - month
    //   - year
    //   - newyear (used for RESTART only!)
    //
    //   MODIFICATIONS:
    //   --------------
    //      Julien (May 2008)
    //      Julien (Sep 2008) : add of newyear variable
    //-------------------------------------------------------------------

    int oldyear = year;

    int jday_step = jday_spinup + (t_count - 1) * deltaT;
    switch (date_mode) {
        case 1: {
            dmy(jday_step, day, month, year);
            jday = jday_step - julday(1, 1, year) + 1;
            break;
        }

        case 2: {
            nlydmy(jday_step, day, month, year);
            jday = jday_step - nlyjulday(1, 1, year) + 1;
            break;
        }

        case 3: {
            clmdmy(jday_step, day, month, year);
            jday = jday_step - clmjulday(1, 1, year) + 1;
            break;
        }
    }

    // Test : are we in a new year ?
    if (oldyear == year)
        newyear = 0;
    else
        newyear = 1;
}

int Date::Update_now_time(int yr, int month, int day) {
    //-------------------------------------------------------------------
    //   Return "now" time as an integer
    //   - now_time
    //   as a fonction of :
    //   - day
    //   - month
    //   - year
    //
    //   Ej : 19990204
    //
    //   MODIFICATIONS:
    //   --------------
    //      Julien (December 2008)
    //-------------------------------------------------------------------
    int now_time = yr * 10000 + month * 100 + day;
    return now_time;
}
string Date::Update_now_time_str(int yr, int month, int day) {
    //-------------------------------------------------------------------
    //   Return "now" time string
    //   - now_time
    //   as a fonction of :
    //   - day
    //   - month
    //   - year
    //
    //   Ej : 19990204
    //
    //   MODIFICATIONS:
    //   --------------
    //      Julien (June 2008)
    //-------------------------------------------------------------------

    int ow_time = yr * 10000 + month * 100 + day;
    ostringstream osdate;
    osdate << ow_time;
    string now_time = osdate.str();
    return now_time;
}

string Date::Update_now_time_str_spinup(int month) {
    //-------------------------------------------------------------------
    //   Return "now" time string for monthly climatological spinup
    //
    //   Ej : M02
    //
    //   MODIFICATIONS:
    //   --------------
    //      Julien (June 2008)
    //-------------------------------------------------------------------

    ostringstream osdate;
    osdate << setfill('0') << setw(2) << month;
    string now_time = "M" + osdate.str();
    return now_time;
}

int Date::dym_startdate_run(
    CParam& param, const dvector zlevel_dym, const int nbstot) {
    float dateini = 1900.0;
    int day, month, year;
    int dt = param.deltaT;
    // cout << param.date_mode << endl;
    // need to convert ndatini to float
    switch (param.date_mode) {
        case 1: {
            idatymd(param.ndatini, year, month, day);
            int jdayini = juldayy(day, month, year);
            // temporal: until the fix of date in GLORYS 7days DYM files
            // (action: Beatriz) dateini = year+(jdayini-dt)/365.0;
            dateini = year + jdayini / 365.0;
            if (leapYear(year)) dateini = year + jdayini / 366.0;
            break;
        }
        case 2: {
            idatymd(param.ndatini, year, month, day);
            int jdayini = nlyjuldayy(day, month, year);
            dateini = year + jdayini / 365.0;
            break;
        }
        case 3: {
            idatymd(param.ndatini, year, month, day);
            int jdayini = clmjuldayy(day, month, year);
            // dateini = year+jdayini/360.0;
            dateini = year + jdayini / 365.0;
            if (leapYear(year)) dateini = year + jdayini / 366.0;
            break;
        }
    }
    int nini = 0;
    while (dateini > zlevel_dym[nini] + 1e-4 * dt) {
        ++nini;
    }
    // nini = 0;
    cout << __FILE__ << " : Time steps to skip in DYM files " << nini
         << endl;  // exit(1);
    return nini;
}

/*
void Date::zlevel_run(CParam& param, const dvector zlevel_dym, const int nbstot,
dvector& zlevel, const int nini)
{
//note: the case with forecast is no longer implemented!
//note: this function is used for writing only.
//note: it uses dates from input DYMs (zlevel_dym)
        int deltaT = param.deltaT;
        int dt = (zlevel_dym[1]-zlevel_dym[0])*365.25;


        int t_series = 0;
        for (int n=0; n<nbstot; n++){
                zlevel[n] = zlevel_dym[n+nini];
        }
//cout << dt << " " << deltaT << " " << zlevel << endl; //exit(1);
}
*/

// Leap year or not ?
int Date::leapYear(int year) {
    return (year & 3) == 0 && (year % 100 != 0 || year % 400 == 0);
}

// Is a day (1-31) within a given month?
int Date::dayWithinMonth(int day, int month, int year) {
    int dh;
    if (day <= 0 || month < 1 || month > 12) return 0;
    dh = daysInMonth[month - 1];
    if (leapYear(year) && month == 2) dh++;
    return day <= dh;
}

// Convert Gregorian calendar date to the corresponding Julian day
// Not valid before sep 14 1752
unsigned long Date::julday(int day, int month, int year) {
    unsigned long c, ya, jd;
    if (!dayWithinMonth(day, month, year)) return 0UL;
    if (month > 2) {
        month -= 3; /* wash out the leap day */
    } else {
        month += 9;
        year--;
    }
    c = year / 100;
    ya = year - 100 * c;
    jd = ((146097 * c) >> 2) + ((1461 * ya) >> 2) + (153 * month + 2) / 5 +
         day + 1721119;
    return jd;
}

// compute jday within a year: standard date
unsigned long Date::juldayy(int day, int month, int year) {
    unsigned int jd;

    jd = day;
    for (int i = 0; i < month - 1; i++) {
        jd += daysInMonth[i];
        if (i == 1 && leapYear(year)) jd++;
    }
    return jd;
}

// Convert Gregorian calendar date to the corresponding Julian day for a 360
// days calendar
unsigned long Date::clmjulday(int day, int month, int year) {
    unsigned long c, ya, jd;
    c = year / 100;
    ya = year - 100 * c;
    jd = ((36000 * c)) + ((360 * ya)) + (30 * month) + day;
    return jd;
}

// compute jday within a year: 360-days calendar
unsigned long Date::clmjuldayy(int day, int month, int year) {
    unsigned long jd;
    jd = (30 * (month - 1)) + day;
    return jd;
}

// Convert Gregorian calendar date to the corresponding Julian day for a noleap
// year calendar
unsigned long Date::nlyjulday(int day, int month, int year) {
    unsigned long c, ya, jd;
    c = year / 100;
    ya = year - 100 * c;
    int nd = 0;
    for (int i = 0; i < (month - 1); i++) nd += daysInMonth[i];
    jd = ((36500 * c)) + ((365 * ya)) + nd + day;
    return jd;
}

// compute jday within a year: a noleap year calendar
unsigned long Date::nlyjuldayy(int day, int month, int year) {
    unsigned long jd = day;
    for (int i = 0; i < (month - 1); i++) jd += daysInMonth[i];
    return jd;
}

// Convert a Julian day number to its corresponding Gregorian calendar
void Date::dmy(unsigned long julnum, int& d, int& m, int& y) {
    unsigned long dh;
    unsigned long jh = julnum - 1721119;
    y = ((jh << 2) - 1) / 146097;
    jh = (jh << 2) - 1 - 146097 * y;
    dh = jh >> 2;
    jh = ((dh << 2) + 3) / 1461;
    dh = (dh << 2) + 3 - 1461 * jh;
    dh = (dh + 4) >> 2;
    m = (5 * dh - 3) / 153;
    dh = 5 * dh - 3 - 153 * m;
    d = (dh + 5) / 5;
    y = 100 * y + jh;
    if (m < 10)
        m += 3;
    else {
        m -= 9;
        y++;
    }
}

// Convert a Julian day number to its corresponding 360 days calendar  //
// Modifications Olga
void Date::clmdmy(unsigned long julnum, int& d, int& m, int& y) {
    unsigned long dh;
    unsigned long jh = julnum;
    y = (jh) / 36000;
    jh = (jh)-36000 * y;
    dh = jh;
    jh = ((dh) / 360);
    dh = (dh)-360 * jh;
    m = (dh) / 30;
    if (m == 0) {
        dh = dh - 30 * (m);
        d = dh;
        m = 12;
        jh = jh - 1;
    } else {
        dh = dh - 30 * (m);
        d = dh;
    }
    y = 100 * y + jh;
}

// Convert a Julian day number to its corresponding noleap year day
void Date::nlydmy(unsigned long julnum, int& d, int& m, int& y) {
    unsigned long dh;
    unsigned long jh = julnum;
    y = (jh) / 36500;
    jh = (jh)-36500 * y;
    dh = jh;
    jh = ((dh) / 365);
    dh = (dh)-365 * jh;
    unsigned int mm = 0;
    unsigned int nm = 0;
    while (nm < dh) {
        nm += daysInMonth[mm];
        mm++;
    }
    m = mm;
    d = daysInMonth[m - 1] - (nm - dh);
    y = 100 * y + jh;
}
/*
        string Date::MakeDate(int yr, int mo, int jr)
        {
                ostringstream osyr;
                ostringstream osjr;
                osyr << yr;
                osjr << jr;
                string date= MonthName(mo)+" "+osjr.str()+", "+osyr.str()+"  ";
                return date;
        }

        string Date::MakeDate(int yr, int mo)
        {
                ostringstream osyr;
                osyr << yr;
                string date= MonthName(mo)+", "+osyr.str();
                return date;
        }

        string Date::MonthName(int mo)
        {
                switch (mo)
                {
                case 1:
                        return "Jan";
                        break;
                case 2:
                        return "Feb";
                        break;
                case 3:
                        return "Mar";
                        break;
                case 4:
                        return "Apr";
                        break;
                case 5:
                        return "May";
                        break;
                case 6:
                        return "Jun";
                        break;
                case 7:
                        return "Jul";
                        break;
                case 8:
                        return "Aug";
                        break;
                case 9:
                        return "Sep";
                        break;
                case 10:
                        return "Oct";
                        break;
                case 11:
                        return "Nov";
                        break;
                case 12:
                        return "Dec";
                        break;
                }
return "";
        }
*/
