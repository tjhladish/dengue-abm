#ifndef __DATE_H
#define __DATE_H
#include "Parameters.h"

enum DayOfWeekType {
    MONDAY,
    TUESDAY,
    WEDNESDAY,
    THURSDAY,
    FRIDAY,
    SATURDAY,
    SUNDAY,
    NUM_OF_DAYS_IN_WEEK
};

static const std::vector<std::string> DAY_NAMES = {"MON", "TUE", "WED", "THU", "FRI", "SAT", "SUN"};
static const std::vector<std::string> MONTH_NAMES = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};
static const std::vector<size_t> COMMON_DAYS_IN_MONTH = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
static const std::vector<size_t> COMMON_END_DAY_OF_MONTH = {31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};
static const std::vector<size_t> LEAP_DAYS_IN_MONTH = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
static const std::vector<size_t> LEAP_END_DAY_OF_MONTH = {31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366};

struct TimeSeriesAnchorPoint {
    string date;
    double value;
};

class Date {
  public:
    Date():_simulation_day(0),_offset(0),_month_ct(0),_year_ct(0),_julian_day(1),_julian_year(1) {};
    Date(const Parameters* par):_simulation_day(0),_offset(par->startDayOfYear-1),_month_ct(0),_year_ct(0),_julian_day((size_t) par->startDayOfYear),_julian_year((size_t) par->startJulianYear) {};
    Date(const Date& o): _offset(o._offset) {
        _simulation_day = o._simulation_day;
        _month_ct       = o._month_ct;
        _year_ct        = o._year_ct;
        _julian_day     = o._julian_day;
        _julian_year    = o._julian_year; 
    };

    // for use when needing to pass a historical sim day as a Date object, e.g. to model pre-existing natural immunity.  Note that month and year counts are not set!
    Date(const Parameters* par, int sim_day):_simulation_day(sim_day),_offset(par->startDayOfYear-1),_month_ct(0),_year_ct(0),_julian_day(par->startDayOfYear),_julian_year(par->startJulianYear) {};

    template <typename T> friend inline bool operator== (const Date &lhs, const T rhs)  { return lhs.to_ymd().compare(rhs) == 0; }
    template <typename T> friend inline bool operator== (const T lhs, const Date &rhs)  { return rhs.to_ymd().compare(lhs) == 0; }

    template <typename T> friend inline bool operator!= (const Date &lhs, const T rhs)  { return lhs.to_ymd().compare(rhs) != 0; }
    template <typename T> friend inline bool operator!= (const T lhs, const Date &rhs)  { return rhs.to_ymd().compare(lhs) != 0; }

    template <typename T> friend inline bool operator<  (const Date &lhs, const T rhs)  { return lhs.to_ymd().compare(rhs) <  0; }
    template <typename T> friend inline bool operator<  (const T lhs, const Date &rhs)  { return rhs.to_ymd().compare(lhs) >  0; }

    template <typename T> friend inline bool operator>  (const Date &lhs, const T rhs)  { return lhs.to_ymd().compare(rhs) >  0; }
    template <typename T> friend inline bool operator>  (const T lhs, const Date &rhs)  { return rhs.to_ymd().compare(lhs) <  0; }

    template <typename T> friend inline bool operator<= (const Date &lhs, const T rhs)  { return lhs.to_ymd().compare(rhs) <= 0; }
    template <typename T> friend inline bool operator<= (const T lhs, const Date &rhs)  { return rhs.to_ymd().compare(lhs) >= 0; }

    template <typename T> friend inline bool operator>= (const Date &lhs, const T rhs)  { return lhs.to_ymd().compare(rhs) >= 0; }
    template <typename T> friend inline bool operator>= (const T lhs, const Date &rhs)  { return rhs.to_ymd().compare(lhs) <= 0; }

    int offset()                const { return _offset; }
    inline int day()            const { return _simulation_day; }                                   // [0, ...]
    size_t julianDay()          const { return _julian_day; }                                       // [1, {365, 366}]
    static size_t dayOfMonth(size_t julian_day, size_t julian_year) {                               // [1, {29,30,31}]
                                        const size_t jm = julianMonth(julian_day, julian_year);
                                        return jm == 1 ? julian_day : julian_day - end_day_of_month(julian_year)[jm - 2];
    }
    size_t dayOfMonth()         const { return dayOfMonth(julianDay(), julianYear()); }

    size_t nDayPeriod(int n)    const { return (int) day()/n; }                                     // [0, ...]
    size_t week()               const { return (int) day()/NUM_OF_DAYS_IN_WEEK; }                   // [0, ...]
    static size_t julianWeek(size_t julian_day) { return ((julian_day-1)/NUM_OF_DAYS_IN_WEEK) + 1; }// [1, 53]
    size_t julianWeek()         const { return julianWeek(julianDay()); }

    size_t month()              const { return _month_ct; }                                         // [0, ...]
    static size_t julianMonth(size_t julian_day, size_t julian_year) {                              // [1, 12]
        vector<size_t>::const_iterator it;
        // find first month that hasn't ended (hint: it's this month)
        // julianDay()-1 because this isn't upper_or_equal_bound, which would be convenient
        const vector<size_t> end_dom = end_day_of_month(julian_year);
        it = upper_bound(end_dom.begin(), end_dom.end(), julian_day-1);
        return it - end_dom.begin() + 1; // +1 b/c [1, 12], not [0, 11]
    }

    size_t julianMonth()     const {                                                                // [1, 12]
        return julianMonth(julianDay(), julianYear());
    }

    static DayOfWeekType dayOfWeek(int y, int m, int d) {
        // John Conway's "Doomsday Algorithm"
        // static int t[] = {0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4}; // week starting Sunday
        static int t[] = {6, 2, 1, 4, 6, 2, 4, 0, 3, 5, 1, 3};   // week starting Monday
        y -= m < 3;
        return (DayOfWeekType) ((y + y/4 - y/100 + y/400 + t[m-1] + d) % 7);
    }

    DayOfWeekType dayOfWeek() const { return dayOfWeek(julianYear(), julianMonth(), dayOfMonth()); }
    bool isWeekend() const { return dayOfWeek() == SATURDAY or dayOfWeek() == SUNDAY; }
    bool isWeekday() const { return not isWeekend(); }

    string dayOfWeekName(int y, int m, int d) const { return DAY_NAMES[dayOfWeek(y, m, d)]; }
    string dayOfWeekName()   const { return DAY_NAMES[dayOfWeek()]; }
    string monthName()       const { return MONTH_NAMES[julianMonth()-1]; }
    size_t year()            const { return _year_ct; }

    bool endOfPeriod(int n)  const { return (day()+1) % n == 0; }
    bool endOfWeek()         const { return (day()+1) % 7 == 0; }
    bool endOfMonth()        const {
        vector<size_t>::const_iterator it;
        // find out if today corresponds to a month-end
        const vector<size_t> end_dom = end_day_of_month();
        it = find(end_dom.begin(), end_dom.end(), julianDay());
        return it != end_dom.end();
    }
    static vector<size_t> end_day_of_month(size_t year) { return isLeap(year) ? LEAP_END_DAY_OF_MONTH : COMMON_END_DAY_OF_MONTH; }
    //vector<size_t> end_day_of_month() const { return isLeap() ? LEAP_END_DAY_OF_MONTH : COMMON_END_DAY_OF_MONTH; }
    vector<size_t> end_day_of_month() const { return end_day_of_month(julianYear()); }
    bool startOfYear()       const { return startOfJulianYear(); }
    bool endOfYear()         const { return endOfJulianYear(); }                               // is it end of {365,366} day period
    bool startOfJulianYear() const { return julianDay() == 1; }                                // is it Jan 1
    bool endOfJulianYear()   const { return julianDay() == num_days_in_year(julianYear()); }   // is it Dec 31

    void setJulianDay (int jd) { _julian_day = jd; }
    void setJulianYear (int julian_year) { _julian_year = julian_year; }

    size_t julianYear () const { return _julian_year; }

    void increment() {
        if(endOfMonth()) _month_ct++;
        if(endOfYear())  _year_ct++;
        _simulation_day++;
        _julian_day++;
        if (_julian_day > num_days_in_year(_julian_year)) {
            _julian_day = 1;
            _julian_year++;
        }
    }

    void print() {
        cerr << day() << "\t" << julianDay() << "\t" << year()
             << "\t(" << monthName() << " " << dayOfMonth() << ")\t" << month() << "\t" << julianMonth();
        if (endOfWeek()) cerr << " EoW";
        if (endOfMonth()) cerr << " EoM";
        if (startOfYear()) cerr << " SoY";
        if (endOfYear()) cerr << " EoY";
        if (endOfJulianYear()) cerr << " EoJY";
        cerr << endl;
    }

    bool isLeap() const {
        return Date::isLeap(julianYear());
    }

    static bool isLeap (size_t julian_year) {
        if (julian_year % 4 != 0) {
            return false;
        } else if (julian_year % 100 != 0) {
            return true;
        } else if (julian_year % 400 != 0) {
            return false;
        } else {
            return true;
        }
    }

    static size_t num_days_in_year(size_t julian_year) {
        if (isLeap(julian_year)) {
            return LEAP_END_DAY_OF_MONTH.back();
        } else {
            return COMMON_END_DAY_OF_MONTH.back();
        }
    }

    // generic function to get date strings that look like 2020-03-15, 03/15/2020, 315 (but why?), etc.
    string to_string(vector<string> format, string sep="/", bool zero_pad = true) const {
        stringstream ss;
        for (size_t i = 0; i<format.size(); ++i) {
            string part = format[i];
            if (zero_pad) { ss << setfill('0'); }
            if (part == "yyyy") {
                ss << setw(4) << julianYear();
            } else if (part == "mm") {
                ss << setw(2) << julianMonth();
            } else if (part == "dd") {
                ss << setw(2) << dayOfMonth();
            } else {
                cerr << "ERROR: Date::date_string() format argument \"" << part << "\" not supported.  Please use yyyy, mm, and/or dd.\n";
                exit (-1);
            }
            if (i < format.size() - 1) { ss << sep; }
        }
        return ss.str();
    }

    // assumes yyyy-mm-dd format
    static vector<size_t> to_numeric_parts(string date_string){
        vector<string> parts = dengue::util::split(date_string, '-');
        assert(parts.size() == 3);
        size_t julian_year  = (size_t) stoi(parts[0]);
        size_t julian_month = (size_t) stoi(parts[1]);
        size_t day_of_month = (size_t) stoi(parts[2]);
        return {julian_year, julian_month, day_of_month};
    }

    // assumes yyyy-mm-dd format
    static size_t to_julian_day(const vector<size_t>& date_parts) {
        const size_t julian_year  = date_parts[0];
        const size_t julian_month = date_parts[1];
        const size_t day_of_month = date_parts[2];
        assert(julian_month >= 1);
        assert(julian_month <= 12);
        if (julian_month == 1) {
            return day_of_month;
        } else if (isLeap(julian_year)) {
            return LEAP_END_DAY_OF_MONTH[julian_month - 2] + day_of_month;
        } else {
            return COMMON_END_DAY_OF_MONTH[julian_month - 2] + day_of_month;
        }
    }

    // assumes yyyy-mm-dd format
    static size_t to_julian_day(string date_string) {
        vector<size_t> parts = to_numeric_parts(date_string);
        return to_julian_day(parts);
    }

    // assumes yyyy-mm-dd format
    static string to_ymd (size_t julian_year, int julian_day) {
        string sep = "-";
        stringstream ss;
        ss << setfill('0')
           << setw(4) << julian_year << sep
           << setw(2) << julianMonth(julian_day, julian_year) << sep
           << setw(2) << dayOfMonth(julian_day, julian_year);
        return ss.str();
    }

    string to_ymd() const {
        return to_ymd(julianYear(), julianDay());
    }

    static int to_sim_day(size_t julian_start_year, size_t julian_start_day, string date_str) {
        // TODO -- currently assumes that requested date_str does not precede the julian start day/year
        const vector<size_t> parts = to_numeric_parts(date_str);
        const size_t date_year = parts[0];
        const size_t julian_day = to_julian_day(parts);
        int day_tally = 0;
        for (size_t y = julian_start_year; y <= date_year; ++y) {
            day_tally += num_days_in_year(y);
        }
        day_tally -= julian_start_day;
        day_tally -= num_days_in_year(date_year) - julian_day;
        return day_tally;
    }

    static vector<double> linInterpolate(double start, double end, size_t n) {
        // "end" is the value immediately after the sequence produced here
        // e.g. linInterpolate(0.0, 1.0, 5) produces {0.0, 0.2, 0.4, 0.6, 0.8}
        assert(n>=1);
        //cout << "Start: " << start << "; End: " << end << endl;

        // initialize the n values to the start value
        double step = (end - start)/n;
        vector<double> v(n, start);
        for (size_t i = 1; i < n; i++) { v[i] += i * step; }

        return v;
    }

    //static vector<double> linInterpolateTimeSeries(vector<TimeSeriesAnchorPoint> ap, bool truncateStart, size_t julian_start_year, size_t julian_start_day) {
    static vector<double> linInterpolateTimeSeries(vector<TimeSeriesAnchorPoint> ap, size_t julian_start_year, size_t julian_start_day) {
        if (ap.size() < 2) {
            cerr << "Must have at least two anchor points." << endl;
            exit(1);
        }

        // sort anchor points based on date
        sort(ap.begin(), ap.end(), [](TimeSeriesAnchorPoint &a, TimeSeriesAnchorPoint &b){return a.date < b.date;});
        vector<double> v;

        for (size_t i = 0; i < ap.size() - 1; i++) {
            const TimeSeriesAnchorPoint A = ap[i];
            const TimeSeriesAnchorPoint B = ap[i+1];
            size_t n = to_sim_day(julian_start_year, julian_start_day, B.date) - to_sim_day(julian_start_year, julian_start_day, A.date);

            vector<double> v_i = linInterpolate(A.value, B.value, n);
            v.insert(v.end(), v_i.begin(), v_i.end());
        }
        v.push_back(ap.back().value); // bc linInterpolate leaves off last value, to prevent repetitions

        int series_julian_start = to_sim_day(julian_start_year, julian_start_day, ap[0].date);

        // TODO - how could the julian start be less than 0?
        //if (truncateStart and series_julian_start < 0) {
        if (series_julian_start < 0) {
            v.erase(v.begin(), v.begin()-series_julian_start);
        }

        return v;
    }


  private:
    int _simulation_day;
    const int _offset;
    size_t _month_ct;
    size_t _year_ct;
    size_t _julian_day;
    size_t _julian_year;
};


#endif
