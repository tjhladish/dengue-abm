#ifndef IMM_GEN_H
#define IMM_GEN_H

#include <vector>
#include <map>
#include <string>

using std::vector;
using std::map;
using std::string;

template <typename T> inline T sum(vector<T> list) { T sum=0; for (unsigned int i=0; i<list.size(); i++) sum += list[i]; return sum;}

const double EPSILON = 10e-15;

const int NUM_SEROTYPES = 4;
const int MAX_CENSUS_AGE = 85;// used if a census age group has only a minimum value, e.g. '85+'

// Reported DF + DHF cases, 1997-2011 (inclusive)

const vector<int> CASES =                                 {4234, // 1979
     4672, 3377, 1412,  643, 5495,  193,   34,   15,  356,    2, // 1980-1989
        8,  352,   22,   29,  680,   69,  650, 5529,   36,   43, // 1990-1999
        0,  287,  946,   26,   57,  162,  627, 1861,  721, 3212, // 2000-2009
     2517, 6132, 5705};                                          // 2010-2012

const float FRACTION_SEROTYPED = 0.056; // based on Casos Históricos.xlsx from Hector
 
const vector< vector<double> > SEROTYPE_WT = { 
    {1.00,  0.00,  0.00,  0.00},   // 1979 // Fig. 1
    {1.00,  0.00,  0.00,  0.00},   // 1980 // Fig. 1
    {1.00,  0.00,  0.00,  0.00},   // 1981 // Fig. 1
    {1.00,  0.00,  0.00,  0.00},   // 1982 // Fig. 1
    {1.00,  0.00,  0.00,  0.00},   // 1983 // Fig. 1
    {0.50,  0.00,  0.00,  0.50},   // 1984 // Fig. 1
    {1.00,  0.00,  0.00,  0.00},   // 1985 // Fig. 1
    {0.00,  1.00,  0.00,  0.00},   // 1986 // Fig. 1
    {1.00,  0.00,  0.00,  0.00},   // 1987 // Fig. 1
    {1.00,  0.00,  0.00,  0.00},   // 1988 // Fig. 1
    {1.00,  0.00,  0.00,  0.00},   // 1989 // Fig. 1
    {1.00,  0.00,  0.00,  0.00},   // 1990 // Fig. 1
    {0.50,  0.50,  0.00,  0.00},   // 1991 // Fig. 1
    {1.00,  0.00,  0.00,  0.00},   // 1992 // Fig. 1
    {1.00,  0.00,  0.00,  0.00},   // 1993 // Fig. 1
    {0.33,  0.33,  0.00,  0.34},   // 1994 // Fig. 1
    {0.50,  0.50,  0.00,  0.00},   // 1995 // Fig. 1
    {0.25,  0.25,  0.25,  0.25},   // 1996 // Fig. 1
    {0.09,  0.00,  0.87,  0.04},   // 1997 // Fig. 4
    {0.09,  0.00,  0.87,  0.04},   // 1998 // extrapolated
    {0.09,  0.00,  0.87,  0.04},   // 1999 // extrapolated
    {0.09,  0.00,  0.87,  0.04},   // 2000 // extrapolated
    {0.00,  0.60,  0.40,  0.00},   // 2001 // Fig. 4
    {0.04,  0.96,  0.00,  0.00},   // 2002 // Fig. 4
    {0.04,  0.96,  0.00,  0.00},   // 2003 // extrapolated
    {0.04,  0.96,  0.00,  0.00},   // 2004 // extrapolated
    {0.11,  0.89,  0.00,  0.00},   // 2005 // Fig. 4
    {0.27,  0.55,  0.18,  0.00},   // 2006 // Fig. 4
    {0.90,  0.04,  0.04,  0.02},   // 2007 // Fig. 4
    {0.85,  0.15,  0.00,  0.00},   // 2008 // Fig. 4
    {0.46,  0.54,  0.00,  0.00},   // 2009 // Fig. 4
    {0.59,  0.41,  0.00,  0.00},   // 2010 // Fig. 4
    {0.32,  0.68,  0.00,  0.00},   // 2011 // Fig. 4
    {0.34,  0.66,  0.00,  0.00}};  // 2012 // Fig. 4


vector<int> bootstrap_cases(const double EF);

vector< vector<double> > bootstrap_serotypes(const vector<int>& cases); 

long seedgen(void);

int sample_serotype(vector<double> dist);

struct AgeTally {
    AgeTally(string a, int t) { age=a; tally=t; }
    string age;
    int tally;
};


void split(const string& s, char c, vector<string>& v);

map<int,vector<AgeTally> > import_census_data(string filename);

bool recently_infected(vector<int> states);

int age_str_to_int(string age_str);

vector< vector< vector<int> > > initialize_full_population(map<int,vector<AgeTally> >& census, int first_year);

void age_immunity(vector< vector< vector<int> > >& pop);

int extract_tally(vector<AgeTally> census_year, string age_str);

void age_full_population(map<int,vector<AgeTally> >& census, vector< vector< vector<int> > >& full_pop, int new_year);

vector<int> tally_pop_by_age(vector<vector<vector<int > > > full_pop);

void output_immunity_file(string filename, const vector<vector<vector<int> > >& full_pop);

vector<vector<vector<int>>> simulate_immune_dynamics(const float EXPANSION_FACTOR, const int LAST_YEAR);

#endif
