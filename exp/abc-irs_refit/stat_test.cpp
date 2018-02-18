#include <iostream>
#include <vector>
#include "AbcUtil.h"
//#include "ranker.h"

using namespace std;
// To compile on hipergator:
// g++ -O2 -std=c++0x -Wall -I../../AbcSmc -I$HPC_GSL_INC ../../AbcSmc/AbcUtil.o stat_test.cpp -o stats -lm -L$HPC_GSL_LIB/ -lgsl -lgslcblas
// or
// mpicxx -O2 -std=c++11 -w0 -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX -D USING_MPI -I/scratch/lfs/thladish/AbcSmc -I$HPC_GSL_INC ../../AbcSmc/AbcUtil.o stat_test.cpp -o stats -L/scratch/lfs/thladish/AbcSmc -lm -L$HPC_GSL_LIB/ -lgsl -lgslcblas

// To compile on dragonfly:
// g++ -O2 -std=c++0x -Wall -I../../AbcSmc ../../AbcSmc/AbcUtil.o stat_test.cpp -o stats -L../../AbcSmc -lm -lgsl -lgslcblas
int main() {
    // yucatan dengue case data, 1979-2013
    vector<int> cases = {                                                       4234, // 1979
                          4672, 3377, 1412,  643, 5495,  193,   34,   15,  356,    2, // 1980-1989
                             8, 352,    22,   29,  680,   69,  650, 5529,   36,   43, // 1990-1999
                             0, 288,   946,   26,   57,  162,  627, 1862,  722, 3222, // 2000-2009
                          2523, 6219, 5723, 2864, 1073, 1505};                        // 2010-2015

    // yucatan population size, interpolated as necessary, 1979-2012
    vector<ABC::float_type> pop = {                                                                            1033195, // 1979
                              1063733, 1093654, 1123574, 1153495, 1183416, 1213337, 1243257, 1273178, 1303099, 1333019, // 1980-1989
                              1362940, 1401676, 1440413, 1479149, 1517886, 1556622, 1576940, 1597257, 1617575, 1637892, // 1990-1999
                              1658210, 1690358, 1722505, 1754653, 1786800, 1818948, 1846274, 1873600, 1900925, 1928251, // 2000-2009
                              1955577, 1982903, 2010229, 2037554, 2064880, 2092206};                                    // 2010-2015




    //for (unsigned int y = 0; y<cases.size() and y<pop.size(); ++y) cout << 1979+y << " " << 1e5*cases[y]/pop[y] << endl;
    vector<ABC::float_type> incidence;
    for (unsigned int y = 0; y<cases.size() and y<pop.size(); ++y) incidence.push_back(1e5*cases[y]/pop[y]);
{
    cerr << "total cases per 100k people\n";
    int year = 1979;
    for (auto val: incidence) cerr << setprecision(1) << fixed << year++ << " " << val << endl;
}

    ABC::Map<ABC::Col> col(&incidence[0], incidence.size()); // copy data from vector into Col
    cout << col.transpose() << endl;
    // vector<float_type> new_vec(col.data(), col.data()+col.size());

    vector<int> severe = {                                                         0,  // 1979
                             0,    0,    0,    0,    9,    0,    0,    0,    0,    0,  // 1980-1989
                             0,    0,    0,    0,    6,    4,   30,  163,    0,    0,  // 1990-1999
                             0,   36,  197,    6,    6,   39,  162,  390,  149, 1113,  // 2000-2009
                           817, 2124, 2824, 1032,  440,  376};                         // 2010-2015

    vector<double> x(severe.size()); for(unsigned int i = 0; i < x.size(); ++i) x[i] = i;

    int y1995_idx = 16;

    //ABC::LogisticFit* fit = ABC::logistic_reg(x, severe, cases);
    //assert(fit->status == GSL_SUCCESS);

    assert(cases.size() == pop.size());
    assert(cases.size() == severe.size());
    assert(cases.size() > y1995_idx);

    double pre_1995_severe = accumulate(severe.begin(),           severe.begin()+y1995_idx, 0.0);
    double modern_severe   = accumulate(severe.begin()+y1995_idx, severe.end(), 0.0);
    double pre_1995_cases  = accumulate(cases.begin(),            cases.begin()+y1995_idx, 0.0);
    double modern_cases    = accumulate(cases.begin()+y1995_idx,  cases.end(), 0.0);
    cerr << pre_1995_severe << " " << modern_severe << endl;
    cerr << pre_1995_cases << " " << modern_cases << endl;
    double mean_pre_1995_severe  = pre_1995_severe / pre_1995_cases;
    double mean_modern_severe    = modern_severe / modern_cases;

    cout << "mean:               " << ABC::mean(col) << endl;
    cout << "0%   quantile:      " << quantile(incidence, 0.0) << endl;
    cout << "25%  quantile:      " << quantile(incidence, 0.25) << endl;
    cout << "50%  quantile:      " << quantile(incidence, 0.5) << endl;
    cout << "75%  quantile:      " << quantile(incidence, 0.75) << endl;
    cout << "100% quantile:      " << quantile(incidence, 1.0) << endl;
    cout << "stdev:              " << sqrt(ABC::variance(col, ABC::mean(col))) << endl;
    cout << "skewness:           " << ABC::skewness(col) << endl;
    cout << "median crossings:   " << ABC::median_crossings(col) << endl;
    cout << "seroprevalence:     " << 0.6 << endl;
    cout << "pre95 severe prev:  " << mean_pre_1995_severe << endl;
    cout << "modern severe prev: " << mean_modern_severe << endl;

//        {"name" : "mean",       "num_type" : "FLOAT",   "value"     : },
//        {"name" : "min",        "num_type" : "FLOAT",   "value"     : },
//        {"name" : "quant25",    "num_type" : "FLOAT",   "value"     : },
//        {"name" : "median",     "num_type" : "FLOAT",   "value"     : },
//        {"name" : "quant75",    "num_type" : "FLOAT",   "value"     : },
//        {"name" : "max",        "num_type" : "FLOAT",   "value"     : },
//        {"name" : "stdev",      "num_type" : "FLOAT",   "value"     : },
//        {"name" : "skewness",   "num_type" : "FLOAT",   "value"     : },
//        {"name" : "med_xing",   "num_type" : "FLOAT",   "value"     : },
//        {"name" : "seroprev",   "num_type" : "FLOAT",   "value"     : },
//        {"name" : "beta0",      "num_type" : "FLOAT",   "value"     : },
//        {"name" : "beta1",      "num_type" : "FLOAT",   "value"     : },


    // test of median_crossings function
    /*
    Col test(4); 
    test << 1,1,1,1; cout << "mc==1: " << median_crossings(test) << endl;
    test << 1,2,2,2; cout << "mc==1: " << median_crossings(test) << endl;
    test << 1,2,2,1; cout << "mc==2: " << median_crossings(test) << endl;
    test << 1,2,1,2; cout << "mc==3: " << median_crossings(test) << endl;
    test.resize(5);
    test << 1,3,3,2,1; cout << "mc==2: " << median_crossings(test) << endl;
    test << 1,3,3,2,3; cout << "mc==2: " << median_crossings(test) << endl;
    test << 1,1,1,1,1; cout << "mc==1: " << median_crossings(test) << endl;
    test << 1,3,1,3,1; cout << "mc==3: " << median_crossings(test) << endl;
    test << 1,3,2,3,1; cout << "mc==3: " << median_crossings(test) << endl;
    test.resize(1);
    test << 1; cout << "mc==0: " << median_crossings(test) << endl;
    test.resize(0);
    cout << "mc==0: " << median_crossings(test) << endl;
    */
}
