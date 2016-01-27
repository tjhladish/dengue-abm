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
                             0, 287,   946,   26,   57,  162,  627, 1861,  721, 3212, // 2000-2009
                          2517, 6132, 13000, 8937};                                   // 2010-2013

    // yucatan population size, interpolated as necessary, 1979-2012
    vector<ABC::float_type> pop = {                                                                                 1033195, // 1979
                              1063733, 1093654, 1123574, 1153495, 1183416, 1213337, 1243257, 1273178, 1303099, 1333019, // 1980-1989
                              1362940, 1401676, 1440413, 1479149, 1517886, 1556622, 1576940, 1597257, 1617575, 1637892, // 1990-1999
                              1658210, 1690358, 1722505, 1754653, 1786800, 1818948, 1846274, 1873600, 1900925, 1928251, // 2000-2009
                              1955577, 1982903, 2010229, 2037554};                                                      // 2010-2013

    //for (unsigned int y = 0; y<cases.size() and y<pop.size(); ++y) cout << 1979+y << " " << 1e5*cases[y]/pop[y] << endl;
    vector<ABC::float_type> incidence;
    for (unsigned int y = 0; y<cases.size() and y<pop.size(); ++y) incidence.push_back(1e5*cases[y]/pop[y]);


    ABC::Map<ABC::Col> col(&incidence[0], incidence.size()); // copy data from vector into Col
    cout << col.transpose() << endl;
    // vector<float_type> new_vec(col.data(), col.data()+col.size());

    vector<int> severe = {                                                         0,  // 1979
                             0,    0,    0,    0,    9,    0,    0,    0,    0,    0,  // 1980-1989
                             0,    0,    0,    0,    6,    4,   30,  163,    0,    0,  // 1990-1999
                             0,   35,  197,    6,    6,   39,  162,  389,  148, 1110,  // 2000-2009
                           810, 2092, 2497,  914};                                     // 2010-2013

    vector<double> x(severe.size()); for(unsigned int i = 0; i < x.size(); ++i) x[i] = i;

    ABC::LogisticFit* fit = ABC::logistic_reg(x, severe, cases);
    assert(fit->status == GSL_SUCCESS);

    double total_severe = accumulate(severe.begin(), severe.end(), 0.0);
    double total_cases  = accumulate(cases.begin(), cases.end(), 0.0);
    double mean_severe  = total_severe / total_cases;

    cout << "mean:             " << ABC::mean(col) << endl;
    cout << "0%   quantile:    " << quantile(incidence, 0.0) << endl;
    cout << "25%  quantile:    " << quantile(incidence, 0.25) << endl;
    cout << "50%  quantile:    " << quantile(incidence, 0.5) << endl;
    cout << "75%  quantile:    " << quantile(incidence, 0.75) << endl;
    cout << "100% quantile:    " << quantile(incidence, 1.0) << endl;
    cout << "stdev:            " << sqrt(ABC::variance(col, ABC::mean(col))) << endl;
    cout << "skewness:         " << ABC::skewness(col) << endl;
    cout << "median crossings: " << ABC::median_crossings(col) << endl;
    cout << "seroprevalence:   " << 0.6 << endl;
    cout << "severe prev:      " << mean_severe << endl;
    cout << "beta0:            " << fit->beta0 << endl;
    cout << "beta1:            " << fit->beta1 << endl;


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
