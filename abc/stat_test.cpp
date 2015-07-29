#include <iostream>
#include <vector>
#include "AbcUtil.h"

using namespace std;
// To compile on hipergator:
// g++ -O2 -std=c++0x -Wall -I/scratch/lfs/thladish/AbcSmc -I$HPC_GSL_INC ../../AbcSmc/AbcUtil.o stat_test.cpp -o stats -L/scratch/lfs/thladish/AbcSmc -lm -L$HPC_GSL_LIB/ -lgsl -lgslcblas
// or
// mpicxx -O2 -std=c++11 -w0 -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX -D USING_MPI -I/scratch/lfs/thladish/AbcSmc -I$HPC_GSL_INC ../../AbcSmc/AbcUtil.o stat_test.cpp -o stats -L/scratch/lfs/thladish/AbcSmc -lm -L$HPC_GSL_LIB/ -lgsl -lgslcblas

int main() {
    // yucatan dengue case data, 1979-2013
    vector<float_type> cases = {                                                      4234, //1979
                                4672, 3377, 1412,  643, 5495,  193,   34,   15,  356,    2, // 1980-1989
                                   8, 352,    22,   29,  680,   69,  650, 5529,   36,   43, // 1990-1999
                                   0, 287,   946,   26,   57,  162,  627, 1861,  721, 3212, // 2000-2009
                                2517, 6132, 13000, 8937};                                   // 2010-2013

    // yucatan population size, interpolated as necessary, 1979-2012
    vector<float_type> pop = {                                                                                 1033195, // 1979
                              1063733, 1093654, 1123574, 1153495, 1183416, 1213337, 1243257, 1273178, 1303099, 1333019, // 1980-1989
                              1362940, 1401676, 1440413, 1479149, 1517886, 1556622, 1576940, 1597257, 1617575, 1637892, // 1990-1999
                              1658210, 1690358, 1722505, 1754653, 1786800, 1818948, 1846274, 1873600, 1900925, 1928251, // 2000-2009
                              1955577, 1982903, 2010229, 2037554};                                                      // 2010-2013

    //for (unsigned int y = 0; y<cases.size() and y<pop.size(); ++y) cout << 1979+y << " " << 1e5*cases[y]/pop[y] << endl;
    vector<float_type> incidence;
    for (unsigned int y = 0; y<cases.size() and y<pop.size(); ++y) incidence.push_back(1e5*cases[y]/pop[y]);


    Map<Col> col(&incidence[0], incidence.size()); // copy data from vector into Col
    cout << col.transpose() << endl;
    // vector<float_type> new_vec(col.data(), col.data()+col.size());

    cout << "mean: " << mean(col) << endl;
    cout << "median: " << median(col) << endl;
    cout << "stdev: " << sqrt(variance(col, mean(col))) << endl;
    cout << "max: " << max(col) << endl;
    cout << "skewness: " << skewness(col) << endl;
    cout << "median crossings: " << median_crossings(col) << endl;
    
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
