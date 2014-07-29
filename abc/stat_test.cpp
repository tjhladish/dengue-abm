#include <iostream>
#include <vector>
#include "AbcUtil.h"

using namespace std;
// To compile on hipergator:
// g++ -O2 -std=c++0x -Wall -I/scratch/lfs/thladish/AbcSmc -I$HPC_GSL_INC ../../AbcSmc/AbcUtil.o stat_test.cpp -o stats -L/scratch/lfs/thladish/AbcSmc -lm -L$HPC_GSL_LIB/ -lgsl -lgslcblas

int main() {
    // yucatan dengue case data, 1979-2012
    vector<float_type> vdata = {4234, 4672, 3377, 1412, 643, 5495, 193, 34, 15, 356, 2, 8, 352, 22, 29, 680, 69, 650, 5529, 36, 43, 0, 287, 946, 26, 57, 162, 627, 1861, 721, 3212, 2517, 6132, 5705};

    Map<Col> col(&vdata[0], vdata.size()); // copy data from vector into Col
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
