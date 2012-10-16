// Mosquito.cpp

#include <cstdlib>
#include <cstring>
#include <climits>
#include <iostream>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Mosquito.h"
#include "Parameters.h"

using namespace dengue::standard;

int Mosquito::_nNextID = 0;

//double Mosquito::_fDeathProbability[MAXMOSQUITOAGE] =
//  {0.0018,0.002069514,0.002378646,0.002732982,0.003138823,0.003603251,0.004134191,0.004740476,0.005431895,0.006219231,0.007114273,0.008129798,0.009279508,0.01057792,0.01204018,0.01368180,0.01551828,0.01756467,0.01983494,0.02234135,0.02509359,0.02809799,0.03135655,0.03486617,0.03861782,0.04259607,0.04677877,0.05113718,0.05563645,0.06023658,0.0648937,0.0695617,0.07419404,0.07874562,0.08317441,0.08744304,0.09151983,0.09537952,0.09900354,0.1023799,0.1055030,0.1083723,0.1109923,0.1133714,0.1155206,0.1174531,0.1191837,0.1207277,0.1221007,0.1233179,0.1243943,0.1253438,0.1261799,0.1269147,0.1275594,0.1281243,0.1286187,0.1290510,0.1294285,0.1297580,0.1300453}; // probability of death each day
double Mosquito::_fAgeDistribution[MAXMOSQUITOAGE] = {                        // cumulative density of mosquito ages
    0.03115129,0.06224651,0.09327738,0.1242344,0.1551069,0.1858824,0.2165471,
    0.247085,0.2774781,0.3077061,0.3377462,0.3675725,0.3971563,0.4264656,0.4554649,
    0.484115,0.5123732,0.5401928,0.5675238,0.5943126,0.620503,0.6460362,0.6708519,
    0.6948895,0.718089,0.7403926,0.7617461,0.7821008,0.8014145,0.8196537,0.8367943,
    0.8528225,0.8677358,0.8815426,0.8942622,0.9059238,0.9165657,0.9262337,0.9349795,
    0.9428595,0.9499327,0.9562597,0.961901,0.9669162,0.9713627,0.9752957,0.9787666,
    0.981824,0.9845121,0.9868721,0.988941,0.9907526,0.9923371,0.9937217,0.9949305,
    0.9959852,0.9969047,0.997706,0.9984038,1.0
};

//double Mosquito::_fDailyBitingProbs[3] = {0.7, 0.3, 1.0}; // probability of biting at 3 different times of day

Mosquito::Mosquito() {
    _nID = _nNextID++;
    _bDead = false;
    //_bInfected = false;
    _nAgeInfected = -1;
    _eSerotype = NULL_SEROTYPE;
    _nInfectedAtID = -1;
    _pLocation = _pOriginLocation = NULL;
}


Mosquito::Mosquito(Location *p, Serotype serotype, int nInfectedAtID) {
    //assert(nSerotype>0 && nSerotype<=4);
    _nID = _nNextID++;
    _bDead = false;
    _eSerotype = serotype;
    _nInfectedAtID = nInfectedAtID;
    double r = gsl_rng_uniform(RNG);
    _nAgeInfected = 0;
    while (_fAgeDistribution[_nAgeInfected]<r && _nAgeInfected<MAXMOSQUITOAGE-1)
        _nAgeInfected++;
    _nAgeDeath = _nAgeInfected;
    r = 1.0-(gsl_rng_uniform(RNG)*(1.0-_fAgeDistribution[_nAgeDeath]));
    while (_fAgeDistribution[_nAgeDeath]<r && _nAgeDeath<MAXMOSQUITOAGE-1)
        _nAgeDeath++;
    //  cerr << "Mosquito " << _nID << ", age=" << _nAgeInfected << ", death=" << _nAgeDeath << endl;
    _pLocation = _pOriginLocation = p;
}


Mosquito::~Mosquito() {
}
