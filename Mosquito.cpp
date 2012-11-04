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
    while (MOSQUITO_AGE_DISTRIBUTION[_nAgeInfected]<r && _nAgeInfected<MAXMOSQUITOAGE-1)
        _nAgeInfected++;
    _nAgeDeath = _nAgeInfected;
    r = 1.0-(gsl_rng_uniform(RNG)*(1.0-MOSQUITO_AGE_DISTRIBUTION[_nAgeDeath]));
    while (MOSQUITO_AGE_DISTRIBUTION[_nAgeDeath]<r && _nAgeDeath<MAXMOSQUITOAGE-1)
        _nAgeDeath++;
    //  cerr << "Mosquito " << _nID << ", age=" << _nAgeInfected << ", death=" << _nAgeDeath << endl;
    _pLocation = _pOriginLocation = p;
}


Mosquito::~Mosquito() {
}
