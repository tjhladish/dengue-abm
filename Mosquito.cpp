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
    _nAgeInfected = -1;
    _nAgeInfectious = -1;
    _eSerotype = NULL_SEROTYPE;
    _nInfectedAtID = -1;
    _pLocation = _pOriginLocation = NULL;
}


Mosquito::Mosquito(Location* p, Serotype serotype, int nInfectedAtID, int nExternalIncubationPeriod) {
    _nID = _nNextID++;
    _bDead = false;
    _eSerotype = serotype;
    _nInfectedAtID = nInfectedAtID;
    _nAgeInfected = Parameters::sampler(MOSQUITO_AGE_CDF, gsl_rng_uniform(RNG));
    _nAgeInfectious = _nAgeInfected + nExternalIncubationPeriod;
    _nAgeDeath = _nAgeInfected; // can't be younger than this
    double r = 1.0-(gsl_rng_uniform(RNG)*(1.0-MOSQUITO_AGE_CDF[_nAgeInfected]));
    _nAgeDeath = Parameters::sampler(MOSQUITO_AGE_CDF, r, _nAgeDeath);

    _pLocation = _pOriginLocation = p;
    _pLocation->addInfectedMosquito();
}

Mosquito::Mosquito(Location* p, Serotype s, int ageInfd, int ageInfs, int ageDead):
    _pLocation(p), _eSerotype(s), _nAgeInfected(ageInfd), _nAgeInfectious(ageInfs), _nAgeDeath(ageDead) {
    _nID = _nNextID++;
    _bDead = false;
    _nInfectedAtID = -1;
    _pOriginLocation = NULL;
    _pLocation->addInfectedMosquito();
}


Mosquito::~Mosquito() {
    _pLocation->removeInfectedMosquito();
}
