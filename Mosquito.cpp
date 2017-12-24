// Mosquito.cpp

#include <cstdlib>
#include <cstring>
#include <climits>
#include <cmath>
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



Mosquito::Mosquito(Location* p, Serotype serotype, int nInfectedAtID, int nExternalIncubationPeriod, double prob_infecting_bite) {
    _nID = _nNextID++;
    _bDead = false;
    _eSerotype = serotype;
    _nInfectedAtID = nInfectedAtID;
    // extract precalculated age CDF given the specified prob_infecting_bite
    vector<double> age_cdf = MOSQUITO_FIRST_BITE_AGE_CDF_MESH[(int) (prob_infecting_bite * (MOSQUITO_FIRST_BITE_AGE_CDF_MESH.size()-1))]; //MOSQUITO_AGE_CDF;
    _nAgeInfected = Parameters::sampler(age_cdf, gsl_rng_uniform(RNG));
    _nAgeInfectious = _nAgeInfected + nExternalIncubationPeriod;
    _nAgeDeath = _nAgeInfected; // can't be younger than this
    double r = 1.0-(gsl_rng_uniform(RNG)*(1.0-MOSQUITO_DEATHAGE_CDF[_nAgeInfected]));
    _nAgeDeath = Parameters::sampler(MOSQUITO_DEATHAGE_CDF, r, _nAgeDeath);
    //cerr << _nAgeInfected << " " << _nAgeDeath << endl;
    _pLocation = _pOriginLocation = p;
    _pLocation->addInfectedMosquito();
}


Mosquito::Mosquito(RestoreMosquitoPars* rp):
    _pLocation(rp->location), _eSerotype(rp->serotype), _nAgeInfected(rp->age_infected), _nAgeInfectious(rp->age_infectious), _nAgeDeath(rp->age_dead) {
    _nID = _nNextID++;
    _bDead = false;
    _nInfectedAtID = -1;
    _pOriginLocation = NULL;
    _pLocation->addInfectedMosquito();
}


Mosquito::~Mosquito() {
    _pLocation->removeInfectedMosquito();
}
