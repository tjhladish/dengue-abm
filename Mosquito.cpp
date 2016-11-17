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


Mosquito::Mosquito(Location* p, Serotype serotype, int nInfectedAtID, int nExternalIncubationPeriod, int time) {
    _nID = _nNextID++;
    _bDead = false;
    _eSerotype = serotype;
    _nInfectedAtID = nInfectedAtID;
    //const double rho = p->getCurrentVectorControlDailyMortality(time);
    vector<double> age_cdf = MOSQUITO_AGE_CDF;
    // no longer adjusting CDF -- mosquitoes at location already have increased mortality in Community::applyVectorControl()
    /*if (rho > 0) {
                                               //mortality_by_age_i* = 1 - (1 - mortality_by_age_i)(1 - rho)^i
        for (unsigned int i = 0; i < age_cdf.size(); ++i) age_cdf[i] = 1.0 - (1.0 - MOSQUITO_AGE_CDF[i]) * pow(1.0 - rho, i);
    }*/
    _nAgeInfected = Parameters::sampler(age_cdf, gsl_rng_uniform(RNG));
    _nAgeInfectious = _nAgeInfected + nExternalIncubationPeriod;
    _nAgeDeath = _nAgeInfected; // can't be younger than this
    double r = 1.0-(gsl_rng_uniform(RNG)*(1.0-age_cdf[_nAgeInfected]));
    _nAgeDeath = Parameters::sampler(age_cdf, r, _nAgeDeath);

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
