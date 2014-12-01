// Person.cpp

#include <cstdlib>
#include <cstring>
#include <climits>

#include <iostream>
#include <string>
#include <math.h>

#include <assert.h>
#include <bitset>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Person.h"
#include "Community.h"
#include "Parameters.h"

using namespace dengue::standard;

int Person::_nNextID = 1;

const Parameters* Person::_par;

Person::Person() {
    _nID = _nNextID++;
    _nAge = -1;
    _nLifespan = -1;
    _nHomeID = -1;
    _nWorkID = -1;
    for(int i=0; i<STEPS_PER_DAY; i++) _pLocation[i] = NULL;
    _bDead = false;
    _bVaccinated = false;
    _bNaiveVaccineProtection = false;
    _bCase = false;
}


Person::~Person() {
    clearInfectionHistory();
}

void Person::clearInfectionHistory() {
    for (unsigned int i = 0; i < infectionHistory.size(); i++) {
        delete infectionHistory[i];
    }
    infectionHistory.clear();
}

Infection& Person::initializeNewInfection(Serotype serotype) {
    setImmunity(serotype);
    Infection* infection = new Infection(serotype);
    infectionHistory.push_back(infection);
    return *infection;
}


// copyImmunity - copy immune status from person *p
void Person::copyImmunity(const Person *p) {
    assert(p!=NULL);
    _nImmunity = p->_nImmunity;
    _bVaccinated = p->_bVaccinated;
    clearInfectionHistory();
    for (int i=0; i < p->getNumInfections(); i++) {
        infectionHistory.push_back( new Infection(p->infectionHistory[i]) );
    }
}


// resetImmunity - reset immune status (infants)
void Person::resetImmunity() {
    _nImmunity.reset();
    clearInfectionHistory();
    _bVaccinated = false;
    _bNaiveVaccineProtection = false;
    _bCase = false;
    _bDead = false;
}


bool Person::naturalDeath(int t) {
    if (_nLifespan<=_nAge+(t/365.0)) {
        _bDead = true;
        return true;
    }
    return false;
}


void Person::kill(int time) {
    _bDead = true;
}


bool Person::isInfectable(Serotype serotype, int time) {
    bool infectable = true;
    const int maxinfectionparity    = _par->nMaxInfectionParity;
    if (!isSusceptible(serotype)) {
        // immune to this serotype via previous infection OR non-leaky vaccine
        infectable = false;                                   
    } else if (getNumInfections() > 0 and infectionHistory.back()->recoveryTime + _par->nDaysImmune > time) {
        // cross-serotype protection from last infection
        infectable = false;
    } else if (getInfectionParity() >= maxinfectionparity) {
        // already infected max number of times
        infectable = false;
    } else if (isVaccinated() && _par->bVaccineLeaky==true) {
        // potentially protected by leaky vaccine
        // (all-or-none vaccines are handled via isSusceptible conditional)
        double r = gsl_rng_uniform(RNG);
        double protectionProbability = 0;
        // Which level of protection does this person have?
        if (_bNaiveVaccineProtection == true) {
            protectionProbability = _par->fVESs_NAIVE[serotype];
        } else {
            protectionProbability = _par->fVESs[serotype];
        }
        if ( r < protectionProbability ) {
            infectable = false;
        }
    }
    return infectable;
}


// infect - infect this individual
// primary symptomatic is a scaling factor for pathogenicity of primary infections.
// secondary_scaling * primarysymptomatic is the scaling factor for pathogenicity of secondary infections.
// returns true if infection occurs
bool Person::infect(int sourceid, Serotype serotype, int time, int sourceloc) {
    // Bail now if this person can not become infected
    if (not isInfectable(serotype, time)) {
        return false;
    }

    Infection& infection = initializeNewInfection(serotype); // Create a new infection record

    infection.infectedByID  = sourceid; // TODO - What kind of ID is this?
    infection.infectedTime  = time;
    infection.infectedPlace = sourceloc;

    // When do they become infectious?
    infection.infectiousTime = 0;
    double r = gsl_rng_uniform(RNG);
    while (infection.infectiousTime < MAX_INCUBATION && INCUBATION_DISTRIBUTION[infection.infectiousTime] < r) {
        infection.infectiousTime++;
    }
    infection.infectiousTime += time;

    if (_nImmunity.any()) {
        infection.recoveryTime = infection.infectiousTime+INFECTIOUS_PERIOD_SEC;        // should draw from distribution!!!!!!
    } else {
        infection.recoveryTime = infection.infectiousTime+INFECTIOUS_PERIOD_PRI;        // should draw from distribution!!!!!!
    }

    // Determine if this person withdraws (stops going to work/school)
    const double primary_symptomatic   = _par->fPrimaryPathogenicity[(int) serotype];
    const double secondary_scaling     = _par->fSecondaryScaling[(int) serotype];
    const double vaccine_protection    = isVaccinated()   ? _par->fVEP        : 0.0;   // reduced symptoms due to vaccine
    const double secondary_symptomatic = _nImmunity.any() ? secondary_scaling : 1.0;

    const double symptomatic_probability = primary_symptomatic * SYMPTOMATIC_BY_AGE[_nAge] * (1.0 - vaccine_protection) * secondary_symptomatic;

    if (gsl_rng_uniform(RNG) < symptomatic_probability) {
        infection.symptomTime = infection.infectiousTime + 1;                   // symptomatic one day after infectious
        double r = gsl_rng_uniform(RNG);
        const int symptomatic_duration = infection.recoveryTime - infection.symptomTime;
        for (int i=0; i<symptomatic_duration; i++) {
            if (r < 1-pow(0.5,1+i) ) {
                infection.withdrawnTime = infection.symptomTime + i;                // withdraws
                break;
            }
        }
        _bCase = true;
    }

    for (int d=infection.infectiousTime; d<infection.recoveryTime; d++) {
        for (int t=0; t<STEPS_PER_DAY; t++) {
            Community::flagInfectedLocation(_pLocation[t], d);
        }
    }
    if (_par->bRetroactiveMatureVaccine) {
        // if the best vaccine-induced immunity can be acquired retroactively,
        // upgrade this person from naive to mature
        _bNaiveVaccineProtection = false;
    }
    return true;
}


bool Person::isNewlyInfected(int time) {
    if (infectionHistory.size() > 0) {
        Infection* infection = infectionHistory.back();
        if (time == infection->infectedTime) {
            return true;
        }
    }
    return false;
}


bool Person::isInfected(int time) {
    if (infectionHistory.size() > 0) {
        Infection* infection = infectionHistory.back();
        if (time >= infection->infectedTime and time < infection->recoveryTime) {
            return true;
        }
    }
    return false;
}


bool Person::isViremic(int time) {
    if (infectionHistory.size() > 0) {
        Infection* infection = infectionHistory.back();
        if (time >= infection->infectiousTime and time < infection->recoveryTime and not _bDead) {
            return true;
        }
    }
    return false;
}


bool Person::isSymptomatic(int time) {
    if (infectionHistory.size() > 0) {
        Infection* infection = infectionHistory.back();
        if (time >= infection->symptomTime and time < infection->recoveryTime and not _bDead) {
            return true;
        }
    }
    return false;
}


bool Person::isWithdrawn(int time) {
    if (infectionHistory.size() > 0) {
        Infection* infection = infectionHistory.back();
        if (time >= infection->withdrawnTime and time < infection->recoveryTime and not _bDead) {
            return true;
        }
    }
    return false;
}


bool Person::isSusceptible(Serotype serotype) {
    return !_bDead && !(_nImmunity[serotype] == 1);  //1<<0=1  1<<1=2  1<<2=4  1<<3=8   1<<4=16
}


// getInfectionParity - number of serotypes with immunity (including vaccination)
int Person::getInfectionParity() {
    return _nImmunity.count(); // count number of 1's in bit string
}


bool Person::fullySusceptible() {
    bool susceptible = true;
    for (int s = 0; s<(int) NUM_OF_SEROTYPES; ++s) {
        if ( not isSusceptible((Serotype) s) ) { susceptible = false; }
    }
    return susceptible;
}


bool Person::vaccinate() {
    if (!_bDead) {
        //vector<double> _fVES = _par->fVESs;
        _bVaccinated = true;
        if ( fullySusceptible() ) {
            _bNaiveVaccineProtection = true;
        } else {
            _bNaiveVaccineProtection = false;
        }
        if ( _par->bVaccineLeaky == false ) { // all-or-none VE_S protection
            if ( fullySusceptible() ) { // naive against all serotypes
                for (int i=0; i<NUM_OF_SEROTYPES; i++) {
                    if (gsl_rng_uniform(RNG)<_par->fVESs_NAIVE[i]) _nImmunity[i] = 1;                                // protect against serotype i
                }
            } else {
                for (int i=0; i<NUM_OF_SEROTYPES; i++) {
                    if (gsl_rng_uniform(RNG)<_par->fVESs[i]) _nImmunity[i] = 1;                                // protect against serotype i
                }
            }
        }
        return true;
    } else {
        return false;
    }
}

