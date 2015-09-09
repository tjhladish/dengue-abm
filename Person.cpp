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
    for(int i=0; i<(int) NUM_OF_TIME_PERIODS; i++) _pLocation[i] = NULL;
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


// copyImmunity - copy immune status from person* p
void Person::copyImmunity(const Person* p) {
    assert(p!=NULL);
    _nImmunity = p->_nImmunity;
    _bVaccinated = p->_bVaccinated;

    vaccineHistory.clear();
    vaccineHistory.assign(p->vaccineHistory.begin(), p->vaccineHistory.end());

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
    vaccineHistory.clear();
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


bool Person::isInfectable(Serotype serotype, int time) const {
    bool infectable = true;
    const int maxinfectionparity = _par->nMaxInfectionParity;
    if (!isSusceptible(serotype)) {
        // immune to this serotype via previous infection OR non-leaky vaccine
        infectable = false;                                   
    } else if (getNumInfections() > 0 and infectionHistory.back()->infectedTime + _par->nDaysImmune > time) {
        // cross-serotype protection from last infection
        infectable = false;
    } else if (getInfectionParity() >= maxinfectionparity) {
        // already infected max number of times
        infectable = false;
    } else if (isVaccinated() && _par->bVaccineLeaky==true) {
        // potentially protected by leaky vaccine
        // (all-or-none vaccines are handled via isSusceptible conditional)
        double r = gsl_rng_uniform(RNG);
        double protectionProbability = vaccineProtection(serotype, time);
        // Which level of protection does this person have?
        if ( r < protectionProbability ) {
            infectable = false;
        }
    } else {
        // Apparently person is infectable
    }
    return infectable;
}


double Person::remainingEfficacy(const int time) const {
    double remainingFraction = 1.0;
    if (_par->linearlyWaningVaccine) {
        // reduce by fraction of immunity duration that has waned
        int time_since_vac = daysSinceVaccination(time);
        remainingFraction -=  ((double) time_since_vac / _par->vaccineImmunityDuration); 
    }
    return remainingFraction;
}


double Person::vaccineProtection(const Serotype serotype, const int time) const {
    double ves;
    if (not isVaccinated()) {
        ves = 0.0;
    } else {
        int time_since_vac = daysSinceVaccination(time);
        if (time_since_vac > _par->vaccineImmunityDuration) {
            ves = 0.0;
        } else {
            if (_bNaiveVaccineProtection == true) {
                ves = _par->fVESs_NAIVE[serotype];
            } else {
                ves = _par->fVESs[serotype];
            }
            ves *= remainingEfficacy(time);
        }
    }
    return ves;
}


// infect - infect this individual
// primary symptomatic is a scaling factor for pathogenicity of primary infections.
// if secondaryPathogenicityOddsRatio > 1, secondary infections are more often symptomatic
// returns true if infection occurs
bool Person::infect(int sourceid, Serotype serotype, int time, int sourceloc) {
    // Bail now if this person can not become infected
    // TODO - clarify this.  why would a person not be infectable in this scope?
    if (not isInfectable(serotype, time)) {
        return false;
    }

    bool maternalAntibodyEnhancement = false;
    if (getAge() == 0 and time >= 0) {                       // this is an infant, and we aren't reloading an infection history
        Person* mom = getLocation(HOME_NIGHT)->findMom();    // find a cohabitating female of reproductive age
        if (mom and mom->getImmunityBitset().any()) {        // if there is one and she has an infection history
            if (gsl_rng_uniform(RNG) < _par->infantImmuneProb) {
                return false;                                // maternal antibody protection
            } else if (gsl_rng_uniform(RNG) < _par->infantSevereProb) {
                maternalAntibodyEnhancement = true;          // maternal antibodies are a problem
            }
        }
    }

    Infection& infection = initializeNewInfection(serotype); // Create a new infection record

    infection.infectedByID  = sourceid; // TODO - What kind of ID is this?
    infection.infectedTime  = time;
    infection.infectedPlace = sourceloc;

    const bool postprimary = getImmunityBitset().any();      // has this person ever been infected before?

    // When do they become infectious?
    infection.infectiousTime = 0;
    double r = gsl_rng_uniform(RNG);
    while (infection.infectiousTime < MAX_INCUBATION && INCUBATION_DISTRIBUTION[infection.infectiousTime] < r) {
        infection.infectiousTime++;
    }
    infection.infectiousTime += time;

    if (postprimary) {
        infection.recoveryTime = infection.infectiousTime+INFECTIOUS_PERIOD_SEC;        // TODO - draw from distribution
    } else {
        infection.recoveryTime = infection.infectiousTime+INFECTIOUS_PERIOD_PRI;        // TODO - draw from distribution
    }

    const double primary_symptomatic_prob = _par->primaryPathogenicity[(int) serotype] * SYMPTOMATIC_BY_AGE[_nAge];
    double symptomatic_probability = primary_symptomatic_prob;

    if (postprimary and primary_symptomatic_prob != 1.0) {
        const double secondary_odds = _par->secondaryPathogenicityOddsRatio[(int) serotype] * primary_symptomatic_prob / (1.0 - primary_symptomatic_prob);
        symptomatic_probability = secondary_odds / (1.0 + secondary_odds);
    }

    const double effective_VEP = isVaccinated() ? _par->fVEP : 0.0;   // reduced symptoms due to vaccine
    symptomatic_probability *= (1.0 - effective_VEP);

    if (gsl_rng_uniform(RNG) < symptomatic_probability or maternalAntibodyEnhancement) {
        // Person develops symptoms --> this is a case
        // Is this a severe case?
        const double severe_rand = gsl_rng_uniform(RNG);
        if ( (!postprimary and severe_rand < _par->primarySevereFraction[(int) serotype])     // primary infection
             or (postprimary and severe_rand < _par->secondarySevereFraction[(int) serotype]) // secondary infection
             or maternalAntibodyEnhancement) {                                                // infant w/ antibodies
            if (not isVaccinated() or gsl_rng_uniform(RNG) > _par->fVEH*remainingEfficacy(time)) { // is this person unvaccinated or unlucky?
                infection.severeDisease = true;
            }
        }

        // Determine if this person withdraws (stops going to work/school)
        // TODO - the 1 in the below line is a parameter and should be moved out
        infection.symptomTime = infection.infectiousTime + 1;                   // symptomatic one day after infectious
        double r = gsl_rng_uniform(RNG);
        const int symptomatic_duration = infection.recoveryTime - infection.symptomTime;
        for (int i=0; i<symptomatic_duration; i++) {
            if (r < 1-pow(0.5,1+i) ) {                                              // TODO - refactor to sample from geometric
                infection.withdrawnTime = infection.symptomTime + i;                // withdraws
                break;
            }
        }
        _bCase = true;
    }

    // TODO - clarify what's going on with d; sometimes this is an old infection, sometimes a current one
    for (int d=infection.infectiousTime; d<infection.recoveryTime; d++) {
        if (d < 0) continue; // we don't need to worry about flagging locations for past infections
        for (int t=0; t<(int) NUM_OF_TIME_PERIODS; t++) {
            Community::flagInfectedLocation(_pLocation[t], d);
        }
    }

    // if the antibody-primed vaccine-induced immunity can be acquired retroactively, upgrade this person from naive to mature
    if (_par->bRetroactiveMatureVaccine) _bNaiveVaccineProtection = false;

    return true;
}


bool Person::isNewlyInfected(int time) const {
    if (infectionHistory.size() > 0) {
        Infection* infection = infectionHistory.back();
        if (time == infection->infectedTime) {
            return true;
        }
    }
    return false;
}


bool Person::isInfected(int time) const {
    if (infectionHistory.size() > 0) {
        Infection* infection = infectionHistory.back();
        if (time >= infection->infectedTime and time < infection->recoveryTime) {
            return true;
        }
    }
    return false;
}


bool Person::isViremic(int time) const {
    if (infectionHistory.size() > 0) {
        Infection* infection = infectionHistory.back();
        if (time >= infection->infectiousTime and time < infection->recoveryTime and not _bDead) {
            return true;
        }
    }
    return false;
}


bool Person::isSymptomatic(int time) const {
    if (infectionHistory.size() > 0) {
        Infection* infection = infectionHistory.back();
        // TODO: this is a mess of conditionals.  Make it not confusing.
        if (infection->isSymptomatic() and time >= infection->symptomTime and time < infection->recoveryTime and not _bDead) {
            return true;
        }
    }
    return false;
}


bool Person::hasSevereDisease(int time) const {
    if (infectionHistory.size() > 0) {
        Infection* infection = infectionHistory.back();
        if (infection->severeDisease and time >= infection->symptomTime and time < infection->recoveryTime and not _bDead) {
            return true;
        }
    }
    return false;
}


bool Person::isWithdrawn(int time) const {
    if (infectionHistory.size() > 0) {
        Infection* infection = infectionHistory.back();
        if (time >= infection->withdrawnTime and time < infection->recoveryTime and not _bDead) {
            return true;
        }
    }
    return false;
}


bool Person::isSusceptible(Serotype serotype) const {
    return !_bDead && !(_nImmunity[serotype] == 1);  //1<<0=1  1<<1=2  1<<2=4  1<<3=8   1<<4=16
}


// getInfectionParity - number of serotypes with immunity (including vaccination)
int Person::getInfectionParity() const {
    return _nImmunity.count(); // count number of 1's in bit string
}


bool Person::fullySusceptible() const {
    bool susceptible = true;
    for (int s = 0; s<(int) NUM_OF_SEROTYPES; ++s) {
        if ( not isSusceptible((Serotype) s) ) { susceptible = false; }
    }
    return susceptible;
}


bool Person::vaccinate(int time) {
    if (!_bDead) {
        //vector<double> _fVES = _par->fVESs;
        _bVaccinated = true;
        vaccineHistory.push_back(time);
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

