// Person.cpp

#include <cstdlib>
#include <cstring>
#include <climits>

#include <iostream>
#include <string>
#include <math.h>
#include <algorithm>

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


Infection& Person::initializeNewInfection(Serotype serotype, int time, int sourceloc, int sourceid) {
    Infection& infection = initializeNewInfection(serotype);
    infection.infectedTime  = time;
    infection.infectedPlace = sourceloc;
    infection.infectedByID  = sourceid; // TODO - What kind of ID is this?
    infection.infectiousTime = Parameters::sampler(INCUBATION_CDF, gsl_rng_uniform(RNG)) + time;
    return infection;
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
    return isSusceptible(serotype) and // is susceptible to this serotype (i.e., not immune to this serotype via previous infection)
           !isCrossProtected(time) and // not cross-serotype protection from last infection
           !isVaccineProtected(serotype, time); // not vaccine protected at this time
}


double Person::remainingEfficacy(const int time) const {
    double remainingFraction = 1.0;
    if (not isVaccinated()) {
        remainingFraction = 0.0;
    } else {
        if (_par->linearlyWaningVaccine) {
            // reduce by fraction of immunity duration that has waned
            int time_since_vac = daysSinceVaccination(time);
            if (time_since_vac > _par->vaccineImmunityDuration) {
                remainingFraction = 0.0;
            } else {
                remainingFraction -= ((double) time_since_vac) / _par->vaccineImmunityDuration;
            }
        }
    }
    return remainingFraction;
}


double Person::vaccineProtection(const Serotype serotype, const int time) const {
    double ves;
    if (not isVaccinated()) {
        ves = 0.0;
    } else {
        if (daysSinceVaccination(time) > _par->vaccineImmunityDuration) {
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

enum MaternalEffect { MATERNAL_PROTECTION, NO_EFFECT, MATERNAL_ENHANCEMENT };

MaternalEffect _maternal_antibody_effect(Person* p, const Parameters* _par, int time) {
    MaternalEffect effect = NO_EFFECT;
    if (p->getAge() == 0 and time >= 0) {                       // this is an infant, and we aren't reloading an infection history
        Person* mom = p->getLocation(HOME_NIGHT)->findMom();    // find a cohabitating female of reproductive age
        if (mom and mom->getImmunityBitset().any()) {           // if there is one and she has an infection history
            if (gsl_rng_uniform(RNG) < _par->infantImmuneProb) {
                effect = MATERNAL_PROTECTION;
            } else if (gsl_rng_uniform(RNG) < _par->infantSevereProb) {
                effect = MATERNAL_ENHANCEMENT;
            }
        }
    }

    return effect;
}


// infect - infect this individual
// primary symptomatic is a scaling factor for pathogenicity of primary infections.
// if secondaryPathogenicityOddsRatio > 1, secondary infections are more often symptomatic
// returns true if infection occurs
bool Person::infect(int sourceid, Serotype serotype, int time, int sourceloc) {
    // Bail now if this person can not become infected
    // TODO - clarify this.  why would a person not be infectable in this scope?
    if (not isInfectable(serotype, time)) return false;

    MaternalEffect maternal_effect = _maternal_antibody_effect(this, _par, time);
    bool maternalAntibodyEnhancement;
    switch( maternal_effect ) {
        case MATERNAL_PROTECTION:
            return false;
            break;
        case NO_EFFECT:
            maternalAntibodyEnhancement = false;
            break;
        case MATERNAL_ENHANCEMENT:
            maternalAntibodyEnhancement = true;
            break;
        default:
            cerr << "ERROR: Unknown maternal effect: " << maternal_effect << endl;
            exit(-837);
            break;
    }

    const int numPrevInfections = getNumInfections(); // needs to be called before initializing new infection
    const double remaining_efficacy = remainingEfficacy(time);  // before initializing new infection

    // Create a new infection record
    Infection& infection = initializeNewInfection(serotype, time, sourceloc, sourceid);

    double symptomatic_probability = _par->serotypePathogenicityRelativeRisks[(int) serotype] * _par->basePathogenicity;
    double severe_given_case = 0.0;
    switch (numPrevInfections) {
        case 0:
            if (_par->primaryPathogenicityModel == CONSTANT_PATHOGENICITY) {
                symptomatic_probability *= _par->primaryRelativeRisk;
            } else if (_par->primaryPathogenicityModel == ORIGINAL_LOGISTIC) {
                symptomatic_probability *= SYMPTOMATIC_BY_AGE[_nAge];
            } else if (_par->primaryPathogenicityModel == GEOMETRIC_PATHOGENICITY) {
                symptomatic_probability *= 1.0 - pow(1.0 - _par->annualFlavivirusAttackRate, getAge());
            }
            severe_given_case        = _par->primarySevereFraction[(int) serotype];
            break;
        case 1:
            //symptomatic_probability -- no change
            severe_given_case        = _par->secondarySevereFraction[(int) serotype];
            break;
        case 2:
            symptomatic_probability *= _par->postSecondaryRelativeRisk;
            severe_given_case        = _par->tertiarySevereFraction[(int) serotype];
            break;
        case 3:
            symptomatic_probability *= _par->postSecondaryRelativeRisk;
            severe_given_case        = _par->quaternarySevereFraction[(int) serotype];
            break;
        default:
            cerr << "ERROR: Unsupported number of previous infections: " << numPrevInfections << endl;
            exit(-838);
    }

    if (symptomatic_probability > 1.0) symptomatic_probability = 1.0;
    const double effective_VEP = isVaccinated() ? _par->fVEP*remaining_efficacy : 0.0;        // reduced symptoms due to vaccine
    symptomatic_probability *= (1.0 - effective_VEP);
    assert(symptomatic_probability >= 0.0);
    assert(symptomatic_probability <= 1.0);

    infection.recoveryTime = infection.infectiousTime + INFECTIOUS_PERIOD_ASYMPTOMATIC;            // may be changed below 

    if ((gsl_rng_uniform(RNG) < symptomatic_probability) or maternalAntibodyEnhancement) {         // Is this a case?
        const double severe_rand = gsl_rng_uniform(RNG);
        infection.recoveryTime = infection.infectiousTime + INFECTIOUS_PERIOD_MILD;                // may yet be changed below 
        if ( severe_rand < severe_given_case or maternalAntibodyEnhancement) {                     // Is this a severe case?
            if (not isVaccinated() or gsl_rng_uniform(RNG) > _par->fVEH*remaining_efficacy) { // Is this person unvaccinated or vaccinated but unlucky?
                infection.recoveryTime = infection.infectiousTime + INFECTIOUS_PERIOD_SEVERE;
                infection.severeDisease = true;
            }
        }

        // Determine if this person withdraws (stops going to work/school)
        infection.symptomTime = infection.infectiousTime + SYMPTOMATIC_DELAY;
        const int symptomatic_duration = infection.recoveryTime - infection.symptomTime;
        const int symptomatic_active_period = gsl_ran_geometric(RNG, 0.5) - 1; // min generator value is 1 trial
        infection.withdrawnTime = symptomatic_active_period < symptomatic_duration ?
                                  infection.symptomTime + symptomatic_active_period :
                                  infection.withdrawnTime;
    }

    // Flag locations with (non-historical) infections, so that we know to look there for human->mosquito transmission
    // Negative days are historical (pre-simulation) events, and thus we don't care about modeling transmission
    for (int day = std::max(infection.infectiousTime, 0); day < infection.recoveryTime; day++) {
        for (int t=0; t<(int) NUM_OF_TIME_PERIODS; t++) {
            Community::flagInfectedLocation(_pLocation[t], day);
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
    return !_bDead && !(_nImmunity[serotype] == 1);
}


bool Person::isCrossProtected(int time) const {
    return (getNumInfections() > 0) and // has any past infection
           (infectionHistory.back()->infectedTime + _par->nDaysImmune > time); // prev. infection w/in crossprotection period 
}


bool Person::isVaccineProtected(Serotype serotype, int time) const {
    return isVaccinated() and
           ( !_par->bVaccineLeaky or // if the vaccine isn't leaky
            (gsl_rng_uniform(RNG) < vaccineProtection(serotype, time)) ); // or it protects (i.e., doesn't leak this time)
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
