// Person.cpp

#include <cstdlib>
#include <cstring>
#include <climits>

#include <iostream>
#include <string>

#include <assert.h>
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
    _nImmunity = 0;
    _nNumInfections = -1;
    pushInfectionHistory();
    for(int i=0; i<STEPSPERDAY; i++) _pLocation[i] = NULL;
    _bDead = false;
    _bVaccinated = false;
    _bCase = false;
}


Person::~Person() {
}


void Person::pushInfectionHistory() {
    _nNumInfections++;
    assert(_nNumInfections<MAXHISTORY);
    for (int i=_nNumInfections; i>0; i--) {
        _nInfectedByID[i] = _nInfectedByID[i-1];
        _nInfectedPlace[i] = _nInfectedPlace[i-1];
        _nInfectedTime[i] = _nInfectedTime[i-1];
        _nInfectiousTime[i] = _nInfectiousTime[i-1];
        _nSymptomTime[i] = _nSymptomTime[i-1];
        _nWithdrawnTime[i] = _nWithdrawnTime[i-1];
        _nRecoveryTime[i] = _nRecoveryTime[i-1];
        _eSerotype[i] = _eSerotype[i-1];
    }
    _eSerotype[0] = NULL_SEROTYPE; 
    _nInfectedByID[0] = -1;
    _nInfectedPlace[0] = -1;
    _nInfectedTime[0] = _nInfectiousTime[0] = _nSymptomTime[0] = _nRecoveryTime[0] = -100000;
    _nWithdrawnTime[0] = 10000;
}


// copyImmunity - copy immune status from person *p
void Person::copyImmunity(const Person *p) {
    assert(p!=NULL);
    _nImmunity = p->_nImmunity;
    _bVaccinated = p->_bVaccinated;
    _nNumInfections = p->_nNumInfections;
    for (int i=0; i<=_nNumInfections; i++) {
        _nInfectedByID[i] = p->_nInfectedByID[i];
        _nInfectedPlace[i] = p->_nInfectedPlace[i];
        _nInfectedTime[i] = p->_nInfectedTime[i];
        _nInfectiousTime[i] = p->_nInfectiousTime[i];
        _nSymptomTime[i] = p->_nSymptomTime[i];
        _nWithdrawnTime[i] = p->_nWithdrawnTime[i];
        _nRecoveryTime[i] = p->_nRecoveryTime[i];
        _eSerotype[i] = p->_eSerotype[i];
    }
}


// resetImmunity - reset immune status (infants)
void Person::resetImmunity() {
    _nImmunity = 0;
    _nNumInfections = -1;
    pushInfectionHistory();
    _bVaccinated = false;
    _bCase = false;
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


// infect - infect this individual
// primary symptomatic is a scaling factor for pathogenicity of primary infections.
// secondaryscaling * primarysymptomatic is the scaling factor for pathogenicity of secondary infections.
// returns true if infection occurs
bool Person::infect(int sourceid, Serotype serotype, int time, int sourceloc) {
    double primarysymptomatic = _par->fPrimaryPathogenicity[(int) serotype];
    double secondaryscaling = _par->fSecondaryScaling[(int) serotype];
    int maxinfectionparity = _par->nMaxInfectionParity;
    //assert(serotype>=1 && serotype<=4);
    if (!isSusceptible(serotype) ||                                   // already infected by serotype
        _nRecoveryTime[0]+_par->nDaysImmune>time ||                        // cross-serotype protection
        getInfectionParity()>=maxinfectionparity)                     // already infected max number of times
        return false;
    pushInfectionHistory();

    _nInfectedTime[0] = time;

    double r = gsl_rng_uniform(RNG);
    _nInfectiousTime[0] = 0;
    while (_nInfectiousTime[0]<MAXINCUBATION && INCUBATION_DISTRIBUTION[_nInfectiousTime[0]]<r)
        _nInfectiousTime[0]++;
    //    cerr << "symp " << _nSymptomTime << "," << r << endl;
    _nInfectiousTime[0] += time;
    _nRecoveryTime[0] = _nInfectiousTime[0]+5;                        // should draw from distribution!!!!!!
    if (_nImmunity>0)
        _nRecoveryTime[0]--;                                          // secondary infections are one day shorter
    _nInfectedByID[0] = sourceid;
    _nInfectedPlace[0] = sourceloc;
    _eSerotype[0] = serotype;

    if (primarysymptomatic>0.0 &&
        gsl_rng_uniform(RNG)<primarysymptomatic*SYMPTOMATIC_BY_AGE[_nAge] * (isVaccinated()?(1.0-_par->fVEP):1.0) *
    (_nImmunity>0?secondaryscaling:1.0)) {                            // scale for primary or secondary infection
        _nSymptomTime[0] = _nInfectiousTime[0] + 1;                   // symptomatic one day before infectious
        double r = gsl_rng_uniform(RNG);
        if (r<0.5) {
            _nWithdrawnTime[0] = _nSymptomTime[0];                    // withdraws (FIX THIS!!!!)
        }
        else if (r<0.75) {
            _nWithdrawnTime[0] = _nSymptomTime[0] + 1;                // withdraws (FIX THIS!!!!)
        }
        else if (r<0.75) {
            _nWithdrawnTime[0] = _nSymptomTime[0] + 2;                // withdraws (FIX THIS!!!!)
        }
        else if (r<0.875) {
            _nWithdrawnTime[0] = _nSymptomTime[0] + 3;                // withdraws (FIX THIS!!!!)
        }
        else if (r<0.945) {
            _nWithdrawnTime[0] = _nSymptomTime[0] + 4;                // withdraws (FIX THIS!!!!)
        }
        if (_nRecoveryTime[0]<_nWithdrawnTime[0])
            _nWithdrawnTime[0] = 100000;
        _bCase = true;
    }
    else {
        _nSymptomTime[0] = -10000;
    }
    _nImmunity |= (1<<((int) serotype));
    //  cerr << _nAge << "," << serotype << ": " <<  primarysymptomatic << "," << secondaryscaling << "," << SYMPTOMATIC_BY_AGE[_nAge] <<  " ; " << _nSymptomTime[0] << endl;

    return true;
}


bool Person::isInfected(int time) {
    if (time>=_nInfectedTime[0] && time<_nRecoveryTime[0])
        return true;
    else
        return false;
}


bool Person::isViremic(int time) {
    if (time>=_nInfectiousTime[0] && time<_nRecoveryTime[0] && !_bDead)
        return true;
    else
        return false;
}


bool Person::isSymptomatic(int time) {
    if (time>=_nSymptomTime[0] && time<_nRecoveryTime[0] && !_bDead)
        return true;
    else
        return false;
}


bool Person::isWithdrawn(int time) {
    if (time>=_nWithdrawnTime[0] && time<_nRecoveryTime[0] && !_bDead)
        return true;
    else
        return false;
}


bool Person::isSusceptible(Serotype serotype) {
    return !_bDead && !(_nImmunity & (1<<((int) serotype)));
}


// getInfectionParity - number of serotypes with immunity (including vaccination)
int Person::getInfectionParity() {
    return ((_nImmunity & (1<<(1-1))) +
        ((_nImmunity & (1<<(2-1)))>>1) +
        ((_nImmunity & (1<<(3-1)))>>2) +
        ((_nImmunity & (1<<(4-1)))>>3));
}


bool Person::vaccinate() {
    vector<double> _fVES = _par->fVESs;
    if (!_bVaccinated & !_bDead) {
        _bVaccinated = true;
        if ((_fVES[0]==_fVES[1]) &&
            (_fVES[1]==_fVES[2]) &&
        (_fVES[2]==_fVES[3])) {                                       // same protection against all 4 serotypes
            if (gsl_rng_uniform(RNG)<_fVES[0]) {
                _nImmunity = 0xff;                                    // this person is protected against all serotypes
            }
        }
        else {
            if (gsl_rng_uniform(RNG)<_fVES[0])
                _nImmunity |= (1<<(1-1));                             // protected against serotype 1
            if (gsl_rng_uniform(RNG)<_fVES[1])
                _nImmunity |= (1<<(2-1));                             // protected against serotype 2
            if (gsl_rng_uniform(RNG)<_fVES[2])
                _nImmunity |= (1<<(3-1));                             // protected against serotype 3
            if (gsl_rng_uniform(RNG)<_fVES[3])
                _nImmunity |= (1<<(4-1));                             // protected against serotype 4
        }
        return true;
    }
    return false;
}


// generateRandomIDs - generate num sorted unique ints between 0 and bound-1
void Person::generateRandomIDs(int num, int bound, int *ids) {
    for (int i=0; i<num; i++) {
        int r = gsl_rng_uniform_int(RNG, bound-i);
        int j;
        for (j=0; j<i; j++) {
            if (r>=ids[j])
                r++;
            else
                break;
        }
        for (int k=i-1; k>=j; k--)
            ids[k+1]=ids[k];                                          // move entry over
        ids[j]=r;
    }
}
