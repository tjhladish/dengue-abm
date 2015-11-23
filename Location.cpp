// Location.cpp

#include <cstdlib>
#include <cstring>
#include <climits>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Mosquito.h"
#include "Person.h"
#include "Location.h"
#include "Parameters.h"

using namespace dengue::standard;

int Location::_nNextSerial = 0;
//int Location::_nDefaultMosquitoCapacity;

Location::Location()
    : _person((int) NUM_OF_TIME_PERIODS, vector<Person*>(0) ) {
    _serial = _nNextSerial++;
    _ID = 0;
    _nBaseMosquitoCapacity = 0;
    _currentInfectedMosquitoes = 0;
    _coord = make_pair(0.0, 0.0);
}


Location::~Location() {
    _person.clear();
    _neighbors.clear();
}


void Location::addPerson(Person* p, int t) {
    //assert((unsigned) t < _person.size());
    _person[t].push_back(p);
}


bool Location::removePerson(Person* p, int t) {
    //assert((unsigned) t < _person.size());
    for (unsigned int i=0; i<_person[t].size(); i++) {
        if (_person[t][i] == p) {
            _person[t][i] = _person[t].back();
            _person[t].pop_back();
            return true;
        }
    }
    return false;
}


// Calling scope must verify that returned person is not nullptr
Person* Location::findMom() {
    vector<Person*> residents = getResidents();
    vector<Person*> potential_moms;
    int minage = 15;
    int maxage = 45;
    for (auto p: residents) {
        if (p->getSex() == FEMALE and p->getAge() >= minage and p->getAge() <= maxage) potential_moms.push_back(p);
    }
    if (potential_moms.size() == 0) return nullptr;
    int r = gsl_rng_uniform_int(RNG, potential_moms.size());
    Person* mom = potential_moms[r];
    return mom;
}


// addNeighbor - adds location p to the location's neighbor list.
// Note that this relationship is one-way.
void Location::addNeighbor(Location* p) {
    for (unsigned int i=0; i<_neighbors.size(); i++)
        if (_neighbors[i]==p) return;                                                   // already a neighbor
    _neighbors.push_back(p);
}
