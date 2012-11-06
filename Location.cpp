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

int Location::_nNextID = 0;
//int Location::_nDefaultMosquitoCapacity;

Location::Location()
    : _person(STEPSPERDAY, vector<Person*>(0) ) {
    _nID = _nNextID++;
    _nBaseMosquitoCapacity = 0;
    _bUndefined = false;
}


Location::~Location() {
    _person.clear();
    _neighbors.clear();
    /*if (_person.size()) {
        for (unsigned int i = 0; i<_person.size(); i++) {
            for (unsigned int j = 0; j<_person[i].size(); j++) {
                delete _person[i][j]; 
            }
        }
    }
    if (_neighbors.size()) {
        for (unsigned int i = 0; i<_neighbors.size(); i++) {
            delete _neighbors[i];
        }
    }*/
}


void Location::addPerson(Person *p, int t) {
    //assert((unsigned) t < _person.size());
    _person[t].push_back(p);
}


bool Location::removePerson(Person *p, int t) {
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




// addNeighbor - adds location p to the location's neighbor list.
// Note that this relationship is one-way.
void Location::addNeighbor(Location *p) {
    for (unsigned int i=0; i<_neighbors.size(); i++)
        if (_neighbors[i]==p) return;                                                   // already a neighbor
    _neighbors.push_back(p);
}
