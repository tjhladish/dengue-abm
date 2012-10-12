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

using namespace std;

int Location::_nNextID = 0;
int Location::_nDefaultMosquitoCapacity = 20;

Location::Location() {
    _nID = _nNextID++;
    _person = NULL;
    _nMaxPerson = 0;
    _nNumPerson[0] = _nNumPerson[1] = _nNumPerson[2] = 0;
    _neighbors = NULL;
    _nNumNeighbors = _nMaxNeighbors = 0;
    _nBaseMosquitoCapacity = _nDefaultMosquitoCapacity;
    _bUndefined = false;
}


Location::~Location() {
    if (_person) {
        delete [] _person[0];
        delete [] _person[1];
        delete [] _person[2];
        delete [] _person;
    }
    if (_neighbors)
        delete [] _neighbors;
}


void Location::addPerson(Person *p, int t) {
    if (!_person) {
        _person = new Person **[3];
        _nMaxPerson = 20;
        _person[0] = new Person*[_nMaxPerson];
        _person[1] = new Person*[_nMaxPerson];
        _person[2] = new Person*[_nMaxPerson];
    }
    if (_nNumPerson[t]>=_nMaxPerson) {
        _nMaxPerson *= 2;
        for (int t=0; t<3; t++) {
            Person **temp = new Person*[_nMaxPerson];
            for (int i=0; i<_nNumPerson[t]; i++) {
                temp[i] = _person[t][i];
            }
            delete [] _person[t];
            _person[t] = temp;
        }
    }
    _person[t][_nNumPerson[t]] = p;
    _nNumPerson[t]++;
}


bool Location::removePerson(Person *p, int t) {
    assert(_person);
    for (int i=0; i<_nNumPerson[t]; i++) {
        if (_person[t][i] == p) {
            _person[t][i] = _person[t][_nNumPerson[t]-1];
            _nNumPerson[t]--;
            return true;
        }
    }
    return false;
}


Person *Location::getPerson(int id, int timeofday) {
    assert(id<_nNumPerson[timeofday]);
    assert(timeofday<3);
    return _person[timeofday][id];
}


// addNeighbor - adds location p to the location's neighbor list.
// Note that this relationship is one-way.
void Location::addNeighbor(Location *p) {
    if (_neighbors==NULL) {
        _nMaxNeighbors = 20;
        _neighbors = new Location *[_nMaxNeighbors];
    }
    for (int i=0; i<_nNumNeighbors; i++)
        if (_neighbors[i]==p)
            return;                                                   // already a neighbor
    if (_nNumNeighbors>=_nMaxNeighbors) {
        Location **temp = new Location *[_nMaxNeighbors*2];
        for (int i=0; i<_nNumNeighbors; i++)
            temp[i] = _neighbors[i];
        delete [] _neighbors;
        _neighbors = temp;
    }
    _neighbors[_nNumNeighbors++] = p;
}
