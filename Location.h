// Location.h
// A location, which is basically a building.
// Maintains lists of which humans are in the building at 3 points in the day
#ifndef __LOCATION_H
#define __LOCATION_H
#include "Parameters.h"

class Person;

class Location
{
    public:
        Location();
        virtual ~Location();
        int getID() { return _nID; }
        //  void setNumMosquitoes(int n) { _nNumMosquitoes = n; }
        //  int getNumMosquitoes(int n) { return _nNumMosquitoes; }
        void addPerson(Person *p, int t);
        bool removePerson(Person *p, int t);
        int getNumPerson(int timeofday) { return _nNumPerson[timeofday]; }
        int getBaseMosquitoCapacity() { return _nBaseMosquitoCapacity; }
        Person *getPerson(int id, int timeofday);
        void addNeighbor(Location *p);
        int getNumNeighbors() { return _nNumNeighbors; }
        Location *getNeighbor(int n) { return _neighbors[n]; }
        void setUndefined() { _bUndefined=true; }
        bool getUndefined() { return _bUndefined; }
        static void setDefaultMosquitoCapacity(int x) { _nDefaultMosquitoCapacity = x; }
        static int getDefaultMosquitoCapacity() { return _nDefaultMosquitoCapacity; }

        // bitmasks for coding _nPersonTimes
        static const unsigned char TIME_MORNING = 0x01;
        static const unsigned char TIME_MIDDAY = 0x02;
        static const unsigned char TIME_AFTERNOON = 0x04;
    protected:
        int _nID;                                                     // unique identifier
        //  int _nNumMosquitoes;    // number of mosquitoes living here
        Person ***_person;                                            // pointers to person who come to this location
        int _nNumPerson[3];                                           // number of people present in the morning, at midday, and in the afternoon
        int _nMaxPerson;                                              // number of people present in the morning, at midday, and in the afternoon
        int _nBaseMosquitoCapacity;                                   // "baseline" carrying capacity for mosquitoes
        Location **_neighbors;
        int _nNumNeighbors;
        int _nMaxNeighbors;
        bool _bUndefined;
        static int _nNextID;                                          // unique ID to assign to the next Location allocated
        static int _nDefaultMosquitoCapacity;
};
#endif
