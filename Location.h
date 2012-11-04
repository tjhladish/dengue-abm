// Location.h
// A location, which is basically a building.
// Maintains lists of which humans are in the building at 3 points in the day
#ifndef __LOCATION_H
#define __LOCATION_H

class Person;

class Location {
    public:
        Location();
        virtual ~Location();
        int getID() { return _nID; }
        void addPerson(Person *p, int t);
        bool removePerson(Person *p, int t);
        int getNumPerson(int timeofday) { return _person[timeofday].size(); } 
        void setBaseMosquitoCapacity(int capacity) { _nBaseMosquitoCapacity = capacity; }
        int getBaseMosquitoCapacity() { return _nBaseMosquitoCapacity; }
        Person *getPerson(int id, int timeofday);
        void addNeighbor(Location *p);
        int getNumNeighbors() { return _neighbors.size(); }
        Location *getNeighbor(int n) { return _neighbors[n]; }
        void setUndefined() { _bUndefined=true; }
        bool getUndefined() { return _bUndefined; }
        //static void setDefaultMosquitoCapacity(int x) { _nDefaultMosquitoCapacity = x; }

    protected:
        int _nID;                                                     // unique identifier
        std::vector< std::vector<Person*> > _person;                  // pointers to person who come to this location
        int _nBaseMosquitoCapacity;                                   // "baseline" carrying capacity for mosquitoes
        std::vector<Location*> _neighbors;
        bool _bUndefined;
        static int _nNextID;                                          // unique ID to assign to the next Location allocated
        //static int _nDefaultMosquitoCapacity;
};
#endif
