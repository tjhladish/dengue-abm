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
        void setID(int id) { _ID = id; }
        int getID() { return _ID; }
        int getSerial() { return _serial; }
        void addPerson(Person *p, int t);
        bool removePerson(Person *p, int t);
        int getNumPerson(int timeofday) { return _person[timeofday].size(); } 
        void setBaseMosquitoCapacity(int capacity) { _nBaseMosquitoCapacity = capacity; }
        int getBaseMosquitoCapacity() { return _nBaseMosquitoCapacity; }
        int getCurrentInfectedMosquitoes() { return _currentInfectedMosquitoes; }
        void addInfectedMosquito() { _currentInfectedMosquitoes++; }
        void removeInfectedMosquito() { _currentInfectedMosquitoes--; }
        void addNeighbor(Location *p);
        int getNumNeighbors() { return _neighbors.size(); }
        Location *getNeighbor(int n) { return _neighbors[n]; }
        void setUndefined() { _bUndefined=true; }
        bool getUndefined() { return _bUndefined; }
        inline Person* getPerson(int idx, int timeofday) { return _person[timeofday][idx]; }
        void setCoordinates(std::pair<double, double> c) { _coord = c; }
        std::pair<double, double> getCoordinates() { return _coord; }
        void setX(double x) { _coord.first = x; }
        void setY(double y) { _coord.second = y; }
        double getX() { return _coord.first; }
        double getY() { return _coord.second; }

    protected:
        int _ID;                                                     // original identifier in location file
        int _serial;                                                  // unique identifier assigned on construction
        std::vector< std::vector<Person*> > _person;                  // pointers to person who come to this location
        int _nBaseMosquitoCapacity;                                   // "baseline" carrying capacity for mosquitoes
        int _currentInfectedMosquitoes;
        std::vector<Location*> _neighbors;
        bool _bUndefined;
        static int _nNextSerial;                                          // unique ID to assign to the next Location allocated
        std::pair<double, double> _coord;                                  // (x,y) coordinates for location
};
#endif
