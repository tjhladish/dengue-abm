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
        int getID() const { return _ID; }
        int getSerial() const { return _serial; }
        void addPerson(Person *p, int t);
        bool removePerson(Person *p, int t);
        int getNumPerson(TimePeriod timeofday) const { return _person[(int) timeofday].size(); } 
        std::vector<Person*> getResidents() { return _person[HOME_NIGHT]; }
        Person* findMom();                                            // Try to find a resident female of reproductive age
        void setBaseMosquitoCapacity(int capacity) { _nBaseMosquitoCapacity = capacity; }
        int getBaseMosquitoCapacity() const { return _nBaseMosquitoCapacity; }
        int getCurrentInfectedMosquitoes() { return _currentInfectedMosquitoes; }
        void addInfectedMosquito() { _currentInfectedMosquitoes++; }
        void addInfectedMosquitoes(int n) { _currentInfectedMosquitoes += n; }
        void removeInfectedMosquito() { _currentInfectedMosquitoes--; }
        void removeInfectedMosquitoes(int n) { _currentInfectedMosquitoes -= n; }
        void clearInfectedMosquitoes() { _currentInfectedMosquitoes = 0; }
        void addNeighbor(Location *p);
        int getNumNeighbors() const { return _neighbors.size(); }
        Location *getNeighbor(int n) { return _neighbors[n]; }
        inline Person* getPerson(int idx, TimePeriod timeofday) { return _person[(int) timeofday][idx]; }
        void setCoordinates(std::pair<double, double> c) { _coord = c; }
        std::pair<double, double> getCoordinates() { return _coord; }
        void setX(double x) { _coord.first = x; }
        void setY(double y) { _coord.second = y; }
        double getX() const { return _coord.first; }
        double getY() const { return _coord.second; }

        bool operator == ( const Location* other ) const { return ( ( _ID == other->_ID ) && ( _serial == other->_serial ) ); }

    protected:
        int _ID;                                                      // original identifier in location file
        int _serial;                                                  // unique identifier assigned on construction
        std::vector< std::vector<Person*> > _person;                  // pointers to person who come to this location
        int _nBaseMosquitoCapacity;                                   // "baseline" carrying capacity for mosquitoes
        int _currentInfectedMosquitoes;
        std::vector<Location*> _neighbors;
        static int _nNextSerial;                                      // unique ID to assign to the next Location allocated
        std::pair<double, double> _coord;                             // (x,y) coordinates for location
};
#endif
