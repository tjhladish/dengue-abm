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
        LocationType getType() const { return _type; }
        void addPerson(Person *p, int t);
        bool removePerson(Person *p, int t);
        int getNumPerson(TimePeriod timeofday) const { return _person[(int) timeofday].size(); } 
        std::vector<Person*> getResidents() { return _person[HOME_NIGHT]; }
        Person* findMom();                                            // Try to find a resident female of reproductive age
        void setBaseMosquitoCapacity(int capacity) { _nBaseMosquitoCapacity = capacity; }
        int getBaseMosquitoCapacity() const { return _nBaseMosquitoCapacity; }
        int getCurrentInfectedMosquitoes() const { return _currentInfectedMosquitoes; }
        void setVectorControlEvent(const double efficacy, const double daily_mortality, const int start, const int duration) {
            vector_control_efficacy = efficacy;
            vector_control_daily_mortality = daily_mortality;
            vector_control_start_day = start;
            vector_control_end_day = start + duration; // efficacy is 0 on end day
        }
        double getVectorControlEfficacy() const { return vector_control_efficacy; }
        double getCurrentVectorControlEfficacy(int time) const { return (time >= vector_control_start_day and time < vector_control_end_day) ? vector_control_efficacy : 0.0; }
        double getCurrentVectorControlDailyMortality(int time) const { return (time >= vector_control_start_day and time < vector_control_end_day) ? vector_control_daily_mortality: 0.0; }
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
        LocationType _type;
        std::vector< std::vector<Person*> > _person;                  // pointers to person who come to this location
        int _nBaseMosquitoCapacity;                                   // "baseline" carrying capacity for mosquitoes
        int _currentInfectedMosquitoes;
        std::vector<Location*> _neighbors;
        static int _nNextSerial;                                      // unique ID to assign to the next Location allocated
        std::pair<double, double> _coord;                             // (x,y) coordinates for location

        double vector_control_efficacy;                                // expected percentage reduction in mosquito pop at equillibrium
        double vector_control_daily_mortality;                        // daily probability of mosquito death due to IRS at this location
        int vector_control_start_day;
        int vector_control_end_day;
};
#endif
