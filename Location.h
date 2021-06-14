// Location.h
// A location, which is basically a building.
// Maintains lists of which humans are in the building at 3 points in the day
#ifndef __LOCATION_H
#define __LOCATION_H

#include <queue>

class Person;

struct InsecticideTreatmentEvent {
    InsecticideTreatmentEvent(double eff, double m, int s, int d) : efficacy(eff), daily_mortality(m), start_day(s), end_day(s+d) {};
    double efficacy;
    double daily_mortality;
    int start_day;
    int end_day; // efficacy is 0 on end day
};


inline bool operator<(const InsecticideTreatmentEvent& lhs, const InsecticideTreatmentEvent& rhs) {
    return lhs.start_day > rhs.start_day;
}


class Location {
    public:
        Location();
        virtual ~Location();
        void setID(int id) { _ID = id; }
        int getID() const { return _ID; }
        int getSerial() const { return _serial; }
        void setType(LocationType t) { _type = t; }
        LocationType getType() const { return _type; }
        void setTrialArm(int trial_arm) { _trial_arm = trial_arm; }
        int getTrialArm() const { return _trial_arm; }
        void setSurveilled(bool surveilled) { _surveilled = surveilled; }
        bool isSurveilled() const { return _surveilled; }

        void addPerson(Person *p, int t);
        bool removePerson(Person *p, int t);
        int getNumPerson(TimePeriod timeofday) const { return _person[(int) timeofday].size(); }
        std::vector<Person*> getResidents() { return _person[HOME_NIGHT]; }
        Person* findMom();                                            // Try to find a resident female of reproductive age
        void setBaseMosquitoCapacity(int capacity) { _nBaseMosquitoCapacity = capacity; }
                                                                                      // not really killing of I and S mosquitoes in the same way . . .
        int getBaseMosquitoCapacity() const { return _nBaseMosquitoCapacity; }
        int getCurrentInfectedMosquitoes() const { return _currentInfectedMosquitoes; }
        void scheduleVectorControlEvent(const double efficacy, const double daily_mortality, const int start, const int duration) {
            ITQ.emplace(efficacy, daily_mortality, start, duration);
        }
        const InsecticideTreatmentEvent* getCurrentVectorControl() const { return &ITQ.top(); }
        void updateVectorControlQueue(int now) {
            while (vectorControlScheduled() and (now >= ITQ.top().end_day)) ITQ.pop(); }
        inline bool vectorControlScheduled() const { return (not ITQ.empty()); }
        // vectorControlActive() assumes updateVectorControlQueue() has been called recently enough that the top element is not out-of-date
        bool vectorControlActive(int now) const {
            assert(ITQ.empty() or now < ITQ.top().end_day);
            return (vectorControlScheduled() and (now >= ITQ.top().start_day) and (now < ITQ.top().end_day));
        }
        double getCurrentVectorControlEfficacy(int now) const { return vectorControlActive(now) ? ITQ.top().efficacy : 0.0; }
        double getCurrentVectorControlDailyMortality(int now) const { return vectorControlActive(now) ? ITQ.top().daily_mortality : 0.0; }
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
        int _trial_arm;
        bool _surveilled;
        std::vector< std::vector<Person*> > _person;                  // pointers to person who come to this location
        int _nBaseMosquitoCapacity;                                   // "baseline" carrying capacity for mosquitoes
        int _currentInfectedMosquitoes;
        std::vector<Location*> _neighbors;
        static int _nNextSerial;                                      // unique ID to assign to the next Location allocated
        std::pair<double, double> _coord;                             // (x,y) coordinates for location

        std::priority_queue<InsecticideTreatmentEvent> ITQ;           // insecticide treatment event priority queue
};
#endif
