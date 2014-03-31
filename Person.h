// Person.h
// A single individual.

#ifndef __PERSON_H
#define __PERSON_H
#include <bitset>
#include <vector>
#include <climits>
#include "Parameters.h"

class Location;

class Infection {
    friend class Person;
    Infection() {
        infectedByID = INT_MIN;
        infectedPlace = INT_MIN;
        infectedTime = INT_MIN;
        infectiousTime = INT_MIN;
        symptomTime = INT_MIN;
        recoveryTime = INT_MIN;

        withdrawnTime = INT_MAX;
        serotype = NULL_SEROTYPE;
    };

    Infection* operator=(const Infection* o) {
        if (this != o) {
            infectedByID   = o->infectedByID;
            infectedPlace  = o->infectedPlace;
            infectedTime   = o->infectedTime;
            infectiousTime = o->infectiousTime;
            symptomTime    = o->symptomTime;
            recoveryTime   = o->recoveryTime;
            withdrawnTime  = o->withdrawnTime;
            serotype       = o->serotype;
        }
        return this;
    }

    int infectedByID;                               // who infected this person
    int infectedPlace;                              // where infected?
    int infectedTime;                               // when infected?
    int infectiousTime;                             // when infectious period starts
    int symptomTime;                                // when symptoms start
    int recoveryTime;                               // when recovered?
    int withdrawnTime;                              // when person withdraws to home
    Serotype serotype; 
};

class Person
{
    public:
        Person();
        ~Person();
        inline int getID() { return _nID; }
        int getAge() { return _nAge; }
        void setAge(int n) { _nAge = n; }
        int getLifespan() { return _nLifespan; }
        void setLifespan(int n) { _nLifespan = n; }
        int getHomeID() { return _nHomeID; }
        void setHomeID(int n) { _nHomeID = n; }
        int getWorkID() { return _nWorkID; }
        void setWorkID(int n) { _nWorkID = n; }
        void setImmunity(Serotype serotype) { _nImmunity[(int) serotype] = 1; }
        void copyImmunity(const Person *p);
        void resetImmunity();
        void appendToSwapProbabilities(std::pair<int, double> p) { _swap_probabilities.push_back(p); }
        std::vector<std::pair<int, double> > getSwapProbabilities() { return _swap_probabilities; }

        bool isSusceptible(Serotype serotype);                        // is susceptible to serotype (and is alive)
        int getInfectionParity();                                     // is immune to n serotypes
        inline Location *getLocation(int timeofday) { return _pLocation[timeofday]; }
        inline void setLocation(Location *p, int timeofday) { _pLocation[timeofday] = p; }

        inline int getInfectedByID(int infectionsago=0)   { return infectionHistory[getNumInfections() - 1 - infectionsago]->infectedByID; }
        inline int getInfectedPlace(int infectionsago=0)  { return infectionHistory[getNumInfections() - 1 - infectionsago]->infectedPlace; }
        inline int getInfectedTime(int infectionsago=0)   { return infectionHistory[getNumInfections() - 1 - infectionsago]->infectedTime; }
        inline int getInfectiousTime(int infectionsago=0) { return infectionHistory[getNumInfections() - 1 - infectionsago]->infectiousTime; }
        inline int getSymptomTime(int infectionsago=0)    { return infectionHistory[getNumInfections() - 1 - infectionsago]->symptomTime; }
        inline int getRecoveryTime(int infectionsago=0)   { return infectionHistory[getNumInfections() - 1 - infectionsago]->recoveryTime; }
        inline int getWithdrawnTime(int infectionsago=0)  { return infectionHistory[getNumInfections() - 1 - infectionsago]->withdrawnTime; }
        inline Serotype getSerotype(int infectionsago=0)  { return infectionHistory[getNumInfections() - 1 - infectionsago]->serotype; }

        inline void setRecoveryTime(int time, int infectionsago=0) { infectionHistory[getNumInfections() - 1 - infectionsago]->recoveryTime = time; }
        bool isWithdrawn(int time);                                   // at home sick?
        inline int getNumInfections() const { return infectionHistory.size(); }

        bool infect(int sourceid, Serotype serotype, int time, int sourceloc);
        bool isViremic(int time);

        void kill(int time);
        bool isDead() { return _bDead; }
        bool naturalDeath(int t);                                     // die of old age check?

        bool isNewlyInfected(int time);                               // became infected today?
        bool isInfected(int time);                                    // is currently infected
        bool isSymptomatic(int time);                                 // has symptoms
        bool isVaccinated() {                                         // has been vaccinated
            return _bVaccinated;
        }
        bool isInfectable(Serotype serotype, int time);               // more complicated than isSusceptible

        bool vaccinate();                                 // vaccinate this person
        static void setPar(const Parameters* par) { _par = par; }

        static const double _fIncubationDistribution[MAX_INCUBATION];

        Infection& initializeNewInfection();

    protected:
        int _nID;                                                     // unique identifier
        int _nHomeID;                                                 // family membership
        int _nWorkID;                                                 // ID of location of work
        bool _bCase;                                                  // ever detected as case?
        Location *_pLocation[STEPS_PER_DAY];                            // where this person is at morning, day, and evening
        int _nAge;                                                    // age in years
        int _nLifespan;                                               // lifespan in years
        bool _bDead;                                                  // is dead
        std::bitset<NUM_OF_SEROTYPES> _nImmunity;                     // bitmask of serotype exposure
        bool _bVaccinated;                                            // has been vaccinated
        std::vector<std::pair<int,double> > _swap_probabilities;      // list of the nearest people one year younger, with distances
        std::vector<Infection*> infectionHistory;
        void clearInfectionHistory();

        static const Parameters* _par;
        static int _nNextID;                                          // unique ID to assign to the next Person allocated
};
#endif
