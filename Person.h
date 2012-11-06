// Person.h
// A single individual.

#ifndef __PERSON_H
#define __PERSON_H
#include <bitset>
#include "Parameters.h"

class Location;

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

        bool isSusceptible(Serotype serotype);                        // is susceptible to serotype (and is alive)
        int getInfectionParity();                                     // is immune to n serotypes
        Serotype getSerotype(int infectionsago=0) {                   // last infected with
            return _eSerotype[infectionsago];
        }
        Location *getLocation(int timeofday) { return _pLocation[timeofday]; }
        void setLocation(Location *p, int timeofday) { _pLocation[timeofday] = p; }
        int getInfectedTime(int infectionsago=0) { return _nInfectedTime[infectionsago]; }
        int getInfectedPlace(int infectionsago=0) { return _nInfectedPlace[infectionsago]; }
        int getInfectiousTime(int infectionsago=0) { return _nInfectiousTime[infectionsago]; }
        int getSymptomTime(int infectionsago=0) { return _nSymptomTime[infectionsago]; }
        int getRecoveryTime(int infectionsago=0) { return _nRecoveryTime[infectionsago]; }
        bool isWithdrawn(int time);                                   // at home sick?
        int getWithdrawnTime(int infectionsago=0) {                   // first day at home sick?
            return _nWithdrawnTime[infectionsago];
        }
        int getInfectedByID(int infectionsago=0) { return _nInfectedByID[infectionsago]; }
        int getNumInfections() { return _nNumInfections; }

        bool infect(int sourceid, Serotype serotype, int time, int sourceloc);
        bool isViremic(int time);

        void kill(int time);
        bool isDead() { return _bDead; }
        bool naturalDeath(int t);                                     // die of old age check?

        bool isInfected(int time);                                    // is currently infected
        bool isSymptomatic(int time);                                 // has symptoms
        bool isVaccinated() {                                         // has been vaccinated
            return _bVaccinated;
        }
        bool vaccinate();                                 // vaccinate this person
        static void setPar(const Parameters* par) { _par = par; }

        static const double _fIncubationDistribution[MAXINCUBATION];

    protected:
        void pushInfectionHistory();

    protected:
        int _nID;                                                     // unique identifier
        int _nHomeID;                                                 // family membership
        int _nWorkID;                                                 // ID of location of work
        bool _bCase;                                                  // ever detected as case?
        Location *_pLocation[STEPSPERDAY];                            // where this person is at morning, day, and evening
        int _nAge;                                                    // age in years
        int _nLifespan;                                               // lifespan in years
        bool _bDead;                                                  // is dead
        std::bitset<NUM_OF_SEROTYPES> _nImmunity;                     // bitmask of serotype exposure
        bool _bVaccinated;                                            // has been vaccinated

        // _nNumInfections counts the number of infections, up to MAXHISTORY infections.
        // arrays that use this index store the most recent infection at index 0 and data from n years ago at index n
        int _nNumInfections;                                          // number of dengue infections, index of the first exposure
        int _nInfectedByID[MAXHISTORY];                               // who infected this person
        int _nInfectedPlace[MAXHISTORY];                              // where infected?
        int _nInfectedTime[MAXHISTORY];                               // when infected?
        int _nInfectiousTime[MAXHISTORY];                             // when infectious period starts
        int _nSymptomTime[MAXHISTORY];                                // when symptoms start
        int _nWithdrawnTime[MAXHISTORY];                              // when person withdraws to home
        int _nRecoveryTime[MAXHISTORY];                               // when recovered?
        Serotype _eSerotype[MAXHISTORY]; 

        static const Parameters* _par;
        static int _nNextID;                                          // unique ID to assign to the next Person allocated
};
#endif
