// Person.h
// A single individual.
// Note that serotypes are numbered 1-4

#ifndef __PERSON_H
#define __PERSON_H
#include "Parameters.h"

class Location;

class Person
{
    public:
        Person();
        ~Person();
        int getID() { return _nID; }
        int getAge() { return _nAge; }
        void setAge(int n) { _nAge = n; }
        int getLifespan() { return _nLifespan; }
        void setLifespan(int n) { _nLifespan = n; }
        int getHomeID() { return _nHomeID; }
        void setHomeID(int n) { _nHomeID = n; }
        int getWorkID() { return _nWorkID; }
        void setWorkID(int n) { _nWorkID = n; }
        void setImmunity(Serotype serotype) { _nImmunity |= 1<<((int) serotype); }
        void copyImmunity(const Person *p);
        void resetImmunity();

        bool isSusceptible(Serotype serotype);                       // is susceptible to serotype (and is alive)
        int getInfectionParity();                                     // is immune to n serotypes
        Serotype getSerotype(int infectionsago=0) {                        // last infected with
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

        bool infect(gsl_rng *rng, int sourceid, Serotype serotype, int time, int sourceloc, double primarysymptomatic, double secondaryscaling, int maxinfectionparity);
        bool isViremic(int time);

        void kill(int time);
        bool isDead() { return _bDead; }
        bool naturalDeath(int t);                                     // die of old age check?

        static void generateRandomIDs(gsl_rng *rng, int num, int bound, int *ids);

        //  bool isCase();                // is detected as a case
        bool isInfected(int time);                                    // is currently infected
        bool isSymptomatic(int time);                                 // has symptoms
        bool isVaccinated() {                                         // has been vaccinated
            return _bVaccinated;
        }
        //  double getSusceptibility();   // 0=not susceptible, 1=max susceptibility
        //  double getInfectiousness();   // 0=not infectious
        bool vaccinate(gsl_rng *rng);                                 // vaccinate this person
        static void setVES(double f) { _fVES[0] = _fVES[1] = _fVES[2] = _fVES[3] = f; }
        static void setVESs(double f1,double f2,double f3,double f4) { _fVES[0] = f1; _fVES[1] = f2; _fVES[2] = f3; _fVES[3] = f4; }
        static void setVEI(double f) { _fVEI = f; }
        static void setVEP(double f) { _fVEP = f; }
        static double *getVES() { return _fVES; }
        static double getVEI() { return _fVEI; }
        static double getVEP() { return _fVEP; }
        static void setDaysImmune(int n) { _nDaysImmune = n; }
        static int getDaysImmune() { return _nDaysImmune; }

        const static int MAXPERSONAGE = 95;                           // maximum age-1 for a person
        static const int MAXINCUBATION = 9;                           // max incubation period in days
        static const double _fIncubationDistribution[MAXINCUBATION];
        static const int MAXHISTORY = 50;                             // length of exposure history in years

    protected:
        void pushInfectionHistory();

    protected:
        int _nID;                                                     // unique identifier
        int _nHomeID;                                                 // family membership
        int _nWorkID;                                                 // ID of location of work
        bool _bCase;                                                  // ever detected as case?
        Location *_pLocation[3];                                      // where this person is at morning, day, and evening
        int _nAge;                                                    // age in years
        int _nLifespan;                                               // lifespan in years
        bool _bDead;                                                  // is dead
        int _nImmunity;                                               // bitmask of serotype exposure
        bool _bVaccinated;                                            // has been vaccinated

        // __nNumInfections counts the number of infections, up to MAXHISTORY infections.
        // arrays that use this index store the most recent infection at index 0 and data from n years ago at index n
        int _nNumInfections;                                          // number of dengue infections, index of the first exposure
        int _nInfectedByID[MAXHISTORY];                               // who infected this person
        int _nInfectedPlace[MAXHISTORY];                              // where infected?
        int _nInfectedTime[MAXHISTORY];                               // when infected?
        int _nInfectiousTime[MAXHISTORY];                             // when infectious period starts
        int _nSymptomTime[MAXHISTORY];                                // when symptoms start
        int _nWithdrawnTime[MAXHISTORY];                              // when person withdraws to home
        int _nRecoveryTime[MAXHISTORY];                               // when recovered?
        Serotype _eSerotype[MAXHISTORY];                              // from 1-4

        static int _nDaysImmune;                                      // length of complete cross-protective immunity in days
        static double _fVES[4];                                       // vaccine protection (all-or-none, different for each serotype)
        static double _fVEI;                                          // vaccine protection
        static double _fVEP;                                          // vaccine protection
        static int _nNextID;                                          // unique ID to assign to the next Person allocated
};
#endif
