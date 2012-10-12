// Community.h
// Manages relationship between people, locations, and mosquitoes
#ifndef __COMMUNITY_H
#define __COMMUNITY_H
#include <string>
#include "Parameters.h"
using namespace std;

class Person;
class Mosquito;
class Location;

class Community
{
    public:
        Community();
        virtual ~Community();
        bool loadPopulation(string szPop,string szImm);
        bool loadLocations(string szLocs,string szNet);
        int getNumPerson() { return _nNumPerson; }
        Person *getPerson(int n) { return _person+n; }
        int getNumInfected(int day);
        int getNumSymptomatic(int day);
        int getNumSusceptible(Serotype serotype);
        void populate(gsl_rng *rng, Person **parray, int targetpop);
        int getNumLocation() { return _nNumLocation; }
        bool infect(gsl_rng *rng, int id, Serotype serotype, int day);
        int addMosquito(gsl_rng *rng, Location *p, Serotype serotype, int nInfectedByID);
        int getDay() {                                                // what day is it?
            return _nDay;
        }
        void tick(gsl_rng *rng);                                      // simulate one day
        void setBetaPM(double f) { _fBetaPM = f; }
        void setBetaMP(double f) { _fBetaMP = f; }
        void setMosquitoMoveProbability(double f) { _fMosquitoMoveProb = f; }
        void setMosquitoTeleportProbability(double f) { _fMosquitoTeleportProb = f; }
        double getMosquitoMoveProbability() { return _fMosquitoMoveProb; }
        double setMosquitoTeleportProbability() { return _fMosquitoTeleportProb; }
        void setNoSecondaryTransmission() { _bNoSecondaryTransmission = true; }
        void setMosquitoMultiplier(double f) {                        // seasonality multiplier for number of mosquitoes
            _fMosquitoCapacityMultiplier = f;
        }
        double getMosquitoMultiplier() { return _fMosquitoCapacityMultiplier; }
        int getNumInfectiousMosquitoes();
        int getNumExposedMosquitoes();
        void vaccinate(gsl_rng *rng, double f, int age=-1);
        void setVES(double f);
        void setVESs(double f1,double f2,double f3,double f4);
        void setVEI(double f);
        void setVEP(double f);
        void setPrimaryPathogenicityScaling(double f1,double f2,double f3,double f4) {
            _fPrimarySymptomaticScaling[0] = f1;
            _fPrimarySymptomaticScaling[1] = f2;
            _fPrimarySymptomaticScaling[2] = f3;
            _fPrimarySymptomaticScaling[3] = f4;
        }
        void setSecondaryPathogenicityScaling(double f1,double f2,double f3,double f4) {
            _fSecondarySymptomaticScaling[0] = f1;
            _fSecondarySymptomaticScaling[1] = f2;
            _fSecondarySymptomaticScaling[2] = f3;
            _fSecondarySymptomaticScaling[3] = f4;
        }
        double getSecondaryPathogenicityScaling(Serotype serotype) { return _fSecondarySymptomaticScaling[(int) serotype]; }
        int getMaxInfectionParity() { return _nMaxInfectionParity; }
        void setMaxInfectionParity(int n) { _nMaxInfectionParity=n; }
        Mosquito *getInfectiousMosquito(int n);
        Mosquito *getExposedMosquito(int n);
        const int *getNumNewlyInfected(Serotype serotype) { return _nNumNewlyInfected[(int) serotype]; }
        const int *getNumNewlySymptomatic(Serotype serotype) { return _nNumNewlySymptomatic[(int) serotype]; }

        //const static int MAXSEROTYPES = 4;                            // maximum number of serotypes
        const static int STEPSPERDAY = 3;                             // number of time steps per day
        const static int MAXINCUBATION = 15;                          // maximum incubation period for humans
        const static int MOSQUITOINCUBATION = 11;                     // number of days for mosquito incubation (extrinsic incubation period)
        const static int MAXRUNTIME = 7400;                           // maximum number of simulation days (+ extra for mosquito lifetime)
        const static double SYMPTOMATICBYAGE[Person::MAXPERSONAGE];   // for some serotypes, the fraction who are symptomatic upon primary infection

    protected:
        string _szPopulationFilename;                                 // population data filename
        string _szImmunityFilename;                                   // immune status data filename
        string _szLocationFilename;                                   // location filename
        string _szNetworkFilename;                                    // location network filename
        //  Mosquito *_mosquito;
        Person *_person;                                              // the array index is equal to the ID
        Person ***_personAgeCohort;                                   // array of pointers to people of the same age
        int _nPersonAgeCohortSizes[Person::MAXPERSONAGE];             // size of each age cohort
        int _nPersonAgeCohortMaxSize;                                 // size of largest cohort
        int _nOriginalPersonAgeCohortSizes[Person::MAXPERSONAGE];     // size of each age cohort at simulation start time
        double *_fMortality;                                          // mortality by year, starting from 0
        Location *_location;                                          // the array index is equal to the ID
        int **_numLocationMosquitoCreated;                            // number of instantiated mosquitoes at this location at time t.  the first 
                                                                      // array index is equal to the location ID, the second to the day
        Person ***_exposedQueue;                                      // queue of people with n days of latency left
        int _nExposedQueueCapacity;                                   // capacity for each day
        Mosquito ***_infectiousMosquitoQueue;                         // queue of infectious mosquitoes with n days left to live
        Mosquito ***_exposedMosquitoQueue;                            // queue of exposed mosquitoes with n days of latency left
        int _nMosquitoQueueCapacity;                                  // capacity for each day
        int _nDay;                                                    // current day
        int _nNumPerson;                                              // number of persons in the simulation
        int _nMaxPerson;                                              // max of persons ever in the simulation
        int _nNumLocation;                                            // number of locations in the simulation
        int _nMaxInfectionParity;                                     // maximum number of infections (serotypes) per person
        double _fBetaPM;                                              // scales person-to-mosquito transmission
        double _fBetaMP;                                              // scales mosquito-to-person transmission (includes bite rate)
        double _fMosquitoMoveProb;                                    // daily probability of mosquito migration
        double _fMosquitoTeleportProb;                                // daily probability of mosquito teleportation (long-range movement)
        bool _bNoSecondaryTransmission;
        double _fMosquitoCapacityMultiplier;                          // seasonality multiplier for mosquito capacity
        double _fPrimarySymptomaticScaling[NUM_OF_SEROTYPES];
        double _fSecondarySymptomaticScaling[NUM_OF_SEROTYPES];
        int _nNumNewlyInfected[NUM_OF_SEROTYPES][MAXRUNTIME];
        int _nNumNewlySymptomatic[NUM_OF_SEROTYPES][MAXRUNTIME];

        void expandExposedQueues();
        void expandMosquitoQueues();
        void moveMosquito(Mosquito *m, gsl_rng *rng);
        double _fDailyBitingPDF[STEPSPERDAY];                         // probability of biting at 3 different times of day (as defined in Location.h)
};
#endif
