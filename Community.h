// Community.h
// Manages relationship between people, locations, and mosquitoes
#ifndef __COMMUNITY_H
#define __COMMUNITY_H
#include <string>
#include <vector>
#include <map>

class Person;
class Mosquito;
class Location;

class Community {
    public:
        Community(const Parameters* parameters);
        virtual ~Community();
        bool loadPopulation(std::string szPop,std::string szImm, std::string szSwap);
        bool loadLocations(std::string szLocs,std::string szNet);
        int getNumPerson() { return _nNumPerson; }
        Person *getPerson(int n) { return _person+n; }
        int getNumInfected(int day);
        int getNumSymptomatic(int day);
        std::vector<int> getNumSusceptible();
        void populate(Person **parray, int targetpop);
        Person* getPersonByID(int id);
        bool infect(int id, Serotype serotype, int day);
        int attemptToAddMosquito(Location *p, Serotype serotype, int nInfectedByID);
        int getDay() {                                                // what day is it?
            return _nDay;
        }
        void swapImmuneStates(); 
        void updateWithdrawnStatus(); 
        void mosquitoToHumanTransmission();
        void humanToMosquitoTransmission();
        void tick(int day);                                                  // simulate one day
        void setNoSecondaryTransmission() { _bNoSecondaryTransmission = true; }
        void setMosquitoMultiplier(double f) {                        // seasonality multiplier for number of mosquitoes
            _fMosquitoCapacityMultiplier = f;
        }
        double getMosquitoMultiplier() { return _fMosquitoCapacityMultiplier; }
        void setExternalIncubation(double n) { assert(n<MAX_MOSQUITO_INCUBATION); _nExternalIncubation=n; }
        int getExternalIncubation() { return _nExternalIncubation; }
        int getNumInfectiousMosquitoes();
        int getNumExposedMosquitoes();
        void vaccinate(double f, int age=-1);
        void setVES(double f);
        void setVESs(std::vector<double> f);
        Mosquito *getInfectiousMosquito(int n);
        Mosquito *getExposedMosquito(int n);
        std::vector< std::vector<int> > getNumNewlyInfected() { return _nNumNewlyInfected; }
        std::vector< std::vector<int> > getNumNewlySymptomatic() { return _nNumNewlySymptomatic; }
        static void flagInfectedLocation(Location* _pLoc, int day) { _isHot[_pLoc][day] = true; }



    protected:
        static const Parameters* _par;
        Person *_person;                                              // the array index is equal to the ID
        std::vector< std::vector<Person*> > _personAgeCohort;         // array of pointers to people of the same age
        int _nPersonAgeCohortSizes[MAX_PERSON_AGE];                   // size of each age cohort
        double *_fMortality;                                          // mortality by year, starting from 0
        std::vector<Location*> _location;                             // the array index is equal to the ID
        std::vector< std::vector<Person*> > _exposedQueue;            // queue of people with n days of latency left
        std::vector< std::vector<Mosquito*> > _infectiousMosquitoQueue;  // queue of infectious mosquitoes with n days
                                                                         // left to live
        std::vector< std::vector<Mosquito*> > _exposedMosquitoQueue;  // queue of exposed mosquitoes with n days of latency left
        int _nDay;                                                    // current day
        int _nNumPerson;                                              // number of persons in the simulation
        int _nMaxInfectionParity;                                     // maximum number of infections (serotypes) per person
        bool _bNoSecondaryTransmission;
        double _fMosquitoCapacityMultiplier;                          // seasonality multiplier for mosquito capacity
        int _nExternalIncubation;                                     // external incubation period in days
        std::vector< std::vector<int> > _nNumNewlyInfected;
        std::vector< std::vector<int> > _nNumNewlySymptomatic;
        static std::map< Location*, std::map<int, bool> > _isHot;
        bool _uniformSwap;                                            // use original swapping (==true); or parse swap file (==false)

        void expandExposedQueues();
        void expandMosquitoQueues();
        void moveMosquito(Mosquito *m);
        void _advanceTimers();
        void _modelMosquitoMovement();
};
#endif
