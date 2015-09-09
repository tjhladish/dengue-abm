// Community.h
// Manages relationship between people, locations, and mosquitoes
#ifndef __COMMUNITY_H
#define __COMMUNITY_H
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <numeric>


class Person;
class Mosquito;
class Location;

class Community {
    public:
        Community(const Parameters* parameters);
        virtual ~Community();
        bool loadPopulation(std::string szPop,std::string szImm, std::string szSwap);
        bool loadLocations(std::string szLocs,std::string szNet);
        bool loadMosquitoes(std::string moslocFilename, std::string mosFilename);
        int getNumPerson() const { return _nNumPerson; }
        Person *getPerson(int n) const { return _person+n; }
        int getNumInfected(int day);
        int getNumSymptomatic(int day);
        std::vector<int> getNumSusceptible();
        void populate(Person **parray, int targetpop);
        Person* getPersonByID(int id);
        bool infect(int id, Serotype serotype, int day);
        void attemptToAddMosquito(Location *p, Serotype serotype, int nInfectedByID);
        int getDay() { return _nDay; }                                // what day is it?
        void swapImmuneStates(); 
        void updateDiseaseStatus();
        void mosquitoToHumanTransmission();
        void humanToMosquitoTransmission();
        void tick(int day);                                           // simulate one day
        void setNoSecondaryTransmission() { _bNoSecondaryTransmission = true; }
        void setMosquitoMultiplier(double f) { _fMosquitoCapacityMultiplier = f; }  // seasonality multiplier for number of mosquitoes
        void applyMosquitoMultiplier(double f);                    // sets multiplier and kills off infectious mosquitoes as necessary
        double getMosquitoMultiplier() const { return _fMosquitoCapacityMultiplier; }
        void setExtrinsicIncubation(int n) { _EIP = n; }
        int getExtrinsicIncubation() const { return _EIP; }
        int getNumInfectiousMosquitoes();
        int getNumExposedMosquitoes();
        void vaccinate(int time, double f, int age=-1);
        void boost(int time, double f);
        void setVES(double f);
        void setVESs(std::vector<double> f);
        Mosquito *getInfectiousMosquito(int n);
        Mosquito *getExposedMosquito(int n);
        std::vector< std::vector<int> > getNumNewlyInfected() { return _nNumNewlyInfected; }
        std::vector< std::vector<int> > getNumNewlySymptomatic() { return _nNumNewlySymptomatic; }
        std::vector< std::vector<int> > getNumVaccinatedCases() { return _nNumVaccinatedCases; }
        std::vector< std::vector<int> > getNumSevereCases() { return _nNumSevereCases; }
        static void flagInfectedLocation(Location* _pLoc, int day);

        int ageIntervalSize(int ageMin, int ageMax) { return std::accumulate(_nPersonAgeCohortSizes+ageMin, _nPersonAgeCohortSizes+ageMax,0); }

        void reset();                                                 // reset the state of the community; experimental! 
        const std::vector<Location*> getLocations() const { return _location; }
        const std::vector< std::vector<Mosquito*> > getInfectiousMosquitoes() const { return _infectiousMosquitoQueue; }
        const std::vector< std::vector<Mosquito*> > getExposedMosquitoes() const { return _exposedMosquitoQueue; }

    protected:
        static const Parameters* _par;
        Person *_person;                                              // the array index is equal to the ID
        std::vector< std::vector<Person*> > _personAgeCohort;         // array of pointers to people of the same age
        int _nPersonAgeCohortSizes[NUM_AGE_CLASSES];                  // size of each age cohort
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
        int _EIP;                                                     // extrinsic incubation period in days
        std::vector< std::vector<int> > _nNumNewlyInfected;
        std::vector< std::vector<int> > _nNumNewlySymptomatic;
        std::vector< std::vector<int> > _nNumVaccinatedCases;
        std::vector< std::vector<int> > _nNumSevereCases;
        static std::vector<std::unordered_set<Location*> > _isHot;
        bool _uniformSwap;                                            // use original swapping (==true); or parse swap file (==false)

        void expandExposedQueues();
        void expandMosquitoQueues();
        void moveMosquito(Mosquito *m);
        void mosquitoFilter(std::vector<Mosquito*>& mosquitoes, const double survival_prob);
        void _advanceTimers();
        void _modelMosquitoMovement();
};
#endif
