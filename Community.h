// Community.h
// Manages relationship between people, locations, and mosquitoes
#ifndef __COMMUNITY_H
#define __COMMUNITY_H
#include <string>
#include <vector>
#include <map>
#include "types.h"

class Person;
class Mosquito;
class Location;

class Community
{
    public:
        Community(const Parameters* parameters);
        virtual ~Community();
        bool loadPopulation(std::string szPop,std::string szImm);
        bool loadLocations(std::string szLocs,std::string szNet);
        int getNumPerson() { return _nNumPerson; }
        Person *getPerson(int n) { return _person+n; }
        int getNumInfected(int day);
        int getNumSymptomatic(int day);
        std::vector<int> getNumSusceptible();
        void populate(Person **parray, int targetpop);
        bool infect(int id, Serotype serotype, int day);
        int addMosquito(Location *p, Serotype serotype, int nInfectedByID);
        int getDay() {                                                // what day is it?
            return _nDay;
        }
        void swapImmuneStates(); 
        void updateWithdrawnStatus(); 
        void mosquitoToHumanTransmission();
        void humanToMosquitoTransmission();
        void tick();                                                  // simulate one day
        void setNoSecondaryTransmission() { _bNoSecondaryTransmission = true; }
        void setMosquitoMultiplier(double f) {                        // seasonality multiplier for number of mosquitoes
            _fMosquitoCapacityMultiplier = f;
        }
        double getMosquitoMultiplier() { return _fMosquitoCapacityMultiplier; }
        int getNumInfectiousMosquitoes();
        int getNumExposedMosquitoes();
        void vaccinate(double f, int age=-1);
        void setVES(double f);
        void setVESs(std::vector<double> f);
        Mosquito *getInfectiousMosquito(int n);
        Mosquito *getExposedMosquito(int n);
        VectorMP< VectorMP<int> > getNumNewlyInfected() { return _nNumNewlyInfected; }
        VectorMP< VectorMP<int> > getNumNewlySymptomatic() { return _nNumNewlySymptomatic; }
        static void flagInfectedLocation(Location* _pLoc, int day) { _isHot[_pLoc][day] = true; }



    protected:
        static const Parameters* _par;
        std::string _szPopulationFilename;                            // population data filename
        std::string _szImmunityFilename;                              // immune status data filename
        std::string _szLocationFilename;                              // location filename
        std::string _szNetworkFilename;                               // location network filename
        Person *_person;                                              // the array index is equal to the ID
        std::vector< std::vector<Person*> > _personAgeCohort;         // array of pointers to people of the same age
        int _nPersonAgeCohortSizes[MAX_PERSON_AGE];                   // size of each age cohort
        double *_fMortality;                                          // mortality by year, starting from 0
        std::vector<Location*> _location;                             // the array index is equal to the ID
        std::vector< std::vector<int> > _numLocationMosquitoCreated;  // number of instantiated mosquitoes at 
                                                                      // this location at time t.  the first 
                                                                      // array index is equal to the location ID,
                                                                      // the second to the day
        VectorMP< VectorMP<Person*> > _exposedQueue;            // queue of people with n days of latency left
        //std::vector< std::vector<Mosquito*> > _infectiousMosquitoQueue;  // queue of infectious mosquitoes with n days
        //                                                                 // left to live
        //std::vector< std::vector<Mosquito*> > _exposedMosquitoQueue;  // queue of exposed mosquitoes with n days of latency left
        VectorMP< VectorMP<Mosquito*> > _infectiousMosquitoQueue;  // queue of infectious mosquitoes with n days
                                                                         // left to live
        VectorMP< VectorMP<Mosquito*> > _exposedMosquitoQueue;  // queue of exposed mosquitoes with n days of latency left
        int _nDay;                                                    // current day
        int _nNumPerson;                                              // number of persons in the simulation
        int _nMaxInfectionParity;                                     // maximum number of infections (serotypes) per person
        bool _bNoSecondaryTransmission;
        double _fMosquitoCapacityMultiplier;                          // seasonality multiplier for mosquito capacity
        VectorMP< VectorMP<int> > _nNumNewlyInfected;
        VectorMP< VectorMP<int> > _nNumNewlySymptomatic;
        //std::vector< std::vector<int> > _nNumNewlyInfected;
        //std::vector< std::vector<int> > _nNumNewlySymptomatic;
        static MapMP< Location*, MapMP<int, bool> > _isHot;

        //void expandExposedQueues();
        //void expandMosquitoQueues();
        void moveMosquito(Mosquito *m);
        void _advanceTimers();
        void _modelMosquitoMovement();
};
#endif
