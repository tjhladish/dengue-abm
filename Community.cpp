// Community.cpp

#include <cstdlib>
#include <cstring>
#include <climits>
#include <iostream>
#include <string>
#include <iterator>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Person.h"
#include "Mosquito.h"
#include "Location.h"
#include "Community.h"
#include "Parameters.h"

using namespace dengue::standard;

const Parameters* Community::_par;
map< Location*, map<int, bool> > Community::_isHot;

// Community
Community::Community(const Parameters* parameters) :
    _exposedQueue(MAX_INCUBATION, vector<Person*>(0)),
    _infectiousMosquitoQueue(MAX_MOSQUITO_AGE-MOSQUITO_INCUBATION, vector<Mosquito*>(0)),
    _exposedMosquitoQueue(MOSQUITO_INCUBATION, vector<Mosquito*>(0)),
    _nNumNewlyInfected(NUM_OF_SEROTYPES, vector<int>(MAX_RUN_TIME)),
    _nNumNewlySymptomatic(NUM_OF_SEROTYPES, vector<int>(MAX_RUN_TIME))
    {
    _par = parameters;
    _nDay = 0;
    _nNumPerson = 0;
    _person = NULL;
    _fMosquitoCapacityMultiplier = 1.0;
    _fMortality = NULL;
    _bNoSecondaryTransmission = false;
}


Community::~Community() {
    if (_person)
        delete [] _person;

    _location.clear();
    _exposedQueue.clear();
    _infectiousMosquitoQueue.clear();
    _exposedMosquitoQueue.clear();
    _personAgeCohort.clear();
    _numLocationMosquitoCreated.clear();
}


bool Community::loadPopulation(string szPop,string szImm) {
    _szPopulationFilename = szPop;
    _szImmunityFilename = szImm;
    ifstream iss(_szPopulationFilename.c_str());
    if (!iss) {
        cerr << "ERROR: " << _szPopulationFilename << " not found." << endl;
        return false;
    }
    // count lines
    int maxPerson = 0;
    while (iss) {
        char buffer[500];
        iss.getline(buffer,500);
        maxPerson++;
    }
    iss.close();
    iss.open(_szPopulationFilename.c_str());
    _nNumPerson=0;
    _person = new Person[maxPerson];
    int agecounts[MAX_PERSON_AGE];
    for (int i=0; i<MAX_PERSON_AGE; i++)
        agecounts[i] = 0;
    while (iss) {
        char buffer[500];
        iss.getline(buffer,500);
        istringstream line(buffer);
        int id, age, house, work;
        int temp1, temp2;
        string gender;
        if (line >> id >> house >> age >> gender >> temp1 >> temp2 >> work) {
            _person[_nNumPerson].setAge(age);
            _person[_nNumPerson].setHomeID(house);
            _person[_nNumPerson].setLocation(_location[house], 0);
            _person[_nNumPerson].setLocation(_location[work], 1);
            _person[_nNumPerson].setLocation(_location[house], 2);
            _location[house]->addPerson(_person+_nNumPerson, 0);
            _location[work]->addPerson(_person+_nNumPerson, 1);
            _location[house]->addPerson(_person+_nNumPerson, 2);
            assert(age<MAX_PERSON_AGE);
            agecounts[age]++;
            _nNumPerson++;
        }
        else {
            //      cerr << "error" << endl;
            //      return false;
            //      cerr << "header" << endl;
        }
    }
    iss.close();

    if (_szImmunityFilename.length()>0) {
        ifstream immiss(_szImmunityFilename.c_str());
        if (!immiss) {
            cerr << "ERROR: " << _szImmunityFilename << " not found." << endl;
            return false;
        }
        char buffer[500];
        int part;
        vector<int> parts;
        while (immiss) {
            immiss.getline(buffer,500);
            istringstream line(buffer);
            while (line >> part) parts.push_back(part);

            if (parts.size() >= 1 + NUM_OF_SEROTYPES) {
                int id = parts[0];
                for (int i=id-1; i<=id; i++) {
                    if (_person[i].getID()==id) {
                        for (unsigned int s=0; s<NUM_OF_SEROTYPES; s++) {
                            if (parts[s+1]>0) _person[i].setImmunity((Serotype) s);
                        }
                        break;
                    }
                }
            }
            parts.clear();
        }
        immiss.close();
    }

    // keep track of all age cohorts for aging and mortality
    _personAgeCohort.clear();
    _personAgeCohort.resize(MAX_PERSON_AGE, vector<Person*>(0));

    for (int i=0; i<_nNumPerson; i++) {
        int age = _person[i].getAge();
        assert(age<MAX_PERSON_AGE);
        _personAgeCohort[age].push_back(_person + i);
        _nPersonAgeCohortSizes[age]++;
    }
    return true;
}


bool Community::loadLocations(string szLocs,string szNet) {
    _szLocationFilename = szLocs;
    _szNetworkFilename = szNet;
    ifstream iss(_szLocationFilename.c_str());
    if (!iss) {
        cerr << "ERROR: " << _szLocationFilename << " not found." << endl;
        return false;
    }
    int maxLocationID = 0;
    while (iss) {
        char buffer[500];
        iss.getline(buffer,500);
        istringstream line(buffer);
        int locID;
        if (line >> locID) {
            if (locID>maxLocationID) maxLocationID = locID;
        }
    }
    iss.close();
    maxLocationID+=1;  // we don't actually use index 0 in the location vector . . .
    cerr << maxLocationID << " locations" << endl;
    _location.clear();
    _location.resize(maxLocationID);
    for (unsigned int i=0; i<_location.size(); i++) {
        _location[i] = new Location();
        _location[i]->setBaseMosquitoCapacity(_par->nDefaultMosquitoCapacity); 
        _isHot[_location[i]] = map<int,bool>();
    } 

    _numLocationMosquitoCreated.clear();
    _numLocationMosquitoCreated.resize(maxLocationID, vector<int>(MAX_RUN_TIME, 0));

    iss.open(_szNetworkFilename.c_str());
    if (!iss) {
        cerr << "ERROR: " << _szNetworkFilename << " not found." << endl;
        return false;
    }
    while (iss) {
        char buffer[500];
        iss.getline(buffer,500);
        istringstream line(buffer);
        int locID1,locID2;
        if (line >> locID1 >> locID2) {
            //      cerr << locID1 << " , " << locID2 << endl;
            _location[locID1]->addNeighbor(_location[locID2]);            // should check for ID
            _location[locID2]->addNeighbor(_location[locID1]);
        }
        else {
            //      cerr << "header" << endl;
        }
    }
    iss.close();
    for (unsigned int i=0; i<_location.size(); i++)
        if (_location[i]->getNumNeighbors()<=0)
            _location[i]->setUndefined();
    return true;
}


// infect - infects person id
bool Community::infect(int id, Serotype serotype, int day) {
    for (int i=0; i<_nNumPerson; i++)
    if (_person[i].getID()==id) {
        bool result =  _person[i].infect(-1, serotype, day, 0);
        if (result)
            _nNumNewlyInfected[(int) serotype][_nDay]++;
        return result;
        //      cerr << "inf " << i << " at "  << _person[i].getLocation(0)->getID() << endl;
    }
    return false;
}


// vaccinate - vaccinate fraction f of the population
// if age>=0, then vaccinate only those who are "age" years old
void Community::vaccinate(double f, int age) {
    if (f<=0.0)
        return;
    if (age<0) {                                                      // vaccinate everyone
        for (int i=0; i<_nNumPerson; i++)
            if (gsl_rng_uniform(RNG)<f)
                _person[i].vaccinate();
    }
    else {
        for (int pnum=0; pnum<_nPersonAgeCohortSizes[age]; pnum++) {
            Person *p = _personAgeCohort[age][pnum];
            assert(p!=NULL);
            if (!p->isVaccinated() && gsl_rng_uniform(RNG)<f)
                p->vaccinate();
        }
    }
}


// returns number of days mosquito has left to live
int Community::addMosquito(Location *p, Serotype serotype, int nInfectedByID) {
    Mosquito *m = new Mosquito(p, serotype, nInfectedByID);
    int daysleft = m->getAgeDeath() - m->getAgeInfected();
    int daysinfectious = daysleft - MOSQUITO_INCUBATION;
    if (daysinfectious<=0) {
        delete m;
        return daysleft;                                              // dies before infectious
    }

    //on the latest possible day
    // add mosquito to latency queue
    _exposedMosquitoQueue.back().push_back(m); 
    return daysleft;
}


int Community::getNumInfectiousMosquitoes() {
    int count = 0;
    for (unsigned int i=0; i<_infectiousMosquitoQueue.size(); i++) {
        count += _infectiousMosquitoQueue[i].size();
    }
    return count;
}


int Community::getNumExposedMosquitoes() {
    int count = 0;
    for (unsigned int i=0; i<_exposedMosquitoQueue.size(); i++) {
        count += _exposedMosquitoQueue[i].size();
    }
    return count;
}


Mosquito *Community::getInfectiousMosquito(int n) {
    for (unsigned int i=0; i<_infectiousMosquitoQueue.size(); i++) {
        int bin_size = _infectiousMosquitoQueue[i].size();
        if (n >= bin_size) {
            n -= bin_size;
        } else {
            return _infectiousMosquitoQueue[i][n]; 
        }
    }
    return NULL;
}


Mosquito *Community::getExposedMosquito(int n) {
    for (unsigned int i=0; i<_exposedMosquitoQueue.size(); i++) {
        int bin_size = _exposedMosquitoQueue[i].size();
        if (n >= bin_size) {
            n -= bin_size;
        } else {
            return _exposedMosquitoQueue[i][n]; 
        }
    }
    return NULL;
}


void Community::moveMosquito(Mosquito *m) {
    double r = gsl_rng_uniform(RNG);
    if (r<_par->fMosquitoMove) {
        if (r<_par->fMosquitoTeleport) {                               // teleport
            int locID;
            do {
                locID = gsl_rng_uniform_int(RNG,_location.size());
            } while (_location[locID]->getUndefined());
            m->setLocation(_location[locID]);
        }                                                             // move to neighbor
        else {
            Location *pLoc = m->getLocation();
            if (pLoc->getNumNeighbors()>0) {
                int n = gsl_rng_uniform_int(RNG,pLoc->getNumNeighbors());
                m->setLocation(pLoc->getNeighbor(n));
            }
        }
    }
}


void Community::swapImmuneStates() {
    map< Location*, map<int,bool> >::iterator it;
    for ( it=_isHot.begin() ; it != _isHot.end(); it++ ) (*it).second.clear();
    //    cerr << "swap!" << endl;
    // big annual swap!
    // For people of age x, copy immune status from people of age x-1
    for (int age=MAX_PERSON_AGE-1; age>0; age--) {
        //      cerr << "age " << age << ", " << _nPersonAgeCohortSizes[age] << " people" << endl;
        int age1 = age-1;
        for (int pnum=0; pnum<_nPersonAgeCohortSizes[age]; pnum++) {
            Person *p = _personAgeCohort[age][pnum];
            assert(p!=NULL);
            int r = gsl_rng_uniform_int(RNG,_nPersonAgeCohortSizes[age1]);
            p->copyImmunity(_personAgeCohort[age1][r]);

            // update map of locations with infectious people
            if (p->getRecoveryTime() > _nDay) {
                for (int d = p->getInfectiousTime(); d < p->getRecoveryTime(); d++) {
                    for (int t=0; t<STEPS_PER_DAY; t++) {
                        flagInfectedLocation(p->getLocation(t), d);
                    }
                }
            }
        }
    }
    // For people of age 0, reset immunity
    for (int pnum=0; pnum<_nPersonAgeCohortSizes[0]; pnum++) {
        Person *p = _personAgeCohort[0][pnum];
        assert(p!=NULL);
        p->resetImmunity();
    }
    return;
}


void Community::updateWithdrawnStatus() {
    for (int i=0; i<_nNumPerson; i++) {
        Person *p = _person+i;
        if (p->getSymptomTime()==_nDay)                               // started showing symptoms
            _nNumNewlySymptomatic[(int) p->getSerotype()][_nDay]++;
        if (p->getWithdrawnTime()==_nDay) {                           // started withdrawing
            p->getLocation(0)->addPerson(p,1);                        // stays at home at mid-day
            p->getLocation(1)->removePerson(p,1);                     // does not go to work
        }
        if (_nDay>0 &&
            p->isWithdrawn(_nDay-1) &&
        p->getRecoveryTime()==_nDay) {                                // just stopped withdrawing
            p->getLocation(1)->addPerson(p,1);                        // goes back to work
            p->getLocation(0)->removePerson(p,1);                     // stops staying at home
        }
    }
    return;
}


void Community::mosquitoToHumanTransmission() {
    for(unsigned int i=0; i<_infectiousMosquitoQueue.size(); i++) {
        for(unsigned int j=0; j<_infectiousMosquitoQueue[i].size(); j++) {
            Mosquito* m = _infectiousMosquitoQueue[i][j];
            Location *pLoc = m->getLocation();
            if (gsl_rng_uniform(RNG)<_par->betaMP) {                      // infectious mosquito bites

                // take sum of people in the location, weighting by time of day
                double exposuretime[STEPS_PER_DAY];
                double totalExposureTime = 0;
                for (int timeofday=0; timeofday<STEPS_PER_DAY; timeofday++) {
                    exposuretime[timeofday] = pLoc->getNumPerson(timeofday) * DAILY_BITING_PDF[timeofday];
                    totalExposureTime += exposuretime[timeofday];
                }
                double r = gsl_rng_uniform(RNG) * totalExposureTime;
                for (int timeofday=0; timeofday<STEPS_PER_DAY; timeofday++) {
                    if (r<exposuretime[timeofday]) {
                        // bite at this time of day
                        int idx = floor(r*pLoc->getNumPerson(timeofday)/exposuretime[timeofday]);
                        Person *p = pLoc->getPerson(idx, timeofday);
//                        cerr << " infect person " << idx << "/" << pLoc->getNumPerson(timeofday) 
//                            << "  " << p->getID() << " with " << m->getSerotype() << " from " << m->getID() << endl;
                        Serotype serotype = m->getSerotype();
                        if (p->infect(m->getID(), serotype, _nDay, pLoc->getID())) {
                            _nNumNewlyInfected[(int) serotype][_nDay]++;
                            if (_bNoSecondaryTransmission) {
                                p->kill(_nDay);                       // kill secondary cases so they do not transmit
                            }
                            else {
                                // NOTE: We are storing the location ID of infection, not person ID!!!
                                //	      cerr << "Person " << p->getID() << " infected" << endl;
                                // add to queue
                                _exposedQueue[p->getInfectiousTime()-_nDay].push_back(p);
                            }
                        }
                        break;
                    }
                    r -= exposuretime[timeofday];
                }
            }
        }
    }
    return;
}


void Community::humanToMosquitoTransmission() {
    for (unsigned int loc=0; loc<_location.size(); loc++) {
        if (!_location[loc]->getUndefined() && _isHot[_location[loc]].count(_nDay)) {
            _isHot[_location[loc]].erase(_nDay);  // works if transmission is modeled daily
            double sumviremic = 0.0;
            double sumnonviremic = 0.0;
            vector<double> sumserotype(NUM_OF_SEROTYPES,0.0);                                    // serotype fractions at location

            // calculate fraction of people who are viremic
            for (int timeofday=0; timeofday<STEPS_PER_DAY; timeofday++) {
                for (int i=_location[loc]->getNumPerson(timeofday)-1; i>=0; i--) {
                    Person *p = _location[loc]->getPerson(i, timeofday);
                    if (p->isViremic(_nDay)) {
                        double vaceffect = (p->isVaccinated()?(1.0-_par->fVEI):1.0);
                        int serotype = (int) p->getSerotype();
                        if (vaceffect==1.0) {
                            sumviremic += DAILY_BITING_PDF[timeofday];
                            sumserotype[serotype] += DAILY_BITING_PDF[timeofday];
                        } else {
                            sumviremic += DAILY_BITING_PDF[timeofday]*vaceffect;
                            sumserotype[serotype] += DAILY_BITING_PDF[timeofday]*vaceffect;
                            // a vaccinated person is treated like a fraction of an infectious person and a fraction of a non-infectious person
                            sumnonviremic += DAILY_BITING_PDF[timeofday]*(1.0-vaceffect);
                        }
                    } else {
                        sumnonviremic += DAILY_BITING_PDF[timeofday];
                    }
                }
            }

            if (sumviremic>0.0) {
                for (int i=0; i<NUM_OF_SEROTYPES; i++)
                    sumserotype[i] /= sumviremic;
                int m;                                                // number of susceptible mosquitoes
                int locid = _location[loc]->getID();                   // location ID
                m = int(_location[loc]->getBaseMosquitoCapacity()*_fMosquitoCapacityMultiplier+0.5);
                //if (_numLocationMosquitoCreated[locid])
                    m -= _numLocationMosquitoCreated[locid][_nDay];
                if (m<0)
                    m=0;
                                                                      // how many susceptible mosquitoes bite viremic hosts in this location?
                int numbites = gsl_ran_binomial(RNG, _par->betaPM*sumviremic/(sumviremic+sumnonviremic), m);
                while (numbites-->0) {
                    int serotype;                                     // which serotype infects mosquito
                    if (sumserotype[0]==1.0) {
                        serotype = 0;
                    } else {
                        double r = gsl_rng_uniform(RNG);
                        for (serotype=0; serotype<NUM_OF_SEROTYPES && r>sumserotype[serotype]; serotype++)
                            r -= sumserotype[serotype];
                    }
                    /*if (serotype==0) {
                        cerr << "serotypes: ";
                        for (int i=0; i<5; i++) cerr << sumserotype[i] << ",";
                        cerr << endl;
                    }*/
                    int daysleft = addMosquito(_location[loc], (Serotype) serotype, locid);
                    assert(_nDay+daysleft<MAX_RUN_TIME);
                    for (int j=0; j<daysleft; j++)
                        _numLocationMosquitoCreated[locid][_nDay+j]++;
                }
            }
            //	  cerr << "Mosquito infected by Person " << _person[i].getID() << " at " 
            //         << _person[i].getLocation(timeofday)->getID() << "." << endl;
        }
    }
    return;
}


void Community::_advanceTimers() {
    // advance incubation in people
    for (unsigned int i=0; i<_exposedQueue.size()-1; i++) {
        _exposedQueue[i] = _exposedQueue[i+1];
    }
    _exposedQueue.back().clear();


    // advance age of infectious mosquitoes
    for (unsigned int i=0; i<_infectiousMosquitoQueue.size()-1; i++) {
        _infectiousMosquitoQueue[i] = _infectiousMosquitoQueue[i+1];
    }
    _infectiousMosquitoQueue.back().clear();

    // advance incubation period of exposed mosquitoes
    for (unsigned int mnum=0; mnum<_exposedMosquitoQueue[0].size(); mnum++) {
        Mosquito *m = _exposedMosquitoQueue[0][mnum];
        // incubation over: some mosquitoes become infectious
        int daysinfectious = m->getAgeDeath() - m->getAgeInfected() - MOSQUITO_INCUBATION;
        assert((unsigned) daysinfectious < _infectiousMosquitoQueue.size());
        _infectiousMosquitoQueue[daysinfectious].push_back(m);
    }

    for (unsigned int i=0; i<_exposedMosquitoQueue.size()-1; i++) {
        _exposedMosquitoQueue[i] = _exposedMosquitoQueue[i+1];
    }
    _exposedMosquitoQueue.back().clear();
    
    return;
}


void Community::_modelMosquitoMovement() {
    // move mosquitoes
    for(unsigned int i=0; i<_infectiousMosquitoQueue.size(); i++) {
        for(unsigned int j=0; j<_infectiousMosquitoQueue[i].size(); j++) {
            Mosquito* m = _infectiousMosquitoQueue[i][j];
            moveMosquito(m);
        }
    }
    for(unsigned int i=0; i<_exposedMosquitoQueue.size(); i++) {
        for(unsigned int j=0; j<_exposedMosquitoQueue[i].size(); j++) {
            Mosquito* m = _exposedMosquitoQueue[i][j];
            moveMosquito(m);
        }
    }
    return;
}


void Community::tick() {
    assert(_nDay<MAX_RUN_TIME);
    if (_nDay%365==364) { swapImmuneStates(); }                       // randomize and advance immune states
    updateWithdrawnStatus();                                          // make people stay home or return to work
    mosquitoToHumanTransmission();                                    // infect people

    // maybe add the occasional random infection of a person to reflect sporadic travel
    humanToMosquitoTransmission();                                    // infect mosquitoes in each location
    _advanceTimers();                                                 // advance H&M incubation periods and M ages
    _modelMosquitoMovement();                                         // probabilistic movement of mosquitos

    _nDay++;
    return;
}


// getNumInfected - counts number of infected residents
int Community::getNumInfected(int day) {
    int count=0;
    for (int i=0; i<_nNumPerson; i++)
        if (_person[i].isInfected(day))
            count++;
    return count;
}


// getNumSymptomatic - counts number of symptomatic residents
int Community::getNumSymptomatic(int day) {
    int count=0;
    for (int i=0; i<_nNumPerson; i++)
        if (_person[i].isSymptomatic(day))
            count++;
    return count;
}


// getNumSusceptible - counts number of susceptible residents
vector<int> Community::getNumSusceptible() {
    vector<int> counts(NUM_OF_SEROTYPES, 0);
    //int count=0;
    for (int i=0; i<_nNumPerson; i++) {
        for (int s=0; s<NUM_OF_SEROTYPES; s++) {
            if (_person[i].isSusceptible((Serotype) s))
                counts[s]++;
        }
    }
    return counts;
}

