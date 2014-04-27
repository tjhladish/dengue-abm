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
    _infectiousMosquitoQueue(MAX_MOSQUITO_AGE, vector<Mosquito*>(0)),
    _exposedMosquitoQueue(MAX_MOSQUITO_INCUBATION, vector<Mosquito*>(0)),
    _nNumNewlyInfected(NUM_OF_SEROTYPES, vector<int>(MAX_RUN_TIME)),
    _nNumNewlySymptomatic(NUM_OF_SEROTYPES, vector<int>(MAX_RUN_TIME))
    {
    _par = parameters;
    _nDay = 0;
    _nNumPerson = 0;
    _person = NULL;
    _fMosquitoCapacityMultiplier = 1.0;
    _nExternalIncubation = 11; // default external incubation period of 11 days (Nishiura & Halstead 2007)
    _fMortality = NULL;
    _bNoSecondaryTransmission = false;
    _uniformSwap = true;
    for (int a = 0; a<MAX_PERSON_AGE; a++) _nPersonAgeCohortSizes[a] = 0;
}


Community::~Community() {
    if (_person)
        delete [] _person;

    _location.clear();
    _exposedQueue.clear();
    _infectiousMosquitoQueue.clear();
    _exposedMosquitoQueue.clear();
    _personAgeCohort.clear();
}


bool Community::loadPopulation(string populationFilename, string immunityFilename, string swapFilename) {
    ifstream iss(populationFilename.c_str());
    if (!iss) {
        cerr << "ERROR: " << populationFilename << " not found." << endl;
        return false;
    }
    // count lines
    int maxPerson = -1; // header line isn't a person
    char buffer[500];
    while (iss) {
        iss.getline(buffer,500);
        if (iss) maxPerson++;
    }
    iss.close();

    iss.open(populationFilename.c_str());
    _nNumPerson=0;
    _person = new Person[maxPerson];

    int agecounts[MAX_PERSON_AGE];
    for (int i=0; i<MAX_PERSON_AGE; i++)
        agecounts[i] = 0;

    istringstream line;
    int id, age, house, work;
    string hh_serial;
    int pernum;
    string gender;  // might be "M", "F", "1", or "2"
    while (iss) {
        iss.getline(buffer,500);
        line.clear(); 
        line.str(buffer);
        /*
        pid hid age sex hh_serial pernum workid
        1 1 31 1 2748179000 1 442670
        2 1 29 2 2748179000 2 395324
        3 1 10 2 2748179000 3 468423
        4 2 32 1 2748114000 1 397104
        5 2 30 2 2748114000 2 396166
        */
        if (line >> id >> house >> age >> gender >> hh_serial >> pernum >> work) {
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
    }
    iss.close();

    if (immunityFilename.length()>0) {
        ifstream immiss(immunityFilename.c_str());
        if (!immiss) {
            cerr << "ERROR: " << immunityFilename << " not found." << endl;
            return false;
        }
        int part;
        vector<int> parts;
        istringstream line;
        while (immiss) {
            immiss.getline(buffer,500);
            line.clear(); 
            line.str(buffer);
            while (line >> part) parts.push_back(part);

            if (parts.size() >= 1 + NUM_OF_SEROTYPES) {
                int id = parts[0];
                for (int i=id-1; i<=id; i++) {
                    if (_person[i].getID()==id) {
                        for (unsigned int s=0; s<NUM_OF_SEROTYPES; s++) {
                            if (parts[s+1]>0) {
                                _person[i].setImmunity((Serotype) s);
                                _person[i].initializeNewInfection();
                                _person[i].setRecoveryTime(-365*parts[s+1]); // last dengue infection was x years ago
                            }
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

    //cerr << "cohort sizes: " << _nPersonAgeCohortSizes[0] << " " << _nPersonAgeCohortSizes[25] << " " << _nPersonAgeCohortSizes[75] << " " <<    _nPersonAgeCohortSizes[100] << endl;

    if (swapFilename == "") {
        _uniformSwap = true;
    } else {
        iss.open(swapFilename.c_str());
        if (!iss) {
            cerr << "ERROR: " << swapFilename << " not found." << endl;
            return false;
        }

        int id1, id2;
        double prob;
        istringstream line;
        while (iss) {
            iss.getline(buffer,500);
            line.clear(); 
            line.str(buffer);

            if (line >> id1 >> id2 >> prob) {
                Person* person = getPersonByID(id1);
                person->appendToSwapProbabilities(make_pair(id2, prob));
            }
        }
        iss.close();
        _uniformSwap = false;
    }

    return true;
}


bool Community::loadLocations(string locationFilename,string networkFilename) {
    ifstream iss(locationFilename.c_str());
    if (!iss) {
        cerr << "ERROR: " << locationFilename << " not found." << endl;
        return false;
    }
    _location.clear();
 
    // This is a hack for backward compatibility.  Indices should start at zero.
    Location* dummy = new Location();
    dummy->setBaseMosquitoCapacity(_par->nDefaultMosquitoCapacity); 
    _isHot[dummy] = map<int,bool>();
    _location.push_back(dummy); // first val is a dummy, for backward compatibility
    // End of hack
    char buffer[500];
    int locID;
    string locType;
    double locX, locY;
    istringstream line(buffer);

    while (iss) {
        iss.getline(buffer,500);
        line.clear();
        line.str(buffer);
        if (line >> locID >> locType >> locX >> locY) {
            if (locID != (signed) _location.size()) {
                cerr << "WARNING: Location ID's must be sequential integers" << endl;
                return false;
            }
            Location* newLoc = new Location();
            newLoc->setID(locID);
            newLoc->setX(locX);
            newLoc->setY(locY);
            if (_par->eMosquitoDistribution==CONSTANT) {
                // all houses have same number of mosquitoes
                newLoc->setBaseMosquitoCapacity(_par->nDefaultMosquitoCapacity); 
            } else if (_par->eMosquitoDistribution==EXPONENTIAL) {
                // exponential distribution of mosquitoes -dlc
                // gsl takes the 1/lambda (== the expected value) as the parameter for the exp RNG
                newLoc->setBaseMosquitoCapacity(gsl_ran_exponential(RNG, _par->nDefaultMosquitoCapacity));
            } else {
                cerr << "ERROR: Invalid mosquito distribution: " << _par->eMosquitoDistribution << endl;
                cerr << "       Valid distributions include CONSTANT and EXPONENTIAL" << endl;
                return false;
            }

            _isHot[newLoc] = map<int,bool>(); // _isHot flags locations with infections
            _location.push_back(newLoc);
        }
    }
    iss.close();
    //cerr << _location.size() << " locations" << endl;

    //_numLocationMosquitoCreated.clear();
    //_numLocationMosquitoCreated.resize(_location.size(), vector<int>(MAX_RUN_TIME, 0));

    iss.open(networkFilename.c_str());
    if (!iss) {
        cerr << "ERROR: " << networkFilename << " not found." << endl;
        return false;
    }
    int locID1, locID2;
    //istringstream line(buffer);
    while (iss) {
        iss.getline(buffer,500);
        line.clear();
        line.str(buffer);
        if (line >> locID1 >> locID2) { // data (non-header) line
            //      cerr << locID1 << " , " << locID2 << endl;
            _location[locID1]->addNeighbor(_location[locID2]);            // should check for ID
            _location[locID2]->addNeighbor(_location[locID1]);
        }
    }
    iss.close();

    for (unsigned int i=0; i<_location.size(); i++)
        if (_location[i]->getNumNeighbors()<=0)
            _location[i]->setUndefined();
    return true;
}


Person* Community::getPersonByID(int id) {
    // This assumes that IDs start at 1, and tries to guess
    // that person with ID id is in position id-1
    int i = 0;
    if (_person[id-1].getID()==id) {
        i = id-1;
    } else {
        for (i=0; i<_nNumPerson; i++) {
            if (_person[i].getID()==id) {
                break;
            }
        }
    }
    return &_person[i];
}


// infect - infects person id
bool Community::infect(int id, Serotype serotype, int day) {
    Person* person = getPersonByID(id);

    bool result =  person->infect(-1, serotype, day, 0);
    if (result) _nNumNewlyInfected[(int) serotype][_nDay]++;
    return result;
}


// vaccinate - vaccinate fraction f of the population
// if age>=0, then vaccinate only those who are "age" years old
void Community::vaccinate(double f, int age) {
    if (f<=0.0)
        return;
    if (age<0) {                                                      // vaccinate everyone
        for (int i=0; i<_nNumPerson; i++) {
            if (gsl_rng_uniform(RNG)<f) {
                _person[i].vaccinate();
            }
        }
    }
    else {
        if (age>MAX_PERSON_AGE)
            age=MAX_PERSON_AGE;
        for (int pnum=0; pnum<_nPersonAgeCohortSizes[age]; pnum++) {
            Person *p = _personAgeCohort[age][pnum];
            assert(p!=NULL);
            if (!p->isVaccinated() && gsl_rng_uniform(RNG)<f) {
                p->vaccinate();
            }
        }
    }
}


// returns number of days mosquito has left to live
int Community::attemptToAddMosquito(Location *p, Serotype serotype, int nInfectedByID) {
    int eip = getExternalIncubation();
    Mosquito *m = new Mosquito(p, serotype, nInfectedByID, eip);
    int daysleft = m->getAgeDeath() - m->getAgeInfected();
    int daysinfectious = daysleft - eip;
    if (daysinfectious<=0) {
        delete m;
        return daysleft;                                              // dies before infectious
    }

    //on the latest possible day
    // add mosquito to latency queue
    //    _exposedMosquitoQueue.back().push_back(m);
    _exposedMosquitoQueue[eip-1].push_back(m);
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
            } while (_location[locID]->getUndefined());               // why would it be undefined?
            m->updateLocation(_location[locID]);
        } else {                                                            // move to neighbor
            Location *pLoc = m->getLocation();
            double x1 = pLoc->getX();
            double y1 = pLoc->getY();

            int degree = pLoc->getNumNeighbors();
            int neighbor=0; // neighbor is an index

            if (_par->szMosquitoMoveModel == "weighted") {
                vector<double> weights(degree, 0);
                double sum_weights = 0.0;

                // Prefer nearby neighbors
                // Calculate distance-based weights to select each of the degree neighbors
                for (int i=0; i<degree; i++) {
                    Location* loc2 = pLoc->getNeighbor(i);
                    double x2 = loc2->getX();
                    double y2 = loc2->getY();
                    double distance_squared = pow(x1-x2,2) + pow(y1-y2,2);
                    double w = 1.0 / distance_squared;
                    sum_weights += w;
                    weights[i] = w;
                }
                double r2 = gsl_rng_uniform(RNG);
                neighbor = degree-1; // neighbor is (still) an index
                int idx;
                for ( idx = 0; idx < degree - 1; idx++ ) {
                    weights[idx] /= sum_weights; // normalize prob
                    if ( r2 < weights[idx] ) {
                        break;
                    } else {
                        r2 -= weights[idx];
                    }
                }
                neighbor = idx; 
            } else {
                // Ignore actual distances
                if (degree>0) {
                    neighbor = gsl_rng_uniform_int(RNG,pLoc->getNumNeighbors());
                }
            }

            m->updateLocation(pLoc->getNeighbor(neighbor));
        }    
    }
}


void Community::swapImmuneStates() {
    map< Location*, map<int,bool> >::iterator it;
    for ( it=_isHot.begin() ; it != _isHot.end(); it++ ) (*it).second.clear();

    // For people of age x, copy immune status from people of age x-1
    for (int age=MAX_PERSON_AGE-1; age>0; age--) {
        for (int pnum=0; pnum<_nPersonAgeCohortSizes[age]; pnum++) {
            //cerr << "age " << age << ": " << pnum << " of " << _nPersonAgeCohortSizes[age] << endl;
            Person *p = _personAgeCohort[age][pnum];
            assert(p!=NULL);
            if (_uniformSwap == true) {
                // For people of age x, copy immune status from people of age x-1
                int r = gsl_rng_uniform_int(RNG,_nPersonAgeCohortSizes[age-1]);
                p->copyImmunity(_personAgeCohort[age-1][r]);
            } else {
                // Same as above, but use weighted sampling based on swap probs from file
                double r = gsl_rng_uniform(RNG);
                const vector<pair<int, double> >& swap_probs = p->getSwapProbabilities();
                int n;
                for (n = 0; n < (signed) swap_probs.size() - 1; n++) {
                    if (r < swap_probs[n].second) {
                        break;
                    } else {
                        r -= swap_probs[n].second;
                    }
                }
                const int id = swap_probs[n].first;
                p->copyImmunity(getPersonByID(id)); 
            }

            // update map of locations with infectious people
            if (p->getNumInfections() > 0 and p->getRecoveryTime() > _nDay) {
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
        if (p->getNumInfections() == 0) continue;
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
                for (int t=0; t<STEPS_PER_DAY; t++) {
                    exposuretime[t] = pLoc->getNumPerson(t) * DAILY_BITING_PDF[t];
                    totalExposureTime += exposuretime[t];
                }
                if ( totalExposureTime > 0 ) {
                    double r = gsl_rng_uniform(RNG) * totalExposureTime;
                    int timeofday;
                    for (timeofday=0; timeofday<STEPS_PER_DAY - 1; timeofday++) {
                        if (r<exposuretime[timeofday]) {
                            // bite at this time of day
                            break;
                        }
                        r -= exposuretime[timeofday];
                    }
                    int idx = floor(r*pLoc->getNumPerson(timeofday)/exposuretime[timeofday]);
                    Person *p = pLoc->getPerson(idx, timeofday);
                    Serotype serotype = m->getSerotype();
                    if (p->infect(m->getID(), serotype, _nDay, pLoc->getID())) {
                        _nNumNewlyInfected[(int) serotype][_nDay]++;
                        if (_bNoSecondaryTransmission) {
                            p->kill(_nDay);                       // kill secondary cases so they do not transmit
                        }
                        else {
                            // NOTE: We are storing the location ID of infection, not person ID!!!
                            // add to queue
                            _exposedQueue[p->getInfectiousTime()-_nDay].push_back(p);
                        }
                    }
                }
            }
        }
    }
    return;
}


void Community::humanToMosquitoTransmission() {
    for (unsigned int loc_idx=0; loc_idx<_location.size(); loc_idx++) {
        Location* loc = _location[loc_idx];
        if (!loc->getUndefined() && _isHot[loc].count(_nDay)) {
            _isHot[loc].erase(_nDay);  // works if transmission is modeled daily
            double sumviremic = 0.0;
            double sumnonviremic = 0.0;
            vector<double> sumserotype(NUM_OF_SEROTYPES,0.0);                                    // serotype fractions at location

            // calculate fraction of people who are viremic
            for (int timeofday=0; timeofday<STEPS_PER_DAY; timeofday++) {
                for (int i=loc->getNumPerson(timeofday)-1; i>=0; i--) {
                    Person *p = loc->getPerson(i, timeofday);
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
                for (int i=0; i<NUM_OF_SEROTYPES; i++) {
                    sumserotype[i] /= sumviremic;
                }
                int locid = loc->getID();                   // location ID
                int m = int(loc->getBaseMosquitoCapacity()*_fMosquitoCapacityMultiplier+0.5);  // number of mosquitoes
                m -= loc->getCurrentInfectedMosquitoes(); // subtract off the number of already-infected mosquitos
                if (m<0) // should never happen, but set minimum number of susceptible mosquitos to 0
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
                    int daysleft = attemptToAddMosquito(loc, (Serotype) serotype, locid);
                    if (daysleft > 0) loc->addInfectedMosquito();
                }
            }
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

    // delete infected mosquitoes that are dying today
    vector<Mosquito*>::iterator itr; 
    for(itr = _infectiousMosquitoQueue.back().begin(); itr != _infectiousMosquitoQueue.back().end(); ++itr ) {
        delete (*itr);
    }
    _infectiousMosquitoQueue.back().clear();

    // advance incubation period of exposed mosquitoes
    for (unsigned int mnum=0; mnum<_exposedMosquitoQueue[0].size(); mnum++) {
        Mosquito *m = _exposedMosquitoQueue[0][mnum];
        // incubation over: some mosquitoes become infectious
        int daysinfectious = m->getAgeDeath() - m->getAgeInfectious(); // - MOSQUITO_INCUBATION;
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


void Community::tick(int day) {
    _nDay = day;
    assert(_nDay<MAX_RUN_TIME);
    if ((_nDay-100)%365==364) { swapImmuneStates(); }                 // randomize and advance immune states
    updateWithdrawnStatus();                                          // make people stay home or return to work
    mosquitoToHumanTransmission();                                    // infect people

    // maybe add the occasional random infection of a person to reflect sporadic travel
    humanToMosquitoTransmission();                                    // infect mosquitoes in each location
    _advanceTimers();                                                 // advance H&M incubation periods and M ages
    _modelMosquitoMovement();                                         // probabilistic movement of mosquitos
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
