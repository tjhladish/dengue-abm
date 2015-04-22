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
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Person.h"
#include "Mosquito.h"
#include "Location.h"
#include "Community.h"
#include "Parameters.h"

using namespace dengue::standard;

const Parameters* Community::_par;
vector< unordered_set<Location*> > Community::_isHot;

Community::Community(const Parameters* parameters) :
    _exposedQueue(MAX_INCUBATION, vector<Person*>(0)),
    _infectiousMosquitoQueue(MAX_MOSQUITO_AGE, vector<Mosquito*>(0)),
    // reserving MAX_MOSQUITO_AGE is simpler than figuring out what the maximum
    // possible EIP is when EIP is variable
    _exposedMosquitoQueue(MAX_MOSQUITO_AGE, vector<Mosquito*>(0)),
    _nNumNewlyInfected(NUM_OF_SEROTYPES, vector<int>(parameters->nRunLength + MAX_MOSQUITO_AGE)),
    _nNumNewlySymptomatic(NUM_OF_SEROTYPES, vector<int>(parameters->nRunLength + MAX_MOSQUITO_AGE)),
    _nNumVaccinatedCases(NUM_OF_SEROTYPES, vector<int>(parameters->nRunLength + MAX_MOSQUITO_AGE))
    {
    _par = parameters;
    _nDay = 0;
    _nNumPerson = 0;
    _person = NULL;
    _fMosquitoCapacityMultiplier = 1.0;
    _EIP = 11; // default external incubation period of 11 days (Nishiura & Halstead 2007)
    _fMortality = NULL;
    _bNoSecondaryTransmission = false;
    _uniformSwap = true;
    for (int a = 0; a<NUM_AGE_CLASSES; a++) _nPersonAgeCohortSizes[a] = 0;
    _isHot.resize(_par->nRunLength);
}


void Community::reset() { // used for r-zero calculations, to reset pop after a single intro
    // reset people
    for (int i=0; i<_nNumPerson; i++) {
        Person* p = _person+i;
        if (p->isWithdrawn(_nDay)) {
            p->getLocation(1)->addPerson(p,1);                        // goes back to work
            p->getLocation(0)->removePerson(p,1);                     // stops staying at home
        }
        p->resetImmunity(); // no past infections, not dead, not vaccinated
    }

    // reset locations
    for (unsigned int i = 0; i < _location.size(); i++ ) _location[i]->clearInfectedMosquitoes();

    for (auto &e: _isHot) e.clear();

    // clear community queues & tallies
    for (unsigned int i = 0; i < _exposedQueue.size(); i++ ) _exposedQueue[i].clear();
    _exposedQueue.clear();

    for (unsigned int i = 0; i < _infectiousMosquitoQueue.size(); i++ ) _infectiousMosquitoQueue[i].clear();
    _infectiousMosquitoQueue.clear();

    for (unsigned int i = 0; i < _exposedMosquitoQueue.size(); i++ ) _exposedMosquitoQueue[i].clear();
    _exposedMosquitoQueue.clear();

    for (unsigned int i = 0; i < _nNumNewlyInfected.size(); i++ ) _nNumNewlyInfected[i].clear();
    _nNumNewlyInfected.clear();

    for (unsigned int i = 0; i < _nNumNewlySymptomatic.size(); i++ ) _nNumNewlySymptomatic[i].clear();
    _nNumNewlySymptomatic.clear();

    for (unsigned int i = 0; i < _nNumVaccinatedCases.size(); i++ ) _nNumVaccinatedCases[i].clear();
    _nNumVaccinatedCases.clear();

    _exposedQueue.resize(MAX_INCUBATION, vector<Person*>(0));
    _infectiousMosquitoQueue.resize(MAX_MOSQUITO_AGE, vector<Mosquito*>(0));
    _exposedMosquitoQueue.resize(MAX_MOSQUITO_AGE, vector<Mosquito*>(0));
    _nNumNewlyInfected.resize(NUM_OF_SEROTYPES, vector<int>(_par->nRunLength + MAX_MOSQUITO_AGE));
    _nNumNewlySymptomatic.resize(NUM_OF_SEROTYPES, vector<int>(_par->nRunLength + MAX_MOSQUITO_AGE));
    _nNumVaccinatedCases.resize(NUM_OF_SEROTYPES, vector<int>(_par->nRunLength + MAX_MOSQUITO_AGE));
}


Community::~Community() {
    if (_person)
        delete [] _person;

    Person::reset_ID_counter();

    for (auto &e: _isHot) e.clear();

    for (unsigned int i = 0; i < _location.size(); i++ ) delete _location[i];
    _location.clear();

    for (unsigned int i = 0; i < _exposedQueue.size(); i++ ) _exposedQueue[i].clear();
    _exposedQueue.clear();

    for (unsigned int i = 0; i < _infectiousMosquitoQueue.size(); i++ ) _infectiousMosquitoQueue[i].clear();
    _infectiousMosquitoQueue.clear();

    for (unsigned int i = 0; i < _exposedMosquitoQueue.size(); i++ ) _exposedMosquitoQueue[i].clear();
    _exposedMosquitoQueue.clear();

    for (unsigned int i = 0; i < _personAgeCohort.size(); i++ ) _personAgeCohort[i].clear();
    _personAgeCohort.clear();
    
    for (unsigned int i = 0; i < _nNumNewlyInfected.size(); i++ ) _nNumNewlyInfected[i].clear();
    _nNumNewlyInfected.clear();

    for (unsigned int i = 0; i < _nNumNewlySymptomatic.size(); i++ ) _nNumNewlySymptomatic[i].clear();
    _nNumNewlySymptomatic.clear();

    for (unsigned int i = 0; i < _nNumVaccinatedCases.size(); i++ ) _nNumVaccinatedCases[i].clear();
    _nNumVaccinatedCases.clear();

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

    int agecounts[NUM_AGE_CLASSES];
    for (int i=0; i<NUM_AGE_CLASSES; i++) agecounts[i] = 0;

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
            assert(age<NUM_AGE_CLASSES);
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
        int line_no = 0;
        while (immiss) {
            line_no++;
            immiss.getline(buffer,500);
            line.clear(); 
            line.str(buffer);
            while (line >> part) parts.push_back(part);

            // 1+ without age, 2+ with age
            if (parts.size() == 1 + NUM_OF_SEROTYPES or parts.size() == 2 + NUM_OF_SEROTYPES) {
                const int id = parts[0];
                Person* person = getPersonByID(id);
                unsigned int offset = parts.size() - NUM_OF_SEROTYPES;
                vector<pair<int,Serotype> > infection_history;
                for (unsigned int f=offset; f<offset+NUM_OF_SEROTYPES; f++) {
                    Serotype s = (Serotype) (f - offset);
                    const int infection_time = parts[f];
                    if (infection_time == 0) {
                        continue; // no infection for this serotype
                    } else if (infection_time<0) {
                        infection_history.push_back(make_pair(infection_time, s));
                    } else {
                        cerr << "ERROR: Found positive-valued infection time in population immunity file:\n\t";
                        cerr << "person " << person->getID() << ", serotype " << s+1 << ", time " << infection_time << "\n\n";
                        cerr << "Infection time should be provided as a negative integer indicated how many days\n";
                        cerr << "before the start of simulation the infection began.";
                        exit(-359);
                    }
                }
                sort(infection_history.begin(), infection_history.end());
                for (auto p: infection_history) person->infect(p.second, p.first + _nDay);
            } else if (parts.size() == 0) {
                continue; // skipping blank line, or line that doesn't start with ints
            } else {
                cerr << "ERROR: Unexpected number of values on one line in population immunity file.\n\t";
                cerr << "line num, line: " << line_no << ", " << buffer << "\n\n";
                cerr << "Expected " << 1+NUM_OF_SEROTYPES << " values (person id followed by infection time for each serotype),\n";
                cerr << "found " << parts.size() << endl;
                exit(-361);
            }
            parts.clear();
        }
        immiss.close();
    }

    // keep track of all age cohorts for aging and mortality
    _personAgeCohort.clear();
    _personAgeCohort.resize(NUM_AGE_CLASSES, vector<Person*>(0));

    for (int i=0; i<_nNumPerson; i++) {
        int age = _person[i].getAge();
        assert(age<NUM_AGE_CLASSES);
        _personAgeCohort[age].push_back(_person + i);
        _nPersonAgeCohortSizes[age]++;
    }

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
                if (person) person->appendToSwapProbabilities(make_pair(id2, prob));
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

            _location.push_back(newLoc);
        }
    }
    iss.close();
    //cerr << _location.size() << " locations" << endl;

    iss.open(networkFilename.c_str());
    if (!iss) {
        cerr << "ERROR: " << networkFilename << " not found." << endl;
        return false;
    }
    int locID1, locID2;
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

    return true;
}

bool Community::loadMosquitoes(string moslocFilename, string mosFilename) {
    if (moslocFilename == "" and mosFilename == "") return true; // nothing to do
    assert(_location.size() > 0); // make sure loadLocations() was already called

    ifstream iss_mosloc(moslocFilename.c_str());
    if (!iss_mosloc) { cerr << "ERROR: " << moslocFilename << " not found." << endl; return false; }

    string buffer;
    int locID, baseMos, infdMos;
    istringstream line(buffer);

    while (iss_mosloc) {
        if (!getline(iss_mosloc, buffer)) break;
        line.clear();
        line.str(buffer);
        if (line >> locID >> baseMos >> infdMos) {
            if (locID >= (signed) _location.size()) {
                cerr << "ERROR: Location ID in mosquito location file greater than largest valid location"
                     << " ID: " << locID << " in file: " << moslocFilename << endl;
                return false;
            }
            Location* loc = _location[locID];
            loc->setBaseMosquitoCapacity(baseMos);
            loc->clearInfectedMosquitoes();
            loc->addInfectedMosquitoes(infdMos);
        }
    }
    iss_mosloc.close();

    ifstream iss_mos(mosFilename.c_str());
    if (!iss_mos) { cerr << "ERROR: " << mosFilename << " not found." << endl; return false; }

    for (unsigned int i = 0; i < _exposedMosquitoQueue.size(); i++ ) _exposedMosquitoQueue[i].clear();
    _exposedMosquitoQueue.clear();
    _exposedMosquitoQueue.resize(MAX_MOSQUITO_AGE, vector<Mosquito*>(0));

    for (unsigned int i = 0; i < _infectiousMosquitoQueue.size(); i++ ) _infectiousMosquitoQueue[i].clear();
    _infectiousMosquitoQueue.clear();
    _infectiousMosquitoQueue.resize(MAX_MOSQUITO_AGE, vector<Mosquito*>(0));

    char queue;
    int sero, idx, ageInfd, ageInfs, ageDead;

    while (iss_mos) {
        if (!getline(iss_mos, buffer)) break;
        line.clear();
        line.str(buffer);
        if (line >> locID >> sero >> queue >> idx >> ageInfd >> ageInfs >> ageDead) {
            if (locID >= (signed) _location.size()) {
                cerr << "ERROR: Location ID in mosquito file greater than largest valid location"
                     << " ID: " << locID << " in file: " << mosFilename << endl;
                return false;
            }
            assert(sero < NUM_OF_SEROTYPES);
            Location* loc = _location[locID];
            Mosquito* m = new Mosquito(loc, (Serotype) sero, ageInfd, ageInfs, ageDead);
            if (queue == 'e') {
                assert(idx < (signed) _exposedMosquitoQueue.size());
                _exposedMosquitoQueue[idx].push_back(m);
            } else if (queue == 'i') {
                assert(idx < (signed) _infectiousMosquitoQueue.size());
                _infectiousMosquitoQueue[idx].push_back(m);
            } else {
                cerr << "ERROR: unknown queue type: " << queue << endl;
                return false;
            }
        }
    }
    iss_mos.close();

    return true;
}


Person* Community::getPersonByID(int id) {
    // This assumes that IDs start at 1, and tries to guess
    // that person with ID id is in position id-1
    // TODO - make that not true (about starting at 1)
    if(id < 1 or id > _nNumPerson) {
        cerr << "ERROR: failed to find person with id " << id << " max: " << _nNumPerson << endl;
        assert(id > 0 and id <= _nNumPerson);
    }

    int i = 0;
    Person* person = NULL;
    if (_person[id-1].getID()==id) {
        i = id-1;
        person = &_person[i];
    } else {
        for (i=0; i<_nNumPerson; i++) {
            if (_person[i].getID()==id) {
                person = &_person[i];
                break;
            }
        }
    }

    if (not person) {
        cerr << "ERROR: failed to find person with id " << id << endl;
        exit(-2001);
    }
    return person;
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
void Community::vaccinate(int time, double f, int age) {
    // This approach to vaccination is somewhat problematic.  Age classes can be vaccinated multiple times,
    // so the probability of an individual being vaccinated becomes 1 - (1 - f)^n, where n is the number
    // of times an age class is specified, either explicitly or implicitly by using a negative value for age
    assert(f>=0.0 and f<=1.0);
    if (age<0) {                                                      // vaccinate everyone
        for (int i=0; i<_nNumPerson; i++) {
            Person p = _person[i];
            if (!p.isVaccinated() && gsl_rng_uniform(RNG)<f) {
                p.vaccinate(time);
            }
        }
    } else {
        // is the specified age a valid index?
        if (age <= NUM_AGE_CLASSES - 1) {
            for (int pnum=0; pnum<_nPersonAgeCohortSizes[age]; pnum++) {
                Person* p = _personAgeCohort[age][pnum];
                assert(p!=NULL);
                if (!p->isVaccinated() && gsl_rng_uniform(RNG)<f) {
                    p->vaccinate(time);
                }
            }
        }
    }
}


void Community::boost(int time, double f) { // re-vaccinate people who have less than threshold immunity left
    assert(f>=0.0 and f<=1.0);
    for (int i=0; i<_nNumPerson; i++) {
        Person* p = _person + i;
        // boost if ~fraction (or more) of immunity has waned
        if (p->isVaccinated() and (p->daysSinceVaccination(time) >= f * _par->vaccineImmunityDuration) ) {
            p->vaccinate(time);
        }
    }
}


// returns number of days mosquito has left to live
int Community::attemptToAddMosquito(Location* p, Serotype serotype, int nInfectedByID) {
    int eip = getExtrinsicIncubation();
    // It doesn't make sense to have an EIP that is greater than the mosquitoes lifespan
    // Truncating also makes vector sizing more straightforward
    eip = eip > MAX_MOSQUITO_AGE ? MAX_MOSQUITO_AGE : eip;
    Mosquito* m = new Mosquito(p, serotype, nInfectedByID, eip);
    int daysleft = m->getAgeDeath() - m->getAgeInfected();
    int daysinfectious = daysleft - eip;
    if (daysinfectious<=0) {
        delete m;
        return daysleft;                                              // dies before infectious
    }

    // add mosquito to latency queue
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


Mosquito* Community::getInfectiousMosquito(int n) {
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


Mosquito* Community::getExposedMosquito(int n) {
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


void Community::moveMosquito(Mosquito* m) {
    double r = gsl_rng_uniform(RNG);
    if (r<_par->fMosquitoMove) {
        if (r<_par->fMosquitoTeleport) {                               // teleport
            int locID = gsl_rng_uniform_int(RNG,_location.size());
            m->updateLocation(_location[locID]);
        } else {                                                            // move to neighbor
            Location* pLoc = m->getLocation();
            double x1 = pLoc->getX();
            double y1 = pLoc->getY();

            int degree = pLoc->getNumNeighbors();
            if (degree == 0) return; // movement isn't possible; no neighbors exist
            int neighbor=0; // neighbor is an index

            if (_par->mosquitoMoveModel == "weighted") {
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
    for (auto &e: _isHot) e.clear();

    // For people of age x, copy immune status from people of age x-1
    for (int age=NUM_AGE_CLASSES-1; age>0; age--) {
        for (int pnum=0; pnum<_nPersonAgeCohortSizes[age]; pnum++) {
            //cerr << "age " << age << ": " << pnum << " of " << _nPersonAgeCohortSizes[age] << endl;
            Person* p = _personAgeCohort[age][pnum];
            assert(p!=NULL);
            if (_uniformSwap == true) {
                // For people of age x, copy immune status from people of age x-1
                // TODO: this may not be safe, if there are age gaps, i.e. people of age N with no one of age N-1
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
        Person* p = _personAgeCohort[0][pnum];
        assert(p!=NULL);
        p->resetImmunity();
    }
    return;
}


void Community::updateWithdrawnStatus() {
    for (int i=0; i<_nNumPerson; i++) {
        Person* p = _person+i;
        if (p->getNumInfections() == 0) continue;
        if (p->getSymptomTime()==_nDay) {                              // started showing symptoms
            _nNumNewlySymptomatic[(int) p->getSerotype()][_nDay]++;
            if (p->isVaccinated()) {
                _nNumVaccinatedCases[(int) p->getSerotype()][_nDay]++;
            }
        }
        if (p->getWithdrawnTime()==_nDay) {                           // started withdrawing
            p->getLocation(0)->addPerson(p,1);                        // stays at home at mid-day
            p->getLocation(1)->removePerson(p,1);                     // does not go to work
        } else if (p->isWithdrawn(_nDay-1) &&
        p->getRecoveryTime()==_nDay) {                                // just stopped withdrawing
            p->getLocation(1)->addPerson(p,1);                        // goes back to work
            p->getLocation(0)->removePerson(p,1);                     // stops staying at home
        }
    }
    return;
}


void Community::flagInfectedLocation(Location* _pLoc, int day) {
    if (day < _par->nRunLength) _isHot[day].insert(_pLoc); 
}


void Community::mosquitoToHumanTransmission() {
    for(unsigned int i=0; i<_infectiousMosquitoQueue.size(); i++) {
        for(unsigned int j=0; j<_infectiousMosquitoQueue[i].size(); j++) {
            Mosquito* m = _infectiousMosquitoQueue[i][j];
            Location* pLoc = m->getLocation();
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
                    Person* p = pLoc->getPerson(idx, timeofday);
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
    for (Location* loc: _isHot[_nDay]) {
        double sumviremic = 0.0;
        double sumnonviremic = 0.0;
        vector<double> sumserotype(NUM_OF_SEROTYPES,0.0);                                    // serotype fractions at location

        // calculate fraction of people who are viremic
        for (int timeofday=0; timeofday<STEPS_PER_DAY; timeofday++) {
            for (int i=loc->getNumPerson(timeofday)-1; i>=0; i--) {
                Person* p = loc->getPerson(i, timeofday);
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
            if (m<0) m=0; // more infected mosquitoes than the base capacity, presumable due to immigration
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
    _isHot[_nDay].clear();
    return;
}


void Community::_advanceTimers() {
    // advance incubation in people
    for (unsigned int i=0; i<_exposedQueue.size()-1; i++) {
        _exposedQueue[i] = _exposedQueue[i+1];
    }
    _exposedQueue.back().clear();

    // delete infected mosquitoes that are dying today
    vector<Mosquito*>::iterator itr; 
    for(itr = _infectiousMosquitoQueue.front().begin(); itr != _infectiousMosquitoQueue.front().end(); ++itr ) {
        delete (*itr);
    }

    // advance age of infectious mosquitoes
    for (unsigned int i=0; i<_infectiousMosquitoQueue.size()-1; i++) {
        _infectiousMosquitoQueue[i] = _infectiousMosquitoQueue[i+1];
#ifndef __INTEL_COMPILER
        _infectiousMosquitoQueue[i].shrink_to_fit();
#endif
    }
    
    _infectiousMosquitoQueue.back().clear();
#ifndef __INTEL_COMPILER
    _infectiousMosquitoQueue.back().shrink_to_fit();
#endif

    assert(_exposedMosquitoQueue.size() > 0);
    // advance incubation period of exposed mosquitoes
    for (unsigned int mnum=0; mnum<_exposedMosquitoQueue[0].size(); mnum++) {
        Mosquito* m = _exposedMosquitoQueue[0][mnum];
        // incubation over: some mosquitoes become infectious
        int daysinfectious = m->getAgeDeath() - m->getAgeInfectious(); // - MOSQUITO_INCUBATION;
        assert((unsigned) daysinfectious < _infectiousMosquitoQueue.size());
        _infectiousMosquitoQueue[daysinfectious].push_back(m);
    }

    for (unsigned int i=0; i<_exposedMosquitoQueue.size()-1; i++) {
        _exposedMosquitoQueue[i] = _exposedMosquitoQueue[i+1];
#ifndef __INTEL_COMPILER
        _exposedMosquitoQueue[i].shrink_to_fit();
#endif
    }
    _exposedMosquitoQueue.back().clear();
#ifndef __INTEL_COMPILER
    _exposedMosquitoQueue.back().shrink_to_fit();
#endif
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

/*{ 
    int w = 0;
    int v = 0;
    int d = 0;
    int ni = 0;
    int i = 0;
    int s = 0;
    int va = 0;
    int in = 0;
    int fs = 0;
    for (int k=0; k<_nNumPerson; k++) {
        Person* p = _person + k;
        if ( p->isWithdrawn(day) )               w++;
        if ( p->isViremic(day) )                 v++;
        if ( p->isDead() )                       d++;

        if ( p->isNewlyInfected(day) )          ni++;
        if ( p->isInfected(day) )                i++; 
        if ( p->isSymptomatic(day) )             s++;

        if ( p->isVaccinated() )                va++;
        if ( p->isInfectable(SEROTYPE_1, day) ) in++;
        if ( p->fullySusceptible() )            fs++;
    }

    cerr << w << " " << v << " " << d << " | " << ni << " " << i << " " << s << " | " << va << " " << in << " " << fs;
    cerr << " || " << _nNumPerson << endl;
}*/

    _nDay = day;
    if ((_nDay+1)%365==0) { swapImmuneStates(); }                     // randomize and advance immune states on
                                                                      // last day of simulator year

    updateWithdrawnStatus();                                          // make people stay home or return to work
    mosquitoToHumanTransmission();                                    // infect people

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
    for (int i=0; i<_nNumPerson; i++) {
        for (int s=0; s<NUM_OF_SEROTYPES; s++) {
            if (_person[i].isSusceptible((Serotype) s))
                counts[s]++;
        }
    }
    return counts;
}
