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
vector< set<Location*, LocPtrComp> > Community::_isHot;
vector<Person*> Community::_peopleByAge;
map<int, set<pair<Person*,Person*> > > Community::_delayedBirthdays;

int mod(int k, int n) { return ((k %= n) < 0) ? k+n : k; } // correct for non-negative n

Community::Community(const Parameters* parameters) :
    _exposedQueue(MAX_INCUBATION, vector<Person*>(0)),
    _infectiousMosquitoQueue(MAX_MOSQUITO_AGE, vector<Mosquito*>(0)),
    // reserving MAX_MOSQUITO_AGE is simpler than figuring out what the maximum
    // possible EIP is when EIP is variable
    _exposedMosquitoQueue(MAX_MOSQUITO_AGE, vector<Mosquito*>(0)),
    _nNumNewlyInfected(NUM_OF_SEROTYPES, vector<int>(parameters->nRunLength + MAX_MOSQUITO_AGE)),
    _nNumNewlySymptomatic(NUM_OF_SEROTYPES, vector<int>(parameters->nRunLength + MAX_MOSQUITO_AGE)),
    _nNumVaccinatedCases(NUM_OF_SEROTYPES, vector<int>(parameters->nRunLength + MAX_MOSQUITO_AGE)),
    _nNumSevereCases(NUM_OF_SEROTYPES, vector<int>(parameters->nRunLength + MAX_MOSQUITO_AGE))
    {
    _par = parameters;
    _nDay = 0;
    _fMosquitoCapacityMultiplier = 1.0;
    _expectedEIP = -1;
    _EIP_emu = -1;
    _fMortality = NULL;
    _bNoSecondaryTransmission = false;
    _uniformSwap = true;
    for (int a = 0; a<NUM_AGE_CLASSES; a++) _nPersonAgeCohortSizes[a] = 0;
    _isHot.resize(_par->nRunLength);
}


void Community::reset() { // used for r-zero calculations, to reset pop after a single intro
    // reset people
    for (Person* p: _people) {
        if (p->isWithdrawn(_nDay)) {
            p->getLocation(WORK_DAY)->addPerson(p,WORK_DAY);                    // goes back to work
            p->getLocation(HOME_MORNING)->removePerson(p,WORK_DAY);             // stops staying at home
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
    if (_people.size() > 0) { for (Person* p: _people) delete p; }

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
    string buffer;
    int agecounts[NUM_AGE_CLASSES];
    for (int i=0; i<NUM_AGE_CLASSES; i++) agecounts[i] = 0;

    istringstream line;
    int id, age, house, work;
    string hh_serial;
    int pernum;
    int sex; // per IPUMS, expecting 1 for male, 2 for female
    while ( getline(iss,buffer) ) {
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
        if (line >> id >> house >> age >> sex >> hh_serial >> pernum >> work) {
            if (house < 1000) {
                Person* p = new Person();
                _people.push_back(p);
                p->setAge(age);
                p->setSex((SexType) sex);
                p->setHomeID(house);
                p->setLocation(_location[house], HOME_MORNING);
                p->setLocation(_location[work], WORK_DAY);
                p->setLocation(_location[house], HOME_NIGHT);
                _location[house]->addPerson(p, HOME_MORNING);
                _location[work]->addPerson(p, WORK_DAY);
                _location[house]->addPerson(p, HOME_NIGHT);
                assert(age<NUM_AGE_CLASSES);
                agecounts[age]++;
                p->infect(SEROTYPE_1,-10);
            }
        }
    }
    iss.close();

    _peopleByAge = _people;
    sort(_peopleByAge.begin(), _peopleByAge.end(), PerPtrComp());

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
        while ( getline(immiss,buffer) ) {
            line_no++;
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

    for (Person* p: _people) {
        int age = p->getAge();
        assert(age<NUM_AGE_CLASSES);
        _personAgeCohort[age].push_back(p);
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

        while ( getline(iss, buffer) ) {
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
    string locTypeStr;
    double locX, locY;
    istringstream line(buffer);

    while (iss) {
        iss.getline(buffer,500);
        line.clear();
        line.str(buffer);
        if (line >> locID >> locTypeStr >> locX >> locY) {
            if (locID != (signed) _location.size()) {
                cerr << "WARNING: Location ID's must be sequential integers" << endl;
                return false;
            }
            const LocationType locType = (locTypeStr == "house") ? HOME : (locTypeStr == "work") ? WORK : (locTypeStr == "school") ? SCHOOL : NUM_OF_LOCATION_TYPES;
            if (locType == NUM_OF_LOCATION_TYPES) {
                cerr << "ERROR: Parsed unknown location type: " << locTypeStr << " from location file: " << locationFilename << endl;
                return false;
            }
            Location* newLoc = new Location();
            newLoc->setID(locID);
            newLoc->setX(locX);
            newLoc->setY(locY);
            newLoc->setType(locType);

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
    int locID, baseMos;
    istringstream line(buffer);

    while (iss_mosloc) {
        if (!getline(iss_mosloc, buffer)) break;
        line.clear();
        line.str(buffer);
        if (line >> locID >> baseMos) { // there may be an infected_mosquito_ct field, but that is handled
                                        // when we call the mos constructor while parsing the mosquito file
            if (locID >= (signed) _location.size()) {
                cerr << "ERROR: Location ID in mosquito location file greater than largest valid location"
                     << " ID: " << locID << " in file: " << moslocFilename << endl;
                return false;
            }
            Location* loc = _location[locID];
            loc->setBaseMosquitoCapacity(baseMos);
            loc->clearInfectedMosquitoes();
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
            RestoreMosquitoPars restorePars(loc, (Serotype) sero, ageInfd, ageInfs, ageDead);
            Mosquito* m = new Mosquito(&restorePars);
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
    if(id < 1 or id > getNumPeople()) {
        cerr << "ERROR: failed to find person with id " << id << " max: " << getNumPeople() << endl;
        assert(id > 0 and id <= getNumPeople());
    }

    int i = 0;
    Person* person = NULL;
    if (_people[id-1]->getID()==id) {
        i = id-1;
        person = _people[i];
    } else {
        for (Person* p: _people) {
            if (p->getID()==id) {
                person = p;
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


void Community::vaccinate(VaccinationEvent ve) {
    // This approach to vaccination is somewhat problematic.  Age classes can be vaccinated multiple times,
    // so the probability of an individual being vaccinated becomes 1 - (1 - ve.coverage)^n, where n is the number
    // of times an age class is specified, either explicitly or implicitly by using a negative value for age

    // Valid coverage and age?
    assert(ve.coverage >= 0.0 and ve.coverage <= 1.0);
    assert(ve.age <= (signed) _personAgeCohort.size());

    for (Person* p: _personAgeCohort[ve.age]) {
        assert(p != NULL);
        if (!p->isVaccinated() && gsl_rng_uniform(RNG) < ve.coverage) p->vaccinate(ve.simDay);
    }
}


void Community::boost(int time, int interval, int maxDoses) { // re-vaccinate people who were last vaccinated interval days ago
    for (Person* p: _people) {
        if (p->isVaccinated()) {
            const int timeSinceLastVaccination = p->daysSinceVaccination(time);
            if (timeSinceLastVaccination >= interval and p->getNumVaccinations() < maxDoses) {
                p->vaccinate(time);
            }
        }
    }
}


// returns number of days mosquito has left to live
void Community::attemptToAddMosquito(Location* p, Serotype serotype, int nInfectedByID) {
    int eip = (int) (getEIP() + 0.5);

    // It doesn't make sense to have an EIP that is greater than the mosquitoes lifespan
    // Truncating also makes vector sizing more straightforward
    eip = eip > MAX_MOSQUITO_AGE ? MAX_MOSQUITO_AGE : eip;
    Mosquito* m = new Mosquito(p, serotype, nInfectedByID, eip);
    int daysleft = m->getAgeDeath() - m->getAgeInfected();
    int daysinfectious = daysleft - eip;
    if (daysinfectious<=0) {
        // dies before infectious
        delete m;
    } else {
        if (eip == 0) {
            // infectious immediately -- unlikely, but supported
            // we don't push onto index 0, because we're at the end of the day already;
            // this mosquito would be destroyed before being allowed to transmit
            _infectiousMosquitoQueue[daysinfectious].push_back(m);
        } else {
            // more typically, add mosquito to latency queue
            _exposedMosquitoQueue[eip].push_back(m);
        }
    }
    return;
}


void Community::mosquitoFilter(vector<Mosquito*>& mosquitoes, const double survival_prob) {
    if (survival_prob >= 1.0) return;
    const unsigned int nmos = mosquitoes.size();
    if (nmos == 0) return;
    gsl_ran_shuffle(RNG, mosquitoes.data(), nmos, sizeof(Mosquito*));
    const int survivors = gsl_ran_binomial(RNG, survival_prob, nmos);
    for (unsigned int m = survivors; m<mosquitoes.size(); ++m) delete mosquitoes[m];
    mosquitoes.resize(survivors);
}


void Community::applyMosquitoMultiplier(double current) {
    const double prev = getMosquitoMultiplier();
    setMosquitoMultiplier(current);
    if (current < prev) {
        const double survival_prob = current/prev;
        for (unsigned int day = 0; day < _exposedMosquitoQueue.size(); ++day) mosquitoFilter(_exposedMosquitoQueue[day], survival_prob);
        for (unsigned int day = 0; day < _infectiousMosquitoQueue.size(); ++day) mosquitoFilter(_infectiousMosquitoQueue[day], survival_prob);
    }
}


void Community::applyVectorControl() {
    for (Location* loc: _location) loc->updateVectorControlQueue(_nDay); // make sure proper VC is active for tomorrow -- must be at end

    vector<Mosquito*> swap;
    for (unsigned int day = 0; day < _exposedMosquitoQueue.size(); ++day) {
        for (Mosquito* m: _exposedMosquitoQueue[day]) {
            const float vc_rho = m->getLocation()->getCurrentVectorControlDailyMortality(_nDay);
            if (vc_rho > 0 and gsl_rng_uniform(RNG) < vc_rho) {
                delete m;
            } else {
                swap.push_back(m);
            }
        }
        _exposedMosquitoQueue[day] = swap;
        swap.clear();
    }

    for (unsigned int day = 0; day < _infectiousMosquitoQueue.size(); ++day) {
        for (Mosquito* m: _infectiousMosquitoQueue[day]) {
            const float vc_rho = m->getLocation()->getCurrentVectorControlDailyMortality(_nDay);
            if (vc_rho > 0 and gsl_rng_uniform(RNG) < vc_rho) {
                delete m;
            } else {
                swap.push_back(m);
            }
        }
        _infectiousMosquitoQueue[day] = swap;
        swap.clear();
    }
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
        if (r<_par->fMosquitoTeleport) {                // teleport
            int locID = gsl_rng_uniform_int(RNG,_location.size());
            m->updateLocation(_location[locID]);
        } else {                                        // move to neighbor
            Location* pLoc = m->getLocation();
            double x1 = pLoc->getX();
            double y1 = pLoc->getY();

            int degree = pLoc->getNumNeighbors();
            if (degree == 0) return;                    // movement isn't possible; no neighbors exist
            int neighbor=0;                             // neighbor is an index

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
                neighbor = degree-1;                    // neighbor is (still) an index
                int idx;
                for ( idx = 0; idx < degree - 1; idx++ ) {
                    weights[idx] /= sum_weights;        // normalize prob
                    if ( r2 < weights[idx] ) {
                        break;
                    } else {
                        r2 -= weights[idx];
                    }
                }
                neighbor = idx;
            } else {                                    // Alternatively, ignore distances when choosing destination
                if (degree>0) {
                    neighbor = gsl_rng_uniform_int(RNG,pLoc->getNumNeighbors());
                }
            }

            m->updateLocation(pLoc->getNeighbor(neighbor));
        }
    }
}


void Community::_processBirthday(Person* p) {
    Person* donor;
    if (p->getAge() == 0) {
        //p->resetImmunity();
        donor = nullptr;
    } else if (_uniformSwap == true) {
        // For people of age x, copy immune status from people of age x-1
        // TODO: this may not be safe, if there are age gaps, i.e. people of age N with no one of age N-1
        const int donor_age = p->getAge() - 1;
        int r = gsl_rng_uniform_int(RNG,_nPersonAgeCohortSizes[donor_age]);
        donor = _personAgeCohort[donor_age][r];
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
        donor = getPersonByID(id);
    }
    if (_par->delayBirthdayIfInfected) {
        _swapIfNeitherInfected(p, donor);
    } else {
        if (donor) {
            p->copyImmunity(donor);
        } else {
            p->resetImmunity();
        }
        // update map of locations with infectious people -- not necessary if birthdays are delayed until after infection resolves
        // if this is a historical infection, person may have neg values for
        // infectious dates that we don't need to deal with--they're in the past
        if (p->isInfected(_nDay)) {
            const int start_date = p->getInfectiousTime() > 0 ? p->getInfectiousTime() : 0;
            for (int d = start_date; d < p->getRecoveryTime(); d++) {
                for (int t=0; t<(int) NUM_OF_TIME_PERIODS; t++) {
                    flagInfectedLocation(p->getLocation((TimePeriod) t), d);
                }
            }
        }
    }
}


void Community::_swapIfNeitherInfected(Person* p, Person* donor) {
    int process_date = p->isInfected(_nDay) ? p->getRecoveryTime() : _nDay; // delay birthday if p is infected
    process_date = donor and donor->isInfected(_nDay) ? std::max(process_date, donor->getRecoveryTime()) : process_date; // and/or delay if donor is infected
    if (process_date == _nDay) {
        if (donor) {
            p->copyImmunity(donor);
        } else {
            p->resetImmunity();
        }
    } else {
        if (_delayedBirthdays.count(process_date) == 0) _delayedBirthdays[process_date] = set<pair< Person*, Person*> >();
        _delayedBirthdays[process_date].insert(make_pair(p, donor));
    }
}


void Community::_processDelayedBirthdays() {
    if (_delayedBirthdays.count(_nDay) > 0) {
        for (pair<Person*, Person*> recip_donor: _delayedBirthdays[_nDay]) {
            Person* p     = recip_donor.first;
            Person* donor = recip_donor.second;
            _swapIfNeitherInfected(p, donor);
        }
        _delayedBirthdays.erase(_nDay);
    }
}

// 1.) Process delayed birthdays every day -- need to handle runs where normal birthdays are handled with interval > 1
// 2.) Test delayed birthdays with having birthdays once per year on Julian day 100
void Community::swapImmuneStates() {
    const int julian = _nDay % 365;
    const int bday_ivl = _par->birthdayInterval;
    // For people of age x, copy immune status from people of age x-1

    for (int y = 0; y < ceil(_peopleByAge.size()/365.0); ++y) {
        const int minval = (365*y + julian + 1 <= bday_ivl ) ? 0 : 365*y + julian - bday_ivl + 1;
        const int maxval = (365*y + julian >= (signed) _peopleByAge.size()) ? _peopleByAge.size() - 1 : 365*y + julian;
        for (int pidx = minval; pidx <= maxval; ++pidx) {
            Person* p = _peopleByAge[pidx];
            assert(p!=NULL);
            _processBirthday(p);
        }
    }
    return;
}


void Community::updateDiseaseStatus() {
    for (Person* p: _people) {
        if (p->getNumNaturalInfections() == 0) continue;
        if (p->getSymptomTime()==_nDay) {                              // started showing symptoms today
            _nNumNewlySymptomatic[(int) p->getSerotype()][_nDay]++;
            if (p->isVaccinated()) {
                _nNumVaccinatedCases[(int) p->getSerotype()][_nDay]++;
            }
            if (p->hasSevereDisease(_nDay)) {                          // symptoms will be severe at onset
                _nNumSevereCases[(int) p->getSerotype()][_nDay]++;     // if they're going to be severe
            }
        }
        if (p->getWithdrawnTime()==_nDay) {                            // started withdrawing
            p->getLocation(HOME_MORNING)->addPerson(p,WORK_DAY);       // stays at home at mid-day
            p->getLocation(WORK_DAY)->removePerson(p,WORK_DAY);        // does not go to work
        } else if (p->isWithdrawn(_nDay-1) and
        p->getRecoveryTime()==_nDay) {                                 // just stopped withdrawing
            p->getLocation(WORK_DAY)->addPerson(p,WORK_DAY);           // goes back to work
            p->getLocation(HOME_MORNING)->removePerson(p,WORK_DAY);    // stops staying at home
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
                double exposuretime[(int) NUM_OF_TIME_PERIODS];
                double totalExposureTime = 0;
                for (int t=0; t<(int) NUM_OF_TIME_PERIODS; t++) {
                    exposuretime[t] = pLoc->getNumPerson((TimePeriod) t) * DAILY_BITING_PDF[t];
                    totalExposureTime += exposuretime[t];
                }
                if ( totalExposureTime > 0 ) {
                    double r = gsl_rng_uniform(RNG) * totalExposureTime;
                    int timeofday;
                    for (timeofday=0; timeofday<(int) NUM_OF_TIME_PERIODS - 1; timeofday++) {
                        if (r<exposuretime[timeofday]) {
                            // bite at this time of day
                            break;
                        }
                        r -= exposuretime[timeofday];
                    }
                    int idx = floor(r*pLoc->getNumPerson((TimePeriod) timeofday)/exposuretime[timeofday]);
                    Person* p = pLoc->getPerson(idx, (TimePeriod) timeofday);
                    Serotype serotype = m->getSerotype();
                    if (p->infect(m->getID(), serotype, _nDay, pLoc->getID())) {
                        _nNumNewlyInfected[(int) serotype][_nDay]++;
                        if (_bNoSecondaryTransmission) {
                            p->kill();                       // kill secondary cases so they do not transmit
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
        for (int timeofday=0; timeofday<(int) NUM_OF_TIME_PERIODS; timeofday++) {
            for (int i=loc->getNumPerson((TimePeriod) timeofday)-1; i>=0; i--) {
                Person* p = loc->getPerson(i, (TimePeriod) timeofday);
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
            int m = int(loc->getBaseMosquitoCapacity() * (1.0-loc->getCurrentVectorControlEfficacy(_nDay)) * getMosquitoMultiplier() + 0.5);  // number of mosquitoes
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
                attemptToAddMosquito(loc, (Serotype) serotype, locid);
            }
        }
    }
    _isHot[_nDay].clear();

    int total_mos_sum = 0;
    int inf_mos_sum = 0;
    for (int i = 1; i<=1000; ++i) {
        Location* loc = _location[i];
        total_mos_sum += int(loc->getBaseMosquitoCapacity() * (1.0-loc->getCurrentVectorControlEfficacy(_nDay)) * getMosquitoMultiplier() + 0.5);  // number of mosquitoes
        inf_mos_sum += loc->getCurrentInfectedMosquitoes();
    }
    cerr << "day, total, infected mos: " <<_nDay << ", " << total_mos_sum/1000.0 << ", " << inf_mos_sum/1000.0 << endl;
    
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
    _nDay = day;
    //if ((_nDay+1)%365==0) { swapImmuneStates(1.0); }                     // randomize and advance immune states on
    //_processDelayedBirthdays();
    //if ((_nDay+1) % _par->birthdayInterval == 0) { swapImmuneStates(); }     // randomize and advance some immune states
    if (_par->vectorControlEvents.size() > 0) applyVectorControl();   // also advances vector control status to next day
    /*{
    const Location* _l = _location[4];
    cerr << "day " << day
         << ", loc " << _l->getID()
         << ", " << (int) _l->getType()
         << ", base: " << _l->getBaseMosquitoCapacity()
         << ", start: " << _l->getCurrentVectorControl()->start_day
         << ", cur: " << (int)(_l->getBaseMosquitoCapacity() * (1.0-_l->getCurrentVectorControlEfficacy(_nDay)) * getMosquitoMultiplier() + 0.5)
         << ", cur inf: " << _l->getCurrentInfectedMosquitoes()
    //     << ", set eff: " << _l->getVectorControlEfficacy()
         << ", cur eff: " << _l->getCurrentVectorControlEfficacy(_nDay)
         << ", mortality: " << _l->getCurrentVectorControlDailyMortality(_nDay) << endl;
    }*/
                                                                      // last day of simulator year

    updateDiseaseStatus();                                            // make people stay home or return to work
    mosquitoToHumanTransmission();                                    // infect people

    humanToMosquitoTransmission();                                    // infect mosquitoes in each location
    _advanceTimers();                                                 // advance H&M incubation periods and M ages
    _modelMosquitoMovement();                                         // probabilistic movement of mosquitos
    return;
}


// getNumInfected - counts number of infected residents
int Community::getNumInfected(int day) {
    int count=0;
    for (Person* p: _people) { if (p->isInfected(day)) count++; }
    return count;
}


// getNumSymptomatic - counts number of symptomatic residents
int Community::getNumSymptomatic(int day) {
    int count=0;
    for (Person* p: _people) { if (p->isSymptomatic(day)) count++; }
    return count;
}

// getNumSusceptible - counts number of susceptible residents
vector<int> Community::getNumSusceptible() {
    vector<int> counts(NUM_OF_SEROTYPES, 0);
    for (Person* p: _people) {
        for (int s=0; s<NUM_OF_SEROTYPES; s++) {
            if (p->isSusceptible((Serotype) s)) counts[s]++;
        }
    }
    return counts;
}
