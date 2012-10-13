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

using namespace std;

const double Community::SYMPTOMATICBYAGE[MAXPERSONAGE] = {
    0.05189621,0.05189621,0.05189621,0.05189621,0.05189621,
    0.1017964,0.1017964,0.1017964,0.1017964,0.1017964,
    0.2774451,0.2774451,0.2774451,0.2774451,0.2774451,
    0.4870259,0.4870259,0.4870259,0.4870259,0.4870259,
    0.4870259,0.4870259,0.4870259,0.4870259,0.4870259,
    0.8522954,0.8522954,0.8522954,0.8522954,0.8522954,
    0.8522954,0.8522954,0.8522954,0.8522954,0.8522954,
    0.9600798,0.9600798,0.9600798,0.9600798,0.9600798,
    0.9600798,0.9600798,0.9600798,0.9600798,0.9600798,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
};

//2005 Thai mortality data by age from Porapakkham 2010
//double thaimortality[Community::MAXPERSONAGE] = {0.0157,0.0009,0.0009,0.0009,0.0009,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005,0.0007,0.0007,0.0007,0.0007,0.0007,0.0009,0.0009,0.0009,0.0009,0.0009,0.0016,0.0016,0.0016,0.0016,0.0016,0.002,0.002,0.002,0.002,0.002,0.0022,0.0022,0.0022,0.0022,0.0022,0.0028,0.0028,0.0028,0.0028,0.0028,0.0038,0.0038,0.0038,0.0038,0.0038,0.005,0.005,0.005,0.005,0.005,0.0077,0.0077,0.0077,0.0077,0.0077,0.012,0.012,0.012,0.012,0.012,0.0185,0.0185,0.0185,0.0185,0.0185,0.0287,0.0287,0.0287,0.0287,0.0287,0.0457,0.0457,0.0457,0.0457,0.0457,0.0767,0.0767,0.0767,0.0767,0.0767,0.1434,0.1434,0.1434,0.1434,0.1434};

// Community
Community::Community() {
    _nDay = 0;
    _nNumPerson = 0;
    _nNumLocation = 0;
    _person = NULL;
    _location = NULL;
    _fBetaPM = 1.0;
    _fBetaMP = 0.05;
    _fMosquitoMoveProb = 0.3;
    _fMosquitoTeleportProb = 0.02;
    _fDailyBitingPDF[0] = 0.08;
    _fDailyBitingPDF[1] = 0.76;
    _fDailyBitingPDF[2] = 0.16;
    _nMaxInfectionParity = NUM_OF_SEROTYPES;
    _exposedQueue = new Person **[MAXINCUBATION];
    _nExposedQueueCapacity = 100;
    _fMosquitoCapacityMultiplier = 1.0;
    for (int i=0; i<MAXINCUBATION; i++) {
        _exposedQueue[i] = new Person *[_nExposedQueueCapacity];
        _exposedQueue[i][0] = NULL;
    }
    _infectiousMosquitoQueue = new Mosquito**[MAXMOSQUITOAGE-MOSQUITOINCUBATION];
    _nMosquitoQueueCapacity = 100;
    for (int i=0; i<MAXMOSQUITOAGE-MOSQUITOINCUBATION; i++) {
        _infectiousMosquitoQueue[i] = new Mosquito *[_nMosquitoQueueCapacity];
        _infectiousMosquitoQueue[i][0] = NULL;
    }
    _exposedMosquitoQueue = new Mosquito**[MOSQUITOINCUBATION];
    for (int i=0; i<MOSQUITOINCUBATION; i++) {
        _exposedMosquitoQueue[i] = new Mosquito *[_nMosquitoQueueCapacity];
        _exposedMosquitoQueue[i][0] = NULL;
    }
    _personAgeCohort = NULL;
    _fMortality = NULL;
    //  _fMortality = thaimortality;
    _numLocationMosquitoCreated = NULL;
    _bNoSecondaryTransmission = false;
    _fPrimarySymptomaticScaling[0] = _fPrimarySymptomaticScaling[1] = _fPrimarySymptomaticScaling[2] = _fPrimarySymptomaticScaling[3] = 1.0;
    _fSecondarySymptomaticScaling[0] = _fSecondarySymptomaticScaling[1] = _fSecondarySymptomaticScaling[2] = _fSecondarySymptomaticScaling[3] = 1.0;
    for (int i=0; i<MAXRUNTIME; i++)
        for (int j=0; j<NUM_OF_SEROTYPES; j++)
            _nNumNewlyInfected[j][i] = _nNumNewlySymptomatic[j][i] = 0;
}


Community::~Community() {
    if (_person)
        delete [] _person;
    if (_location)
        delete [] _location;
    if (_exposedQueue) {
        for (int i=0; i<MAXINCUBATION; i++)
            if (_exposedQueue[i])
                delete [] _exposedQueue[i];
        delete [] _exposedQueue;
    }
    if (_infectiousMosquitoQueue) {
        for (int i=0; i<MAXMOSQUITOAGE-MOSQUITOINCUBATION; i++)
            if (_infectiousMosquitoQueue[i])
                delete [] _infectiousMosquitoQueue[i];
        delete [] _infectiousMosquitoQueue;
    }
    if (_exposedMosquitoQueue) {
        for (int i=0; i<MOSQUITOINCUBATION; i++)
            if (_exposedMosquitoQueue[i])
                delete [] _exposedMosquitoQueue[i];
        delete [] _exposedMosquitoQueue;
    }
    if (_personAgeCohort) {
        for (int i=0; i<MAXPERSONAGE; i++)
            if (_personAgeCohort[i])
                delete [] _personAgeCohort[i];
        delete [] _personAgeCohort;
    }
    if (_numLocationMosquitoCreated) {
        for (int i=0; i<_nNumLocation; i++)
            if (_numLocationMosquitoCreated[i])
                delete [] _numLocationMosquitoCreated[i];
        delete []_numLocationMosquitoCreated;
    }
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
    int numlines = 0;
    while (iss) {
        char temp[500];
        iss.getline(temp,500);
        numlines++;
    }
    iss.close();
    iss.open(_szPopulationFilename.c_str());
    _nNumPerson=0;
    _nMaxPerson = numlines*4;
    _person = new Person[_nMaxPerson];
    int agecounts[MAXPERSONAGE];
    for (int i=0; i<MAXPERSONAGE; i++)
        agecounts[i] = 0;
    while (iss) {
        char temp[500];
        iss.getline(temp,500);
        istringstream line(temp);
        int id, age, house, work;//, imm1, imm2, imm3, imm4;
        //imm1 = imm2 = imm3 = imm4 = 0;
        int temp1, temp2;
        string gender;
        //    if (line >> id >> age >> gender >> house >> work >> imm1 >> imm2 >> imm3 >> imm4) {
                                                                      // >> imm1 >> imm2 >> imm3 >> imm4) {
        if (line >> id >> house >> age >> gender >> temp1 >> temp2 >> work) {
            _person[_nNumPerson].setAge(age);
            _person[_nNumPerson].setHomeID(house);
            _person[_nNumPerson].setLocation(_location+house, 0);
            _person[_nNumPerson].setLocation(_location+work, 1);
            _person[_nNumPerson].setLocation(_location+house, 2);
            _location[house].addPerson(_person+_nNumPerson, 0);
            _location[work].addPerson(_person+_nNumPerson, 1);
            _location[house].addPerson(_person+_nNumPerson, 2);
            assert(age<MAXPERSONAGE);
            agecounts[age]++;
            /*if (imm1>0)
                _person[_nNumPerson].setImmunity(1);
            if (imm2>0)
                _person[_nNumPerson].setImmunity(2);
            if (imm3>0)
                _person[_nNumPerson].setImmunity(3);
            if (imm4>0)
                _person[_nNumPerson].setImmunity(4);*/
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
        while (immiss) {
            char temp[500];
            immiss.getline(temp,500);
            istringstream line(temp);
            int id, imm1, imm2, imm3, imm4;
            string gender;
            if (line >> id >> imm1 >> imm2 >> imm3 >> imm4) {
                for (int i=id-1; i<=id; i++) {
                    if (_person[i].getID()==id) {
                        if (imm1>0)
                            _person[i].setImmunity(SEROTYPE_1);
                        if (imm2>0)
                            _person[i].setImmunity(SEROTYPE_2);
                        if (imm3>0)
                            _person[i].setImmunity(SEROTYPE_3);
                        if (imm4>0)
                            _person[i].setImmunity(SEROTYPE_4);
                        break;
                    }
                }
            }
        }
        immiss.close();
    }

    // keep track of all age cohorts for aging and mortality
    _nPersonAgeCohortMaxSize = 0;
    for (int i=0; i<MAXPERSONAGE; i++)
        if (_nPersonAgeCohortMaxSize<agecounts[i])
            _nPersonAgeCohortMaxSize=agecounts[i];
    _nPersonAgeCohortMaxSize *= 1.25;
    _personAgeCohort = new Person **[MAXPERSONAGE];
    for (int i=0; i<MAXPERSONAGE; i++) {
        _nPersonAgeCohortSizes[i] = 0;
        _personAgeCohort[i] = new Person *[_nPersonAgeCohortMaxSize];
    }
    for (int i=0; i<_nNumPerson; i++) {
        int age = _person[i].getAge();
        assert(age<MAXPERSONAGE);
        _personAgeCohort[age][_nPersonAgeCohortSizes[age]++]=_person+i;
    }
    for (int i=0; i<MAXPERSONAGE; i++)
        _nOriginalPersonAgeCohortSizes[i] = _nPersonAgeCohortSizes[i];
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
    _nNumLocation=0;
    while (iss) {
        char temp[500];
        iss.getline(temp,500);
        istringstream line(temp);
        int temp1;
        if (line >> temp1) {
            if (temp1>_nNumLocation)
                _nNumLocation = temp1;
        }
    }
    iss.close();
    _nNumLocation+=1;                                                 // the ID of the highest location + 1
    cerr << _nNumLocation << " locations" << endl;
    _location = new Location[_nNumLocation];

    _numLocationMosquitoCreated = new int *[_nNumLocation];
    for (int i=0; i<_nNumLocation; i++)
        _numLocationMosquitoCreated[i] = NULL;

    iss.open(_szNetworkFilename.c_str());
    if (!iss) {
        cerr << "ERROR: " << _szNetworkFilename << " not found." << endl;
        return false;
    }
    while (iss) {
        char temp[500];
        iss.getline(temp,500);
        istringstream line(temp);
        int temp1,temp2;
        if (line >> temp1 >> temp2) {
            //      cerr << temp1 << " , " << temp2 << endl;
            _location[temp1].addNeighbor(_location+temp2);            // should check for ID
            _location[temp2].addNeighbor(_location+temp1);
        }
        else {
            //      cerr << "header" << endl;
        }
    }
    iss.close();
    for (int i=0; i<_nNumLocation; i++)
        if (_location[i].getNumNeighbors()<=0)
            _location[i].setUndefined();
    return true;
}


// infect - infects person id
bool Community::infect(gsl_rng *rng, int id, Serotype serotype, int day) {
    for (int i=0; i<_nNumPerson; i++)
    if (_person[i].getID()==id) {
        bool result =  _person[i].infect(rng, -1, serotype, day, 0,
            _fPrimarySymptomaticScaling[(int) serotype],
            _fSecondarySymptomaticScaling[(int) serotype],
            _nMaxInfectionParity);
        if (result)
            _nNumNewlyInfected[(int) serotype][_nDay]++;
        return result;
        //      cerr << "inf " << i << " at "  << _person[i].getLocation(0)->getID() << endl;
    }
    return false;
}


// vaccinate - vaccinate fraction f of the population
// if age>=0, then vaccinate only those who are "age" years old
void Community::vaccinate(gsl_rng *rng, double f, int age) {
    if (f<=0.0)
        return;
    if (age<0) {                                                      // vaccinate everyone
        for (int i=0; i<_nNumPerson; i++)
            if (gsl_rng_uniform(rng)<f)
                _person[i].vaccinate(rng);
    }
    else {
        for (int pnum=0; pnum<_nPersonAgeCohortSizes[age]; pnum++) {
            Person *p = _personAgeCohort[age][pnum];
            assert(p!=NULL);
            if (!p->isVaccinated() && gsl_rng_uniform(rng)<f)
                p->vaccinate(rng);
        }
    }
}


void Community::setVES(double f) {
    Person::setVES(f);
}


void Community::setVESs(double f1,double f2,double f3,double f4) {
    Person::setVESs(f1,f2,f3,f4);
}


void Community::setVEI(double f) {
    Person::setVEI(f);
}


void Community::setVEP(double f) {
    Person::setVEP(f);
}


// expandMosquitoQueues - increases the capacity of the queues
// that keep track of the mosquito populations
void Community::expandMosquitoQueues() {
    for (int j=0; j<MAXMOSQUITOAGE-MOSQUITOINCUBATION; j++) {
        Mosquito **temp = new Mosquito *[_nMosquitoQueueCapacity*2];
        for (int k=0; k<_nMosquitoQueueCapacity; k++)
            temp[k] = _infectiousMosquitoQueue[j][k];
        delete [] _infectiousMosquitoQueue[j];
        _infectiousMosquitoQueue[j] = temp;
    }
    for (int j=0; j<MOSQUITOINCUBATION; j++) {
        Mosquito **temp = new Mosquito *[_nMosquitoQueueCapacity*2];
        for (int k=0; k<_nMosquitoQueueCapacity; k++)
            temp[k] = _exposedMosquitoQueue[j][k];
        delete [] _exposedMosquitoQueue[j];
        _exposedMosquitoQueue[j] = temp;
    }
    _nMosquitoQueueCapacity*=2;
}


// returns number of days mosquito has left to live
int Community::addMosquito(gsl_rng *rng, Location *p, Serotype serotype, int nInfectedByID) {
    assert(serotype>=0 && serotype<4);
    Mosquito *m = new Mosquito(rng, p, serotype, nInfectedByID);
    int daysleft = m->getAgeDeath()-m->getAgeInfected();
    int daysinfectious = daysleft-MOSQUITOINCUBATION;
    if (daysinfectious<=0) {
        delete m;
        return daysleft;                                              // dies before infectious
    }
    // add mosquito to latency queue
    for (int i=0; i<_nMosquitoQueueCapacity; i++) {
        if (_exposedMosquitoQueue[MOSQUITOINCUBATION-1][i]==NULL) {
            if (i+2>_nMosquitoQueueCapacity)
                expandMosquitoQueues();
            // add mosquito to latency queue
            _exposedMosquitoQueue[MOSQUITOINCUBATION-1][i] = m;
            _exposedMosquitoQueue[MOSQUITOINCUBATION-1][i+1]=NULL;
            break;
        }
    }
    return daysleft;
}


int Community::getNumInfectiousMosquitoes() {
    int count = 0;
    for (int mosAge=0; mosAge<MAXMOSQUITOAGE-MOSQUITOINCUBATION-1; mosAge++)
        for (Mosquito **m=_infectiousMosquitoQueue[mosAge]; *m!=NULL; m++)
            count++;
    return count;
}


int Community::getNumExposedMosquitoes() {
    int count = 0;
    for (int i=0; i<MOSQUITOINCUBATION; i++)
        for (Mosquito **m=_exposedMosquitoQueue[i]; *m!=NULL; m++)
            count++;
    return count;
}


Mosquito *Community::getInfectiousMosquito(int n) {
    int count = 0;
    for (int mosAge=0; mosAge<MAXMOSQUITOAGE-MOSQUITOINCUBATION-1; mosAge++)
    for (Mosquito **m=_infectiousMosquitoQueue[mosAge]; *m!=NULL; m++) {
        if (count==n)
            return *m;
        count++;
    }
    return NULL;
}


Mosquito *Community::getExposedMosquito(int n) {
    int count = 0;
    for (int i=0; i<MOSQUITOINCUBATION; i++)
    for (Mosquito **m=_exposedMosquitoQueue[i]; *m!=NULL; m++) {
        if (count==n)
            return *m;
        count++;
    }
    return NULL;
}


// expandExposedQueues - increases the capacity of the queues
// that keep track of exposed humans
void Community::expandExposedQueues() {
    for (int i=0; i<MAXINCUBATION; i++) {
        Person **temp = new Person *[_nExposedQueueCapacity*2];
        for (int j=0; j<_nExposedQueueCapacity; j++)
            temp[j] = _exposedQueue[i][j];
        delete [] _exposedQueue[i];
        _exposedQueue[i] = temp;
    }
    _nExposedQueueCapacity *= 2;
}


void Community::moveMosquito(Mosquito *m, gsl_rng *rng) {
    double r = gsl_rng_uniform(rng);
    if (r<_fMosquitoMoveProb) {
        if (r<_fMosquitoTeleportProb) {                               // teleport
            int location;
            do {
                location = gsl_rng_uniform_int(rng,_nNumLocation);
            } while (_location[location].getUndefined());
            m->setLocation(_location+location);
        }                                                             // move to neighbor
        else {
            Location *pLoc = m->getLocation();
            if (pLoc->getNumNeighbors()>0) {
                int n = gsl_rng_uniform_int(rng,pLoc->getNumNeighbors());
                m->setLocation(pLoc->getNeighbor(n));
            }
        }
    }
}


void Community::tick(gsl_rng *rng) {
    assert(_nDay<MAXRUNTIME);

    //  if (false) {
    if (_nDay%365==364) {
        //    cerr << "swap!" << endl;
        // big annual swap!
        // For people of age x, copy immune status from people of age x-1
        for (int age=MAXPERSONAGE-1; age>0; age--) {
            //      cerr << "age " << age << ", " << _nPersonAgeCohortSizes[age] << " people" << endl;
            int age1 = age-1;
            for (int pnum=0; pnum<_nPersonAgeCohortSizes[age]; pnum++) {
                Person *p = _personAgeCohort[age][pnum];
                assert(p!=NULL);
                int r = gsl_rng_uniform_int(rng,_nPersonAgeCohortSizes[age1]);
                p->copyImmunity(_personAgeCohort[age1][r]);
            }
        }
        // For people of age 0, reset immunity
        for (int pnum=0; pnum<_nPersonAgeCohortSizes[0]; pnum++) {
            Person *p = _personAgeCohort[0][pnum];
            assert(p!=NULL);
            p->resetImmunity();
        }
    }

    // check for withdrawn people
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

    // infect people
    for (int mosAge=0; mosAge<MAXMOSQUITOAGE-MOSQUITOINCUBATION-1; mosAge++) {
        for (Mosquito **m=_infectiousMosquitoQueue[mosAge]; *m!=NULL; m++) {
            Location *pLoc = (**m).getLocation();
            if (gsl_rng_uniform(rng)<_fBetaMP) {                      // infectious mosquito bites
                // take sum of people in the location, weighting by time of day
                double exposuretime[STEPSPERDAY];
                for (int timeofday=0; timeofday<STEPSPERDAY; timeofday++)
                    exposuretime[timeofday] = pLoc->getNumPerson(timeofday) * _fDailyBitingPDF[timeofday];
                double r = gsl_rng_uniform(rng)*(exposuretime[0]+exposuretime[1]+exposuretime[2]);
                for (int timeofday=0; timeofday<STEPSPERDAY; timeofday++) {
                    if (r<exposuretime[timeofday]) {
                        // bite at this time of day
                        int id = floor(r*pLoc->getNumPerson(timeofday)/exposuretime[timeofday]);
                        Person *p = pLoc->getPerson(id, timeofday);
                        //	    cerr << " infect person " << id << "/" << pLoc->getNumPerson(timeofday) << "  " << p->getID() << " with " << (**m).getSerotype() << endl;
                        Serotype serotype = (**m).getSerotype();
                        if (p->infect(rng, (**m).getID(), serotype, _nDay, pLoc->getID(),
                            _fPrimarySymptomaticScaling[(int) serotype],
                            _fSecondarySymptomaticScaling[(int) serotype],
                        _nMaxInfectionParity)) {
                            _nNumNewlyInfected[(int) p->getSerotype()][_nDay]++;
                            if (_bNoSecondaryTransmission) {
                                p->kill(_nDay);                       // kill secondary cases so they do not transmit
                            }
                            else {
                                // NOTE: We are storing the location ID of infection, not person ID!!!
                                //	      cerr << "Person " << p->getID() << " infected" << endl;
                                // add to queue
                                Person **q = _exposedQueue[p->getInfectiousTime()-_nDay];
                                int count=0;
                                while (*q) {
                                    q++;
                                    count++;
                                }
                                *q = p;
                                q++;
                                *q = NULL;
                                if (count+2>=_nExposedQueueCapacity)
                                    expandExposedQueues();
                            }
                        }
                        break;
                    }
                    r -= exposuretime[timeofday];
                }
            }
        }
    }

    // maybe add the occasional random infection of a person to reflect sporadic travel

    // infect mosquitoes in each location
    for (int loc=0; loc<_nNumLocation; loc++) {
        if (!_location[loc].getUndefined()) {
            double sumviremic = 0.0;
            double sumnonviremic = 0.0;
            double sumserotype[5];                                    // serotype fractions at location
            sumserotype[0] = sumserotype[1] = sumserotype[2] = sumserotype[3] = sumserotype[4] = 0.0;

            // calculate fraction of people who are viremic
            for (int timeofday=0; timeofday<STEPSPERDAY; timeofday++) {
                for (int i=_location[loc].getNumPerson(timeofday)-1; i>=0; i--) {
                    Person *p = _location[loc].getPerson(i, timeofday);
                    if (p->isViremic(_nDay)) {
                        double vaceffect = (p->isVaccinated()?(1.0-Person::getVEI()):1.0);
                        int serotype = (int) p->getSerotype();
                        //assert(serotype<4);
                        if (vaceffect==1.0) {
                            sumviremic += _fDailyBitingPDF[timeofday];
                            sumserotype[serotype] += _fDailyBitingPDF[timeofday];
                        } else {
                            sumviremic += _fDailyBitingPDF[timeofday]*vaceffect;
                            sumserotype[serotype] += _fDailyBitingPDF[timeofday]*vaceffect;
                            // a vaccinated person is treated like a fraction of an infectious person and a fraction of a non-infectious person
                            sumnonviremic += _fDailyBitingPDF[timeofday]*(1.0-vaceffect);
                        }
                    } else {
                        sumnonviremic += _fDailyBitingPDF[timeofday];
                    }
                }
            }

            if (sumviremic>0.0) {
                for (int i=0; i<5; i++)
                    sumserotype[i] /= sumviremic;
                int m;                                                // number of susceptible mosquitoes
                int locid = _location[loc].getID();                   // location ID
                m = int(_location[loc].getBaseMosquitoCapacity()*_fMosquitoCapacityMultiplier+0.5);
                if (_numLocationMosquitoCreated[locid])
                    m -= _numLocationMosquitoCreated[locid][_nDay];
                if (m<0)
                    m=0;
                                                                      // how many susceptible mosquitoes bite viremic hosts in this location?
                int numbites = gsl_ran_binomial(rng, _fBetaPM*sumviremic/(sumviremic+sumnonviremic), m);
                while (numbites-->0) {
                    int serotype;                                     // which serotype infects mosquito
                    if (sumserotype[0]==1.0) {
                        serotype = 0;
                    } else {
                        double r = gsl_rng_uniform(rng);
                        for (serotype=0; serotype<5 && r>sumserotype[serotype]; serotype++)
                            r -= sumserotype[serotype];
                    }
                    /*if (serotype==0) {
                        cerr << "serotypes: ";
                        for (int i=0; i<5; i++)
                            cerr << sumserotype[i] << ",";
                        cerr << endl;
                    }*/
                    assert (serotype >= 0 && serotype<4);
                    int daysleft = addMosquito(rng, _location+loc, (Serotype) serotype, locid);
                    if (!_numLocationMosquitoCreated[locid]) {
                        _numLocationMosquitoCreated[locid] = new int[MAXRUNTIME];
                        for (int j=0; j<MAXRUNTIME; j++)
                            _numLocationMosquitoCreated[locid][j] = 0;//_location[loc].getMosquitoCapacity();
                    }
                    assert(_nDay+daysleft<MAXRUNTIME);
                    for (int j=0; j<daysleft; j++)
                        _numLocationMosquitoCreated[locid][_nDay+j]++;
                }
            }
            //	  cerr << "Mosquito infected by Person " << _person[i].getID() << " at " << _person[i].getLocation(timeofday)->getID() << "." << endl;
        }
    }

    // advance incubation in people
    Person **p = _exposedQueue[0];
    for (int i=0; i<MAXINCUBATION-1; i++)
        _exposedQueue[i] = _exposedQueue[i+1];
    _exposedQueue[MAXINCUBATION-1] = p;
    _exposedQueue[MAXINCUBATION-1][0] = NULL;

    // advance age of infectious mosquitoes
    Mosquito **pm = _infectiousMosquitoQueue[0];
    pm[0] = NULL;
    for (int i=0; i<MAXMOSQUITOAGE-MOSQUITOINCUBATION-1; i++)
        _infectiousMosquitoQueue[i] = _infectiousMosquitoQueue[i+1];
    _infectiousMosquitoQueue[MAXMOSQUITOAGE-MOSQUITOINCUBATION-1] = pm;

    // advance incubation period of exposed mosquitoes
    for (int mnum=0; ; mnum++) {
        // incubation over: some mosquitoes become infectious
        Mosquito *m = _exposedMosquitoQueue[0][mnum];
        if (m==NULL)
            break;
        int daysinfectious = m->getAgeDeath()-m->getAgeInfected()-MOSQUITOINCUBATION;
        for (int i=0; i<_nMosquitoQueueCapacity; i++) {
            if (_infectiousMosquitoQueue[daysinfectious][i]==NULL) {
                if (i+2>_nMosquitoQueueCapacity)
                    expandMosquitoQueues();
                // move mosquito to infectious queue
                _infectiousMosquitoQueue[daysinfectious][i] = m;
                _infectiousMosquitoQueue[daysinfectious][i+1]=NULL;
                break;
            }
        }
    }
    pm = _exposedMosquitoQueue[0];
    pm[0] = NULL;
    for (int i=0; i<MOSQUITOINCUBATION-1; i++)
        _exposedMosquitoQueue[i] = _exposedMosquitoQueue[i+1];
    _exposedMosquitoQueue[MOSQUITOINCUBATION-1] = pm;

    // move mosquitoes
    for (int mosAge=0; mosAge<MAXMOSQUITOAGE-MOSQUITOINCUBATION-1; mosAge++)
        for (Mosquito **m=_infectiousMosquitoQueue[mosAge]; *m!=NULL; m++)
            moveMosquito(*m, rng);
    for (int mosdays=0; mosdays<MOSQUITOINCUBATION; mosdays++)
        for (Mosquito **m=_exposedMosquitoQueue[mosdays]; *m!=NULL; m++)
            moveMosquito(*m, rng);

    _nDay++;
    //  if (_nDay%365==1) {
    if (false) {
        for (int blah=0; blah<50; blah++) {
            cerr << "Day " <<  _nDay << endl;

            int newborncount = 0;                                     // number of people who die/are born this year

            // kill off oldest age bracket
            // use _personAgeCohort[MAXPERSONAGE-1] to hold newborns
            for (int i=0; i<_nPersonAgeCohortSizes[MAXPERSONAGE-1]; i++) {
                assert (_nNumPerson+1<_nMaxPerson);
                Person *p = _personAgeCohort[MAXPERSONAGE-1][i];
                p->kill(_nDay);
                // copy attributes to newborn
                Person *newborn = _person + _nNumPerson;
                newborn->setHomeID(p->getHomeID());
                newborn->setWorkID(p->getWorkID());                   // infants work at old place?
                newborn->setAge(0);
                _personAgeCohort[MAXPERSONAGE-1][i] = _person + _nNumPerson;
                newborncount++;
                assert (newborncount<_nPersonAgeCohortMaxSize);
                _nNumPerson++;
            }

            for (int age=1; age<MAXPERSONAGE-1; age++) {
                int numkill = gsl_ran_binomial(rng, _fMortality[age], _nPersonAgeCohortSizes[age]);
                //      cerr << _nDay << "," << age << ", kill " << numkill << "/" << _nPersonAgeCohortSizes[age] << endl;
                for (int i=0; i<numkill; i++) {
                    int r = gsl_rng_uniform_int(rng,_nPersonAgeCohortSizes[age]);
                    Person *p = _personAgeCohort[age][r];
                    assert(p!=NULL);
                    //	cerr << "select " << r << ", age " << age << ", size=" << _nPersonAgeCohortSizes[age] << ", numperson=" << _nNumPerson << endl;
                    if (!p->isDead()) {
                        //	    cerr << ". kill " << p->getID() << ", age=" << p->getAge();
                        p->kill(_nDay);
                        // copy attributes to newborn
                        Person *newborn = _person + _nNumPerson;
                        newborn->setHomeID(p->getHomeID());
                        newborn->setWorkID(p->getWorkID());           // infants work at old place?
                        newborn->setAge(0);
                        _personAgeCohort[MAXPERSONAGE-1][newborncount] = newborn;
                        newborncount++;
                        assert (newborncount<_nPersonAgeCohortMaxSize);
                        _nNumPerson++;
                        assert(_nNumPerson<_nMaxPerson);
                        _nPersonAgeCohortSizes[age]--;
                        // move the last person in the list to this spot
                        if (_nPersonAgeCohortSizes[age]>0)
                            _personAgeCohort[age][r] =
                                _personAgeCohort[age][_nPersonAgeCohortSizes[age]];
                        else
                            _personAgeCohort[age][r] = NULL;
                    }
                }
            }

            /* 
            // kill younger people based on population pyramid
            // do this in 5-year age groups
            int diff[MAXPERSONAGE];
            for (int i=1; i<MAXPERSONAGE; i++)
              diff[i-1] = _nPersonAgeCohortSizes[i-1] - _nOriginalPersonAgeCohortSizes[i];
            for (int group=0; group<(MAXPERSONAGE-2)/5; group++) {
              int groupdiff = diff[group*5] + diff[group*5+1] + diff[group*5+2] + diff[group*5+3] + diff[group*5+4];
              //      cerr << "Group " << group << " diff= " << groupdiff << endl;
              if (groupdiff>0) {
            // kill excess people in this age group
            int groupsize = 0;
            for (int i=1; i<5; i++)
              groupsize += _nPersonAgeCohortSizes[group*5-1+i];
            if (group>0)
              groupsize += _nPersonAgeCohortSizes[group*5-1];
            assert(groupsize>0);
            while (groupdiff>0) {
              int r = gsl_rng_uniform_int(rng,groupsize);
              int age=group*5-1;
              if (age<1)
                age=1;
              while(_nPersonAgeCohortSizes[age]<=r) {
                r-=_nPersonAgeCohortSizes[age];
                age++;
                assert(age<MAXPERSONAGE);
              }
              assert(r<_nPersonAgeCohortSizes[age]);
              Person *p = _personAgeCohort[age][r];
              assert(p!=NULL);
              cerr << "select " << r << ", age " << age << ", size=" << _nPersonAgeCohortSizes[age] << ", numperson=" << _nNumPerson << endl;
              if (!p->isDead()) {
                //	    cerr << ". kill " << p->getID() << ", age=" << p->getAge();
                p->kill();
                // copy attributes to newborn
                Person *newborn = _person + _nNumPerson;
                newborn->setHomeID(p->getHomeID());
                newborn->setWorkID(p->getWorkID()); // infants work at old place?
                newborn->setAge(0);
                _personAgeCohort[MAXPERSONAGE-1][r] = _person + _nNumPerson;
                newborncount++;
                assert (newborncount<_nPersonAgeCohortMaxSize);
                _nNumPerson++;
                assert(_nNumPerson<_nMaxPerson);
                groupdiff--;
                _nPersonAgeCohortSizes[age]--;
                // move the last person in the list to this spot
                _personAgeCohort[age][r] =
                  _personAgeCohort[age][--_nPersonAgeCohortSizes[age]];
                groupsize--;
              }
            }
              }
            }
            */
            // age population 1 year
            Person **temp = _personAgeCohort[MAXPERSONAGE-1];
            for (int i=MAXPERSONAGE-1; i>0; i--) {
                _personAgeCohort[i] = _personAgeCohort[i-1];
                _nPersonAgeCohortSizes[i] = _nPersonAgeCohortSizes[i-1];
            }
            _nPersonAgeCohortSizes[0] = newborncount;
            cerr << "newborncount " << newborncount << endl;
            _personAgeCohort[0] = temp;
            for (int i=0; i<_nNumPerson; i++)
                if (!_person[i].isDead())
                    _person[i].setAge(_person[i].getAge()+1);

            // children age out and start working?

            for (int i=0; i<MAXPERSONAGE; i++) {
                cout << _nDay << " " << i << " " << _nPersonAgeCohortSizes[i] << endl;
            }
            _nDay++;
        }
    }
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
int Community::getNumSusceptible(Serotype serotype) {
    int count=0;
    for (int i=0; i<_nNumPerson; i++)
        if (_person[i].isSusceptible(serotype))
            count++;
    return count;
}
