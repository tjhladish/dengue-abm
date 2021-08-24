// Person.h
// A single individual.

#ifndef __PERSON_H
#define __PERSON_H
#include <bitset>
#include <vector>
#include <climits>
#include "Parameters.h"
#include "Location.h"

class Location;
class Mosquito;

class Infection {
    friend class Person;
    Infection() {
        infectedBy     = nullptr;
        infectedLoc    = nullptr;
        infectedTime   = INT_MIN;
        infectiousTime = INT_MIN;
        symptomTime    = INT_MIN;
        recoveryTime   = INT_MIN;

        withdrawnTime  = INT_MAX;
        _serotype      = NULL_SEROTYPE;
        severeDisease  = false;
    };

    Infection(const Serotype sero) {
        infectedBy     = nullptr;
        infectedLoc    = nullptr;
        infectedTime   = INT_MIN;
        infectiousTime = INT_MIN;
        symptomTime    = INT_MIN;
        recoveryTime   = INT_MIN;

        withdrawnTime  = INT_MAX;
        _serotype      = sero;
        severeDisease  = false;
    };

    Infection(const Infection* o) {
        infectedBy     = o->infectedBy;
        infectedLoc    = o->infectedLoc;
        infectedTime   = o->infectedTime;
        infectiousTime = o->infectiousTime;
        symptomTime    = o->symptomTime;
        recoveryTime   = o->recoveryTime;
        withdrawnTime  = o->withdrawnTime;
        _serotype      = o->_serotype;
        severeDisease  = o->severeDisease;
    }

    Mosquito* infectedBy;                           // who infected this person
    Location* infectedLoc;                          // where infected?
    int infectedTime;                               // when infected?
    int infectiousTime;                             // when infectious period starts
    int symptomTime;                                // when symptoms start
    int recoveryTime;                               // when recovered?
    int withdrawnTime;                              // when person withdraws to home
    Serotype _serotype;
    bool severeDisease;

  public:

    bool isLocallyAcquired() const { return (bool) infectedBy; }
    int getInfectedTime() const { return infectedTime; }
    bool isSymptomatic() const { return symptomTime > infectedTime; }
    bool isSevere()      const { return severeDisease; }
    Serotype serotype()  const { return _serotype; }

    Location* getInfectedLoc() const { return infectedLoc; }
};

class Person {
    public:
        Person();
        ~Person();
        inline int getID() const { return _nID; }
        int getAge() const { return _nAge; }
        void setAge(int n) { _nAge = n; }
        SexType getSex() const { return _sex; }
        void setSex(SexType sex) { _sex = sex; }
        int getLifespan() const { return _nLifespan; }
        void setLifespan(int n) { _nLifespan = n; }
        Location* getHomeLoc() const { return getLocation(HOME_MORNING); }
        Location* getDayLoc() const { return getLocation(WORK_DAY); }
        void setImmunity(Serotype serotype) { _nImmunity[(int) serotype] = 1; }
        const std::bitset<NUM_OF_SEROTYPES> getImmunityBitset() const { return _nImmunity; }
        const std::string getImmunityString() const { return _nImmunity.to_string(); }
        void copyImmunity(const Person *p);
        void resetImmunity();
        void appendToSwapProbabilities(std::pair<int, double> p) { _swap_probabilities.push_back(p); }
        std::vector<std::pair<int, double> > getSwapProbabilities() const { return _swap_probabilities; }

        bool isSusceptible(Serotype serotype) const;                  // is susceptible to serotype (and is alive)
        bool isCrossProtected(int time) const;
        bool isVaccineProtected(Serotype serotype, int time) const;

        inline Location* getLocation(TimePeriod timeofday) const { return _pLocation[(int) timeofday]; }
        inline void setLocation(Location* p, TimePeriod timeofday) { _pLocation[(int) timeofday] = p; }

        inline Mosquito* getInfectedBy(int infectionsago=0) const  { return getInfection(infectionsago)->infectedBy; }
        inline Location* getInfectedLoc(int infectionsago=0) const { return getInfection(infectionsago)->infectedLoc; }
        inline int getInfectedTime(int infectionsago=0) const      { return getInfection(infectionsago)->infectedTime; }
        inline int getInfectiousTime(int infectionsago=0) const    { return getInfection(infectionsago)->infectiousTime; }
        inline int getSymptomTime(int infectionsago=0) const       { return getInfection(infectionsago)->symptomTime; }
        inline int getRecoveryTime(int infectionsago=0) const      { return getInfection(infectionsago)->recoveryTime; }
        inline int getWithdrawnTime(int infectionsago=0) const     { return getInfection(infectionsago)->withdrawnTime; }
        inline Serotype getSerotype(int infectionsago=0) const     { return getInfection(infectionsago)->serotype(); }
        const Infection* getInfection(int infectionsago=0) const   { return infectionHistory[getNumNaturalInfections() - 1 - infectionsago]; }

        inline void setRecoveryTime(int time, int infectionsago=0) { infectionHistory[getNumNaturalInfections() - 1 - infectionsago]->recoveryTime = time; }
        bool isWithdrawn(int time) const;                             // at home sick?
        inline int getNumNaturalInfections() const { return infectionHistory.size(); }
        inline int getEffectiveNumInfections() const {
            int order = getNumNaturalInfections();
            if (isVaccinated()) {
                if (_par->whoDiseaseOutcome == INC_NUM_INFECTIONS or (order == 0 and _par->whoDiseaseOutcome == INC_INFECTIONS_NAIVE)) {
                    order += 1;
                }
            };
            return order;
        }

        int getNumVaccinations() const { return vaccineHistory.size(); }
        const std::vector<int>& getVaccinationHistory() const { return vaccineHistory; }
        const std::vector<Infection*>& getInfectionHistory() const { return infectionHistory; }
        int daysSinceVaccination(int time) const { assert( vaccineHistory.size() > 0); return time - vaccineHistory.back(); } // isVaccinated() should be called first
        double vaccineProtection(const Serotype serotype, const int time) const;

        bool infect(Mosquito* mos, int time, Location* loc, Serotype serotype);
        inline bool infect(int time, Serotype serotype) {return infect(nullptr, time, nullptr, serotype);}
        bool isViremic(int time) const;

        void kill();
        bool isDead() const { return _bDead; }
        bool naturalDeath(int t);                                     // die of old age check?

        bool isNewlyInfected(int time) const;                         // became infected today?
        bool isInfected(int time) const;                              // is currently infected
        bool isSymptomatic(int time) const;                           // has symptoms
        bool hasSevereDisease(int time) const;                        // used for estimating hospitalizations
        bool isVaccinated() const {                                   // has been vaccinated
            return _bVaccinated;
        }
        bool isInfectable(Serotype serotype, int time) const;         // more complicated than isSusceptible
        double remainingEfficacy(const int time) const;

        bool fullySusceptible() const;
                                                                      // does this person's immune state permit vaccination?
                                                                      // NB: inaccurate test results are possible
        bool isSeroEligible(VaccineSeroConstraint vsc, double falsePos, double falseNeg) const;
        bool vaccinate(int time);                                     // vaccinate this person
        static void setPar(const Parameters* par) { _par = par; }

        static const double _fIncubationDistribution[MAX_INCUBATION];

        Infection& initializeNewInfection(Serotype serotype);
        Infection& initializeNewInfection(Mosquito* mos, int time, Location* loc, Serotype serotype);

        static void reset_ID_counter() { _nNextID = 1; }

    protected:
        int _nID;                                                     // unique identifier
        int _nHomeID;                                                 // family membership
        int _nWorkID;                                                 // ID of location of work
        Location *_pLocation[(int) NUM_OF_TIME_PERIODS];              // where this person is at morning, day, and evening
        int _nAge;                                                    // age in years
        SexType _sex;                                                 // sex (gender)
        int _nLifespan;                                               // lifespan in years
        bool _bDead;                                                  // is dead
        std::bitset<NUM_OF_SEROTYPES> _nImmunity;                     // bitmask of serotype infection
        bool _bVaccinated;                                            // has been vaccinated
        bool _bNaiveVaccineProtection; // if vaccinated, do we use the naive or non-naive VE_S?

        std::vector<std::pair<int,double> > _swap_probabilities;      // list of the nearest people one year younger, with distances
        std::vector<Infection*> infectionHistory;
        std::vector<int> vaccineHistory;
        void clearInfectionHistory();

        static const Parameters* _par;
        static int _nNextID;                                          // unique ID to assign to the next Person allocated
};
#endif
