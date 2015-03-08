#include <cstdlib>
#include <cstring>
#include <climits>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Parameters.h"
#include "Person.h"
#include "Mosquito.h"
#include "Location.h"
#include "Community.h"
#include "Utility.h"

using namespace dengue::standard;
using namespace dengue::util;

class Date {
  public:
    Date():_offset(0),_simulation_day(0) {};
    Date(const Parameters* par):_offset(par->startDayOfYear-1),_simulation_day(0) {};

    int offset()             const { return _offset; }
    inline int day()         const { return _simulation_day; }                // [0, ...]
    int julianDay()          const { return ((day() + offset()) % 365) + 1; } // [1, 365]
    int dayOfMonth()         const { return julianMonth() == 1 ?              // [1, {29,30,31}]
                                            julianDay() :
                                            julianDay() - END_DAY_OF_MONTH[julianMonth()-2]; }
    int week()               const { return (int) day()/7 ; }                 // [0, ...]
    int julianWeek()         const { return (int) ((julianDay()-1)/7) + 1; }  // [1, 53]
    int month()              const { return _month_ct; }                      // [0, ...]
    int julianMonth()        const {                                          // [1, 12]
        vector<int>::const_iterator it;
        // find first month that hasn't ended (hint: it's this month)
        // julianDay()-1 because this isn't upper_or_equal_bound, which would be convenient
        it = upper_bound(END_DAY_OF_MONTH.begin(), END_DAY_OF_MONTH.end(), julianDay()-1);
        return it - END_DAY_OF_MONTH.begin() + 1; // +1 b/c [1, 12], not [0, 11]
    }
    string monthName()       const { return MONTH_NAMES[julianMonth()-1]; }
    int year()               const { return (int) (day()/365); }

    bool endOfWeek()         const { return (day()+1) % 7 == 0; }
    bool endOfMonth()        const {
        vector<int>::const_iterator it;
        // find out if today corresponds to a month-end
        it = find(END_DAY_OF_MONTH.begin(), END_DAY_OF_MONTH.end(), julianDay());
        return it != END_DAY_OF_MONTH.end();
    }
    bool startOfYear()       const { return day() % 365 == 0; }     // is it beginning of 365 day period
    bool endOfYear()         const { return (day()+1) % 365 == 0; } // is it end of 365 day period
    bool startOfJulianYear() const { return julianDay() == 1; }     // is it Jan 1
    bool endOfJulianYear()   const { return julianDay() == 365; }   // is it Dec 31

    void increment() {
        if(endOfMonth()) _month_ct++;
        _simulation_day++;
    }

    void print() {
        cerr << day() << "\t" << julianDay() << "\t" << year()
             << "\t(" << monthName() << " " << dayOfMonth() << ")\t" << month() << "\t" << julianMonth();
        if (endOfWeek()) cerr << " EoW";
        if (endOfMonth()) cerr << " EoM";
        if (startOfYear()) cerr << " SoY";
        if (endOfYear()) cerr << " EoY";
        if (endOfJulianYear()) cerr << " EoJY";
        cerr << endl;
    }

  private:
    const int _offset;
    int _simulation_day;
    int _month_ct;
};

const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);

// Predeclare local functions
Community* build_community(const Parameters* par);
void seed_epidemic(const Parameters* par, Community* community);
vector<int> simulate_epidemic(const Parameters* par, Community* community, const int process_id=0);
void write_immunity_file(const Parameters* par, const Community* community, const int int_label, string filename="");
void write_immunity_by_age_file(const Parameters* par, const Community* community, const int year, string filename="");
void write_output(const Parameters* par, Community* community, vector<int> initial_susceptibles);
void write_daily_buffer( vector<string>& buffer, const int process_id);

Community* build_community(const Parameters* par) {
    Community* community = new Community(par);
    Person::setPar(par);

    if (!community->loadLocations(par->locationFilename, par->networkFilename)) {
        cerr << "ERROR: Could not load locations" << endl;
        exit(-1);
    }
    if (!community->loadPopulation(par->populationFilename, par->immunityFilename, par->swapProbFilename)) {
        cerr << "ERROR: Could not load population" << endl;
        exit(-1);
    }

    if (!par->abcVerbose) {
        cerr << community->getNumPerson() << " people" << endl;
    }

    if (!par->bSecondaryTransmission) {
        community->setNoSecondaryTransmission();
    }

    if (par->fPreVaccinateFraction>0.0) {
        community->vaccinate(par->fPreVaccinateFraction);
    }

    if (par->nSizePrevaccinateAge>0) {
        for (int j=0; j<par->nSizePrevaccinateAge; j++) {
            for (int k=par->nPrevaccinateAgeMin[j]; k<=par->nPrevaccinateAgeMax[j]; k++) {
                community->vaccinate(par->fPrevaccinateAgeFraction[j],k);
            }
        }
    }
    
    //for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
    //    par->nNumInitialSusceptible[serotype] = community->getNumSusceptible((Serotype) serotype);
    //}
    return community;
}


void seed_epidemic(const Parameters* par, Community* community) {
    // epidemic may be seeded with initial exposure OR initial infection 
    bool attempt_initial_infection = true;
    for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
        // Normal usage, to simulate epidemic
        if (par->nInitialExposed[serotype] > 0) {
            attempt_initial_infection = false;
            for (int i=0; i<par->nInitialExposed[serotype]; i++)
                community->infect(gsl_rng_uniform_int(RNG, community->getNumPerson()), (Serotype) serotype,0);
        }
    }
    if (attempt_initial_infection) {
        for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
            // Useful for estimating R0
            if(par->nInitialInfected[serotype] > 0) {
                int count = community->getNumInfected(0);

                // must infect nInitialInfected persons -- this bit is mysterious
                while (community->getNumInfected(0) < count + par->nInitialInfected[serotype]) {
                    community->infect(gsl_rng_uniform_int(RNG, community->getNumPerson()), (Serotype) serotype,0);
                }
            }
        }
    }
    return;
}


void write_yearly_people_file(const Parameters* par, const Community* community, int time) {
    ofstream yearlyPeopleOutputFile;
    ostringstream ssFilename;
    ssFilename << par->yearlyPeopleOutputFilename << ((int)(time/365)) << ".csv";
    cerr << "outputing yearly people information to " << ssFilename.str() << endl;
    yearlyPeopleOutputFile.open(ssFilename.str().c_str());
    if(yearlyPeopleOutputFile.fail()) {
        cerr << "ERROR: People file '" << par->yearlyPeopleOutputFilename << "' cannot be open for writing." << endl;
        exit(-1);
    }
    yearlyPeopleOutputFile << "pid,serotype,infectiontime,symptomtime,withdrawtime,recoverytime,immdenv1,immdenv2,immdenv3,immdenv4" << endl;
    for (int i=0; i<community->getNumPerson(); i++) {
        Person *p = community->getPerson(i);
        for (int j=p->getNumInfections()-1; j>=0; j--) {
            yearlyPeopleOutputFile << p->getID() << "," 
                << 1 + (int) p->getSerotype(j) << "," 
                << p->getInfectedTime(j) << "," 
                << p->getSymptomTime(j) << "," 
                << p->getWithdrawnTime(j) << "," 
                << p->getRecoveryTime(j) << ",";
            for (int s = 0; s < NUM_OF_SEROTYPES - 1; ++s) {
                yearlyPeopleOutputFile << (p->isSusceptible((Serotype) s)?0:1) << ","; 
            }
                yearlyPeopleOutputFile << (p->isSusceptible((Serotype) (NUM_OF_SEROTYPES - 1))?0:1) << endl;
        }
    }
    yearlyPeopleOutputFile.close();
    return;
}

void initialize_seasonality(const Parameters* par, Community* community, int& nextMosquitoMultiplierIndex, int& nextEIPindex, Date& date) {
    const int mosquitoMultiplierTotalDuration = par->getMosquitoMultiplierTotalDuration();
    int currentDayOfYearOffset = 0;
    if (mosquitoMultiplierTotalDuration > 0) {
        nextMosquitoMultiplierIndex = 0;
        while (currentDayOfYearOffset <= date.offset()) {
            currentDayOfYearOffset += par->mosquitoMultipliers[nextMosquitoMultiplierIndex].duration;
            if (currentDayOfYearOffset > date.offset()) {
                const double mm = par->mosquitoMultipliers[nextMosquitoMultiplierIndex].value;
                community->setMosquitoMultiplier(mm);
            }
            nextMosquitoMultiplierIndex = (nextMosquitoMultiplierIndex+1)%par->mosquitoMultipliers.size();
        }
    }

    const int EIPtotalDuration = par->getEIPtotalDuration();
    currentDayOfYearOffset = 0;
    if (EIPtotalDuration > 0) {
        nextEIPindex = 0;
        while (currentDayOfYearOffset <= date.offset()) {
            currentDayOfYearOffset += par->extrinsicIncubationPeriods[nextEIPindex].duration;
            if (currentDayOfYearOffset > date.offset()) {
                const double eip = par->extrinsicIncubationPeriods[nextEIPindex].value;
                community->setExtrinsicIncubation(eip);
            }
            nextEIPindex = (nextEIPindex+1)%par->extrinsicIncubationPeriods.size();
        }
    }
}


void periodic_output(const Parameters* par, const Community* community, map<string, vector<int> >& periodic_incidence, const Date& date, const int process_id, vector<int>& epi_sizes) {
        // local transmission          = total                          - introductions
        periodic_incidence["daily"][1] = periodic_incidence["daily"][2] - periodic_incidence["daily"][0];
        if (par->dailyOutput) {
            cerr << "day: " << date.day() << " "; for (auto v: periodic_incidence["daily"]) cerr << v << " ";
            cerr << community->getExtrinsicIncubation() << " " << community->getMosquitoMultiplier()*par->nDefaultMosquitoCapacity << endl;
        }

        if (par->weeklyOutput) {
            for (unsigned int i = 0; i < periodic_incidence["daily"].size(); ++i) periodic_incidence["weekly"][i] += periodic_incidence["daily"][i];
            if (date.endOfWeek()) {
                cerr << "week: " << date.week() << " "; for (auto v: periodic_incidence["weekly"]) cerr << v << " "; cerr << endl;
                periodic_incidence["weekly"] = {0,0,0};
            }
        }

        if (par->monthlyOutput) {
            for (unsigned int i = 0; i < periodic_incidence["daily"].size(); ++i) periodic_incidence["monthly"][i] += periodic_incidence["daily"][i];
            if (date.endOfMonth()) {
                cerr << "month: " << date.julianMonth() << " "; for (auto v: periodic_incidence["monthly"]) cerr << v << " "; cerr << endl;
                periodic_incidence["monthly"] = {0,0,0};
            }
        }

        // handle several things that happen yearly

        for (unsigned int i = 0; i < periodic_incidence["daily"].size(); ++i) periodic_incidence["yearly"][i] += periodic_incidence["daily"][i];

        if (date.endOfYear()) {
            if (par->abcVerbose) cout << hex << process_id << dec << " T: " << date.day() << " annual: " << periodic_incidence["yearly"][2] << endl;

            epi_sizes.push_back(periodic_incidence["yearly"][2]);

            if (par->yearlyPeopleOutputFilename.length() > 0) write_yearly_people_file(par, community, date.day());

            if (par->yearlyOutput) {
                cerr << "year: " << date.year() + 1 << " "; 
                for (auto v: periodic_incidence["yearly"]) cerr << v << " "; cerr << endl;
            }
            periodic_incidence["yearly"] = {0,0,0};
        }

        periodic_incidence["daily"] = {0,0,0};
}

vector<int> simulate_epidemic(const Parameters* par, Community* community, const int process_id, vector<int> &serotested_ids, double &seropos_87) {
    assert(serotested_ids.size() > 0);
    vector<int> epi_sizes;
    Date date(par);
    int nextMosquitoMultiplierIndex = 0;
    const int mosquitoMultiplierTotalDuration = par->getMosquitoMultiplierTotalDuration();
    int nextEIPindex = 0;
    const int EIPtotalDuration = par->getEIPtotalDuration();

    initialize_seasonality(par, community, nextMosquitoMultiplierIndex, nextEIPindex, date);
    vector<string> daily_output_buffer;
    if (par->bSecondaryTransmission and not par->abcVerbose) {
        daily_output_buffer.push_back("time,type,id,location,serotype,symptomatic,withdrawn,new_infection");
        //cout << "time,type,id,location,serotype,symptomatic,withdrawn,new_infection" << endl;
    }

    map<string, vector<int> > periodic_incidence { {"daily", vector<int>(3,0)},  // { introductions, local transmission, total }
                                                   {"weekly", vector<int>(3,0)},
                                                   {"monthly", vector<int>(3,0)},
                                                   {"yearly", vector<int>(3,0)} };
    for (; date.day() < par->nRunLength; date.increment()) {
        if ( date.endOfYear() and date.year() == 127 ) { // This should correspond to April 9 (day 99) of 1987
                                                         // for a 155 year simulation starting on day 100
            // calculate seroprevalence among 8-14 year old merida residents
            for (int id: serotested_ids) {
                const double seropos = community->getPersonByID(id)->getNumInfections() > 0 ? 1.0 : 0.0;
                seropos_87 += seropos;
            }
            seropos_87 /= serotested_ids.size();
        } 
        //date.print();

        // phased vaccination
        if (date.startOfYear()) {
            for (int i=0; i<par->nSizeVaccinate; i++) {
                if (date.year()==par->nVaccinateYear[i]) {
                    community->vaccinate(par->fVaccinateFraction[i],par->nVaccinateAge[i]);
                    if (not par->abcVerbose) {
                        cerr << "vaccinating " << par->fVaccinateFraction[i]*100 << "% of age " << par->nVaccinateAge[i] << endl;
                    }
                }
            }
        }


        {   // seed epidemic
            int numperson = community->getNumPerson();
            for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
                const int ai_year_lookup = date.year() % par->annualIntroductions.size();
                const double intros = par->annualIntroductions[ai_year_lookup];
                const int de_year_lookup = date.year() % par->nDailyExposed.size();
                const double serotype_weight = par->nDailyExposed[de_year_lookup][serotype];
                const double annual_intros_weight = par->annualIntroductionsCoef;
                const double expected_num_exposed = serotype_weight * annual_intros_weight * intros; 
                if (expected_num_exposed <= 0) continue;
                const int num_exposed = gsl_ran_poisson(RNG, expected_num_exposed);
                for (int i=0; i<num_exposed; i++) {
                    // gsl_rng_uniform_int returns on [0, numperson-1]
                    int transmit_to_id = gsl_rng_uniform_int(RNG, numperson) + 1; 
                    if (community->infect(transmit_to_id, (Serotype) serotype, date.day())) {
                        periodic_incidence["daily"][0]++;
                    }
                }
            }
        }

        // should the mosquito population change?
        if (par->mosquitoMultipliers.size() > 0) {
            const int nextMosquitoStart = par->mosquitoMultipliers[nextMosquitoMultiplierIndex].start;
            if ( ((date.day()+date.offset())%mosquitoMultiplierTotalDuration) == nextMosquitoStart) {
                //cerr << "updating mosquitoes on day " << date.day() << ", which is day " << date.julianDay() 
                //     << " of the year. Using index " << nextMosquitoMultiplierIndex << endl;
                community->setMosquitoMultiplier(par->mosquitoMultipliers[nextMosquitoMultiplierIndex].value);
                nextMosquitoMultiplierIndex = (nextMosquitoMultiplierIndex+1)%par->mosquitoMultipliers.size();
            }
        }

        // should the EIP change?
        if (par->extrinsicIncubationPeriods.size() > 0) {
            const int nextEIPstart = par->extrinsicIncubationPeriods[nextEIPindex].start;
            if ( ((date.day()+date.offset())%EIPtotalDuration) == nextEIPstart) {
                community->setExtrinsicIncubation(par->extrinsicIncubationPeriods[nextEIPindex].value);
                nextEIPindex = (nextEIPindex+1)%par->extrinsicIncubationPeriods.size();
            }
        }

        community->tick(date.day());

        for (int i=community->getNumPerson()-1; i>=0; i--) {
            Person *p = community->getPerson(i);
            if (p->isInfected(date.day()) and p->isNewlyInfected(date.day())) {
                ++periodic_incidence["daily"][2];
                /*
                // this is commented out because we don't usually output daily data,
                // and there's no point wasting ~100 mb of ram per process on long runs
                stringstream ss;
                ss << date.day() << ","
                   << p->getID() << ","
                   << p->getLocation(0)->getID() << ","
                   << 1 + (int) p->getSerotype() << ","
                   << (p->isSymptomatic(date.day())?1:0);
                daily_output_buffer.push_back(ss.str());
                */
            }
        }

        periodic_output(par, community, periodic_incidence, date, process_id, epi_sizes);

    }

    //write_daily_buffer(daily_output_buffer, process_id);
    return epi_sizes;
}


vector<int> simulate_epidemic(const Parameters* par, Community* community, const int process_id) {
    vector<int> dummy1 (1,1);
    double dummy2 = 0.0;
    return simulate_epidemic(par, community, process_id, dummy1, dummy2);
}


void write_immunity_file() {
    /*{
      epi_sizes.push_back(periodic_incidence["monthly"][2]);
      string imm_filename;
      stringstream ss_imm_filename;
      //we're starting in april, want to start numbering at 1
      ss_imm_filename << "immunity." << process_id << ".mon" << date.month();
      //ss_imm_filename << "immunity." << process_id << ".year" << date.year() << ".mon" << date.month();
      imm_filename = ss_imm_filename.str();
      int dummy = 0;
      write_immunity_by_age_file(par, community, dummy, imm_filename);
    }*/
}


void write_daily_buffer( vector<string>& buffer, const int process_id) {
    stringstream ss_filename;
    ss_filename << "daily_output." << process_id;
    ofstream file;
    file.open(ss_filename.str());

    for (string s: buffer) file << s << endl;
    file.close();
}


void write_immunity_by_age_file(const Parameters* par, const Community* community, const int year, string filename) {
    if (filename == "") {
        stringstream ss_filename;
        ss_filename << "imm_vs_age.year" << year;
        filename = ss_filename.str();
    }
    // total count, denv1, denv2, denv3, denv4, imm_to_one, imm_to_all
    vector< vector<int> > tally(NUM_AGE_CLASSES, vector<int>(NUM_OF_SEROTYPES+3, 0)); 
    for (int i = 0; i<community->getNumPerson(); ++i) {
        Person* p = community->getPerson(i);
        const int age = p->getAge();
        tally[age][0]++;
        const int numInfections = p->getNumInfections();
        for (int k = 0; k<numInfections; ++k) {
            const int s = (int) p->getSerotype(k);
            tally[age][s+1]++;
        }
        if (numInfections > 0) tally[age][NUM_OF_SEROTYPES+1]++;
        if (numInfections == NUM_OF_SEROTYPES) tally[age][NUM_OF_SEROTYPES+2]++;
    }

    ofstream file;
    file.open(filename);
    for (int a = 0; a<NUM_AGE_CLASSES; ++a) {
        file << a;
        for (int i: tally[a]) file << " " << i;
        file << endl;
    }
    file.close();
}


void write_immunity_file(const Parameters* par, const Community* community, const int int_label, string filename) {
    if (filename == "") {
        stringstream ss_filename;
        ss_filename << "immunity." << int_label;
        filename = ss_filename.str();
    }
    ofstream file;
    file.open(filename);
    file << "pid age imm1 imm2 imm3 imm4\n";
    for (int i = 0; i<community->getNumPerson(); ++i) {
        Person* p = community->getPerson(i);
        vector<int> infection_history(NUM_OF_SEROTYPES, 0); // 0 is no infection; 1 means last year, 2 means 2 years ago ...
//if (i == 0) cerr << "person 0 num infections: " << p->getNumInfections()  << " in year " << year << endl;
        for (int k = 0; k<p->getNumInfections(); ++k) {
            int s = (int) p->getSerotype(k);
            infection_history[s] = p->getInfectedTime(k) - par->nRunLength;
            // how many years ago did this person get infected?  within last 365 days = 1 year ago
            // infection time is relative to end of simulation of par->nRunLength days
            //infection_history[s] = 1 + (int) (par->nRunLength - p->getInfectedTime(k))/365;
//if (i == 0) cerr << "infection time: " << p->getInfectedTime()  << " output as: " << infection_history[s] << endl;
        }
        file << p->getID() << " " << p->getAge() << " ";
        for (auto sero: infection_history) file << sero << " ";
        file << endl;
    }
    file.close();
}


void daily_detailed_output(Community* community, int t) {
    // print out infectious mosquitoes
/*    for (int i=community->getNumInfectiousMosquitoes()-1; i>=0; i--) {
        Mosquito *p = community->getInfectiousMosquito(i);
        cout << t << ",mi," << p->getID() << "," << p->getLocation()->getID() << "," << "," << "," << endl;
    }
    // print out exposed mosquitoes
    for (int i=community->getNumExposedMosquitoes()-1; i>=0; i--) {
        Mosquito *p = community->getExposedMosquito(i);
        // "current" location
        cout << t << ",me," << p->getID() << "," << p->getLocation()->getID() << "," << 1 + (int) p->getSerotype() << "," << "," << "," << endl;
    }*/
    // print out infected people
    for (int i=community->getNumPerson()-1; i>=0; i--) {
        Person *p = community->getPerson(i);
        if (p->isInfected(t)) {
            // home location
            cout << t 
                 << ",p,"
                 << p->getID() << "," 
                 << p->getLocation(0)->getID() << "," 
                 << 1 + (int) p->getSerotype() << "," 
                 << (p->isSymptomatic(t)?1:0) << "," 
                 << (p->isWithdrawn(t)?1:0) << ","
                 << (p->isNewlyInfected(t)?1:0) << endl;
        }
    }
}

void write_output(const Parameters* par, Community* community, vector<int> numInitialSusceptible) {
    if (!par->bSecondaryTransmission) {
        // outputs
        //   number of secondary infections by serotype (4)
        //   number of households infected
        //   age of index case
        //   age(s) of secondary cases
        vector<int> numCurrentSusceptible = community->getNumSusceptible();
        for (int s=0; s<NUM_OF_SEROTYPES; s++) {
            cout << (numInitialSusceptible[s] - numCurrentSusceptible[s]) << " ";
        }
        //    cout << "secondary infections" << endl;

        int ages[100];
        int times[100];
        int numages=0;
        int indexage=-1;
        int homeids[100];
        int numhomes=0;
        // ages of infected people
        for (int i=community->getNumPerson()-1; i>=0; i--) {
            Person *p = community->getPerson(i);
            int t = p->getInfectedTime();
            if (t>=0) {
                if (t==0)
                    indexage = p->getAge();
                else {
                    ages[numages] = p->getAge();
                    times[numages] = t;
                    numages++;
                }
                int homeid = p->getHomeID();
                bool bFound = false;
                for (int j=0; j<numhomes; j++)
                    if (homeids[j]==homeid)
                        bFound = true;
                if (!bFound)
                    homeids[numhomes++] = homeid;
            }
        }
        cout << indexage << " " << numhomes << " " << numages;
        for (int i=0; i<numages; i++)
            cout << " " << ages[i];
        for (int i=0; i<numages; i++)
            cout << " " << times[i];
        cout << endl;
    }

    // output daily infected/symptomatic file
    if (par->dailyOutputFilename.length()>0) {
        cerr << "outputing daily infected/symptomatic information to " << par->dailyOutputFilename << endl;
        ofstream dailyOutputFile;
        dailyOutputFile.open(par->dailyOutputFilename.c_str());
        if(dailyOutputFile.fail()) {
            cerr << "ERROR: Daily file '" << par->dailyOutputFilename << "' cannot be open for writing." << endl;
            exit(-1);
        }
        dailyOutputFile << "day,newly infected DENV1,newly infected DENV2,newly infected DENV3,newly infected DENV4,"
                  << "newly symptomatic DENV1,newly symptomatic DENV2,newly symptomatic DENV3,newly symptomatic DENV4" << endl;
        vector< vector<int> > infected =    community->getNumNewlyInfected();
        vector< vector<int> > symptomatic = community->getNumNewlySymptomatic();
        for (int t=0; t<par->nRunLength; t++) {
            dailyOutputFile << t << ",";
            for (int i=0; i<NUM_OF_SEROTYPES; i++)   dailyOutputFile << infected[i][t] << ",";
            for (int i=0; i<NUM_OF_SEROTYPES-1; i++) dailyOutputFile << symptomatic[i][t] << ","; 
            dailyOutputFile << symptomatic[NUM_OF_SEROTYPES-1][t] << endl;
        }
        dailyOutputFile.close();
    }

    // output people file
    if (par->peopleOutputFilename.length()>0) {
        cerr << "outputing people information to " << par->peopleOutputFilename << endl;
        ofstream peopleOutputFile;
        peopleOutputFile.open(par->peopleOutputFilename.c_str());
        if(peopleOutputFile.fail()) {
            cerr << "ERROR: People file '" << par->peopleOutputFilename << "' cannot be open for writing." << endl;
            exit(-1);
        }
        peopleOutputFile << "pid,serotype,infectiontime,symptomtime,withdrawtime,recoverytime,immdenv1,immdenv2,immdenv3,immdenv4,vaccinated" << endl;
        for (int i=0; i<community->getNumPerson(); i++) {
            Person *p = community->getPerson(i);
            for (int j=p->getNumInfections()-1; j>=0; j--) {
                peopleOutputFile << p->getID() << "," 
                    << 1 + (int) p->getSerotype(j) << "," 
                    << p->getInfectedTime(j) << "," 
                    << p->getSymptomTime(j) << "," 
                    << p->getWithdrawnTime(j) << "," 
                    << p->getRecoveryTime(j) << ","; 
                for (int s = 0; s < NUM_OF_SEROTYPES; ++s) {
                    peopleOutputFile << (p->isSusceptible((Serotype) s)?0:1) << ","; 
                }
                    peopleOutputFile << (p->isVaccinated()?1:0) << endl;
            }
        }
        peopleOutputFile.close();
    }
}
