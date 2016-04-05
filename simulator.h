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
#include "sys/stat.h"

using namespace dengue::standard;
using namespace dengue::util;

enum ReportingType {
    INTRO,
    LOCAL,
    INFECTION,
    CASE,
    DSS,
    NUM_OF_REPORTING_TYPES
};

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
vector<int> simulate_epidemic(const Parameters* par, Community* community, const string process_id = "0");
void write_immunity_file(const Parameters* par, const Community* community, const string label, string filename, int runLength);
void write_immunity_by_age_file(const Parameters* par, const Community* community, const int year, string filename="");
void write_output(const Parameters* par, Community* community, vector<int> initial_susceptibles);
void write_daily_buffer( vector<string>& buffer, const int process_id, string filename);

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
                community->setExpectedExtrinsicIncubation(eip);
            }
            nextEIPindex = (nextEIPindex+1)%par->extrinsicIncubationPeriods.size();
        }
    }
}


void _aggregator(map<string, vector<int> >& periodic_incidence, string key) {
    for (unsigned int i = 0; i < periodic_incidence["daily"].size(); ++i) periodic_incidence[key][i] += periodic_incidence["daily"][i];
}


void _reporter(stringstream& ss, map<string, vector<int> > &periodic_incidence, const Parameters* par, const string process_id, const string label, const int value, string key) {
        ss << process_id << dec << " " << par->serial << label << value << " ";
        for (auto v: periodic_incidence[key]) ss << v << " ";
        for (auto v: par->reportedFraction) ss << v << " ";
}


void periodic_output(const Parameters* par, const Community* community, map<string, vector<int> >& periodic_incidence, const Date& date, const string process_id, vector<int>& epi_sizes) {
    stringstream ss;
    // local transmission              = total                                  - introductions
    periodic_incidence["daily"][LOCAL] = periodic_incidence["daily"][INFECTION] - periodic_incidence["daily"][INTRO];
    if (par->dailyOutput) {
        _reporter(ss, periodic_incidence, par, process_id, " day: ", date.day(), "daily");
        ss << community->getExpectedExtrinsicIncubation() << " " << community->getMosquitoMultiplier()*par->nDefaultMosquitoCapacity << endl;
    }

    if (par->weeklyOutput) {
        _aggregator(periodic_incidence, "weekly");
        if (date.endOfWeek()) {
            _reporter(ss, periodic_incidence, par, process_id, " week: ", date.week(), "weekly"); ss << endl;
            periodic_incidence["weekly"] = vector<int>(NUM_OF_REPORTING_TYPES, 0);
        }
    }

    if (par->monthlyOutput) {
        _aggregator(periodic_incidence, "monthly");
        if (date.endOfMonth()) {
            _reporter(ss, periodic_incidence, par, process_id, " month: ", date.julianMonth(), "monthly"); ss << endl;
            periodic_incidence["monthly"] = vector<int>(NUM_OF_REPORTING_TYPES, 0);
        }
    }

    // handle several things that happen yearly
    _aggregator(periodic_incidence, "yearly");
    if (date.endOfYear()) {
        if (par->abcVerbose) {
            cout << process_id << dec << " " << par->serial << " T: " << date.day() << " annual: "; 
            for (auto v: periodic_incidence["yearly"]) cout << v << " "; cout << endl;
        }

        epi_sizes.push_back(periodic_incidence["yearly"][2]);

        if (par->yearlyPeopleOutputFilename.length() > 0) write_yearly_people_file(par, community, date.day());
        if (par->yearlyOutput) _reporter(ss, periodic_incidence, par, process_id, " year: ", date.year(), "yearly"); ss << endl;
        periodic_incidence["yearly"] = vector<int>(NUM_OF_REPORTING_TYPES, 0);
    }

    periodic_incidence["daily"] = vector<int>(NUM_OF_REPORTING_TYPES, 0);
    string output = ss.str();
    fputs(output.c_str(), stderr);
}

void update_vaccinations(const Parameters* par, Community* community, const Date &date) {
    for (VaccinationEvent ve: par->vaccinationEvents) {
        // Normal, initial vaccination
        if (date.day() == ve.simDay) {
            if (not par->abcVerbose) cerr << "vaccinating " << ve.coverage*100 << "% of age " << ve.age << " on day " << ve.simDay << endl;
            community->vaccinate(ve);
        } else if (date.day() > ve.simDay) {
            // Re-vaccination via ...
            int timeSinceVac = date.day() - ve.simDay;
            const int doseInterval = par->vaccineDoseInterval;
            const int boostInterval = par->vaccineBoostingInterval;
            if (timeSinceVac % doseInterval == 0 and timeSinceVac / doseInterval <= par->numVaccineDoses) {
                // Multi-dose vaccine
                community->boost(date.day(), doseInterval, par->numVaccineDoses);
            } else if (par->linearlyWaningVaccine and par->vaccineBoosting and timeSinceVac % boostInterval == 0) {
                // Boosting
                community->boost(date.day(), boostInterval);
            }
        }
    }
}


int seed_epidemic(const Parameters* par, Community* community, const Date &date) {
    int introduced_infection_ct = 0;
    const int numperson = community->getNumPerson();
    for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
        const int ai_year_lookup = date.year() % par->annualIntroductions.size();
        const double intros = par->annualIntroductions[ai_year_lookup];
        const int de_year_lookup = date.year() % par->nDailyExposed.size();
        const double serotype_weight = par->nDailyExposed[de_year_lookup][serotype];
        const double annual_intros_weight = par->annualIntroductionsCoef;
        const double expected_num_exposed = serotype_weight * annual_intros_weight * intros;
        if (expected_num_exposed <= 0) continue;
        assert(expected_num_exposed <= numperson);
        const int num_exposed = gsl_ran_poisson(RNG, expected_num_exposed);
        for (int i=0; i<num_exposed; i++) {
            // gsl_rng_uniform_int returns on [0, numperson-1]
            int transmit_to_id = gsl_rng_uniform_int(RNG, numperson) + 1;
            if (community->infect(transmit_to_id, (Serotype) serotype, date.day())) {
                introduced_infection_ct++;
            }
        }
    }
    return introduced_infection_ct;
}


void update_mosquito_population(const Parameters* par, Community* community, const Date &date, int& nextMosquitoMultiplierIndex) {
    const int mosquitoMultiplierTotalDuration = par->getMosquitoMultiplierTotalDuration();
    // should the mosquito population change?
    if (par->mosquitoMultipliers.size() > 0) {
        const int nextMosquitoStart = par->mosquitoMultipliers[nextMosquitoMultiplierIndex].start;
        if ( ((date.day()+date.offset())%mosquitoMultiplierTotalDuration) == nextMosquitoStart) {
            //cerr << "updating mosquitoes on day " << date.day() << ", which is day " << date.julianDay()
            //     << " of the year. Using index " << nextMosquitoMultiplierIndex << endl;
            community->applyMosquitoMultiplier(par->mosquitoMultipliers[nextMosquitoMultiplierIndex].value);
            nextMosquitoMultiplierIndex = (nextMosquitoMultiplierIndex+1)%par->mosquitoMultipliers.size();
        }
    }
}


void update_extrinsic_incubation_period(const Parameters* par, Community* community, const Date &date, int& nextEIPindex) {
    // should the EIP change?
    const int EIPtotalDuration = par->getEIPtotalDuration();
    if (par->extrinsicIncubationPeriods.size() > 0) {
        const int nextEIPstart = par->extrinsicIncubationPeriods[nextEIPindex].start;
        if ( ((date.day()+date.offset())%EIPtotalDuration) == nextEIPstart) {
            community->setExpectedExtrinsicIncubation(par->extrinsicIncubationPeriods[nextEIPindex].value);
            nextEIPindex = (nextEIPindex+1)%par->extrinsicIncubationPeriods.size();
        }
    }
}


void advance_simulator(const Parameters* par, Community* community, Date &date, const string process_id, map<string, vector<int> > &periodic_incidence, int &nextMosquitoMultiplierIndex, int &nextEIPindex, vector<int> &epi_sizes) {
    periodic_incidence["daily"][INTRO] += seed_epidemic(par, community, date);
    update_mosquito_population(par, community, date, nextMosquitoMultiplierIndex);
    update_extrinsic_incubation_period(par, community, date, nextEIPindex);

    community->tick(date.day());

    // TODO - make it cleaner to iterate through pop
    for (int i=community->getNumPerson()-1; i>=0; i--) {
        Person *p = community->getPerson(i);
        // TODO - it should be sufficient to only check second conditional
        if (p->isInfected(date.day()) and p->isNewlyInfected(date.day())) {
            ++periodic_incidence["daily"][INFECTION];
            const Infection* infec = p->getInfection();
            if (infec->isSymptomatic()) ++periodic_incidence["daily"][CASE];
            if (infec->isSevere())      ++periodic_incidence["daily"][DSS];
        }
    }

    periodic_output(par, community, periodic_incidence, date, process_id, epi_sizes);
    return;
}


map<string, vector<int> > construct_tally() {
    // { introductions, local transmission, total, case, severe}
    map<string, vector<int> > periodic_incidence { {"daily", vector<int>(NUM_OF_REPORTING_TYPES,0)},
                                                   {"weekly", vector<int>(NUM_OF_REPORTING_TYPES,0)},
                                                   {"monthly", vector<int>(NUM_OF_REPORTING_TYPES,0)},
                                                   {"yearly", vector<int>(NUM_OF_REPORTING_TYPES,0)} };
    return periodic_incidence;
}


vector<int> simulate_epidemic(const Parameters* par, Community* community, const string process_id) {
    vector<int> epi_sizes;
    Date date(par);
    int nextMosquitoMultiplierIndex = 0;
    int nextEIPindex = 0;

    initialize_seasonality(par, community, nextMosquitoMultiplierIndex, nextEIPindex, date);
    vector<string> daily_output_buffer;

    if (par->bSecondaryTransmission and not par->abcVerbose) {
        //daily_output_buffer.push_back("time,type,id,location,serotype,symptomatic,withdrawn,new_infection");
        //cout << "time,type,id,location,serotype,symptomatic,withdrawn,new_infection" << endl;
    }

    map<string, vector<int> > periodic_incidence = construct_tally();

    for (; date.day() < par->nRunLength; date.increment()) {
        update_vaccinations(par, community, date); 
        advance_simulator(par, community, date, process_id, periodic_incidence, nextMosquitoMultiplierIndex, nextEIPindex, epi_sizes);
    }

    // write_daily_buffer(daily_output_buffer, process_id, dailyfilename);
    return epi_sizes;
}


vector<long double> simulate_who_fitting(const Parameters* par, Community* community, const string process_id, vector<int> &serotested_ids) {
    assert(serotested_ids.size() > 0);
    vector<long double> metrics;
    vector<int> epi_sizes;
    Date date(par);
    int nextMosquitoMultiplierIndex = 0;
    int nextEIPindex = 0;

    initialize_seasonality(par, community, nextMosquitoMultiplierIndex, nextEIPindex, date);
    vector<string> daily_output_buffer;

    if (par->bSecondaryTransmission and not par->abcVerbose) {
        daily_output_buffer.push_back("time,type,id,location,serotype,symptomatic,withdrawn,new_infection");
        //cout << "time,type,id,location,serotype,symptomatic,withdrawn,new_infection" << endl;
    }

    map<string, vector<int> > periodic_incidence = construct_tally();

    for (; date.day() < par->nRunLength; date.increment()) {
        if ( date.julianDay() == 99 ) {
                                                               // for a 135 year simulation
            // calculate seroprevalence among 9 year old merida residents
            double seropos_9yo = 0.0;
            for (int id: serotested_ids) {
                const double seropos = community->getPersonByID(id)->getNumInfections() > 0 ? 1.0 : 0.0;
                seropos_9yo += seropos;
            }
            seropos_9yo /= serotested_ids.size();
            metrics.push_back(seropos_9yo);
        }

        advance_simulator(par, community, date, process_id, periodic_incidence, nextMosquitoMultiplierIndex, nextEIPindex, epi_sizes);
    }

    return metrics;
}


vector<int> simulate_abc(const Parameters* par, Community* community, const string process_id, vector<int> &serotested_ids, double &seropos_87) {
    assert(serotested_ids.size() > 0);
    vector<int> epi_sizes;
    Date date(par);
    int nextMosquitoMultiplierIndex = 0;
    int nextEIPindex = 0;

    initialize_seasonality(par, community, nextMosquitoMultiplierIndex, nextEIPindex, date);
    vector<string> daily_output_buffer;

    if (par->bSecondaryTransmission and not par->abcVerbose) {
        daily_output_buffer.push_back("time,type,id,location,serotype,symptomatic,withdrawn,new_infection");
        //cout << "time,type,id,location,serotype,symptomatic,withdrawn,new_infection" << endl;
    }

    map<string, vector<int> > periodic_incidence = construct_tally();

    for (; date.day() < par->nRunLength; date.increment()) {
        if ( date.julianDay() == 99 and date.year() == 108 ) { // This should correspond to April 9 (day 99) of 1987
                                                               // for a 135 year simulation
            // calculate seroprevalence among 8-14 year old merida residents
            for (int id: serotested_ids) {
                const double seropos = community->getPersonByID(id)->getNumInfections() > 0 ? 1.0 : 0.0;
                seropos_87 += seropos;
            }
            seropos_87 /= serotested_ids.size();
        }

//        if (date.day() == 125*365) {
//            string imm_filename = "/scratch/lfs/thladish/imm_1000_yucatan/immunity2003." + process_id;
//            write_immunity_file(par, community, process_id, imm_filename, date.day());
//        }

        advance_simulator(par, community, date, process_id, periodic_incidence, nextMosquitoMultiplierIndex, nextEIPindex, epi_sizes);
    }

    return epi_sizes;
}


bool fileExists(const std::string& filename) {
    struct stat buf;
    return stat(filename.c_str(), &buf) != -1;
}


void write_daily_buffer( vector<string>& buffer, const int process_id, string filename = "" ) {
    if (filename == "") {
        stringstream ss_filename;
        ss_filename << "daily_output." << process_id;
        filename = ss_filename.str();
    }

    string all_output;
    for (const auto &line : buffer) all_output += (line + "\n");

    if (fileExists(filename)) {
        cerr << "WARNING: Daily output file already exists: " << filename << endl << "WARNING: Aborting write.\n";
        return;
    }

    ofstream file;
    file.open(filename);

    if (file.is_open()) {  // TODO - add this check everywhere a file is opened
        file << all_output;
        file.close();
    } else {
        cerr << "ERROR: Could not open daily buffer file for output: " << filename << endl;
        exit(-842);
    }
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


void write_immunity_file(const Parameters* par, const Community* community, const string label, string filename, int runLength) {
    if (filename == "") {
        stringstream ss_filename;
        ss_filename << "immunity." << label;
        filename = ss_filename.str();
    }
    ofstream file;
    file.open(filename);
    file << "pid age imm1 imm2 imm3 imm4\n";
    for (int i = 0; i<community->getNumPerson(); ++i) {
        Person* p = community->getPerson(i);
        vector<int> infection_history(NUM_OF_SEROTYPES, 0); // 0 is no infection; -1 means yesterday, -2 means 2 days ago ...
        for (int k = 0; k<p->getNumInfections(); ++k) {
            int s = (int) p->getSerotype(k);
            infection_history[s] = p->getInfectedTime(k) - runLength;
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
                 << p->getLocation(HOME_MORNING)->getID() << "," 
                 << 1 + (int) p->getSerotype() << "," 
                 << (p->isSymptomatic(t)?1:0) << "," 
                 << (p->isWithdrawn(t)?1:0) << ","
                 << (p->isNewlyInfected(t)?1:0) << endl;
        }
    }
}


void write_mosquito_location_data(const Community* community, string mos_filename, string loc_filename) {
    ofstream mos_file;
    mos_file.open(mos_filename);
    mos_file << "locID sero queue idx ageInfd ageInfs ageDead\n";
    const vector< vector<Mosquito*> > exposed = community->getExposedMosquitoes();
    // Exposed mosquitoes, by incubation days left
    for (unsigned int i = 0; i < exposed.size(); ++i) {
        const vector<Mosquito*>& mosquitoes = exposed[i];
        for (const Mosquito* m: mosquitoes) {
            mos_file << m->getLocation()->getID() << " " << m->getSerotype()    << " "
                     << "e " << i                 << " " << m->getAgeInfected() << " "
                     << m->getAgeInfectious()     << " " << m->getAgeDeath()    << endl;
        }
    }

    const vector< vector<Mosquito*> > infectious = community->getInfectiousMosquitoes();
    // Infectious mosquitoes, by days left to live
    for (unsigned int i = 0; i < infectious.size(); ++i) {
        const vector<Mosquito*>& mosquitoes = infectious[i];
        for (const Mosquito* m: mosquitoes) {
            mos_file << m->getLocation()->getID() << " " << m->getSerotype()    << " "
                     << "i " << i                 << " " << m->getAgeInfected() << " "
                     << m->getAgeInfectious()     << " " << m->getAgeDeath()    << endl;
        }
    }
    mos_file.close();

    ofstream loc_file;
    loc_file.open(loc_filename);
    loc_file << "locID baseMos infdMos\n";
    // Mosquitoes by location
    const vector<Location*> locations = community->getLocations();
    for (Location* loc: locations) {
        loc_file << loc->getID() << " ";
        loc_file << loc->getBaseMosquitoCapacity() << " ";
        loc_file << loc->getCurrentInfectedMosquitoes() << endl;
    }
    loc_file.close();
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
