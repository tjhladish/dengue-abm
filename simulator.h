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
#include "Date.h"

using namespace dengue::standard;
using namespace dengue::util;

enum IncidenceReportingType {
    INTRO_INF,
    TOTAL_INF,
    INTRO_CASE,
    TOTAL_CASE,
    INTRO_DSS,
    TOTAL_DSS,
    NUM_OF_INCIDENCE_REPORTING_TYPES
};

enum PrevalenceReportingType {
    INTRO_INF_PREV,
    TOTAL_INF_PREV,
    INTRO_CASE_PREV,
    TOTAL_CASE_PREV,
    INTRO_DSS_PREV,
    TOTAL_DSS_PREV,
    NUM_OF_PREVALENCE_REPORTING_TYPES
};

// metrics order: inf arm 1, sym arm 1, dss arm 1, inf arm 2, sym arm 2, dss arm 2, p10k aar
enum ProtoMetricType {
    ARM_1_INF,
    ARM_1_SYM,
    ARM_1_DSS,
    ARM_2_INF,
    ARM_2_SYM,
    ARM_2_DSS,
    AAR,
    NUM_OF_PROTO_METRIC_TYPES
};

// class Date {
//   public:
//     Date():_offset(0),_simulation_day(0) {};
//     Date(const Parameters* par):_offset(par->startDayOfYear-1),_simulation_day(0) {};
//
//     int offset()             const { return _offset; }
//     inline int day()         const { return _simulation_day; }                // [0, ...]
//     int julianDay()          const { return ((day() + offset()) % 365) + 1; } // [1, 365]
//     int dayOfMonth()         const { return julianMonth() == 1 ?              // [1, {29,30,31}]
//                                             julianDay() :
//                                             julianDay() - END_DAY_OF_MONTH[julianMonth()-2]; }
//     int nDayPeriod(int n)    const { return (int) day()/n; }                  // [0, ...]
//     int week()               const { return (int) day()/7; }                  // [0, ...]
//     int julianWeek()         const { return (int) ((julianDay()-1)/7) + 1; }  // [1, 53]
//     int month()              const { return _month_ct; }                      // [0, ...]
//     int julianMonth()        const {                                          // [1, 12]
//         vector<int>::const_iterator it;
//         // find first month that hasn't ended (hint: it's this month)
//         // julianDay()-1 because this isn't upper_or_equal_bound, which would be convenient
//         it = upper_bound(END_DAY_OF_MONTH.begin(), END_DAY_OF_MONTH.end(), julianDay()-1);
//         return it - END_DAY_OF_MONTH.begin() + 1; // +1 b/c [1, 12], not [0, 11]
//     }
//     string monthName()       const { return MONTH_NAMES[julianMonth()-1]; }
//     int year()               const { return (int) (day()/365); }
//
//     bool endOfPeriod(int n)  const { return (day()+1) % n == 0; }
//     bool endOfWeek()         const { return (day()+1) % 7 == 0; }
//     bool endOfMonth()        const {
//         vector<int>::const_iterator it;
//         // find out if today corresponds to a month-end
//         it = find(END_DAY_OF_MONTH.begin(), END_DAY_OF_MONTH.end(), julianDay());
//         return it != END_DAY_OF_MONTH.end();
//     }
//     bool startOfYear()       const { return day() % 365 == 0; }     // is it beginning of 365 day period
//     bool endOfYear()         const { return (day()+1) % 365 == 0; } // is it end of 365 day period
//     bool startOfJulianYear() const { return julianDay() == 1; }     // is it Jan 1
//     bool endOfJulianYear()   const { return julianDay() == 365; }   // is it Dec 31
//
//     void increment() {
//         if(endOfMonth()) _month_ct++;
//         _simulation_day++;
//     }
//
//     void print() {
//         cerr << day() << "\t" << julianDay() << "\t" << year()
//              << "\t(" << monthName() << " " << dayOfMonth() << ")\t" << month() << "\t" << julianMonth();
//         if (endOfWeek()) cerr << " EoW";
//         if (endOfMonth()) cerr << " EoM";
//         if (startOfYear()) cerr << " SoY";
//         if (endOfYear()) cerr << " EoY";
//         if (endOfJulianYear()) cerr << " EoJY";
//         cerr << endl;
//     }
//
//   private:
//     const int _offset;
//     int _simulation_day;
//     int _month_ct;
// };

const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);

// Predeclare local functions
Community* build_community(const Parameters* par);
void seed_epidemic(const Parameters* par, Community* community);
map<ProtoMetricType, vector<int>> simulate_epidemic(const Parameters* par, Community* community, const string process_id = "0");
void write_immunity_file(const Community* community, const string label, string filename, int runLength);
void write_immunity_by_age_file(const Community* community, const int year, string filename="");
void write_daily_buffer( vector<string>& buffer, const string process_id, string filename);

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
        cerr << community->getNumPeople() << " people" << endl;
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
                community->infect(gsl_rng_uniform_int(RNG, community->getNumPeople()), (Serotype) serotype,0);
        }
    }
    if (attempt_initial_infection) {
        for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
            // Useful for estimating R0
            if(par->nInitialInfected[serotype] > 0) {
                int count = community->getNumInfected(0);

                // must infect nInitialInfected persons -- this bit is mysterious
                while (community->getNumInfected(0) < count + par->nInitialInfected[serotype]) {
                    community->infect(gsl_rng_uniform_int(RNG, community->getNumPeople()), (Serotype) serotype,0);
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
    for (Person* p: community->getPeople()) {
        for (int j=p->getNumNaturalInfections()-1; j>=0; j--) {
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


void _aggregator(map<string, vector<int> >& periodic_incidence, string key, string daily_key = "daily") {
    for (unsigned int i = 0; i < periodic_incidence[daily_key].size(); ++i) periodic_incidence[key][i] += periodic_incidence[daily_key][i];
}


void _reporter(stringstream& ss, map<string, vector<int> > &periodic_incidence, vector<int> &periodic_prevalence, const Parameters* par, const string process_id, const string label, const int value, string key) {
        ss << process_id << dec << " " << par->serial << label << value << " ";
        for (auto v: periodic_incidence[key]) ss << v << " ";
        if(key=="daily") {
            ss << "| ";
            for (auto v: periodic_prevalence) { ss <<  v << " "; }
        }
        ss << "| ";
        for (auto v: par->reportedFraction) { ss << v << " "; }
}


void periodic_output(const Parameters* par, const Community* community, map<string, vector<int> > &periodic_incidence, vector<int> &periodic_prevalence, const Date& date,
                     const string process_id, map<ProtoMetricType, vector<int>> &proto_metrics) {
    stringstream ss;
//if (date.day() >= 25*365 and date.day() < 36*365) {
//if (date.day() >= 116*365) {                         // daily output starting in 1995, assuming Jan 1, 1879 simulation start
//if (date.day() >= 99*365 and date.day() < 105*365) { // daily output for summer/winter IRS comparison
    if (par->dailyOutput) {
        _reporter(ss, periodic_incidence, periodic_prevalence, par, process_id, " day: ", date.day(), "daily");
        ss << community->getExpectedExtrinsicIncubation() << " " << community->getMosquitoMultiplier()*par->nDefaultMosquitoCapacity << endl;
    }
//}
    vector<int> dummy;
     if (par->periodicOutput) {
        _aggregator(periodic_incidence, "n_day");
        const int n = par->periodicOutputInterval;
        if (date.endOfPeriod(n)) {
            _reporter(ss, periodic_incidence, dummy, par, process_id, " " + to_string(n) + "_day: ", date.nDayPeriod(n), "n_day"); ss << endl;
            periodic_incidence["n_day"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);
        }
    }

    if (par->weeklyOutput) {
        _aggregator(periodic_incidence, "weekly");
        if (date.endOfWeek()) {
            _reporter(ss, periodic_incidence, dummy, par, process_id, " week: ", date.week(), "weekly"); ss << endl;
            periodic_incidence["weekly"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);
        }
    }

    if (par->monthlyOutput) {
        _aggregator(periodic_incidence, "monthly");
        if (date.endOfMonth()) {
            _reporter(ss, periodic_incidence, dummy, par, process_id, " month: ", date.julianMonth(), "monthly"); ss << endl;
            periodic_incidence["monthly"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);
        }
    }

/*    if (par->simulateTrial) {
        _reporter(ss, periodic_incidence, periodic_prevalence, par, process_id, " day (arm 1 ): ", date.day(), "daily-arm1");
        ss << endl;
        _reporter(ss, periodic_incidence, periodic_prevalence, par, process_id, " day (arm 2 ): ", date.day(), "daily-arm2");
        ss << endl;
    }*/

    // handle several things that happen yearly
    _aggregator(periodic_incidence, "yearly");
    if (par->simulateTrial) {
        _aggregator(periodic_incidence, "yearly-arm1", "daily-arm1");
        _aggregator(periodic_incidence, "yearly-arm2", "daily-arm2");
    }

    if (date.endOfYear()) {
        if (par->abcVerbose) {
            cout << process_id << dec << " " << par->serial << " T: " << date.day() << " annual: ";
            for (auto v: periodic_incidence["yearly"]) { cout << v << " "; } cout << endl;
        }

        if (par->yearlyOutput) {
            _reporter(ss, periodic_incidence, dummy, par, process_id, " year ( total ): ", date.year(), "yearly"); ss << endl;
        }
        // periodic_incidence["yearly"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);

        if (par->simulateTrial) {
            proto_metrics[ARM_1_INF].push_back(periodic_incidence["yearly-arm1"][TOTAL_INF]);
            proto_metrics[ARM_1_SYM].push_back(periodic_incidence["yearly-arm1"][TOTAL_CASE]);
            proto_metrics[ARM_1_DSS].push_back(periodic_incidence["yearly-arm1"][TOTAL_DSS]);

            proto_metrics[ARM_2_INF].push_back(periodic_incidence["yearly-arm2"][TOTAL_INF]);
            proto_metrics[ARM_2_SYM].push_back(periodic_incidence["yearly-arm2"][TOTAL_CASE]);
            proto_metrics[ARM_2_DSS].push_back(periodic_incidence["yearly-arm2"][TOTAL_DSS]);

            if (par->yearlyOutput) {
                _reporter(ss, periodic_incidence, dummy, par, process_id, " year ( arm_1 ): ", date.year(), "yearly-arm1"); ss << endl;
                _reporter(ss, periodic_incidence, dummy, par, process_id, " year ( arm_2 ): ", date.year(), "yearly-arm2"); ss << endl;
            }

            periodic_incidence["yearly-arm1"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);
            periodic_incidence["yearly-arm2"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);
        } //else {
            // proto_metrics.push_back(periodic_incidence["yearly"][2]);
        //}

        proto_metrics[AAR].push_back(periodic_incidence["yearly"][TOTAL_INF]);

        if (par->yearlyPeopleOutputFilename.length() > 0) write_yearly_people_file(par, community, date.day());
        periodic_incidence["yearly"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);

    }

    periodic_incidence["daily"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);
    if (par->simulateTrial) {
        periodic_incidence["daily-arm1"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);
        periodic_incidence["daily-arm2"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);
    }

    periodic_prevalence         = vector<int>(NUM_OF_PREVALENCE_REPORTING_TYPES, 0);

    string output = ss.str();
    fputs(output.c_str(), stderr);
    //fputs(output.c_str(), stdout);
}

void update_vaccinations(const Parameters* par, Community* community, const Date &date) {
    const int doseInterval = par->vaccineDoseInterval;
    assert(doseInterval > 0); // neg is nonsensical, 0 is disallowed due to mod operation
    //const int boostInterval = par->vaccineBoostingInterval;
    for (CatchupVaccinationEvent cve: par->catchupVaccinationEvents) {
        // Normal, initial vaccination -- boosting, multiple doses handled in Community::tick()
        if (date.day() == cve.simDay) {
            if (not par->abcVerbose) cerr << "vaccinating " << cve.coverage*100 << "% of age " << cve.age << " on day " << cve.simDay << endl;
            community->vaccinate(cve);
        }
    }
}

void schedule_vector_control(const Parameters* par, Community* community) {
    for (VectorControlEvent vce: par->vectorControlEvents) {
        const string loc_label = vce.locationType == 0 ? "houses" : vce.locationType == 1 ? "workplaces" : vce.locationType == 2 ? "schools" : "unknown location type";
        if (not par->abcVerbose) cerr << "will start treating " << vce.coverage*100 << "% of " << loc_label << " on day " << vce.campaignStart << endl;

        double rho = par->calculate_daily_vector_control_mortality(vce.efficacy);
        if (vce.strategy == UNIFORM_STRATEGY) {
            for (Location* loc: community->getLocations() ) {
                if (loc->getType() == vce.locationType and  gsl_rng_uniform(RNG) < vce.coverage) {
                   // location will be treated
                   const int loc_treatment_date = vce.campaignStart + gsl_rng_uniform_int(RNG, vce.campaignDuration);
                   loc->scheduleVectorControlEvent(vce.efficacy, rho, loc_treatment_date, vce.efficacyDuration);
                }
            }
        } else if (vce.strategy == TIRS_STUDY_STRATEGY) {
            for (Location* loc: community->getLocations() ) {
                if (loc->getType() == vce.locationType and loc->getTrialArm() == 2) {
                   // location will be treated
                   const int loc_treatment_date = vce.campaignStart + gsl_rng_uniform_int(RNG, vce.campaignDuration);
                   loc->scheduleVectorControlEvent(vce.efficacy, rho, loc_treatment_date, vce.efficacyDuration);
                }
            }
        } else if (vce.strategy == MAX_MOSQUITOES_STRATEGY) { //  TODO - implement me
            cerr << "ERROR: MAX_MOSQUITOES_STRATEGY is not yet implemented\n";
            exit(-831);
        } else {
            cerr << "ERROR: Unsupported vector control strategy\n";
            exit(-832);
        }
    }
}


int seed_epidemic(const Parameters* par, Community* community, const Date &date) {
    int introduced_infection_ct = 0;
    const int numperson = community->getNumPeople();
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
            int transmit_to_id = gsl_rng_uniform_int(RNG, numperson);
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


void advance_simulator(const Parameters* par, Community* community, Date &date, const string process_id, map<string, vector<int> > &periodic_incidence, vector<int> &periodic_prevalence,
                       int &nextMosquitoMultiplierIndex, int &nextEIPindex, map<ProtoMetricType, vector<int>> &proto_metrics) {
    update_mosquito_population(par, community, date, nextMosquitoMultiplierIndex);
    update_extrinsic_incubation_period(par, community, date, nextEIPindex);
    community->tick(date);

    seed_epidemic(par, community, date);

    for (Person* p: community->getPeople()) {
        if (p->isInfected(date.day())) {
            const Infection* infec = p->getInfection();
            bool intro  = not infec->isLocallyAcquired();
            bool symp   = infec->isSymptomatic();
            bool severe = infec->isSevere();

            // Prevalence
            periodic_prevalence[INTRO_INF_PREV]  += intro;
            periodic_prevalence[INTRO_CASE_PREV] += intro and symp;
            periodic_prevalence[INTRO_DSS_PREV]  += intro and severe;

            periodic_prevalence[TOTAL_INF_PREV]  += 1;
            periodic_prevalence[TOTAL_CASE_PREV] += symp;
            periodic_prevalence[TOTAL_DSS_PREV]  += severe;

            // Incidence
            if (p->isNewlyInfected(date.day())) {
                periodic_incidence["daily"][INTRO_INF]  += intro;
                periodic_incidence["daily"][INTRO_CASE] += intro and symp;
                periodic_incidence["daily"][INTRO_DSS]  += intro and severe;

                periodic_incidence["daily"][TOTAL_INF]  += 1;
                periodic_incidence["daily"][TOTAL_CASE] += symp;
                periodic_incidence["daily"][TOTAL_DSS]  += severe;

                const Location* home = p->getLocation(HOME_MORNING);
                if (par->simulateTrial) { // TODO check whether p is a surveilled person (e.g. a child, for TIRS study)
                    const int age = p->getAge();
                    const TrialArmState arm = home->isSurveilled() and age >= 2 and age <= 15 ? home->getTrialArm() : NOT_IN_TRIAL;

                    if (arm == TRIAL_ARM_1) {
                        periodic_incidence["daily-arm1"][INTRO_INF]  += intro;
                        periodic_incidence["daily-arm1"][INTRO_CASE] += intro and symp;
                        periodic_incidence["daily-arm1"][INTRO_DSS]  += intro and severe;

                        periodic_incidence["daily-arm1"][TOTAL_INF]  += 1;
                        periodic_incidence["daily-arm1"][TOTAL_CASE] += symp;
                        periodic_incidence["daily-arm1"][TOTAL_DSS]  += severe;
                    } else if (arm == TRIAL_ARM_2) {
                        periodic_incidence["daily-arm2"][INTRO_INF]  += intro;
                        periodic_incidence["daily-arm2"][INTRO_CASE] += intro and symp;
                        periodic_incidence["daily-arm2"][INTRO_DSS]  += intro and severe;

                        periodic_incidence["daily-arm2"][TOTAL_INF]  += 1;
                        periodic_incidence["daily-arm2"][TOTAL_CASE] += symp;
                        periodic_incidence["daily-arm2"][TOTAL_DSS]  += severe;
                    }
                }
            }
        }
    }

    periodic_output(par, community, periodic_incidence, periodic_prevalence, date, process_id, proto_metrics);
    return;
}



map<string, vector<int> > construct_tally() {
    // { introductions, local transmission, total, case, severe}
    map<string, vector<int> > periodic_incidence { {"daily", vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES,0)},
                                                   {"n_day", vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES,0)},
                                                   {"weekly", vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES,0)},
                                                   {"monthly", vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES,0)},
                                                   {"yearly", vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES,0)},
                                                   {"yearly-arm1", vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES,0)},
                                                   {"yearly-arm2", vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES,0)},
                                                   {"daily-arm1", vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES,0)},
                                                   {"daily-arm2", vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES,0)}};
    return periodic_incidence;
}


map<ProtoMetricType, vector<int>> simulate_epidemic_with_seroprev(const Parameters* par, Community* community, const string process_id, bool capture_sero_prev, vector< vector<double> > &sero_prev, int sero_prev_aggregation_julian_start=0) {
    sero_prev = vector< vector<double> > (5, vector<double>(par->nRunLength/365, 0.0)); // rows are infection history: 0, 1, and 2+ infections
    map<ProtoMetricType, vector<int>> proto_metrics;
    Date date(par);
    int nextMosquitoMultiplierIndex = 0;
    int nextEIPindex = 0;

    initialize_seasonality(par, community, nextMosquitoMultiplierIndex, nextEIPindex, date);
    schedule_vector_control(par, community);
    vector<string> daily_output_buffer;

    if (par->bSecondaryTransmission and not par->abcVerbose) {
        daily_output_buffer.push_back("day,year,id,age,location,vaccinated,serotype,symptomatic,severity");
    }

    map<string, vector<int> > periodic_incidence = construct_tally();
    vector<int> periodic_prevalence(NUM_OF_PREVALENCE_REPORTING_TYPES, 0);

    for (; date.day() < par->nRunLength; date.increment()) {
        update_vaccinations(par, community, date);
        advance_simulator(par, community, date, process_id, periodic_incidence, periodic_prevalence, nextMosquitoMultiplierIndex, nextEIPindex, proto_metrics);
        if (capture_sero_prev and ((int) date.julianDay() == ((sero_prev_aggregation_julian_start+364) % 365 ) + 1)) { // +1 because julianDay is [1,365])), avg(avg(interventions are specified on [0,364]
            // tally current seroprevalence stats
            int vaccinated_tally = 0;
            vector<double> inf_ct_tally(5,0.0);         // [0,4] past infections
            for (Person* p: community->getPeople()) {
                // ignoring possible seroconversion due to vaccine
                vaccinated_tally += p->isVaccinated();
                const int ct = p->getNumNaturalInfections();
                ++inf_ct_tally[ct];
            }
            double N = community->getNumPeople();
            for (unsigned int i = 0; i < inf_ct_tally.size(); ++i) sero_prev[i][date.year()] = inf_ct_tally[i] / N;
            cerr << "vaccinated N, fraction: " << vaccinated_tally << " " << (double) vaccinated_tally / N << endl;
        }
    }
/*
    stringstream ss_filename;
    ss_filename << "/scratch/lfs/thladish/who-feb-2016/daily." << process_id << "." << par->randomseed;
    string dailyfilename = ss_filename.str();
    write_daily_buffer(daily_output_buffer, process_id, dailyfilename);
*/
    return proto_metrics;
}


map<ProtoMetricType, vector<int>> simulate_epidemic(const Parameters* par, Community* community, const string process_id) {
    vector< vector<double> > sero_prev;
    bool capture_sero_prev = false;
    return simulate_epidemic_with_seroprev(par, community, process_id, capture_sero_prev, sero_prev);
}


// vector<long double> simulate_who_fitting(const Parameters* par, Community* community, const string process_id, vector<int> &serotested_ids) {
//     assert(serotested_ids.size() > 0);
//     vector<long double> metrics;
//     vector<int> proto_metrics;
//     Date date(par);
//     int nextMosquitoMultiplierIndex = 0;
//     int nextEIPindex = 0;
//
//     initialize_seasonality(par, community, nextMosquitoMultiplierIndex, nextEIPindex, date);
//     vector<string> daily_output_buffer;
//
//     if (par->bSecondaryTransmission and not par->abcVerbose) {
//         daily_output_buffer.push_back("time,type,id,location,serotype,symptomatic,withdrawn,new_infection");
//         //cout << "time,type,id,location,serotype,symptomatic,withdrawn,new_infection" << endl;
//     }
//
//     map<string, vector<int> > periodic_incidence = construct_tally();
//     vector<int> periodic_prevalence(NUM_OF_PREVALENCE_REPORTING_TYPES, 0);
//
//     for (; date.day() < par->nRunLength; date.increment()) {
//         if ( date.julianDay() == 99 ) {
//                                                                // for a 135 year simulation
//             // calculate seroprevalence among 9 year old merida residents
//             double seropos_9yo = 0.0;
//             for (int id: serotested_ids) {
//                 const double seropos = community->getPersonByID(id)->getNumNaturalInfections() > 0 ? 1.0 : 0.0;
//                 seropos_9yo += seropos;
//             }
//             seropos_9yo /= serotested_ids.size();
//             metrics.push_back(seropos_9yo);
//         }
//
//         advance_simulator(par, community, date, process_id, periodic_incidence, periodic_prevalence, nextMosquitoMultiplierIndex, nextEIPindex, proto_metrics);
//     }
//
//     return metrics;
// }
//
// //c('0-4', '5-9', '10-14', '15-19', '20-29', '30-39', '40-49', '50-59', '60+')
//
// vector<int> simulate_abc(const Parameters* par, Community* community, const string process_id, vector<int> &serotested_ids_87, double &seropos_87, vector<int> &serotested_ids_14, vector<double> &seropos_14_by_age) {
//     assert(serotested_ids_87.size() > 0);
//     vector<int> proto_metrics;
//     Date date(par);
//     int nextMosquitoMultiplierIndex = 0;
//     int nextEIPindex = 0;
//
//     initialize_seasonality(par, community, nextMosquitoMultiplierIndex, nextEIPindex, date);
//     vector<string> daily_output_buffer;
//
//     if (par->bSecondaryTransmission and not par->abcVerbose) {
//         daily_output_buffer.push_back("day,year,id,age,location,vaccinated,serotype,symptomatic,severe");
//     }
//
//     map<string, vector<int> > periodic_incidence = construct_tally();
//     vector<int> periodic_prevalence(NUM_OF_PREVALENCE_REPORTING_TYPES, 0);
//
//     const vector<int> upper_age_bound_14 = {4, 9, 14, 19, 29, 39, 49, 59, INT_MAX};
//     vector<int> seropos_14_sample_size(upper_age_bound_14.size(), 0);
//     assert(upper_age_bound_14.size() == seropos_14_by_age.size());
//
//     for (; date.day() < par->nRunLength; date.increment()) {
//         if ( date.julianDay() == 99 and date.year() == 108 ) { // This corresponds to April 9 (day 99) of 1987
//                                                                // for a simulation starting Jan 1, 1879
//             cerr << "1987 serosurvey\n";
//             // calculate seroprevalence among 8-14 year old merida residents
//             for (int id: serotested_ids_87) {
//                 const double seropos = community->getPersonByID(id)->getNumNaturalInfections() > 0 ? 1.0 : 0.0;
//                 seropos_87 += seropos;
//             }
//             seropos_87 /= serotested_ids_87.size();
//         } else if ( date.julianDay() == 99 and date.year() == 135 ) { // This corresponds to April 9 (day 99) of 2014
//             cerr << "2014 serosurvey\n";
//             // calculate seroprevalence among all merida residents
//             for (int id: serotested_ids_14) {
//                 const Person* p = community->getPersonByID(id);
//                 const int age = p->getAge();
//                 assert(age >= 0);
//                 unsigned int age_cat;
//                 for (age_cat = 0; age_cat<upper_age_bound_14.size() and age>upper_age_bound_14[age_cat]; ++age_cat) {/*this space intentionally left blank*/}
//                 const double seropos = p->getNumNaturalInfections() > 0 ? 1.0 : 0.0;
//                 seropos_14_by_age[age_cat] += seropos;
//                 seropos_14_sample_size[age_cat]++;
//             }
//             for (unsigned int age_cat = 0; age_cat < seropos_14_by_age.size(); ++age_cat) {
//                 seropos_14_by_age[age_cat] /= seropos_14_sample_size[age_cat];
//             }
//         }
//         advance_simulator(par, community, date, process_id, periodic_incidence, periodic_prevalence, nextMosquitoMultiplierIndex, nextEIPindex, proto_metrics);
//
// /*        if ( date.julianDay() == 365 and date.year() == 121 ) { // December 31 (day 365) of 2000
//             string imm_filename = "/ufrc/longini/tjhladish/imm_1000_yucatan-irs_refit/immunity2000." + process_id;
//             write_immunity_file(community, process_id, imm_filename, date.day());
//         }
//
//         if ( date.julianDay() == 99 and date.year() == 135 ) { // April 9 (day 99) of 2014
//             string imm_filename = "/ufrc/longini/tjhladish/imm_1000_yucatan-irs_refit/immunity2014_04_09." + process_id;
//             write_immunity_file(community, process_id, imm_filename, date.day());
//         }*/
//     }
//
//     //string dailyfilename = "";
//     //write_daily_buffer(daily_output_buffer, process_id, dailyfilename);
//
//     return proto_metrics;
// }


bool fileExists(const std::string& filename) {
    struct stat buf;
    return stat(filename.c_str(), &buf) != -1;
}


void write_daily_buffer( vector<string>& buffer, const string process_id, string filename = "" ) {
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


void write_immunity_by_age_file(const Community* community, const int year, string filename) {
    if (filename == "") {
        stringstream ss_filename;
        ss_filename << "imm_vs_age.year" << year;
        filename = ss_filename.str();
    }
    // total count, denv1, denv2, denv3, denv4, imm_to_one, imm_to_all
    vector< vector<int> > tally(NUM_AGE_CLASSES, vector<int>(NUM_OF_SEROTYPES+3, 0));
    for (Person* p: community->getPeople()) {
        const int age = p->getAge();
        tally[age][0]++;
        const int numInfections = p->getNumNaturalInfections();
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


void write_immunity_file(const Community* community, const string label, string filename, int runLength) {
    if (filename == "") {
        stringstream ss_filename;
        ss_filename << "immunity." << label;
        filename = ss_filename.str();
    }
    ofstream file;
    file.open(filename);
    file << "pid age imm1 imm2 imm3 imm4\n";
    for (Person* p: community->getPeople()) {
        vector<int> infection_history(NUM_OF_SEROTYPES, 0); // 0 is no infection; -1 means yesterday, -2 means 2 days ago ...
        for (int k = 0; k<p->getNumNaturalInfections(); ++k) {
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
    for (Person *p: community->getPeople()) {
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

void import_csv_to_db(string filename, string table, string db) {
    stringstream ss;
    ss << "sqlite3 " << db << " '.mode csv' '.import " << filename << ' ' << table << "'";
    string cmd_str = ss.str();
    int retval = system(cmd_str.c_str());
    if (retval == -1) { cerr << "System failed to import " << table << " data to db\n"; }
    return;
}

void generate_sim_data_db(const Parameters* /*par*/, const Community* community, const unsigned long int serial, vector<string> tables) {
    vector<stringstream> filenames(tables.size());
    for (size_t i = 0; i < tables.size(); ++i) {
        filenames[i] << "./" << tables[i] << "_" << serial << ".csv";
    }

    map<string, ofstream> ofiles;
    for (size_t i = 0; i < tables.size(); ++i) {
        ofiles[tables[i]] = ofstream(filenames[i].str(), std::ios::trunc);
    }

    bool all_files_open = true;
    for (const auto& [table, ofile] : ofiles) {
        all_files_open = all_files_open and (bool)ofile;
    }
    if (not all_files_open) { cerr << "FILES FAILED TO OPEN" << endl; exit(-1); }

    for (Person* p : community->getPeople()) {
        for (Infection* inf : p->getInfectionHistory()) {
            if (not inf) { continue; }
            int inf_place_id = inf->getInfectedLoc() ? inf->getInfectedLoc()->getID() : -1;
            int inf_by_id    = inf->getInfectedBy() ? inf->getInfectedBy()->getID() : -1;
            int inf_owner_id = inf->getInfectionOwner() ? inf->getInfectionOwner()->getID() : -1;

            ofiles["infection_history"] << inf << ','
                                        << inf_place_id << ','
                                        << inf_by_id << ','
                                        << inf_owner_id << ','
                                        << inf->getInfectedTime() << ','
                                        << inf->getInfectiousTime() << ','
                                        << inf->getSymptomTime() << ','
                                        << inf->getRecoveryTime() << ','
                                        << inf->getWithdrawnTime() << ','
                                        << inf->isSevere() << ','
                                        << inf->serotype() << endl;
        }
    }
    for (auto& [table, ofile] : ofiles) { ofile.close(); }

    stringstream ss, db;
    db << "sim_data_" << serial << ".sqlite";
    ss << "sqlite3 " << db.str() << " '.read gen_sim_db.sql'";

    string cmd_str = ss.str();
    int retval = system(cmd_str.c_str());
    if (retval == -1) { cerr << "System call to `sqlite3 " << db.str() << " '.read gen_sim_db.sql'` failed\n"; }

    for (size_t i = 0; i < tables.size(); ++i) {
        import_csv_to_db(filenames[i].str(), tables[i], db.str());
    }

    ss.str(string());
    ss << "rm";
    for (size_t i = 0; i < filenames.size(); ++i) { ss << ' ' << filenames[i].str(); }

    cmd_str = ss.str();
    retval = system(cmd_str.c_str());
    if (retval == -1) { cerr << "System call to delete infection csv files failed\n"; }
}
