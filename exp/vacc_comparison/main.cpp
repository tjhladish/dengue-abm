#include <unistd.h>
#include <AbcSmc/AbcSmc.h>
#include "simulator.h"
#include <cstdlib>
#include "CCRC32.h"
#include "Utility.h"
#include "ImmunityGenerator.h"
#include <unordered_set>

using namespace std;

using dengue::util::to_string;
using dengue::util::mean;
using dengue::util::stdev;
using dengue::util::max_element;

time_t GLOBAL_START_TIME;

string calculate_process_id(vector< long double> &args, string &argstring);
const string SIM_POP = "merida";
const string HOME(std::getenv("HOME"));
const string pop_dir = HOME + "/work/dengue-abm/pop-" + SIM_POP;
const string output_dir("/ufrc/longini/tjhladish");
//const string output_dir("/scratch/lfs/thladish");

const int FIRST_YEAR          = 1879;                                 // inclusive
const int FIRST_OBSERVED_YEAR = 1979;
const int LAST_OBSERVED_YEAR  = 2013;                                 // inclusive
const int DDT_START           = 1956 - FIRST_YEAR;                    // simulator year 77, counting from year 1
const int DDT_DURATION        = 23;                                   // 23 years long
const int FITTED_DURATION     = LAST_OBSERVED_YEAR - FIRST_OBSERVED_YEAR + 1;  // 35 for 1979-2013 inclusive

//const int RESTART_BURNIN    = 80;
const int FORECAST_DURATION = 5;
//const int TOTAL_DURATION    = RESTART_BURNIN + FORECAST_DURATION;

Parameters* define_simulator_parameters(vector<long double> args, const unsigned long int rng_seed) {
    Parameters* par = new Parameters();
    par->define_defaults();

//  0 mildEF real,
//  1 severeEF real,
//  2 sec_sev real,
//  3 base_path real,
//  3 pss_ratio real,
//  4 exp_coef real,
//  5 num_mos real,
//  6 beta real,
//
//  Pseudo pars
//  7 beta_multiplier real,

    // The following FOI targets are the needed product of # mosquitoes * beta in order to achieve 9 yo equillibrium seroprevalences
    // of 10, 30, 50, 70, and 90%
    //    vector<float> foi_targets = {6.5, 7.5, 9, 12, 20}; // NEEDS TO BE UPDATED IF TRANSMISSION MODEL CHANGES
    //const vector<float> foi_targets = {6.2, 8.4, 10.4, 12.5, 20}; // NEEDS TO BE UPDATED IF TRANSMISSION MODEL CHANGES
    const vector<float> foi_targets = {5.7, 7.8, 10.0, 13.5, 21}; // NEEDS TO BE UPDATED IF TRANSMISSION MODEL CHANGES

//0.11739814169|0.0
//0.33525855234|1.0
//0.5314984042|2.0
//0.646959106000001|3.0
//0.868876936|4.0


//0.111192782416|0.0
//0.20122664603|1.0
//0.3716522091|2.0
//0.668049612|3.0
//0.894214538|4.0

    double _mild_EF      = args[0];
    double _severe_EF    = args[1];
    double _base_path    = args[2];
    double _sec_severity = args[3];
    double _pss_ratio    = args[4];
    double _exp_coef     = args[5];
    double _nmos         = args[6];
    double _betamp       = foi_targets[(int) args[7]] / _nmos;   // mp and pm and not separately
    double _betapm       = _betamp;                              // identifiable, so they're the same
    par->reportedFraction = {0.0, 1.0/_mild_EF, 1.0/_severe_EF}; // no asymptomatic infections are reported

    vector<long double> abc_args(&args[0], &args[6]);
    string argstring;
    const string process_id = calculate_process_id(abc_args, argstring);

    par->randomseed              = rng_seed;
    par->dailyOutput             = false;
    par->monthlyOutput           = false;
    par->yearlyOutput            = true;
    par->abcVerbose              = false; // needs to be false to get WHO daily output
    //const int runLengthYears     = RESTART_BURNIN + FORECAST_DURATION;
    const int runLengthYears     = DDT_START + DDT_DURATION + FITTED_DURATION + FORECAST_DURATION;
    par->nRunLength              = runLengthYears*365;
    par->startDayOfYear          = 100;
    par->annualIntroductionsCoef = pow(10, _exp_coef);

    // pathogenicity values fitted in
    // Reich et al, Interactions between serotypes of dengue highlight epidemiological impact of cross-immunity, Interface, 2013
    // Normalized from Fc values in supplement table 2, available at
    // http://rsif.royalsocietypublishing.org/content/10/86/20130414/suppl/DC1
//par->defineSerotypeRelativeRisks();
par->serotypePathogenicityRelativeRisks = vector<double>(NUM_OF_SEROTYPES, 1.0);
    par->useAgeStructuredPrimaryPathogenicity = false;
    par->basePathogenicity = _base_path;
    par->primaryRelativeRisk = 0.45/_base_path;
    par->postSecondaryRelativeRisk = 1.0/6.0;

    par->primarySevereFraction    = vector<double>(NUM_OF_SEROTYPES, _sec_severity*_pss_ratio);
    par->secondarySevereFraction  = vector<double>(NUM_OF_SEROTYPES, _sec_severity);
    par->tertiarySevereFraction   = vector<double>(NUM_OF_SEROTYPES, _sec_severity/4.0);
    par->quaternarySevereFraction = vector<double>(NUM_OF_SEROTYPES, _sec_severity/4.0);

    par->betaPM = _betapm;
    par->betaMP = _betamp;
    par->fMosquitoMove = 0.15;
    par->mosquitoMoveModel = "weighted";
    par->fMosquitoTeleport = 0.0;
    par->nDefaultMosquitoCapacity = (int) _nmos;
    par->eMosquitoDistribution = EXPONENTIAL;

    par->nDaysImmune = 730;
    par->fVESs = vector<double>(NUM_OF_SEROTYPES, 0.0);

    //par->hospitalizedFraction = 0.25; // fraction of cases assumed to be hospitalized

    par->simulateAnnualSerotypes = true;
    //par->normalizeSerotypeIntros = true;
    bool use_exact_known_seros = true;
    bool assume_future_hyperendemicity = true;
    par->nDailyExposed = generate_serotype_sequences(RNG, FIRST_YEAR, FIRST_OBSERVED_YEAR, LAST_OBSERVED_YEAR + 100, LAST_OBSERVED_YEAR, use_exact_known_seros, assume_future_hyperendemicity, TRANSPOSE);
    if (par->simulateAnnualSerotypes) par->generateAnnualSerotypes(runLengthYears + 100);
    //par->nDailyExposed = {{0.25, 0.25, 0.25, 0.25}};

    par->annualIntroductions = {1.0};

    // load daily EIP
    // EIPs calculated by ../raw_data/weather/calculate_daily_eip.R, based on reconstructed yucatan temps
    par->loadDailyEIP(pop_dir + "/seasonal_EIP_24hr.out");

    // we add the startDayOfYear offset because we will end up discarding that many values from the beginning
    // mosquitoMultipliers are indexed to start on Jan 1
    par->loadDailyMosquitoMultipliers(pop_dir + "/mosquito_seasonality.out", par->nRunLength + par->startDayOfYear);

    par->populationFilename       = pop_dir    + "/population-"         + SIM_POP + ".txt";
    par->immunityFilename         = "";
    par->locationFilename         = pop_dir    + "/locations-"          + SIM_POP + ".txt";
    par->networkFilename          = pop_dir    + "/network-"            + SIM_POP + ".txt";
    par->swapProbFilename         = pop_dir    + "/swap_probabilities-" + SIM_POP + ".txt";
    //par->mosquitoFilename         = output_dir + "/mos_mer_who/mos."       + to_string(process_id);
    //par->mosquitoLocationFilename = output_dir + "/mosloc_mer_who/mosloc." + to_string(process_id);

    return par;
}

// Take a list of values, return original indices sorted by value
vector<int> ordered(vector<int> const& values) {

    vector<pair<int,int> > pairs(values.size());
    for(unsigned int pos=0; pos<values.size(); pos++) {
        pairs[pos] = make_pair(values[pos],pos);
    }

    //bool comparator ( const mypair& l, const mypair& r) { return l.first < r.first; }
    std::sort( pairs.rbegin(), pairs.rend() ); // sort greatest to least
    vector<int> indices(values.size());
    for(unsigned int i=0; i < pairs.size(); i++) indices[i] = pairs[i].second;

    return indices;
}


string calculate_process_id(vector< long double> &args, string &argstring) {
    // CCRC32 checksum based on string version of argument values
    CCRC32 crc32;
    crc32.Initialize();

    for (unsigned int i = 0; i < args.size(); i++) argstring += to_string((long double) args[i]) + " ";

    const unsigned char* argchars = reinterpret_cast<const unsigned char*> (argstring.c_str());
    const int len = argstring.length();
    const int process_id = crc32.FullCRC(argchars, len);

    return to_string(process_id);
}


string report_process_id (vector<long double> &args, const unsigned long int serial, const ABC::MPI_par* mp, const time_t start_time) {
    double dif = difftime (start_time, GLOBAL_START_TIME);

    string argstring;
    const string process_id = calculate_process_id(args, argstring);

    stringstream ss;
    ss << mp->mpi_rank << " begin " << process_id << " " << dec << serial << " " << dif << " " << argstring << endl;
    string output = ss.str();
    fputs(output.c_str(), stderr);

    return to_string(process_id);
}


void append_if_finite(vector<long double> &vec, double val) {
    if (isfinite(val)) { 
        vec.push_back((long double) val);
    } else {
        vec.push_back(0);
    }
}


vector<long double> tally_counts(const Parameters* par, Community* community) {
    const int discard_days = 365*RESTART_BURNIN;
    vector< vector<int> > infected    = community->getNumNewlyInfected();
    const int num_years = FORECAST_DURATION;
    vector<vector<int> > i_tally(NUM_OF_SEROTYPES, vector<int>(num_years+1, 0));

    vector<long double> metrics(num_years, 0.0);
    for (int t=discard_days; t<par->nRunLength; t++) {
        // use epidemic years, instead of calendar years
        const int y = (t-discard_days)/365;
        for (int s=0; s<NUM_OF_SEROTYPES; s++) {
            i_tally[s][y] += infected[s][t];
        }
    }
    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = 0; y<num_years; ++y) {
            metrics[y] += i_tally[s][y];
        }
    }

    vector<int> pop_sizes  = vector<int>(10,0);
    vector<double> seropos = vector<double>(10,0.0); // seronegative at vaccination time

    for (int i=community->getNumPerson()-1; i>=0; i--) {
        Person *p = community->getPerson(i);
        const int final_age  = p->getAge();
        const int age_at_vac_start = final_age - FORECAST_DURATION;
        if (age_at_vac_start >= 9 and age_at_vac_start <= 18) {
            int age_class = age_at_vac_start - 9;

            const int serosurvey_day = (RESTART_BURNIN - age_class)*365;
            const bool ever_infected = p->getNumInfections() > 0 ? true : false;
            const bool isSeropos = ever_infected and p->getInfectionHistory().front()->getInfectedTime() < serosurvey_day;
            ++pop_sizes[age_class];
            if (isSeropos) ++seropos[age_class];
        }
    }

    long double seroprev_mean = 0.0;
    for (unsigned int i = 0; i < seropos.size(); ++i) seroprev_mean += seropos[i] / pop_sizes[i];
    assert(seropos.size() > 0);
    seroprev_mean /= seropos.size();

    metrics.push_back(seroprev_mean);
    return metrics;
}


vector<long double> simulator(vector<long double> args, const unsigned long int rng_seed, const unsigned long int serial, const ABC::MPI_par* mp) {
    gsl_rng_set(RNG, rng_seed); // seed the rng using sys time and the process id

    //initialize bookkeeping for run
    time_t start, end;
    time (&start);
    const string process_id = report_process_id(args, serial, mp, start);

    cerr << "SCENARIO " << rng_seed;
    for (auto _p: args) cerr << " " << _p; cerr << endl;

    Parameters* par = define_simulator_parameters(args, rng_seed);

    //vector<int> vaccine_durations = {-1, 10*365};
//  Posterior pars
//  0 mildEF real,
//  1 severeEF real,
//  2 base_path real
//  3 sec_sev real,
//  4 pss_ratio real
//  5 exp_coef real,
//  6 num_mos real,
//
//  Pseudo pars
//  7 foi_target real,
//  8 catchup int,
//  9 vaccine_mechanism int,
//  10 vaccine int,
//  11 target int,
//  12 catchup_to int,
//  13 coverage int 

    //vector<int> target_ages = {9, 16};         // default 9
    vector<int> catchup_ages = {17, 30};         // default 17
    vector<double> coverage_levels = {0.5, 0.8}; // default 0.8

    bool catchup           = (bool) args[8];
    int vaccine_mechanism  = (int) args[9];

    bool vaccine           = (bool) args[10];
    int target             = (int) args[11];
    int catchup_to         = catchup_ages[(int) args[12]];
    double coverage        = coverage_levels[(int) args[13]];

    bool nonsensical_parameters = false;
    // only run a non-vaccination campaign if all the vaccine parameters are 0
    // TODO - this should be reworked with new "catchup-to" parameter
    if (not vaccine and catchup) { nonsensical_parameters = true; }
    if (nonsensical_parameters) {
        vector<long double> dummy(FORECAST_DURATION + 1, 0.0);
        delete par;
        return dummy;
    }

    Community* community = build_community(par);
    if (vaccine) {
        double target_coverage  = coverage;
        double catchup_coverage = coverage;

        par->fVEH                    = 0.803; // fraction of hospitalized cases prevented by vaccine
        par->bVaccineLeaky           = true;
        par->numVaccineDoses         = 3;
        par->vaccineDoseInterval     = 182;
        par->vaccineBoosting         = false;
        if (catchup) {
            for (int catchup_age = target + 1; catchup_age <= catchup_to; catchup_age++) {
                par->vaccinationEvents.emplace_back(catchup_age, RESTART_BURNIN*365, catchup_coverage);
            }
        } 

        for (int vacc_year = RESTART_BURNIN; vacc_year < TOTAL_DURATION; vacc_year++) {
            par->vaccinationEvents.emplace_back(target, vacc_year*365, target_coverage);
        }
    }

    par->fVESs.clear();
    par->fVESs_NAIVE.clear();
    if (vaccine_mechanism == 0) {                                  // WHO "baseline" scenario: A2b + B2 + C3a
        // perfect efficacy that wanes rapidly
        par->fVESs             = vector<double>(NUM_OF_SEROTYPES, 1.0);
        par->fVESs_NAIVE       = vector<double>(NUM_OF_SEROTYPES, 1.0);
        par->whoDiseaseOutcome = INC_NUM_INFECTIONS;
        par->whoBreakthrough   = BREAKTHROUGH_SEROCONVERSION;
        par->linearlyWaningVaccine   = true;
        par->vaccineImmunityDuration = 2*365;
        par->whoWaning         = UNIVERSAL_WANING;
    } else if (vaccine_mechanism == 1 or vaccine_mechanism == 2) { // PLoS NTDs mechanisms
        par->fVESs             = {0.6, 0.54, 0.9, 0.95};
        par->fVESs_NAIVE       = {0.3, 0.27, 0.45, 0.48};
        par->whoDiseaseOutcome = VAC_ISNT_INFECTION;
        par->whoBreakthrough   = BREAKTHROUGH_SEROCONVERSION;
        par->whoWaning         = UNIVERSAL_WANING;
        if (vaccine_mechanism == 1) {
            par->linearlyWaningVaccine   = false;                  // main-text mechanism
        } else if (vaccine_mechanism == 2){
            par->linearlyWaningVaccine   = true;                   // a waning model, albeit faster than those in the SI
            par->vaccineImmunityDuration = 2*365;
        } else {
            cerr << "You shall not pass! (Unsupported vaccine mechanism: " << vaccine_mechanism << ")\n";
            exit(-153);
        }
    } else {
        cerr << "Unsupported vaccine mechanism: " << vaccine_mechanism << endl;
        exit(-152);
    }

    seed_epidemic(par, community);
    simulate_epidemic(par, community, process_id);

    string imm_filename = "/ufrc/longini/tjhladish/imm_baseline-merida/immunity." + process_id;
    write_immunity_file(par, community, process_id, imm_filename, par->nRunLength);

    //vector<int> epi_sizes = simulate_epidemic(par, community, process_id);

    time (&end);
    double dif = difftime (end,start);

    vector<long double> metrics = tally_counts(par, community);

    stringstream ss;
    ss << mp->mpi_rank << " end " << hex << process_id << " " << dec << dif << " ";

    for (auto i: args) ss << i << " ";
    for (auto i: metrics) ss << i << " ";
    ss << endl;

    string output = ss.str();
    fputs(output.c_str(), stderr);

    /*const int pop_size = community->getNumPerson();
    for (unsigned int i = 0; i < epi_sizes.size(); i++) { 
        // convert infections to cases per 1,000
        metrics[i] = ((long double) 1e3*epi_sizes[i])/pop_size; 
    }*/

    delete par;
    delete community;
    return metrics;
}

void usage() {
    cerr << "\n\tUsage: ./abc_sql abc_config_sql.json --process\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate -n <number of simulations per database write>\n\n";

}


int main(int argc, char* argv[]) {

    if (not (argc == 3 or argc == 5 or argc == 6) ) {
        usage();
        exit(100);
    }

    bool process_db = false;
    bool simulate_db = false;
    int buffer_size = -1;

    for (int i=2; i < argc;  i++ ) {
        if ( strcmp(argv[i], "--process") == 0  ) {
            process_db = true;
        } else if ( strcmp(argv[i], "--simulate") == 0  ) {
            simulate_db = true;
            buffer_size = buffer_size == -1 ? 1 : buffer_size;
        } else if ( strcmp(argv[i], "-n" ) == 0 ) {
            buffer_size = atoi(argv[++i]);
        } else {
            usage();
            exit(101);
        }
    }

    AbcSmc* abc = new AbcSmc();
    abc->parse_config(string(argv[1]));
    if (process_db) {
        gsl_rng_set(RNG, time(NULL) * getpid()); // seed the rng using sys time and the process id
        abc->process_database(RNG);
    }

    if (simulate_db) {
        time(&GLOBAL_START_TIME);
        abc->set_simulator(simulator);
        abc->simulate_next_particles(buffer_size);
    }

    return 0;
}
