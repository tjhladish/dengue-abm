#include <unistd.h>
#include "AbcSmc.h"
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

const unsigned int calculate_process_id(vector< long double> &args, string &argstring);
const string SIM_POP = "merida";
const string HOME(std::getenv("HOME"));
const string pop_dir = HOME + "/work/dengue/pop-" + SIM_POP;
const string output_dir("/scratch/lfs/thladish");

const int RESTART_BURNIN    = 50;
const int FORECAST_DURATION = 30;
const int TOTAL_DURATION    = RESTART_BURNIN + FORECAST_DURATION;

Parameters* define_simulator_parameters(vector<long double> args, const unsigned long int rng_seed) {
    Parameters* par = new Parameters();
    par->define_defaults();

//  0 mildEF real,
//  1 severeEF real,
//  2 sec_sev real,
//  3 pss_ratio real
//  4 exp_coef real,
//  5 num_mos real,
//  6 beta real,
//
//  Pseudo pars
//  7 beta_multiplier real,

    //const vector<float> beta_multipliers = {0.52, 0.70, 0.92, 1.25, 2.25};
    const vector<float> beta_multipliers = {0.48, 0.725, 0.975, 1.375, 2.8};

    double _mild_EF      = args[0];
    double _severe_EF    = args[1];
    double _sec_severity = args[2];
    double _pss_ratio    = args[3];
    double _exp_coef     = args[4];
    double _nmos         = args[5];
    double _betamp       = args[6] * beta_multipliers[(int) args[7]]; // mp and pm and not separately
    double _betapm       = args[6] * beta_multipliers[(int) args[7]]; // identifiable, so they're the same

    vector<long double> abc_args(&args[0], &args[7]);
    string argstring;
    const string process_id = to_string(calculate_process_id(abc_args, argstring));

    par->randomseed              = rng_seed;
    par->dailyOutput             = false;
    par->monthlyOutput           = false;
    par->yearlyOutput            = true;
    par->abcVerbose              = false; // needs to be false to get WHO daily output
    const int runLengthYears     = RESTART_BURNIN + FORECAST_DURATION;
    par->nRunLength              = runLengthYears*365;
    par->startDayOfYear          = 100;
    par->annualIntroductionsCoef = pow(10, _exp_coef);

    // pathogenicity values fitted in
    // Reich et al, Interactions between serotypes of dengue highlight epidemiological impact of cross-immunity, Interface, 2013
    // Normalized from Fc values in supplement table 2, available at
    // http://rsif.royalsocietypublishing.org/content/10/86/20130414/suppl/DC1
    par->primaryPathogenicity    = {1.000, 0.825, 0.833, 0.317};
    par->secondaryPathogenicity  = {1.000, 0.825, 0.833, 0.317};
    par->tertiaryPathogenicity   = vector<double>(NUM_OF_SEROTYPES, 0.0);
    par->quaternaryPathogenicity = vector<double>(NUM_OF_SEROTYPES, 0.0);
    par->reportedFraction = {0.0, 1.0/_mild_EF, 1.0/_severe_EF}; // no asymptomatic infections are reported

    par->primarySevereFraction.clear();
    par->primarySevereFraction.resize(NUM_OF_SEROTYPES, _sec_severity*_pss_ratio);
    par->secondarySevereFraction.clear();
    par->secondarySevereFraction.resize(NUM_OF_SEROTYPES, _sec_severity);
    par->tertiarySevereFraction   = {0,0,0,0};
    par->quaternarySevereFraction = {0,0,0,0};

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
    //par->fVEH = 0.803;                // fraction of hospitalized cases prevented by vaccine

    par->simulateAnnualSerotypes = true;
    par->normalizeSerotypeIntros = true;
    if (par->simulateAnnualSerotypes) par->generateAnnualSerotypes();

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


const unsigned int calculate_process_id(vector< long double> &args, string &argstring) {
    // CCRC32 checksum based on string version of argument values
    CCRC32 crc32;
    crc32.Initialize();

    for (unsigned int i = 0; i < args.size(); i++) argstring += to_string((long double) args[i]) + " ";

    const unsigned char* argchars = reinterpret_cast<const unsigned char*> (argstring.c_str());
    const int len = argstring.length();
    const int process_id = crc32.FullCRC(argchars, len);

    return process_id;
}


const unsigned int report_process_id (vector<long double> &args, const ABC::MPI_par* mp, const time_t start_time) {
    double dif = difftime (start_time, GLOBAL_START_TIME);

    string argstring;
    const unsigned int process_id = calculate_process_id(args, argstring);

    stringstream ss;
    ss << mp->mpi_rank << " begin " << hex << process_id << " " << dif << " " << argstring << endl;
    string output = ss.str();
    fputs(output.c_str(), stderr);

    return process_id;
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
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                                  IMPORTANT:                                           //
    // Update dummy metrics vector in calling function if number of metrics is changed here! //
    ///////////////////////////////////////////////////////////////////////////////////////////
    //vector< vector<int> > severe      = community->getNumSevereCases();
    //vector< vector<int> > symptomatic = community->getNumNewlySymptomatic();
    vector< vector<int> > infected    = community->getNumNewlyInfected();
    const int num_years = FORECAST_DURATION;
    //vector<vector<int> > h_tally(NUM_OF_SEROTYPES, vector<int>(num_years+1, 0)); // +1 to handle run lengths of a non-integral number of years
    //vector<vector<int> > s_tally(NUM_OF_SEROTYPES, vector<int>(num_years+1, 0)); // (any extra fraction of a year will be discarded)
    vector<vector<int> > i_tally(NUM_OF_SEROTYPES, vector<int>(num_years+1, 0));

    vector<long double> metrics;
    for (int t=discard_days; t<par->nRunLength; t++) {
        // use epidemic years, instead of calendar years
        const int y = (t-discard_days)/365;
        for (int s=0; s<NUM_OF_SEROTYPES; s++) {
            //h_tally[s][y] += severe[s][t];
            //s_tally[s][y] += symptomatic[s][t];
            i_tally[s][y] += infected[s][t];
        }
    }
    // flatten data structures into the metrics vector
    // this could be tightened up using the right stride
//    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
//        for (int y = 0; y<num_years; ++y) {
//            metrics.push_back(h_tally[s][y]);
//        }
//    }
//    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
//        for (int y = 0; y<num_years; ++y) {
//            metrics.push_back(s_tally[s][y]);
//        }
//    }
    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = 0; y<num_years; ++y) {
            metrics.push_back(i_tally[s][y]);
        }
    }
    return metrics;
}


vector<long double> simulator(vector<long double> args, const unsigned long int rng_seed, const ABC::MPI_par* mp) {
    gsl_rng_set(RNG, rng_seed); // seed the rng using sys time and the process id

    //initialize bookkeeping for run
    time_t start, end;
    time (&start);
    const unsigned int process_id = report_process_id(args, mp, start);

    cerr << "SCENARIO " << rng_seed;
    for (auto _p: args) cerr << " " << _p; cerr << endl;

    Parameters* par = define_simulator_parameters(args, rng_seed);

    //vector<int> vaccine_durations = {-1, 10*365};
//  Posterior pars
//  0 mildEF real,
//  1 severeEF real,
//  2 sec_sev real,
//  3 pss_ratio real
//  4 exp_coef real,
//  5 num_mos real,
//  6 beta real,
//
//  Pseudo pars
//  7 beta_multiplier real,
//  8 catchup int,
//  9 vaccine_mechanism int,
//  10 vaccine int,
//  11 target int,
//  12 catchup_to int,
//  13 coverage int 

    vector<int> target_ages = {9, 16};        // default 9
    vector<int> catchup_ages = {17, 30};      // default 17
    vector<double> coverage_levels = {0.5, 0.8}; // default 0.8

    bool catchup           = (bool) args[8];
    int vaccine_mechanism  = (int) args[9];

    bool vaccine           = (bool) args[10];
    int vaccine_duration   = 2*365; 
    bool boosting          = false;
    int target             = target_ages[(int) args[11]];
    int catchup_to         = catchup_ages[(int) args[12]];
    double coverage        = coverage_levels[(int) args[13]];
    par->fVESs       = {0.6, 0.54, 0.9, 0.95};
    par->fVESs_NAIVE = {0.3, 0.27, 0.45, 0.48};
    // previous numbers used for WHO brazil:
    //par->fVESs       = {0.592, 0.498, 0.871, 0.914};
    //par->fVESs_NAIVE = {0.296, 0.249, 0.435, 0.457};

    if (vaccine_duration == -1) {
        par->linearlyWaningVaccine = false;
        par->vaccineImmunityDuration = INT_MAX;
    } else {
        par->linearlyWaningVaccine = true;
        par->vaccineImmunityDuration = vaccine_duration;
        par->vaccineBoosting = boosting;
    }

    bool nonsensical_parameters = false;
    // only run a non-vaccination campaign if all the vaccine parameters are 0
    // TODO - this should be reworked with new "catchup-to" parameter
    if (not vaccine and (catchup or boosting)) { nonsensical_parameters = true; }
    if ((vaccine_duration == -1) and boosting) { nonsensical_parameters = true; }
    if (nonsensical_parameters) {
        // 3 is because vaccinated cases, total cases, and infections are reported
        vector<long double> dummy(TOTAL_DURATION*NUM_OF_SEROTYPES, 0.0);
        delete par;
        return dummy;
    }

    Community* community = build_community(par);
    //community->loadMosquitoes(par->mosquitoLocationFilename, par->mosquitoFilename);

    if (vaccine) {
        double target_coverage  = coverage;
        double catchup_coverage = coverage;

        par->bVaccineLeaky = true;

        if (catchup) {
            for (int catchup_age = target + 1; catchup_age <= catchup_to; catchup_age++) {
                par->nVaccinateYear.push_back(RESTART_BURNIN);
                par->nVaccinateAge.push_back(catchup_age);
                par->fVaccinateFraction.push_back(catchup_coverage);
                par->nSizeVaccinate++;
            }
        } 

        for (int vacc_year = RESTART_BURNIN; vacc_year < TOTAL_DURATION; vacc_year++) {
            par->nVaccinateYear.push_back(vacc_year);
            par->nVaccinateAge.push_back(target);
            par->fVaccinateFraction.push_back(target_coverage);
            par->nSizeVaccinate++;
        }
    }

    if (vaccine_mechanism == 0) {        // "baseline" scenario: A2b + B2 + C3a
        par->whoDiseaseOutcome = INC_INFECTIONS_NAIVE;
        par->whoBreakthrough   = BREAKTHROUGH_SEROCONVERSION;
        par->whoWaning         = NAIVE_WANING_ONLY;
    } else if (vaccine_mechanism == 1) { // also want A1 + B2 + C3a
        par->whoDiseaseOutcome = VAC_ISNT_INFECTION;
        par->whoBreakthrough   = BREAKTHROUGH_SEROCONVERSION;
        par->whoWaning         = NAIVE_WANING_ONLY;
    } else {
        cerr << "Unsupported vaccine mechanism: " << vaccine_mechanism << endl;
        exit(-152);
    }

    seed_epidemic(par, community);
    simulate_epidemic(par, community, process_id);
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
