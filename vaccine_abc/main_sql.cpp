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

const int DDT_START          = 77;  // simulator year 77
const int DDT_DURATION       = 23;  // 23 years long
const int FITTED_DURATION    = 30;  // 35 for 1979-2013 inclusive
const int RESTART_BURNIN     =  5;
const int FORECAST_DURATION  = 20;

Parameters* define_simulator_parameters(vector<long double> args, const unsigned long int rng_seed) {
    Parameters* par = new Parameters();
    par->define_defaults();

    double _mild_EF      = args[0];
    double _severe_EF    = args[1];
    double _sec_severity = args[2];
    double _exp_coef     = args[3];
    double _nmos         = args[4];
    double _betamp       = args[5]; // mp and pm and not separately
    double _betapm       = args[5]; // identifiable, so they're the same

    string HOME(std::getenv("HOME"));
    string pop_dir(HOME + "/work/dengue/pop-yucatan");
    string output_dir("/scratch/lfs/thladish");
    string imm_dir(output_dir + "/imm_tmp-5");
    string sero_dir(output_dir + "/sero_tmp-5");
    string mos_dir(output_dir + "/mos_tmp-5");
    string mosloc_dir(output_dir + "/mosloc_tmp-5");
    vector<long double> abc_args(&args[0], &args[6]);
    string argstring;
    const string process_id = to_string(calculate_process_id(abc_args, argstring));

    par->randomseed              = rng_seed;
    par->dailyOutput             = true;
    par->abcVerbose              = true;
    const int runLengthYears     = RESTART_BURNIN + FORECAST_DURATION;
    par->nRunLength              = runLengthYears*365;
    par->startDayOfYear          = 100;
    par->annualIntroductionsCoef = pow(10, _exp_coef);

    // pathogenicity values fitted in
    // Reich et al, Interactions between serotypes of dengue highlight epidemiological impact of cross-immunity, Interface, 2013
    // Normalized from Fc values in supplement table 2, available at
    // http://rsif.royalsocietypublishing.org/content/10/86/20130414/suppl/DC1
    par->primaryPathogenicity = {1.000, 0.825, 0.833, 0.317};
    par->secondaryPathogenicityOddsRatio = {1.0, 1.0, 1.0, 1.0};
    par->reportedFraction = {0.0, 1.0/_mild_EF, 1.0/_severe_EF}; // no asymptomatic infections are reported
    par->primarySevereFraction.clear();
    par->primarySevereFraction.resize(NUM_OF_SEROTYPES, 0.0);
    par->secondarySevereFraction.clear();
    par->secondarySevereFraction.resize(NUM_OF_SEROTYPES, _sec_severity);

    par->betaPM = _betapm;
    par->betaMP = _betamp;
    par->fMosquitoMove = 0.15;
    par->mosquitoMoveModel = "weighted";
    par->fMosquitoTeleport = 0.0;
    par->nDefaultMosquitoCapacity = (int) _nmos;
    par->eMosquitoDistribution = EXPONENTIAL;

    par->nDaysImmune = 730;
    par->fVESs.clear();
    par->fVESs.resize(NUM_OF_SEROTYPES, 0);

    //par->hospitalizedFraction = 0.25; // fraction of cases assumed to be hospitalized
    //par->fVEH = 0.803;                // fraction of hospitalized cases prevented by vaccine

    par->simulateAnnualSerotypes = false;
    par->normalizeSerotypeIntros = false;
    if (par->simulateAnnualSerotypes) par->generateAnnualSerotypes();
    par->annualSerotypeFilename = sero_dir + "/annual_serotypes." + to_string(process_id);
    par->loadAnnualSerotypes();

    // we need to throw out the serotype years that were already used during ABC
    int abc_duration = DDT_START + DDT_DURATION + FITTED_DURATION - RESTART_BURNIN;
    vector<vector<float> >(par->nDailyExposed.begin()+abc_duration, par->nDailyExposed.end()).swap(par->nDailyExposed);

    par->annualIntroductions = {1.0};

    // load daily EIP
    // EIPs calculated by ../raw_data/weather/calculate_daily_eip.R, based on reconstructed yucatan temps
    par->loadDailyEIP(pop_dir + "/seasonal_avg_eip.out");

    // we add the startDayOfYear offset because we will end up discarding that many values from the beginning
    // mosquitoMultipliers are indexed to start on Jan 1
    par->loadDailyMosquitoMultipliers(pop_dir + "/mosquito_seasonality.out", par->nRunLength + par->startDayOfYear);

    par->populationFilename       = pop_dir    + "/population-yucatan.txt";
    par->immunityFilename         = imm_dir    + "/immunity." + to_string(process_id);
    par->locationFilename         = pop_dir    + "/locations-yucatan.txt";
    par->networkFilename          = pop_dir    + "/network-yucatan.txt";
    par->swapProbFilename         = pop_dir    + "/swap_probabilities-yucatan.txt";
    par->mosquitoFilename         = mos_dir    + "/mos."    + to_string(process_id);
    par->mosquitoLocationFilename = mosloc_dir + "/mosloc." + to_string(process_id);
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


const unsigned int report_process_id (vector<long double> &args, const MPI_par* mp, const time_t start_time) {
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


vector<long double> tally_counts(const Parameters* par, Community* community, const int discard_days) {
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                                  IMPORTANT:                                           //
    // Update dummy metrics vector in calling function if number of metrics is changed here! //
    ///////////////////////////////////////////////////////////////////////////////////////////
    vector< vector<int> > severe      = community->getNumSevereCases();
    vector< vector<int> > symptomatic = community->getNumNewlySymptomatic();
    vector< vector<int> > infected    = community->getNumNewlyInfected();
    const int num_years = (int) par->nRunLength/365;
    vector<vector<int> > severe_tally(NUM_OF_SEROTYPES, vector<int>(num_years+1, 0)); // +1 to handle run lengths of a non-integral number of years
    vector<vector<int> > symptomatic_tally(NUM_OF_SEROTYPES, vector<int>(num_years+1, 0)); // (any extra fraction of a year will be discarded)
    vector<vector<int> > infected_tally(NUM_OF_SEROTYPES, vector<int>(num_years+1, 0));

    vector<long double> metrics;
    for (int t=discard_days; t<par->nRunLength; t++) {
        // use epidemic years, instead of calendar years
        const int y = (t-discard_days)/365;
        for (int s=0; s<NUM_OF_SEROTYPES; s++) {
            severe_tally[s][y]      += severe[s][t];
            symptomatic_tally[s][y] += symptomatic[s][t];
            infected_tally[s][y]    += infected[s][t];
        }
    }

    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = 0; y<num_years; ++y) {
            metrics.push_back(infected_tally[s][y]);
        }
    }
    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = 0; y<num_years; ++y) {
            metrics.push_back(symptomatic_tally[s][y] - severe_tally[s][y]);
        }
    }
    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = 0; y<num_years; ++y) {
            metrics.push_back(severe_tally[s][y]);
        }
    }
    return metrics;
}


vector<long double> simulator(vector<long double> args, const unsigned long int rng_seed, const MPI_par* mp) {
    gsl_rng_set(RNG, rng_seed); // seed the rng using sys time and the process id

    Parameters* par = define_simulator_parameters(args, rng_seed);

    const int discard_days        = 0;
    const int vac_start_year      = 5;
    const int years_simulated     = 25;
    vector<float> vector_controls = {0.0, 0.1, 0.25, 0.5};
    vector<int> vaccine_durations = {-1, 4*365, 10*365, 20*365};

//  Posterior pars
//  0 caseEF real,
//  1 mos_mov real,
//  2 exp_coef real,
//  3 num_mos real,
//  4 beta real,
//
//  Pseudo pars
//  5 vac real,
//  6 catchup real,
//  7 target real,
//  8 catchup_to real,
//  9 vec_control real,
// 10 vac_waning real,
// 11 vac_boosting real

    bool vaccine                  = (bool) args[6];
    bool catchup                  = (bool) args[7];
    int target                    = (int) args[8];
    int catchup_to                = (int) args[9];
    float vector_reduction        = vector_controls[(unsigned int) args[10]];
    int vaccine_duration          = vaccine_durations[(unsigned int) args[11]];
    bool boosting                 = (bool) args[12];

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
        vector<long double> dummy(years_simulated*NUM_OF_SEROTYPES*3, 0.0);
        delete par;
        return dummy;
    }

    Community* community = build_community(par);
    community->loadMosquitoes(par->mosquitoLocationFilename, par->mosquitoFilename);

    if (vaccine) {
        double default_target_coverage = 0.8;
        double default_catchup_coverage = 0.6;
        int num_target  = community->ageIntervalSize(9,10); // default target
        int num_catchup = community->ageIntervalSize(10,31);// default catchup

        double target_coverage  = default_target_coverage*num_target/community->ageIntervalSize(target, target+1);
        double catchup_coverage = default_catchup_coverage*num_catchup/community->ageIntervalSize(target+1, catchup_to+1);

        target_coverage = target_coverage > 1.0 ? 1.0 : target_coverage;
        catchup_coverage = catchup_coverage > 1.0 ? 1.0 : catchup_coverage;

        par->bVaccineLeaky = true;

        par->fVESs.clear();
        par->fVESs = {0.6, 0.54, 0.9, 0.95};

        par->fVESs_NAIVE.clear();
        par->fVESs_NAIVE = {0.3, 0.27, 0.45, 0.48};

        if (catchup) {
            for (int catchup_age = target + 1; catchup_age <= catchup_to; catchup_age++) {
                par->nVaccinateYear.push_back(vac_start_year);
                par->nVaccinateAge.push_back(catchup_age);
                par->fVaccinateFraction.push_back(catchup_coverage);  // imperial used 0.5
                par->nSizeVaccinate++;
            }
        } 

        for (int vacc_year = vac_start_year; vacc_year < years_simulated; vacc_year++) {
            par->nVaccinateYear.push_back(vacc_year);
            par->nVaccinateAge.push_back(target);
            par->fVaccinateFraction.push_back(target_coverage);      // imperial used 0.9
            par->nSizeVaccinate++;
        }
    }

/*
    const vector<float> MOSQUITO_MULTIPLIER_DEFAULTS = {0.179,0.128,0.123,0.0956,0.195,0.777,0.940,0.901,1.0,0.491,0.301,0.199};
    mosquitoMultipliers.clear();
    mosquitoMultipliers.resize(DAYS_IN_MONTH.size());
    int running_sum = 0;
    for (unsigned int j=0; j<mosquitoMultipliers.size(); j++) {
        mosquitoMultipliers[j].start = running_sum;
        mosquitoMultipliers[j].duration = DAYS_IN_MONTH[j];
        mosquitoMultipliers[j].value = MOSQUITO_MULTIPLIER_DEFAULTS[j];
        running_sum += DAYS_IN_MONTH[j];
    }
*/

    if (vector_reduction > 0.0) {
        //cerr << "ERROR: vector reduction not supported yet.  need to work out how 100 day + 20 yr simulation should work\n";
        //exit(-1834);
        assert( vector_reduction <= 1.0); 
        assert( par->mosquitoMultipliers.size() == 12); // expecting values for 12 months
        for (unsigned int i = 0; i < par->mosquitoMultipliers.size(); ++i) {
            float vector_coeff = 1.0 - vector_reduction; // normally, there's no reduction in vectors
            par->mosquitoMultipliers[i].value *= vector_coeff;
        }
    }

    //initialize bookkeeping for run
    time_t start ,end;
    time (&start);
    const unsigned int process_id = report_process_id(args, mp, start);

    seed_epidemic(par, community);
    simulate_epidemic(par, community, process_id);
    //vector<int> epi_sizes = simulate_epidemic(par, community, process_id);

    time (&end);
    double dif = difftime (end,start);

    vector<long double> metrics = tally_counts(par, community, discard_days);

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
