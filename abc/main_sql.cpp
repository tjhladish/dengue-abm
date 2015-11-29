#include <unistd.h>
#include "AbcSmc.h"
#include "simulator.h"
#include <cstdlib>
#include "CCRC32.h"
#include "Utility.h"
#include "ImmunityGenerator.h"
#include "yucatan_serotype_generator.h"
#include <unordered_set>

using namespace std;

using dengue::util::to_string;
using dengue::util::mean;
using dengue::util::stdev;
using dengue::util::max_element;

time_t GLOBAL_START_TIME;

const unsigned int calculate_process_id(vector< long double> &args, string &argstring);
const string SIM_POP = "merida";


const int FIRST_YEAR          = 1879;                                 // inclusive
const int FIRST_OBSERVED_YEAR = 1979;
const int LAST_YEAR           = 2013;                                 // inclusive
const int DDT_START           = 1956 - FIRST_YEAR;                    // simulator year 77, counting from year 1
const int DDT_DURATION        = 23;                                   // 23 years long
const int FITTED_DURATION     = LAST_YEAR - FIRST_OBSERVED_YEAR + 1;  // 35 for 1979-2013 inclusive

Parameters* define_simulator_parameters(vector<long double> args, const unsigned long int rng_seed) {
    Parameters* par = new Parameters();
    par->define_defaults();

    double _mild_EF      = args[0];
    double _severe_EF    = args[1];
    double _sec_severity = args[2];
    double _pss_ratio    = args[3];
    double _exp_coef     = args[4];
    double _nmos         = args[5];
    double _betamp       = args[6]; // mp and pm and not separately
    double _betapm       = args[6]; // identifiable, so they're the same

    string HOME(std::getenv("HOME"));
    string pop_dir = HOME + "/work/dengue/pop-" + SIM_POP;
    string output_dir = "/scratch/lfs/thladish";

    vector<long double> abc_args(&args[0], &args[6]);
    string argstring;
    const string process_id = to_string(calculate_process_id(abc_args, argstring));

    par->randomseed              = rng_seed;
    par->dailyOutput             = false;
    par->abcVerbose              = true;
    int runLengthYears           = DDT_START + DDT_DURATION + FITTED_DURATION;
    par->nRunLength              = runLengthYears*365;
    par->startDayOfYear          = 100;
    par->annualIntroductionsCoef = pow(10, _exp_coef);

    // pathogenicity values fitted in
    // Reich et al, Interactions between serotypes of dengue highlight epidemiological impact of cross-immunity, Interface, 2013
    // Normalized from Fc values in supplement table 2, available at
    // http://rsif.royalsocietypublishing.org/content/10/86/20130414/suppl/DC1
    par->primaryPathogenicity    = {1.000, 0.825, 0.833, 0.317};
    par->secondaryPathogenicity  = par->primaryPathogenicity;
    par->tertiaryPathogenicity   = {0,0,0,0};
    par->quaternaryPathogenicity = {0,0,0,0};
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
    par->fVESs.clear();
    par->fVESs.resize(NUM_OF_SEROTYPES, 0);

    par->nDailyExposed = generate_serotype_sequences(RNG, FIRST_YEAR, FIRST_OBSERVED_YEAR, LAST_YEAR, TRANS_AND_NORM);

    //par->simulateAnnualSerotypes = false;
    //par->normalizeSerotypeIntros = true;
    // generate some extra years of serotypes, for subsequent intervention modeling
    //if (par->simulateAnnualSerotypes) par->generateAnnualSerotypes(runLengthYears+50);
    // 77 year burn-in, 23 years of no dengue, then re-introduction
    // annualIntros is indexed in terms of simulator (not calendar) years
    par->annualIntroductions = vector<double>(DDT_START, 1.0);
    par->annualIntroductions.resize(DDT_START+DDT_DURATION, 0.1); // 90% reduction in intros for 20 years
    par->annualIntroductions.resize(runLengthYears, 1.0);

    // load daily EIP
    // EIPs calculated by ../raw_data/weather/calculate_daily_eip.R, based on reconstructed Merida temps
    par->loadDailyEIP(pop_dir + "/seasonal_avg_eip.out");

    // we add the startDayOfYear offset because we will end up discarding that many values from the beginning
    // mosquitoMultipliers are indexed to start on Jan 1
    par->loadDailyMosquitoMultipliers(pop_dir + "/mosquito_seasonality.out", par->nRunLength + par->startDayOfYear);
    assert(par->mosquitoMultipliers.size() >= (DDT_START+DDT_DURATION)*365);
    for (int i = DDT_START*365; i < (DDT_START+DDT_DURATION)*365; ++i) par->mosquitoMultipliers[i].value *= 0.23; // 77% reduction in mosquitoes

    par->populationFilename       = pop_dir    + "/population-"         + SIM_POP + ".txt";
    par->immunityFilename         = "";
    par->locationFilename         = pop_dir    + "/locations-"          + SIM_POP + ".txt";
    par->networkFilename          = pop_dir    + "/network-"            + SIM_POP + ".txt";
    par->swapProbFilename         = pop_dir    + "/swap_probabilities-" + SIM_POP + ".txt";
    par->mosquitoFilename         = output_dir + "/mos_tmp/mos."       + to_string(process_id);
    par->mosquitoLocationFilename = output_dir + "/mosloc_tmp/mosloc." + to_string(process_id);
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

vector<int> read_pop_ids (string filename) {
    ifstream is(filename);
    istream_iterator<double> start(is), end;
    vector<int> ids(start, end);
    return ids;
}

void tally_counts(const Parameters* par, Community* community, vector<int>& all_cases_annual, vector<int>& severe_cases_annual) {
    vector< vector<int> > all_cases_daily = community->getNumNewlySymptomatic();
    vector< vector<int> > severe_cases_daily    = community->getNumSevereCases();

    const int num_years = (int) par->nRunLength/365;

    all_cases_annual.clear();
    all_cases_annual.resize(num_years, 0);

    severe_cases_annual.clear();
    severe_cases_annual.resize(num_years, 0);

    for (int t=0; t<par->nRunLength; t++) {
        const int y = t/365;
        for (int s=0; s<NUM_OF_SEROTYPES; s++) {
            all_cases_annual[y]    += all_cases_daily[s][t];
            severe_cases_annual[y] += severe_cases_daily[s][t];
        }
    }

    return;
}

vector<long double> simulator(vector<long double> args, const unsigned long int rng_seed, const MPI_par* mp) {
    gsl_rng_set(RNG, rng_seed); // seed the rng using sys time and the process id
    // initialize bookkeeping for run
    time_t start ,end;
    time (&start);
    const unsigned int process_id = report_process_id(args, mp, start);

    // initialize & run simulator
    const Parameters* par = define_simulator_parameters(args, rng_seed);

//    string sero_filename = "/scratch/lfs/thladish/sero_tmp/annual_serotypes." + to_string(process_id);
//    par->writeAnnualSerotypes(sero_filename);

    gsl_rng_set(RNG, rng_seed);
    Community* community = build_community(par);
    //seed_epidemic(par, community);
    double seropos_87 = 0.0;
    vector<int> serotested_ids = read_pop_ids("../pop-" + SIM_POP + "/8-14_merida_ids.txt");
    simulate_abc(par, community, process_id, serotested_ids, seropos_87);

    // We might want to write the immunity and mosquito files if this is the real posterior
//    string imm_filename = "/scratch/lfs/thladish/imm_tmp/immunity." + to_string(process_id);
//    write_immunity_file(par, community, process_id, imm_filename, par->nRunLength);
//    write_mosquito_location_data(community, par->mosquitoFilename, par->mosquitoLocationFilename);

    vector<int> all_cases;
    vector<int> severe_cases;
    tally_counts(par, community, all_cases, severe_cases);

    // throw out everything before end of DDT period
    vector<int>(all_cases.begin()+DDT_START+DDT_DURATION, all_cases.end()).swap(all_cases);
    vector<int>(severe_cases.begin()+DDT_START+DDT_DURATION, severe_cases.end()).swap(severe_cases);

    const double total_severe = accumulate(severe_cases.begin(), severe_cases.end(), 0.0);
    const double total_cases  = accumulate(all_cases.begin(), all_cases.end(), 0.0);
    const double mean_severe  = total_severe / total_cases;

    time (&end);
    double dif = difftime (end,start);

    const double mild_reporting = par->reportedFraction[(int) MILD];
    const double severe_reporting = par->reportedFraction[(int) SEVERE];

    // convert all cases and severe fraction to reported cases
    Col reported(all_cases.size());
    const int pop_size = community->getNumPerson();
    for (unsigned int i = 0; i < all_cases.size(); i++) {
        const float_type mild_reported   = (all_cases[i] - severe_cases[i]) * mild_reporting;
        const float_type severe_reported = severe_cases[i] * severe_reporting;
        // convert to reported cases per 100,000
        reported[i] = 1e5 * (mild_reported + severe_reported) / pop_size;
    }

    stringstream ss;
    ss << mp->mpi_rank << " end " << hex << process_id << " " << dec << dif << " ";
    // parameters
    ss << 1.0/par->reportedFraction[(int) MILD] << " " 
       << 1.0/par->reportedFraction[(int) SEVERE] << " " 
       << par->secondarySevereFraction[0] << " " 
       << log(par->annualIntroductionsCoef)/log(10) << " "
       << par->nDefaultMosquitoCapacity << " " 
       << par->betaMP << " | ";

    // metrics
    float_type _mean             = mean(reported);
    float_type _min              = quantile(reported, 0.00);
    float_type _quant25          = quantile(reported, 0.25);
    float_type _median           = quantile(reported, 0.50);
    float_type _quant75          = quantile(reported, 0.75);
    float_type _max              = quantile(reported, 1.00);
    float_type _stdev            = sqrt(variance(reported, _mean));
    float_type _skewness         = skewness(reported);
    float_type _median_crossings = median_crossings(reported);
    float_type _seropos          = seropos_87;
    float_type _severe_prev      = mean_severe;

    // logistic regression requires x values that are, in this case, just sequential ints
    vector<double> x(all_cases.size()); for(unsigned int i = 0; i < x.size(); ++i) x[i] = i;

    LogisticFit* fit = logistic_reg(x, severe_cases, all_cases);
    float_type _beta0, _beta1;
    if (fit->status == GSL_SUCCESS) {
        _beta0 = fit->beta0;
        _beta1 = fit->beta1;
    } else {
        _beta0 = 100;
        _beta1 = 100;
    }

    ss << _mean << " "
       << _min << " "
       << _quant25 << " "
       << _median << " "
       << _quant75 << " "
       << _max << " "
       << _stdev << " "
       << _skewness << " "
       << _median_crossings << " "
       << _seropos << " "
       << _severe_prev << " "
       << _beta0 << " "
       << _beta1 << endl;

    string output = ss.str();
    fputs(output.c_str(), stderr);

    vector<long double> metrics;
    append_if_finite(metrics, _mean);
    append_if_finite(metrics, _min);
    append_if_finite(metrics, _quant25);
    append_if_finite(metrics, _median);
    append_if_finite(metrics, _quant75);
    append_if_finite(metrics, _max);
    append_if_finite(metrics, _stdev);
    append_if_finite(metrics, _skewness);
    append_if_finite(metrics, _median_crossings);
    append_if_finite(metrics, _seropos);
    append_if_finite(metrics, _severe_prev);
    append_if_finite(metrics, _beta0);
    append_if_finite(metrics, _beta1);

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
