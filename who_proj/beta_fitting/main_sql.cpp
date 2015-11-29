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

const int RESTART_BURNIN     =  0;
const int FORECAST_DURATION  = 100;

Parameters* define_simulator_parameters(vector<long double> args, const unsigned long int rng_seed) {
    Parameters* par = new Parameters();
    par->define_defaults();

    //const vector<float> beta_multipliers = {0.4, 0.8, 1.0, 1.2, 1.6};
    const vector<float> beta_multipliers = {0.52, 0.70, 0.92, 1.25, 2.25};

    double _mild_EF      = args[0];
    double _severe_EF    = args[1];
    double _sec_severity = args[2];
    double _exp_coef     = args[3];
    double _nmos         = args[4];
    double _betamp       = args[5] * beta_multipliers[(int) args[6]]; // mp and pm and not separately
    double _betapm       = args[5] * beta_multipliers[(int) args[6]]; // identifiable, so they're the same

    vector<long double> abc_args(&args[0], &args[6]);
    string argstring;
    const string process_id = to_string(calculate_process_id(abc_args, argstring));

    par->randomseed              = rng_seed;
    par->dailyOutput             = false;
    par->yearlyOutput            = true;
    par->abcVerbose              = true;
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

    par->primarySevereFraction    = vector<double>(NUM_OF_SEROTYPES, _sec_severity/20); // ala Neil
    par->secondarySevereFraction  = vector<double>(NUM_OF_SEROTYPES, _sec_severity);
    par->tertiarySevereFraction   = vector<double>(NUM_OF_SEROTYPES, 0.0);
    par->quaternarySevereFraction = vector<double>(NUM_OF_SEROTYPES, 0.0);

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

    par->simulateAnnualSerotypes = true;
    par->normalizeSerotypeIntros = true;
    if (par->simulateAnnualSerotypes) par->generateAnnualSerotypes();

    par->annualIntroductions = {1.0};

    // load daily EIP
    // EIPs calculated by ../raw_data/weather/calculate_daily_eip.R, based on reconstructed yucatan temps
    par->loadDailyEIP(pop_dir + "/seasonal_avg_eip.out");

    // we add the startDayOfYear offset because we will end up discarding that many values from the beginning
    // mosquitoMultipliers are indexed to start on Jan 1
    par->loadDailyMosquitoMultipliers(pop_dir + "/mosquito_seasonality.out", par->nRunLength + par->startDayOfYear);

    par->populationFilename       = pop_dir    + "/population-"         + SIM_POP + ".txt";
    par->immunityFilename         = "";
    par->locationFilename         = pop_dir    + "/locations-"          + SIM_POP + ".txt";
    par->networkFilename          = pop_dir    + "/network-"            + SIM_POP + ".txt";
    par->swapProbFilename         = pop_dir    + "/swap_probabilities-" + SIM_POP + ".txt";
    par->mosquitoFilename         = output_dir + "/mos_mer_who/mos."       + to_string(process_id);
    par->mosquitoLocationFilename = output_dir + "/mosloc_mer_who/mosloc." + to_string(process_id);

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


vector<int> read_pop_ids (string filename) {
cerr << filename << endl;
    ifstream is(filename);
    istream_iterator<double> start(is), end;
    vector<int> ids(start, end);
    return ids;
}


vector<long double> simulator(vector<long double> args, const unsigned long int rng_seed, const MPI_par* mp) {
    gsl_rng_set(RNG, rng_seed); // seed the rng using sys time and the process id

    //initialize bookkeeping for run
    time_t start ,end;
    time (&start);
    const unsigned int process_id = report_process_id(args, mp, start);

    Parameters* par = define_simulator_parameters(args, rng_seed);

    string sero_filename = "/scratch/lfs/thladish/sero_mer_who/annual_serotypes." + to_string(process_id);
    par->writeAnnualSerotypes(sero_filename);

//  Posterior pars
//  0 mildEF real,
//  1 severeEF real,
//  2 sec_sev real,
//  3 exp_coef real,
//  4 num_mos real,
//  5 beta real,
//
//  Pseudo pars
//  6 beta_multiplier real,

    Community* community = build_community(par);
    //community->loadMosquitoes(par->mosquitoLocationFilename, par->mosquitoFilename);

    vector<int> serotested_ids = read_pop_ids(pop_dir + "/9_merida_ids.txt");

    seed_epidemic(par, community);
    vector<long double> seropos_9yo = simulate_who_fitting(par, community, process_id, serotested_ids);

    // We might want to write the immunity and mosquito files
    string imm_filename = "/scratch/lfs/thladish/imm_mer_who/immunity." + to_string(process_id);
    write_immunity_file(par, community, process_id, imm_filename, par->nRunLength);
    write_mosquito_location_data(community, par->mosquitoFilename, par->mosquitoLocationFilename);

    time (&end);
    double dif = difftime (end,start);

    stringstream ss;
    ss << mp->mpi_rank << " end " << hex << process_id << " " << dec << dif << " ";

    for (auto i: args) ss << i << " ";
    for (auto i: seropos_9yo) ss << i << " ";
    ss << endl;

    string output = ss.str();
    fputs(output.c_str(), stderr);

    delete par;
    delete community;

    return seropos_9yo;
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
