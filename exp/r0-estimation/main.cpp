#include <unistd.h>
#include <AbcSmc/AbcSmc.h>
#include "simulator.h"
#include <cstdlib>
#include "CCRC32.h"
#include "Utility.h"
#include <unordered_set>

using namespace std;

using dengue::util::to_string;
using dengue::util::mean;
using dengue::util::stdev;
using dengue::util::max_element;
using ABC::float_type;

time_t GLOBAL_START_TIME;

//unsigned int calculate_process_id(vector<double> &args, string &argstring);
const string SIM_POP = "yucatan";

Parameters* define_simulator_parameters(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial) {
    Parameters* par = new Parameters();
    par->define_defaults();
    par->serial = serial;

    double _mild_RF      = args[0];
  //double _p95_mild_RF  = args[1]; // not needed in this scope
    double _severe_RF    = args[2];
    double _base_path    = args[3];
    double _sec_severity = args[4];
    double _pss_ratio    = args[5];
  //double _exp_coef     = args[6];
    double _nmos         = args[7];
  //double _eip          = args[8]; // for doing R0 sensitivity analysis
    double _betamp       = 0.25; // beta values from chao et al
    double _betapm       = 0.10; //

    par->reportedFraction = {0.0, _mild_RF, _severe_RF}; // no asymptomatic infections are reported

    string HOME(std::getenv("HOME"));
    string pop_dir = HOME + "/work/dengue/pop-" + SIM_POP;
    vector<double> abc_args(&args[0], &args[7]); // args[8] if passing in EIP for R0 sensitivity analysis
    string argstring;
    //const string process_id = to_string(calculate_process_id(abc_args, argstring));

    par->randomseed              = rng_seed;
    par->periodicOutput          = true;
    par->periodicOutputInterval  = 25;
    par->dailyOutput             = false;
    par->weeklyOutput            = false;
    par->monthlyOutput           = false;
    par->yearlyOutput            = false;
    par->abcVerbose              = true;
    par->nRunLength              = 100;
    par->annualIntroductionsCoef = 0.0;
    par->birthdayInterval        = 365;

    // pathogenicity values fitted in
    // Reich et al, Interactions between serotypes of dengue highlight epidemiological impact of cross-immunity, Interface, 2013
    // Normalized from Fc values in supplement table 2, available at
    // http://rsif.royalsocietypublishing.org/content/10/86/20130414/suppl/DC1
    par->defineSerotypeRelativeRisks();
    par->basePathogenicity = _base_path;
    par->primaryPathogenicityModel = ORIGINAL_LOGISTIC;
    par->postSecondaryRelativeRisk = 0.1;

    par->primarySevereFraction    = vector<double>(NUM_OF_SEROTYPES, _sec_severity*_pss_ratio);
    par->secondarySevereFraction  = vector<double>(NUM_OF_SEROTYPES, _sec_severity);
    par->tertiarySevereFraction   = vector<double>(NUM_OF_SEROTYPES, _sec_severity/5.0);
    par->quaternarySevereFraction = vector<double>(NUM_OF_SEROTYPES, _sec_severity/5.0);

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
    par->loadDailyEIP(pop_dir + "/seasonal_EIP_24hr.out");
    par->loadDailyMosquitoMultipliers(pop_dir + "/mosquito_seasonality.out"); // do not use 2-argument version for R0 estimation
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // USING FIXED EIP AND FIXED MOSQUITO POP FOR SENSITIVITY ANALYSIS
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // par->extrinsicIncubationPeriods.clear();
    // par->extrinsicIncubationPeriods.emplace_back(0, 365, _eip);
    // par->mosquitoMultipliers.clear();
    // par->mosquitoMultipliers.emplace_back(0, 365, 1.0);

    par->populationFilename       = pop_dir    + "/population-"         + SIM_POP + ".txt";
    par->immunityFilename         = "";
    par->locationFilename         = pop_dir    + "/locations-"          + SIM_POP + ".txt";
    par->networkFilename          = pop_dir    + "/network-"            + SIM_POP + ".txt";
    par->swapProbFilename         = pop_dir    + "/swap_probabilities-" + SIM_POP + ".txt";

    par->nInitialInfected = {1,0,0,0};
    par->bSecondaryTransmission = false;
    return par;
}

string calculate_process_id(vector<double> &args, string &argstring) {
    // CCRC32 checksum based on string version of argument values
    CCRC32 crc32;
    crc32.Initialize();

    for (unsigned int i = 0; i < args.size(); i++) argstring += to_string((double) args[i]) + " ";

    const unsigned char* argchars = reinterpret_cast<const unsigned char*> (argstring.c_str());
    const int len = argstring.length();
    const int process_id = crc32.FullCRC(argchars, len);

    return to_string(process_id);
}


string report_process_id (vector<double> &args, const unsigned long int serial, const ABC::MPI_par* mp, const time_t start_time) {
    double dif = difftime (start_time, GLOBAL_START_TIME);

    string argstring;
    const string process_id = calculate_process_id(args, argstring);

    cerr << "pid in report_process_id (num args = " << args.size() << "): " << process_id << endl;
    stringstream ss;
    ss << mp->mpi_rank << " begin " << process_id << " " << dec << serial << " " << dif << " " << argstring << endl;
    string output = ss.str();
    fputs(output.c_str(), stderr);

    return to_string(process_id);
}



double tally_counts(const Parameters* par, Community* community) {
    vector< vector<int> > infected    = community->getNumNewlyInfected();
/*
cerr << "numNewlyInfected:\n";
for (auto v: infected) {
    for (int t=0; t<par->nRunLength; t++) cerr << v[t] << " ";
}
cerr << endl;
*/
                               // we're estimating R-zero, so
    double metric = -1.0; // subtract one for patient zero
    for (int t=0; t<par->nRunLength; t++) {
        for (int s=0; s<NUM_OF_SEROTYPES; s++) {
            metric += infected[s][t];
        }
    }
    cout << metric << endl;
    return metric;
}


vector<double> simulator(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const ABC::MPI_par* mp) {
    gsl_rng_set(RNG, rng_seed); // seed the rng using sys time and the process id

    time_t start ,end;
    time (&start);

    const unsigned int realization = (int) args[8]; // args[9] if passing in EIP for R0 sensitivity analysis
    const string process_id = report_process_id(args, serial, mp, start) + "." + to_string(realization);

    Parameters* par = define_simulator_parameters(args, rng_seed, serial);
    Community* community = build_community(par);
    //const vector<int> MONTH_START = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};

    vector<double> metrics;

    //for (unsigned int month = 0; month < MONTH_START.size(); ++month) { // to measure R0 on first of each month only
    for (unsigned int day = 0; day < 365; ++day) {
    //for (unsigned int day = 0; day < 1; ++day) { // when running code for R0 sensitivity analysis
        par->startDayOfYear = day;

        seed_epidemic(par, community);
        simulate_epidemic(par, community, process_id);

        metrics.push_back( tally_counts(par, community) );

        community->reset(); // reset immune states, remove infected mosquitoes, etc.
    }

    time (&end);
    double dif = difftime (end,start);

    stringstream ss;
    ss << mp->mpi_rank << " end " << process_id << " " << dec << serial << " " << dif << " ";
    for (auto i: args) ss << i << " ";
    ss << "| ";
    for (auto i: metrics) ss << i << " ";
    ss << endl;

    string output = ss.str();
    fputs(output.c_str(), stderr);

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
