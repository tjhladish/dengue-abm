#include <unistd.h>
#include "AbcSmc.h"
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

//unsigned int calculate_process_id(vector< long double> &args, string &argstring);
const string SIM_POP = "yucatan";

Parameters* define_simulator_parameters(vector<long double> args, const unsigned long int rng_seed, const unsigned long int serial) {
    Parameters* par = new Parameters();
    par->define_defaults();
    par->serial = serial;

    double _mild_EF      = args[0];
    double _severe_EF    = args[1];
    double _base_path    = args[2];
    double _sec_severity = args[3];
    double _pss_ratio    = args[4];
    double _exp_coef     = args[5];
    double _nmos         = args[6];
    double _betamp       = 0.25; // beta values from chao et al
    double _betapm       = 0.10; //
    par->reportedFraction = {0.0, 1.0/_mild_EF, 1.0/_severe_EF}; // no asymptomatic infections are reported

    string HOME(std::getenv("HOME"));
    string pop_dir = HOME + "/work/dengue/pop-" + SIM_POP;
    vector<long double> abc_args(&args[0], &args[7]);
    string argstring;
    //const string process_id = to_string(calculate_process_id(abc_args, argstring));

    par->randomseed              = rng_seed;
    par->dailyOutput             = false;
    par->weeklyOutput            = true;
    par->monthlyOutput           = true;
    par->yearlyOutput            = true;
    par->abcVerbose              = true;
    par->nRunLength              = 100;
    par->annualIntroductionsCoef = 0.0;
    //if ( SIM_POP == "merida") {
        par->annualIntroductionsCoef = 0.0;//pow(10, _exp_coef);
    //} else if ( SIM_POP == "yucatan" ) {
    //    // assuming parameter fitting was done on merida population 
    //    par->annualIntroductionsCoef = pow(10, _exp_coef)*1819498.0/839660.0;
    //} else {
    //    cerr << "ERROR: Unknown simulation population (SIM_POP): " << SIM_POP << endl;
    //    exit(523);
    //}

    // pathogenicity values fitted in
    // Reich et al, Interactions between serotypes of dengue highlight epidemiological impact of cross-immunity, Interface, 2013
    // Normalized from Fc values in supplement table 2, available at
    // http://rsif.royalsocietypublishing.org/content/10/86/20130414/suppl/DC1
    par->defineSerotypeRelativeRisks();
    par->basePathogenicity = _base_path;
    par->primaryPathogenicityModel = ORIGINAL_LOGISTIC;
    //par->annualFlavivirusAttackRate = _flav_ar;
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
    par->populationFilename       = pop_dir    + "/population-"         + SIM_POP + ".txt";
    par->immunityFilename         = "";
    par->locationFilename         = pop_dir    + "/locations-"          + SIM_POP + ".txt";
    par->networkFilename          = pop_dir    + "/network-"            + SIM_POP + ".txt";
    par->swapProbFilename         = pop_dir    + "/swap_probabilities-" + SIM_POP + ".txt";

    par->nInitialInfected = {1,0,0,0};
    par->bSecondaryTransmission = false;
    return par;
}


long double tally_counts(const Parameters* par, Community* community) {
    vector< vector<int> > infected    = community->getNumNewlyInfected();
/*
cerr << "numNewlyInfected:\n";
for (auto v: infected) {
    for (int t=0; t<par->nRunLength; t++) cerr << v[t] << " ";
}
cerr << endl;
*/
                               // we're estimating R-zero, so
    long double metric = -1.0; // subtract one for patient zero
    for (int t=0; t<par->nRunLength; t++) {
        for (int s=0; s<NUM_OF_SEROTYPES; s++) {
            metric += infected[s][t];
        }
    }
    cout << metric << " ";
    return metric;
}


vector<long double> simulator(vector<long double> args, const unsigned long int rng_seed, const unsigned long int serial, const ABC::MPI_par* mp) {
    gsl_rng_set(RNG, rng_seed); // seed the rng using sys time and the process id

    time_t start ,end;
    time (&start);

    Parameters* par = define_simulator_parameters(args, rng_seed, serial);
    Community* community = build_community(par);
    //const vector<int> MONTH_START = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};

    vector<long double> metrics;

    //for (unsigned int month = 0; month < MONTH_START.size(); ++month) {
    for (unsigned int day = 0; day < 365; ++day) {
        par->startDayOfYear = day;

        seed_epidemic(par, community);
        simulate_epidemic(par, community);

        metrics.push_back( tally_counts(par, community) );

        community->reset(); // reset immune states, remove infected mosquitoes, etc.
    }

    time (&end);
    double dif = difftime (end,start);

    stringstream ss;
    ss << " end " << dif << " ";
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
