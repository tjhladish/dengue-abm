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

time_t GLOBAL_START_TIME;

//const unsigned int calculate_process_id(vector< long double> &args, string &argstring);

Parameters* define_simulator_parameters(vector<long double> args, const unsigned long int rng_seed) {
    Parameters* par = new Parameters();
    par->define_defaults();

    double _mild_EF      = args[0];
    double _severe_EF    = args[1];
    //double _sec_severity = args[2]; // not used
    //double _exp_coef     = args[3]; // not used
    double _nmos         = args[4];
    double _betamp       = args[5]; // mp and pm and not separately
    double _betapm       = args[5]; // identifiable, so they're the same

    string HOME(std::getenv("HOME"));
    string pop_dir = HOME + "/work/dengue/pop-yucatan";

    vector<long double> abc_args(&args[0], &args[6]);
    string argstring;
    //const string process_id = to_string(calculate_process_id(abc_args, argstring));

    par->abcVerbose = true;
    par->nRunLength = 100;
    par->annualIntroductionsCoef = 0.0;
    par->randomseed = rng_seed;

    par->primaryPathogenicity = {1.000, 0.825, 0.833, 0.317};
    par->secondaryPathogenicityOddsRatio = {1.0, 1.0, 1.0, 1.0};
    par->reportedFraction = {0.0, 1.0/_mild_EF, 1.0/_severe_EF};
    par->betaPM = _betapm;
    par->betaMP = _betamp;
    par->fMosquitoMove = 0.15;
    par->mosquitoMoveModel = "weighted";
    par->fMosquitoTeleport = 0.0;
    par->nDefaultMosquitoCapacity = (int) _nmos;
    par->eMosquitoDistribution = EXPONENTIAL;

    par->nDaysImmune = 730;

    par->loadDailyEIP(pop_dir + "/seasonal_avg_eip.out");
    par->loadDailyMosquitoMultipliers(pop_dir + "/mosquito_seasonality.out");
    par->populationFilename = pop_dir + "/population-yucatan.txt";
    par->immunityFilename   = "";
    par->locationFilename   = pop_dir + "/locations-yucatan.txt";
    par->networkFilename    = pop_dir + "/network-yucatan.txt";
    par->swapProbFilename   = pop_dir + "/swap_probabilities-yucatan.txt";

    par->nInitialInfected = {1};
    par->bSecondaryTransmission = false;
    par->weeklyOutput = true;
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


vector<long double> simulator(vector<long double> args, const unsigned long int rng_seed, const MPI_par* mp) {
    time_t start ,end;
    time (&start);

    Parameters* par = define_simulator_parameters(args, rng_seed);
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
