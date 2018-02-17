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
bool RESUME;

const unsigned int calculate_process_id(vector< long double> &args, string &argstring);

Parameters* define_simulator_parameters(vector<long double> args, const unsigned long int rng_seed) {
    Parameters* par = new Parameters();
    par->define_defaults();

    double _EF       = args[0];
    double _mos_move = args[1];
    double _exp_coef = args[2];
    double _nmos     = args[3];
    
    double _betamp   = args[4]; // mp and pm and not separately
    double _betapm   = args[4]; // identifiable, so they're the same
    string HOME(std::getenv("HOME"));
    string pop_dir = HOME + "/work/dengue/pop-toy"; 
    vector<long double> abc_args(&args[0], &args[5]);
    string argstring;
    const string process_id = to_string(calculate_process_id(abc_args, argstring));

    par->randomseed = rng_seed;
    par->abcVerbose = true;
    par->nRunLength = 5*365;
    par->startDayOfYear = 1;
    par->annualIntroductionsCoef = pow(10,_exp_coef);

    // pathogenicity values fitted in
    // Reich et al, Interactions between serotypes of dengue highlight epidemiological impact of cross-immunity, Interface, 2013
    // Normalized from Fc values in supplement table 2, available at
    // http://rsif.royalsocietypublishing.org/content/10/86/20130414/suppl/DC1
    par->fPrimaryPathogenicity = {1.000, 0.825, 0.833, 0.317};
    par->fSecondaryScaling = {1.0, 1.0, 1.0, 1.0};

    par->betaPM = _betapm;
    par->betaMP = _betamp;
    par->expansionFactor = _EF;
    par->fMosquitoMove = _mos_move;
    par->mosquitoMoveModel = "weighted";
    par->fMosquitoTeleport = 0.0;
    par->nDefaultMosquitoCapacity = (int) _nmos;
    par->eMosquitoDistribution = EXPONENTIAL;

    par->nDaysImmune = 730;
    par->fVESs.clear();
    par->fVESs.resize(4, 0);

    par->nDailyExposed = {{1,1,1,1}};
    //par->simulateAnnualSerotypes = true;
    //par->normalizeSerotypeIntros = true;
    //if (par->simulateAnnualSerotypes) par->generateAnnualSerotypes();

    par->annualIntroductions = {1.0};

    par->populationFilename = pop_dir + "/population-toy.txt";
    par->locationFilename   = pop_dir + "/locations-toy.txt";
    par->networkFilename    = pop_dir + "/network-toy.txt";
    par->swapProbFilename   = pop_dir + "/swap_probabilities-toy.txt";
    if (RESUME) {
        // if resuming this particular run
        par->immunityFilename = "./debug_immunity/immunity." + to_string(process_id);
    } else {
        // otherwise use a vanilla immunity file, or "" if fully susceptible
        par->immunityFilename = pop_dir + "/immunity-toy.txt";
    }

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


vector<int> tally_counts(const Parameters* par, Community* community) {
    vector< vector<int> > symptomatic = community->getNumNewlySymptomatic();
    //vector< vector<int> > infected    = community->getNumNewlyInfected();
    const int num_years = (int) par->nRunLength/365;
    vector<int> s_tally(num_years, 0);
    //vector<vector<int> > i_tally(NUM_OF_SEROTYPES, vector<int>(num_years+1, 0)); // +1 to handle run lengths of a non-integral number of years

    for (int t=0; t<par->nRunLength; t++) {
        const int y = t/365;
        for (int s=0; s<NUM_OF_SEROTYPES; s++) {
            s_tally[y] += symptomatic[s][t];
    //        i_tally[s][y] += infected[s][t];
        }
    }
    /*for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = discard_years; y<num_years; ++y) {
            metrics.push_back(s_tally[s][y]);
        }
    }*/
    /*for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = discard_years; y<num_years; ++y) {
            metrics.push_back(i_tally[s][y]);
        }
    }*/
    return s_tally;
}

vector<long double> simulator(vector<long double> args, const unsigned long int rng_seed, const MPI_par* mp) {
    gsl_rng_set(RNG, rng_seed); // seed the rng using sys time and the process id
    // initialize bookkeeping for run
    time_t start ,end;
    time (&start);
    const unsigned int process_id = report_process_id(args, mp, start);

    // initialize & run simulator 
    Parameters* par = define_simulator_parameters(args, rng_seed); 
    par->mosquitoFilename         = "./debug_mos/mos." + to_string(process_id);
    par->mosquitoLocationFilename = "./debug_mosloc/mosloc." + to_string(process_id);

    Community* community = build_community(par);
    //seed_epidemic(par, community);
    double seropos_87 = 0.0;
    vector<int> serotested_ids = read_pop_ids("8-14_merida_ids.txt");
    simulate_abc(par, community, process_id, serotested_ids, seropos_87);

    if (RESUME)  {
        community->loadMosquitoes(par->mosquitoLocationFilename, par->mosquitoFilename);
    } else {
        write_mosquito_location_data(community, par->mosquitoFilename, par->mosquitoLocationFilename);

        // We might want to write the immunity file if this is the real posterior
        string imm_filename = "./debug_immunity/immunity." + to_string(process_id);
        write_immunity_file(par, community, process_id, imm_filename, par->nRunLength);
    }

    vector<int> case_sizes = tally_counts(par, community);
    vector<long double> metrics(case_sizes.begin(), case_sizes.end());

    time (&end);
    double dif = difftime (end,start);

    stringstream ss;
    ss << hex << process_id << dec << " ";
    for (auto v: metrics) ss << v << " ";
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
    RESUME = false;

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
        } else if ( strcmp(argv[i], "--resume") == 0 ) {
            RESUME = true;
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
