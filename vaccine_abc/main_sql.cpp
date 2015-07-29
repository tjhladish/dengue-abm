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

unsigned int calculate_process_id(vector< long double> &args, string &argstring);

Parameters* define_simulator_parameters(vector<long double> args, const unsigned long int rng_seed) {
    Parameters* par = new Parameters();
    par->define_defaults();

    //double _caseEF   = args[0]; // don't think we're doing anything with this
    double _mos_move = args[1];
    double _exp_coef = args[2];
    double _nmos     = args[3];
    
    double _betamp   = args[4]; // these two aren't independently identifiable
    double _betapm   = args[4];
    string HOME(std::getenv("HOME"));
    string pop_dir = HOME + "/work/dengue/pop-yucatan"; 

    string imm_dir("/scratch/lfs/thladish/imm_posterior_files");
    string sero_dir("/scratch/lfs/thladish/annual_serotype_files");
    string mos_dir("/scratch/lfs/thladish/mosquito_files");
    string mosloc_dir("/scratch/lfs/thladish/mosquito_location_files");
    vector<long double> abc_args(&args[0], &args[5]);
    string argstring;
    string particle_id = to_string(calculate_process_id(abc_args, argstring));

    par->randomseed = rng_seed;
    par->abcVerbose = true;
    par->nRunLength = 100 + (20*365);
    par->startDayOfYear = 1;
    par->annualIntroductionsCoef = pow(10,_exp_coef);

    // pathogenicity values fitted in
    // Reich et al, Interactions between serotypes of dengue highlight epidemiological impact of cross-immunity, Interface, 2013
    // Normalized from Fc values in supplement table 2, available at
    // http://rsif.royalsocietypublishing.org/content/10/86/20130414/suppl/DC1
    par->fPrimaryPathogenicity[0] = 1.000;
    par->fPrimaryPathogenicity[1] = 0.825;
    par->fPrimaryPathogenicity[2] = 0.833;
    par->fPrimaryPathogenicity[3] = 0.317;

    par->fSecondaryScaling[0] = 1.0;
    par->fSecondaryScaling[1] = 1.0;
    par->fSecondaryScaling[2] = 1.0;
    par->fSecondaryScaling[3] = 1.0;
    par->betaPM = _betapm;
    par->betaMP = _betamp;
    //par->expansionFactor = _caseEF;
    par->fMosquitoMove = _mos_move;
    par->mosquitoMoveModel = "weighted";
    par->fMosquitoTeleport = 0.0;
    par->nDefaultMosquitoCapacity = (int) _nmos;
    par->eMosquitoDistribution = EXPONENTIAL;

    par->nDaysImmune = 730;
    par->fVESs.clear();
    par->fVESs.resize(4, 0);

    par->hospitalizedFraction = 0.25; // fraction of cases assumed to be hospitalized
    par->fVEH = 0.803;                // fraction of hospitalized cases prevented by vaccine

    par->simulateAnnualSerotypes = false;
    par->normalizeSerotypeIntros = false;
    if (par->simulateAnnualSerotypes) par->generateAnnualSerotypes();
    par->annualSerotypeFilename = sero_dir + "/annual_serotypes." + particle_id;
    par->loadAnnualSerotypes();

    // we need to throw out the serotype years that were already used during ABC
    int abc_duration = 155;
    vector<vector<float> >(par->nDailyExposed.begin()+abc_duration, par->nDailyExposed.end()).swap(par->nDailyExposed);

    par->annualIntroductions = vector<double>(1, 1.0);

    par->populationFilename       = pop_dir + "/population-yucatan.txt";
    par->immunityFilename         = imm_dir + "/immunity." + particle_id;
//cerr << "argstring: " << argstring << " " << par->immunityFilename << endl;
    par->locationFilename         = pop_dir + "/locations-yucatan.txt";
    par->networkFilename          = pop_dir + "/network-yucatan.txt";
    par->swapProbFilename         = pop_dir + "/swap_probabilities-yucatan.txt";
    par->mosquitoFilename         = mos_dir + "/mos." + particle_id;
    par->mosquitoLocationFilename = mosloc_dir + "/mosloc." + particle_id;

    par->monthlyOutput = true;

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


unsigned int calculate_process_id(vector< long double> &args, string &argstring) {
    // CCRC32 checksum based on string version of argument values
    CCRC32 crc32;
    crc32.Initialize();

    for (unsigned int i = 0; i < args.size(); i++) argstring += to_string((long double) args[i]) + " ";

    const unsigned char* argchars = reinterpret_cast<const unsigned char*> (argstring.c_str());
    const int len = argstring.length();
    const int process_id = crc32.FullCRC(argchars, len);
    return process_id;
}


unsigned int report_process_id (vector<long double> &args, const MPI_par* mp, const time_t start_time) {
    double dif = difftime (start_time, GLOBAL_START_TIME);

    string argstring;
    unsigned int process_id = calculate_process_id(args, argstring);

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
    vector<vector<int> > h_tally(NUM_OF_SEROTYPES, vector<int>(num_years+1, 0)); // +1 to handle run lengths of a non-integral number of years
    vector<vector<int> > s_tally(NUM_OF_SEROTYPES, vector<int>(num_years+1, 0)); // (any extra fraction of a year will be discarded)
    vector<vector<int> > i_tally(NUM_OF_SEROTYPES, vector<int>(num_years+1, 0));

    vector<long double> metrics;
    for (int t=discard_days; t<par->nRunLength; t++) {
        const int y = t/365;
        for (int s=0; s<NUM_OF_SEROTYPES; s++) {
            h_tally[s][y] += severe[s][t];
            s_tally[s][y] += symptomatic[s][t];
            i_tally[s][y] += infected[s][t];
        }
    }
    // flatten data structures into the metrics vector
    // this could be tightened up using the right stride
    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = 0; y<num_years; ++y) {
            metrics.push_back(h_tally[s][y]);
        }
    }
    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = 0; y<num_years; ++y) {
            metrics.push_back(s_tally[s][y]);
        }
    }
    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = 0; y<num_years; ++y) {
            metrics.push_back(i_tally[s][y]);
        }
    }
    return metrics;
}


vector<long double> simulator(vector<long double> args, const unsigned long int rng_seed, const MPI_par* mp) {
    gsl_rng_set(RNG, rng_seed); // seed the rng using sys time and the process id

    Parameters* par = define_simulator_parameters(args, rng_seed);

    int discard_days = 100;
    int years_simulated = 20;
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

    bool vaccine                  = (bool) args[5];
    bool catchup                  = (bool) args[6];
    int target                    = (int) args[7];
    int catchup_to                = (int) args[8];
    float vector_reduction        = vector_controls[(unsigned int) args[9]];
    int vaccine_duration          = vaccine_durations[(unsigned int) args[10]];
    bool boosting                 = (bool) args[11];

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
        double default_coverage = 0.7;
        int num_target  = community->ageIntervalSize(9,10); // default target
        int num_catchup = community->ageIntervalSize(10,31);// default catchup

        double target_coverage  = default_coverage*num_target/community->ageIntervalSize(target, target+1);
        double catchup_coverage = default_coverage*num_catchup/community->ageIntervalSize(target+1, catchup_to+1);

        target_coverage = target_coverage > 1.0 ? 1.0 : target_coverage;
        catchup_coverage = catchup_coverage > 1.0 ? 1.0 : catchup_coverage;

        par->bVaccineLeaky = true;

        par->fVESs.clear();
        par->fVESs = {0.6, 0.54, 0.9, 0.95};

        par->fVESs_NAIVE.clear();
        par->fVESs_NAIVE = {0.3, 0.27, 0.45, 0.48};

        if (catchup) {
            for (int catchup_age = target + 1; catchup_age <= catchup_to; catchup_age++) {
                par->nVaccinateYear.push_back(0);
                par->nVaccinateAge.push_back(catchup_age);
                par->fVaccinateFraction.push_back(catchup_coverage);  // imperial used 0.5
                par->nSizeVaccinate++;
            }
        } 

        for (int vacc_year = 0; vacc_year < years_simulated; vacc_year++) {
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
