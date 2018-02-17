#include "mpi.h"
#include "AbcSmc.h"
#include "mpi_simulator.h"
#include <cstdlib>
#include "CCRC32.h"
#include "Utility.h"
#include "ImmunityGenerator.h"

using namespace std;

using dengue::util::to_string;
using dengue::util::mean;
using dengue::util::stdev;
using dengue::util::max_element;

time_t GLOBAL_START_TIME;

void setup_mpi(MPI_par &m, int &argc, char **argv) {
    /* MPI variables */
    m.comm  = MPI_COMM_WORLD;
    m.info  = MPI_INFO_NULL;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(m.comm, &m.mpi_size);
    MPI_Comm_rank(m.comm, &m.mpi_rank);  
}

Parameters* define_simulator_parameters(vector<long double> args) {
    Parameters* par = new Parameters();
    par->define_defaults();

    double _EF       = args[0];
    double _mos_move = args[1];
    double _exp_coef = args[2];
    double _betamp   = args[3];
    double _betapm   = args[4];
    double _nmos     = args[5];
    
    string HOME(std::getenv("HOME"));
    string pop_dir = HOME + "/dengue/pop-yucatan"; 
//    string WORK(std::getenv("WORK"));
//    string imm_dir = WORK + "/initial_immunity";
//    string imm_dir = pop_dir + "/immunity";

    par->randomseed = 5500;
    par->nRunLength = 4500;
    par->annualIntroductionsCoef = pow(10,_exp_coef);

    par->nDailyExposed[0] = 1.0; 
    par->nDailyExposed[1] = 1.0;
    par->nDailyExposed[2] = 1.0;
    par->nDailyExposed[3] = 0.0;
    par->fPrimaryPathogenicity[0] = 1.0;
    par->fPrimaryPathogenicity[1] = 0.25;
    par->fPrimaryPathogenicity[2] = 1.0;
    par->fPrimaryPathogenicity[3] = 0.25;
    par->fSecondaryScaling[0] = 1.0;
    par->fSecondaryScaling[1] = 1.0;
    par->fSecondaryScaling[2] = 1.0;
    par->fSecondaryScaling[3] = 1.0;
    par->betaPM = _betapm;
    par->betaMP = _betamp;
    par->expansionFactor = _EF;
    par->fMosquitoMove = _mos_move;
    par->szMosquitoMoveModel = "weighted";
    par->fMosquitoTeleport = 0.0;
    par->nDefaultMosquitoCapacity = (int) _nmos;
    par->eMosquitoDistribution = EXPONENTIAL;

    {
        static const int _time_periods[] = {7,7,7,7,7,7,7,7,7,7,
                                          7,7,7,7,7,7,7,7,7,7,
                                          7,7,7,7,7,7,7,7,7,7,
                                          7,7,7,7,7,7,7,7,7,7,
                                          7,7,7,7,7,7,7,7,7,7,
                                          7,8};  
        static const double _multipliers[] = {0.05,0.04,0.05,0.04,0.03,0.04,0.05,0.03,0.02,0.02,
                                            0.03,0.03,0.05,0.04,0.05,0.04,0.06,0.07,0.08,0.09,
                                            0.11,0.15,0.18,0.19,0.24,0.28,0.38,0.46,0.45,0.61,
                                            0.75,0.97,0.91,1.00,0.94,0.85,0.79,0.71,0.65,0.65,
                                            0.42,0.30,0.26,0.27,0.11,0.10,0.11,0.12,0.09,0.08,
                                            0.04,0.07};
        int num_periods = sizeof(_multipliers) / sizeof(_multipliers[0]);
        //assert(_time_periods.size() == _multipliers.size());
        par->nSizeMosquitoMultipliers = num_periods;
        par->nMosquitoMultiplierCumulativeDays[0] = 0;

        for (unsigned int j=0; j<num_periods; j++) {
            par->nMosquitoMultiplierDays[j] = _time_periods[j];
            par->nMosquitoMultiplierCumulativeDays[j+1] =
                par->nMosquitoMultiplierCumulativeDays[j]+par->nMosquitoMultiplierDays[j];
            par->fMosquitoMultipliers[j] = _multipliers[j];
        }
    }

    par->nDaysImmune = 730;
    par->fVESs.clear();
    par->fVESs.resize(4, 0.7);

    par->annualIntroductions = {  400519,  // 2000 total cases in Americas, reported by PAHO
                                  652212,  // 2001
                                 1015420,  // 2002
                                  517617,  // 2003
                                  267050,  // 2004
                                  427627,  // 2005
                                  552141,  // 2006
                                  900782,  // 2007
                                  908926,  // 2008
                                 1134001,  // 2009
                                 1663276,  // 2010
                                 1093252,  // 2011
                                 1120902,  // 2012
                                 2386836}; // 2013

    par->szPopulationFile = pop_dir + "/population-yucatan.txt";
    par->szImmunityFile   = "";
    //par->szImmunityFile   = imm_dir + "/" + to_string((int) _EF) + ".txt"; // stupidest cast in the world.  what is wrong with icc?!
    par->szLocationFile   = pop_dir + "/locations-yucatan.txt";
    par->szNetworkFile    = pop_dir + "/network-yucatan.txt";
    par->szSwapProbFile   = pop_dir + "/swap_probabilities-yucatan.txt";

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


void sample_immune_history(Community* community, const Parameters* par) {
    const int last_immunity_year = 1999;
    vector<vector<vector<int>>> full_pop = simulate_immune_dynamics(par->expansionFactor, last_immunity_year);

    int N = community->getNumPerson();
    for(int i=0; i<N; i++) {
        Person* person = community->getPerson(i);

        int age = person->getAge();
        int maxAge = MAX_CENSUS_AGE;

        age = age <= maxAge ? age : maxAge;
        int r = gsl_rng_uniform_int(RNG, full_pop[age].size());
        vector<int> states = full_pop[age][r];
        // we need to go through the states, greatest to least
        // Serotype 0, 1,  2, 3   -->   1, 3, 2, 0
        //        { 0, 12, 1, 2 } --> { 1, 3, 2, 0 }
        vector<int> indices = ordered(states);
        for (int serotype: indices) {
            if (states[serotype]>0) {
                person->setImmunity((Serotype) serotype);
                person->initializeNewInfection();
                person->setRecoveryTime(-365*states[serotype]); // last dengue infection was x years ago
            }
        }
    }
}

unsigned int report_process_id (vector<long double> &args, const MPI_par* mp, const time_t start_time) {
    // CCRC32 checksum based on string version of argument values
    CCRC32 crc32;
    crc32.Initialize();

    string argstring;
    for (unsigned int i = 0; i < args.size(); i++) argstring += to_string((long double) args[i]) + " ";

    const unsigned char* argchars = reinterpret_cast<const unsigned char*> (argstring.c_str());
    const int len = argstring.length();
    const int process_id = crc32.FullCRC(argchars, len);
    //fprintf(stderr, "%Xbegin\n", process_id);

    double dif = difftime (start_time, GLOBAL_START_TIME);

    stringstream ss;
    ss << mp->mpi_rank << " begin " << hex << process_id << " " << dif << " " << argstring << endl;
    string output = ss.str();
    fprintf(stderr, output.c_str());

    return process_id;
}

void append_if_finite(vector<long double> &vec, double val) {
    if (isfinite(val)) { 
        vec.push_back((long double) val);
    } else {
        vec.push_back(0);
    }
}

// wrapper for simulator
// must take vector of doubles (ABC paramters) 
// and return vector of doubles (ABC metrics)
vector<long double> simulator(vector<long double> args, const MPI_par* mp) {
    // initialize bookkeeping for run
    time_t start ,end;
    time (&start);
    const unsigned int process_id = report_process_id(args, mp, start);

    // initialize & run simulator 
    const Parameters* par = define_simulator_parameters(args); 
    gsl_rng_set(RNG, par->randomseed);
    Community* community = build_community(par);
    sample_immune_history(community, par);
    //vector<int> initial_susceptibles = community->getNumSusceptible();
    seed_epidemic(par, community);
    vector<int> epi_sizes = simulate_epidemic(par, community, process_id);

    time (&end);
    double dif = difftime (end,start);

    // calculate linear regression based on estimated reported cases
    const double ef = par->expansionFactor;
    vector<double> x(epi_sizes.size());
    vector<double> y(epi_sizes.size());
    for (unsigned int i = 0; i < epi_sizes.size(); i++) { 
        // convert infections to cases
        y[i] = ((double) epi_sizes[i])/ef; 
        // year index, for regression
        x[i] = i+1.0;
    }
    Fit* fit = lin_reg(x, y);

    stringstream ss;
    ss << mp->mpi_rank << " end " << hex << process_id << " " << dec << dif << " ";
    // parameters
    ss << par->expansionFactor << " " << par->fMosquitoMove << " " << par->annualIntroductionsCoef << " "
       << par->betaMP << " " << par->betaPM << " " << par->nDefaultMosquitoCapacity << " ";
    // metrics
    ss << mean(y) << " " << stdev(y) << " " << max_element(y) << " "
       << fit->m << " " << fit->b << " " << fit->rsq << endl;

    string output = ss.str();
    fprintf(stderr, output.c_str());
    
    vector<long double> metrics;
    append_if_finite(metrics, mean(y) );
    append_if_finite(metrics, stdev(y) );
    append_if_finite(metrics, max_element(y) );
    append_if_finite(metrics, fit->m );
    append_if_finite(metrics, fit->b );
    append_if_finite(metrics, fit->rsq );

    delete par;
    delete community;
    return metrics;
}


int main(int argc, char* argv[]) {
    MPI_par mp;
    setup_mpi(mp, argc, argv);

    if (argc != 2) {
        cerr << "\n\tUsage: ./abc_mpi abc_config_file.json\n\n";
        return 100;
    }
    
    const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);
    gsl_rng_set(RNG, time (NULL) * getpid()); // seed the rng using sys time and the process id

    AbcSmc* abc = new AbcSmc(mp);
    abc->set_simulator(simulator);
    abc->parse_config(string(argv[1]));
    time(&GLOBAL_START_TIME);
    abc->run(RNG);

    MPI_Finalize();
    return 0;
}
