#include "mpi.h"
#include "AbcSmc.h"
#include "mpi_simulator.h"
#include <cstdlib>
#include "CCRC32.h"
#include "Utility.h"

using namespace std;

template <typename T> inline T sum(vector<T> list) { T sum=0; for (unsigned int i=0; i<list.size(); i++) sum += list[i]; return sum;}

template <typename T> inline double mean(vector<T> list) { return (double) sum(list) / list.size(); }

template <typename T>
double variance(vector<T> & numbers) {
    double x = mean(numbers);
    double var_num = 0;
    int N = numbers.size();
    if (N == 1) return 0;
    for (int i=0; i<N; i++) var_num += pow(numbers[i] - x, 2);
    double var = var_num/(N-1);
    return var;
}

template <typename T>
double stdev(vector<T> & numbers) { return sqrt( variance(numbers) ); }

template <typename T> inline
T max_element(vector<T> list) {
    T element = list[0];
    for (unsigned int i = 0; i < list.size(); i++) {
        element = max(element, list[i]);
    }
    return element;
}


template <typename T>
inline std::string to_string (const T& t) {
    std::stringstream ss;
    ss << t;
    return ss.str();
}


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

    double _EF        = args[0];
    double _mos_move  = args[1];
    double _daily_exp = args[2];
    double _betamp    = args[3];
    double _betapm    = args[4];
    
    string HOME(std::getenv("HOME"));
    string pop_dir = HOME + "/dengue/pop-yucatan"; 
    string WORK(std::getenv("WORK"));
    string imm_dir = WORK + "/initial_immunity";

    par->randomseed = 5500;
    par->nRunLength = 4500;
    par->nDailyExposed[0] = _daily_exp;
    par->nDailyExposed[1] = _daily_exp;
    par->nDailyExposed[2] = _daily_exp;
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
    par->nDefaultMosquitoCapacity = 50;
    par->eMosquitoDistribution = CONSTANT;

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

    par->szPopulationFile = pop_dir + "/population-yucatan.txt";
    par->szImmunityFile   = imm_dir + "/" + to_string((int) _EF) + ".txt"; // stupidest cast in the world.  what is wrong with icc?!
    par->szLocationFile   = pop_dir + "/locations-yucatan.txt";
    par->szNetworkFile    = pop_dir + "/network-yucatan.txt";
    par->szSwapProbFile   = pop_dir + "/swap_probabilities-yucatan.txt";

    return par;
}

unsigned int report_process_id (vector<long double> &args) {
    // CCRC32 checksum based on string version of argument values
    CCRC32 crc32;
    crc32.Initialize();

    string argstring;
    for (unsigned int i = 0; i < args.size(); i++) argstring += to_string((long double) args[i]) + " ";

    const unsigned char* argchars = reinterpret_cast<const unsigned char*> (argstring.c_str());
    const int len = argstring.length();
    const int process_id = crc32.FullCRC(argchars, len);
    fprintf(stderr, "%Xbegin\n", process_id);
    return process_id;
}

// wrapper for simulator
// must take vector of doubles (ABC paramters) 
// and return vector of doubles (ABC metrics)
vector<long double> simulator(vector<long double> args) {
    // initialize bookkeeping for run
    const unsigned int proccess_id = report_process_id(args);
    time_t start ,end;
    time (&start);

    // initialize & run simulator 
    const Parameters* par = define_simulator_parameters(args); 
    gsl_rng_set(RNG, par->randomseed);
    Community* community = build_community(par);
    vector<int> initial_susceptibles = community->getNumSusceptible();
    seed_epidemic(par, community);
    vector<int> epi_sizes = simulate_epidemic(par, community);

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
    ss << hex << proccess_id << "end " << dec << dif << " ";
    // parameters
    ss << par->expansionFactor << " " << par->fMosquitoMove << " " << par->nDailyExposed[0] << " "
       << par->betaMP << " " << par->betaPM << " ";
    // metrics
    ss << mean(y) << " " << stdev(y) << " " << max_element(y) << " "
       << fit->m << " " << fit->b << " " << fit->rsq << endl;

    string output = ss.str();
    fprintf(stderr, output.c_str());
    
    //vector<long double> metrics = {(long double) mean(y), (long double) stdev(y), 
    //(long double) max_element(y), (long double) fit->m, (long double) fit->b, (long double) fit->rsq};
    vector<long double> metrics;  // I hate this compiler
    metrics.push_back( (long double) mean(y) );
    metrics.push_back( (long double) stdev(y) );
    metrics.push_back( (long double) max_element(y) );
    metrics.push_back( (long double) fit->m );
    metrics.push_back( (long double) fit->b );
    metrics.push_back( (long double) fit->rsq );

    return metrics;
}


int main(int argc, char* argv[]) {
    MPI_par mp;
    setup_mpi(mp, argc, argv);

    if (argc != 2) {
        cerr << "\n\tUsage: ./abc abc_config_file.json\n\n";
        return 100;
    }
    
    const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);
    gsl_rng_set(RNG, time (NULL) * getpid()); // seed the rng using sys time and the process id

    AbcSmc* abc = new AbcSmc(mp);
    abc->set_simulator(simulator);
    abc->parse_config(string(argv[1]));
    abc->run(RNG);

    MPI_Finalize();
    return 0;
}
