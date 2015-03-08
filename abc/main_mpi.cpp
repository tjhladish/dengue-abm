#include "mpi.h"
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
    double _nmos     = args[3];
    
    double _betamp   = args[4]; // mp and pm and not separately
    double _betapm   = args[4]; // identifiable, so they're the same
    //double _betamp   = 0.25;
    //double _betapm   = 0.1;
    string HOME(std::getenv("HOME"));
    string pop_dir = HOME + "/work/dengue/pop-yucatan"; 
//    string WORK(std::getenv("WORK"));
//    string imm_dir = WORK + "/initial_immunity";
//    string imm_dir = pop_dir + "/immunity";

    //par->randomseed = 5500;
    par->abcVerbose = true;
    par->nRunLength = 155*365;
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
    par->expansionFactor = _EF;
    par->fMosquitoMove = _mos_move;
    par->mosquitoMoveModel = "weighted";
    par->fMosquitoTeleport = 0.0;
    par->nDefaultMosquitoCapacity = (int) _nmos;
    par->eMosquitoDistribution = EXPONENTIAL;

    par->nDaysImmune = 730;
    par->fVESs.clear();
    par->fVESs.resize(4, 0);

    par->simulateAnnualSerotypes = true;
    par->normalizeSerotypeIntros = true;
    if (par->simulateAnnualSerotypes) par->generateAnnualSerotypes();

    // 100 year burn-in, 20 years of no dengue, then re-introduction
    par->annualIntroductions = vector<double>(100, 1.0);
    par->annualIntroductions.resize(120, 0.0);
    par->annualIntroductions.resize(155, 1.0);

    par->populationFilename = pop_dir + "/population-yucatan.txt";
    par->immunityFilename   = "";
    par->locationFilename   = pop_dir + "/locations-yucatan.txt";
    par->networkFilename    = pop_dir + "/network-yucatan.txt";
    par->swapProbFilename   = pop_dir + "/swap_probabilities-yucatan.txt";
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


/*
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
}*/

unsigned int report_process_id (vector<long double> &args, const MPI_par* mp, const time_t start_time) {
    // CCRC32 checksum based on string version of argument values
    CCRC32 crc32;
    crc32.Initialize();

    string argstring;
    for (unsigned int i = 0; i < args.size(); i++) argstring += to_string((long double) args[i]) + " ";

    const unsigned char* argchars = reinterpret_cast<const unsigned char*> (argstring.c_str());
    const int len = argstring.length();
    const int process_id = crc32.FullCRC(argchars, len);
    double dif = difftime (start_time, GLOBAL_START_TIME);

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

vector<long double> simulator(vector<long double> args, const MPI_par* mp) {
    // initialize bookkeeping for run
    time_t start ,end;
    time (&start);
    const unsigned int process_id = report_process_id(args, mp, start);

    // initialize & run simulator 
    const Parameters* par = define_simulator_parameters(args); 
    Community* community = build_community(par);
    //seed_epidemic(par, community);
    double seropos_87 = 0.0;
    vector<int> serotested_ids = read_pop_ids("8-14_merida_ids.txt");
    vector<int> epi_sizes = simulate_epidemic(par, community, process_id, serotested_ids, seropos_87);
    vector<int>(epi_sizes.begin()+120, epi_sizes.end()).swap(epi_sizes); // throw out first 120 values

    time (&end);
    double dif = difftime (end,start);

    // calculate linear regression based on estimated reported cases
    const double ef = par->expansionFactor;
//    vector<double> x(epi_sizes.size());
//    vector<double> y(epi_sizes.size());
    Col y(epi_sizes.size());
    const int pop_size = community->getNumPerson();
    for (unsigned int i = 0; i < epi_sizes.size(); i++) { 
        // convert infections to cases per 100,000
        y[i] = ((float_type) 1e5*epi_sizes[i])/(ef*pop_size); 
        // year index, for regression
 //       x[i] = i+1.0;
    }
    //Fit* fit = lin_reg(x, y);

    stringstream ss;
    ss << mp->mpi_rank << " end " << hex << process_id << " " << dec << dif << " ";
    // parameters
    ss << par->expansionFactor << " " << par->fMosquitoMove << " " << par->annualIntroductionsCoef << " "
       << par->nDefaultMosquitoCapacity << " " << par->betaMP << " ";
    // metrics
    float_type _mean             = mean(y);
    float_type _median           = median(y);
    float_type _stdev            = sqrt(variance(y, _mean));
    float_type _max              = max(y);
    float_type _skewness         = skewness(y);
    float_type _median_crossings = median_crossings(y);
    float_type _seropos          = seropos_87;
    ss << _mean << " " << _median << " " << _stdev << " " << _max << " " << _skewness << " " << _median_crossings << " " << _seropos << endl;
      
//       << fit->m << " " << fit->b << " " << fit->rsq << endl;

    string output = ss.str();
    fputs(output.c_str(), stderr);
    
    vector<long double> metrics;
    append_if_finite(metrics, _mean);
    append_if_finite(metrics, _median);
    append_if_finite(metrics, _stdev);
    append_if_finite(metrics, _max);
    append_if_finite(metrics, _skewness);
    append_if_finite(metrics, _median_crossings);
    append_if_finite(metrics, _seropos);

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
    
    //const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);
    gsl_rng_set(RNG, time (NULL) * getpid()); // seed the rng using sys time and the process id

    AbcSmc* abc = new AbcSmc(mp);
    abc->set_simulator(simulator);
    abc->parse_config(string(argv[1]));
    time(&GLOBAL_START_TIME);
    abc->run(RNG);

    MPI_Finalize();
    return 0;
}
