#include "mpi.h"
#include <unistd.h>
#include "AbcSmc.h"
#include "mpi_simulator.h"
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


Parameters* define_default_parameters(const int years_simulated) {
    Parameters* par = new Parameters();
    par->define_defaults();

    double _EF       = 40;
    double _mos_move = 0.5;
    //double _exp_coef = 1;
    double _nmos     = 50;
    //double _nmos     = 60;
    //double _nmos     = 80;
    
    double _betamp   = 0.25;
    double _betapm   = 0.1;
    string HOME(std::getenv("HOME"));
    string pop_dir = HOME + "/work/dengue/pop-yucatan"; 
//    string WORK(std::getenv("WORK"));
//    string imm_dir = WORK + "/initial_immunity";
//    string imm_dir = pop_dir + "/immunity";

//    par->randomseed = 5500;
    par->nRunLength = years_simulated*365 + 100;
    par->annualIntroductionsCoef = 1;
    par->nDailyExposed = {{1.0, 1.0, 1.0, 1.0}};

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

        for (int j=0; j<num_periods; j++) {
            par->nMosquitoMultiplierDays[j] = _time_periods[j];
            par->nMosquitoMultiplierCumulativeDays[j+1] =
                par->nMosquitoMultiplierCumulativeDays[j]+par->nMosquitoMultiplierDays[j];
            par->fMosquitoMultipliers[j] = _multipliers[j];
        }
    }

    par->nDaysImmune = 365;

    //par->annualSerotypeFile = pop_dir + "/serotypes-yucatan.txt";
    //par->loadAnnualSerotypes(par->annualSerotypeFile);

    par->szPopulationFile = pop_dir + "/population-yucatan.txt";
    par->szImmunityFile   = "";
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


void generate_homogeneous_immune_history(Community* community, const Parameters* par, float attack_rate) {
    int N = community->getNumPerson();
    for(int i=0; i<N; i++) {
        Person* person = community->getPerson(i);

        int age = person->getAge();
        age = age <= MAX_CENSUS_AGE ? age : MAX_CENSUS_AGE;

        // determine order in which person will be exposed to serotypes
        int sero_order[NUM_OF_SEROTYPES];
        for (i = 0; i < NUM_OF_SEROTYPES; i++) { sero_order[i] = i; }
        gsl_ran_shuffle (RNG, sero_order, NUM_OF_SEROTYPES, sizeof (int));

        // don't allow multiple infections in one year
        unordered_set<int> infection_years;
        for (int s=0; s<NUM_OF_SEROTYPES; ++s) {
            Serotype serotype = (Serotype) sero_order[s];
            // sample year of life in which infection occurred
            // gsl_ran_geometric returns integers >= 1
            int year = gsl_ran_geometric(RNG, attack_rate)-1;
            // person is old enough to have been infected in 'year'
            // AND they were not infected with a different serotype
            // already in that year
            if (age > year and infection_years.count(year) == 0) {
                infection_years.insert(year);
                person->setImmunity((Serotype) serotype);
                person->initializeNewInfection();
                person->setRecoveryTime(-365*(age-year)); // last dengue infection was x years ago
            }
        }
    }
}

void sample_immune_history(Community* community, const Parameters* par) {
    const int last_immunity_year = 2002;
    vector<vector<vector<int>>> full_pop = simulate_immune_dynamics(par->expansionFactor, last_immunity_year);

    int N = community->getNumPerson();
    for(int i=0; i<N; i++) {
        Person* person = community->getPerson(i);

        int age = person->getAge();

        age = age <= MAX_CENSUS_AGE ? age : MAX_CENSUS_AGE;
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


vector<long double> tally_counts(const Parameters* par, Community* community) {
    vector< vector<int> > infected    = community->getNumNewlyInfected();
    vector< vector<int> > symptomatic = community->getNumNewlySymptomatic();
    const int num_years = (int) par->nRunLength/365;
    vector<vector<int> > i_tally(NUM_OF_SEROTYPES, vector<int>(num_years+1, 0)); // +1 to handle run lengths of a non-integral number of years
    vector<vector<int> > s_tally(NUM_OF_SEROTYPES, vector<int>(num_years+1, 0)); // (any extra fraction of a year will be discarded)

    vector<long double> metrics;
    for (int t=0; t<par->nRunLength; t++) {
        const int y = t/365;
        for (int s=0; s<NUM_OF_SEROTYPES; s++) {
            i_tally[s][y] += infected[s][t];
            s_tally[s][y] += symptomatic[s][t];
        }
    }
    // flatten data structures into the metrics vector
    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = 0; y<num_years; ++y) {
            metrics.push_back(i_tally[s][y]);
        }
    }
    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = 0; y<num_years; ++y) {
            metrics.push_back(s_tally[s][y]);
        }
    }
    return metrics;
}


vector<long double> simulator(vector<long double> args, const MPI_par* mp) {

    const int years_simulated = 40;
    const int burnin = 20; // e.g., 0 means start vaccinating in the first simulated year

    vector<int> target_ages = {2,6,10,14};
    bool vaccine      = (bool) args[0];
    bool retro        = (bool) args[1];
    bool catchup      = (bool) args[2];
    int target        = target_ages[(unsigned int) args[3]];
    bool full_catchup = (bool) args[4];
    bool all_mature   = (bool) args[5];

    const int max_catchup_age = full_catchup ? 100 : 46;

    Parameters* par = define_default_parameters(years_simulated); 
    if (vaccine) {
        par->bVaccineLeaky = true;

        par->fVESs.clear();
        par->fVESs = {0.7, 0.35, 0.7, 0.7};

        par->fVESs_NAIVE.clear();
        if (all_mature) {
            par->fVESs_NAIVE = par->fVESs;
        } else {
            par->fVESs_NAIVE = {0.35, 0.0, 0.35, 0.35};
        }

        if (catchup) {
            for (int catchup_age = target + 1; catchup_age <= max_catchup_age; catchup_age++) {
                par->nVaccinateYear.push_back(burnin);
                par->nVaccinateAge.push_back(catchup_age);
                par->fVaccinateFraction.push_back(0.7);
                par->nSizeVaccinate++;
            }
        } 

        for (int vacc_year = burnin; vacc_year < years_simulated; vacc_year++) {
            par->nVaccinateYear.push_back(vacc_year);
            par->nVaccinateAge.push_back(target);
            par->fVaccinateFraction.push_back(0.7);
            par->nSizeVaccinate++;
        }

        if (retro) par->bRetroactiveMatureVaccine = true;
    }


    //initialize bookkeeping for run
    time_t start ,end;
    time (&start);
    const unsigned int process_id = report_process_id(args, mp, start);

    // initialize & run simulator 
//    gsl_rng_set(RNG, par->randomseed);
    Community* community = build_community(par);
    //sample_immune_history(community, par); // use historical data
    // initialize population immunity to make burn-in go faster
    float approx_attack_rate = 0.05; // approximate probability someone will be infected in a given year
    generate_homogeneous_immune_history(community, par, approx_attack_rate);

    seed_epidemic(par, community);
    simulate_epidemic(par, community, process_id);
    //vector<int> epi_sizes = simulate_epidemic(par, community, process_id);

    time (&end);
    double dif = difftime (end,start);

    vector<long double> metrics = tally_counts(par, community);

    stringstream ss;
    ss << mp->mpi_rank << " end " << hex << process_id << " " << dec << dif << " ";

    for (auto i: args) ss << i << " ";
    for (auto i: metrics) ss << i << " ";
    ss << endl;

    string output = ss.str();
    fprintf(stderr, output.c_str());

    /*const int pop_size = community->getNumPerson();
    for (unsigned int i = 0; i < epi_sizes.size(); i++) { 
        // convert infections to cases per 1,000
        metrics[i] = ((long double) 1e3*epi_sizes[i])/pop_size; 
    }*/

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

//    const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);
    gsl_rng_set(RNG, time (NULL) * getpid()); // seed the rng using sys time and the process id

    AbcSmc* abc = new AbcSmc(mp);
    abc->set_simulator(simulator);
    abc->parse_config(string(argv[1]));
    time(&GLOBAL_START_TIME);
    abc->run(RNG);

    //caseFile << "vaccine retro catchup target year inf1 inf2 inf3 inf4 sym1 sym2 sym3 sym4" << endl;

    MPI_Finalize();
    return 0;
}
