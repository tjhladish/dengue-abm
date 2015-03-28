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

Parameters* define_default_parameters(const int years_simulated) {
    Parameters* par = new Parameters();
    par->define_defaults();

    double _EF       = 38;
    double _mos_move = 0.45;
    double _exp_coef = -0.763;
    
    double _betamp   = 0.14;
    double _betapm   = 0.14;
    string HOME(std::getenv("HOME"));
    string pop_dir = HOME + "/work/dengue/pop-yucatan"; 
//    string WORK(std::getenv("WORK"));
//    string imm_dir = WORK + "/initial_immunity";
//    string imm_dir = pop_dir + "/immunity";

    par->abcVerbose = true;
    par->nRunLength = years_simulated*365;
    par->startDayOfYear = 100;
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
    par->eMosquitoDistribution = CONSTANT;

    par->nDaysImmune = 730;

    par->simulateAnnualSerotypes = true;
    par->normalizeSerotypeIntros = true;
    if (par->simulateAnnualSerotypes) par->generateAnnualSerotypes();

//    par->dailyEIPfilename   = pop_dir + "/merida_eip.out";
//    par->loadDailyEIP( par->dailyEIPfilename );
//    int period_length = 365;
//    par->extrinsicIncubationPeriods = shuffle_periods( RNG, par->extrinsicIncubationPeriods, period_length );

    par->populationFilename = pop_dir + "/population-yucatan.txt";
    par->immunityFilename   = "";
    par->locationFilename   = pop_dir + "/locations-yucatan.txt";
    par->networkFilename    = pop_dir + "/network-yucatan.txt";
    par->swapProbFilename   = pop_dir + "/swap_probabilities-yucatan.txt";

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


void generate_homogeneous_immune_history(Community* community, float attack_rate) {
    assert(attack_rate <= 1);
    assert(attack_rate >= 0);
    if (attack_rate == 0) return;

    int N = community->getNumPerson();
    for(int i=0; i<N; i++) {
        Person* person = community->getPerson(i);

        int age = person->getAge();
        age = age <= MAX_CENSUS_AGE ? age : MAX_CENSUS_AGE;

        // determine order in which person will be exposed to serotypes
        int sero_order[NUM_OF_SEROTYPES];
        for (int j = 0; j < NUM_OF_SEROTYPES; j++) { sero_order[j] = j; }
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
            //
            // gsl_rng_gemetric can return -1, for example if attack_rate == 0
            if (age > year and infection_years.count(year) == 0 and year > -1) { 
                infection_years.insert(year);
                // simulation should start in interepidemic period, so "historical" infections should be
                // offset by ~6 months
                int infectionTime = -365*(age-year)+182; // last dengue infection was x.5 years ago
                person->infect((Serotype) serotype, infectionTime);
            }
        }
    }
}

/*void sample_immune_history(Community* community, const Parameters* par) {
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
                // Ideally, simulate_immune_dynamics() would return days ago, instead of years ago,
                // so that we could handle infection times more consistently
                person->infect((Serotype) serotype, -365*states[serotype]+182);
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
    ss << mp->mpi_rank << " BEGIN " << hex << process_id << " " << dif << " " << argstring << endl;
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


vector<long double> tally_counts(const Parameters* par, Community* community, const int discard_years) {
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                                  IMPORTANT:                                           //
    // Update dummy metrics vector in calling function if number of metrics is changed here! //
    ///////////////////////////////////////////////////////////////////////////////////////////
    vector< vector<int> > vac_symptomatic = community->getNumVaccinatedCases();
    vector< vector<int> > symptomatic = community->getNumNewlySymptomatic();
    vector< vector<int> > infected    = community->getNumNewlyInfected();
    const int num_years = (int) par->nRunLength/365;
    vector<vector<int> > vc_tally(NUM_OF_SEROTYPES, vector<int>(num_years+1, 0)); // +1 to handle run lengths of a non-integral number of years
    vector<vector<int> > s_tally(NUM_OF_SEROTYPES, vector<int>(num_years+1, 0)); // (any extra fraction of a year will be discarded)
    vector<vector<int> > i_tally(NUM_OF_SEROTYPES, vector<int>(num_years+1, 0)); // +1 to handle run lengths of a non-integral number of years

    vector<long double> metrics;
    for (int t=0; t<par->nRunLength; t++) {
        const int y = t/365;
        for (int s=0; s<NUM_OF_SEROTYPES; s++) {
            vc_tally[s][y] += vac_symptomatic[s][t];
            s_tally[s][y] += symptomatic[s][t];
            i_tally[s][y] += infected[s][t];
        }
    }
    // flatten data structures into the metrics vector
    // this could be tightened up using the right stride
    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = discard_years; y<num_years; ++y) {
            metrics.push_back(vc_tally[s][y]);
        }
    }
    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = discard_years; y<num_years; ++y) {
            metrics.push_back(s_tally[s][y]);
        }
    }
    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = discard_years; y<num_years; ++y) {
            metrics.push_back(i_tally[s][y]);
        }
    }
    return metrics;
}


vector<long double> simulator(vector<long double> args, const unsigned long int rng_seed, const MPI_par* mp) {
    gsl_rng_set(RNG, 12345); // seed the rng using sys time and the process id
//    gsl_rng_set(RNG, rng_seed); // seed the rng using sys time and the process id
    const int years_simulated = 70;
    const int burnin = 50; // e.g., 0 means start vaccinating in the first simulated year
    const int discard_years = 30; // initial years of burn-in not to report

    Parameters* par = define_default_parameters(years_simulated); 

    vector<int> target_ages = {2,6,10,14};
    vector<float> vector_controls = {0.0, 0.1, 0.25, 0.5};
    vector<int> vaccine_durations = {-1, 4*365, 10*365, 20*365};

    bool vaccine                  = true;
    bool retro                    = false;
    bool catchup                  = true;
    int target                    = 2;
    int catchup_to                = 46;
    bool all_mature               = false;
    par->nDefaultMosquitoCapacity = 80;
    float vector_reduction        = 0.0;
    int vaccine_duration          = 20*365;
    bool boosting                 = true;
/*
    bool vaccine                  = (bool) args[0];
    bool retro                    = (bool) args[1];
    bool catchup                  = (bool) args[2];
    int target                    = target_ages[(unsigned int) args[3]];
    int catchup_to                = (int) args[4];
    bool all_mature               = (bool) args[5];
    par->nDefaultMosquitoCapacity = (int) args[6];
    float vector_reduction        = vector_controls[(unsigned int) args[8]];
    int vaccine_duration          = vaccine_durations[(unsigned int) args[9]]; 
    bool boosting                 = (bool) args[10];
*/
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
    if (not vaccine and (retro or catchup or all_mature or target > target_ages[0])) { nonsensical_parameters = true; } 
    if (nonsensical_parameters) {
        // 3 is because vaccinated cases, total cases, and infections are reported
        vector<long double> dummy((years_simulated-discard_years)*NUM_OF_SEROTYPES*3, 0.0);
        delete par;
        return dummy;
    }

    if (vaccine) {
        par->bVaccineLeaky = true;

        par->fVESs.clear();
        //par->fVESs = {0.7, 0.35, 0.7, 0.7};
        par->fVESs = {0.6, 0.54, 0.9, 0.95};

        par->fVESs_NAIVE.clear();
        if (all_mature) {
            par->fVESs_NAIVE = par->fVESs;
        } else {
            //par->fVESs_NAIVE = {0.35, 0.0, 0.35, 0.35};
            par->fVESs_NAIVE = {0.3, 0.27, 0.45, 0.48};
        }

        if (catchup) {
            for (int catchup_age = target + 1; catchup_age <= catchup_to; catchup_age++) {
                par->nVaccinateYear.push_back(burnin);
                par->nVaccinateAge.push_back(catchup_age);
                par->fVaccinateFraction.push_back(0.7);  // imperial used 0.5
                par->nSizeVaccinate++;
            }
        } 

        for (int vacc_year = burnin; vacc_year < years_simulated; vacc_year++) {
            par->nVaccinateYear.push_back(vacc_year);
            par->nVaccinateAge.push_back(target);
            par->fVaccinateFraction.push_back(0.7);      // imperial used 0.9
            par->nSizeVaccinate++;
        }

        if (retro) par->bRetroactiveMatureVaccine = true;
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
//    cerr << "using vector reduction: " << vector_reduction << endl;
        assert( vector_reduction <= 1.0); 
        assert( par->mosquitoMultipliers.size() == 12); // expecting values for 12 months
        int running_sum = par->mosquitoMultipliers.back().start + DAYS_IN_MONTH[0]; // start on january 1 of year two
        const int months_simulated = years_simulated*12;
        par->mosquitoMultipliers.resize(months_simulated);
        for (unsigned int i = 12; i < months_simulated; ++i) {
            const int month_of_year = i % 12;
            par->mosquitoMultipliers[i].start = running_sum;
            par->mosquitoMultipliers[i].duration = DAYS_IN_MONTH[month_of_year];
            float vector_coeff = i/12 >= burnin ? 1.0 - vector_reduction : 1.0; // normally, there's no reduction in vectors
            par->mosquitoMultipliers[i].value = vector_coeff * par->mosquitoMultipliers[month_of_year].value;
//cerr << par->mosquitoMultipliers[i].value << " ";
            running_sum += DAYS_IN_MONTH[month_of_year];
        }
//cerr << endl;
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
    float approx_attack_rate = (float) args[7]/100.0; // approximate per-serotype probability someone will be infected in a given year
    generate_homogeneous_immune_history(community, approx_attack_rate);

    seed_epidemic(par, community);
    simulate_epidemic(par, community, process_id);
    //vector<int> epi_sizes = simulate_epidemic(par, community, process_id);

    time (&end);
    double dif = difftime (end,start);

    vector<long double> metrics = tally_counts(par, community, discard_years);

    stringstream ss;
    ss << mp->mpi_rank << " END " << hex << process_id << " " << dec << dif << " ";

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
