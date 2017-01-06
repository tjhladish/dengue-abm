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

string calculate_process_id(vector<double> &args, string &argstring);
const string SIM_POP = "merida";
const string HOME_DIR(std::getenv("HOME"));
const string pop_dir = HOME_DIR + "/work/dengue/pop-" + SIM_POP;
const string output_dir("/ufrc/longini/tjhladish/");

const int RESTART_BURNIN    = 25;
const int FORECAST_DURATION = 51;
const bool RUN_FORECAST     = true;
const int TOTAL_DURATION    = RUN_FORECAST ? RESTART_BURNIN + FORECAST_DURATION : RESTART_BURNIN;
const vector<float> SP9_TARGETS = {0.1, 0.3, 0.5, 0.7, 0.9};

Parameters* define_simulator_parameters(vector<double> args, const unsigned long int rng_seed) {
    Parameters* par = new Parameters();
    par->define_defaults();

//  0 mildEF real,
//  1 severeEF real,
//  2 sec_sev real,
//  3 base_path real,
//  3 pss_ratio real,
//  4 exp_coef real,
//  5 num_mos real,
//  6 beta real,
//
//  Pseudo pars
//  7 beta_multiplier real,

    // The following FOI targets are the needed product of # mosquitoes * beta in order to achieve 9 yo equillibrium seroprevalences
    // of 10, 30, 50, 70, and 90%
    //    vector<float> foi_targets = {6.5, 7.5, 9, 12, 20}; // NEEDS TO BE UPDATED IF TRANSMISSION MODEL CHANGES
//0.111192782416|0.0
//0.20122664603|1.0
//0.3716522091|2.0
//0.668049612|3.0
//0.894214538|4.0
    //const vector<float> foi_targets = {6.2, 8.4, 10.4, 12.5, 20}; // NEEDS TO BE UPDATED IF TRANSMISSION MODEL CHANGES
//0.11739814169|0.0
//0.33525855234|1.0
//0.5314984042|2.0
//0.646959106000001|3.0
//0.868876936|4.0

    const vector<float> foi_targets = {6.07, 8.15, 10.35, 13.53, 21.21}; // try 6
    /* 0.0|0.0960594171 // I don't think we can do better than this without a larger sample size
     * 1.0|0.305756767
     * 2.0|0.51102584
     * 3.0|0.69006929
     * 4.0|0.89070218
     */

    double _mild_EF      = args[0];
    double _severe_EF    = args[1];
    double _base_path    = args[2];
    double _sec_severity = args[3];
    double _pss_ratio    = args[4];
    double _exp_coef     = args[5];
    double _nmos         = args[6];
    double _betamp       = foi_targets[(int) args[7]] / _nmos;   // mp and pm and not separately
    double _betapm       = _betamp;                              // identifiable, so they're the same
    par->reportedFraction = {0.0, 1.0/_mild_EF, 1.0/_severe_EF}; // no asymptomatic infections are reported

    vector<double> pid_args(&args[0], &args[14]);  // really fucking ugly, but I have to hack the args to match the pid for the right immunity file
    pid_args[8] = 0;
    pid_args[9] = 0;
    pid_args[10] = 0;
    pid_args[11] = 9;
    pid_args[12] = 0;
    pid_args[13] = 0;

    string argstring;
    const string imm_file_pid = calculate_process_id(pid_args, argstring);

    cerr << "pid in define_simulator_parameters: " << imm_file_pid << endl;

    par->randomseed              = rng_seed;
    par->dailyOutput             = false;
    par->weekyOutput             = true;
    par->monthlyOutput           = false;
    par->yearlyOutput            = true;
    par->abcVerbose              = false; // needs to be false to get WHO daily output
    const int runLengthYears     = TOTAL_DURATION;
    par->nRunLength              = runLengthYears*365;
    par->startDayOfYear          = 100;
    par->annualIntroductionsCoef = pow(10, _exp_coef);

    // pathogenicity values fitted in
    // Reich et al, Interactions between serotypes of dengue highlight epidemiological impact of cross-immunity, Interface, 2013
    // Normalized from Fc values in supplement table 2, available at
    // http://rsif.royalsocietypublishing.org/content/10/86/20130414/suppl/DC1
//par->defineSerotypeRelativeRisks();
par->serotypePathogenicityRelativeRisks = vector<double>(NUM_OF_SEROTYPES, 1.0);
    //par->useAgeStructuredPrimaryPathogenicity = false;
    par->primaryPathogenicityModel = CONSTANT_PATHOGENICITY;
    par->basePathogenicity = _base_path;
    par->primaryRelativeRisk = 0.45/_base_path;
    par->postSecondaryRelativeRisk = 1.0/6.0;

    par->primarySevereFraction    = vector<double>(NUM_OF_SEROTYPES, _sec_severity*_pss_ratio);
    par->secondarySevereFraction  = vector<double>(NUM_OF_SEROTYPES, _sec_severity);
    par->tertiarySevereFraction   = vector<double>(NUM_OF_SEROTYPES, _sec_severity/4.0);
    par->quaternarySevereFraction = vector<double>(NUM_OF_SEROTYPES, _sec_severity/4.0);

    par->betaPM = _betapm;
    par->betaMP = _betamp;
    par->fMosquitoMove = 0.15;
    par->mosquitoMoveModel = "weighted";
    par->fMosquitoTeleport = 0.0;
    par->nDefaultMosquitoCapacity = (int) _nmos;
    par->eMosquitoDistribution = EXPONENTIAL;

    par->nDaysImmune = 730;
    par->fVESs = vector<double>(NUM_OF_SEROTYPES, 0.0);

    //par->hospitalizedFraction = 0.25; // fraction of cases assumed to be hospitalized

    //par->simulateAnnualSerotypes = true;
    //par->normalizeSerotypeIntros = true;
    //if (par->simulateAnnualSerotypes) par->generateAnnualSerotypes();
    par->nDailyExposed = {{0.25, 0.25, 0.25, 0.25}};

    par->annualIntroductions = {1.0};

    // load daily EIP
    // EIPs calculated by ../raw_data/weather/calculate_daily_eip.R, based on reconstructed yucatan temps
    par->loadDailyEIP(pop_dir + "/seasonal_EIP_24hr.out");

    // we add the startDayOfYear offset because we will end up discarding that many values from the beginning
    // mosquitoMultipliers are indexed to start on Jan 1
    par->loadDailyMosquitoMultipliers(pop_dir + "/mosquito_seasonality.out", par->nRunLength + par->startDayOfYear);

    par->populationFilename       = pop_dir    + "/population-"         + SIM_POP + ".txt";
    par->immunityFilename         = "/ufrc/longini/tjhladish/imm_who-baseline-seroprev-july2016/immunity." + imm_file_pid;
    //par->immunityFilename         = "";
    par->locationFilename         = pop_dir    + "/locations-"          + SIM_POP + ".txt";
    par->networkFilename          = pop_dir    + "/network-"            + SIM_POP + ".txt";
    par->swapProbFilename         = pop_dir    + "/swap_probabilities-" + SIM_POP + ".txt";
    //par->mosquitoFilename         = output_dir + "/mos_mer_who/mos."       + to_string(process_id);
    //par->mosquitoLocationFilename = output_dir + "/mosloc_mer_who/mosloc." + to_string(process_id);

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

    cerr << "pid in report_process_id: " << process_id << endl;
    stringstream ss;
    ss << mp->mpi_rank << " begin " << process_id << " " << dec << serial << " " << dif << " " << argstring << endl;
    string output = ss.str();
    fputs(output.c_str(), stderr);

    return to_string(process_id);
}


void append_if_finite(vector<double> &vec, double val) {
    if (isfinite(val)) { 
        vec.push_back((double) val);
    } else {
        vec.push_back(0);
    }
}


void prime_population(Community* community, const gsl_rng* RNG, double p_prime) { // p_prime is per-sero, yearly probability of infection
    for (int i = 0; i < community->getNumPerson(); ++i) {
        Person* p = community->getPerson(i);
        const int age = p->getAge();
        const int birthday = gsl_rng_uniform_int(RNG, 365);
        const int age_days = age*365 + birthday;
        if (age_days == 0) continue;
        const double p_infec = 1.0 - pow(1.0-p_prime, age);
        for (int s = 0; s < (int) NUM_OF_SEROTYPES; ++s) {
            if (p_infec < gsl_rng_uniform(RNG)) {
                const int day = -1*gsl_rng_uniform_int(RNG, age_days); // negative day == before "now"/start of simulation
                p->infect((Serotype) s, day);
            }
        }
    }
}


vector<double> tally_counts(const Parameters* par, Community* community, const int vc_timing) {
    // aggregate based on the timing of the annual start of vector control
    int discard_days = INT_MAX;
    for (VectorControlEvent vce: par->vectorControlEvents) {
        discard_days = discard_days < vce->campaignStart ? discard_days : vce->campaignStart;
    }
    const int pre_intervention_output = 5; // years
    discard_days = discard_days==INT_MAX ? 365*(RESTART_BURNIN-pre_intervention_output) : discard_days - (365*pre_intervention_output);

    //vector< vector<int> > severe      = community->getNumSevereCases();
    vector< vector<int> > symptomatic = community->getNumNewlySymptomatic();

    //vector< vector<int> > infected    = community->getNumNewlyInfected();
    //const int num_years = FORECAST_DURATION;
    const int num_years = 55;
    vector<vector<int> > s_tally(NUM_OF_SEROTYPES, vector<int>(num_years+1, 0));

    vector<double> metrics(num_years, 0.0);
    for (int t=discard_days; t<par->nRunLength; t++) {
        // use epidemic years, instead of calendar years
        const int y = (t-discard_days)/365;
        for (int s=0; s<NUM_OF_SEROTYPES; s++) {
            s_tally[s][y] += symptomatic[s][t];
        }
        //cout << "d,i:" << t << "," << infected[0][t] + infected[1][t] + infected[2][t] + infected[3][t] << endl;
    }
    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = 0; y<num_years; ++y) {
            metrics[y] += s_tally[s][y];
        }
    }

/*
    vector<int> pop_sizes  = vector<int>(10,0);
    vector<double> seropos = vector<double>(10,0.0); // seronegative at vaccination time

    for (int i=community->getNumPerson()-1; i>=0; i--) {
        Person *p = community->getPerson(i);
        const int final_age  = p->getAge();
        const int age_at_vac_start = RUN_FORECAST ? final_age - FORECAST_DURATION : final_age;
        if (age_at_vac_start >= 9 and age_at_vac_start <= 18) {
            int age_class = age_at_vac_start - 9;

            const int serosurvey_day = (RESTART_BURNIN - age_class)*365;
            const bool ever_infected = p->getNumNaturalInfections() > 0 ? true : false;
            const bool isSeropos = ever_infected and p->getInfectionHistory().front()->getInfectedTime() < serosurvey_day;
            ++pop_sizes[age_class];
            if (isSeropos) ++seropos[age_class];
        }
    }
    double seroprev_mean = 0.0;
    for (unsigned int i = 0; i < seropos.size(); ++i) seroprev_mean += seropos[i] / pop_sizes[i];
    assert(seropos.size() > 0);
    seroprev_mean /= seropos.size();

    metrics.push_back(seroprev_mean);
*/
    return metrics;
}


vector<double> simulator(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const ABC::MPI_par* mp) {
    gsl_rng_set(RNG, rng_seed); // seed the rng using sys time and the process id

    //initialize bookkeeping for run
    time_t start, end;
    time (&start);
    const string process_id = report_process_id(args, serial, mp, start);

    cerr << "SCENARIO " << rng_seed;
    for (auto _p: args) cerr << " " << _p; cerr << endl;

    Parameters* par = define_simulator_parameters(args, rng_seed);

    //vector<int> vaccine_durations = {-1, 10*365};
//  Posterior pars
//  0 mildEF real,
//  1 severeEF real,
//  2 base_path real
//  3 sec_sev real,
//  4 pss_ratio real
//  5 exp_coef real,
//  6 num_mos real,
//
//  Pseudo pars
//  7 foi_target real,
//  8 catchup int,
//  9 vaccine_mechanism int,
//  10 vaccine int,
//  11 target int,
//  12 catchup_to int,
//  13 coverage int 
//
//  14 vector_control
//  15 campaign_duration
//  16 timing
//  17 vc_coverage
//  18 vc_efficacy

    //vector<int> target_ages = {9, 16};         // default 9
    vector<int> catchup_ages = {17, 30};         // default 17
    vector<double> coverage_levels = {0.5, 0.8}; // default 0.8
    vector<int> vc_campaign_duration_levels = {1, 90, 365};

    bool catchup           = (bool) args[8];
    int vaccine_mechanism  = (int) args[9];

    bool vaccine           = (bool) args[10];
    int target             = (int) args[11];
    int catchup_to         = catchup_ages[(int) args[12]];
    double coverage        = coverage_levels[(int) args[13]];

    const bool vector_control     = (bool) args[14];
    //const int vc_campaignDuration = (int) args[15]; // number of days to achieve coverage
    const int vc_campaignDuration = vc_campaign_duration_levels[(int) args[15]]; // number of days to achieve coverage
    const int vc_timing           = (int) args[16]; 
    const double vc_coverage      = args[17];
    const double vc_efficacy       = args[18];       // expected % reduction in equillibrium mosquito population in treated houses

    bool nonsensical_parameters = false;
    // only run a non-vaccination campaign if all the vaccine parameters are 0
    // TODO - this should be reworked with new "catchup-to" parameter
    if (not vaccine and catchup) { nonsensical_parameters = true; }
    if (nonsensical_parameters) {
        vector<double> dummy(FORECAST_DURATION + 1, 0.0);
        delete par;
        return dummy;
    }

    Community* community = build_community(par);

    if (vector_control) {
        const int efficacyDuration = 90;       // number of days efficacy is maintained
        const LocationType locType = HOME;
        const LocationSelectionStrategy lss = UNIFORM_STRATEGY;
        for (int vec_cont_year = RESTART_BURNIN; vec_cont_year < TOTAL_DURATION; vec_cont_year++) {
        //for (int vec_cont_year = 0; vec_cont_year < TOTAL_DURATION; vec_cont_year++) {
            // TODO - address situation where startDate could be negative if startDayOfYear is larger than vc_timing
            //const int startDate = (vec_cont_year*365) + (vc_timing - par->startDayOfYear)%365; // 151 days after Jan 1 is June 1, offset by julian startDay
            int startDate = 0;
            if (vc_timing >= par->startDayOfYear) { // start before Jan 1
                startDate = vec_cont_year*365 + vc_timing - par->startDayOfYear; // 151 days after Jan 1 is June 1, offset by julian startDay
            } else {                                // start after Jan 1
                startDate = (vec_cont_year+1)*365 + vc_timing - par->startDayOfYear;
            }
            par->vectorControlEvents.emplace_back(startDate, vc_campaignDuration, vc_coverage, vc_efficacy, efficacyDuration, locType, lss);
        }
    }

    if (vaccine) {
        double target_coverage  = coverage;
        double catchup_coverage = coverage;

        if (catchup) {
            for (int catchup_age = target + 1; catchup_age <= catchup_to; catchup_age++) {
                par->vaccinationEvents.emplace_back(catchup_age, RESTART_BURNIN*365, catchup_coverage);
            }
        } 

        for (int vacc_year = RESTART_BURNIN; vacc_year < TOTAL_DURATION; vacc_year++) {
            par->vaccinationEvents.emplace_back(target, vacc_year*365, target_coverage);
        }
    }

    if (vaccine_mechanism == 0) {        // "baseline" scenario: A2b + B2 + C3a
        // perfect efficacy that wanes rapidly -- most benefit comes from vaccine-as-infection assumption
        par->fVESs                   = vector<double>(NUM_OF_SEROTYPES, 1.0);
        par->fVESs_NAIVE             = vector<double>(NUM_OF_SEROTYPES, 1.0);
        par->fVEH                    = 0.803; // fraction of hospitalized cases prevented by vaccine
        par->linearlyWaningVaccine   = true;
        par->vaccineImmunityDuration = 2*365;
        par->bVaccineLeaky           = true;
        par->numVaccineDoses         = 3;
        par->vaccineDoseInterval     = 182;
        par->vaccineBoosting         = false;

        par->whoDiseaseOutcome = INC_NUM_INFECTIONS;
        par->whoBreakthrough   = BREAKTHROUGH_SEROCONVERSION;
        par->whoWaning         = UNIVERSAL_WANING; // Naive-only was specified, but seems unrealistic
    } else if (vaccine_mechanism == 1) { // also want A1 + B2 + C3a
        // imperfect efficacy that does not wane
        par->fVESs = {0.6, 0.54, 0.9, 0.95};
        par->fVESs_NAIVE = {0.3, 0.27, 0.45, 0.48};

        par->fVEH                    = 0.803; // fraction of hospitalized cases prevented by vaccine
        par->linearlyWaningVaccine   = false;
        par->vaccineImmunityDuration = INT_MAX;
        par->bVaccineLeaky           = true;
        par->numVaccineDoses         = 3;
        par->vaccineDoseInterval     = 182;
        par->vaccineBoosting         = false;

        par->whoDiseaseOutcome = VAC_ISNT_INFECTION;
        par->whoBreakthrough   = BREAKTHROUGH_SEROCONVERSION;
        par->whoWaning         = UNIVERSAL_WANING; // Naive-only was specified, but seems unrealistic
    } else {
        cerr << "Unsupported vaccine mechanism: " << vaccine_mechanism << endl;
        exit(-152);
    }

    //const double AR   = 1.0 - pow(1.0 - SP9_TARGETS[(int) args[7]], 1.0/9.0); // Overall target attack rate
    //const double AR_S = 1.0 - pow(1.0 - AR, 1.0/4.0);         // Per-serotype target AR

    //prime_population(community, RNG, AR_S);

    seed_epidemic(par, community);
    simulate_epidemic(par, community, process_id);

    //string imm_filename = "/ufrc/longini/tjhladish/imm_who-baseline-seroprev-july2016/immunity." + process_id;
    //write_immunity_file(community, process_id, imm_filename, par->nRunLength);

    time (&end);
    double dif = difftime (end,start);

    vector<double> metrics = tally_counts(par, community, vc_timing);

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
        metrics[i] = ((double) 1e3*epi_sizes[i])/pop_size; 
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
        cerr << "simulate db\n";
        time(&GLOBAL_START_TIME);
        abc->set_simulator(simulator);
        abc->simulate_next_particles(buffer_size);
    }

    return 0;
}
