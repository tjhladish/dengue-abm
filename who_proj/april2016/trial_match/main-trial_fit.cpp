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
using ABC::float_type;

time_t GLOBAL_START_TIME;

string calculate_process_id(vector< double> &args, string &argstring);
const string SIM_POP = "yucatan";
const string HOME(std::getenv("HOME"));
const string pop_dir = HOME + "/work/dengue/pop-" + SIM_POP;
const string output_dir("/scratch/lfs/thladish");

const int FIRST_YEAR          = 0;                                  // inclusive
const int LAST_YEAR           = 50;                                 // inclusive
const int HISTORICAL_DURATION = LAST_YEAR - FIRST_YEAR + 1; // HISTORICAL == before intervention period
//const int BURNIN_YEARS        = 125;                              // should be zero if we're running from FIRST_YEAR; 125 should correspond to resuming after 2003
//const int RESTART_BURNIN    = 80;
const int FORECAST_DURATION = 6;
const int TOTAL_DURATION    = HISTORICAL_DURATION + FORECAST_DURATION;
const int NUM_OF_ABC_PARS   = 2;

Parameters* define_simulator_parameters(vector<double> args, const unsigned long int rng_seed) {
    Parameters* par = new Parameters();
    par->define_defaults();

//  Fitted pars
//  0 foi_multiplier real,
//  1 vw_rate real,
//  2 initial_efficacy real,

    //const vector<float> foi_targets = {8,9,10,11,12,13,14,15,16,17}; // NEEDS TO BE UPDATED IF TRANSMISSION MODEL CHANGES
    //const vector<float> beta_multipliers = {1.50, 1.60, 1.70, 1.80, 1.90};
    //const vector<float> beta_multipliers = {0.52, 0.70, 0.92, 1.25, 2.25};   // from dec fit

//  sqlite> .schema parameters
//  CREATE TABLE parameters( serial int primary key, seed blob, mild_EF real, severe_EF real, base_path real, sec_sev real, pss_ratio real, exp_coef real, num_mos real );
//  sqlite> select * from parameters where serial = 49050; // serial corresponds to best-ranking particle
//  49050|152353960|10.628|4.0888|0.595058|0.00557066|0.675647|-0.416732|82.8077


    double _mild_EF      = 10.628;
    double _severe_EF    = 4.0888;
    double _base_path    = 0.595058;
    double _sec_severity = 0.00557066;
    double _pss_ratio    = 0.675647;
    double _exp_coef     = -0.416732;
    double _nmos         = 82.8077*2.3049;
    double _betamp       = 0.25;//foi_targets[(int) args[7]] / _nmos;
    double _betapm       = 0.10;//_betamp;
    par->reportedFraction = {0.0, 1.0/_mild_EF, 1.0/_severe_EF}; // no asymptomatic infections are reported

    vector<double> abc_args(&args[0], &args[NUM_OF_ABC_PARS]);
    string argstring;
    const string process_id = to_string(calculate_process_id(abc_args, argstring));

    par->randomseed              = rng_seed;
    //par->dailyOutput             = true;
    par->monthlyOutput           = true;
    par->yearlyOutput            = true;
    par->abcVerbose              = true; // needs to be false to get WHO daily output
    par->nRunLength              = TOTAL_DURATION*365;
    par->startDayOfYear          = 100;
    par->annualIntroductionsCoef = pow(10, _exp_coef);

    // pathogenicity values fitted in
    // Reich et al, Interactions between serotypes of dengue highlight epidemiological impact of cross-immunity, Interface, 2013
    // Normalized from Fc values in supplement table 2, available at
    // http://rsif.royalsocietypublishing.org/content/10/86/20130414/suppl/DC1
    par->defineSerotypeRelativeRisks();
//par->serotypePathogenicityRelativeRisks = vector<double>(NUM_OF_SEROTYPES, 1.0);
    par->basePathogenicity = _base_path;
    par->primaryPathogenicityModel = ORIGINAL_LOGISTIC;
    //par->annualFlavivirusAttackRate = _flav_ar;
    //par->primaryRelativeRisk = 0.45/_base_path;
    //par->primaryPathogenicityModel = ORIGINAL_LOGISTIC;
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
    par->immunityFilename         = "";
    par->locationFilename         = pop_dir    + "/locations-"          + SIM_POP + ".txt";
    par->networkFilename          = pop_dir    + "/network-"            + SIM_POP + ".txt";
    par->swapProbFilename         = pop_dir    + "/swap_probabilities-" + SIM_POP + ".txt";
    //par->mosquitoFilename         = output_dir + "/mos_mer_who/mos."       + to_string(process_id);
    //par->mosquitoLocationFilename = output_dir + "/mosloc_mer_who/mosloc." + to_string(process_id);

    return par;
}


void prime_population(Community* community, const gsl_rng* RNG, const double target_SP9) { // target_SP9 = target seroprev in 9 year olds at equillibrium
    const double AR   = 1.0 - pow(1.0 - target_SP9, 1.0/9.0); // Overall target attack rate
    const double AR_S = 1.0 - pow(1.0 - AR, 1.0/4.0);         // Per-serotype target AR
    for (int i = 0; i < community->getNumPerson(); ++i) {
        Person* p = community->getPerson(i);
        int age = p->getAge();
        const int birthday = gsl_rng_uniform_int(RNG, 365);
        const int age_days = age*365 + birthday;
        if (age_days == 0) continue;
        const double p_infec = 1.0 - pow(1.0-AR_S, age_days/365.0);
        for (int s = 0; s < (int) NUM_OF_SEROTYPES; ++s) {
            if (p_infec < gsl_rng_uniform(RNG)) {
                const int day = -1*gsl_rng_uniform_int(RNG, age_days); // negative day == before "now"/start of simulation
                p->infect((Serotype) s, day);
            }
        }
    }
}



string calculate_process_id(vector< double> &args, string &argstring) {
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

    stringstream ss;
    ss << mp->mpi_rank << " begin " << process_id << " " << dec << serial << " " << dif << " " << argstring << endl;
    string output = ss.str();
    fputs(output.c_str(), stderr);

    return to_string(process_id);
}
    

vector<double> tally_counts(Parameters* par, Community* community) {
    const int VAC_DAY = HISTORICAL_DURATION*365;
    vector<double> metrics;

    vector< vector<int> > age_bins = {{2,5},{6,11},{12,14}};
    vector< vector<double> > infection_attack_rate_age  = {{0,0},{0,0},{0,0}}; 
    vector< vector<double> > infection_attack_rate_sero = {{0,0},{0,0}}; 

    vector<int> pop_sizes                            = {0,0,0};
    vector< vector<int> > arm_pop_sizes              = {{0,0},{0,0},{0,0}}; // age; placebo = 0, vaccine = 1
    vector< vector<int> > sero_pop_sizes             = {{0,0},{0,0}};       // seroneg = 0, seropos = 1; placebo = 0, vaccine = 1
    vector<double> seroneg                           = {0,0,0}; // seronegative at vaccination time
    vector< vector<double> > attack_rate_age         = {{0,0},{0,0},{0,0}};
    vector< vector<double> > attack_rate_sero        = {{0,0},{0,0}};
    // Dimensions are [age_class][vaccination_status][hospital_phase_year]
    // Only first year of hospital phase is used for fitting; others only reporting for model validation
    vector< vector< vector<double> > > mild_attack_rate_hosp   = {{{0,0,0},{0,0,0}},{{0,0,0},{0,0,0}},{{0,0,0},{0,0,0}}};
    vector< vector< vector<double> > > severe_attack_rate_hosp = {{{0,0,0},{0,0,0}},{{0,0,0},{0,0,0}},{{0,0,0},{0,0,0}}};

    for (int i=community->getNumPerson()-1; i>=0; i--) {
        Person *p = community->getPerson(i);
        const int final_age  = p->getAge();
        const int age_at_vac = final_age - FORECAST_DURATION;
        if (age_at_vac >= age_bins.front()[0] and age_at_vac <= age_bins.back()[1]) {
            int age_class = 0;
            for (auto bin: age_bins) {
                if (age_at_vac > bin[1]) ++age_class;
            }
            const bool ever_infected = p->getNumNaturalInfections() > 0 ? true : false;
            const bool seropos = ever_infected and p->getInfectionHistory().front()->getInfectedTime() < VAC_DAY;
            const int vaccinated = (int) p->isVaccinated();
            ++pop_sizes[age_class];
            ++arm_pop_sizes[age_class][vaccinated];
            ++sero_pop_sizes[(int) seropos][vaccinated];
            if (not seropos) ++seroneg[age_class];

            for (const Infection* infec: p->getInfectionHistory()) {
                const int itime = infec->getInfectedTime();
                // sanity check -- tally all infections
                if (itime >= VAC_DAY + 365 + 28 and itime < VAC_DAY + (365*2) + 28) {
                    ++infection_attack_rate_age[age_class][vaccinated];
                    ++infection_attack_rate_sero[(int) seropos][vaccinated];
                }
                // We only care about symptomatic infections for trial data comparison
                if (not infec->isSymptomatic()) continue;
                if (itime >= VAC_DAY + 365 + 28 and itime < VAC_DAY + (365*2) + 28) {
                    // active surveillance phase
                    ++attack_rate_age[age_class][vaccinated];
                    ++attack_rate_sero[(int) seropos][vaccinated];
                } else if (itime >= VAC_DAY + (365*2) + 28){// and itime < VAC_DAY + (365*3) + 28) {
                    // hospital phase
                    const int hosp_year = (itime - (VAC_DAY + (365*2) + 28)) / 365;
                    if ((const bool) infec->isSevere()) {
                        if (hosp_year < (signed) severe_attack_rate_hosp[age_class][vaccinated].size()) { 
                            ++severe_attack_rate_hosp[age_class][vaccinated][hosp_year];
                        }
                    } else {
                        if (hosp_year < (signed) mild_attack_rate_hosp[age_class][vaccinated].size()) { 
                            ++mild_attack_rate_hosp[age_class][vaccinated][hosp_year];
                        }
                    }
                }
            }
        }
    }

    //vector<long double> metrics;
    //vector< vector<int> > age_bins = {{2,5},{6,11},{12,14}};

    //vector<int> pop_sizes                            = {0,0,0};
    //vector<int> arm_pop_sizes                        = {{0,0},{0,0},{0,0}}; // age; placebo = 0, vaccine = 1
    //vector<int> sero_pop_sizes                       = {{0,0},{0,0}};       // seroneg = 0, seropos = 1; placebo = 0, vaccine = 1
    //
    //vector<double> seroneg                           = {0,0,0}; // seronegative at vaccination time
    //vector< vector<double> > attack_rate_age         = {{0,0},{0,0},{0,0}};
    //vector< vector<double> > attack_rate_sero        = {{0,0},{0,0}};
    //vector< vector< vector<double> > > mild_attack_rate_hosp   = {{{0,0,0},{0,0,0}},{{0,0,0},{0,0,0}},{{0,0,0},{0,0,0}}};
    //vector< vector< vector<double> > > severe_attack_rate_hosp = {{{0,0,0},{0,0,0}},{{0,0,0},{0,0,0}},{{0,0,0},{0,0,0}}};

//        {"name" : "seroneg_2_5",        "num_type" : "FLOAT", "value" : 0.488338192},
//        {"name" : "seroneg_6_11",       "num_type" : "FLOAT", "value" : 0.279266573},
//        {"name" : "seroneg_12_14",      "num_type" : "FLOAT", "value" : 0.188647746},
    for (unsigned int i = 0; i < seroneg.size(); ++i) metrics.push_back(seroneg[i] / pop_sizes[i]);

//        {"name" : "ar_placebo_2_5",     "num_type" : "FLOAT", "value" : 0.056186869},
//        {"name" : "ar_vaccine_2_5",     "num_type" : "FLOAT", "value" : 0.037278658},
//        {"name" : "ar_placebo_6_11",    "num_type" : "FLOAT", "value" : 0.046821793},
//        {"name" : "ar_vaccine_6_11",    "num_type" : "FLOAT", "value" : 0.018951446},
//        {"name" : "ar_placebo_12_14",   "num_type" : "FLOAT", "value" : 0.03630363},
//        {"name" : "ar_vaccine_12_14",   "num_type" : "FLOAT", "value" : 0.009285943},
    for (unsigned int i = 0; i < attack_rate_age.size(); ++i) {
        for (unsigned int j = 0; j < attack_rate_age[i].size(); ++j) metrics.push_back(attack_rate_age[i][j] / arm_pop_sizes[i][j]);
    }

//        {"name" : "ar_placebo_neg",     "num_type" : "FLOAT", "value" : 0.042553191},
//        {"name" : "ar_vaccine_neg",     "num_type" : "FLOAT", "value" : 0.027446301},
//        {"name" : "ar_placebo_pos",     "num_type" : "FLOAT", "value" : 0.038636364},
//        {"name" : "ar_vaccine_pos",     "num_type" : "FLOAT", "value" : 0.00993926},
    for (unsigned int i = 0; i < attack_rate_sero.size(); ++i) {
        for (unsigned int j = 0; j < attack_rate_sero[i].size(); ++j) metrics.push_back(attack_rate_sero[i][j] / sero_pop_sizes[i][j]);
    }

//        {"name" : "hosp_placebo_2_5",   "num_type" : "FLOAT", "value" : 0.001230012},
//        {"name" : "hosp_vaccine_2_5",   "num_type" : "FLOAT", "value" : 0.009168704},
//        {"name" : "hosp_placebo_6_11",  "num_type" : "FLOAT", "value" : 0.004429679},
//        {"name" : "hosp_vaccine_6_11",  "num_type" : "FLOAT", "value" : 0.002779322},
//        {"name" : "hosp_placebo_12_14", "num_type" : "FLOAT", "value" : 0.005208333},
//        {"name" : "hosp_vaccine_12_14", "num_type" : "FLOAT", "value" : 0.001295337}
    assert(mild_attack_rate_hosp.size() == severe_attack_rate_hosp.size());
    //cerr << "Calculating hospital attack rates\n";
    for (unsigned int i = 0; i < mild_attack_rate_hosp.size(); ++i) {
        assert(mild_attack_rate_hosp[i].size() == severe_attack_rate_hosp[i].size());
        for (unsigned int j = 0; j < mild_attack_rate_hosp[i].size(); ++j) {
            assert(mild_attack_rate_hosp[i][j].size() == severe_attack_rate_hosp[i][j].size());
            for (unsigned int y = 0; y < mild_attack_rate_hosp[i][j].size(); ++y) {
                const double hosp_ar = severe_attack_rate_hosp[i][j][y] + (mild_attack_rate_hosp[i][j][y] * par->reportedFraction[1]);
                //cerr << "i, j, y | severe, mild, reported_frac, hosp_ar: " << i << " " << j << " " << y << " | " 
                //     << severe_attack_rate_hosp[i][j][y] << " + " << mild_attack_rate_hosp[i][j][y] << " * " << par->reportedFraction[1] << " = " << hosp_ar << endl;
                metrics.push_back(hosp_ar / arm_pop_sizes[i][j]); // we're pushing three years of hosp results, but last two will be discarded after logging
            }
        }
    }

    // 10 more metrics, to be logged but not reported to ABC in calling function
    for (unsigned int i = 0; i < infection_attack_rate_age.size(); ++i) {
        for (unsigned int j = 0; j < infection_attack_rate_age[i].size(); ++j) metrics.push_back(infection_attack_rate_age[i][j] / arm_pop_sizes[i][j]);
    }
    for (unsigned int i = 0; i < infection_attack_rate_sero.size(); ++i) {
        for (unsigned int j = 0; j < infection_attack_rate_sero[i].size(); ++j) metrics.push_back(infection_attack_rate_sero[i][j] / sero_pop_sizes[i][j]);
    }

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

//  Fitted pars
//  0 foi_multiplier real,
//  1 vw_rate real,
//  2 initial_efficacy real,

    int youngest    = 2;
    int oldest      = 14;
    double coverage = 0.07; // approximately what is necessary to get CYD14-like numbers of vaccinees
    const double target_SP9 = 0.7;

    Community* community = build_community(par);
    prime_population(community, RNG, target_SP9);

    // perfect efficacy that wanes rapidly
    const double EFF = NUM_OF_ABC_PARS == 2 ? args[1] : 1.0;
    par->fVESs                   = vector<double>(NUM_OF_SEROTYPES, EFF);
    par->fVESs_NAIVE             = vector<double>(NUM_OF_SEROTYPES, EFF);
    par->fVEH                    = 0.803; // fraction of hospitalized cases prevented by vaccine
    par->linearlyWaningVaccine   = true;
    par->vaccineImmunityDuration = NUM_OF_ABC_PARS > 1 ? (int) (1.0/args[0]) : 100000;
    par->bVaccineLeaky           = true;
    par->numVaccineDoses         = 3;
    par->vaccineDoseInterval     = 182;
    par->vaccineBoosting         = false;

    for (int vac_age = youngest; vac_age <= oldest; ++vac_age) {
        par->vaccinationEvents.emplace_back(vac_age, HISTORICAL_DURATION*365, coverage);
    }

    par->whoDiseaseOutcome = INC_NUM_INFECTIONS;
    par->whoBreakthrough   = BREAKTHROUGH_SEROCONVERSION;
    par->whoWaning         = UNIVERSAL_WANING; // Naive-only was specified, but seems unrealistic

    seed_epidemic(par, community);
    simulate_epidemic(par, community, process_id);

    time (&end);
    double dif = difftime (end,start);

    vector<double> metrics = tally_counts(par, community);

    stringstream ss;
    ss << mp->mpi_rank << " end " << hex << process_id << " " << dec << dif << " ";

    for (auto i: args) ss << i << " ";
    for (auto i: metrics) ss << i << " ";
    ss << endl;

    string output = ss.str();
    fputs(output.c_str(), stderr);

    const int NUM_METRICS_NEEDED = NUM_OF_ABC_PARS == 1 ? 3 : 19;
    metrics.resize(NUM_METRICS_NEEDED); // Must equal expected number of metrics.  We are passing some extra values for logging/validation purposes

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
