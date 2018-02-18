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

const unsigned int calculate_process_id(vector< long double> &args, string &argstring);
const string SIM_POP = "yucatan";
const string HOME(std::getenv("HOME"));
const string pop_dir = HOME + "/work/dengue/pop-" + SIM_POP;
const string output_dir("/scratch/lfs/thladish");

const int FIRST_YEAR          = 2003;                                 // inclusive
const int LAST_YEAR           = 2032;                                 // inclusive
const int HISTORICAL_DURATION = LAST_YEAR - FIRST_YEAR + 1; // HISTORICAL == before intervention period
//const int BURNIN_YEARS        = 125;                                  // should be zero if we're running from FIRST_YEAR; 125 should correspond to resuming after 2003
//const int RESTART_BURNIN    = 80;
const int FORECAST_DURATION = 6;
const int TOTAL_DURATION    = HISTORICAL_DURATION + FORECAST_DURATION;

Parameters* define_simulator_parameters(vector<long double> args, const unsigned long int rng_seed) {
    Parameters* par = new Parameters();
    par->define_defaults();

//  Posterior pars
//  0 mildEF real,
//  1 severeEF real,
//  2 base_path real,
//  3 sec_sev real,
//  4 pss_ratio real,
//  5 exp_coef real,
//  6 num_mos real,
//
//  Replicate par
//  7 realization int,
//
//  Fitted pars
//  8 foi_multiplier real,
//  9 vac_duration int,
//  10 initial_efficacy real,
//  11 CYD pseudo

    //const vector<float> foi_targets = {8,9,10,11,12,13,14,15,16,17}; // NEEDS TO BE UPDATED IF TRANSMISSION MODEL CHANGES
    //const vector<float> beta_multipliers = {1.50, 1.60, 1.70, 1.80, 1.90};
    //const vector<float> beta_multipliers = {0.52, 0.70, 0.92, 1.25, 2.25};   // from dec fit

    double _mild_EF      = args[0];
    double _severe_EF    = args[1];
    double _base_path    = args[2];
    double _sec_severity = args[3];
    double _pss_ratio    = args[4];
    double _exp_coef     = args[5];
    double _nmos         = args[6];
    double _betamp       = 0.25;//foi_targets[(int) args[7]] / _nmos;
    double _betapm       = 0.10;//_betamp;
    par->reportedFraction = {0.0, 1.0/_mild_EF, 1.0/_severe_EF}; // no asymptomatic infections are reported

    vector<long double> abc_args(&args[0], &args[6]);
    string argstring;
    const string process_id = to_string(calculate_process_id(abc_args, argstring));

    par->randomseed              = rng_seed;
    par->dailyOutput             = false;
    par->monthlyOutput           = false;
    par->yearlyOutput            = true;
    par->abcVerbose              = false; // needs to be false to get WHO daily output
    const int runLengthYears     = RESTART_BURNIN + FORECAST_DURATION;
    par->nRunLength              = runLengthYears*365;
    par->startDayOfYear          = 100;
    par->annualIntroductionsCoef = pow(10, _exp_coef);

    // pathogenicity values fitted in
    // Reich et al, Interactions between serotypes of dengue highlight epidemiological impact of cross-immunity, Interface, 2013
    // Normalized from Fc values in supplement table 2, available at
    // http://rsif.royalsocietypublishing.org/content/10/86/20130414/suppl/DC1
//par->defineSerotypeRelativeRisks();
par->serotypePathogenicityRelativeRisks = vector<double>(NUM_OF_SEROTYPES, 1.0);
    par->useAgeStructuredPrimaryPathogenicity = false;
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
    par->immunityFilename         = "";
    par->locationFilename         = pop_dir    + "/locations-"          + SIM_POP + ".txt";
    par->networkFilename          = pop_dir    + "/network-"            + SIM_POP + ".txt";
    par->swapProbFilename         = pop_dir    + "/swap_probabilities-" + SIM_POP + ".txt";
    //par->mosquitoFilename         = output_dir + "/mos_mer_who/mos."       + to_string(process_id);
    //par->mosquitoLocationFilename = output_dir + "/mosloc_mer_who/mosloc." + to_string(process_id);

    return par;
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


const unsigned int report_process_id (vector<long double> &args, const ABC::MPI_par* mp, const time_t start_time) {
    double dif = difftime (start_time, GLOBAL_START_TIME);

    string argstring;
    const unsigned int process_id = calculate_process_id(args, argstring);

    stringstream ss;
    ss << mp->mpi_rank << " begin " << hex << process_id << " " << dec << dif << " " << start_time << " " << argstring << endl;
    string output = ss.str();
    fputs(output.c_str(), stderr);

    return process_id;
}
    

vector<long double> tally_counts(Parameters* par, Community* community) {
    const int VAC_DAY = RESTART_BURNIN*365;
    vector<long double> metrics;

    vector< vector<int> > age_bins = {{2,5},{6,11},{12,14}};
    vector< vector<double> > infection_attack_rate_age  = {{0,0},{0,0},{0,0}}; 
    vector< vector<double> > infection_attack_rate_sero = {{0,0},{0,0}}; 

    vector<int> pop_sizes                            = {0,0,0};
    vector< vector<int> > arm_pop_sizes              = {{0,0},{0,0},{0,0}}; // age; placebo = 0, vaccine = 1
    vector< vector<int> > sero_pop_sizes             = {{0,0},{0,0}};       // seroneg = 0, seropos = 1; placebo = 0, vaccine = 1
    vector<double> seroneg                           = {0,0,0}; // seronegative at vaccination time
    vector< vector<double> > attack_rate_age         = {{0,0},{0,0},{0,0}};
    vector< vector<double> > attack_rate_sero        = {{0,0},{0,0}};
    vector< vector<double> > mild_attack_rate_hosp   = {{0,0},{0,0},{0,0}};
    vector< vector<double> > severe_attack_rate_hosp = {{0,0},{0,0},{0,0}};

    for (int i=community->getNumPerson()-1; i>=0; i--) {
        Person *p = community->getPerson(i);
        const int final_age  = p->getAge();
        const int age_at_vac = final_age - FORECAST_DURATION;
        if (age_at_vac >= age_bins.front()[0] and age_at_vac <= age_bins.back()[1]) {
            int age_class = 0;
            for (auto bin: age_bins) {
                if (age_at_vac > bin[1]) ++age_class;
            }
            const bool ever_infected = p->getNumInfections() > 0 ? true : false;
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
                } else if (itime >= VAC_DAY + (365*2) + 28 and itime < VAC_DAY + (365*3) + 28) {
                    // hospital phase
                    if ((const bool) infec->isSevere()) {
                        ++severe_attack_rate_hosp[age_class][vaccinated]; 
                    } else {
                        ++mild_attack_rate_hosp[age_class][vaccinated]; 
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
    //vector< vector<double> > mild_attack_rate_hosp   = {{0,0},{0,0},{0,0}};
    //vector< vector<double> > severe_attack_rate_hosp = {{0,0},{0,0},{0,0}};

    for (unsigned int i = 0; i < seroneg.size(); ++i) metrics.push_back(seroneg[i] / pop_sizes[i]);
    for (unsigned int i = 0; i < attack_rate_age.size(); ++i) {
        for (unsigned int j = 0; j < attack_rate_age[i].size(); ++j) metrics.push_back(attack_rate_age[i][j] / arm_pop_sizes[i][j]);
    }
    for (unsigned int i = 0; i < attack_rate_sero.size(); ++i) {
        for (unsigned int j = 0; j < attack_rate_sero[i].size(); ++j) metrics.push_back(attack_rate_sero[i][j] / sero_pop_sizes[i][j]);
    }
    for (unsigned int i = 0; i < mild_attack_rate_hosp.size(); ++i) {
        for (unsigned int j = 0; j < mild_attack_rate_hosp[i].size(); ++j) metrics.push_back(mild_attack_rate_hosp[i][j] / arm_pop_sizes[i][j]);
    }
    for (unsigned int i = 0; i < severe_attack_rate_hosp.size(); ++i) {
        for (unsigned int j = 0; j < severe_attack_rate_hosp[i].size(); ++j) metrics.push_back(severe_attack_rate_hosp[i][j] / arm_pop_sizes[i][j]);
    }

    for (unsigned int i = 0; i < infection_attack_rate_age.size(); ++i) {
        for (unsigned int j = 0; j < infection_attack_rate_age[i].size(); ++j) metrics.push_back(infection_attack_rate_age[i][j] / arm_pop_sizes[i][j]);
    }
    for (unsigned int i = 0; i < infection_attack_rate_sero.size(); ++i) {
        for (unsigned int j = 0; j < infection_attack_rate_sero[i].size(); ++j) metrics.push_back(infection_attack_rate_sero[i][j] / sero_pop_sizes[i][j]);
    }

    return metrics;
}


vector<long double> simulator(vector<long double> args, const unsigned long int rng_seed, const ABC::MPI_par* mp) {
    gsl_rng_set(RNG, rng_seed); // seed the rng using sys time and the process id

    //initialize bookkeeping for run
    time_t start, end;
    time (&start);
    const unsigned int process_id = report_process_id(args, mp, start);

    cerr << "SCENARIO " << rng_seed;
    for (auto _p: args) cerr << " " << _p; cerr << endl;

    Parameters* par = define_simulator_parameters(args, rng_seed);

    //vector<int> vaccine_durations = {-1, 10*365};
//  Posterior pars
//  0 mildEF real,
//  1 severeEF real,
//  2 base_path real,
//  3 sec_sev real,
//  4 pss_ratio real,
//  5 exp_coef real,
//  6 num_mos real,
//
//  Replicate par
//  7 realization int,
//
//  Fitted pars
//  8 foi_multiplier real,
//  9 vac_duration int,
//  10 initial_efficacy real,
//  11 CYD pseudo


    int youngest    = 2;
    int oldest      = 14;
    double coverage = 0.07; // approximately what is necessary to get CYD14-like numbers of vaccinees

    Community* community = build_community(par);

    // perfect efficacy that wanes rapidly
    par->fVESs                   = vector<double>(NUM_OF_SEROTYPES, args[10]);
    par->fVESs_NAIVE             = vector<double>(NUM_OF_SEROTYPES, args[10]);
    par->fVEH                    = 0.803; // fraction of hospitalized cases prevented by vaccine
    par->linearlyWaningVaccine   = true;
    par->vaccineImmunityDuration = (int) args[9];
    par->bVaccineLeaky           = true;
    par->numVaccineDoses         = 3;
    par->vaccineDoseInterval     = 182;
    par->vaccineBoosting         = false;

    for (int vac_age = youngest; vac_age <= oldest; ++vac_age) {
        par->vaccinationEvents.emplace_back(vac_age, RESTART_BURNIN*365, coverage);
    }

    par->whoDiseaseOutcome = INC_NUM_INFECTIONS;
    par->whoBreakthrough   = BREAKTHROUGH_SEROCONVERSION;
    par->whoWaning         = UNIVERSAL_WANING; // Naive-only was specified, but seems unrealistic

    seed_epidemic(par, community);
    simulate_epidemic(par, community, process_id);

    time (&end);
    double dif = difftime (end,start);

    vector<long double> metrics = tally_counts(par, community);

    stringstream ss;
    ss << mp->mpi_rank << " end " << hex << process_id << " " << dec << dif << " ";

    for (auto i: args) ss << i << " ";
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
