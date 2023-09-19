#include <unistd.h>
#include <AbcSmc/AbcSmc.h>
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
const string SIM_POP = "yucatan";
const string HOME_DIR(std::getenv("HOME"));
const string pop_dir = HOME_DIR + "/work/dengue/pop-" + SIM_POP;
const string output_dir("/ufrc/longini/tjhladish/");
const string imm_dir(output_dir + "imm_1000_yucatan-irs_refit3");

const int RESTART_BURNIN    = 10; // 10 for normal results, 100 for foi-effect analysis if changing foi from imm file
const int FORECAST_DURATION = 21; // normally 41; using 11 for efficacy duration sensitivity analysis
const bool RUN_FORECAST     = true;
const int TOTAL_DURATION    = RUN_FORECAST ? RESTART_BURNIN + FORECAST_DURATION : RESTART_BURNIN;
const int JULIAN_TALLY_DATE = 146; // intervention julian date - 1

//Parameters* define_simulator_parameters(vector<double> args, const unsigned long int rng_seed) {
Parameters* define_simulator_parameters(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const string process_id) {
    Parameters* par = new Parameters();
    par->define_defaults();
    par->serial = serial;

    //00   "mild_rf",
    //01   "p95_mrf",
    //02   "severe_rf",
    //03   "sec_path",
    //04   "sec_sev",
    //05   "pss_ratio",
    //06   "exp_coef",
    //07   "num_mos",

    //08   "realization",
    //09   "vector_control",
    //10   "campaign_duration",
    //11   "timing",
    //12   "vc_coverage",
    //13   "vc_efficacy",

    double _mild_RF      = args[0];
  //double _p95_mild_RF  = args[1]; // not needed in this scope
    double _severe_RF    = args[2];
    double _base_path    = args[3];
    double _sec_severity = args[4];
    double _pss_ratio    = args[5];
    double _exp_coef     = args[6];
    double _nmos         = args[7];
    double _foi_mult     = args[14];
    double _betamp       = 0.25; // beta values from chao et al
    double _betapm       = 0.10; //

    par->reportedFraction = {0.0, _mild_RF, _severe_RF}; // no asymptomatic infections are reported

    par->randomseed              = rng_seed;
    par->dailyOutput             = false; // turn on for daily prevalence figure, probably uncomment filter in simulator.h for daily output to get only rel. days
    par->periodicOutput          = false;
    par->periodicOutputInterval  = 5;
    par->weeklyOutput            = false;
    par->monthlyOutput           = false;
    par->yearlyOutput            = true;
    par->abcVerbose              = false; // needs to be false to get WHO daily output
    const int runLengthYears     = TOTAL_DURATION;
    par->nRunLength              = runLengthYears*365;
    par->startDayOfYear          = 1;
    par->birthdayInterval        = 1;
    par->delayBirthdayIfInfected = false;
    par->annualIntroductionsCoef = _exp_coef;
    //par->annualIntroductionsCoef = pow(10, _exp_coef);

    // pathogenicity values fitted in
    // Reich et al, Interactions between serotypes of dengue highlight epidemiological impact of cross-immunity, Interface, 2013
    // Normalized from Fc values in supplement table 2, available at
    // http://rsif.royalsocietypublishing.org/content/10/86/20130414/suppl/DC1
    par->defineSerotypeRelativeRisks();
//par->serotypePathogenicityRelativeRisks = vector<double>(NUM_OF_SEROTYPES, 1.0);
    //par->useAgeStructuredPrimaryPathogenicity = false;
    par->basePathogenicity = _base_path;
    par->primaryPathogenicityModel = ORIGINAL_LOGISTIC;
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
    par->nDefaultMosquitoCapacity = (int) (_nmos * _foi_mult);
    par->eMosquitoDistribution = EXPONENTIAL;

    par->nDaysImmune = 730;
    par->fVESs = vector<double>(NUM_OF_SEROTYPES, 0.0);

    //par->hospitalizedFraction = 0.25; // fraction of cases assumed to be hospitalized

    par->simulateAnnualSerotypes = false;
    //par->normalizeSerotypeIntros = true;
    //if (par->simulateAnnualSerotypes) par->generateAnnualSerotypes();
    par->nDailyExposed = {{1.0, 1.0, 1.0, 1.0}};
    //par->nDailyExposed = {{0.25, 0.25, 0.25, 0.25}};

    par->annualIntroductions = {1.0};

    // load daily EIP
    // EIPs calculated by ../raw_data/weather/calculate_daily_eip.R, based on reconstructed yucatan temps
    par->loadDailyEIP(pop_dir + "/seasonal_EIP_24hr.out");

    // we add the startDayOfYear offset because we will end up discarding that many values from the beginning
    // mosquitoMultipliers are indexed to start on Jan 1
    //par->loadDailyMosquitoMultipliers(pop_dir + "/mosquito_seasonality.out", par->nRunLength + par->startDayOfYear);
    par->loadDailyMosquitoMultipliers(pop_dir + "/mosquito_seasonality.out");

    par->populationFilename       = pop_dir    + "/population-"         + SIM_POP + ".txt";
    //par->immunityFilename         = "/ufrc/longini/tjhladish/imm_who-baseline-seroprev-july2016/immunity." + imm_file_pid;
    par->immunityFilename         = imm_dir    + "/immunity2030."       + process_id;
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

    cerr << "pid in report_process_id (num args = " << args.size() << "): " << process_id << endl;
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

int julian_to_sim_day (const Parameters* par, const int julian, const int intervention_year) {
    int startDate = intervention_year*365 + julian - par->startDayOfYear;
    if (julian < par->startDayOfYear) { // start intervention in following year
        startDate += 365;
    }
    return startDate;
}


vector<double> tally_counts(const Parameters* par, Community* community, int pre_intervention_output) {
    const int discard_days = julian_to_sim_day(par, JULIAN_TALLY_DATE, RESTART_BURNIN-pre_intervention_output);

    //vector< vector<int> > severe      = community->getNumSevereCases();
    vector< vector<int> > symptomatic = community->getNumNewlySymptomatic();

    //vector< vector<int> > infected    = community->getNumNewlyInfected();
    const int num_years = FORECAST_DURATION + pre_intervention_output - 1; // typically 55
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

    // initialize bookkeeping for run
    time_t start, end;
    time (&start);
    //const string process_id = report_process_id(args, serial, mp, start);
    vector<double> abc_args(&args[0], &args[8]);
    const unsigned int realization = 0; //(int) args[9];

    const string process_id = report_process_id(abc_args, serial, mp, start) + "." + to_string(realization);
//cerr << "IDS " << process_id << " " << serial;
//for (auto _p: args) cerr << " " << _p; cerr << endl;
//exit(-10);
    report_process_id(args, serial, mp, start);

    cerr << "SCENARIO " << rng_seed;
    for (auto _p: args) { cerr << " " << _p; } cerr << endl;

    Parameters* par = define_simulator_parameters(args, rng_seed, serial, process_id);
    //00   "mild_rf",
    //01   "p95_mrf",
    //02   "severe_rf",
    //03   "sec_path",
    //04   "sec_sev",
    //05   "pss_ratio",
    //06   "exp_coef",
    //07   "num_mos",

    //08   "realization",
    //09   "vector_control",
    //10   "campaign_duration",
    //11   "timing",
    //12   "vc_coverage",
    //13   "vc_efficacy",
    //14   "foi_multiplier"
    //15   "vc_strategy_duration"
    //16   "efficacy_duration"

    vector<int> vc_campaign_duration_levels = {1, 90, 365};

    const bool vector_control     = (bool) args[9];
    const int vc_campaignDuration = vc_campaign_duration_levels[(int) args[10]]; // number of days to achieve coverage
    const int vc_timing           = (int) args[11];
    const double vc_coverage      = args[12];
    const double vc_efficacy      = args[13];       // expected % reduction in equillibrium mosquito population in treated houses
    const double foi_mult         = args[14];
    const int vc_years            = (int) args[15];
    const int efficacyDuration    = (int) args[16]; // default was 90; number of days efficacy is maintained
    const bool vaccine            = (bool) args[17];
    const int vaccine_mechanism   = (int) args[18];
    const double coverage         = args[19]; 
    const int target              = (int) args[20];
    const bool catchup            = (bool) args[21];
    const int catchup_to          = (int) args[22];
    const double seroTestFalsePos = args[23];
    const double seroTestFalseNeg = args[24];

    if (foi_mult != 1.0) {
        par->immunityFilename = "/ufrc/longini/tjhladish/imm_1000_yucatan-alt_foi/immunity2130." + process_id + ".foi" + to_string(foi_mult);
    }
    Community* community = build_community(par);

    if (vector_control) {
        assert(vc_timing - 1 == JULIAN_TALLY_DATE);
        const LocationType locType = HOME;
        const LocationSelectionStrategy lss = UNIFORM_STRATEGY;
        for (int vec_cont_year = RESTART_BURNIN; vec_cont_year < RESTART_BURNIN + vc_years; vec_cont_year++) {
            const int startDate = julian_to_sim_day(par, vc_timing, vec_cont_year);
            par->vectorControlEvents.emplace_back(startDate, vc_campaignDuration, vc_coverage, vc_efficacy, efficacyDuration, locType, lss);
        }
    }
    if (vaccine) {
        double target_coverage  = coverage;
        double catchup_coverage = coverage;

        if (catchup) {
            for (int catchup_age = target; catchup_age <= catchup_to; catchup_age++) {
                const int vacDate = julian_to_sim_day(par, JULIAN_TALLY_DATE + 1, RESTART_BURNIN);
                par->catchupVaccinationEvents.emplace_back(catchup_age, vacDate, catchup_coverage);
            }
        } 

        par->vaccineTargetAge = target;
        par->vaccineTargetCoverage = target_coverage;
        par->vaccineTargetStartDate = julian_to_sim_day(par, JULIAN_TALLY_DATE + 1, RESTART_BURNIN);
        par->seroTestFalsePos = seroTestFalsePos;
        par->seroTestFalseNeg = seroTestFalseNeg;
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
        par->vaccineSeroConstraint   = VACCINATE_ALL_SERO_STATUSES;
        par->whoDiseaseOutcome       = INC_NUM_INFECTIONS;
    } else if (vaccine_mechanism == 1) { // also want A1 + B2 + C3a
        // imperfect efficacy that does not wane
        par->fVESs                   = vector<double>(NUM_OF_SEROTYPES, 0.7);
        par->fVESs_NAIVE             = vector<double>(NUM_OF_SEROTYPES, 0.7);
        //par->fVESs = {0.6, 0.54, 0.9, 0.95};
        //par->fVESs_NAIVE = {0.3, 0.27, 0.45, 0.48};

        par->fVEH                    = 0.803; // fraction of hospitalized cases prevented by vaccine
        par->linearlyWaningVaccine   = false;
        par->vaccineImmunityDuration = INT_MAX;
        par->bVaccineLeaky           = true;
        par->numVaccineDoses         = 1;
        par->vaccineDoseInterval     = INT_MAX;
        par->vaccineBoosting         = false;
        par->vaccineSeroConstraint   = VACCINATE_ALL_SERO_STATUSES;
        par->whoDiseaseOutcome       = VAC_ISNT_INFECTION;
    } else if (vaccine_mechanism == 2) {        // "baseline" scenario: A2b + B2 + C3a
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
        par->vaccineSeroConstraint   = VACCINATE_SEROPOSITIVE_ONLY;
        par->whoDiseaseOutcome       = INC_NUM_INFECTIONS;
    } else {
        cerr << "Unsupported vaccine mechanism: " << vaccine_mechanism << endl;
        exit(-152);
    }


    seed_epidemic(par, community);
    //simulate_epidemic(par, community, process_id);
    vector< vector<double> > sero_prev;
    bool capture_sero_prev = true;
    int sero_prev_aggregation_start_date = vc_timing; // aggregation interval is [sero_prev_aggregation_start_date, (sero_prev_aggregation_start_date + 364) % 365] on a [0,364] calendar
    simulate_epidemic_with_seroprev(par, community, process_id, capture_sero_prev, sero_prev, sero_prev_aggregation_start_date);

    time (&end);
    double dif = difftime (end,start);

    const int pre_intervention_output = 5; // years
    const int desired_intervention_output = FORECAST_DURATION - 1;
    vector<double> metrics = tally_counts(par, community, pre_intervention_output);

    assert(sero_prev.size() == 5);
    for (unsigned int i = 1; i < sero_prev.size(); ++i) assert(sero_prev[0].size() == sero_prev[i].size());
    assert(sero_prev[0].size() >= pre_intervention_output + desired_intervention_output);

    // flatten sero_prev
    for (auto sero_prev_class: sero_prev) {
        for (int year = RESTART_BURNIN-pre_intervention_output; year < RESTART_BURNIN + desired_intervention_output; ++year) {
            metrics.push_back(sero_prev_class[year]);
        }
    }

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
        time(&GLOBAL_START_TIME);
        abc->set_simulator(simulator);
        abc->simulate_next_particles(buffer_size);
    }

    return 0;
}
