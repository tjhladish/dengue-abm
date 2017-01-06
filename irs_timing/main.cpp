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
const string SIM_POP = "yucatan";
const string HOME_DIR(std::getenv("HOME"));
const string pop_dir = HOME_DIR + "/work/dengue/pop-" + SIM_POP;
const string output_dir("/ufrc/longini/tjhladish/");
const string imm_dir(output_dir + "/imm_1000_yucatan");

const int RESTART_BURNIN    = 25;
const int FORECAST_DURATION = 51;
const bool RUN_FORECAST     = true;
const int TOTAL_DURATION    = RUN_FORECAST ? RESTART_BURNIN + FORECAST_DURATION : RESTART_BURNIN;

//Parameters* define_simulator_parameters(vector<double> args, const unsigned long int rng_seed) {
Parameters* define_simulator_parameters(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const string process_id) {
    Parameters* par = new Parameters();
    par->define_defaults();
    par->serial = serial;

    //00  "mild_expansion_factor",
    //01  "severe_expansion_factor",
    //02  "base_pathogenicity",
    //03  "secondary_severity",
    //04  "primary_secondary_severity_ratio",
    //05  "exposures_coefficient",
    //06  "mosquito_density",
    //07  "realization",
    //08  "vector_control",
    //09  "campaign_duration",
    //10  "timing",
    //11  "vc_coverage",
    //12  "vc_efficacy",

    double _mild_EF      = args[0];
    double _severe_EF    = args[1];
    double _base_path    = args[2];
    double _sec_severity = args[3];
    double _pss_ratio    = args[4];
    double _exp_coef     = args[5];
    double _nmos         = args[6];
    double _betamp       = 0.25; // beta values from chao et al
    double _betapm       = 0.10; //
    par->reportedFraction = {0.0, 1.0/_mild_EF, 1.0/_severe_EF}; // no asymptomatic infections are reported

    par->randomseed              = rng_seed;
    par->dailyOutput             = false;
    par->periodicOutput          = false;
    par->periodicOutputInterval  = 5;
    par->weeklyOutput            = true;
    par->monthlyOutput           = false;
    par->yearlyOutput            = true;
    par->abcVerbose              = false; // needs to be false to get WHO daily output
    const int runLengthYears     = TOTAL_DURATION;
    par->nRunLength              = runLengthYears*365;
    par->startDayOfYear          = 100;
    par->birthdayInterval        = 365; // process birthdays daily
    par->delayBirthdayIfInfected = false;
    par->annualIntroductionsCoef = pow(10, _exp_coef);

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
    par->nDefaultMosquitoCapacity = (int) _nmos;
    par->eMosquitoDistribution = EXPONENTIAL;

    par->nDaysImmune = 730;
    par->fVESs = vector<double>(NUM_OF_SEROTYPES, 0.0);

    //par->hospitalizedFraction = 0.25; // fraction of cases assumed to be hospitalized

    par->simulateAnnualSerotypes = false;
    //par->normalizeSerotypeIntros = true;
    //if (par->simulateAnnualSerotypes) par->generateAnnualSerotypes();
    par->nDailyExposed = {{0.25, 0.25, 0.25, 0.25}};

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
    par->immunityFilename         = imm_dir    + "/immunity2003."       + process_id;
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


vector<double> tally_counts(const Parameters* par, Community* community, const int vc_timing) {
    // aggregate based on the timing of the annual start of vector control
    int discard_days = INT_MAX;
    for (VectorControlEvent vce: par->vectorControlEvents) {
        discard_days = discard_days < vce.campaignStart ? discard_days : vce.campaignStart;
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

    // initialize bookkeeping for run
    time_t start, end;
    time (&start);
    //const string process_id = report_process_id(args, serial, mp, start);
    vector<double> abc_args(&args[0], &args[7]);
    const unsigned int realization = (int) args[7];

    const string process_id = report_process_id(abc_args, serial, mp, start) + "." + to_string(realization);
    report_process_id(args, serial, mp, start);

    cerr << "SCENARIO " << rng_seed;
    for (auto _p: args) cerr << " " << _p; cerr << endl;

    Parameters* par = define_simulator_parameters(args, rng_seed, serial, process_id);
//  Posterior pars
    //00  "mild_expansion_factor",
    //01  "severe_expansion_factor",
    //02  "base_pathogenicity",
    //03  "secondary_severity",
    //04  "primary_secondary_severity_ratio",
    //05  "exposures_coefficient",
    //06  "mosquito_density",

//  Pseudo pars
    //07  "realization",
    //08  "vector_control",
    //09  "campaign_duration",
    //10  "timing",
    //11  "vc_coverage",
    //12  "vc_efficacy",

    vector<int> vc_campaign_duration_levels = {1, 90, 365};

    const bool vector_control     = (bool) args[8];
    const int vc_campaignDuration = vc_campaign_duration_levels[(int) args[9]]; // number of days to achieve coverage
    const int vc_timing           = (int) args[10]; 
    const double vc_coverage      = args[11];
    const double vc_efficacy      = args[12];       // expected % reduction in equillibrium mosquito population in treated houses

    Community* community = build_community(par);

    if (vector_control) {
        const int efficacyDuration = 90;       // number of days efficacy is maintained
        const LocationType locType = HOME;
        const LocationSelectionStrategy lss = UNIFORM_STRATEGY;
        for (int vec_cont_year = RESTART_BURNIN; vec_cont_year < TOTAL_DURATION; vec_cont_year++) {
        //for (int vec_cont_year = 0; vec_cont_year < TOTAL_DURATION; vec_cont_year++) {
            // TODO - address situation where startDate could be negative if startDayOfYear is larger than vc_timing
            int startDate = 0;
            if (vc_timing >= par->startDayOfYear) { // start before Jan 1
                startDate = vec_cont_year*365 + vc_timing - par->startDayOfYear; // 151 days after Jan 1 is June 1, offset by julian startDay
            } else {                                // start after Jan 1
                startDate = (vec_cont_year+1)*365 + vc_timing - par->startDayOfYear;
            }
            par->vectorControlEvents.emplace_back(startDate, vc_campaignDuration, vc_coverage, vc_efficacy, efficacyDuration, locType, lss);
        }
    }

    seed_epidemic(par, community);
    simulate_epidemic(par, community, process_id);

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
