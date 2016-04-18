#include <unistd.h>
#include "AbcSmc.h"
#include "simulator.h"
#include <cstdlib>
#include "CCRC32.h"
#include "Utility.h"
#include "yucatan_serotype_generator.h"
#include <unordered_set>

using namespace std;

using dengue::util::to_string;
using dengue::util::mean;
using dengue::util::stdev;
using dengue::util::max_element;
using ABC::float_type;

time_t GLOBAL_START_TIME;
/* Notes
 *
 * - Revert code to resuming simulation circa 2003, rather than doing entire time period
 * - Add realization pseudo parameter to follow posterior parameters, with values 0-9
 * - Load serotype sequences, throw out initial, unsimulated years
 * - load immunity
 * - filenames for both should be constructed using proccess_id based on epi pars, and the realization #
 *
 *
 */

string calculate_process_id(vector< long double> &args, string &argstring);
const string SIM_POP = "yucatan";

const int FIRST_YEAR          = 1879;                                 // inclusive
const int FIRST_OBSERVED_YEAR = 1979;
const int LAST_YEAR           = 2013;                                 // inclusive
const int DDT_START           = 1956 - FIRST_YEAR;                    // simulator year 77, counting from year 1
const int DDT_DURATION        = 23;                                   // 23 years long
const int HISTORICAL_DURATION = LAST_YEAR - FIRST_YEAR + 1;
const int FITTED_DURATION     = LAST_YEAR - FIRST_OBSERVED_YEAR + 1;  // 35 for 1979-2013 inclusive
const int FORECAST_DURATION   = 20;
const int BURNIN_YEARS        = 125;                                  // should be zero if we're running from FIRST_YEAR; 125 should correspond to resuming after 2003
const int TOTAL_DURATION      = HISTORICAL_DURATION + FORECAST_DURATION - BURNIN_YEARS;

const bool CLIMATE_CHANGE     = true;

Parameters* define_simulator_parameters(vector<long double> args, const unsigned long int rng_seed, const unsigned long int serial, const string process_id) {
    Parameters* par = new Parameters();
    par->define_defaults();
    par->serial = serial;

    double _mild_EF      = args[0];
    double _severe_EF    = args[1];
    double _base_path    = args[2];
    double _sec_severity = args[3];
    double _pss_ratio    = args[4];
    double _exp_coef     = args[5];
    double _nmos         = args[6];
    double _betamp       = 0.25; // beta values from chao et al
    double _betapm       = 0.10; //
    //double _flav_ar      = args[8];

    par->reportedFraction = {0.0, 1.0/_mild_EF, 1.0/_severe_EF}; // no asymptomatic infections are reported

    string HOME(std::getenv("HOME"));
    string pop_dir = HOME + "/work/dengue/pop-" + SIM_POP;
    string output_dir = "/scratch/lfs/thladish";
    string imm_dir(output_dir + "/imm_1000_yucatan");

    par->randomseed              = rng_seed;
    par->dailyOutput             = false;
    par->monthlyOutput           = true;
    par->yearlyOutput            = true;
    par->abcVerbose              = true;
    //const int runLengthYears     = TOTAL_DURATION;
    par->nRunLength              = TOTAL_DURATION*365;
    par->startDayOfYear          = 100;
    //if ( SIM_POP == "merida") {
    par->annualIntroductionsCoef = pow(10, _exp_coef);
    //} else if ( SIM_POP == "yucatan" ) {
    //    // assuming parameter fitting was done on merida population 
    //    par->annualIntroductionsCoef = pow(10, _exp_coef)*1819498.0/839660.0;
    //} else {
    //    cerr << "ERROR: Unknown simulation population (SIM_POP): " << SIM_POP << endl;
    //    exit(523);
    //}

    // pathogenicity values fitted in
    // Reich et al, Interactions between serotypes of dengue highlight epidemiological impact of cross-immunity, Interface, 2013
    // Normalized from Fc values in supplement table 2, available at
    // http://rsif.royalsocietypublishing.org/content/10/86/20130414/suppl/DC1
    par->defineSerotypeRelativeRisks();
    par->basePathogenicity = _base_path;
    par->primaryPathogenicityModel = ORIGINAL_LOGISTIC;
    //par->annualFlavivirusAttackRate = _flav_ar;
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
    par->fVESs.clear();
    par->fVESs.resize(NUM_OF_SEROTYPES, 0);

    // generate introductions of serotypes based on when they were first observed in Yucatan
//    par->nDailyExposed = generate_serotype_sequences(RNG, FIRST_YEAR, FIRST_OBSERVED_YEAR, LAST_YEAR+FORECAST_DURATION, TRANS_AND_NORM);

    par->simulateAnnualSerotypes = false;
    //par->normalizeSerotypeIntros = true;
    // generate some extra years of serotypes, for subsequent intervention modeling
    //if (par->simulateAnnualSerotypes) par->generateAnnualSerotypes(runLengthYears+50);
    // 77 year burn-in, 23 years of no dengue, then re-introduction
    // annualIntros is indexed in terms of simulator (not calendar) years
    string sero_filename = output_dir + "/sero/annual_serotypes." + process_id;
    par->loadAnnualSerotypes(sero_filename);
    vector<vector<float> >(par->nDailyExposed.begin()+BURNIN_YEARS, par->nDailyExposed.end()).swap(par->nDailyExposed);

    par->annualIntroductions = vector<double>(1, 1.0);
    //par->annualIntroductions.resize(DDT_START+DDT_DURATION, 0.1); // 90% reduction in intros for 20 years
    //par->annualIntroductions.resize(HISTORICAL_DURATION+FORECAST_DURATION, 1.0);

    // load daily EIP
    // EIPs calculated by ../raw_data/weather/calculate_daily_eip.R, based on reconstructed Merida temps
    // we add the startDayOfYear offset because we will end up discarding that many values from the beginning
    // mosquitoMultipliers are indexed to start on Jan 1
    if (CLIMATE_CHANGE) {
        par->loadDailyEIP(pop_dir + "/seasonal_EIP_24hr.out", (HISTORICAL_DURATION-BURNIN_YEARS)*365 + par->startDayOfYear);
        string dailyEIPfilename = pop_dir + "/forecast_EIP_0.02degC_per_year.out";
        vector<string> fitted_EIPs = dengue::util::read_vector_file(dailyEIPfilename); // slurp first column from file

        const int duration = 1;
        int start = par->extrinsicIncubationPeriods.back().start;
        for(unsigned int i = par->startDayOfYear; i < fitted_EIPs.size(); ++i) {
            string val_str = fitted_EIPs[i];
            const double value = dengue::util::to_double(val_str);
            start += duration;
            if (value <= 0) {
                cerr << "An EIP <= 0 was read from " << dailyEIPfilename << "." << endl;
                cerr << "Value read: " << value << endl;
                cerr << "This is nonsensical and indicates a non-numerical value in the first column or an actual bad value." << endl;
                exit(113);
            }
            par->extrinsicIncubationPeriods.emplace_back(start, duration, value);
        }
    } else {
        par->loadDailyEIP(pop_dir + "/seasonal_EIP_24hr.out");
    }
    par->loadDailyMosquitoMultipliers(pop_dir + "/mosquito_seasonality.out");
    //par->loadDailyMosquitoMultipliers(pop_dir + "/mosquito_seasonality.out", (365*TOTAL_DURATION) + par->startDayOfYear);
    //assert(par->mosquitoMultipliers.size() >= (DDT_START+DDT_DURATION)*365);
    //for (int i = DDT_START*365; i < (DDT_START+DDT_DURATION)*365; ++i) par->mosquitoMultipliers[i].value *= 0.23; // 77% reduction in mosquitoes
    //vector<DynamicParameter>(par->mosquitoMultipliers.begin()+BURNIN_YEARS, par->mosquitoMultipliers.end()).swap(par->mosquitoMultipliers);

    par->populationFilename       = pop_dir    + "/population-"         + SIM_POP + ".txt";
    // immunity file will be loaded when build community is called if string != ""
    par->immunityFilename         = imm_dir    + "/immunity2003."       + process_id;
    //par->immunityFilename         = "";
    par->locationFilename         = pop_dir    + "/locations-"          + SIM_POP + ".txt";
    par->networkFilename          = pop_dir    + "/network-"            + SIM_POP + ".txt";
    par->swapProbFilename         = pop_dir    + "/swap_probabilities-" + SIM_POP + ".txt";
    //par->mosquitoFilename         = output_dir + "/mos_tmp/mos."       + to_string(process_id);
    //par->mosquitoLocationFilename = output_dir + "/mosloc_tmp/mosloc." + to_string(process_id);
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


string calculate_process_id(vector< long double> &args, string &argstring) {
    // CCRC32 checksum based on string version of argument values
    CCRC32 crc32;
    crc32.Initialize();

    for (unsigned int i = 0; i < args.size(); i++) argstring += to_string((long double) args[i]) + " ";

    const unsigned char* argchars = reinterpret_cast<const unsigned char*> (argstring.c_str());
    const int len = argstring.length();
    const int process_id = crc32.FullCRC(argchars, len);

    return to_string(process_id);
}


string report_process_id (vector<long double> &args, const unsigned long int serial, const ABC::MPI_par* mp, const time_t start_time) {
    double dif = difftime (start_time, GLOBAL_START_TIME);

    string argstring;
    const string process_id = calculate_process_id(args, argstring);

    stringstream ss;
    ss << mp->mpi_rank << " begin " << process_id << " " << dec << serial << " " << dif << " " << argstring << endl;
    string output = ss.str();
    fputs(output.c_str(), stderr);

    return to_string(process_id);
}


void append_if_finite(vector<long double> &vec, double val) {
    if (isfinite(val)) { 
        vec.push_back((long double) val);
    } else {
        vec.push_back(0);
    }
}


vector<long double> tally_counts(const Parameters* par, Community* community, const int discard_days, const int output_years) {
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                                  IMPORTANT:                                           //
    // Update dummy metrics vector in calling function if number of metrics is changed here! //
    ///////////////////////////////////////////////////////////////////////////////////////////
    vector< vector<int> > severe      = community->getNumSevereCases();
    vector< vector<int> > symptomatic = community->getNumNewlySymptomatic();
    vector< vector<int> > infected    = community->getNumNewlyInfected();
    //const int num_years = (int) par->nRunLength/365;
    vector<vector<int> > severe_tally(NUM_OF_SEROTYPES, vector<int>(output_years, 0)); // +1 to handle run lengths of a non-integral number of years
    vector<vector<int> > symptomatic_tally(NUM_OF_SEROTYPES, vector<int>(output_years, 0)); // (any extra fraction of a year will be discarded)
    vector<vector<int> > infected_tally(NUM_OF_SEROTYPES, vector<int>(output_years, 0));

    vector<long double> metrics;
    for (int t=discard_days; t<par->nRunLength; t++) {
        // use epidemic years, instead of calendar years
        const int y = (t-discard_days)/365;
        for (int s=0; s<NUM_OF_SEROTYPES; s++) {
            severe_tally[s][y]      += severe[s][t];
            symptomatic_tally[s][y] += symptomatic[s][t];
            infected_tally[s][y]    += infected[s][t];
        }
    }

    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = 0; y<output_years; ++y) {
            metrics.push_back(infected_tally[s][y]);
        }
    }
    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = 0; y<output_years; ++y) {
            metrics.push_back(symptomatic_tally[s][y] - severe_tally[s][y]);
        }
    }
    for (int s=0; s<NUM_OF_SEROTYPES; s++) {
        for (int y = 0; y<output_years; ++y) {
            metrics.push_back(severe_tally[s][y]);
        }
    }
    return metrics;
}

vector<long double> simulator(vector<long double> args, const unsigned long int rng_seed, const unsigned long int serial, const ABC::MPI_par* mp) {
    gsl_rng_set(RNG, rng_seed); // seed the rng using sys time and the process id
    // initialize bookkeeping for run
    time_t start ,end;
    time (&start);

    vector<long double> abc_args(&args[0], &args[7]);
    const unsigned int realization = (int) args[7];
    const string process_id = report_process_id(abc_args, serial, mp, start) + "." + to_string(realization);

    // initialize & run simulator
    Parameters* par = define_simulator_parameters(args, rng_seed, serial, process_id);

    const int output_years        = FORECAST_DURATION + 5;           // needs to agree with json file
    const int burnin_years        = 10;
    const int discard_days        = 5*365; //(HISTORICAL_DURATION - 5) * 365;
    const int vac_start_year      = burnin_years; //HISTORICAL_DURATION;
    vector<float> vector_controls = {0.0, 0.1, 0.25, 0.5};
    vector<int> vaccine_durations = {-1, 4*365, 10*365, 20*365};

//  0 mild_expansion_factor
//  1 severe_expansion_factor
//  2 base_pathogenicity
//  3 secondary_severity
//  4 primary_secondary_severity_ratio
//  5 exposures_coefficient
//  6 mosquito_density
//  7 realization
//
//  8 vac
//  9 catchup
// 10 target
// 11 catchup_to
// 12 vec_control
// 13 vac_waning
// 14 vac_boosting
// 15 vec_scenario

    bool vaccine           = (bool) args[8];
    bool catchup           = (bool) args[9];
    int target             = (int) args[10];
    int catchup_to         = (int) args[11];
    float vector_reduction = vector_controls[(unsigned int) args[12]];
    int vaccine_duration   = vaccine_durations[(unsigned int) args[13]];
    bool boosting          = (bool) args[14];
    int vector_scenario    = (int) args[15];

    // given current par combinations, these are no-intervention runs
    // Daily output can be used for e.g. avg epi-curve plot
    //par->dailyOutput = ((serial - 4) % 1152 == 0);
    par->dailyOutput = not vaccine;

    if (vaccine_duration == -1) {
        par->linearlyWaningVaccine = false;
        par->vaccineImmunityDuration = INT_MAX;
    } else {
        par->linearlyWaningVaccine = true;
        par->vaccineImmunityDuration = vaccine_duration;
    }

    bool nonsensical_parameters = false;
    // only run a non-vaccination campaign if all the vaccine parameters are 0
    // TODO - this should be reworked with new "catchup-to" parameter
    if (not vaccine and (catchup or boosting)) { nonsensical_parameters = true; }
    if ((vaccine_duration == -1) and boosting) { nonsensical_parameters = true; }
    if (nonsensical_parameters) {
        // 3 is because vaccinated cases, total cases, and infections are reported
        vector<long double> dummy(output_years*NUM_OF_SEROTYPES*3, 0.0);
        delete par;
        return dummy;
    }

    Community* community = build_community(par);
    //community->loadMosquitoes(par->mosquitoLocationFilename, par->mosquitoFilename);

    if (vaccine) {
        double default_target_coverage = 0.8;
        double default_catchup_coverage = 0.6;
        int num_target  = community->ageIntervalSize(9,10); // default target
        int num_catchup = community->ageIntervalSize(10,31);// default catchup

        double target_coverage  = default_target_coverage*num_target/community->ageIntervalSize(target, target+1);
        double catchup_coverage = default_catchup_coverage*num_catchup/community->ageIntervalSize(target+1, catchup_to+1);

        target_coverage = target_coverage > 1.0 ? 1.0 : target_coverage;
        catchup_coverage = catchup_coverage > 1.0 ? 1.0 : catchup_coverage;

        par->bVaccineLeaky       = true;
        par->fVEH                = 0.803; // fraction of hospitalized cases prevented by vaccine
        par->numVaccineDoses     = 3;
        par->vaccineDoseInterval = 182;
        par->vaccineBoosting     = boosting;
        par->fVESs               = {0.6, 0.54, 0.9, 0.95};
        par->fVESs_NAIVE         = {0.3, 0.27, 0.45, 0.48};

        if (catchup) {
            for (int catchup_age = target + 1; catchup_age <= catchup_to; catchup_age++) {
                par->vaccinationEvents.emplace_back(catchup_age, vac_start_year*365, catchup_coverage);
            }
        } 

        for (int vacc_year = vac_start_year; vacc_year < TOTAL_DURATION; vacc_year++) {
            par->vaccinationEvents.emplace_back(target, vacc_year*365, target_coverage);
        }
    }

    int vr_start_day = 0;
    int vr_end_day = 0;

    switch (vector_scenario) {
        case 0:
            vr_start_day = vac_start_year*365;
            vr_end_day   = TOTAL_DURATION*365;
            break;
        case 1:
            vr_start_day = vac_start_year*365;
            vr_end_day   = (vac_start_year+10)*365;
            break;
        case 2:
            vr_start_day = (vac_start_year-5)*365;
            vr_end_day   = vac_start_year*365;
            break;
        default:
            cerr << "Unsupported vector scenario: " << vector_scenario << endl;
            exit(298);
    }

    if (vector_reduction > 0.0) {
        assert( vector_reduction <= 1.0); 
        assert( vr_end_day <= (signed) par->mosquitoMultipliers.size() );
        const float vector_coeff = 1.0 - vector_reduction; // normally, there's no reduction in vectors
        for (int i = vr_start_day; i < vr_end_day; ++i) {
            par->mosquitoMultipliers[i].value *= vector_coeff;
        }
    }

    seed_epidemic(par, community);
    simulate_epidemic(par, community, process_id);
    //vector<int> epi_sizes = simulate_epidemic(par, community, process_id);

    time (&end);
    double dif = difftime (end,start);

    vector<long double> metrics = tally_counts(par, community, discard_days, output_years);

    stringstream ss;
    ss << mp->mpi_rank << " end " << hex << process_id << " " << dec << serial << " " << dif << " ";

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
