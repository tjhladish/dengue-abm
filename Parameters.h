#ifndef __PARAMETERS_H
#define __PARAMETERS_H

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <bitset>
#include <assert.h>
#include <gsl/gsl_rng.h>

static const int VERSION_NUMBER_MAJOR = 1;
static const int VERSION_NUMBER_MINOR = 3;

enum Serotype {
    SEROTYPE_1,
    SEROTYPE_2,
    SEROTYPE_3,
    SEROTYPE_4,
    NUM_OF_SEROTYPES,  // Make sure this is second to last
    NULL_SEROTYPE      // Make sure this is last
};

enum MosquitoDistribution {
    CONSTANT,
    EXPONENTIAL,
    NUM_OF_DISTRIBUTIONS
};

enum TimePeriod {
    HOME_MORNING,
    WORK_DAY,
    HOME_NIGHT,
    NUM_OF_TIME_PERIODS
};

enum SexType {
    UNKNOWN,
    MALE,
    FEMALE,
    NUM_OF_SEX_TYPES
};

enum InfectionOutcome {
    ASYMPTOMATIC,
    MILD,
    SEVERE,
    NUM_OF_INFECTION_OUTCOMES
};

enum PrimaryPathogenicityModel {
    CONSTANT_PATHOGENICITY,
    ORIGINAL_LOGISTIC,
    GEOMETRIC_PATHOGENICITY,
    NUM_OF_PRIMARY_PATHOGENICITY_MODELS
};

extern const gsl_rng* RNG;// = gsl_rng_alloc (gsl_rng_taus2);

static const std::vector<std::string> MONTH_NAMES = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};
static const std::vector<int> DAYS_IN_MONTH = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
static const std::vector<int> END_DAY_OF_MONTH = {31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};

static const int SYMPTOMATIC_DELAY = 1;                       // delay of symptoms after infectious period starts
//static const int INFECTIOUS_PERIOD_PRI = 5;                 // number of days until recovery from primary infection
//static const int INFECTIOUS_PERIOD_POST_PRI = 4;            // number of days until recovery from post-primary infection
static const int INFECTIOUS_PERIOD_ASYMPTOMATIC = 2;          // number of days until recovery for asymptomatic infections
static const int INFECTIOUS_PERIOD_MILD         = 4;          // number of days until recovery for mild cases
static const int INFECTIOUS_PERIOD_SEVERE       = 6;          // number of days until recovery for severe cases

// from Person.h
static const int NUM_AGE_CLASSES = 101;                       // maximum age+1 for a person
static const int MAX_INCUBATION = 9;                          // max incubation period in days
static const int MAX_HISTORY = 50;                            // length of exposure history in years

// cdf of incubation period, starting from day 1 (from Nishiura 2007)
const std::vector<double> INCUBATION_CDF = { 0,0,0.03590193,0.5070053,0.8248687,0.9124343,0.949212,0.974606,1 };

// from Community
                                                              // probability of biting at 3 different times of day (as defined in Location.h)
static const float DAILY_BITING_PDF[(int) NUM_OF_TIME_PERIODS] = {0.08, 0.76, 0.16};

// from Mosquito
static const int MAX_MOSQUITO_AGE = 60;                       // maximum age of mosquito in days

//double MOSQUITO_DEATH_PROBABILITY[MAX_MOSQUITO_AGE] = { // probability of death each day
//    0.0018,0.002069514,0.002378646,0.002732982,0.003138823,0.003603251,0.004134191,
//    0.004740476,0.005431895,0.006219231,0.007114273,0.008129798,0.009279508,
//    0.01057792,0.01204018,0.01368180,0.01551828,0.01756467,0.01983494,0.02234135,
//    0.02509359,0.02809799,0.03135655,0.03486617,0.03861782,0.04259607,0.04677877,
//    0.05113718,0.05563645,0.06023658,0.0648937,0.0695617,0.07419404,0.07874562,
//    0.08317441,0.08744304,0.09151983,0.09537952,0.09900354,0.1023799,0.1055030,
//    0.1083723,0.1109923,0.1133714,0.1155206,0.1174531,0.1191837,0.1207277,0.1221007,
//    0.1233179,0.1243943,0.1253438,0.1261799,0.1269147,0.1275594,0.1281243,0.1286187,
//    0.1290510,0.1294285,0.1297580,0.1300453};

// cumulative density of mosquito ages
static const std::vector<double> MOSQUITO_AGE_CDF = {
    0.03115129,0.06224651,0.09327738,0.1242344,0.1551069,0.1858824,0.2165471,
    0.247085,0.2774781,0.3077061,0.3377462,0.3675725,0.3971563,0.4264656,0.4554649,
    0.484115,0.5123732,0.5401928,0.5675238,0.5943126,0.620503,0.6460362,0.6708519,
    0.6948895,0.718089,0.7403926,0.7617461,0.7821008,0.8014145,0.8196537,0.8367943,
    0.8528225,0.8677358,0.8815426,0.8942622,0.9059238,0.9165657,0.9262337,0.9349795,
    0.9428595,0.9499327,0.9562597,0.961901,0.9669162,0.9713627,0.9752957,0.9787666,
    0.981824,0.9845121,0.9868721,0.988941,0.9907526,0.9923371,0.9937217,0.9949305,
    0.9959852,0.9969047,0.997706,0.9984038,1.0
};

// for some serotypes, the fraction who are symptomatic upon primary infection
static const std::vector<double> SYMPTOMATIC_BY_AGE = {
    0.05189621,0.05189621,0.05189621,0.05189621,0.05189621,
    0.1017964,0.1017964,0.1017964,0.1017964,0.1017964,
    0.2774451,0.2774451,0.2774451,0.2774451,0.2774451,
    0.4870259,0.4870259,0.4870259,0.4870259,0.4870259,
    0.4870259,0.4870259,0.4870259,0.4870259,0.4870259,
    0.8522954,0.8522954,0.8522954,0.8522954,0.8522954,
    0.8522954,0.8522954,0.8522954,0.8522954,0.8522954,
    0.9600798,0.9600798,0.9600798,0.9600798,0.9600798,
    0.9600798,0.9600798,0.9600798,0.9600798,0.9600798,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
};

//       RUNS 1979-2014
// DENV1:  7  11  1  10
// DENV2:  1  1  6  14
// DENV3:  2  1  2  1
// DENV4:  1  1  2  1  3
//
//       GAPS 1979-2014
// DENV1:  1  4  2
// DENV2:  7  4  2  1
// DENV3: 17  3  4  5  1
// DENV4:  5  9  1  9  4

// fitted using ABC to Yucatan serotype data, 1979-2014
static const std::vector<double> MEAN_RUN_LENGTH = { 13.03, 8.99, 1.90, 2.17 };
static const std::vector<double> MEAN_GAP_LENGTH = { 3.33, 6.36, 11.10, 9.43 };

// Fraction of days with precipitation in each month, aggregated over 1979-2013
// Derived from NOAA data for airport in Merida
// Jan        Feb        Mar        Apr        May        Jun
// 0.14774282 0.10505319 0.10095012 0.07866667 0.16069057 0.63969171
// Jul        Aug        Sep        Oct        Nov        Dec
// 0.77393075 0.74148297 0.82327586 0.40394089 0.24817518 0.16411683

//2005 Thai mortality data by age from Porapakkham 2010
//double thaimortality[NUM_AGE_CLASSES] = {
//    0.0157,0.0009,0.0009,0.0009,0.0009,0.0005,0.0005,0.0005,0.0005,0.0005,
//    0.0005,0.0005,0.0005,0.0005,0.0005,0.0007,0.0007,0.0007,0.0007,0.0007,
//    0.0009,0.0009,0.0009,0.0009,0.0009,0.0016,0.0016,0.0016,0.0016,0.0016,
//    0.002,0.002,0.002,0.002,0.002,0.0022,0.0022,0.0022,0.0022,0.0022,
//    0.0028,0.0028,0.0028,0.0028,0.0028,0.0038,0.0038,0.0038,0.0038,0.0038,
//    0.005,0.005,0.005,0.005,0.005,0.0077,0.0077,0.0077,0.0077,0.0077,
//    0.012,0.012,0.012,0.012,0.012,0.0185,0.0185,0.0185,0.0185,0.0185,
//    0.0287,0.0287,0.0287,0.0287,0.0287,0.0457,0.0457,0.0457,0.0457,0.0457,
//    0.0767,0.0767,0.0767,0.0767,0.0767,0.1434,0.1434,0.1434,0.1434,0.1434};

namespace dengue {
    namespace standard {
        using std::cout;
        using std::cerr;
        using std::endl;
        using std::string;
        using std::vector;
        using std::set;
        using std::map;
        using std::pair;
        using std::make_pair;
        using std::ifstream;
        using std::istringstream;
        using std::ofstream;
        using std::ostringstream;
        using std::strcmp;
        using std::strtol;
        using std::strtod;
        using std::bitset;
    }
}


struct DynamicParameter {
    DynamicParameter(){};
    DynamicParameter(int s, int d, double v) : start(s), duration(d), value(v) {};
    int start;
    int duration;
    double value;
};


struct VaccinationEvent {
    VaccinationEvent(){};
    VaccinationEvent(int a, int s, double c): age(a), simDay(s), coverage(c) {};
    int age;
    int simDay;
    double coverage;
};


class Parameters {
public:

    Parameters() { define_defaults(); }
    Parameters(int argc, char *argv[]) { define_defaults(); readParameters(argc, argv); }

    void define_defaults();
    void readParameters(int argc, char *argv[]);
    void validate_parameters();
    void loadAnnualIntroductions(std::string annualIntrosFilename);
    void loadAnnualSerotypes() { loadAnnualSerotypes(annualSerotypeFilename); };
    void loadAnnualSerotypes(std::string annualSerotypeFilename);
    void writeAnnualSerotypes(std::string filename) const;
    void loadDailyEIP(std::string dailyEIPFilename, int desired_size = 0);
    void loadDailyMosquitoMultipliers(std::string mosquitoMultiplierFilename, int desired_size = 0);
    void generateAnnualSerotypes(int total_num_years = -1);
    bool simulateAnnualSerotypes;

    static int sampler (const std::vector<double> CDF, const double rand, unsigned int index = 0) {
        while (index < CDF.size() and CDF[index] < rand) index++;
        return index;
    };

    unsigned long int randomseed;
    int nRunLength;
    double betaPM;                                          // scales person-to-mosquito transmission
    double betaMP;                                          // scales mosquito-to-person transmission (includes bite rate)
    double fMosquitoMove;                                   // daily probability of mosquito migration
    std::string mosquitoMoveModel;                          // weighted or uniform mosquito movement to adj. buildings
    double fMosquitoTeleport;                               // daily probability of mosquito teleportation (long-range movement)
    std::vector<double> fVESs;                              // vaccine efficacy for susceptibility (can be leaky or all-or-none)
    std::vector<double> fVESs_NAIVE;                        // VES for initially immunologically naive people
    double fVEI;                                            // vaccine efficacy to reduce infectiousness
    double fVEP;                                            // vaccine efficacy for pathogenicity
    double fVEH;                                            // vaccine efficacy against hospitalization, given disease
    PrimaryPathogenicityModel primaryPathogenicityModel;    // use age-specific values, or constant?
    double annualFlavivirusAttackRate;                      // used for geometric primary pathogenicity model
    std::vector<double> primarySevereFraction;              // fraction of primary cases (symptomatic infections) that are severe
    std::vector<double> secondarySevereFraction;            // fraction of post-primary cases (symptomatic infections) that are severe
    std::vector<double> tertiarySevereFraction;             // fraction of post-primary cases (symptomatic infections) that are severe
    std::vector<double> quaternarySevereFraction;           // fraction of post-primary cases (symptomatic infections) that are severe
    std::vector<double> hospitalizedFraction;               // Probability of being hospitalized, given asymptomatic, mild, and severe infection
    std::vector<double> reportedFraction;                   // Probability of being reported, given asymptomatic, mild, and severe infection
    double infantImmuneProb;                                // Probability that age 0 person is uninfectable, given maternal immunity
    double infantSevereProb;                                // Probability that age 0 person experiences severe disease, given maternal immunity
    bool bVaccineLeaky;                                     // if false, vaccine is all-or-none
    bool bRetroactiveMatureVaccine;                         // if true, infection causes leaky vaccine to jump from naive to mature protection
    std::vector<int> nInitialExposed;                       // serotypes
    std::vector<std::vector<float> > nDailyExposed;         // dimensions are [year][serotype]
    std::vector<int> nInitialInfected;                      // serotypes
    double basePathogenicity;                               // weighted average (over serotypes) post-primary pathogenicity (Pr{infection->symptomatic})
    std::vector<double> serotypePathogenicityRelativeRisks; // Relative risks of symptoms, normalized to have max of 1.0
    double primaryRelativeRisk;                             // how much less pathogenic are primary infections relative to secondary
    double postSecondaryRelativeRisk;                       // risk of symptoms in post-secondary infections relative to secondary infections
    void defineSerotypeRelativeRisks();                     // should be called after reportedFractions (1/expansion factors) are set, if they're going to be

    int nDefaultMosquitoCapacity;
    MosquitoDistribution eMosquitoDistribution;
    std::vector<DynamicParameter> mosquitoMultipliers;
    int getMosquitoMultiplierTotalDuration() const { return mosquitoMultipliers.size() > 0 ?
                                                      mosquitoMultipliers.back().start + mosquitoMultipliers.back().duration
                                                      : 0; }
    std::vector<DynamicParameter> extrinsicIncubationPeriods;
    int getEIPtotalDuration() const { return extrinsicIncubationPeriods.size() > 0 ?
                                       extrinsicIncubationPeriods.back().start + extrinsicIncubationPeriods.back().duration
                                       : 0; }
    bool bSecondaryTransmission;
    std::string populationFilename;
    std::string immunityFilename;
    std::string networkFilename;
    std::string locationFilename;
    std::string peopleOutputFilename;
    std::string yearlyPeopleOutputFilename;
    std::string dailyOutputFilename;
    std::string swapProbFilename;
    std::string annualIntroductionsFilename;                // time series of some external factor determining introduction rate
    std::string annualSerotypeFilename;                     // time series of some external factor determining introduction rate
    std::string dailyEIPfilename;
    std::string mosquitoFilename;
    std::string mosquitoLocationFilename;
    std::vector<double> annualIntroductions;
    double annualIntroductionsCoef;                         // multiplier to rescale external introductions to something sensible
    bool normalizeSerotypeIntros;                           // is expected # of intros held constant, regardless of serotypes # (>0)
    bool simpleEIP;                                         // do all mosquitoes infected on day X have the same EIP? (default=F, e.g. sampled)
    int nDaysImmune;
    bool linearlyWaningVaccine;
    int vaccineImmunityDuration;
    bool vaccineBoosting;                                   // Are we re-vaccinated, either because of waning or because of multi-dose vaccine
    int numVaccineDoses;                                    // Number of times to boost; default is INT_MAX
    int vaccineDoseInterval;                                // How often to we re-vaccinate for initial vaccine course, in days
    int vaccineBoostingInterval;                            // How often to we re-vaccinate for boosting, in days
    std::vector<VaccinationEvent> vaccinationEvents;
    int startDayOfYear;
    bool dailyOutput;
    bool weeklyOutput;
    bool monthlyOutput;
    bool yearlyOutput;
    bool abcVerbose;
    unsigned long int serial;
};

#endif
