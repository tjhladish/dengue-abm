#ifndef __PARAMETERS_H
#define __PARAMETERS_H

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <assert.h>
#include <bitset>
#include <gsl/gsl_rng.h>

static const int VERSION_NUMBER_MAJOR = 1;
static const int VERSION_NUMBER_MINOR = 1;

enum Serotype {
    SEROTYPE_1,
    SEROTYPE_2,
    SEROTYPE_3,
    SEROTYPE_4,
    NUM_OF_SEROTYPES,  // Make sure this is second to last
    NULL_SEROTYPE      // Make sure this is last
};

enum MosquitoDistribution {
    UNIFORM,
    EXPONENTIAL,
    NUM_OF_DISTRIBUTIONS
};

extern const gsl_rng* RNG;// = gsl_rng_alloc (gsl_rng_taus2);

static const int INFECTIOUS_PERIOD_PRI = 5;                   // number of days until recovery from primary infection
static const int INFECTIOUS_PERIOD_SEC = 4;                   // number of days until recovery from secondary infection

// from Person.h
static const int MAX_PERSON_AGE = 95;                         // maximum age-1 for a person
static const int MAX_INCUBATION = 9;                          // max incubation period in days
static const int MAX_HISTORY = 50;                            // length of exposure history in years

// cdf of incubation period, starting from day 1 (from Nishiura 2007)
const double INCUBATION_DISTRIBUTION[MAX_INCUBATION] = {
    0,0,0.03590193,0.5070053,0.8248687,0.9124343,0.949212,0.974606,1
};

// from Community
static const int STEPS_PER_DAY = 3;                           // number of time steps per day
static const int MOSQUITO_INCUBATION = 11;                    // number of days for mosquito incubation (extrinsic incubation period)
static const int MAX_RUN_TIME = 7400;                         // maximum number of simulation days (+ extra for mosquito lifetime)
static const float DAILY_BITING_PDF[STEPS_PER_DAY] = {0.08, 0.76, 0.16};  // probability of biting at 3 different times of day (as defined in Location.h)
 
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
static const double MOSQUITO_AGE_DISTRIBUTION[MAX_MOSQUITO_AGE] = {
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
static const double SYMPTOMATIC_BY_AGE[MAX_PERSON_AGE] = {
    0.05189621,0.05189621,0.05189621,0.05189621,0.05189621,
    0.1017964,0.1017964,0.1017964,0.1017964,0.1017964,
    0.2774451,0.2774451,0.2774451,0.2774451,0.2774451,
    0.4870259,0.4870259,0.4870259,0.4870259,0.4870259,
    0.4870259,0.4870259,0.4870259,0.4870259,0.4870259,
    0.8522954,0.8522954,0.8522954,0.8522954,0.8522954,
    0.8522954,0.8522954,0.8522954,0.8522954,0.8522954,
    0.9600798,0.9600798,0.9600798,0.9600798,0.9600798,
    0.9600798,0.9600798,0.9600798,0.9600798,0.9600798,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
};

//2005 Thai mortality data by age from Porapakkham 2010
//double thaimortality[MAX_PERSON_AGE] = {
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


class Parameters {
public:

    Parameters(int argc, char *argv[]) {
        readParameters(argc,argv);
    }

    void readParameters(int argc, char *argv[]);

    int randomseed;
    int nRunLength;
    double betaPM;                                          // scales person-to-mosquito transmission
    double betaMP;                                          // scales mosquito-to-person transmission (includes bite rate)
    double fMosquitoMove;                                   // daily probability of mosquito migration
    std::string szMosquitoMoveModel;                          // weighted or uniform mosquito movement to adj. buildings
    double fMosquitoTeleport;                               // daily probability of mosquito teleportation (long-range movement)
    double fVEI;
    double fVEP;
    std::vector<double> fVESs;
    double fPreVaccinateFraction;
    bool bVaccineLeaky; // if false, vaccine is all-or-none
    int nInitialExposed[NUM_OF_SEROTYPES];                  // serotypes
    int nDailyExposed[NUM_OF_SEROTYPES];                    // serotypes
    int nInitialInfected[NUM_OF_SEROTYPES];                 // serotypes
    std::vector<double> fPrimaryPathogenicity;              // serotypes
    std::vector<double> fSecondaryScaling;                  //
    int nDefaultMosquitoCapacity;
    MosquitoDistribution eMosquitoDistribution;
    double fMosquitoMultipliers[54];                        // number of mosquitoes per week (up to 53), conforming to a 365-day cycle
    double nMosquitoMultiplierDays[54];                     // when to change the mosquito population
    double nMosquitoMultiplierCumulativeDays[54];           // when to change the mosquito population
    int nSizeMosquitoMultipliers;
    bool bSecondaryTransmission;
    std::string szPopulationFile;
    std::string szImmunityFile;
    std::string szNetworkFile;
    std::string szLocationFile;
    int nNumInitialSusceptible[NUM_OF_SEROTYPES];
    std::string szPeopleFile;
    std::string szYearlyPeopleFile;
    std::string szDailyFile;
    std::string szSwapProbFile;
    int nDaysImmune;
    int nSizeVaccinate;
    int nVaccinateYear[500];                                // when to vaccinate
    int nVaccinateAge[500];                                 // who to vaccinate
    double fVaccinateFraction[500];                         // fraction of age group to vaccinate
    int nSizePrevaccinateAge;
    int nPrevaccinateAgeMin[100];
    int nPrevaccinateAgeMax[100];
    double fPrevaccinateAgeFraction[100];
    int nMaxInfectionParity;
};

#endif
