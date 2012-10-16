#ifndef __PARAMETERS_H
#define __PARAMETERS_H

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <assert.h>
#include <gsl/gsl_rng.h>

static const int VERSIONNUMBERMAJOR = 1;
static const int VERSIONNUMBERMINOR = 0;

enum Serotype {
    SEROTYPE_1,
    SEROTYPE_2,
    SEROTYPE_3,
    SEROTYPE_4,
    NUM_OF_SEROTYPES,  // Make sure this is second to last
    NULL_SEROTYPE      // Make sure this is last
};

extern const gsl_rng* RNG;// = gsl_rng_alloc (gsl_rng_taus2);

// from Person.h
static const int MAXPERSONAGE = 95;                           // maximum age-1 for a person
static const int MAXINCUBATION = 9;                           // max incubation period in days
static const int MAXHISTORY = 50;                             // length of exposure history in years

// from Community.h
static const int STEPSPERDAY = 3;                             // number of time steps per day
static const int MOSQUITOINCUBATION = 11;                     // number of days for mosquito incubation (extrinsic incubation period)
static const int MAXRUNTIME = 7400;                           // maximum number of simulation days (+ extra for mosquito lifetime)
//static const int MAXINCUBATION = 15;                          // maximum incubation period for humans

// from Mosquito.h
static const int MAXMOSQUITOAGE = 60;                                 // maximum age of mosquito in days


namespace dengue {
    namespace standard {
        using std::cout;
        using std::cerr;
        using std::endl;
        using std::string;
        using std::vector;
        using std::ifstream;
        using std::istringstream;
        using std::ofstream;
        using std::ostringstream;
        using std::strcmp;
        using std::strtol;
        using std::strtod;
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
    double betaPM;
    double betaMP;
    double fMosquitoMove;
    double fMosquitoTeleport;
    double fVES;
    double fVEI;
    double fVEP;
    double fVESs[4];
    double fPreVaccinateFraction;
    int nInitialExposed[NUM_OF_SEROTYPES];                  // serotypes
    int nDailyExposed[NUM_OF_SEROTYPES];                    // serotypes
    int nInitialInfected[NUM_OF_SEROTYPES];                 // serotypes
    double fPrimaryPathogenicity[NUM_OF_SEROTYPES];         // serotypes
    double fSecondaryScaling[NUM_OF_SEROTYPES];             //
    int nDefaultMosquitoCapacity;
    double fMosquitoMultipliers[54];                        // number of mosquitoes per week (up to 53), conforming to a 365-day cycle
    double nMosquitoMultiplierDays[54];                     // when to change the mosquito population
    double nMosquitoMultiplierCumulativeDays[54];           // when to change the mosquito population
    int nSizeMosquitoMultipliers;
    bool bSecondaryTransmission;
    std::string szPopulationFile;
    std::string szImmunityFile;
    std::string szNetworkFile;
    std::string szLocationFile;
    int nNumInitialSusceptible[4];
    std::string szPeopleFile;
    std::string szYearlyPeopleFile;
    std::string szDailyFile;
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
