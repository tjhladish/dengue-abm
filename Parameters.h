#ifndef __PARAMETERS_H
#define __PARAMETERS_H

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
const static int MAXPERSONAGE = 95;                           // maximum age-1 for a person
static const int MAXINCUBATION = 9;                           // max incubation period in days
static const int MAXHISTORY = 50;                             // length of exposure history in years

// from Community.h
const static int STEPSPERDAY = 3;                             // number of time steps per day
//const static int MAXINCUBATION = 15;                          // maximum incubation period for humans
const static int MOSQUITOINCUBATION = 11;                     // number of days for mosquito incubation (extrinsic incubation period)
const static int MAXRUNTIME = 7400;                           // maximum number of simulation days (+ extra for mosquito lifetime)

// from Mosquito.h
static const int MAXMOSQUITOAGE = 60;                                 // maximum age of mosquito in days

#endif
