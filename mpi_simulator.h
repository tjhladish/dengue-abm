#include <cstdlib>
#include <cstring>
#include <climits>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <assert.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Parameters.h"
#include "Person.h"
#include "Mosquito.h"
#include "Location.h"
#include "Community.h"
#include "Utility.h"

using namespace dengue::standard;
const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);

// Predeclare local functions
Community* build_community(const Parameters* par);
void seed_epidemic(const Parameters* par, Community* community);
vector<int> simulate_epidemic(const Parameters* par, Community* community);

Community* build_community(const Parameters* par) {
    Community* community = new Community(par);
    
    if (!community->loadLocations(par->szLocationFile, par->szNetworkFile)) {
        cerr << "ERROR: Could not load locations" << endl;
        exit(-1);
    }
    if (!community->loadPopulation(par->szPopulationFile, par->szImmunityFile, par->szSwapProbFile)) {
        cerr << "ERROR: Could not load population" << endl;
        exit(-1);
    }

    Person::setPar(par);
    //cerr << community->getNumPerson() << " people" << endl;

    if (!par->bSecondaryTransmission) {
        community->setNoSecondaryTransmission();
    }

    if (par->fPreVaccinateFraction>0.0) {
        community->vaccinate(par->fPreVaccinateFraction);
    }

    if (par->nSizePrevaccinateAge>0) {
        for (int j=0; j<par->nSizePrevaccinateAge; j++) {
            for (int k=par->nPrevaccinateAgeMin[j]; k<=par->nPrevaccinateAgeMax[j]; k++) {
                community->vaccinate(par->fPrevaccinateAgeFraction[j],k);
            }
        }
    }
    
    //for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
    //    par->nNumInitialSusceptible[serotype] = community->getNumSusceptible((Serotype) serotype);
    //}
    return community;
}


void seed_epidemic(const Parameters* par, Community* community) {
    // epidemic may be seeded with initial exposure OR initial infection -- not sure why
    bool attempt_initial_infection = true;
    for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
        if (par->nInitialExposed[serotype] > 0) {
            attempt_initial_infection = false;
            for (int i=0; i<par->nInitialExposed[serotype]; i++)
                community->infect(gsl_rng_uniform_int(RNG, community->getNumPerson()), (Serotype) serotype,0);
        }
    }
    if (attempt_initial_infection) {
        for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
            if(par->nInitialInfected[serotype] > 0) {
                int count = community->getNumInfected(0);

                // must infect nInitialInfected persons -- this bit is mysterious
                while (community->getNumInfected(0) < count + par->nInitialInfected[serotype]) {
                    community->infect(gsl_rng_uniform_int(RNG, community->getNumPerson()), (Serotype) serotype,0);
                }
            }
        }
    }
    return;
}


vector<int> simulate_epidemic(const Parameters* par, Community* community) {
    vector<int> epi_sizes;
    int nNextMosquitoMultiplier = 0;
    int nNextExternalIncubation = 0;
    int epi_ctr = 0;
    int daily_ctr = 0;

    for (int t=0; t<par->nRunLength; t++) {
        if (t%1000==0) cerr << "T: " << t << " daily: " << daily_ctr << " annual: " << epi_ctr << endl;
        if (t%365==0) {
        //if ((t-100)%365==0) {
            if (t >= 365) {
                epi_sizes.push_back(epi_ctr);
            }
            epi_ctr = 0;
        }

        // phased vaccination
        if ((t%365)==0) {
            int year = (int)(t/365);
            for (int i=0; i<par->nSizeVaccinate; i++) {
                if (year==par->nVaccinateYear[i]) {
                    community->vaccinate(par->fVaccinateFraction[i],par->nVaccinateAge[i]);
                    //cerr << "vaccinating " << par->fVaccinateFraction[i]*100 << "% of age " << par->nVaccinateAge[i] << endl;
                }
            }
        }

        // seed epidemic
        int intro_count = 0;
        {
            int numperson = community->getNumPerson();
            for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
                if (par->nDailyExposed[serotype] <= 0) continue;
                const int num_exposed = gsl_ran_poisson(RNG, par->nDailyExposed[serotype]);
                for (int i=0; i<num_exposed; i++) {
                    // gsl_rng_uniform_int returns on [0, numperson-1]
                    int transmit_to_id = gsl_rng_uniform_int(RNG, numperson) + 1; 
                    if (community->infect(transmit_to_id, (Serotype) serotype, t)) {
                        intro_count++;
                    }
                }
            }
        }

        // mosquito population seasonality?
        if (par->nSizeMosquitoMultipliers>0 && (t%365)==par->nMosquitoMultiplierCumulativeDays[nNextMosquitoMultiplier]) {
            community->setMosquitoMultiplier(par->fMosquitoMultipliers[nNextMosquitoMultiplier]);
            nNextMosquitoMultiplier = (nNextMosquitoMultiplier+1)%par->nSizeMosquitoMultipliers;
        }
        if (par->nSizeExternalIncubation>0 && (t%365)==par->nExternalIncubationCumulativeDays[nNextExternalIncubation]) {
            community->setExternalIncubation(par->nExternalIncubation[nNextExternalIncubation]);
            nNextExternalIncubation = (nNextExternalIncubation+1)%par->nSizeExternalIncubation;
        }

        int daily_infection_ctr = 0;
        community->tick(t);

        for (int i=community->getNumPerson()-1; i>=0; i--) {
            Person *p = community->getPerson(i);
            if (p->isInfected(t) and p->isNewlyInfected(t)) {
                daily_infection_ctr++;
            }
        }
        epi_ctr += daily_infection_ctr;
        daily_ctr = daily_infection_ctr;
    }
    return epi_sizes;
}
