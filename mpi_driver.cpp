// driver.cpp
// driver for dengue epidemic model

#include <cstdlib>
#include <cstring>
#include <climits>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <assert.h>
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
void write_output(const Parameters* par, Community* community, vector<int> initial_susceptibles);

int main(int argc, char *argv[]) {
    const Parameters* par = new Parameters(argc, argv);

    gsl_rng_set(RNG, par->randomseed);

    Community* community = build_community(par);
    vector<int> initial_susceptibles = community->getNumSusceptible();
    seed_epidemic(par, community);
    vector<int> epi_sizes = simulate_epidemic(par, community);
    cout << mean(epi_sizes) << " " << stdev(epi_sizes) << " " << max_element(epi_sizes) << endl;
    write_output(par, community, initial_susceptibles);
    
    return 0;
}
 
Community* build_community(const Parameters* par) {
    Community* community = new Community(par);

    if (!community->loadLocations(par->szLocationFile, par->szNetworkFile)) {
        cerr << "Could not load locations" << endl;
        exit(-1);
    }
    if (!community->loadPopulation(par->szPopulationFile, par->szImmunityFile, par->szSwapProbFile)) {
        cerr << "Could not load population" << endl;
        exit(-1);
    }

    Person::setPar(par);
    cerr << community->getNumPerson() << " people" << endl;

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
    for (int t=0; t<par->nRunLength; t++) {
        //if (t%10==0) cerr << "Time " << t << " epi_size " << epi_ctr << endl;
        if ((t-100)%365==0) {
            if (t > 365) {
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
                // home location
                //cout << t << ",p," << p->getID() << "," << p->getLocation(0)->getID() << "," << 1 + (int) p->getSerotype() << "," << (p->isSymptomatic(t)?1:0) << "," << (p->isWithdrawn(t)?1:0) << "," << (p->isNewlyInfected(t)?1:0) << endl;
            }
        }
        epi_ctr += daily_infection_ctr;
    }
    return epi_sizes;
}



void write_output(const Parameters* par, Community* community, vector<int> numInitialSusceptible) {
    // output daily infected/symptomatic file
    if (par->szDailyFile.length()>0) {
        cerr << "outputing daily infected/symptomatic information to " << par->szDailyFile << endl;
        ofstream dailyFile;
        dailyFile.open(par->szDailyFile.c_str());
        if(dailyFile.fail()) {
            cerr << "ERROR: Daily file '" << par->szDailyFile << "' cannot be open for writing." << endl;
            exit(-1);
        }
        dailyFile << "day,newly infected DENV1,newly infected DENV2,newly infected DENV3,newly infected DENV4,"
                  << "newly symptomatic DENV1,newly symptomatic DENV2,newly symptomatic DENV3,newly symptomatic DENV4" << endl;
        vector< vector<int> > infected =    community->getNumNewlyInfected();
        vector< vector<int> > symptomatic = community->getNumNewlySymptomatic();
        for (int t=0; t<par->nRunLength; t++) {
            dailyFile << t << ",";
            for (int i=0; i<NUM_OF_SEROTYPES; i++)   dailyFile << infected[i][t] << ",";
            for (int i=0; i<NUM_OF_SEROTYPES-1; i++) dailyFile << symptomatic[i][t] << ","; 
            dailyFile << symptomatic[NUM_OF_SEROTYPES-1][t] << endl;
        }
        dailyFile.close();
    }
}

