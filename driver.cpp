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


using namespace dengue::standard;
const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);

// Predeclare local functions
Community* build_community(Parameters* par);
void seed_epidemic(Parameters* par, Community* community);
void simulate_epidemic(Parameters* par, Community* community);
void write_output(Parameters* par, Community* community);

int main(int argc, char *argv[]) {
    Parameters* par = new Parameters(argc, argv);

    gsl_rng_set(RNG, par->randomseed);

    Community* community = build_community(par);
    seed_epidemic(par, community);
    simulate_epidemic(par, community);
    write_output(par, community);

    return 0;
}
 
Community* build_community(Parameters* par) {
    Community* community = new Community();
    if (!community->loadLocations(par->szLocationFile.c_str(),par->szNetworkFile.c_str())) {
        cerr << "Could not load locations" << endl;
        exit(-1);
    }
    if (!community->loadPopulation(par->szPopulationFile.c_str(),par->szImmunityFile.c_str())) {
        cerr << "Could not load population" << endl;
        exit(-1);
    }
    community->setBetaPM(par->betaPM);
    community->setBetaMP(par->betaMP);
    if (par->fVES<0.0) {
        community->setVESs(par->fVESs);
    } else {
        community->setVES(par->fVES);
    }
    community->setVEI(par->fVEI);
    community->setVEP(par->fVEP);
    community->setMaxInfectionParity(par->nMaxInfectionParity);
    Person::setDaysImmune(par->nDaysImmune);
    community->setMosquitoMoveProbability(par->fMosquitoMove);
    community->setMosquitoTeleportProbability(par->fMosquitoTeleport);
    community->setPrimaryPathogenicityScaling(par->fPrimaryPathogenicity);
    community->setSecondaryPathogenicityScaling(par->fSecondaryScaling);

    cerr << community->getNumPerson() << " people" << endl;

    if (!par->bSecondaryTransmission) {
        community->setNoSecondaryTransmission();
        //    cerr << "No secondary transmission, 1 person infected" << endl;
    }

    if (par->fPreVaccinateFraction>0.0)
        community->vaccinate(par->fPreVaccinateFraction);

    if (par->nSizePrevaccinateAge>0) {
        for (int j=0; j<par->nSizePrevaccinateAge; j++)
            for (int k=par->nPrevaccinateAgeMin[j]; k<=par->nPrevaccinateAgeMax[j]; k++)
                community->vaccinate(par->fPrevaccinateAgeFraction[j],k);
    }

    for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
        par->nNumInitialSusceptible[serotype] = community->getNumSusceptible((Serotype) serotype);
    }
    return community;
}


void seed_epidemic(Parameters* par, Community* community) {
    // epidemic may be seeded with initial exposure OR initial infection -- not sure why
    bool attempt_initial_infection = true;
    for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
        if (par->nInitialExposed[serotype] > 0) {
            attempt_initial_infection = false;
            cerr << "Initial serotype " << serotype+1 << " exposed = " << par->nInitialExposed[serotype] << endl;
            for (int i=0; i<par->nInitialExposed[serotype]; i++)
                community->infect(gsl_rng_uniform_int(RNG, community->getNumPerson()), (Serotype) serotype,0);
        }
    }
    if (attempt_initial_infection) {
        for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
            if(par->nInitialInfected[serotype] > 0) {
                cerr << "Initial serotype " << serotype+1 << " infected = " << par->nInitialInfected[serotype] << endl;
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


void write_people_file(Parameters* par, Community* community, int time) {
    if (time%365==364 && par->szYearlyPeopleFile.length()>0) {
        ofstream peopleFile;
        ostringstream ssFilename;
        ssFilename << par->szYearlyPeopleFile << ((int)(time/365)) << ".csv";
        cerr << "outputing yearly people information to " << ssFilename.str() << endl;
        peopleFile.open(ssFilename.str().c_str());
        if(peopleFile.fail()) {
            cerr << "ERROR: People file '" << par->szPeopleFile << "' cannot be open for writing." << endl;
            exit(-1);
        }
        peopleFile << "pid,serotype,infectiontime,symptomtime,withdrawtime,recoverytime,immdenv1,immdenv2,immdenv3,immdenv4" << endl;
        for (int i=0; i<community->getNumPerson(); i++) {
            Person *p = community->getPerson(i);
            for (int j=p->getNumInfections()-1; j>=0; j--) {
                peopleFile << p->getID() << "," 
                    << 1 + (int) p->getSerotype(j) << "," 
                    << p->getInfectedTime(j) << "," 
                    << p->getSymptomTime(j) << "," 
                    << p->getWithdrawnTime(j) << "," 
                    << p->getRecoveryTime(j) << "," 
                    << (p->isSusceptible(SEROTYPE_1)?0:1) << "," 
                    << (p->isSusceptible(SEROTYPE_2)?0:1) << "," 
                    << (p->isSusceptible(SEROTYPE_3)?0:1) << "," 
                    << (p->isSusceptible(SEROTYPE_4)?0:1) << endl;
            }
        }
        peopleFile.close();
    }
    return;
}

void simulate_epidemic(Parameters* par, Community* community) {
    int nNextMosquitoMultiplier = 0;
    if (par->bSecondaryTransmission) cout << "time,type,id,location,serotype,symptomatic,withdrawn" << endl;
    for (int t=0; t<par->nRunLength; t++) {
        // phased vaccination
        if ((t%365)==0) {
            int year = (int)(t/365);
            for (int i=0; i<par->nSizeVaccinate; i++) {
                if (year==par->nVaccinateYear[i]) {
                    community->vaccinate(par->fVaccinateFraction[i],par->nVaccinateAge[i]);
                    cerr << "vaccinating " << par->fVaccinateFraction[i]*100 << "% of age " << par->nVaccinateAge[i] << endl;
                }
            }
        }

        // seed epidemic
        if (par->nDailyExposed[0]>0 || par->nDailyExposed[1]>0 || par->nDailyExposed[2]>0 || par->nDailyExposed[3]>0) {
            int count = 0;
            int numperson = community->getNumPerson();
            for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++)
                for (int i=0; i<par->nDailyExposed[serotype]; i++)
                    if (community->infect(gsl_rng_uniform_int(RNG, numperson),(Serotype) serotype,t))
                        count++;
        }

        // mosquito population seasonality?
        if (par->nSizeMosquitoMultipliers>0 && (t%365)==par->nMosquitoMultiplierCumulativeDays[nNextMosquitoMultiplier]) {
            community->setMosquitoMultiplier(par->fMosquitoMultipliers[nNextMosquitoMultiplier]);
            nNextMosquitoMultiplier = (nNextMosquitoMultiplier+1)%par->nSizeMosquitoMultipliers;
        }

        if (par->bSecondaryTransmission) {
            // print out infectious mosquitoes
            for (int i=community->getNumInfectiousMosquitoes()-1; i>=0; i--) {
                Mosquito *p = community->getInfectiousMosquito(i);
                cout << t << ",mi," << p->getID() << "," << p->getLocation()->getID() << "," << "," << endl;
            }
            // print out exposed mosquitoes
            for (int i=community->getNumExposedMosquitoes()-1; i>=0; i--) {
                Mosquito *p = community->getExposedMosquito(i);
                // "current" location
                cout << t << ",me," << p->getID() << "," << p->getLocation()->getID() << "," << 1 + (int) p->getSerotype() << "," << "," << endl;
            }
            // print out infected people
            for (int i=community->getNumPerson()-1; i>=0; i--) {
                Person *p = community->getPerson(i);
                if (p->isInfected(t))
                    // home location
                    cout << t 
                         << ",p,"
                         << p->getID() << "," 
                         << p->getLocation(0)->getID() << "," 
                         << 1 + (int) p->getSerotype() << "," 
                         << (p->isSymptomatic(t)?1:0) << "," 
                         << (p->isWithdrawn(t)?1:0) << endl;
                // printing out the home location of each infected person is not useful, but we don't keep track of where they get infected
            }
        }
        write_people_file(par, community, t);

        if (t%100==0) cerr << "Time " << t << endl;
        community->tick();
    }
    return;
}



void write_output(Parameters* par, Community* community) {
    if (!par->bSecondaryTransmission) {
        // outputs
        //   number of secondary infections by serotype (4)
        //   number of households infected
        //   age of index case
        //   age(s) of secondary cases
        for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
            cout << (par->nNumInitialSusceptible[serotype] - community->getNumSusceptible((Serotype) serotype)) << " ";
        }
        //    cout << "secondary infections" << endl;

        int ages[100];
        int times[100];
        int numages=0;
        int indexage=-1;
        int homeids[100];
        int numhomes=0;
        // ages of infected people
        for (int i=community->getNumPerson()-1; i>=0; i--) {
            Person *p = community->getPerson(i);
            int t = p->getInfectedTime();
            if (t>=0) {
                if (t==0)
                    indexage = p->getAge();
                else {
                    ages[numages] = p->getAge();
                    times[numages] = t;
                    numages++;
                }
                int homeid = p->getHomeID();
                bool bFound = false;
                for (int j=0; j<numhomes; j++)
                    if (homeids[j]==homeid)
                        bFound = true;
                if (!bFound)
                    homeids[numhomes++] = homeid;
            }
        }
        cout << indexage << " " << numhomes << " " << numages;
        for (int i=0; i<numages; i++)
            cout << " " << ages[i];
        for (int i=0; i<numages; i++)
            cout << " " << times[i];
        cout << endl;
    }

    // output daily infected/symptomatic file
    if (par->szDailyFile.length()>0) {
        cerr << "outputing daily infected/symptomatic information to " << par->szDailyFile << endl;
        ofstream dailyFile;
        dailyFile.open(par->szDailyFile.c_str());
        if(dailyFile.fail()) {
            cerr << "ERROR: Daily file '" << par->szDailyFile << "' cannot be open for writing." << endl;
            exit(-1);
        }
        dailyFile << "day,newly infected DENV1,newly infected DENV2,newly infected DENV3,newly infected DENV4,newly symptomatic DENV1,newly symptomatic DENV2,newly symptomatic DENV3,newly symptomatic DENV4" << endl;
        const int *pInfected1 = community->getNumNewlyInfected(SEROTYPE_1);
        const int *pInfected2 = community->getNumNewlyInfected(SEROTYPE_2);
        const int *pInfected3 = community->getNumNewlyInfected(SEROTYPE_3);
        const int *pInfected4 = community->getNumNewlyInfected(SEROTYPE_4);
        const int *pSymptomatic1 = community->getNumNewlySymptomatic(SEROTYPE_1);
        const int *pSymptomatic2 = community->getNumNewlySymptomatic(SEROTYPE_2);
        const int *pSymptomatic3 = community->getNumNewlySymptomatic(SEROTYPE_3);
        const int *pSymptomatic4 = community->getNumNewlySymptomatic(SEROTYPE_4);
        for (int t=0; t<par->nRunLength; t++) {
            dailyFile << t << ","
                      << pInfected1[t] << "," << pInfected2[t] << "," << pInfected3[t] << "," << pInfected4[t] << "," 
                      << pSymptomatic1[t] << "," << pSymptomatic2[t] << "," << pSymptomatic3[t] << "," << pSymptomatic4[t] << endl;
        }
        dailyFile.close();
    }

    // output people file
    if (par->szPeopleFile.length()>0) {
        cerr << "outputing people information to " << par->szPeopleFile << endl;
        ofstream peopleFile;
        peopleFile.open(par->szPeopleFile.c_str());
        if(peopleFile.fail()) {
            cerr << "ERROR: People file '" << par->szPeopleFile << "' cannot be open for writing." << endl;
            exit(-1);
        }
        peopleFile << "pid,serotype,infectiontime,symptomtime,withdrawtime,recoverytime,immdenv1,immdenv2,immdenv3,immdenv4,vaccinated" << endl;
        for (int i=0; i<community->getNumPerson(); i++) {
            Person *p = community->getPerson(i);
            for (int j=p->getNumInfections()-1; j>=0; j--) {
                peopleFile << p->getID() << "," 
                    << 1 + (int) p->getSerotype(j) << "," 
                    << p->getInfectedTime(j) << "," 
                    << p->getSymptomTime(j) << "," 
                    << p->getWithdrawnTime(j) << "," 
                    << p->getRecoveryTime(j) << "," 
                    << (p->isSusceptible(SEROTYPE_1)?0:1) << "," 
                    << (p->isSusceptible(SEROTYPE_2)?0:1) << "," 
                    << (p->isSusceptible(SEROTYPE_3)?0:1) << "," 
                    << (p->isSusceptible(SEROTYPE_4)?0:1) << "," 
                    << (p->isVaccinated()?1:0) << endl;
            }
        }
        peopleFile.close();
    }
}

