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

using namespace std;

const int VERSIONNUMBERMAJOR = 1;
const int VERSIONNUMBERMINOR = 0;

int main(int argc, char *argv[]) {
    cerr << "main rng: " << RNG << endl;


    int randomseed = 5489;
    int nRunLength = 100;
    double betaPM = 0.2;
    double betaMP = 0.1;
    double fMosquitoMove = 0.2;
    double fMosquitoTeleport = 0.01;
    double fVES = 0.95;
    double fVEI = 0.0;
    double fVEP = 0.0;
    double fVESs[4] = {0.0,0.0,0.0,0.0};
    double fPreVaccinateFraction = 0.0;
    //const int NUM_OF_SEROTYPES=4;
    int nInitialExposed[NUM_OF_SEROTYPES];                                // serotypes
    int nDailyExposed[NUM_OF_SEROTYPES];                                  // serotypes
    int nInitialInfected[NUM_OF_SEROTYPES];                               // serotypes
    double fPrimaryPathogenicity[NUM_OF_SEROTYPES];                       // serotypes
    double fSecondaryScaling[NUM_OF_SEROTYPES];                           //
    int nDefaultMosquitoCapacity = 20;                                // mosquitoes per location
    double fMosquitoMultipliers[54];                                  // number of mosquitoes per week (up to 53), conforming to a 365-day cycle
    double nMosquitoMultiplierDays[54];                               // when to change the mosquito population
    double nMosquitoMultiplierCumulativeDays[54];                     // when to change the mosquito population
    int nSizeMosquitoMultipliers = 0;
    bool bSecondaryTransmission = true;
    string szPopulationFile("population-64.txt");
    string szImmunityFile("");
    string szNetworkFile("locations-network-64.txt");
    string szLocationFile("locations-64.txt");
    int nNumInitialSusceptible[4];
    string szPeopleFile="";
    string szYearlyPeopleFile="";
    string szDailyFile="";
    int nDaysImmune = 365;
    int nSizeVaccinate = 0;                                           // number of parts in phased vaccination
    int nVaccinateYear[500];                                          // when to vaccinate
    int nVaccinateAge[500];                                           // who to vaccinate
    double fVaccinateFraction[500];                                   // fraction of age group to vaccinate
    int nSizePrevaccinateAge = 0;
    int nPrevaccinateAgeMin[100];
    int nPrevaccinateAgeMax[100];
    double fPrevaccinateAgeFraction[100];
    int nMaxInfectionParity = NUM_OF_SEROTYPES;//Community::MAXSEROTYPES;

    for (int i=0; i<NUM_OF_SEROTYPES; i++) {
        nInitialExposed[i]=0;
        nInitialInfected[i]=0;
        nDailyExposed[i]=0;
    }
    fPrimaryPathogenicity[0] = fPrimaryPathogenicity[2] = 1.0;
    fPrimaryPathogenicity[1] = fPrimaryPathogenicity[3] = 0.25;
    fSecondaryScaling[0] = fSecondaryScaling[1] = fSecondaryScaling[2] = fSecondaryScaling[3] = 1.0;

    //cerr << "Dengue model, Version " << VERSIONNUMBERMAJOR << "." << VERSIONNUMBERMINOR << endl;
    //cerr << "written by Dennis Chao in 2012" << endl;

    if (argc>1) {
        for (int i=1; i<argc; i++) {
            char **end = NULL;
            if (strcmp(argv[i], "-randomseed")==0) {
                randomseed = strtol(argv[i+1],end,10);
                i++;
            }
            else if (strcmp(argv[i], "-runlength")==0) {
                nRunLength = strtol(argv[i+1],end,10);
                i++;
            }
            else if (strcmp(argv[i], "-initialinfected")==0) {
                for (int j=0; j<NUM_OF_SEROTYPES; j++)
                    nInitialInfected[j]=strtol(argv[i+1+j],end,10);
                i+=NUM_OF_SEROTYPES;
            }
            else if (strcmp(argv[i], "-initialexposed")==0) {
                for (int j=0; j<NUM_OF_SEROTYPES; j++)
                    nInitialExposed[j]=strtol(argv[i+1+j],end,10);
                i+=NUM_OF_SEROTYPES;
            }
            else if (strcmp(argv[i], "-dailyexposed")==0) {
                for (int j=0; j<NUM_OF_SEROTYPES; j++)
                    nDailyExposed[j]=strtol(argv[i+1+j],end,10);
                i+=NUM_OF_SEROTYPES;
            }
            else if (strcmp(argv[i], "-primarypathogenicity")==0) {
                for (int j=0; j<NUM_OF_SEROTYPES; j++)
                    fPrimaryPathogenicity[j]=strtod(argv[i+1+j],end);
                i+=NUM_OF_SEROTYPES;
            }
            else if (strcmp(argv[i], "-secondaryscaling")==0) {
                for (int j=0; j<NUM_OF_SEROTYPES; j++)
                    fSecondaryScaling[j]=strtod(argv[i+1+j],end);
                i+=NUM_OF_SEROTYPES;
            }
            else if (strcmp(argv[i], "-betapm")==0) {
                betaPM = strtod(argv[i+1],end);
                i++;
            }
            else if (strcmp(argv[i], "-betamp")==0) {
                betaMP = strtod(argv[i+1],end);
                i++;
            }
            else if (strcmp(argv[i], "-mosquitomove")==0) {
                fMosquitoMove = strtod(argv[i+1],end);
                i++;
            }
            else if (strcmp(argv[i], "-mosquitoteleport")==0) {
                fMosquitoTeleport = strtod(argv[i+1],end);
                i++;
            }
            else if (strcmp(argv[i], "-mosquitocapacity")==0) {
                nDefaultMosquitoCapacity = strtol(argv[i+1],end,10);
                i++;
            }
            else if (strcmp(argv[i], "-mosquitomultipliers")==0) {
                nSizeMosquitoMultipliers = strtol(argv[i+1],end,10);
                i++;
                nMosquitoMultiplierCumulativeDays[0] = 0;
                for (int j=0; j<nSizeMosquitoMultipliers; j++) {
                    nMosquitoMultiplierDays[j] = strtol(argv[i+1],end,10);
                    nMosquitoMultiplierCumulativeDays[j+1] =
                        nMosquitoMultiplierCumulativeDays[j]+nMosquitoMultiplierDays[j];
                    i++;
                    fMosquitoMultipliers[j] = strtod(argv[i+1],end);
                    i++;
                }
            }
            else if (strcmp(argv[i], "-vaccinatephased")==0) {
                nSizeVaccinate = strtol(argv[i+1],end,10);
                i++;
                for (int j=0; j<nSizeVaccinate; j++) {
                    nVaccinateYear[j] = strtol(argv[i+1],end,10);
                    i++;
                    nVaccinateAge[j] = strtol(argv[i+1],end,10);
                    i++;
                    fVaccinateFraction[j] = strtod(argv[i+1],end);
                    i++;
                }
            }
            else if (strcmp(argv[i], "-daysimmune")==0) {
                nDaysImmune = strtol(argv[i+1],end,10);
                i++;
            }
            else if (strcmp(argv[i], "-maxinfectionparity")==0) {
                nMaxInfectionParity = strtol(argv[i+1],end,10);
                i++;
                assert(nMaxInfectionParity>0 && nMaxInfectionParity<=NUM_OF_SEROTYPES);
            }
            else if (strcmp(argv[i], "-VES")==0 || strcmp(argv[i], "-ves")==0) {
                fVES = strtod(argv[i+1],end);
                i++;
            }
            else if (strcmp(argv[i], "-VESs")==0 || strcmp(argv[i], "-vess")==0) {
                // different VES for each serotype
                fVESs[0] = strtod(argv[i+1],end);
                fVESs[1] = strtod(argv[i+2],end);
                fVESs[2] = strtod(argv[i+3],end);
                fVESs[3] = strtod(argv[i+4],end);
                fVES=-1.0;                                            // make fVES parameter invalid
                i+=4;
            }
            else if (strcmp(argv[i], "-VEI")==0 || strcmp(argv[i], "-vei")==0) {
                fVEI = strtod(argv[i+1],end);
                i++;
            }
            else if (strcmp(argv[i], "-VEP")==0 || strcmp(argv[i], "-vep")==0) {
                fVEP = strtod(argv[i+1],end);
                i++;
            }
            else if (strcmp(argv[i], "-prevaccinate")==0) {
                fPreVaccinateFraction = strtod(argv[i+1],end);
                i++;
            }
            else if (strcmp(argv[i], "-prevaccinateage")==0) {
                nSizePrevaccinateAge = strtol(argv[i+1],end,10);
                i++;
                for (int j=0; j<nSizePrevaccinateAge; j++) {
                    nPrevaccinateAgeMin[j] = strtol(argv[i+1],end,10);
                    i++;
                    nPrevaccinateAgeMax[j] = strtol(argv[i+1],end,10);
                    i++;
                    fPrevaccinateAgeFraction[j] = strtod(argv[i+1],end);
                    i++;
                }
            }
            else if (strcmp(argv[i], "-nosecondary")==0) {
                bSecondaryTransmission = false;
            }
            else if (strcmp(argv[i], "-popfile")==0) {
                szPopulationFile = argv[i+1];
                i++;
            }
            else if (strcmp(argv[i], "-immfile")==0) {
                szImmunityFile = argv[i+1];
                i++;
            }
            else if (strcmp(argv[i], "-locfile")==0) {
                szLocationFile = argv[i+1];
                i++;
            }
            else if (strcmp(argv[i], "-netfile")==0) {
                szNetworkFile = argv[i+1];
                i++;
            }
            else if (strcmp(argv[i], "-peoplefile")==0) {
                szPeopleFile = argv[i+1];
                i++;
            }
            else if (strcmp(argv[i], "-yearlypeoplefile")==0) {
                szYearlyPeopleFile = argv[i+1];
                i++;
            }
            else if (strcmp(argv[i], "-dailyfile")==0) {
                szDailyFile = argv[i+1];
                i++;
            }
            else {
                //cerr << "Unknown option: " << argv[i] << endl;
                exit(-1);
            }
        }
    }
    gsl_rng_set(RNG, randomseed);
    //gsl_rng* rng_clone = gsl_rng_clone(RNG);
    ////cerr << "RANDOM NUM ---------------------------------------> " << gsl_rng_uniform(RNG) << endl;
    ////cerr << "RANDOM NUM clone----------------------------------> " << gsl_rng_uniform(rng_clone) << endl;
    ////cerr << "RANDOM NUM2 --------------------------------------> " << gsl_rng_uniform(RNG) << endl;
    //gsl_rng_set(RNG, randomseed+1);
    cerr << "main rng: " << RNG << endl;

    //cerr << "population file = " << szPopulationFile << endl;
    //cerr << "immunity file = " << szImmunityFile << endl;
    //cerr << "location file = " << szLocationFile << endl;
    //cerr << "network file = " << szNetworkFile << endl;
    //cerr << "runlength = " << nRunLength << endl;
    if (nRunLength>MAXRUNTIME) {
        //cerr << "ERROR: runlength is too long: " << nRunLength << endl;
        //cerr << " change Community.h and recompile." << endl;
        exit(-1);
    }
    if (nRunLength==365) {
        //cerr << "ERROR: you probably want runlength to be 364, not 365" << endl;
        exit(-1);
    }
    //cerr << "random seed = " << randomseed << endl;
    //cerr << "beta_PM = " << betaPM << endl;
    //cerr << "beta_MP = " << betaMP << endl;
    //cerr << "days of complete cross protection = " << nDaysImmune << endl;
    //cerr << "maximum infection parity = " << nMaxInfectionParity << endl;
    //cerr << "pathogenicity of primary infection =";
    for (int i=0; i<NUM_OF_SEROTYPES; i++)
        //cerr << " " << fPrimaryPathogenicity[i];
    //cerr << endl;
    //cerr << "pathogenicity scaling for secondary infection =";
    for (int i=0; i<NUM_OF_SEROTYPES; i++)
        //cerr << " " << fSecondaryScaling[i];
    //cerr << endl;
    //cerr << "mosquito move prob = " << fMosquitoMove << endl;
    //cerr << "mosquito teleport prob = " << fMosquitoTeleport << endl;
    //cerr << "default mosquito capacity per building = " << nDefaultMosquitoCapacity << endl;
    Location::setDefaultMosquitoCapacity(nDefaultMosquitoCapacity);
    if (nSizeMosquitoMultipliers>0) {
        cerr << "mosquito seasonal multipliers (days,mult) =";
        for (int j=0; j<nSizeMosquitoMultipliers; j++)
            cerr << " (" << nMosquitoMultiplierDays[j] << "," << fMosquitoMultipliers[j] << ")";
        cerr << endl;
    }
    //cerr << "Pre-vaccinate fraction = " << fPreVaccinateFraction << endl;
    if (nSizeVaccinate>0) {
        //cerr << "Phased vaccinate (year, age, frac) = ";
        for (int j=0; j<nSizeVaccinate; j++) {
            //cerr << " (" << nVaccinateYear[j] << "," << nVaccinateAge[j]  << "," << fVaccinateFraction[j] << ")";
        }
        //cerr << endl;
    }
    if (nSizePrevaccinateAge>0) {
        //cerr << "Pre-vaccinate by age (min, max, frac) = ";
        for (int j=0; j<nSizePrevaccinateAge; j++) {
            //cerr << " (" << nPrevaccinateAgeMin[j] << "," << nPrevaccinateAgeMax[j]  << "," << fPrevaccinateAgeFraction[j] << ")";
        }
        //cerr << endl;
    }
    //cerr << "VE_S = " << fVES << endl;
    //cerr << "VE_I = " << fVEI << endl;
    //cerr << "VE_P = " << fVEP << endl;
    if (fVES>1.0 || fVEI>1.0 || fVEP>1.0) {
        //cerr << "ERROR: VE_S, VE_I, and VE_P must be between 0 and 1" << endl;
        exit(-1);
    }
    if (fVES<0.0) {
        //cerr << "VE_Ss = " << fVESs[0] << "," << fVESs[1] << "," << fVESs[2] << "," << fVESs[3] << endl;
    }

    if (szPeopleFile.length()>0) {
        //cerr << "people output file = " << szPeopleFile << endl;
    } else {
        //cerr << "no people output file" << endl;
    }
    if (szYearlyPeopleFile.length()>0) {
        //cerr << "yearly people output file = " << szYearlyPeopleFile << endl;
    }
    if (szDailyFile.length()>0) {
        //cerr << "daily output file = " << szDailyFile << endl;
    } else {
        //cerr << "no daily output file" << endl;
    }

    Community community;
    if (!community.loadLocations(szLocationFile.c_str(),szNetworkFile.c_str())) {
        //cerr << "Could not load locations" << endl;
        exit(-1);
    }
    if (!community.loadPopulation(szPopulationFile.c_str(),szImmunityFile.c_str())) {
        //cerr << "Could not load population" << endl;
        exit(-1);
    }
    community.setBetaPM(betaPM);
    community.setBetaMP(betaMP);
    if (fVES<0.0) {
        community.setVESs(fVESs[0],fVESs[1],fVESs[2],fVESs[3]);
    } else {
        community.setVES(fVES);
    }
    community.setVEI(fVEI);
    community.setVEP(fVEP);
    community.setMaxInfectionParity(nMaxInfectionParity);
    Person::setDaysImmune(nDaysImmune);
    community.setMosquitoMoveProbability(fMosquitoMove);
    community.setMosquitoTeleportProbability(fMosquitoTeleport);
    community.setPrimaryPathogenicityScaling(fPrimaryPathogenicity[0],fPrimaryPathogenicity[1],fPrimaryPathogenicity[2],fPrimaryPathogenicity[3]);
    community.setSecondaryPathogenicityScaling(fSecondaryScaling[0],fSecondaryScaling[1],fSecondaryScaling[2],fSecondaryScaling[3]);

    //cerr << community.getNumPerson() << " people" << endl;

    // seed epidemic
    if (nInitialExposed[0]>0 ||
        nInitialExposed[1]>0 ||
        nInitialExposed[2]>0 ||
    nInitialExposed[3]>0) {

        for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
            //cerr << "Initial serotype " << serotype+1 << " exposed = " << nInitialExposed[serotype] << endl;
            for (int i=0; i<nInitialExposed[serotype]; i++)
                community.infect(gsl_rng_uniform_int(RNG, community.getNumPerson()), (Serotype) serotype,0);
        }
    } else if (nInitialInfected[0]>0 ||
        nInitialInfected[1]>0 ||
        nInitialInfected[2]>0 ||
    nInitialInfected[3]>0) {
        for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
            //cerr << "Initial serotype " << serotype+1 << " infected = " << nInitialInfected[serotype] << endl;
            int count = community.getNumInfected(0);
                                                                      
            // must infect nInitialInfected persons
            while (community.getNumInfected(0)<count+nInitialInfected[serotype]) {
                community.infect(gsl_rng_uniform_int(RNG, community.getNumPerson()), (Serotype) serotype,0);
            }
        }
    }
    if (!bSecondaryTransmission) {
        community.setNoSecondaryTransmission();
        //    //cerr << "No secondary transmission, 1 person infected" << endl;
    }

    if (fPreVaccinateFraction>0.0)
        community.vaccinate(fPreVaccinateFraction);

    if (nSizePrevaccinateAge>0) {
        for (int j=0; j<nSizePrevaccinateAge; j++)
            for (int k=nPrevaccinateAgeMin[j]; k<=nPrevaccinateAgeMax[j]; k++)
                community.vaccinate(fPrevaccinateAgeFraction[j],k);
    }

    for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
        nNumInitialSusceptible[serotype] = community.getNumSusceptible((Serotype) serotype);
    }

    // run epidemic
    int nNextMosquitoMultiplier = 0;
    if (bSecondaryTransmission)
        cout << "time,type,id,location,serotype,symptomatic,withdrawn" << endl;
    for (int t=0; t<nRunLength; t++) {
        // phased vaccination
        if ((t%365)==0) {
            int year = (int)(t/365);
            for (int i=0; i<nSizeVaccinate; i++) {
                if (year==nVaccinateYear[i]) {
                    community.vaccinate(fVaccinateFraction[i],nVaccinateAge[i]);
                    //cerr << "vaccinating " << fVaccinateFraction[i]*100 << "% of age " << nVaccinateAge[i] << endl;
                }
            }
        }

        // seed epidemic
        if (nDailyExposed[0]>0 ||
            nDailyExposed[1]>0 ||
            nDailyExposed[2]>0 ||
        nDailyExposed[3]>0) {
            int count = 0;
            int numperson = community.getNumPerson();
            for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++)
                for (int i=0; i<nDailyExposed[serotype]; i++)
                    if (community.infect(gsl_rng_uniform_int(RNG, numperson),(Serotype) serotype,t))
                        count++;
        }

        // mosquito population seasonality?
        if (nSizeMosquitoMultipliers>0 && (t%365)==nMosquitoMultiplierCumulativeDays[nNextMosquitoMultiplier]) {
            community.setMosquitoMultiplier(fMosquitoMultipliers[nNextMosquitoMultiplier]);
            nNextMosquitoMultiplier = (nNextMosquitoMultiplier+1)%nSizeMosquitoMultipliers;
        }

        if (bSecondaryTransmission) {
            // print out infectious mosquitoes
            for (int i=community.getNumInfectiousMosquitoes()-1; i>=0; i--) {
                Mosquito *p = community.getInfectiousMosquito(i);
                cout << t << ",mi," << p->getID() << "," << p->getLocation()->getID() << "," << "," << endl;
            }
            // print out exposed mosquitoes
            for (int i=community.getNumExposedMosquitoes()-1; i>=0; i--) {
                Mosquito *p = community.getExposedMosquito(i);
                // "current" location
                cout << t << ",me," << p->getID() << "," << p->getLocation()->getID() << "," << 1 + (int) p->getSerotype() << "," << "," << endl;
            }
            // print out infected people
            for (int i=community.getNumPerson()-1; i>=0; i--) {
                Person *p = community.getPerson(i);
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

        if (t%100==0)
            //cerr << "Time " << t << endl;

        // output people file
        if (t%365==364 && szYearlyPeopleFile.length()>0) {
            ofstream peopleFile;
            ostringstream ssFilename;
            ssFilename << szYearlyPeopleFile << ((int)(t/365)) << ".csv";
            //cerr << "outputing yearly people information to " << ssFilename.str() << endl;
            peopleFile.open(ssFilename.str().c_str());
            if(peopleFile.fail()) {
                //cerr << "ERROR: People file '" << szPeopleFile << "' cannot be open for writing." << endl;
                return false;
            }
            peopleFile << "pid,serotype,infectiontime,symptomtime,withdrawtime,recoverytime,immdenv1,immdenv2,immdenv3,immdenv4" << endl;
            for (int i=0; i<community.getNumPerson(); i++) {
                Person *p = community.getPerson(i);
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

        community.tick();
    }

    //gsl_rng_free(RNG);

    if (!bSecondaryTransmission) {
        // outputs
        //   number of secondary infections by serotype (4)
        //   number of households infected
        //   age of index case
        //   age(s) of secondary cases
        for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
            cout << (nNumInitialSusceptible[serotype] - community.getNumSusceptible((Serotype) serotype)) << " ";
        }
        //    cout << "secondary infections" << endl;

        int ages[100];
        int times[100];
        int numages=0;
        int indexage=-1;
        int homeids[100];
        int numhomes=0;
        // ages of infected people
        for (int i=community.getNumPerson()-1; i>=0; i--) {
            Person *p = community.getPerson(i);
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
    if (szDailyFile.length()>0) {
        //cerr << "outputing daily infected/symptomatic information to " << szDailyFile << endl;
        ofstream dailyFile;
        dailyFile.open(szDailyFile.c_str());
        if(dailyFile.fail()) {
            //cerr << "ERROR: Daily file '" << szDailyFile << "' cannot be open for writing." << endl;
            return false;
        }
        dailyFile << "day,newly infected DENV1,newly infected DENV2,newly infected DENV3,newly infected DENV4,newly symptomatic DENV1,newly symptomatic DENV2,newly symptomatic DENV3,newly symptomatic DENV4" << endl;
        const int *pInfected1 = community.getNumNewlyInfected(SEROTYPE_1);
        const int *pSymptomatic1 = community.getNumNewlySymptomatic(SEROTYPE_1);
        const int *pInfected2 = community.getNumNewlyInfected(SEROTYPE_2);
        const int *pSymptomatic2 = community.getNumNewlySymptomatic(SEROTYPE_2);
        const int *pInfected3 = community.getNumNewlyInfected(SEROTYPE_3);
        const int *pSymptomatic3 = community.getNumNewlySymptomatic(SEROTYPE_3);
        const int *pInfected4 = community.getNumNewlyInfected(SEROTYPE_4);
        const int *pSymptomatic4 = community.getNumNewlySymptomatic(SEROTYPE_4);
        for (int t=0; t<nRunLength; t++) {
            dailyFile << t << ","
                      << pInfected1[t] << "," 
                      << pInfected2[t] << "," 
                      << pInfected3[t] << "," 
                      << pInfected4[t] << "," 
                      << pSymptomatic1[t] << "," 
                      << pSymptomatic2[t] << "," 
                      << pSymptomatic3[t] << "," 
                      << pSymptomatic4[t] << endl;
        }
        dailyFile.close();
    }

    // output people file
    if (szPeopleFile.length()>0) {
        //cerr << "outputing people information to " << szPeopleFile << endl;
        ofstream peopleFile;
        peopleFile.open(szPeopleFile.c_str());
        if(peopleFile.fail()) {
            //cerr << "ERROR: People file '" << szPeopleFile << "' cannot be open for writing." << endl;
            return false;
        }
        peopleFile << "pid,serotype,infectiontime,symptomtime,withdrawtime,recoverytime,immdenv1,immdenv2,immdenv3,immdenv4,vaccinated" << endl;
        for (int i=0; i<community.getNumPerson(); i++) {
            Person *p = community.getPerson(i);
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

    return 0;
}
