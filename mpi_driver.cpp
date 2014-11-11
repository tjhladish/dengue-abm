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
using namespace dengue::util;
const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);

// Predeclare local functions
Community* build_community(const Parameters* par);
void seed_epidemic(const Parameters* par, Community* community);
vector<int> simulate_epidemic(const Parameters* par, Community* community, const int process_id=0);
void write_output(const Parameters* par, Community* community, vector<int> initial_susceptibles);

int main(int argc, char *argv[]) {
    time_t start ,end;
    time (&start);

    srand(time(NULL));
    int proccess_id = rand();
    fprintf(stderr, "%dbegin\n", proccess_id);

    const Parameters* par = new Parameters(argc, argv);

    gsl_rng_set(RNG, par->randomseed);
    Community* community = build_community(par);
    vector<int> initial_susceptibles = community->getNumSusceptible();
    seed_epidemic(par, community);
    vector<int> epi_sizes = simulate_epidemic(par, community);

    const double ef = par->expansionFactor;
    vector<double> x(epi_sizes.size());
    vector<double> y(epi_sizes.size());
    // convert infections to cases
    for (unsigned int i = 0; i < epi_sizes.size(); i++) { 
        y[i] = ((double) epi_sizes[i])/ef; 
        x[i] = i+1.0;
    }

    Fit* fit = lin_reg(x, y);
    
    time (&end);
    double dif = difftime (end,start);

    stringstream ss;
    // wallclock time (seconds)
    ss << proccess_id << "end " << dif << " ";
    // parameters
    ss << par->expansionFactor << " " << par->fMosquitoMove << " " << par->nDailyExposed[0][0] << " "
       << par->betaMP << " " << par->betaPM << " ";
    // metrics
    ss << mean(y) << " " << stdev(y) << " " << max_element(y) << " "
       << fit->m << " " << fit->b << " " << fit->rsq;

    ss << endl;
    string output = ss.str();
    fputs(output.c_str(), stderr);
    
    // metrics to stdout
    cout << mean(y) << " " << stdev(y) << " " << max_element(y) << " "
       << fit->m << " " << fit->b << " " << fit->rsq;


    write_output(par, community, initial_susceptibles);
   
    return 0;
}
 
Community* build_community(const Parameters* par) {
    Community* community = new Community(par);

    if (!community->loadLocations(par->locationFilename, par->networkFilename)) {
        cerr << "ERROR: Could not load locations" << endl;
        exit(-1);
    }
    if (!community->loadPopulation(par->populationFilename, par->immunityFilename, par->swapProbFilename)) {
        cerr << "ERROR: Could not load population" << endl;
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


void write_yearly_people_file(const Parameters* par, Community* community, int time) {
    ofstream yearlyPeopleOutputFile;
    ostringstream ssFilename;
    ssFilename << par->yearlyPeopleOutputFilename << ((int)(time/365)) << ".csv";
    cerr << "outputing yearly people information to " << ssFilename.str() << endl;
    yearlyPeopleOutputFile.open(ssFilename.str().c_str());
    if(yearlyPeopleOutputFile.fail()) {
        cerr << "ERROR: People file '" << par->yearlyPeopleOutputFilename << "' cannot be open for writing." << endl;
        exit(-1);
    }
    yearlyPeopleOutputFile << "pid,serotype,infectiontime,symptomtime,withdrawtime,recoverytime,immdenv1,immdenv2,immdenv3,immdenv4" << endl;
    for (int i=0; i<community->getNumPerson(); i++) {
        Person *p = community->getPerson(i);
        for (int j=p->getNumInfections()-1; j>=0; j--) {
            yearlyPeopleOutputFile << p->getID() << "," 
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
    yearlyPeopleOutputFile.close();
    return;
}

vector<int> simulate_epidemic(const Parameters* par, Community* community, const int process_id) {
    vector<int> epi_sizes;
    int nextMosquitoMultiplierIndex = 0;
    const int mosquitoMultiplierTotalDuration = par->mosquitoMultipliers.size() > 0 ?
                                                par->mosquitoMultipliers.back().start + par->mosquitoMultipliers.back().duration : 0;
    int nextEIPindex = 0;
    const int EIPtotalDuration = par->extrinsicIncubationPeriods.size() > 0 ?
                                 par->extrinsicIncubationPeriods.back().start + par->extrinsicIncubationPeriods.back().duration : 0;
    if (par->bSecondaryTransmission and not par->abcVerbose) cout << "time,type,id,location,serotype,symptomatic,withdrawn,new_infection" << endl;

    vector<int> daily_incidence(3,0); // { introductions, local transmission, total }
    vector<int> weekly_incidence(3,0);
    vector<int> monthly_incidence(3,0);
    vector<int> yearly_incidence(3,0);
    pair<int, int> current_month(0, DAYS_IN_MONTH_CUM[0]);
    for (int t=0; t<par->nRunLength; t++) {
        int year = (int)(t/365);

        // phased vaccination
        if ((t%365)==0) {
            for (int i=0; i<par->nSizeVaccinate; i++) {
                if (year==par->nVaccinateYear[i]) {
                    community->vaccinate(par->fVaccinateFraction[i],par->nVaccinateAge[i]);
                    if (not par->abcVerbose) {
                        cerr << "vaccinating " << par->fVaccinateFraction[i]*100 << "% of age " << par->nVaccinateAge[i] << endl;
                    }
                }
            }
        }

        // seed epidemic
        {
            int numperson = community->getNumPerson();
            for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
                const int ai_year_lookup = year % par->annualIntroductions.size();
                const double intros = par->annualIntroductions[ai_year_lookup];
                const int de_year_lookup = year % par->nDailyExposed.size();
                const double serotype_weight = par->nDailyExposed[de_year_lookup][serotype];
                const double annual_intros_weight = par->annualIntroductionsCoef;
                const double expected_num_exposed = serotype_weight * annual_intros_weight * intros; 
                if (expected_num_exposed <= 0) continue;
                const int num_exposed = gsl_ran_poisson(RNG, expected_num_exposed);
                for (int i=0; i<num_exposed; i++) {
                    // gsl_rng_uniform_int returns on [0, numperson-1]
                    int transmit_to_id = gsl_rng_uniform_int(RNG, numperson) + 1; 
                    if (community->infect(transmit_to_id, (Serotype) serotype, t)) {
                        daily_incidence[0]++;
                    }
                }
            }
        }

        // should the mosquito population change?
        if (par->mosquitoMultipliers.size()>0 && (t%mosquitoMultiplierTotalDuration)==par->mosquitoMultipliers[nextMosquitoMultiplierIndex].start) {
            community->setMosquitoMultiplier(par->mosquitoMultipliers[nextMosquitoMultiplierIndex].value);
            nextMosquitoMultiplierIndex = (nextMosquitoMultiplierIndex+1)%par->mosquitoMultipliers.size();
        }

        // should the EIP change?
        if (par->extrinsicIncubationPeriods.size()>0 && (t%EIPtotalDuration)==par->extrinsicIncubationPeriods[nextEIPindex].start) {
            community->setExtrinsicIncubation(par->extrinsicIncubationPeriods[nextEIPindex].value);
            nextEIPindex = (nextEIPindex+1)%par->extrinsicIncubationPeriods.size();
        }

        community->tick(t);

        for (int i=community->getNumPerson()-1; i>=0; i--) {
            Person *p = community->getPerson(i);
            if (p->isInfected(t) and p->isNewlyInfected(t)) {
                ++daily_incidence[2];
            }
        }

        daily_incidence[1] = daily_incidence[2] - daily_incidence[0]; // local transmission = total - introductions
        if (par->dailyOutput) {
            cerr << "day: " << t << " "; for (auto v: daily_incidence) cerr << v << " "; 
            cerr << community->getExtrinsicIncubation() << " " << community->getMosquitoMultiplier()*par->nDefaultMosquitoCapacity << endl;
            //cerr << "day: " << t << " "; for (auto v: daily_incidence) cerr << v << " "; cerr << endl;
        }

        if (par->weeklyOutput) {
            for (unsigned int i = 0; i < daily_incidence.size(); ++i) weekly_incidence[i] += daily_incidence[i];
            if ((t+1)%7==0) {
                cerr << "week: " << (t+1)/7 << " "; for (auto v: weekly_incidence) cerr << v << " "; cerr << endl;
                weekly_incidence = {0,0,0};
            }
        }

        if (par->monthlyOutput) {
            for (unsigned int i = 0; i < daily_incidence.size(); ++i) monthly_incidence[i] += daily_incidence[i];
            if ((t+1)%365==current_month.second) {
                cerr << "month: " << current_month.first+1 << " "; for (auto v: monthly_incidence) cerr << v << " "; cerr << endl;
                current_month.first = current_month.first == 11 ? 0 : current_month.first + 1; // loop DEC -> JAN
                current_month.second = DAYS_IN_MONTH_CUM[current_month.first];
                monthly_incidence = {0,0,0};
            }
        }

        // handle several things that happen yearly
        for (unsigned int i = 0; i < daily_incidence.size(); ++i) yearly_incidence[i] += daily_incidence[i];
        if ((t+1)%365 == 0) {
            epi_sizes.push_back(yearly_incidence[2]);
            if (par->yearlyPeopleOutputFilename.length() > 0) {
                write_yearly_people_file(par, community, t);
            }
            if (par->yearlyOutput) {
                cerr << "year: " << (t+1)/365 << " "; for (auto v: yearly_incidence) cerr << v << " "; cerr << endl;
                yearly_incidence = {0,0,0};
            }
        }

        if (t%1000==0 and par->abcVerbose) {
            cout << hex << process_id << dec << " T: " << t << " daily: " << daily_incidence[2] << " annual: " << yearly_incidence[2] << endl;
        }

        daily_incidence = {0,0,0};
    }
    return epi_sizes;
}

void daily_detailed_output(Community* community, int t) {
    // print out infectious mosquitoes
    for (int i=community->getNumInfectiousMosquitoes()-1; i>=0; i--) {
        Mosquito *p = community->getInfectiousMosquito(i);
        cout << t << ",mi," << p->getID() << "," << p->getLocation()->getID() << "," << "," << "," << endl;
    }
    // print out exposed mosquitoes
    for (int i=community->getNumExposedMosquitoes()-1; i>=0; i--) {
        Mosquito *p = community->getExposedMosquito(i);
        // "current" location
        cout << t << ",me," << p->getID() << "," << p->getLocation()->getID() << "," << 1 + (int) p->getSerotype() << "," << "," << "," << endl;
    }
    // print out infected people
    for (int i=community->getNumPerson()-1; i>=0; i--) {
        Person *p = community->getPerson(i);
        if (p->isInfected(t)) {
            // home location
            cout << t 
                 << ",p,"
                 << p->getID() << "," 
                 << p->getLocation(0)->getID() << "," 
                 << 1 + (int) p->getSerotype() << "," 
                 << (p->isSymptomatic(t)?1:0) << "," 
                 << (p->isWithdrawn(t)?1:0) << ","
                 << (p->isNewlyInfected(t)?1:0) << endl;
        }
    }
}

void write_output(const Parameters* par, Community* community, vector<int> numInitialSusceptible) {
    if (!par->bSecondaryTransmission) {
        // outputs
        //   number of secondary infections by serotype (4)
        //   number of households infected
        //   age of index case
        //   age(s) of secondary cases
        vector<int> numCurrentSusceptible = community->getNumSusceptible();
        for (int s=0; s<NUM_OF_SEROTYPES; s++) {
            cout << (numInitialSusceptible[s] - numCurrentSusceptible[s]) << " ";
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
    if (par->dailyOutputFilename.length()>0) {
        cerr << "outputing daily infected/symptomatic information to " << par->dailyOutputFilename << endl;
        ofstream dailyOutputFile;
        dailyOutputFile.open(par->dailyOutputFilename.c_str());
        if(dailyOutputFile.fail()) {
            cerr << "ERROR: Daily file '" << par->dailyOutputFilename << "' cannot be open for writing." << endl;
            exit(-1);
        }
        dailyOutputFile << "day,newly infected DENV1,newly infected DENV2,newly infected DENV3,newly infected DENV4,"
                  << "newly symptomatic DENV1,newly symptomatic DENV2,newly symptomatic DENV3,newly symptomatic DENV4" << endl;
        vector< vector<int> > infected =    community->getNumNewlyInfected();
        vector< vector<int> > symptomatic = community->getNumNewlySymptomatic();
        for (int t=0; t<par->nRunLength; t++) {
            dailyOutputFile << t << ",";
            for (int i=0; i<NUM_OF_SEROTYPES; i++)   dailyOutputFile << infected[i][t] << ",";
            for (int i=0; i<NUM_OF_SEROTYPES-1; i++) dailyOutputFile << symptomatic[i][t] << ","; 
            dailyOutputFile << symptomatic[NUM_OF_SEROTYPES-1][t] << endl;
        }
        dailyOutputFile.close();
    }

    // output people file
    if (par->peopleOutputFilename.length()>0) {
        cerr << "outputing people information to " << par->peopleOutputFilename << endl;
        ofstream peopleOutputFile;
        peopleOutputFile.open(par->peopleOutputFilename.c_str());
        if(peopleOutputFile.fail()) {
            cerr << "ERROR: People file '" << par->peopleOutputFilename << "' cannot be open for writing." << endl;
            exit(-1);
        }
        peopleOutputFile << "pid,serotype,infectiontime,symptomtime,withdrawtime,recoverytime,immdenv1,immdenv2,immdenv3,immdenv4,vaccinated" << endl;
        for (int i=0; i<community->getNumPerson(); i++) {
            Person *p = community->getPerson(i);
            for (int j=p->getNumInfections()-1; j>=0; j--) {
                peopleOutputFile << p->getID() << "," 
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
        peopleOutputFile.close();
    }
}
