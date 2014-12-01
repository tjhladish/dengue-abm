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
using namespace dengue::util;
const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);

// Predeclare local functions
Community* build_community(const Parameters* par);
void seed_epidemic(const Parameters* par, Community* community);
vector<int> simulate_epidemic(const Parameters* par, Community* community, const int process_id=0);
void write_immunity_file(const Parameters* par, const Community* community, const int int_label, string filename="");
void write_immunity_by_age_file(const Parameters* par, const Community* community, const int year, string filename="");
void write_output(const Parameters* par, Community* community, vector<int> initial_susceptibles);
void write_daily_buffer( vector<string>& buffer, const int process_id);

Community* build_community(const Parameters* par) {
    Community* community = new Community(par);
    Person::setPar(par);

    if (!community->loadLocations(par->locationFilename, par->networkFilename)) {
        cerr << "ERROR: Could not load locations" << endl;
        exit(-1);
    }
    if (!community->loadPopulation(par->populationFilename, par->immunityFilename, par->swapProbFilename)) {
        cerr << "ERROR: Could not load population" << endl;
        exit(-1);
    }

    if (!par->abcVerbose) {
        cerr << community->getNumPerson() << " people" << endl;
    }

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
    // epidemic may be seeded with initial exposure OR initial infection 
    bool attempt_initial_infection = true;
    for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
        // Normal usage, to simulate epidemic
        if (par->nInitialExposed[serotype] > 0) {
            attempt_initial_infection = false;
            for (int i=0; i<par->nInitialExposed[serotype]; i++)
                community->infect(gsl_rng_uniform_int(RNG, community->getNumPerson()), (Serotype) serotype,0);
        }
    }
    if (attempt_initial_infection) {
        for (int serotype=0; serotype<NUM_OF_SEROTYPES; serotype++) {
            // Useful for estimating R0
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
                << p->getRecoveryTime(j) << ",";
            for (int s = 0; s < NUM_OF_SEROTYPES - 1; ++s) {
                yearlyPeopleOutputFile << (p->isSusceptible((Serotype) s)?0:1) << ","; 
            }
                yearlyPeopleOutputFile << (p->isSusceptible((Serotype) (NUM_OF_SEROTYPES - 1))?0:1) << endl;
        }
    }
    yearlyPeopleOutputFile.close();
    return;
}

void initialize_seasonality(const Parameters* par, Community* community, int& nextMosquitoMultiplierIndex, int& nextEIPindex, pair<int,int>& current_month) {
    const int desiredDayOfYearOffset = par->startDayOfYear - 1; // yields 0 for Jan 1

    const int mosquitoMultiplierTotalDuration = par->getMosquitoMultiplierTotalDuration();
    int currentDayOfYearOffset = 0;
    if (mosquitoMultiplierTotalDuration > 0) {
        nextMosquitoMultiplierIndex = 0;
        while (currentDayOfYearOffset <= desiredDayOfYearOffset) {
            currentDayOfYearOffset += par->mosquitoMultipliers[nextMosquitoMultiplierIndex].duration;
            if (currentDayOfYearOffset > desiredDayOfYearOffset) {
                const double mm = par->mosquitoMultipliers[nextMosquitoMultiplierIndex].value;
                community->setMosquitoMultiplier(mm);
            }
            nextMosquitoMultiplierIndex = (nextMosquitoMultiplierIndex+1)%par->mosquitoMultipliers.size();
        }
    }

    const int EIPtotalDuration = par->getEIPtotalDuration();
    currentDayOfYearOffset = 0;
    if (EIPtotalDuration > 0) {
        nextEIPindex = 0;
        while (currentDayOfYearOffset <= desiredDayOfYearOffset) {
            currentDayOfYearOffset += par->extrinsicIncubationPeriods[nextEIPindex].duration;
            if (currentDayOfYearOffset > desiredDayOfYearOffset) {
                const double eip = par->extrinsicIncubationPeriods[nextEIPindex].value;
                community->setExtrinsicIncubation(eip);
            }
            nextEIPindex = (nextEIPindex+1)%par->extrinsicIncubationPeriods.size();
        }
    }

    currentDayOfYearOffset = 0;
    for (int month=0; month<12; ++month) {
        if (DAYS_IN_MONTH_CUM[month] > desiredDayOfYearOffset) {
            current_month.first = month;
            current_month.second = DAYS_IN_MONTH_CUM[month];
            //int day_of_month = month==0 ? desiredDayOfYearOffset + 1 : desiredDayOfYearOffset - DAYS_IN_MONTH_CUM[month-1] + 1;
            //cerr << "Start date: " << MONTH_NAMES[month] << " " << day_of_month << endl;
            break;
        }
    }
}

vector<int> simulate_epidemic(const Parameters* par, Community* community, const int process_id) {
    vector<int> epi_sizes;
    const int dayOfYearOffset = par->startDayOfYear - 1; // yields 0 for Jan 1
    int nextMosquitoMultiplierIndex = 0;
    const int mosquitoMultiplierTotalDuration = par->getMosquitoMultiplierTotalDuration();
    int nextEIPindex = 0;
    const int EIPtotalDuration = par->getEIPtotalDuration();
    pair<int, int> current_month(0, DAYS_IN_MONTH_CUM[0]);
    initialize_seasonality(par, community, nextMosquitoMultiplierIndex, nextEIPindex, current_month);
    vector<string> daily_output_buffer;
    if (par->bSecondaryTransmission and not par->abcVerbose) {
        daily_output_buffer.push_back("time,type,id,location,serotype,symptomatic,withdrawn,new_infection");
        //cout << "time,type,id,location,serotype,symptomatic,withdrawn,new_infection" << endl;
    }

    vector<int> daily_incidence(3,0); // { introductions, local transmission, total }
    vector<int> weekly_incidence(3,0);
    vector<int> monthly_incidence(3,0);
    vector<int> yearly_incidence(3,0);
    for (int t=0; t<par->nRunLength; t++) {
        int year = (int)(t/365);
        //cerr << "Sim day: " << t << endl;

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
        if (par->mosquitoMultipliers.size() > 0) {
            const int nextMosquitoStart = par->mosquitoMultipliers[nextMosquitoMultiplierIndex].start;
            if ( ((t+dayOfYearOffset)%mosquitoMultiplierTotalDuration) == nextMosquitoStart) {
                //cerr << "updating mosquitoes on day " << t << ", which is day " << t+dayOfYearOffset+1 
                //     << " of the year. Using index " << nextMosquitoMultiplierIndex << endl;
                community->setMosquitoMultiplier(par->mosquitoMultipliers[nextMosquitoMultiplierIndex].value);
                nextMosquitoMultiplierIndex = (nextMosquitoMultiplierIndex+1)%par->mosquitoMultipliers.size();
            }
        }

        // should the EIP change?
        if (par->extrinsicIncubationPeriods.size() > 0) {
            const int nextEIPstart = par->extrinsicIncubationPeriods[nextEIPindex].start;
            if ( ((t+dayOfYearOffset)%EIPtotalDuration) == nextEIPstart) {
                community->setExtrinsicIncubation(par->extrinsicIncubationPeriods[nextEIPindex].value);
                nextEIPindex = (nextEIPindex+1)%par->extrinsicIncubationPeriods.size();
            }
        }

        community->tick(t);

        for (int i=community->getNumPerson()-1; i>=0; i--) {
            Person *p = community->getPerson(i);
            if (p->isInfected(t) and p->isNewlyInfected(t)) {
                ++daily_incidence[2];
                /*
                // this is commented out because we don't usually output daily data,
                // and there's no point wasting ~100 mb of ram per process on long runs
                stringstream ss;
                ss << t << ","
                   << p->getID() << ","
                   << p->getLocation(0)->getID() << ","
                   << 1 + (int) p->getSerotype() << ","
                   << (p->isSymptomatic(t)?1:0);
                daily_output_buffer.push_back(ss.str());
                */
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
            if ((t+dayOfYearOffset+1)%365==current_month.second) {
/*{
  epi_sizes.push_back(monthly_incidence[2]);
  string imm_filename;
  stringstream ss_imm_filename;
  //we're starting in april, want to start numbering at 1
  int month = (current_month.first+8)%12 + 1;
  ss_imm_filename << "immunity." << process_id << ".year" << (t+1)/365 << ".mon" << month;
  imm_filename = ss_imm_filename.str();
  int dummy = 0;
  write_immunity_by_age_file(par, community, dummy, imm_filename);
}*/
                cerr << "month: " << current_month.first+1 << " "; for (auto v: monthly_incidence) cerr << v << " "; cerr << endl;
                current_month.first = current_month.first == 11 ? 0 : current_month.first + 1; // loop DEC -> JAN
                current_month.second = DAYS_IN_MONTH_CUM[current_month.first];
                monthly_incidence = {0,0,0};
            }
        }

        // handle several things that happen yearly
        for (unsigned int i = 0; i < daily_incidence.size(); ++i) yearly_incidence[i] += daily_incidence[i];

        if ((t+1)%365==0 and par->abcVerbose) {
            cout << hex << process_id << dec << " T: " << t << " daily: " << daily_incidence[2] << " annual: " << yearly_incidence[2] << endl;
        }

        if ((t+1)%365 == 0) {
            const int year = (t+1)/365;
            //epi_sizes.push_back(yearly_incidence[2]);
            if (par->yearlyPeopleOutputFilename.length() > 0) {
                write_yearly_people_file(par, community, t);
            }
            if (par->yearlyOutput) {
                cerr << "year: " << year << " "; for (auto v: yearly_incidence) cerr << v << " "; cerr << endl;
            }
            yearly_incidence = {0,0,0};
        }

        daily_incidence = {0,0,0};
    }

    //write_daily_buffer(daily_output_buffer, process_id);
    return epi_sizes;
}


void write_daily_buffer( vector<string>& buffer, const int process_id) {
    stringstream ss_filename;
    ss_filename << "daily_output." << process_id;
    ofstream file;
    file.open(ss_filename.str());

    for (string s: buffer) file << s << endl;
    file.close();
}


void write_immunity_by_age_file(const Parameters* par, const Community* community, const int year, string filename) {
    if (filename == "") {
        stringstream ss_filename;
        ss_filename << "imm_vs_age.year" << year;
        filename = ss_filename.str();
    }
    // total count, denv1, denv2, denv3, denv4, imm_to_one, imm_to_all
    vector< vector<int> > tally(NUM_AGE_CLASSES, vector<int>(NUM_OF_SEROTYPES+3, 0)); 
    for (int i = 0; i<community->getNumPerson(); ++i) {
        Person* p = community->getPerson(i);
        const int age = p->getAge();
        tally[age][0]++;
        const int numInfections = p->getNumInfections();
        for (int k = 0; k<numInfections; ++k) {
            const int s = (int) p->getSerotype(k);
            tally[age][s+1]++;
        }
        if (numInfections > 0) tally[age][NUM_OF_SEROTYPES+1]++;
        if (numInfections == NUM_OF_SEROTYPES) tally[age][NUM_OF_SEROTYPES+2]++;
    }

    ofstream file;
    file.open(filename);
    for (int a = 0; a<NUM_AGE_CLASSES; ++a) {
        file << a;
        for (int i: tally[a]) file << " " << i;
        file << endl;
    }
    file.close();
}


void write_immunity_file(const Parameters* par, const Community* community, const int int_label, string filename) {
    if (filename == "") {
        stringstream ss_filename;
        ss_filename << "immunity." << int_label;
        filename = ss_filename.str();
    }
    ofstream file;
    file.open(filename);
    file << "pid age imm1 imm2 imm3 imm4\n";
    for (int i = 0; i<community->getNumPerson(); ++i) {
        Person* p = community->getPerson(i);
        vector<int> infection_history(NUM_OF_SEROTYPES, 0); // 0 is no infection; 1 means last year, 2 means 2 years ago ...
//if (i == 0) cerr << "person 0 num infections: " << p->getNumInfections()  << " in year " << year << endl;
        for (int k = 0; k<p->getNumInfections(); ++k) {
            int s = (int) p->getSerotype(k);
            infection_history[s] = p->getInfectedTime(k) - par->nRunLength;
            // how many years ago did this person get infected?  within last 365 days = 1 year ago
            // infection time is relative to end of simulation of par->nRunLength days
            //infection_history[s] = 1 + (int) (par->nRunLength - p->getInfectedTime(k))/365;
//if (i == 0) cerr << "infection time: " << p->getInfectedTime()  << " output as: " << infection_history[s] << endl;
        }
        file << p->getID() << " " << p->getAge() << " ";
        for (auto sero: infection_history) file << sero << " ";
        file << endl;
    }
    file.close();
}


void daily_detailed_output(Community* community, int t) {
    // print out infectious mosquitoes
/*    for (int i=community->getNumInfectiousMosquitoes()-1; i>=0; i--) {
        Mosquito *p = community->getInfectiousMosquito(i);
        cout << t << ",mi," << p->getID() << "," << p->getLocation()->getID() << "," << "," << "," << endl;
    }
    // print out exposed mosquitoes
    for (int i=community->getNumExposedMosquitoes()-1; i>=0; i--) {
        Mosquito *p = community->getExposedMosquito(i);
        // "current" location
        cout << t << ",me," << p->getID() << "," << p->getLocation()->getID() << "," << 1 + (int) p->getSerotype() << "," << "," << "," << endl;
    }*/
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
                    << p->getRecoveryTime(j) << ","; 
                for (int s = 0; s < NUM_OF_SEROTYPES; ++s) {
                    peopleOutputFile << (p->isSusceptible((Serotype) s)?0:1) << ","; 
                }
                    peopleOutputFile << (p->isVaccinated()?1:0) << endl;
            }
        }
        peopleOutputFile.close();
    }
}
