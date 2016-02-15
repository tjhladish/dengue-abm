#include "Parameters.h"
#include "Location.h"
#include "Utility.h"
#include <fstream>
#include <sstream>
#include <numeric>
#include <gsl/gsl_randist.h>
#include <climits> // INT_MAX

void Parameters::define_defaults() {
    serial = 0;
    randomseed = 5489;
    nRunLength = 100;
    betaPM = 0.2;
    betaMP = 0.1;
    fMosquitoMove = 0.2;
    mosquitoMoveModel = "weighted";
    fMosquitoTeleport = 0.0;
    fVESs = vector<double>(NUM_OF_SEROTYPES, 0.7);
    fVESs_NAIVE.clear();
    fVEI = 0.0;
    fVEP = 0.0;
    fVEH = 0.803;
    primarySevereFraction    = vector<double>(NUM_OF_SEROTYPES, 0.0);  // never severe, except possibly for infants
    secondarySevereFraction  = vector<double>(NUM_OF_SEROTYPES, 0.5);
    tertiarySevereFraction   = vector<double>(NUM_OF_SEROTYPES, 0.0);
    quaternarySevereFraction = vector<double>(NUM_OF_SEROTYPES, 0.0);
    hospitalizedFraction = {0.0, 0.15, 0.9};            // rough estimates in mex
    bVaccineLeaky = false;
    nDefaultMosquitoCapacity = 50;                      // mosquitoes per location
    eMosquitoDistribution = CONSTANT;
    bSecondaryTransmission = true;
    populationFilename = "population.txt";
    immunityFilename = "";
    networkFilename = "network.txt";
    locationFilename = "locations.txt";
    peopleOutputFilename = "";
    yearlyPeopleOutputFilename = "";
    dailyOutputFilename = "";
    swapProbFilename = "";
    annualIntroductionsFilename = "";                   // time series of some external factor determining introduction rate
    annualIntroductionsCoef = 1;                        // multiplier to rescale external introductions to something sensible
    normalizeSerotypeIntros = false;
    annualIntroductions = {1.0};
    nDaysImmune = 365;
    reportedFraction = {0.0, 0.05, 1.0};                // fraction of asymptomatic, mild, and severe cases reported
    nDailyExposed.push_back(vector<float>(NUM_OF_SEROTYPES, 0.0)); // default: no introductions
    annualSerotypeFilename = "";

    extrinsicIncubationPeriods.clear();
    extrinsicIncubationPeriods.emplace_back(0, 1, 11);  // default: 11 days (Nishiura & Halstead 2007), starting on day 0 with (reused) duration 1 day
    dailyEIPfilename = "";
    simpleEIP = false;                                  // default: sample EIPs from a log-normal distribution, using expected incubation periods (Chan & Johanson 2012)
                                                        // 'true' means use EIPs literally as provided (all mosquitoes infected on day X have same EIP)
    nInitialExposed  = vector<int>(NUM_OF_SEROTYPES, 0);
    nInitialInfected = vector<int>(NUM_OF_SEROTYPES, 0);

    basePathogenicity = 1.0;                            // preferably this value is fit; default interpretation is Pr{symptomatic | secondary denv1 infection}
    postSecondaryRelativeRisk = 0.1;                    // risk of symptoms in post-secondary infections relative to secondary infections

    primaryPathogenicityModel = ORIGINAL_LOGISTIC;
    primaryRelativeRisk= 0.5;
    annualFlavivirusAttackRate = 0.05;                  // used for geometric pathogenicity model

    // Probabilities GIVEN maternal antibodies (from random cohabitating female of reproductive age)
    // Values estimated based on Fig. 5 of Halstead et al, Dengue hemorrhagic fever in infants: research opportunities ignored, EID, 2002
    // http://stacks.cdc.gov/view/cdc/13972/cdc_13972_DS1.pdf
    infantImmuneProb = 3.0/12.0;    // given exposure and maternal antibodies, chance of resisting infection (first 3 months)
    infantSevereProb = 6.0/9.0;     // given infection and maternal antibodies, chance of severe disease (next 6 months)

    vaccinationEvents.clear();
    numVaccineDoses = 3;
    vaccineDoseInterval = 182;

    linearlyWaningVaccine = false;
    vaccineImmunityDuration = INT_MAX;
    vaccineBoosting = false;
    vaccineBoostingInterval = 730;
    bRetroactiveMatureVaccine = false;

    const vector<float> MOSQUITO_MULTIPLIER_DEFAULTS = {0.179, 0.128, 0.123, 0.0956, 0.195, 0.777, 0.940, 0.901, 1.0, 0.491, 0.301, 0.199};
    mosquitoMultipliers.clear();
    mosquitoMultipliers.resize(DAYS_IN_MONTH.size());
    int running_sum = 0;
    for (unsigned int j=0; j<mosquitoMultipliers.size(); j++) {
        mosquitoMultipliers[j].start = running_sum;
        mosquitoMultipliers[j].duration = DAYS_IN_MONTH[j];
        mosquitoMultipliers[j].value = MOSQUITO_MULTIPLIER_DEFAULTS[j];
        running_sum += DAYS_IN_MONTH[j];
    }

    startDayOfYear = 1;

    dailyOutput = false;
    weeklyOutput = false;
    monthlyOutput = false;
    yearlyOutput = false;
    abcVerbose = false;
}

void Parameters::readParameters(int argc, char* argv[]) {
    cerr << "Dengue model, Version " << VERSION_NUMBER_MAJOR << "." << VERSION_NUMBER_MINOR << endl;
    cerr << "written by Dennis Chao and Thomas Hladish in 2012-2014" << endl;

    if (argc>1) {
        for (int i=1; i<argc; i++) {
            char** end = NULL;
            if (strcmp(argv[i], "-randomseed")==0) {
                randomseed = strtol(argv[++i],end,10);
            }
            else if (strcmp(argv[i], "-runlength")==0) {
                nRunLength = strtol(argv[++i],end,10);
            }
            else if (strcmp(argv[i], "-initialinfected")==0) {
                for (int j=0; j<NUM_OF_SEROTYPES; j++) nInitialInfected[j]=strtol(argv[++i],end,10);
            }
            else if (strcmp(argv[i], "-initialexposed")==0) {
                for (int j=0; j<NUM_OF_SEROTYPES; j++) nInitialExposed[j]=strtol(argv[++i],end,10);
            }
            else if (strcmp(argv[i], "-dailyexposed")==0) {
                if (annualSerotypeFilename == "") {
                    for (int j=0; j<NUM_OF_SEROTYPES; j++) nDailyExposed[0][j]=strtod(argv[++i],end);
                } else {
                    // This warning doesn't get thrown if -dailyexposed comes before the -annualserotypefile
                    cerr << "WARNING: Annual serotype file specified.  Ignoring daily exposed parameter.\n";
                }
            }
            else if (strcmp(argv[i], "-annualintroscoef")==0) {
                annualIntroductionsCoef = strtod(argv[++i],end);
            }
            else if (strcmp(argv[i], "-betapm")==0) {
                betaPM = strtod(argv[++i],end);
            }
            else if (strcmp(argv[i], "-betamp")==0) {
                betaMP = strtod(argv[++i],end);
            }
            else if (strcmp(argv[i], "-mosquitomove")==0) {
                fMosquitoMove = strtod(argv[++i],end);
            }
            else if (strcmp(argv[i], "-mosquitomovemodel")==0) {
                mosquitoMoveModel = argv[++i];
            }
            else if (strcmp(argv[i], "-mosquitoteleport")==0) {
                fMosquitoTeleport = strtod(argv[++i],end);
            }
            else if (strcmp(argv[i], "-mosquitocapacity")==0) {
                nDefaultMosquitoCapacity = strtol(argv[++i],end,10);
            }
            else if (strcmp(argv[i], "-mosquitodistribution")==0) {
                const char* argstr = {argv[++i]};
                if (strcmp(argstr, "constant")==0) {
                    eMosquitoDistribution = CONSTANT;
                } else if (strcmp(argstr, "exponential")==0) {
                    eMosquitoDistribution = EXPONENTIAL;
                } else {
                    cerr << "ERROR: Invalid mosquito distribution specified." << endl;
                    exit(-1);
                }
            }
            else if (strcmp(argv[i], "-mosquitomultipliers")==0) {
                mosquitoMultipliers.clear();
                mosquitoMultipliers.resize( strtol(argv[++i],end,10) );
                assert(mosquitoMultipliers.size() > 0);
                int running_sum = 0;
                for (unsigned int j=0; j<mosquitoMultipliers.size(); j++) {
                    mosquitoMultipliers[j].start = running_sum;
                    mosquitoMultipliers[j].duration = strtol(argv[++i],end,10);
                    mosquitoMultipliers[j].value = strtod(argv[++i],end);
                    running_sum += mosquitoMultipliers[j].duration;
                }
            }
            else if (strcmp(argv[i], "-extrinsicincubations")==0) {
                extrinsicIncubationPeriods.clear();
                extrinsicIncubationPeriods.resize( strtol(argv[++i],end,10) );
                assert(extrinsicIncubationPeriods.size() > 0);
                int running_sum = 0;
                for (unsigned int j=0; j<extrinsicIncubationPeriods.size(); j++) {
                    extrinsicIncubationPeriods[j].start = running_sum;
                    extrinsicIncubationPeriods[j].duration = strtol(argv[++i],end,10);
                    extrinsicIncubationPeriods[j].value = strtod(argv[++i],end);
                    running_sum += extrinsicIncubationPeriods[j].duration;
                }
            }
            else if (strcmp(argv[i], "-daysimmune")==0) {
                nDaysImmune = strtol(argv[++i],end,10);
            }
            else if (strcmp(argv[i], "-startdayofyear")==0) {
                startDayOfYear = strtol(argv[++i],end,10);
            }
            else if (strcmp(argv[i], "-VES")==0 || strcmp(argv[i], "-ves")==0) {
                fVESs.clear();
                fVESs.resize(NUM_OF_SEROTYPES, strtod(argv[++i],end));
            }
            else if (strcmp(argv[i], "-VESs")==0 || strcmp(argv[i], "-vess")==0) {
                fVESs.clear();
                fVESs.resize(NUM_OF_SEROTYPES);
                // different VES for each serotype
                for (int j=0; j<NUM_OF_SEROTYPES; j++) fVESs[j] = strtod(argv[++i],end);
            }
            else if (strcmp(argv[i], "-VESsnaive")==0 || strcmp(argv[i], "-vessnaive")==0) {
                fVESs_NAIVE.clear(); fVESs_NAIVE.resize(NUM_OF_SEROTYPES, 0);
                for (int j=0; j<NUM_OF_SEROTYPES; j++) fVESs_NAIVE[j] = strtod(argv[++i],end);
            }
            else if (strcmp(argv[i], "-VEI")==0 || strcmp(argv[i], "-vei")==0) {
                fVEI = strtod(argv[++i],end);
            }
            else if (strcmp(argv[i], "-VEP")==0 || strcmp(argv[i], "-vep")==0) {
                fVEP = strtod(argv[++i],end);
            }
            else if (strcmp(argv[i], "-vaccineleaky")==0) {
                bVaccineLeaky=true;
            }
            else if (strcmp(argv[i], "-retroactivematurevaccine")==0) {
                bRetroactiveMatureVaccine=true;
            }
            else if (strcmp(argv[i], "-nosecondary")==0) {
                bSecondaryTransmission = false;
            }
            else if (strcmp(argv[i], "-popfile")==0) {
                populationFilename = argv[++i];
            }
            else if (strcmp(argv[i], "-immfile")==0) {
                immunityFilename = argv[++i];
            }
            else if (strcmp(argv[i], "-locfile")==0) {
                locationFilename = argv[++i];
            }
            else if (strcmp(argv[i], "-netfile")==0) {
                networkFilename = argv[++i];
            }
            else if (strcmp(argv[i], "-peopleoutputfile")==0) {
                peopleOutputFilename = argv[++i];
            }
            else if (strcmp(argv[i], "-yearlypeopleoutputfile")==0) {
                yearlyPeopleOutputFilename = argv[++i];
            }
            else if (strcmp(argv[i], "-dailyoutputfile")==0) {
                dailyOutputFilename = argv[++i];
            }
            else if (strcmp(argv[i], "-probfile")==0) {
                swapProbFilename = argv[++i];
            }
            else if (strcmp(argv[i], "-annualintrosfile")==0) {
                annualIntroductionsFilename = argv[++i];
                loadAnnualIntroductions(annualIntroductionsFilename);
            }
            else if (strcmp(argv[i], "-annualserotypefile")==0) {
                annualSerotypeFilename = argv[++i];
                loadAnnualSerotypes(annualSerotypeFilename);
            }
            else if (strcmp(argv[i], "-simulateannualserotypes")==0) {
                simulateAnnualSerotypes = true;
            }
            else if (strcmp(argv[i], "-normalizeserotypeintros")==0) {
                normalizeSerotypeIntros = true;
            }
            else if (strcmp(argv[i], "-dailyeipfile")==0) {
                dailyEIPfilename = argv[++i];
                loadDailyEIP(dailyEIPfilename);
            }
            else if (strcmp(argv[i], "-dailyoutput")==0) {
                dailyOutput = true;
            }
            else if (strcmp(argv[i], "-weeklyoutput")==0) {
                weeklyOutput = true;
            }
            else if (strcmp(argv[i], "-monthlyoutput")==0) {
                monthlyOutput = true;
            }
            else if (strcmp(argv[i], "-yearlyoutput")==0) {
                yearlyOutput = true;
            }
            else if (strcmp(argv[i], "-abcverbose")==0) {
                abcVerbose = true;
            }
            else {
                cerr << "Unknown option: " << argv[i] << endl;
                cerr << "Check arguments and formatting." << endl;
                exit(-1);
            }
        }
    }

    gsl_rng_set(RNG, randomseed);
    // runlength and randomseed need to be set before calling generateAnnualSerotypes()
    if (simulateAnnualSerotypes) generateAnnualSerotypes();
    validate_parameters();
}

void Parameters::validate_parameters() {
    cerr << "population file = " << populationFilename << endl;
    cerr << "immunity file = " << immunityFilename << endl;
    cerr << "location file = " << locationFilename << endl;
    cerr << "network file = " << networkFilename << endl;
    cerr << "swap probabilities file = " << swapProbFilename << endl;
    cerr << "runlength = " << nRunLength << endl;
    cerr << "start day of year (1 is Jan 1st) = " << startDayOfYear << endl;
    cerr << "random seed = " << randomseed << endl;
    cerr << "beta_PM = " << betaPM << endl;
    cerr << "beta_MP = " << betaMP << endl;
    cerr << "days of complete cross protection = " << nDaysImmune << endl;
    cerr << "mosquito move prob = " << fMosquitoMove << endl;
    cerr << "mosquito move model = " << mosquitoMoveModel << endl;
    if ( mosquitoMoveModel != "uniform" and mosquitoMoveModel != "weighted" ) {
        cerr << "ERROR: invalid mosquito movement model requested:" << endl;
        cerr << " -mosquitomovemodel may be uniform or weighted " << nRunLength << endl;
        exit(-1);
    }
    cerr << "mosquito teleport prob = " << fMosquitoTeleport << endl;
    cerr << "default mosquito capacity per building = " << nDefaultMosquitoCapacity << endl;
    if (annualSerotypeFilename == "") {
        cerr << "number of daily exposures / serotype weights =";
        for (int i=0; i<NUM_OF_SEROTYPES; i++) {
            cerr << " " << nDailyExposed[0][i];
        }
        cerr << endl;
    } else {
        if (simulateAnnualSerotypes) {
            cerr << "ERROR: -simulateannualserotypes and -annualserotypefile cannot both be used" << endl;
        } else {
            cerr << "annual serotype file = " << annualSerotypeFilename << endl;
        }
    }
    if (simulateAnnualSerotypes) {
        cerr << "simulating annual serotypes" << endl;
    }
    if (dailyEIPfilename != "") {
        cerr << "daily EIP file = " << dailyEIPfilename << endl;
    }
    if (eMosquitoDistribution==CONSTANT) {
        cerr << "mosquito capacity distribution is constant" << endl;
    } else if (eMosquitoDistribution==EXPONENTIAL)
        cerr << "mosquito capacity distribution is exponential" << endl;
        if (mosquitoMultipliers.size()>0) {
            cerr << "mosquito seasonal multipliers (days,mult) =";
            for (unsigned int j=0; j<mosquitoMultipliers.size(); j++) {
                cerr << " (" << mosquitoMultipliers[j].duration << "," << mosquitoMultipliers[j].value << ")";
            }
            cerr << endl;
        }
    if (extrinsicIncubationPeriods.size()>0) {
        cerr << "extrinsic incubation periods (days,EIP) =";
        for (unsigned int j=0; j<extrinsicIncubationPeriods.size(); j++) {
            cerr << " (" << extrinsicIncubationPeriods[j].duration << "," << extrinsicIncubationPeriods[j].value << ")";
            if (j >= 12) {
                cerr << " . . . " << endl << "\t" << extrinsicIncubationPeriods.size() - j << " more extrinsic incubation periods not displayed." << endl;
                break;
            }
        }
        cerr << endl;
    }
    cerr << "VE_I = " << fVEI << endl;
    cerr << "VE_P = " << fVEP << endl;
    if (fVEI>1.0 || fVEP>1.0) {
        cerr << "ERROR: VE_S, VE_I, and VE_P must be between 0 and 1" << endl;
        exit(-1);
    }
    cerr << "VE_Ss = " << fVESs[0] << "," << fVESs[1] << "," << fVESs[2] << "," << fVESs[3] << endl;

    if (fVESs_NAIVE.size()==0) {
      fVESs_NAIVE.clear();
      fVESs_NAIVE = fVESs; // naive people have the same VE_S as non-naive
      cerr << "Vaccine protection is independent of infection history (naive == non-naive)" << endl;
    } else {
      cerr << "VE_Ss naive = " << fVESs_NAIVE[0] << "," << fVESs_NAIVE[1] << "," << fVESs_NAIVE[2] << "," << fVESs_NAIVE[3] << endl;
    }

    if (bVaccineLeaky) {
        cerr << "VE_S is leaky" << endl;
    } else {
        cerr << "VE_S is all-or-none" << endl;
    }

    if (bRetroactiveMatureVaccine) {
        if (bVaccineLeaky) {
            cerr << "Vaccine protection is upgraded from naive to non-naive upon infection (retroactive maturity)" << endl;
        } else {
            cerr << "-retroactivematurevaccine is only defined for leaky vaccines (-vaccineleaky)" << endl;
            exit(-1);
        }
    }

    if (peopleOutputFilename.length()>0) {
        cerr << "people output file = " << peopleOutputFilename << endl;
    } else {
        cerr << "no people output file" << endl;
    }
    if (annualIntroductionsFilename.length()>0) {
        cerr << "annual introductions file = " << annualIntroductionsFilename << endl;
    }
    if (yearlyPeopleOutputFilename.length()>0) {
        cerr << "yearly people output file = " << yearlyPeopleOutputFilename << endl;
    }
    if (dailyOutputFilename.length()>0) {
        cerr << "daily output file = " << dailyOutputFilename << endl;
    } else {
        cerr << "no daily output file" << endl;
    }
}


void Parameters::loadAnnualIntroductions(string annualIntrosFilename) {
    ifstream iss(annualIntrosFilename.c_str());
    if (!iss) {
        cerr << "ERROR: " << annualIntrosFilename << " not found." << endl;
        exit(114);
    }
    annualIntroductions.clear();

    char buffer[500];
    double intros;
    istringstream line(buffer);

    while (iss) {
        iss.getline(buffer,500);
        line.clear();
        line.str(buffer);
        if (line >> intros) {
            annualIntroductions.push_back(intros);
        }
    }

    iss.close();
    return;
}


void Parameters::loadAnnualSerotypes(string annualSerotypeFilename) {
    ifstream iss(annualSerotypeFilename.c_str());
    if (!iss) {
        cerr << "ERROR: " << annualSerotypeFilename << " not found." << endl;
        exit(115);
    }

    // get rid of anything there now
    for (auto &v: nDailyExposed) v.clear();
    nDailyExposed.clear();

    char sep = ' ';
    string line;

    while ( getline(iss,line) ) {
        // expecting four serotype intro rates per line
        // each line corresponds to one year
        vector<string> fields = dengue::util::split(line, sep);

        if (fields.size() == NUM_OF_SEROTYPES) {
            vector<float>row(fields.size());
            for( unsigned int i=0; i < fields.size(); i++ ) {
                row[i] = dengue::util::string2double(fields[i]);
            }

            nDailyExposed.push_back(row);
        } else {
            cerr << "WARNING: Found line with unexpected number of values in annual serotypes file" << endl;
            cerr << "\tFile: " << annualSerotypeFilename << endl;
            cerr << "\tLine: " << line << endl;
        }
    }

    return;
}


void Parameters::writeAnnualSerotypes(string filename) const {
    char sep = ' ';

    ofstream file;
    file.open(filename);
    for (auto &year: nDailyExposed) {
        for (auto val: year) file << val << sep;
        file << endl;
    }

    file.close();
}


void Parameters::generateAnnualSerotypes(int total_num_years) { // default arg value is -1
    enum State {GAP, RUN};
    // get rid of anything there now
    for (auto &v: nDailyExposed) v.clear();
    nDailyExposed.clear();

    if (total_num_years == -1) {
        total_num_years = (int) ceil((double) nRunLength / 365.0);
    }
    nDailyExposed.resize(total_num_years, vector<float>(NUM_OF_SEROTYPES));

    for (int s = 0; s<NUM_OF_SEROTYPES; ++s) {
        const float p_gap = 1.0/MEAN_GAP_LENGTH[s];
        const float p_run = 1.0/MEAN_RUN_LENGTH[s];
        const float p_gap_start = MEAN_GAP_LENGTH[s]/(MEAN_GAP_LENGTH[s] + MEAN_RUN_LENGTH[s]);

        int years_so_far = 0;
        State state = p_gap_start > gsl_rng_uniform(RNG) ? GAP : RUN;
        while (years_so_far < total_num_years) {
            if (state == GAP) {
                unsigned int gap = gsl_ran_geometric(RNG, p_gap);
                // We may not be starting at the beginning of a gap
                // Also: gsl_rng_uniform_int() returns ints on [0,n-1]
                if (years_so_far == 0) gap = gsl_rng_uniform_int(RNG, gap) + 1;
                while (gap > 0 and (unsigned) years_so_far < nDailyExposed.size()) {
                    nDailyExposed[years_so_far++][s] = 0.0;
                    --gap;
                }
                state = RUN;
            } else {
                unsigned int run = gsl_ran_geometric(RNG, p_run);
                if (years_so_far == 0) run = gsl_rng_uniform_int(RNG, run) + 1;
                while (run > 0 and (unsigned) years_so_far < nDailyExposed.size()) {
                    nDailyExposed[years_so_far++][s] = 1.0;
                    --run;
                }
                state = GAP;
            }
        }
    }
    if (not abcVerbose) {
        cerr << "Serotype runs:" << endl;
        for (auto &y: nDailyExposed) {
            for (auto v: y)  cerr << v << " "; cerr << endl;
        }
    }
    if (normalizeSerotypeIntros) {
        for (unsigned int i = 0; i < nDailyExposed.size(); ++i) {
            float total = accumulate(nDailyExposed[i].begin(), nDailyExposed[i].end(), 0.0F);
            if (total > 0) {
                for (unsigned int s = 0; s < nDailyExposed[i].size(); ++s) nDailyExposed[i][s] /= total;
            }
        }
    }

    if (not abcVerbose) {
        cerr << "Serotype runs (normalized):" << endl;
        for (auto &y: nDailyExposed) {
            for (auto v: y)  cerr << v << " "; cerr << endl;
        }
    }
    return;
}


void Parameters::loadDailyEIP(string dailyEIPfilename, int desired_size) { // default desired size == 0
    // expecting first value on each line to be the EIP mu
    // other, subsequent values are permitted, but ignored
    vector<string> first_column = dengue::util::read_vector_file(dailyEIPfilename);

    extrinsicIncubationPeriods.clear();

    const int duration = 1;
    int start = 0;
    double value = 0;
    for(string val_str: first_column) {
        value = dengue::util::to_double(val_str);
        if (value <= 0) {
            cerr << "An EIP <= 0 was read from " << dailyEIPfilename << "." << endl;
            cerr << "Value read: " << value << endl;
            cerr << "This is nonsensical and indicates a non-numerical value in the first column or an actual bad value." << endl;
            exit(113);
        }
        extrinsicIncubationPeriods.emplace_back(start, duration, value);
        start += duration;
    }

    // this allows us to repeat a sequence many times and then modify it elsewhere,
    // e.g. for vector control modeling
    const int orig_len = extrinsicIncubationPeriods.size();
    while ((signed) extrinsicIncubationPeriods.size() < desired_size) {
        for (int i = 0; i < orig_len; ++i) {
            auto value = extrinsicIncubationPeriods[i].value;
            extrinsicIncubationPeriods.emplace_back(start, duration, value);
            start += duration;
        }
    }

    if (desired_size > 0 and (signed) extrinsicIncubationPeriods.size() > desired_size) extrinsicIncubationPeriods.resize(desired_size);

    return;
}


void Parameters::loadDailyMosquitoMultipliers(string mosquitoMultiplierFilename, int desired_size) { // default desired size == 0
    // expecting first value on each line to be the mosquito multiplier (floats on [0,1])
    // other, subsequent values are permitted, but ignored
    vector<string> first_column = dengue::util::read_vector_file(mosquitoMultiplierFilename);

    mosquitoMultipliers.clear();

    const int duration = 1;  // this assumption is what makes it daily
    int start = 0;
    double value = 0;
    for(string val_str: first_column) {

        value = dengue::util::to_double(val_str);
        if (value < 0.0 or value > 1.0) {
            cerr << "A mosquito multiplier < 0 or > 1 was read from " << mosquitoMultiplierFilename << "." << endl;
            cerr << "Value read: " << value << endl;
            cerr << "This is nonsensical and indicates a non-numerical value in the first column or an actual bad value." << endl;
            exit(119);
        }
        mosquitoMultipliers.emplace_back(start, duration, value);
        start += duration;
    }

    // this allows us to repeat a sequence many times and then modify it elsewhere,
    // e.g. for vector control modeling
    const int orig_len = mosquitoMultipliers.size();
    while ((signed) mosquitoMultipliers.size() < desired_size) {
        for (int i = 0; i < orig_len; ++i) {
            auto value = mosquitoMultipliers[i].value;
            mosquitoMultipliers.emplace_back(start, duration, value);
            start += duration;
        }
    }

    if (desired_size > 0 and (signed) mosquitoMultipliers.size() > desired_size) mosquitoMultipliers.resize(desired_size);

    return;
}


void Parameters::defineSerotypeRelativeRisks() { // should be called after reportedFractions (1/expansion factors) are set, if they're going to be
    // values fitted in
    // Reich et al, Interactions between serotypes of dengue highlight epidemiological impact of cross-immunity, Interface, 2013
    // Normalized from Fc values in supplement table 2, available at
    // http://rsif.royalsocietypublishing.org/content/10/86/20130414/suppl/DC1
    vector<double> reported_by_sero = {1.20, 0.99, 1.00, 0.38}; // Reich et al Table 2 Fc medians; total cases per 1,000 infections
    vector<double> df_by_sero = {353.0/1223.0, 176.0/1477.0, 397.0/1542.0, 95.0/480.0}; // Nisalak et al http://www.ajtmh.org/content/68/2/191.full.pdf+html Table 4, p. 199
    vector<double> dhf_by_sero;
    for (double val: df_by_sero) dhf_by_sero.push_back(1.0 - val);

    vector<double> implied_cases_by_sero;
    for (unsigned int i = 0; i < reported_by_sero.size(); ++i) {
        implied_cases_by_sero.push_back(reported_by_sero[i]*(df_by_sero[i]/reportedFraction[(int) MILD] + dhf_by_sero[i]/reportedFraction[(int) SEVERE]));
    }

    assert(implied_cases_by_sero[0] > 0);
    vector<double> relative_risk_by_sero = implied_cases_by_sero;
    for (double& val: relative_risk_by_sero) val /= implied_cases_by_sero[0];
    serotypePathogenicityRelativeRisks = relative_risk_by_sero;
}
