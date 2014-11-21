#include "Parameters.h"
#include "Location.h"
#include "Utility.h"
#include <fstream>
#include <sstream>
#include <numeric>
#include <gsl/gsl_randist.h>

void Parameters::define_defaults() {
    randomseed = 5489;
    nRunLength = 100;
    betaPM = 0.2;
    betaMP = 0.1;
    fMosquitoMove = 0.2;
    mosquitoMoveModel = "weighted";
    fMosquitoTeleport = 0.01;
    fVEI = 0.0;
    fVEP = 0.0;
    fVESs.clear(); fVESs.resize(NUM_OF_SEROTYPES, 0.95);
    fVESs_NAIVE.clear();
    bVaccineLeaky = false;
    fPreVaccinateFraction = 0.0;
    nDefaultMosquitoCapacity = 20;                      // mosquitoes per location
    eMosquitoDistribution = CONSTANT;
    bSecondaryTransmission = true;
    populationFilename = "population-64.txt";
    immunityFilename = "";
    networkFilename = "locations-network-64.txt";
    locationFilename = "locations-64.txt";
    peopleOutputFilename = "";
    yearlyPeopleOutputFilename = "";
    dailyOutputFilename = "";
    swapProbFilename = "";
    annualIntroductionsFilename = "";  // time series of some external factor determining introduction rate
    annualIntroductionsCoef = 1;     // multiplier to rescale external introductions to something sensible
    normalizeSerotypeIntros = false;
    annualIntroductions.clear();
    annualIntroductions.push_back(1.0);
    nDaysImmune = 365;
    nSizeVaccinate = 0;                                 // number of parts in phased vaccination
    nSizePrevaccinateAge = 0;
    nMaxInfectionParity = NUM_OF_SEROTYPES;
    expansionFactor = 1;
    nDailyExposed.push_back(vector<float>(NUM_OF_SEROTYPES, 0.0)); // default is no introductions
    annualSerotypeFilename = "";
    dailyEIPfilename = "";

    nInitialExposed.clear();
    nInitialExposed.resize(NUM_OF_SEROTYPES, 0);
    nInitialInfected.clear();
    nInitialInfected.resize(NUM_OF_SEROTYPES, 0);

    fPrimaryPathogenicity.clear();
    fPrimaryPathogenicity.resize(NUM_OF_SEROTYPES, 1.0);
    if (NUM_OF_SEROTYPES == 4) {
        // values fitted in
        // Reich et al, Interactions between serotypes of dengue highlight epidemiological impact of cross-immunity, Interface, 2013
        // Normalized from Fc values in supplement table 2, available at
        // http://rsif.royalsocietypublishing.org/content/10/86/20130414/suppl/DC1
        fPrimaryPathogenicity[0] = 1.000;
        fPrimaryPathogenicity[1] = 0.825;
        fPrimaryPathogenicity[2] = 0.833;
        fPrimaryPathogenicity[3] = 0.317;
    }

    fSecondaryScaling.clear();
    fSecondaryScaling.resize(NUM_OF_SEROTYPES, 1.0);

    nVaccinateYear.clear();
    nVaccinateAge.clear();
    fVaccinateFraction.clear();

    const vector<float> MOSQUITO_MULTIPLIER_DEFAULTS = {0.179,0.128,0.123,0.0956,0.195,0.777,0.940,0.901,1.0,0.491,0.301,0.199};
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

void Parameters::readParameters(int argc, char *argv[]) {
    cerr << "Dengue model, Version " << VERSION_NUMBER_MAJOR << "." << VERSION_NUMBER_MINOR << endl;
    cerr << "written by Dennis Chao and Thomas Hladish in 2012-2014" << endl;

    if (argc>1) {
        for (int i=1; i<argc; i++) {
            char **end = NULL;
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
            else if (strcmp(argv[i], "-primarypathogenicity")==0) {
                for (int j=0; j<NUM_OF_SEROTYPES; j++) fPrimaryPathogenicity[j]=strtod(argv[++i],end);
            }
            else if (strcmp(argv[i], "-secondaryscaling")==0) {
                for (int j=0; j<NUM_OF_SEROTYPES; j++) fSecondaryScaling[j]=strtod(argv[++i],end);
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
            else if (strcmp(argv[i], "-expansionfactor")==0) {
                expansionFactor = strtod(argv[++i],end);
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
            else if (strcmp(argv[i], "-vaccinatephased")==0) {
                nSizeVaccinate = strtol(argv[++i],end,10);
                for (int j=0; j<nSizeVaccinate; j++) {
                    nVaccinateYear.push_back( strtol(argv[++i],end,10) );
                    nVaccinateAge.push_back( strtol(argv[++i],end,10) );
                    fVaccinateFraction.push_back( strtod(argv[++i],end) );
                }
            }
            else if (strcmp(argv[i], "-daysimmune")==0) {
                nDaysImmune = strtol(argv[++i],end,10);
            }
            else if (strcmp(argv[i], "-startdayofyear")==0) {
                startDayOfYear = strtol(argv[++i],end,10);
            }
            else if (strcmp(argv[i], "-maxinfectionparity")==0) {
                nMaxInfectionParity = strtol(argv[++i],end,10);
                assert(nMaxInfectionParity>0 && nMaxInfectionParity<=NUM_OF_SEROTYPES);
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
            else if (strcmp(argv[i], "-prevaccinate")==0) {
                fPreVaccinateFraction = strtod(argv[++i],end);
            }
            else if (strcmp(argv[i], "-prevaccinateage")==0) {
                nSizePrevaccinateAge = strtol(argv[++i],end,10);
                for (int j=0; j<nSizePrevaccinateAge; j++) {
                    nPrevaccinateAgeMin[j] = strtol(argv[++i],end,10);
                    nPrevaccinateAgeMax[j] = strtol(argv[++i],end,10);
                    fPrevaccinateAgeFraction[j] = strtod(argv[++i],end);
                }
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
    cerr << "maximum infection parity = " << nMaxInfectionParity << endl;
    cerr << "pathogenicity of primary infection =";
    for (int i=0; i<NUM_OF_SEROTYPES; i++) {
        cerr << " " << fPrimaryPathogenicity[i];
    }
    cerr << endl;
    cerr << "pathogenicity scaling for secondary infection =";
    for (int i=0; i<NUM_OF_SEROTYPES; i++) {
        cerr << " " << fSecondaryScaling[i];
    }
    cerr << endl;
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
    cerr << "Pre-vaccinate fraction = " << fPreVaccinateFraction << endl;
    if (nSizeVaccinate>0) {
        cerr << "Phased vaccinate (year, age, frac) = ";
        for (int j=0; j<nSizeVaccinate; j++) {
            cerr << " (" << nVaccinateYear[j] << "," << nVaccinateAge[j]  << "," << fVaccinateFraction[j] << ")";
        }
        cerr << endl;
    }
    if (nSizePrevaccinateAge>0) {
        cerr << "Pre-vaccinate by age (min, max, frac) = ";
        for (int j=0; j<nSizePrevaccinateAge; j++) {
            cerr << " (" << nPrevaccinateAgeMin[j] << "," << nPrevaccinateAgeMax[j]  << "," << fPrevaccinateAgeFraction[j] << ")";
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
    for (auto v: nDailyExposed) v.clear();
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


void Parameters::generateAnnualSerotypes() {
    enum State {GAP, RUN};
    // get rid of anything there now
    for (auto v: nDailyExposed) v.clear();
    nDailyExposed.clear();

    int run_length_years = (int) ceil((double) nRunLength / 365.0);
    nDailyExposed.resize(run_length_years, vector<float>(NUM_OF_SEROTYPES));

    const float p_gap = 1.0/MEAN_GAP_LENGTH;
    const float p_run = 1.0/MEAN_RUN_LENGTH;
    const float p_gap_start = MEAN_GAP_LENGTH/(MEAN_GAP_LENGTH + MEAN_RUN_LENGTH);

    for (int s = 0; s<NUM_OF_SEROTYPES; ++s) {
        int years_so_far = 0;
        State state = p_gap_start > gsl_rng_uniform(RNG) ? GAP : RUN;
        while (years_so_far < run_length_years) {
            if (state == GAP) {
                unsigned int gap = gsl_ran_geometric(RNG, p_gap);
                while (gap > 0 and (unsigned) years_so_far < nDailyExposed.size()) {
                    nDailyExposed[years_so_far++][s] = 0.0;
                    --gap;
                }
                state = RUN;
            } else {
                unsigned int run = gsl_ran_geometric(RNG, p_run);
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
        for (auto y: nDailyExposed) {
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
        for (auto y: nDailyExposed) {
            for (auto v: y)  cerr << v << " "; cerr << endl;
        }
    }
    return;
}


void Parameters::loadDailyEIP(string dailyEIPfilename) {
    ifstream iss(dailyEIPfilename.c_str());
    if (!iss) {
        cerr << "ERROR: " << dailyEIPfilename << " not found." << endl;
        exit(116);
    }

    // get rid of anything there now
    extrinsicIncubationPeriods.clear();

    char **end = NULL;
    char sep = ' ';
    string line;

    const int duration = 1;
    int start = 0;
    double value = 0;
    while ( getline(iss,line) ) {
        // expecting first value on each line to be the EIP for any mosquito infected on that day
        // other, subsequent values are permitted, but ignored
        vector<string> fields = dengue::util::split(line, sep);

        if (fields.size() >= 1) {
            value = strtod(fields[0].c_str(), end);
            if (value <= 0) {
                cerr << "An EIP <= 0 was read from " << dailyEIPfilename << "." << endl;
                cerr << "This is nonsensical and indicates a non-numerical value in the first column or an actual bad value." << endl;
                exit(113);
            }
            extrinsicIncubationPeriods.emplace_back(start, duration, value);
            start += duration;
        } else {
            cerr << "WARNING: Found line with unexpected number of values in daily EIP file" << endl;
            cerr << "\tFile: " << dailyEIPfilename << endl;
            cerr << "\tLine: " << line << endl;
        }
    }

    return;
}

