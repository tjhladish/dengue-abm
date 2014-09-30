#include "Parameters.h"
#include "Location.h"
#include "Utility.h"
#include <fstream>
#include <sstream>

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
    annualIntroductions.clear();
    annualIntroductions.push_back(1.0);
    nDaysImmune = 365;
    nSizeVaccinate = 0;                                 // number of parts in phased vaccination
    nSizePrevaccinateAge = 0;
    nMaxInfectionParity = NUM_OF_SEROTYPES;
    expansionFactor = 1;
    nDailyExposed.push_back(std::vector<float>(NUM_OF_SEROTYPES, 0.0)); // default is no introductions
    annualSerotypeFilename = "";
    dailyEIPfilename = "";

    for (int i=0; i<NUM_OF_SEROTYPES; i++) {
        nInitialExposed[i]=0;
        nInitialInfected[i]=0;
    }
    fPrimaryPathogenicity.clear();
    fPrimaryPathogenicity.resize(NUM_OF_SEROTYPES, 1.0);
    fPrimaryPathogenicity[1] = fPrimaryPathogenicity[3] = 0.25;

    fSecondaryScaling.clear();
    fSecondaryScaling.resize(NUM_OF_SEROTYPES, 1.0);

    nVaccinateYear.clear();
    nVaccinateAge.clear();
    fVaccinateFraction.clear();
}

void Parameters::readParameters(int argc, char *argv[]) {
    std::cerr << "Dengue model, Version " << VERSION_NUMBER_MAJOR << "." << VERSION_NUMBER_MINOR << std::endl;
    std::cerr << "written by Dennis Chao and Thomas Hladish in 2012-2014" << std::endl;

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
                    std::cerr << "WARNING: Annual serotype file specified.  Ignoring daily exposed parameter.\n"; 
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
                    std::cerr << "ERROR: Invalid mosquito distribution specified." << std::endl;
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
            else if (strcmp(argv[i], "-dailyeipfile")==0) {
                dailyEIPfilename = argv[++i];
                loadDailyEIP(dailyEIPfilename);
            }
            else {
                std::cerr << "Unknown option: " << argv[i] << std::endl;
                std::cerr << "Check arguments and formatting." << std::endl;
                exit(-1);
            }
        }
    }
    validate_parameters();
}

void Parameters::validate_parameters() {
    std::cerr << "population file = " << populationFilename << std::endl;
    std::cerr << "immunity file = " << immunityFilename << std::endl;
    std::cerr << "location file = " << locationFilename << std::endl;
    std::cerr << "network file = " << networkFilename << std::endl;
    std::cerr << "swap probabilities file = " << swapProbFilename << std::endl;
    std::cerr << "runlength = " << nRunLength << std::endl;
    if (nRunLength>MAX_RUN_TIME) {
        std::cerr << "ERROR: runlength is too long: " << nRunLength << std::endl;
        std::cerr << " change Parameters.h and recompile." << std::endl;
        exit(-1);
    }
    if (nRunLength==365) {
        std::cerr << "WARNING: you probably want runlength to be 364, not 365" << std::endl;
    }
    std::cerr << "random seed = " << randomseed << std::endl;
    std::cerr << "beta_PM = " << betaPM << std::endl;
    std::cerr << "beta_MP = " << betaMP << std::endl;
    std::cerr << "days of complete cross protection = " << nDaysImmune << std::endl;
    std::cerr << "maximum infection parity = " << nMaxInfectionParity << std::endl;
    std::cerr << "pathogenicity of primary infection =";
    for (int i=0; i<NUM_OF_SEROTYPES; i++) {
        std::cerr << " " << fPrimaryPathogenicity[i];
    }
    std::cerr << std::endl;
    std::cerr << "pathogenicity scaling for secondary infection =";
    for (int i=0; i<NUM_OF_SEROTYPES; i++) {
        std::cerr << " " << fSecondaryScaling[i];
    }
    std::cerr << std::endl;
    std::cerr << "mosquito move prob = " << fMosquitoMove << std::endl;
    std::cerr << "mosquito move model = " << mosquitoMoveModel << std::endl;
    if ( mosquitoMoveModel != "uniform" and mosquitoMoveModel != "weighted" ) {
        std::cerr << "ERROR: invalid mosquito movement model requested:" << std::endl;
        std::cerr << " -mosquitomovemodel may be uniform or weighted " << nRunLength << std::endl;
        exit(-1);
    }
    std::cerr << "mosquito teleport prob = " << fMosquitoTeleport << std::endl;
    std::cerr << "default mosquito capacity per building = " << nDefaultMosquitoCapacity << std::endl;
    if (annualSerotypeFilename == "") {
        std::cerr << "number of daily exposures / serotype weights =";
        for (int i=0; i<NUM_OF_SEROTYPES; i++) {
            std::cerr << " " << nDailyExposed[0][i];
        }
        std::cerr << std::endl;
    } else {
        std::cerr << "annual serotype file = " << annualSerotypeFilename << std::endl;
    }
    if (dailyEIPfilename != "") {
        std::cerr << "daily EIP file = " << dailyEIPfilename << std::endl;
    }
    if (eMosquitoDistribution==CONSTANT) {
        std::cerr << "mosquito capacity distribution is constant" << std::endl;
    } else if (eMosquitoDistribution==EXPONENTIAL)
        std::cerr << "mosquito capacity distribution is exponential" << std::endl;
        if (mosquitoMultipliers.size()>0) {
            std::cerr << "mosquito seasonal multipliers (days,mult) =";
            for (unsigned int j=0; j<mosquitoMultipliers.size(); j++) {
                std::cerr << " (" << mosquitoMultipliers[j].duration << "," << mosquitoMultipliers[j].value << ")";
            }
            std::cerr << std::endl;
        }
    if (extrinsicIncubationPeriods.size()>0) {
        std::cerr << "extrinsic incubation periods (days,EIP) =";
        for (unsigned int j=0; j<extrinsicIncubationPeriods.size(); j++) {
            std::cerr << " (" << extrinsicIncubationPeriods[j].duration << "," << extrinsicIncubationPeriods[j].value << ")";
        }
        std::cerr << std::endl;
    }
    std::cerr << "Pre-vaccinate fraction = " << fPreVaccinateFraction << std::endl;
    if (nSizeVaccinate>0) {
        std::cerr << "Phased vaccinate (year, age, frac) = ";
        for (int j=0; j<nSizeVaccinate; j++) {
            std::cerr << " (" << nVaccinateYear[j] << "," << nVaccinateAge[j]  << "," << fVaccinateFraction[j] << ")";
        }
        std::cerr << std::endl;
    }
    if (nSizePrevaccinateAge>0) {
        std::cerr << "Pre-vaccinate by age (min, max, frac) = ";
        for (int j=0; j<nSizePrevaccinateAge; j++) {
            std::cerr << " (" << nPrevaccinateAgeMin[j] << "," << nPrevaccinateAgeMax[j]  << "," << fPrevaccinateAgeFraction[j] << ")";
        }
        std::cerr << std::endl;
    }
    std::cerr << "VE_I = " << fVEI << std::endl;
    std::cerr << "VE_P = " << fVEP << std::endl;
    if (fVEI>1.0 || fVEP>1.0) {
        std::cerr << "ERROR: VE_S, VE_I, and VE_P must be between 0 and 1" << std::endl;
        exit(-1);
    }
    std::cerr << "VE_Ss = " << fVESs[0] << "," << fVESs[1] << "," << fVESs[2] << "," << fVESs[3] << std::endl;

    if (fVESs_NAIVE.size()==0) {
      fVESs_NAIVE.clear();
      fVESs_NAIVE = fVESs; // naive people have the same VE_S as non-naive
      std::cerr << "Vaccine protection is independent of infection history (naive == non-naive)" << std::endl;
    } else {
      std::cerr << "VE_Ss naive = " << fVESs_NAIVE[0] << "," << fVESs_NAIVE[1] << "," << fVESs_NAIVE[2] << "," << fVESs_NAIVE[3] << std::endl;
    }

    if (bVaccineLeaky) {
        std::cerr << "VE_S is leaky" << std::endl;
    } else {
        std::cerr << "VE_S is all-or-none" << std::endl;
    }

    if (bRetroactiveMatureVaccine) {
        if (bVaccineLeaky) {
            std::cerr << "Vaccine protection is upgraded from naive to non-naive upon infection (retroactive maturity)" << std::endl;
        } else {
            std::cerr << "-retroactivematurevaccine is only defined for leaky vaccines (-vaccineleaky)" << std::endl;
            exit(-1); 
        } 
    }

    if (peopleOutputFilename.length()>0) {
        std::cerr << "people output file = " << peopleOutputFilename << std::endl;
    } else {
        std::cerr << "no people output file" << std::endl;
    }
    if (annualIntroductionsFilename.length()>0) {
        std::cerr << "annual introductions file = " << annualIntroductionsFilename << std::endl;
    }
    if (yearlyPeopleOutputFilename.length()>0) {
        std::cerr << "yearly people output file = " << yearlyPeopleOutputFilename << std::endl;
    }
    if (dailyOutputFilename.length()>0) {
        std::cerr << "daily output file = " << dailyOutputFilename << std::endl;
    } else {
        std::cerr << "no daily output file" << std::endl;
    }
}


void Parameters::loadAnnualIntroductions(std::string annualIntrosFilename) {
    std::ifstream iss(annualIntrosFilename.c_str());
    if (!iss) {
        std::cerr << "ERROR: " << annualIntrosFilename << " not found." << std::endl;
        exit(114);
    }
    annualIntroductions.clear();
 
    char buffer[500];
    double intros;
    std::istringstream line(buffer);

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


void Parameters::loadAnnualSerotypes(std::string annualSerotypeFilename) {
    std::ifstream iss(annualSerotypeFilename.c_str());
    if (!iss) {
        std::cerr << "ERROR: " << annualSerotypeFilename << " not found." << std::endl;
        exit(115);
    }

    // get rid of anything there now
    for (auto v: nDailyExposed) v.clear();
    nDailyExposed.clear();

    char sep = ' ';
    std::string line;

    while ( getline(iss,line) ) {
        // expecting four serotype intro rates per line
        // each line correspons to one year
        std::vector<std::string> fields = dengue::util::split(line, sep);

        if (fields.size() == NUM_OF_SEROTYPES) {
            std::vector<float>row(fields.size());
            for( unsigned int i=0; i < fields.size(); i++ ) {
                row[i] = dengue::util::string2double(fields[i]);
            }

            nDailyExposed.push_back(row);
        } else {
            std::cerr << "WARNING: Found line with unexpected number of values in annual serotypes file" << std::endl;
            std::cerr << "\tFile: " << annualSerotypeFilename << std::endl;
            std::cerr << "\tLine: " << line << std::endl;
        }
    }

    return;
}


void Parameters::loadDailyEIP(std::string dailyEIPfilename) {
    std::ifstream iss(dailyEIPfilename.c_str());
    if (!iss) {
        std::cerr << "ERROR: " << dailyEIPfilename << " not found." << std::endl;
        exit(116);
    }

    // get rid of anything there now
    extrinsicIncubationPeriods.clear();

    char **end = NULL;
    char sep = ' ';
    std::string line;

    const int duration = 1;
    int start = 0;
    double value = 0;
    while ( getline(iss,line) ) {
        // expecting first value on each line to be the EIP for any mosquito infected on that day
        // other, subsequent values are permitted, but ignored
        std::vector<std::string> fields = dengue::util::split(line, sep);

        if (fields.size() >= 1) {
            value = strtod(fields[0].c_str(), end);
            if (value <= 0) {
                std::cerr << "An EIP <= 0 was read from " << dailyEIPfilename << "." << std::endl;
                std::cerr << "This is nonsensical and indicates a non-numerical value in the first column or an actual bad value." << std::endl;
                exit(113);
            }
            extrinsicIncubationPeriods.emplace_back(start, duration, value);
            start += duration;
        } else {
            std::cerr << "WARNING: Found line with unexpected number of values in daily EIP file" << std::endl;
            std::cerr << "\tFile: " << dailyEIPfilename << std::endl;
            std::cerr << "\tLine: " << line << std::endl;
        }
    }

    return;
}

