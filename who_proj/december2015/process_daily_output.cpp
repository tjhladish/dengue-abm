#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <unordered_set>
#include <string>
#include <algorithm>
#include <iterator>
#include <assert.h>
#include <glob.h>
#include <sys/stat.h>
#include "sqdb.h"

/*
asymtomatic, mild, severe (3)
by
age (100)
by
year (30)
by
serotype (4)
by

vaccination(2)
by
mechanism (2)
by
catchup (2)
by
intensity (5)
*/


using namespace std;

typedef map< string, vector< vector< vector< vector<double> > > > > AgType;
const int MAX_AGE = 100;
const int BURNIN  = 20;
const int INTERVENTION_DURATION = 30;

inline vector<string> glob(const string& pat){
    glob_t glob_result;
    glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
    vector<string> ret;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        ret.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return ret;
}


string slurp(ifstream& in) {
    return static_cast<stringstream const&>(stringstream() << in.rdbuf()).str();
}


bool fileExists(const std::string& filename) {
    struct stat buf;
    return stat(filename.c_str(), &buf) != -1;
}


inline string extract_checksum_from_filename(const string& s) {
    string::size_type p  = s.find('.');
    string::size_type pp = s.find('.', p + 2); 
    return s.substr(p + 1, pp - p - 1);
    //return stoul(s.substr(p + 1, pp - p - 1));
}

inline string extract_seed_from_filename(const string& s) {
    return s.substr(s.rfind('.') + 1);
    //return stoul(s.substr(p + 1, pp - p - 1));
}

enum SerotypeState {SERO1, SERO2, SERO3, SERO4, NUM_OF_SEROTYPES};
enum SeverityState {ASYMPTOMATIC, MILD, SEVERE, NUM_OF_SEVERITY_TYPES};
//enum HistoryState {PRIMARY, SECONDARY, POST_SECONDARY, NUM_OF_HISTORY_TYPES};

enum VaccineScenario {VACCINE, NO_VACCINE, NUM_OF_VACCINE_TYPES};
enum CatchupScenario {CATCHUP, NO_CATCHUP, NUM_OF_CATCHUP_TYPES};

enum MechanismScenario {BASELINE, ALTERNATE, NUM_OF_MECHANISM_TYPES};
enum IntensityScenario {SP10, SP30, SP50, SP70, SP90, NUM_OF_INTENSITY_TYPES};

//enum MechanismScenario {WANING, NO_WANING, NUM_OF_WANING_TYPES};
//enum BoostingScenario {BOOSTING, NO_BOOSTING, NUM_OF_BOOSTING_TYPES};

struct Scenario {
    Scenario(int v, int c, int w, int b) {
        assert(v < NUM_OF_VACCINE_TYPES);
        assert(c < NUM_OF_CATCHUP_TYPES);
        assert(w < NUM_OF_MECHANISM_TYPES);
        assert(b < NUM_OF_INTENSITY_TYPES);

        vaccine_scenario = (VaccineScenario) v; 
        catchup_scenario = (CatchupScenario) c;
        mechanism_scenario = (MechanismScenario) w;
        intensity_scenario = (IntensityScenario) b;
    }

    string asKey() {
        stringstream ss;
        ss << vaccine_scenario
           << catchup_scenario
           << mechanism_scenario 
           << intensity_scenario;
        return ss.str();
    }

    VaccineScenario vaccine_scenario;
    CatchupScenario catchup_scenario;
    MechanismScenario mechanism_scenario;
    IntensityScenario intensity_scenario;

    friend ostream& operator<<(ostream& os, const Scenario& s) {
        os << s.vaccine_scenario << " ";
        os << s.catchup_scenario << " ";
        os << s.mechanism_scenario << " ";
        os << s.intensity_scenario << " ";
        return os;
    }
};


bool read_scenarios_from_database (string database_filename, map<string, Scenario*> &scenarios, int vac, int catchup, int vac_mech, int beta_mult) {
    sqdb::Db db(database_filename.c_str());

    // make sure database looks intact
    if ( !db.TableExists("jobs") or !db.TableExists("parameters") or !db.TableExists("metrics") ) {
        cerr << "ERROR: Failed to read SMC set from database because one or more tables are missing.\n";
        return false;
    }
//need to select the 1000 rows for a particular scenario
    stringstream select_ss;
    select_ss << "select seed from parameters where vac = " << vac << " and catchup = " << catchup << " and vac_mech = " << vac_mech << " and beta_multiplier = " << beta_mult << ";";
    sqdb::Statement s = db.Query( select_ss.str().c_str() );

    while (s.Next()) {
        const string seed = s.GetField(0);
        assert(scenarios.find(seed) == scenarios.end()); // assert: this is not a seed collision
        scenarios[seed] = new Scenario(vac, catchup, vac_mech, beta_mult);
    }

    return true;
}


void initialize_aggregate_datastructure(map<string, Scenario*>& scenarios, AgType& ag, map<string, int>& N) {
    unordered_set<string> uniqe_scenarios;
    for (auto scen: scenarios) {
        string scenario = (scen.second)->asKey(); 
        uniqe_scenarios.insert(scenario); 
    } 

    // dimensions: scenario, serotype, age, outcome, year
    for (auto scenario: uniqe_scenarios) {
        N[scenario] = 0;
        ag[scenario] = vector< vector< vector< vector<double> > > >
                       ((int) NUM_OF_SEROTYPES, vector< vector< vector<double> > >(
                                       MAX_AGE+1, vector< vector<double> >(
                                               (int) NUM_OF_SEVERITY_TYPES, vector<double>(
                                                       INTERVENTION_DURATION, 0.0))));
    }
}

/*
asymtomatic, mild, severe (3)
by
age (100)
by
year (30)
by
serotype (4)
by

mechanism (2)
by
catchup (2)
by
intensity (5)
*/

void report_average_counts(AgType& ag, map<string, int>& N, string filename) {
    string all_output  = "serotype,age,outcome,y00,y01,y02,y03,y04,y05,";
           all_output += "y06,y07,y08,y09,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29\n";
    for (auto &scenario: ag) {
        const int divisor = N[scenario.first];
        int sero_ct = 0;
        for (auto &serotype: scenario.second) {
            int age_ct = 0;
            for (auto &age: serotype) {
                int out_ct = 0;
                for (auto &outcome: age) {
                    string line = scenario.first + ","
                                  + to_string(sero_ct) + ","
                                  + to_string(age_ct) + ","
                                  + to_string(out_ct);
                    for (auto &val: outcome) { val /= divisor; line += ("," + to_string(val)); }
                    all_output += (line + "\n");
                    ++out_ct;
                }
                ++age_ct;
            } 
            ++sero_ct;
        }
    }

    if (fileExists(filename)) {
        cerr << "ERROR: Digest output file already exists: " << filename << endl << "ERROR: Aborting write.\n";
        return;
    }

    ofstream file;
    file.open(filename);

    if (file.is_open()) {
        file << all_output;
        file.close();
    } else {
        cerr << "ERROR: Could not open daily buffer file for output: " << filename << endl;
        exit(-842);
    }
}


void process_daily_files(map<string, Scenario*> scenarios, string daily_dir, string output_dir) {
    int file_ctr = 0;
    vector<string> daily_filenames = glob(daily_dir + "/daily.*");
    AgType ag;
    map<string, int> N;
    initialize_aggregate_datastructure(scenarios, ag, N);

    cerr << "globbed " << daily_filenames.size() << " filenames\n";
    #pragma omp parallel for num_threads(10) 
    for (unsigned int i = 0; i<daily_filenames.size(); ++i) {
        string daily_filename = daily_filenames[i];
        const string seed = extract_seed_from_filename(daily_filename);

        if (scenarios.count(seed) == 0) {
            //if (i < 10) cerr << "no match: <" << seed << ">\n";
            continue;
        } else {
            cerr << file_ctr << " " << daily_filename << endl;
            Scenario* scenario = scenarios[seed];
            const string scenarioKey = scenario->asKey();
            #pragma omp atomic
            ++N[scenarioKey];

            // Read entire file into a stringstream
            ifstream fin(daily_filename.c_str());
            if (!fin) { cerr << "ERROR: Could not open " << daily_filename << endl; exit(114); }
            stringstream iss;
            copy(istreambuf_iterator<char>(fin), istreambuf_iterator<char>(), ostreambuf_iterator<char>(iss));
            fin.close();

            char buffer[500];
            istringstream line(buffer);

            // day,year,id,age,location,vaccinated,serotype,symptomatic,severity
            // 30,0,134867,90,29256,0,1,1,0
            int day;
            int year;
            int id;
            int age;
            int location;
            int vaccinated;
            int serotype;
            int symptomatic;
            int severe;
            while (iss) {
                iss.getline(buffer,500);
                line.clear();
                string line_str(buffer);
                replace(line_str.begin(), line_str.end(), ',', ' ');
                line.str(line_str);
//                           11832    32    69242   15     14685        1              2            0          0
                if (line >> day >> year >> id >> age >> location >> vaccinated >> serotype >> symptomatic >> severe) {
                    if (year < BURNIN) continue;
//                    cerr << line_str << endl;
                    const int sero = serotype - 1;
                    const int y = year - BURNIN;
                    const int outcome = severe==1 ? 2 : symptomatic==1 ? 1 : 0;
//                    cerr << "a";
//                    assert(sero < (int) NUM_OF_SEROTYPES); 
//                    assert(age <= MAX_AGE);
//                    assert(outcome < (int) NUM_OF_SEVERITY_TYPES);
//                    assert(year < INTERVENTION_DURATION);

                    #pragma omp atomic
                    ++ag[scenarioKey][sero][age][outcome][y];
//                    cerr << "b";
                } else {
                    continue; // Didn't process line, normal for header or EOF
                    //cerr << "WARNING: Could not parse line: " << line.str() << endl;
                }
 //               cerr << "c";
            }
            #pragma omp atomic
            ++file_ctr;
        }

        stringstream output_filename_ss;
        output_filename_ss << output_dir << "/daily_agg." << scenarios.begin()->second->asKey();
        report_average_counts(ag, N, output_filename_ss.str());
        cout << "read files: " << file_ctr << endl;
    }
}


void usage() {
    cerr << "\n\tUsage: ./process_who database_filename path_to_daily_files path_for_output vaccine_par catchup_par vaccine_mechanism_par beta_multiplier_par\n\n";
}


int main (int argc, char* argv[]) {

    /*if (not (argc == 4) ) {
        usage();
        exit(100);
    }*/

    string db_filename = string(argv[1]);
    string daily_dir   = string(argv[2]);
    string output_dir  = string(argv[3]);
    int vac            = atoi(argv[4]);
    int catchup        = atoi(argv[5]);
    int vac_mech       = atoi(argv[6]);
    int beta_mult      = atoi(argv[7]);

    /*if (fileExists(digest_filename)) {
        cerr << "ERROR: Digest output file already exists: " << digest_filename << endl << "ERROR: Aborting write.\n";
        exit(101);
    }*/


    map<string, Scenario*> scenarios;
    
    if (read_scenarios_from_database(db_filename, scenarios, vac, catchup, vac_mech, beta_mult)) {
        process_daily_files(scenarios, daily_dir, output_dir); 
    } else {
        cerr << "Failed to read scenarios\n"; 
    }
    int counter = 0;
    for (auto row: scenarios) {
        if (++counter > 10) break;
        cout << row.first << " " << *(row.second) << endl; // database read sanity check
    }
}
