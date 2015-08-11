#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <unordered_set>
#include <string>
#include <algorithm>
#include <assert.h>
#include <glob.h>
#include <sys/stat.h>
#include "sqdb.h"

using namespace std;

typedef map< string, vector< vector< vector< vector< vector<double> > > > > > AgType;
const int MAX_AGE = 100;
const int YEARS_SIMULATED = 20;

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


inline unsigned int extract_seed_from_filename(const string& s) {
    string::size_type p  = s.find('.');
    string::size_type pp = s.find('.', p + 2); 
    return stoul(s.substr(p + 1, pp - p - 1));
}

enum SerotypeState {SERO1, SERO2, SERO3, SERO4, NUM_OF_SEROTYPES};
enum SeverityState {ASYMPTOMATIC, NO_HOSPITAL, HOSPITAL, NUM_OF_SEVERITY_TYPES};
enum HistoryState {PRIMARY, SECONDARY, POST_SECONDARY, NUM_OF_HISTORY_TYPES};

enum VaccineScenario {VACCINE, NO_VACCINE, NUM_OF_VACCINE_TYPES};
enum CatchupScenario {CATCHUP, NO_CATCHUP, NUM_OF_CATCHUP_TYPES};
enum WaningScenario {WANING, NO_WANING, NUM_OF_WANING_TYPES};
enum BoostingScenario {BOOSTING, NO_BOOSTING, NUM_OF_BOOSTING_TYPES};

struct Scenario {
    Scenario(int v, int c, int w, int b) {
        assert(v < NUM_OF_VACCINE_TYPES);
        assert(c < NUM_OF_CATCHUP_TYPES);
        assert(w < NUM_OF_WANING_TYPES);
        assert(b < NUM_OF_BOOSTING_TYPES);

        vaccine_scenario = (VaccineScenario) v; 
        catchup_scenario = (CatchupScenario) c;
        waning_scenario = (WaningScenario) w;
        boosting_scenario = (BoostingScenario) b;
    }

    string asKey() {
        stringstream ss;
        ss << vaccine_scenario
           << catchup_scenario
           << waning_scenario
           << boosting_scenario;
        return ss.str();
    }

    VaccineScenario vaccine_scenario;
    CatchupScenario catchup_scenario;
    WaningScenario waning_scenario;
    BoostingScenario boosting_scenario;

    friend ostream& operator<<(ostream& os, const Scenario& s) {
        os << s.vaccine_scenario << " ";
        os << s.catchup_scenario << " ";
        os << s.waning_scenario << " ";
        os << s.boosting_scenario << " ";
        return os;
    }
};

bool read_scenarios_from_database (string database_filename, map<int, Scenario*> &scenarios) {
    sqdb::Db db(database_filename.c_str());

    // make sure database looks intact
    if ( !db.TableExists("jobs") or !db.TableExists("parameters") or !db.TableExists("metrics") ) {
        cerr << "ERROR: Failed to read SMC set from database because one or more tables are missing.\n";
        return false;
    }

    // join all three tables for rows with smcSet = t, slurp and store values
    string select_str = "select seed, vac, catchup, vac_waning, vac_boosting from parameters;";

    sqdb::Statement s = db.Query( select_str.c_str() );

    while (s.Next()) {
        const int seed      = s.GetField(0);
        const int vac       = s.GetField(1);
        const int catchup   = s.GetField(2);
        const int waning    = s.GetField(3);
        const int boosting  = s.GetField(4);

        assert(scenarios.find(seed) == scenarios.end()); // not a seed collision
        
        scenarios[seed] = new Scenario(vac, catchup, waning, boosting);
    }

    return true;
}


void initialize_aggregate_datastructure(map<int, Scenario*>& scenarios, AgType& ag, map<string, int>& N) {
    unordered_set<string> uniqe_scenarios;
    for (auto scen: scenarios) {
        string scenario = (scen.second)->asKey(); 
        uniqe_scenarios.insert(scenario); 
    } 

         //  scenario      serotype  history  age  outcome year
//typedef map< string, map< int, map< int, map< int, map< int, vector<int> > > > > > AgType;
    for (auto scenario: uniqe_scenarios) {
        N[scenario] = 0;
        ag[scenario] = vector< vector< vector< vector< vector<double> > > > >
                       ((int) NUM_OF_SEROTYPES, vector< vector< vector< vector<double> > > >(
                               (int) NUM_OF_HISTORY_TYPES, vector< vector< vector <double> > >(
                                       MAX_AGE+1, vector< vector<double> >(
                                               (int) NUM_OF_SEVERITY_TYPES, vector<double>(
                                                       YEARS_SIMULATED, 0.0)))));
    }
}


void report_average_counts(AgType& ag, map<string, int>& N, string filename) {
    string all_output  = "scenario,serotype,history,age,outcome,y00,y01,y02,y03,y04,y05,";
           all_output += "y06,y07,y08,y09,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19\n";
    for (auto &scenario: ag) {
        const int divisor = N[scenario.first];
        int sero_ct = 0;
        for (auto &serotype: scenario.second) {
            int hist_ct = 0;
            for (auto &history: serotype) {
                int age_ct = 0;
                for (auto &age: history) {
                    int out_ct = 0;
                    for (auto &outcome: age) {
                        //for (auto val: outcome) cout << val << " "; cout << endl;
                        string line = scenario.first + ","
                                      + to_string(sero_ct) + ","
                                      + to_string(hist_ct) + ","
                                      + to_string(age_ct) + ","
                                      + to_string(out_ct);
                        for (auto &val: outcome) { val /= divisor; line += ("," + to_string(val)); }
                        //for (auto val: outcome) cout << val << " "; cout << endl;
                        all_output += (line + "\n");
                        ++out_ct;
                    }
                    ++age_ct;
                }
                ++hist_ct;
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

    if (file.is_open()) {  // TODO - add this check everywhere a file is opened
        file << all_output;
        file.close();
    } else {
        cerr << "ERROR: Could not open daily buffer file for output: " << filename << endl;
        exit(-842);
    }
}


void process_daily_files(map<int, Scenario*> scenarios, string daily_path, string digest_filename) {
    int file_ctr = 0;
    //vector<string> daily_filenames = glob(daily_path);
    vector<string> daily_filenames = glob(daily_path + "/daily.*");
    AgType ag;
    map<string, int> N;
    initialize_aggregate_datastructure(scenarios, ag, N);

    #pragma omp parallel for num_threads(3) 
    for (unsigned int i = 0; i<daily_filenames.size(); ++i) {
        string daily_filename = daily_filenames[i];
   // for (string daily_filename: daily_filenames) {
        cerr << file_ctr << " " << daily_filename << endl;
        Scenario* scenario = scenarios[ extract_seed_from_filename(daily_filename) ];
        const string scenarioKey = scenario->asKey();
        #pragma omp atomic
        ++N[scenarioKey];

        ifstream iss(daily_filename.c_str());
        if (!iss) {
            cerr << "ERROR: Could not open " << daily_filename << endl;
            exit(114);
        }
        // day,id,age,sero,case,hosp,hist,lastvac
        // 0,1497628,6,0,F,F,---+,-1
        char buffer[500];
        istringstream line(buffer);
        
        int day;
        int id;
        int age;
        int sero;
        char case_char;
        char hosp_char;
        string hist;
        int lastvac;

        while (iss) {
            iss.getline(buffer,500);
            line.clear();
            string line_str(buffer);
            replace(line_str.begin(), line_str.end(), ',', ' ');
            line.str(line_str);
            if (line >> day >> id >> age >> sero >> case_char >> hosp_char >> hist >> lastvac) {
                const int year = day / 365;
                const bool isCase = case_char == 'T' ? true : false;
                const bool isHosp = hosp_char == 'T' ? true : false;
                int num_prev_infections = count(hist.begin(), hist.end(), '+') - 1; // -1 b/c current one is being reported
                num_prev_infections = num_prev_infections > 2 ? 2 : num_prev_infections;
                const int outcome = isHosp ? 2 : isCase ? 1 : 0; 
                //  scenario      serotype  history  age  outcome year
                //cerr << scenario->asKey() << " " << sero << " " << num_prev_infections << " " << age << " " << outcome << " " << year << endl;;
                #pragma omp atomic
                ++ag[scenarioKey][sero][num_prev_infections][age][outcome][year];
                
//line_ctr++;
                //cerr << line_str << " | " << year << " " << age << " " << isCase << " " << isHosp << " " << num_infections << endl;
                //continue;
            } else {
                continue;
                //cerr << "WARNING: Could not parse line: " << line.str() << endl;
            }
        }


        iss.close();
        //exit(-1); 
        #pragma omp atomic
        ++file_ctr;
//if (++file_ctr >= max_file_ct) break;
    }

    report_average_counts(ag, N, digest_filename);
cout << "read files: " << file_ctr << endl;
}


void usage() {
    cerr << "\n\tUsage: ./process_who database_filename path_to_daily_files digest_filename\n\n";
}



int main (int argc, char* argv[]) {

    if (not (argc == 4) ) {
        usage();
        exit(100);
    }

    string db_filename     = string(argv[1]);
    string daily_path= string(argv[2]);
    string digest_filename = string(argv[3]);

    if (fileExists(digest_filename)) {
        cerr << "ERROR: Digest output file already exists: " << digest_filename << endl << "ERROR: Aborting write.\n";
        exit(101);
    }


    map<int, Scenario*> scenarios;
    
    if (read_scenarios_from_database(db_filename, scenarios)) {
        process_daily_files(scenarios, daily_path, digest_filename); 
    }
    //for (auto row: scenarios) cout << row.first << " " << *(row.second) << endl; // database read sanity check
}
