#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <assert.h>
#include <glob.h>
#include "sqdb.h"

using namespace std;

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


inline unsigned int extract_seed_from_filename(const string& s) {
    string::size_type p  = s.find('.');
    string::size_type pp = s.find('.', p + 2); 
    return stoul(s.substr(p + 1, pp - p - 1));
}

enum SerotypeState {SERO1, SERO2, SERO3, SERO4, NUM_OF_SEROTYPES};
enum SeverityState {ASYMPTOMATIC, NO_HOSPITAL, HOSPITAL, NUM_OF_SEVERITY_TYPES};
enum ImmunityState {PRIMARY, SECONDARY, POST_SECONDARY, NUM_OF_IMMUNITY_TYPES};

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


void process_daily_files(map<int, Scenario*> scenarios, string daily_path, string digest_filename) {
const int max_file_ct = 100;
int file_ctr = 0;
int line_ctr = 0;
    vector<string> daily_filenames = glob(daily_path + "/daily.*");
    
    for (string daily_filename: daily_filenames) {
        //cerr << daily_filename << " " << extract_seed_from_filename(daily_filename) << endl;

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
                int year = day / 365;
                bool isCase = case_char == 'T' ? true : false;
                bool isHosp = hosp_char == 'T' ? true : false;
                int num_infections = count(hist.begin(), hist.end(), '+');
line_ctr++;
                //cerr << line_str << " | " << year << " " << age << " " << isCase << " " << isHosp << " " << num_infections << endl;
                //continue;
            } else {
                cerr << "WARNING: Could not parse line: " << line.str() << endl;
            }
        }


        iss.close();
        //exit(-1); 
if (++file_ctr >= max_file_ct) break;
    }
cerr << "read files, lines: " << file_ctr << ", " << line_ctr << endl;
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
    string daily_path      = string(argv[2]);
    string digest_filename = string(argv[3]);

    map<int, Scenario*> scenarios;
    
    if (read_scenarios_from_database(db_filename, scenarios)) {
        process_daily_files(scenarios, daily_path, digest_filename); 
    }
    //for (auto row: scenarios) cout << row.first << " " << *(row.second) << endl; // database read sanity check
}
