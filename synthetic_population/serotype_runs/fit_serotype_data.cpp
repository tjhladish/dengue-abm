#include <iostream>
#include <vector>
#include <numeric>
#include <unistd.h>
#include "AbcSmc.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


using namespace std;

const int SEQUENCE_LENGTH = 36;
const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);
const vector<int> INTRO_YEARS = {0, 7, 17, 5}; // First year that each serotype appears, respectively
const vector<int> UNOBSERVED_YEARS = {19, 20, 21, 24};
int serotype = 1;

vector<long double> revisedSerotypeGenerator(vector<long double> args, const unsigned long int rng_seed, const MPI_par* mp) {
    gsl_rng_set(RNG, rng_seed);
    //const int serotype = (int) args[0] - 1;
    const long double geometric_mean_run = args[0];
    const long double geometric_mean_gap = args[1];
    const float p_gap = 1.0/geometric_mean_gap;
    const float p_run = 1.0/geometric_mean_run;

    vector<long double> obs_run_and_gap_means(2, 0.0);

    enum State {GAP, RUN, NEITHER};
    vector<State> series(INTRO_YEARS[serotype], NEITHER);

    while (series.size() < SEQUENCE_LENGTH) {
        unsigned int run = gsl_ran_geometric(RNG, p_run);
        series.resize(series.size() + run, RUN);
        unsigned int gap = gsl_ran_geometric(RNG, p_gap);
        series.resize(series.size() + gap, GAP);
    }
    series.resize(SEQUENCE_LENGTH);

    for (auto year: UNOBSERVED_YEARS) series[year] = NEITHER;

    vector<int> runs;
    vector<int> gaps;

    State last_state = NEITHER;
    int tally = 0;
    for (auto s: series) {
        if (s == last_state) {          // continuing
            ++tally;
        } else if (last_state == RUN) { // end of a run
            runs.push_back(tally);
            tally = 1;
        } else if (last_state == GAP) { // end of a gap
            gaps.push_back(tally);
            tally = 1;
        } else {
            tally = 1;
        }
        last_state = s;
    }

    if (last_state == RUN) { // ends on a run
            runs.push_back(tally);
    } else if (last_state == GAP) {
            gaps.push_back(tally);
    }

    for (auto i: series) cerr << i << " "; cerr << endl;
    for (auto i: runs) cerr << i << " "; cerr << endl;
    for (auto i: gaps) cerr << i << " "; cerr << endl;

    if (runs.size() > 0) {
        obs_run_and_gap_means[0] = accumulate(runs.begin(), runs.end(), 0.0) / runs.size();
    } else {
        obs_run_and_gap_means[0] = 1; // minimum possible value
    }

    if (gaps.size() > 0) {
        obs_run_and_gap_means[1] = accumulate(gaps.begin(), gaps.end(), 0.0) / gaps.size();
    } else {
        obs_run_and_gap_means[1] = 1;
    }

    return obs_run_and_gap_means;
}


void usage() {
    cerr << "\n\tUsage: ./abc_sql abc_config_sql.json --process\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate -n <number of simulations per database write>\n\n";
}


int main(int argc, char* argv[]) {

    if (not (argc >= 3) ) {
        usage();
        exit(100);
    }

    bool process_db  = false;
    bool simulate_db = false;
    int buffer_size  = -1;

    for (int i=2; i < argc;  i++ ) {
        if ( strcmp(argv[i], "--process") == 0  ) {
            process_db = true;
        } else if ( strcmp(argv[i], "--simulate") == 0  ) {
            simulate_db = true;
            buffer_size = buffer_size == -1 ? 1 : buffer_size;
        } else if ( strcmp(argv[i], "--serotype" ) == 0 ) {
            serotype = atoi(argv[++i]) - 1;
        } else if ( strcmp(argv[i], "-n" ) == 0 ) {
            buffer_size = atoi(argv[++i]);
        } else {
            usage();
            exit(101);
        }
    }

    AbcSmc* abc = new AbcSmc();
    abc->parse_config(string(argv[1]));
    if (process_db) {
        gsl_rng_set(RNG, time(NULL) * getpid()); // seed the rng using sys time and the process id
        abc->process_database(RNG);
    }

    if (simulate_db) {
        abc->set_simulator(revisedSerotypeGenerator);
        abc->simulate_next_particles(buffer_size);
    }

    return 0;
}
