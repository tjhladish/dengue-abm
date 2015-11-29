#include "yucatan_serotype_generator.h"

using namespace std;

int main(int argc, char* argv[]) {
    const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);

    int first_year = 1969;
    int first_observed_year = 1979;
    int last_year  = 2013;
    vector< vector<StreakState> > all_series = generate_serotype_sequences(RNG, first_year, first_observed_year, last_year, NO_TRANSFORM);

    for(auto s: all_series) {
        for (auto v: s) cerr << v << " "; cerr << endl;
    }

    return 0;
}
