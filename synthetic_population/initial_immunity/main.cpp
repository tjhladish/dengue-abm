#include "immunity_generator.h"
#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main(int argc, char** argv) {
    const bool write_immunity_file = false;
    if (argc != 3) {
        cerr << "\n\tUsage: ./immgen <expansion_factor> <last_year_to_simulate>\n\n";
        exit(100);
    }

    const double EXPANSION_FACTOR = atof(argv[1]);
    vector<vector<vector<int>>> full_pop = simulate_immune_dynamics(EXPANSION_FACTOR, atoi(argv[2]));

    //string immunity_filename = "/work/01856/thladish/initial_immunity/" + to_string(EXPANSION_FACTOR) + ".txt";
    if (write_immunity_file) {
        string immunity_filename = to_string(static_cast<long double>(EXPANSION_FACTOR)) + ".txt";
        output_immunity_file(immunity_filename, full_pop);
    }

    cout << "\t\t\tDone.\n";

    return 0;
}

