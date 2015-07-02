#include "simulator.h"

int main(int argc, char* argv[]) {
    const Parameters* par = new Parameters(argc, argv);

    Community* community = build_community(par);
    vector<int> initial_susceptibles = community->getNumSusceptible();
    seed_epidemic(par, community);
    simulate_epidemic(par, community);
    write_output(par, community, initial_susceptibles);
    
    return 0;
}
