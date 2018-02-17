#include <unistd.h>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>

const gsl_rng * RNG = gsl_rng_alloc (gsl_rng_taus2);

using namespace std;

int main(int argc, char* argv[]) {

    unsigned long int seed = time (NULL) * getpid();
    cerr << "seed1: " << seed << endl;
    gsl_rng_set(RNG, seed); // seed the rng using sys time and the process id

    unsigned long int r = gsl_rng_get(RNG);
    cerr << gsl_rng_get(RNG) << endl;
    cerr << gsl_rng_get(RNG) << endl;
    cerr << "reseeding\n";// with: " << r << "\n";
    //gsl_rng_set(RNG, seed); // seed the rng using sys time and the process id
    gsl_rng_set(RNG, r); // seed the rng using sys time and the process id
    cerr << gsl_rng_get(RNG) << endl;
    cerr << gsl_rng_get(RNG) << endl;
    cerr << gsl_rng_get(RNG) << endl;
    return 0;
}
