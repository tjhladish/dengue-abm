#include "simulator.h"
#include <time.h>

int main(int argc, char* argv[]) {
    time_t start ,end;
    time (&start);

    srand(time(NULL));
    int proccess_id = rand();
    fprintf(stderr, "%dbegin\n", proccess_id);

    const Parameters* par = new Parameters(argc, argv);

    Community* community = build_community(par);
    vector<int> initial_susceptibles = community->getNumSusceptible();
    seed_epidemic(par, community);
    vector<int> epi_sizes = simulate_epidemic(par, community);

    const double ef = par->expansionFactor;
    vector<double> x(epi_sizes.size());
    vector<double> y(epi_sizes.size());
    // convert infections to cases
    for (unsigned int i = 0; i < epi_sizes.size(); i++) { 
        y[i] = ((double) epi_sizes[i])/ef; 
        x[i] = i+1.0;
    }

    Fit* fit = lin_reg(x, y);
    
    time (&end);
    double dif = difftime (end,start);

    stringstream ss;
    // wallclock time (seconds)
    ss << proccess_id << "end " << dif << " ";
    // parameters
    ss << par->expansionFactor << " " << par->fMosquitoMove << " " << par->nDailyExposed[0][0] << " "
       << par->betaMP << " " << par->betaPM << " ";
    // metrics
    ss << mean(y) << " " << stdev(y) << " " << max_element(y) << " "
       << fit->m << " " << fit->b << " " << fit->rsq;

    ss << endl;
    string output = ss.str();
    fputs(output.c_str(), stderr);
    
    // metrics to stdout
    cout << mean(y) << " " << stdev(y) << " " << max_element(y) << " "
       << fit->m << " " << fit->b << " " << fit->rsq;


    write_output(par, community, initial_susceptibles);
   
    return 0;
}
