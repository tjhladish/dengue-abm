#include <iostream>
#include <vector>
#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_roots.h"

using namespace std;


vector<double> mu =  {0.03109522, 0.03103087, 0.03095702, 0.03087250, 0.03077550, 0.03066470,
                      0.03053790, 0.03039310, 0.03022800, 0.03004010, 0.02982630, 0.02958380,
                      0.02930930, 0.02899930, 0.02865010, 0.02825820, 0.02781960, 0.02733100,
                      0.02678880, 0.02619040, 0.02553320, 0.02481570, 0.02403760, 0.02319950,
                      0.02230360, 0.02135350, 0.02035470, 0.01931370, 0.01823920, 0.01714060,
                      0.01602820, 0.01491330, 0.01380680, 0.01271960, 0.01166160, 0.01064190,
                      0.00966800, 0.00874580, 0.00788000, 0.00707320, 0.00632700, 0.00564130,
                      0.00501520, 0.00444650, 0.00393300, 0.00347090, 0.00305740, 0.00268810,
                      0.00236000, 0.00206890, 0.00181160, 0.00158450, 0.00138460, 0.00120880,
                      0.00105470, 0.00091950, 0.00080130, 0.00069780, 0.00159620, 0.0};

double alpha_irs_expr (double rho, void *params) {
    double alpha_irs = *(double *) params;

    double alpha_irs_guess = 1.0;
    for (unsigned int k = 1; k < mu.size(); ++k) { // k in notation goes from 1..N, but C++ goes from 0..N-1, no?
        double partial = 1.0;
        for (unsigned int j = 0; j < k; ++j) {       // I think indexing should match outer loop, but unsure
            partial *= 1.0 - mu[j];
        }
        
        // alpha_irs ==: 1+\sum_{k=2}^{N}(1-rho)^{k-1}\prod_{j=1}^{j=k-1} (1-mu_j)
        alpha_irs_guess += pow(1.0 - rho, k) * partial;
    }
    cerr << rho << "\t" << alpha_irs_guess << endl;
    return alpha_irs_guess - alpha_irs;
}

int find_rho (double alpha_irs) {
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double rho = 0;
    double rho_lo = 0.0, rho_hi = 1.0;
    gsl_function F;

    F.function = &alpha_irs_expr;
    F.params = &alpha_irs;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, rho_lo, rho_hi);

    printf ("using %s method\n", gsl_root_fsolver_name (s));

    printf ("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "root", "err(est)");

    do {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        rho = gsl_root_fsolver_root (s);
        rho_lo = gsl_root_fsolver_x_lower (s);
        rho_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (rho_lo, rho_hi, 0, 0.001);

        if (status == GSL_SUCCESS) { printf ("Converged:\n"); }

        printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n", iter, rho_lo, rho_hi, rho, rho_hi - rho_lo);
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);

    return status;
}

//    p_i = asymptotic mosquito population in age (day) cohort i
//    p_i* = ibid, but w/ IRS
//    mu_i = mortality prob in age cohort i // extracted from current CDF
//    rho = daily IRS probability (fit parameter)
//    p_0 = emerging population
//    P = asymptotic population at location
//    P* = ibid, w/ IRS
//    N = max age // length(mu)
//
//    p_{i+1} = (1-mu_{i})*p_{i}
//
//    p_1 = p_1* = p_0
//    => p_i = p_0\prod_{j=1}^{j=i-1} (1-mu_j)
//
//    P = \sum_k p_k = p_0(1+\sum_k\prod_{j=1}^{j=k-1} (1-mu_j)) = p_0 \alpha // k goes from 1..N
//    (if mu const, then \alpha is a simple geo series, but that's not the case for us)
//
//    p_{i+1}* = (1-mu_{i})(1-rho)p_{i}*
//    => p_i* = p_0(1-rho)^{i-1}\prod_{j=1}^{j=i-1}(1-mu_j)
//    P* = p_0(1+\sum_{k=2}^{N}(1-rho)^{k-1}\prod_{j=1}^{j=k-1} (1-mu_j)) = p_0 \alpha* = P(\alpha* / \alpha)
//
//    P(age == i)* = p_i*/P* =  (1-rho)^{i-1}\prod_{j=1}^{j=i-1}(1-mu_j) / \alpha*

int main() {
    const float efficacy = 0.9;
    // calculate alpha
    double alpha = 1.0; // accounts for k == 0 term
    //cerr << mu.size() << endl;
    for (unsigned int k = 1; k < mu.size(); ++k) { // k in notation goes from 1..N, but C++ goes from 0..N-1, no?
        double partial = 1.0;
        for (unsigned int j = 0; j < k; ++j) {       // I think indexing should match outer loop, but unsure
            partial *= 1.0 - mu[j];
        }
        alpha += partial;
    }
    
    cerr << "alpha: " << alpha << endl;
    const double alpha_irs = (1.0 - efficacy) * alpha;
    // calculate rho
    // alpha_irs ==: 1+\sum_{k=2}^{N}(1-rho)^{k-1}\prod_{j=1}^{j=k-1} (1-mu_j)
    find_rho(alpha_irs);

}
