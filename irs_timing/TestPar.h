#ifndef __TESTPAR_H
#define __TESTPAR_H

#include <algorithm>
#include <numeric>
#include <functional>
#include <vector>
#include <gsl/gsl_roots.h>

// probability of death at an age, given survival to that age
static const std::vector<double> MOSQUITO_DAILY_DEATH_PROBABILITY = {
  0.0018,0.002069514,0.002378646,0.002732982,0.003138823,0.003603251,0.004134191,
  0.004740476,0.005431895,0.006219231,0.007114273,0.008129798,0.009279508,
  0.01057792,0.01204018,0.01368180,0.01551828,0.01756467,0.01983494,0.02234135,
  0.02509359,0.02809799,0.03135655,0.03486617,0.03861782,0.04259607,0.04677877,
  0.05113718,0.05563645,0.06023658,0.0648937,0.0695617,0.07419404,0.07874562,
  0.08317441,0.08744304,0.09151983,0.09537952,0.09900354,0.1023799,0.1055030,
  0.1083723,0.1109923,0.1133714,0.1155206,0.1174531,0.1191837,0.1207277,0.1221007,
  0.1233179,0.1243943,0.1253438,0.1261799,0.1269147,0.1275594,0.1281243,0.1286187,
  0.1290510,0.1294285,0.1297580,0.1300453, 1
};

std::vector<double> norm_weights_to_pdf(std::vector<double> weights) {
  std::vector<double> pdf(weights.size());
  auto total = std::accumulate(
    weights.begin(), weights.end(), // take all of weights
    0.0 // starting with an initial value of 0.0
    // add everything up (default accumulate op is +)
  ); // i.e., the sum of weights
  std::transform(
    weights.begin(), weights.end(), // for all weights
    pdf.begin(), // into pdf
    [total](double rel){ return rel/total; }
    // the probability of particular index
    // is the weight of that index, divided by total
  );
  return pdf;
};

std::vector<double> cdf_from_pdf(std::vector<double> pdf) {
  std::vector<double> cdf(pdf.size());
  std::partial_sum(
    pdf.begin(), pdf.end(),
    cdf.begin()
  );
  return cdf;
};

std::vector<double> complement(std::vector<double> ps) {
  std::vector<double> res(ps.size());
  std::transform(
    ps.begin(), ps.end(), // for all ps
    res.begin(), // put into res
    [](double p) { return 1-p; } // after taking complement
  );
  return res;
};

// probability of survival at an age, given survival to that age
static const std::vector<double> MOSQUITO_DAILY_SURVIVAL_PROBABILITY =
  complement(MOSQUITO_DAILY_DEATH_PROBABILITY);

std::vector<double> cumprod(std::vector<double> ps) {
  std::vector<double> res(ps.size());
  std::partial_sum(
    ps.begin(), ps.end(), // for all p in ps
    res.begin(), // put into res
    std::multiplies<double>()
    // the cumulative product at each p (cp_i = prod from 0 to i p_i)
  );
  return res;
};

// probability of surviving an age from birth
static const std::vector<double> MOSQUITO_SURVIVE_CUMPROB =
  cumprod(MOSQUITO_DAILY_SURVIVAL_PROBABILITY);

std::vector<double> relative_fraction(std::vector<double> ps) {
  std::vector<double> res(ps.size(), 1.0);
  std::copy(
    ps.begin(), ps.end()-1,
    // survival fraction == the fraction that survived previous day
    // so day 1 fraction = 1, day 2 fraction = day 1 survival prob from birth, etc
    res.begin()+1
  );
  return res;
};

// the relative fraction of mosquitos from a birth cohort alive at an age
// also, given no migration etc, the unnormalized steady-state age distribution
static const std::vector<double> MOSQUITO_AGE_RELFRAC =
  relative_fraction(MOSQUITO_SURVIVE_CUMPROB);

// the pdf of mosquito age
static const std::vector<double> MOSQUITO_AGE_PDF = norm_weights_to_pdf(MOSQUITO_AGE_RELFRAC);

// the cdf of mosquito age
static const std::vector<double> MOSQUITO_AGE_CDF = cdf_from_pdf(MOSQUITO_AGE_PDF);

std::vector<double> death_age_cdf(std::vector<double> survive_age_prob, std::vector<double> die_age_prob) {
  std::vector<double> pdf(die_age_prob);

  std::transform(
    survive_age_prob.begin(), survive_age_prob.end()-1,
    // using survival probs for previous day
    die_age_prob.begin()+1, // deaths probs for same day
    pdf.begin()+1, // overwriting from day 2 on
    std::multiplies<double>() // combine the probabilities
  );

  return cdf_from_pdf(pdf);
};


static const std::vector<double> MOSQUITO_DEATHAGE_CDF =
  death_age_cdf(MOSQUITO_SURVIVE_CUMPROB, MOSQUITO_DAILY_DEATH_PROBABILITY);

static const double Mref =
  std::accumulate(MOSQUITO_AGE_RELFRAC.begin(), MOSQUITO_AGE_RELFRAC.end(), 0.0);

struct RhoParams { double eff };

double mstar (double mu, vector<double> phat) {
  vector<int> pows(phat.size());
  std::iota(pows.begin(), pows.end(), 0);
  vector<double> survives(pows.size());
  std::transform(
    pows.begin(), pows.end(),
    phat.begin(),
    survives.begin(),
    [mu](int power, double phat){ return phat*pow(1-mu, power); }
  );
  return std::accumulate(survives.begin(), survives.end(), 0.0);
}

double _irs_expr (double rho_par, void *params) {
    // Logit transformation, for potentially more efficient search of par space.
    // In testing it was 5% slower, however, because of the extra calculations for
    // the transformation.  Convergence did not happen in substantially fewer iterations.
    //const double rho = 1.0 / (1.0 + exp(-rho_par));

    RhoParams* rho_params = (RhoParams *) params;

    //cerr << rho_par << "\t" << rho << "\t" << alpha_irs_guess << "\t" << alpha_irs << " diff: " << alpha_irs_guess - alpha_irs << endl;
    return mstar(rho_par, MOSQUITO_AGE_RELFRAC) - Mref*(1-rho_params->eff);
}


  double __find_rho (RhoParams* rho_params) {
      int status = GSL_CONTINUE;
      int iter = 0, max_iter = 100;
      const gsl_root_fsolver_type *T;
      gsl_root_fsolver *s;
      //double rho_lo = -5; // for logit search
      //double rho_hi = 5;
      double rho_lo = 0.0;
      double rho_hi = 1.0;
      double rho = (rho_hi + rho_lo)/2;
      gsl_function F;

      F.function = &_irs_expr;
      F.params = rho_params;

      T = gsl_root_fsolver_brent;
      s = gsl_root_fsolver_alloc (T);
      gsl_root_fsolver_set (s, &F, rho_lo, rho_hi);

      //printf ("using %s method\n", gsl_root_fsolver_name (s));
      //printf ("%5s [%9s, %9s] %9s %9s\n", "iter", "lower", "upper", "root", "err(est)");
      while (status == GSL_CONTINUE and iter++ < max_iter) {
          status = gsl_root_fsolver_iterate (s);
          rho = gsl_root_fsolver_root (s);
          rho_lo = gsl_root_fsolver_x_lower (s);
          rho_hi = gsl_root_fsolver_x_upper (s);
          status = gsl_root_test_interval (rho_lo, rho_hi, 0, 0.001);

       //   if (status == GSL_SUCCESS) { printf ("Converged:\n"); }

       //   const double rho = 1.0 / (1.0 + exp(-rho_par));
       //   printf ("%5d [%.7f, %.7f] %.7f %+.7f\n", iter, rho_lo, rho_hi, rho, rho_hi - rho_lo);
      }

      gsl_root_fsolver_free (s);
      if (not (status == GSL_SUCCESS)) {
          cerr << "ERROR: Root-finding (for daily mosquito mortality due to IRS) did not converge\n";
          exit(-732);
      }
      return rho;
  }


double calculate_daily_vector_control_mortality(const float efficacy) const {

    RhoParams* rho_params = new RhoParams(efficacy);

    //cerr << "max supported efficacy given irs model: " << 1.0 - (1.0/Mref) << endl;
    if ((1.0 - efficacy) * Mref < 1.0) {
        cerr << "ERROR: Requested efficacy (" << efficacy << ") is too high\n";
        return -100;
    }

    const double rho = __find_rho(rho_params);
    cerr << "IRS efficacy requested, daily mortality found: " << efficacy << ", " << rho << endl;
    delete rho_params;
    return rho;
}

// e.g. today_biting_age_pdf = weight_biting_age_pdf(MOSQUITO_AGE_PDF, today_inf_p);

std::vector<double> weight_biting_age_pdf(
  std::vector<double> raw_age_pdf,
  double prob_infecting_bite
) {
  std::vector<double> pdf(raw_age_pdf);
  std::vector<double> biting_prob(raw_age_pdf.size(), prob_infecting_bite);
  biting_prob[1] = 1.0;

  std::partial_sum(
    biting_prob.begin(), biting_prob.end(),
    biting_prob.begin(), // just overwrite vals
    [](double prev, double p) { return prev*(1-p); }
  ); // TODO: this needs some checking and is more complex than necessary for
  // uniform probability.  however: sets up having historical mosquito FOI

  std::transform(
    pdf.begin(), pdf.end(),
    biting_prob.begin(),
    pdf.begin(),
    // alt could use a series of biting probabilities
    [](double age_p, double prev_non_inf_p){ return prev_non_inf_p * age_p; }
  );

  return pdf;
};


#endif
