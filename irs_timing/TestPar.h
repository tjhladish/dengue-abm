#ifndef __TESTPAR_H
#define __TESTPAR_H

#include <algorithm>
#include <numeric>
#include <functional>
#include <vector>

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

std::vector<double> cdf_from_relfrac(std::vector<double> rf) {
  return cdf_from_pdf(norm_weights_to_pdf(rf));
};

// the cdf of mosquito age
static const std::vector<double> MOSQUITO_AGE_CDF = cdf_from_relfrac(MOSQUITO_AGE_RELFRAC);

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

#endif
