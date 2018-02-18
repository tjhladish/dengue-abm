#include <algorithm>
#include <numeric>
#include <functional>
#include <vector>
#include <assert.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <chrono>
#include <thread>

using namespace std;

// probability of death at an age, given survival to that age
static const vector<double> MOSQUITO_DAILY_DEATH_PROBABILITY = { // probability of death each day
    0.0018,     0.002069514,0.002378646,0.002732982,0.003138823,0.003603251,0.004134191,0.004740476, //  0- 7
    0.005431895,0.006219231,0.007114273,0.008129798,0.009279508,0.01057792, 0.01204018, 0.01368180,  //  8-15
    0.01551828, 0.01756467, 0.01983494, 0.02234135, 0.02509359, 0.02809799, 0.03135655, 0.03486617,  // 16-23
    0.03861782, 0.04259607, 0.04677877, 0.05113718, 0.05563645, 0.06023658, 0.0648937,  0.0695617,   // 24-31
    0.07419404, 0.07874562, 0.08317441, 0.08744304, 0.09151983, 0.09537952, 0.09900354, 0.1023799,   // 32-39
    0.1055030,  0.1083723,  0.1109923,  0.1133714,  0.1155206,  0.1174531,  0.1191837,  0.1207277,   // 40-47
    0.1221007,  0.1233179,  0.1243943,  0.1253438,  0.1261799,  0.1269147,  0.1275594,  0.1281243,   // 48-55
    0.1286187,  0.1290510,  0.1294285,  1.0};//0.1297580,  0.1300453};                               // 56-59


vector<double> norm_weights_to_pdf(vector<double> weights) {
  vector<double> pdf(weights.size());
  auto total = accumulate(
    weights.begin(), weights.end(), // take all of weights
    0.0 // starting with an initial value of 0.0
    // add everything up (default accumulate op is +)
  ); // i.e., the sum of weights
  transform(
    weights.begin(), weights.end(), // for all weights
    pdf.begin(), // into pdf
    [total](double rel){ return rel/total; }
    // the probability of particular index
    // is the weight of that index, divided by total
  );
  return pdf;
}

vector<double> cdf_from_pdf(vector<double> pdf) {
  vector<double> cdf(pdf.size());
  partial_sum(
    pdf.begin(), pdf.end(),
    cdf.begin()
  );
  return cdf;
}

vector<double> complement(vector<double> ps) {
  vector<double> res(ps.size());
  transform(
    ps.begin(), ps.end(), // for all ps
    res.begin(), // put into res
    [](double p) { return 1-p; } // after taking complement
  );
  return res;
}

// probability of survival at an age, given survival to that age
static const vector<double> MOSQUITO_DAILY_SURVIVAL_PROBABILITY =
  complement(MOSQUITO_DAILY_DEATH_PROBABILITY);

vector<double> cumprod(vector<double> ps) {
  vector<double> res(ps.size());
  partial_sum(
    ps.begin(), ps.end(), // for all p in ps
    res.begin(), // put into res
    multiplies<double>()
    // the cumulative product at each p (cp_i = prod from 0 to i p_i)
  );
  return res;
}

// probability of surviving an age from birth
static const vector<double> MOSQUITO_SURVIVE_CUMPROB =
  cumprod(MOSQUITO_DAILY_SURVIVAL_PROBABILITY);

vector<double> relative_fraction(vector<double> ps) {
  vector<double> res(ps.size(), 1.0);
  copy(
    ps.begin(), ps.end()-1,
    // survival fraction == the fraction that survived previous day
    // so day 1 fraction = 1, day 2 fraction = day 1 survival prob from birth, etc
    res.begin()+1
  );
  return res;
}

// the relative fraction of mosquitos from a birth cohort alive at an age
// also, given no migration etc, the unnormalized steady-state age distribution
static const vector<double> MOSQUITO_AGE_RELFRAC =
  relative_fraction(MOSQUITO_SURVIVE_CUMPROB);

// the pdf of mosquito age
static const vector<double> MOSQUITO_AGE_PDF = norm_weights_to_pdf(MOSQUITO_AGE_RELFRAC);

// the cdf of mosquito age
static const vector<double> MOSQUITO_AGE_CDF = cdf_from_pdf(MOSQUITO_AGE_PDF);

vector<double> death_age_cdf(vector<double> survive_age_prob, vector<double> die_age_prob) {
  vector<double> pdf(die_age_prob);

  transform(
    survive_age_prob.begin(), survive_age_prob.end()-1,
    // using survival probs for previous day
    die_age_prob.begin()+1, // deaths probs for same day
    pdf.begin()+1, // overwriting from day 2 on
    multiplies<double>() // combine the probabilities
  );

  return cdf_from_pdf(pdf);
}


static const vector<double> MOSQUITO_DEATHAGE_CDF = death_age_cdf(MOSQUITO_SURVIVE_CUMPROB, MOSQUITO_DAILY_DEATH_PROBABILITY);

// e.g. today_biting_age_pdf = weight_biting_age_pdf(MOSQUITO_AGE_PDF, today_inf_p);

vector<double> weight_biting_age_pdf1( vector<double> raw_age_pdf, double prob_infecting_bite) {
    vector<double> pdf(raw_age_pdf);
    vector<double> biting_prob(raw_age_pdf.size(), prob_infecting_bite);
    biting_prob[1] = 1.0;

    partial_sum(
        biting_prob.begin(), biting_prob.end(),
        biting_prob.begin(), // just overwrite vals
        [](double prev, double p) { return prev*(1-p); }
        ); // TODO: this needs some checking and is more complex than necessary for
    // uniform probability.  however: sets up having historical mosquito FOI

    transform(
        pdf.begin(), pdf.end(),
        biting_prob.begin(),
        pdf.begin(),
        // alt could use a series of biting probabilities
        [](double age_p, double prev_non_inf_p){ return prev_non_inf_p * age_p; }
        );

    return pdf;
}

vector<double> weight_biting_age_pdf2(const double prob_infecting_bite) {
    vector<double> pdf(MOSQUITO_AGE_PDF.size(), 0.0); // prob of biting at age 0 is 0
    vector<double> biting_prob(MOSQUITO_AGE_PDF.size(), prob_infecting_bite);

    double denominator = 0.0;
    for (unsigned int age = 1; age < MOSQUITO_AGE_PDF.size(); ++age) {
        const double prob_inf_bite_given_age = pow(1.0 - prob_infecting_bite, age - 1) * prob_infecting_bite;
        pdf[age] = prob_inf_bite_given_age * MOSQUITO_AGE_PDF[age];
        denominator += pdf[age];
    }
   //cerr << denominator << " "; 
    for (unsigned int age = 1; age < MOSQUITO_AGE_PDF.size(); ++age) { pdf[age] /= denominator; }

    return pdf;
}

int main() {
    const int N = 101;
    const int NUM_MOS_AGES = 60;
    vector< vector<double> > biting_age_pdfs(N, vector<double>(NUM_MOS_AGES));
    assert(MOSQUITO_AGE_PDF.size() == NUM_MOS_AGES);
    for (unsigned int i = 0; i < biting_age_pdfs.size(); ++i) {
        biting_age_pdfs[i] = weight_biting_age_pdf2((double) i / (N-1));
        cerr << (double) i / (N-1);
        cerr << " " << accumulate( biting_age_pdfs[i].begin(), biting_age_pdfs[i].end(), 0.0 );
        //for (auto v: biting_age_pdfs[i]) cerr << " " << v;
        cerr << endl;
    }
    this_thread::sleep_for(chrono::milliseconds(100000));
}

