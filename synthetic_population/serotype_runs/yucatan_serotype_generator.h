#ifndef SERO_GENERATOR_H
#define SERO_GENERATOR_H
#include <iostream>
#include <vector>
#include <numeric>
#include <unistd.h>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//enum StreakState {GAP, RUN};
typedef float StreakState;
const float GAP = 0.0;
const float RUN = 1.0;

enum Transformations {NO_TRANSFORM, TRANSPOSE, TRANS_AND_NORM};

/*
fit values
tjhladish@dragonfly:~/work/dengue/synthetic_population/serotype_runs/serotype_refit$ 

sqlite3 denv1.sqlite 'select avg(geo_run), avg(geo_gap) from parameters p, jobs j where p.serial=j.serial and posterior > -1 and smcSet = 9'
13.8267781189999|2.695982357

sqlite3 denv2.sqlite 'select avg(geo_run), avg(geo_gap) from parameters p, jobs j where p.serial=j.serial and posterior > -1 and smcSet = 9'
9.53394009100002|8.31885120000003

sqlite3 denv3.sqlite 'select avg(geo_run), avg(geo_gap) from parameters p, jobs j where p.serial=j.serial and posterior > -1 and smcSet = 9'
3.08501077|4.91723286599999

sqlite3 denv4.sqlite 'select avg(geo_run), avg(geo_gap) from parameters p, jobs j where p.serial=j.serial and posterior > -1 and smcSet = 9'
2.54990213|8.49890623599998
*/

const std::vector< std::vector<StreakState> > fitting_period = {
//  1980's               1990's               2000's               2010's
{1, 1,1,1,1,1,1,0,1,1,1, 1,1,1,1,1,1,1,1,0,0, 0,0,1,0,0,1,1,1,1,1, 1,1,1,1},  // sero 1
{0, 0,0,0,0,0,0,1,0,0,0, 0,1,0,0,1,1,1,1,0,0, 0,1,1,0,1,1,1,1,1,1, 1,1,1,1},  // sero 2
{0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,1,1,0,0, 0,1,0,0,0,0,1,1,0,0, 0,0,0,1},  // sero 3
{0, 0,0,0,0,1,0,0,0,0,0, 0,0,0,0,1,0,1,1,0,0, 0,0,0,0,0,0,0,1,0,0, 0,0,1,1}}; // sero 4

void append_streak(const gsl_rng* RNG, std::vector<StreakState>& series, double prob, StreakState type, bool include_zeroes = false) {
    unsigned int streak = gsl_ran_geometric(RNG, prob) - (int) include_zeroes;
    series.resize(series.size() + streak, type);
}


std::vector< std::vector<StreakState> > generate_serotype_sequences(const gsl_rng* RNG, const int first_year, const int first_observed_year, const int last_year, bool use_exact_known_seros, Transformations xform = TRANS_AND_NORM) {
                                              // Runs        Gaps
    const std::vector< std::vector<double> > geo_means = {{13.8267781, 2.6959824},
                                                {9.5339401,  8.3188512},
                                                {3.0850108,  4.9172329},
                                                {2.5499021,  8.4989062}};
    const int NUM_SEROTYPES = 4;
    std::vector<bool> PRE_1979_CIRCULATION = {true, false, false, false};
    if (use_exact_known_seros) PRE_1979_CIRCULATION = {true, true, true, true}; // this isn't known, just a different model
    const std::vector<int> INTRO_YEARS = {0, 7, 17, 5}; // First year that each serotype appears, respectively, relative to first_observed_year

    const int SEQUENCE_LENGTH = last_year - first_year + 1;
    const int PREHISTORY_LENGTH = first_observed_year - first_year;
    //gsl_rng_set(RNG, time(NULL) * getpid()); // seed the rng using sys time and the process id

    const std::vector<double> p_run = {1/geo_means[0][0], 1/geo_means[1][0], 1/geo_means[2][0], 1/geo_means[3][0]};
    const std::vector<double> p_gap = {1/geo_means[0][1], 1/geo_means[1][1], 1/geo_means[2][1], 1/geo_means[3][1]};

    std::vector< std::vector<StreakState> > all_series(NUM_SEROTYPES);

    for (int i = 0; i < NUM_SEROTYPES; ++i) {
        std::vector<StreakState> series;
        if (PRE_1979_CIRCULATION[i]) {
            while ((signed) series.size() < PREHISTORY_LENGTH) {
                append_streak(RNG, series, p_gap[i], GAP);
                append_streak(RNG, series, p_run[i], RUN);
            }
            series.resize(PREHISTORY_LENGTH);
            reverse(series.begin(), series.end()); 
        } else {
            series.resize(PREHISTORY_LENGTH, GAP);
        }
        
        if (use_exact_known_seros) {
            series.insert(series.end(), fitting_period[i].begin(), fitting_period[i].end());
            bool include_zeroes = true;
            append_streak(RNG, series, p_run[i], RUN, include_zeroes);
        } else {
            series.resize(series.size() + INTRO_YEARS[i], GAP);
        }

        while ((signed) series.size() < SEQUENCE_LENGTH) {
            append_streak(RNG, series, p_run[i], RUN);
            append_streak(RNG, series, p_gap[i], GAP);
        }

        series.resize(SEQUENCE_LENGTH);
        all_series[i] = series;
    }
 
    std::vector< std::vector<StreakState> > all_series_t(SEQUENCE_LENGTH, std::vector<StreakState>(NUM_SEROTYPES));
    if (xform == TRANSPOSE or xform == TRANS_AND_NORM) {
        // transpose matrix, because that's what epidemic simulator is expecting
        for (unsigned int s = 0; s < all_series.size(); ++s) {
            for (unsigned int y = 0; y < all_series[s].size(); ++y) {
                all_series_t[y][s] = all_series[s][y]; 
            }
        }
    }

    if (xform == TRANS_AND_NORM) {
        for (unsigned int y = 0; y < all_series_t.size(); ++y) {
            float total = accumulate(all_series_t[y].begin(), all_series_t[y].end(), 0.0F);
            if (total > 0) {
                for (unsigned int s = 0; s < all_series_t[y].size(); ++s) all_series_t[y][s] /= total;
            }
        }
    }
 
    if (xform == TRANSPOSE or xform == TRANS_AND_NORM) {
        return all_series_t;
    } else {
        return all_series;
    }
}




#endif
