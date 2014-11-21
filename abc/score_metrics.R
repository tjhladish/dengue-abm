prior_means = c('exp_fac'  = 35,
                'mos_mov'  = 0.5,
                'exp_coef' = -1,
                'num_mos'  = 62.5,
                'beta'     = 0.15)

empirical = c('mean'     = 101.956,
              'median'   = 30.6399,
              'stdev'    = 145.52,
              'max'      = 464.334,
              'skewness' = 1.31321,
              'mc'       = 0.242424)

#d0 = read.table("dengue_predictive_prior-varEIP-small.00", header=T)
#d1 = read.table("dengue_predictive_prior-varEIP-small.01", header=T)
#mean1 = apply(d1, 2, mean)

diff_sim = function(d, func) {
    sims = apply(d, 2, func)

    total_error = 0
    rms_error = 0

    for (i in 1:length(sims)) {
        n = names(sims[i])
        sim_val = as.numeric(sims[n])
        if (n %in% names(empirical)) {
            emp_val = as.numeric(empirical[n])
            error = 100*(sim_val - emp_val) / emp_val
            total_error = total_error + abs(error)/length(empirical)
            rms_error = rms_error + error**2/length(empirical)
            print(sprintf("%10s %10.3f %10.3f, %10.1f%%", n, sim_val, emp_val, error))
        } else if (n %in% names(prior_means)) {
            print(sprintf("%10s %10.3f %10.3f", n, sim_val, prior_means[n]))
        }
    }
    print(sprintf("Total mean error: %10.3f", total_error))
    print(sprintf("RMS error: %10.3f", sqrt(rms_error)))
}
