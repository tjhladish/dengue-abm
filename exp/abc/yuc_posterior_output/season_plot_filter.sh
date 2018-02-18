awk '{if ($3 >= 42240 && $3 < 48080) print $0}' *.err > posterior_daily_1995-2011.out 
