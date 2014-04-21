#!/usr/bin/perl
use strict;

my ($EF, $mos_move, $mos_cap, $daily_exp, $beta, $mp, $pm) = @ARGV;

my $betamp = $beta * $mp;
my $betapm = $beta * $pm;


my $command = "./model -randomseed 5500 \\
-locfile ./pop-yucatan/locations-yucatan.txt \\
-netfile ./pop-yucatan/network-yucatan.txt \\
-popfile ./pop-yucatan/population-yucatan.txt \\
-immfile ./pop-yucatan/immunity-yucatan/$EF.txt \\
-probfile ./pop-yucatan/swap_probabilities-yucatan.txt \\
-mosquitomove $mos_move \\
-mosquitomovemodel weighted \\
-mosquitoteleport 0.0 \\
-betapm $betapm \\
-betamp $betamp \\
-mosquitomultipliers 52 7 0.05 7 0.04 7 0.05 7 0.04 7 0.03 7 0.04 7 0.05 7 0.03 7 0.02 7 0.02 7 0.03 7 0.03 7 0.05 7 0.04 7 0.05 7 0.04 7 0.06 7 0.07 7 0.08 7 0.09 7 0.11 7 0.15 7 0.18 7 0.19 7 0.24 7 0.28 7 0.38 7 0.46 7 0.45 7 0.61 7 0.75 7 0.97 7 0.91 7 1.00 7 0.94 7 0.85 7 0.79 7 0.71 7 0.65 7 0.65 7 0.42 7 0.30 7 0.26 7 0.27 7 0.11 7 0.10 7 0.11 7 0.12 7 0.09 7 0.08 7 0.04 8 0.07 \\
-primarypathogenicity 1.0 0.25 1.0 0.25 \\
-secondaryscaling 1.0 1.0 1.0 1.0 \\
-mosquitocapacity $mos_cap \\
-daysimmune 730 \\
-runlength 470 \\
-dailyexposed $daily_exp $daily_exp $daily_exp $daily_exp \\
-ves 0.7\n";

print $command;

#./run_dengue.py 3 3.25 10 10
#<mean_epi_size> <sd_epi_size> <max_epi_size> <prob_epidemic>
#
#close(STDERR);

my $output = `$command`;
print "output: $output";
