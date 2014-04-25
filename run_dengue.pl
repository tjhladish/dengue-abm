#!/usr/bin/perl
use strict;

my ($EF, $mos_move, $daily_exp, $betamp, $betapm) = @ARGV;

my $dengue_dir = "/scratch/lfs/thladish/dengue";

my $command = "$dengue_dir/mpi_model -randomseed 5500 \\
-locfile $dengue_dir/pop-yucatan/locations-yucatan.txt \\
-netfile $dengue_dir/pop-yucatan/network-yucatan.txt \\
-popfile $dengue_dir/pop-yucatan/population-yucatan.txt \\
-immfile $dengue_dir/pop-yucatan/immunity/$EF.txt \\
-expansionfactor $EF \\
-probfile $dengue_dir/pop-yucatan/swap_probabilities-yucatan.txt \\
-mosquitomove $mos_move \\
-mosquitomovemodel weighted \\
-mosquitoteleport 0.0 \\
-betapm $betapm \\
-betamp $betamp \\
-mosquitomultipliers 52 7 0.05 7 0.04 7 0.05 7 0.04 7 0.03 7 0.04 7 0.05 7 0.03 7 0.02 7 0.02 7 0.03 7 0.03 7 0.05 7 0.04 7 0.05 7 0.04 7 0.06 7 0.07 7 0.08 7 0.09 7 0.11 7 0.15 7 0.18 7 0.19 7 0.24 7 0.28 7 0.38 7 0.46 7 0.45 7 0.61 7 0.75 7 0.97 7 0.91 7 1.00 7 0.94 7 0.85 7 0.79 7 0.71 7 0.65 7 0.65 7 0.42 7 0.30 7 0.26 7 0.27 7 0.11 7 0.10 7 0.11 7 0.12 7 0.09 7 0.08 7 0.04 8 0.07 \\
-primarypathogenicity 1.0 0.25 1.0 0.25 \\
-secondaryscaling 1.0 1.0 1.0 1.0 \\
-mosquitocapacity 100 \\
-daysimmune 730 \\
-runlength 4500 \\
-dailyexposed $daily_exp $daily_exp $daily_exp 0.0 \\
-ves 0.7\n";

#print $command; # debugging

#my $output_str = `$command`;
exec $command;
