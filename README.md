# dengue
---
This stochastic, discrete-time, agent-based model of dengue transmission explicitly represents a synthetic population of humans and mosquitoes with which transmission of all four Dengue serotypes can be simulated.

## Daily Transmission Cycle
This ABM is structured around a daily cycle of human-to-mosquito and mosquito-to-human transmission.
Uninfected mosqutioes can acquire dengue after biting infectious humans and can further spread the virus to susceptible humans.
Mosquito infection follows an SEI model of infection progressing between **S**usceptible, **E**xposed, and **I**nfectious.
Human infection similarly follows an SEIR model of infection progressing between **S**usceptible, **E**xposed, **I**nfectious, and **R**esistant. Moreover, human infections can develop varying levels of severity: asymptomatic (**I<sub>A</sub>**), symptomatic (**I<sub>S</sub>**), and withdrawn (**I<sub>W</sub>**).
The model also uses fitted seasonal rain and temperature patterns to reliably determine the synthetic mosquito population size and the extrinsic incubation period (EIP) for mosquito dengue infections.

## Human Synthetic Population
The human synthetic population for this ABM is assumed to have a constant size and age structure in order to reduce computational complexity. Individuals have fixed age, gender, and household assignment.
Within the ABM's daily simulation cycle, humans move between their designated home and a daytime location (either a workplace or a school).
Instead of modeling an aging population, immune histories are annually shifted from younger to older people in order to capture the accumulation of immunity over time which is vital for modeling long-term transmission dynamics.

## Mosquito Synthetic Population
The overall mosquito population size is calculated from location-specific sizes that vary based on seasonal rainfall patterns. Upon infection, a mosquito will be assigned an current age, EIP, and age of death sanmpled from respective distributions.
Mosquito movement can occur between nearby locations weighted by the inverse square distance to a specific location from the origin.

## Instructions For Use
Brief instructions for the dengue transmission model code version 1.0.
by Dennis Chao and Thomas Hladish (2012)

We have released the source code for our dengue model (described in
the PLoS Neglected Tropical Diseases paper "Controlling dengue with
vaccines in Thailand") under the GPLv3 (see the files LICENSE and
gpl.txt).

This bundle includes C++ source code and data files, all placed in the
"denguemodelcode/" directory. There is no user interface, and the
documentation is minimal and poor. There is a Makefile that works in
Linux, and it produces an executable called "model". Run it on the
command line with appropriate flags to specify model parameters, the
locations of the necessary data files, and the desired output file names.
The model sends a list of all exposed and infectious mosquitoes as well
as all infected people to stdout.

### Data files for the Bangphae model:

- locations-bangphae.txt: a list of all locations in Bangphae
- population-bangphae.txt: a list of all people in the synthetic population for Bangphae
- immunity-bangphae.txt: prior exposure to each of the four serotypes for each of the people in population-bangphae.txt
- network-bangphae.txt: a list of all "adjacent" locations corresponding to ids listed in locations-bangphae.txt. This is used for mosquito movement.

### Command-line options:

  - `randomseed [seed]`: supply a random number seed to the GSL generator
  - `runlength [days]`: length of the simulation run in days. run the model for 364 days for a one-year simulation, otherwise the population will get shuffled on day 365.
  - `initialinfected [num]`: number of randomly selected individuals to infect before the simulation
  - `initialexposed [num]`: number of randomly selected individuals to expose (and possibly infect) before the simulation
  - `dailyexposed [num]`: number of randomly selected individuals to expose (and possibly infect) at the start of each day
  - `vaccinatephased [n] [y] [a] [f] [y] [a] [f]...`: used to vaccinate different age groups in different years. the first argument is the number of vaccination groups desired. this is followed by triplets specifying: the year the group should be vaccinated (occurs on January 1, and the first year is 0), the age cohort to be vaccinated (in years), and the fraction of this cohort to vaccinate (from 0.0 to 1.0).
  - `primarypathogenicity [f1] [f2] [f3] [f4]`: fraction of primary infections that result in illness, by serotype
  - `secondaryscaling [f1] [f2] [f3] [f4]`: fraction of non-primary infections (relative to the first) that result in illness, by serotype. the actual pathogenicity fraction for secondary infections is therefore the product of the primary (see the "primarypathogenicity" argument) and secondary fraction
  - `betapm [x]`: the probability that a mosquito is infected when biting an infectious person
  - `betamp [x]`: the probability that a person is infected when bitten by infectious mosquito
  - `mosquitomove [p]`: daily probability of mosquito movement to an adjacent location
  - `mosquitomovemodel [s]`: mosquitoes move to any neighbor with equal probability if "uniform" or weighted by inverse distance squared if "weighted"
  - `mosquitoteleport [p]`: daily probability of mosquito "teleportation" (to anywhere in the synthetic population)
  - `mosquitocapacity [n]`: mean number of mosquitoes per location
  - `mosquitodistribution [s]`: distribution of mosquitos per location. Set to "constant" for all locations to have the same number of mosquitoes or "exponential" for the number to be exponentially distributed.
  - `mosquitomultipliers [n] [d] [f] [d] [f]...`: relative number of mosquitoes for seasonality. the first argument is the number of pairs of numbers coming up. each pair consists of an integer that specifies a number of days followed by a floating point number that is a multiplier for the mosquito capacity to set the number of mosquitoes per location for this number of days. the number of days should sum to 365, unless you are trying to be funny and make dengue season fall out of sync with the calendar year.
  - `externalincubations [n] [d1] [d2] [d3] [d4]...`: external incubation periods. the first argument is the number of pairs of numbers coming up. each pair consists of an integer that specifies a number of days followed by an integer that is the external incubation period for this number of days. the number of days should sum to 365.
  - `daysimmune`: number of days after recovery that a person has perfect cross-protective immunity to all other serotypes
  - `VES [n]`: reduction in susceptibility of vaccinees, assuming all-or-none protection (0.0-1.0)
  - `VESs [n1] [n2] [n3] [n4]`: reduction in susceptibility (0.0-1.0) of vaccinees to each of 4 serotypes
  - `VESsnaive [n1] [n2] [n3] [n4]`: reduction in susceptibility (0.0-1.0) of vaccinees who had no prior exposure to any serotype. Only works for all-or-none vaccines.
  - `VEI [n]`: reduction in infectiousness of vaccinees (0.0-1.0)
  - `VEP [n]`: reduction in probability of becoming ill upon infection of vaccinees (0.0-1.0)
  - `vaccineleaky`: Setting this flag makes VES leaky. Default is all-or-none.
  - `prevaccinate [f]`: fraction of the population to pre-vaccinate
  - `prevaccinateage [n] [a1] [a2] [f] [a1] [a2] [f]...`: prevaccinate by age group. the first argument is the number of age groups, followed by triplets specifying the groups. fraction `[f]` of those who are ages `[a1]` to `[a2]` (in years) are pre-vaccinated.
  - `nosecondary`: no secondary transmission allowed. this is used for R0 estimation
  - `maxinfectionparity [n]`: specifies the maximum number of serotypes that can (serially) infect a single individual. default is 4.
  - `popfile [filename]`: location of the input file that contains the synthetic population
  - `immfile [filename]`: location of the input file that contains the prior immunity information for the synthetic population
  - `locfile [filename]`: location of the input file that contains the locations for the model (i.e., houses, classrooms, workplaces)
  - `netfile [filename]`: location of the input file that lists every pair of adjacent locations corresponding to the information in "locfile"
  - `probfile [filename]`: location of the (optional) input file that contains information for swapping immune statuses at the end of each year
  - `peoplefile [filename]`: specifies the name of the output file that will contain the information for every infection in a simulation run
  - `yearlypeoplefile [filename]`: specifies the filename prefix of the output file that will contain the information for every infection each year in a simulation run. the output filenames will have the year and ".csv" appended (e.g., filename5.csv)
  - `dailyfile [filename]`: specifies the name of the output file that will contain the number of people infected and symptomatic each day by serotype

### Instructions:

A reasonable way to run the model of Bangphae from the command line is:

`./model -locfile locations-bangphae.txt -netfile network-bangphae.txt -popfile population-bangphae.txt -immfile immunity-bangphae.txt -mosquitomove 0.15 -mosquitoteleport 0.01 -betapm 0.1 -betamp 0.25 -mosquitomultipliers 12 31 0.11 28 0.09 31 0.20 30 0.30 31 0.83 30 1.00 31 0.50 31 0.37 30 0.38 31 0.26 30 0.30 31 0.19 -primarypathogenicity 1.0 0.25 1.0 0.25 -secondaryscaling 1.0 1.0 1.0 1.0 -mosquitocapacity 42 -runlength 364 -dailyexposed 2 2 2 2 -daysimmune 365  -prevaccinate 0.30 -ves 0.7 -peoplefile people-output-bangphaeseasonal-vac30-ves70.csv -dailyfile daily-output-bangphaeseasonal-vac30-ves70.csv > output-bangphaeseasonal-vac30-ves70.csv`

And an example for a 10-year simulation:

`./model -randomseed 5489 -locfile locations-bangphae.txt -netfile network-bangphae.txt -popfile population-bangphae.txt  -immfile immunity-bangphae.txt -mosquitomove 0.15 -mosquitoteleport 0.01 -betapm 0.1 -betamp 0.25 -mosquitomultipliers 12 31 0.11 28 0.09 31 0.20 30 0.30 31 0.83 30 1.00 31 0.50 31 0.37 30 0.38 31 0.26 30 0.30 31 0.19 -primarypathogenicity 1.0 0.25 1.0 0.25 -secondaryscaling 1.0 1.0 1.0 1.0 -mosquitocapacity 42 -daysimmune 120 -runlength 3650 -dailyexposed 2 2 2 2 -ves 0.7 -yearlypeoplefile people-output-bangphae-multiseason-randomseed5489-y -dailyfile daily-output-bangphae-multiseason-randomseed5489.csv > output-bangphae-multiseason-randomseed5489.csv`
