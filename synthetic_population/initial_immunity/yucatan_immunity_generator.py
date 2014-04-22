#!/usr/bin/python
from random import random, shuffle, sample, choice, randint
from collections import OrderedDict
from sys import exit
from copy import deepcopy
from sys import argv, exit

if len(argv) != 3:
    print "\n\tUsage: ./yucatan_immunity_generator.py <expansion_factor> <last_year>" 
    print "\tNB: <last_year> is the last year of data to include when simulating immunity.\n"
    exit()

NUM_SEROTYPES = 4
MAX_CENSUS_AGE = 85 # used if a census age group has only a minimum value, e.g. '85+'

# Reported DF + DHF cases, 1997-2011 (inclusive)

cases = [4234,	#1979
         4672,	#1980
         3377,	#1981
         1412,	#1982
          643,	#1983
         5495,	#1984
          193,	#1985
           34,	#1986
           15,	#1987
          356,	#1988
            2,	#1989
            8,	#1990
          352,	#1991
           22,	#1992
           29,	#1993
          680,	#1994
           69,	#1995
          650,	#1996
         5529,	#1997
           36,	#1998
           43,	#1999
            0,	#2000
          287,	#2001
          946,	#2002
           26,	#2003
           57,	#2004
          162,	#2005
          627,	#2006
         1861,	#2007
          721,	#2008
         3212,	#2009
         2517,	#2010
         6132,	#2011
         5705]	#2012
 
serotype_wt = [ [1.00,	0.00,	0.00,	0.00],   # 1979 # Fig. 1
                [1.00,	0.00,	0.00,	0.00],   # 1980 # Fig. 1
                [1.00,	0.00,	0.00,	0.00],   # 1981 # Fig. 1
                [1.00,	0.00,	0.00,	0.00],   # 1982 # Fig. 1
                [1.00,	0.00,	0.00,	0.00],   # 1983 # Fig. 1
                [0.50,	0.00,	0.00,	0.50],   # 1984 # Fig. 1
                [1.00,	0.00,	0.00,	0.00],   # 1985 # Fig. 1
                [0.00,	1.00,	0.00,	0.00],   # 1986 # Fig. 1
                [1.00,	0.00,	0.00,	0.00],   # 1987 # Fig. 1
                [1.00,	0.00,	0.00,	0.00],   # 1988 # Fig. 1
                [1.00,	0.00,	0.00,	0.00],   # 1989 # Fig. 1
                [1.00,	0.00,	0.00,	0.00],   # 1990 # Fig. 1
                [0.50,	0.50,	0.00,	0.00],   # 1991 # Fig. 1
                [1.00,	0.00,	0.00,	0.00],   # 1992 # Fig. 1
                [1.00,	0.00,	0.00,	0.00],   # 1993 # Fig. 1
                [0.33,	0.33,	0.00,	0.34],   # 1994 # Fig. 1
                [0.50,	0.50,	0.00,	0.00],   # 1995 # Fig. 1
                [0.25,	0.25,	0.25,	0.25],   # 1996 # Fig. 1
                [0.09,	0.00,	0.87,	0.04],   # 1997 # Fig. 4
                [0.09,	0.00,	0.87,	0.04],   # 1998 # extrapolated
                [0.09,	0.00,	0.87,	0.04],   # 1999 # extrapolated
                [0.09,	0.00,	0.87,	0.04],   # 2000 # extrapolated
                [0.00,	0.60,	0.40,	0.00],   # 2001 # Fig. 4
                [0.04,	0.96,	0.00,	0.00],   # 2002 # Fig. 4
                [0.04,	0.96,	0.00,	0.00],   # 2003 # extrapolated
                [0.04,	0.96,	0.00,	0.00],   # 2004 # extrapolated
                [0.11,	0.89,	0.00,	0.00],   # 2005 # Fig. 4
                [0.27,	0.55,	0.18,	0.00],   # 2006 # Fig. 4
                [0.90,	0.04,	0.04,	0.02],   # 2007 # Fig. 4
                [0.85,	0.15,	0.00,	0.00],   # 2008 # Fig. 4
                [0.46,	0.54,	0.00,	0.00],   # 2009 # Fig. 4
                [0.59,	0.41,	0.00,	0.00],   # 2010 # Fig. 4
                [0.34,	0.66,	0.00,	0.00] ]  # 2011 # Fig. 4

def sample_serotype(wts):
    r = random()
    for i, w in enumerate(wts):
        if r < w:
            return i
        else:
            r -= w
    return len(wts) - 1

def import_census_data(filename):
    '''
year 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85+
1970 26222 25943 25663 25384 25105 24825 24546 24267 23433 22599 21765 20931 20098 19254 18411 17568 16725 15882 15239 14597 13955 13312 12670 12185 11700 11215 10730 10245 9827 9410 8993 8575 8158 8105 8051 7998 7945 7892 7529 7165 6802 6438 6075 5889 5702 5516 5330 5144 4864 4584 4305 4025 3746 3632 3519 3406 3293 3179 3120 3061 3002 2943 2884 2749 2613 2478 2343 2208 2073 1938 1804 1669 1534 1386 1238 1090 942 794 749 704 659 614 569 523 478 2623
1971 26661 26406 26150 25895 25639 25384 25128 24873 24070 23267 22464 21662 20859 20030 19201 18372 17542 16713 16041 15369 14697 14025 13353 12844 12335 11826 11317 10808 10377 9946 9514 9083 8652 8576 8499 8423 8346 8270 7888 7507 7125 6743 6361 6167 5972 5777 5583 5388 5102 4817 4531 4245 3960 3836 3712 3588 3464 3340 3274 3209 3143 3078 3012 2869 2727 2584 2441 2298 2158 2018 1878 1738 1598 1449 1300 1151 1002 852 803 754 705 656 607 558 509 2804
    '''
    census = OrderedDict() # Will be 2D (year, age) OrderedDict containing (absolute) number 
                           # of people in each age group in each year.
                           # Last bin may be aggregated, i.e 85+ for everyone 85 and over.
    ages = []
    for line in file(filename):
        p = line.strip().split()
        if p[0] == 'year':
            ages = p[1:]
            continue
        else:
            year = int(p[0])
            counts = map(int, p[1:])
            #print year, counts
            census[year] = OrderedDict(zip(ages, counts))
    return census


has_immunity = lambda s: int(s>0)


def recently_infected(states):
    # "recently infected" means in the last one or two years
    if 1 in states or 2 in states:
        return True
    else:
        return False


def age_immunity(pop):
    #print "Before aging:"
    #z = 0
    #for i in range(len(pop[10])):
    #    if sum(pop[10][i]) > 0:
    #        print i, pop[10][i]
    #        z += 1
    #        if z == 10:
    #            break

    for age in range(len(pop)):
        for i in range(len(pop[age])):
            #if sum(pop[age][i]) > 0:
            #    print "\nBefore aging:"
            #    print i, age, pop[age][i], id(pop[age][i]) 
            for s in range(NUM_SEROTYPES):
                if pop[age][i][s] >= 1:
                   pop[age][i][s] += 1
            #if sum(pop[age][i]) > 0:
            #    print "After aging:"
            #    print i, age, pop[age][i], id(pop[age][i]) 

    return pop


def age_str_to_int(age_str):
    try:
        age = int(age_str)
    except ValueError:
        age = MAX_CENSUS_AGE
    return age


def initialize_full_population(census, first_year):
    full_pop = [[] for i in range(MAX_CENSUS_AGE + 1)]
    for age_str in census[first_year]:
        age = age_str_to_int(age_str)
        for i in range(census[first_year][age_str]):
            full_pop[age].append([0,0,0,0])
        #print first_year, age, census[first_year][age_str]
    return full_pop

def age_full_population(census, full_pop, new_year):
    '''
    "Age" the population by sampling from younger age classes as necessary to make up
    the census numbers for the new year, e.g. if there are ten 30-year-olds in 1980,
    and we need eight 31-year-olds in 1981, we will randomly select 8 of the 10.  If
    we need MORE than the number available, we will sample with replacement to get the
    required number.

    Age 0 individuals are special because they don't have a younger age class to
    inherit from.  MAX_CENSUS_AGE individuals are likely special because they may be
    in a binned age class with no specified upper bound.  This is tricky because 85+
    individuals in 1980 who survive to 1981 are still 85+, but some of the age 84
    individuals from 1980 have joined their ranks.  Because of this, for the eldest
    age class we sample from the two eldest groups.  This does not account for
    differing mortality rates among the group, but since it's a relatively small
    fraction of the population that is likely to have relatively homogenous immunity
    (due to their long exposure to epidemics), this assumption is probably not
    consequential.

    Note: This will fail if there is a year with 0 individuals in an age class other
    than the last.
    '''
    full_pop = age_immunity(full_pop)
    new_pop = []
    age_keys = census[new_year].keys()
    age_keys.reverse() # census is an OrderedDict
    for age_str in age_keys: # going from oldest to youngest
        age = age_str_to_int(age_str)
        if age == 0:
            # nothing to inherit for infants
            full_pop[age] = [[0,0,0,0] for i in range(census[new_year][age_str])]
            continue
        sampling_pop = full_pop[age - 1]
        if age == MAX_CENSUS_AGE: # We are assuming here that the last age class is a multi-year bin, eg. 85+
            sampling_pop += full_pop[MAX_CENSUS_AGE]

        if len(sampling_pop) >= census[new_year][age_str]:
            full_pop[age] = sample(sampling_pop, census[new_year][age_str])
        else:
            #print "year, age, sample size, census size: ", new_year, age, len(sampling_pop), census[new_year][age_str]
            full_pop[age] = [deepcopy(choice(sampling_pop)) for i in xrange(census[new_year][age_str])]

    return full_pop
   
EXPANSION_FACTOR = int(argv[1])
LAST_YEAR = int(argv[2])

YEARS = range(1979,LAST_YEAR + 1)
first_year = YEARS[0]

census = import_census_data('interpolated_ages-yucatan.out')

full_pop = initialize_full_population(census, first_year)

# year number, 1979 = 0
year = 0 

kids_in_1987 = 0 # in the 8-14 range
# should be ~ 60% of the population
seropositive_kids_in_1987 = 0

for ya in range(len(YEARS),0,-1):
    if year > 0:
        full_pop = age_full_population(census, full_pop, YEARS[year])  # delete after debugging
    num_of_infections = cases[year] * EXPANSION_FACTOR
    pop_size_by_age = [len(full_pop[age]) for age in range(len(full_pop))]
    pop_size = sum(pop_size_by_age)
    print "year, years ago, pop, cases, infections, expansion factor"
    print YEARS[year], ya, pop_size, cases[year], num_of_infections, EXPANSION_FACTOR
    print "expected serotype distribution: ", serotype_wt[year]
    serotype_tally = [0 for i in range(NUM_SEROTYPES)]
    for i in range(num_of_infections):
        s = sample_serotype(serotype_wt[year])
        infection_occurred = False
        while not infection_occurred: 
            test = randint(0, pop_size - 1) # randint includes endpoints
            age = 0
            while test >= pop_size_by_age[age]:
                test -= pop_size_by_age[age]
                age += 1
            #if sum(full_pop[age][test]) > 0:
            #    print "age, state: ", age, full_pop[age][test]
            if full_pop[age][test][s] == 0 and not recently_infected(full_pop[age][test]):
                # infect test person
                full_pop[age][test][s] = 1
                serotype_tally[s] += 1
                infection_occurred = True

    print "actual serotype counts: ", serotype_tally
    if num_of_infections > 0:
        print "actual serotype frequency: ", [float(c)/sum(serotype_tally) for c in serotype_tally]
    else:
        print "actual serotype frequency: [NA, NA, NA, NA]"
    
    if YEARS[year] == 1987:
        for age in range(8, 15):
            for kid in full_pop[age]:
                kids_in_1987 += 1.0
                seropos = sum(kid)
                if seropos > 0:
                    seropositive_kids_in_1987 += 1.0
        print 'Number of kids 8-14 in 1987:', kids_in_1987
        print 'Number of seropositive kids 8-14 in 1987:', seropositive_kids_in_1987
        print 'Fraction seropositive:', seropositive_kids_in_1987/kids_in_1987
        print 'Empirical fraction seropositive: ~0.6'
        print 'Expansion factor:', 0.6/(seropositive_kids_in_1987/kids_in_1987)

    year += 1
    print

#for age in range(len(full_pop)):
#    for person in full_pop[age]:
#        print age, ' '.join(map(str, person))
#
#exit()
print "\t\t\tDone."

# read in population data
header = True
fo = open('/work/01856/thladish/initial_immunity/' + str(EXPANSION_FACTOR) + '.txt', 'w')
#fo = open(str(EXPANSION_FACTOR) + '.txt', 'w')
fo.write('pid age imm1 imm2 imm3 imm4\n')

print "Building immunity file . . . "
counter = 0
for line in file('../../pop-yucatan/population-yucatan.txt'):
    if header:
        header = False
        continue
    if counter % 100000 == 0:
        print counter
    counter += 1
    p = line.split()
    pid = p[0]
    age = p[2]
    
    age_int = int(age)
    age_int = age_int if age_int <= MAX_CENSUS_AGE else MAX_CENSUS_AGE
    states = choice(full_pop[age_int])
    fo.write(' '.join([pid, age] + map(str,states)) + '\n')

fo.close()

print "\t\t\tDone."
