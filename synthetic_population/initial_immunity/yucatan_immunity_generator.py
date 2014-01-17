#!/usr/bin/python
from random import random
from collections import OrderedDict

NUM_SEROTYPES = 4
census = OrderedDict() # Will be 2D (year, age) OrderedDict containing (absolute) number of people in each age group in each year.
            # Last bin may be aggregated, i.e 85+ for everyone 85 and over.

# Reported DF + DHF cases, 1997-2011 (inclusive)
# Values are estimated from Fig. 1 of Hector's baseline studies manuscript
# http://arohatgi.info/WebPlotDigitizer/ was used to extract values

cases = [4220,	#1979
         4690,	#1980
         3360,	#1981
         1400,	#1982
          630,	#1983
         5490,	#1984
          190,	#1985
            0,	#1986
            0,	#1987
          330,	#1988
            0,	#1989
            0,	#1990
          330,	#1991
            0,	#1992
           20,	#1993
          650,	#1994
           50,	#1995
          610,	#1996
         5530,	#1997
           20,	#1998
           30,	#1999
            0,	#2000
          280,	#2001
          930,	#2002
            0,	#2003
           50,	#2004
          160,	#2005
          590,	#2006
         1840,	#2007
          700,	#2008
         3220,	#2009
         2520,	#2010
         6140]	#2011

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
                [0.00,	0.61,	0.39,	0.00],   # 2001 # Fig. 4
                [0.04,	0.96,	0.00,	0.00],   # 2002 # Fig. 4
                [0.04,	0.96,	0.00,	0.00],   # 2003 # extrapolated
                [0.04,	0.96,	0.00,	0.00],   # 2004 # extrapolated
                [0.11,	0.89,	0.00,	0.00],   # 2005 # Fig. 4
                [0.27,	0.54,	0.19,	0.00],   # 2006 # Fig. 4
                [0.91,	0.03,	0.04,	0.02],   # 2007 # Fig. 4
                [0.85,	0.15,	0.00,	0.00],   # 2008 # Fig. 4
                [0.45,	0.55,	0.00,	0.00],   # 2009 # Fig. 4
                [0.59,	0.41,	0.00,	0.00],   # 2010 # Fig. 4
                [0.32,	0.68,	0.00,	0.00] ]  # 2011 # Fig. 4

def import_census_data(filename):
    '''
year 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85+
1970 26222 25943 25663 25384 25105 24825 24546 24267 23433 22599 21765 20931 20098 19254 18411 17568 16725 15882 15239 14597 13955 13312 12670 12185 11700 11215 10730 10245 9827 9410 8993 8575 8158 8105 8051 7998 7945 7892 7529 7165 6802 6438 6075 5889 5702 5516 5330 5144 4864 4584 4305 4025 3746 3632 3519 3406 3293 3179 3120 3061 3002 2943 2884 2749 2613 2478 2343 2208 2073 1938 1804 1669 1534 1386 1238 1090 942 794 749 704 659 614 569 523 478 2623
1971 26661 26406 26150 25895 25639 25384 25128 24873 24070 23267 22464 21662 20859 20030 19201 18372 17542 16713 16041 15369 14697 14025 13353 12844 12335 11826 11317 10808 10377 9946 9514 9083 8652 8576 8499 8423 8346 8270 7888 7507 7125 6743 6361 6167 5972 5777 5583 5388 5102 4817 4531 4245 3960 3836 3712 3588 3464 3340 3274 3209 3143 3078 3012 2869 2727 2584 2441 2298 2158 2018 1878 1738 1598 1449 1300 1151 1002 852 803 754 705 656 607 558 509 2804
    '''
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


has_immunity = lambda s: int(s>0)


def recently_infected(states):
    if 1 in states or 2 in states:
        return True
    else:
        return False


def age_immunity(pop):
    for i in range(len(pop)):
        for s in range(NUM_SEROTYPES):
            if pop[i]['state'][s] == 1 or pop[i]['state'][s] == 2:
               pop[i]['state'][s] += 1
    return pop


#print cases
   
EXPANSION_FACTOR = 1

import_census_data('interpolated_ages-yucatan.out')
# read in population data
pop = []
header = True

print "Reading in population . . . "
for line in file('../../pop-yucatan/population-yucatan_final.txt'):
    if header:
        header = False
        continue
    p = line.split()
    pop.append({'pid':p[0], 'age':int(p[2]), 'state':[0,0,0,0]})

print "\t\t\tDone.\nFiltering by age . . . "
N     = len(pop)
NOW   = 2012
YEARS = range(1979,2012)

surviving_N = [[] for i in range(len(YEARS))] # number of people alive today who were also alive each year in 1979-2011

# years ago, starting in 1979
i = -1

kids_in_1987 = 0 # in the 8-14 range
# should be ~ 60% of the population
seropositive_kids_in_1987 = 0

for ya in range(len(YEARS),0,-1):
    i += 1
    pop = age_immunity(pop)
    for idx in range(len(pop)):
        person = pop[idx]
        if person['age'] >= ya:
            surviving_N[i].append(idx)
    surviving_N_i = len(surviving_N[i])
    total_N_i = sum(census[YEARS[i]].values()) 

    expected_num_of_cases = round(cases[i]*float(surviving_N_i)/total_N_i)
    expected_num_of_infections = expected_num_of_cases * EXPANSION_FACTOR
    print "year, years ago, total_pop, surviving_pop, total_cases, surviving_cases"
    print YEARS[i], ya, total_N_i, surviving_N_i, cases[i], expected_num_of_cases
    for s in range(4):
        print "\t" + str(s+1) + ":", serotype_wt[i][s] * expected_num_of_cases
        if serotype_wt[i][s] == 0.0:
            continue
        vulnerable_people = []
        for idx in surviving_N[i]:
            if pop[idx]['state'][s] == 0 and not recently_infected(pop[idx]['state']):
                vulnerable_people.append(idx)
        print "\t\tvulnerable ct:", len(vulnerable_people)
        hazard = expected_num_of_infections/len(vulnerable_people)
        for idx in vulnerable_people:
            if random() < hazard:
                pop[idx]['state'][s] = 1
    
    if YEARS[i] == 1987:
        for idx in surviving_N[i]:
            if pop[idx]['age']-ya >= 8 and pop[idx]['age']-ya <= 14:
                kids_in_1987 += 1.0
                seropos = sum(pop[idx]['state'])
                if seropos > 0:
                    seropositive_kids_in_1987 += 1.0
        print 'Number of kids 8-14 in 1987:', kids_in_1987
        print 'Number of seropositive kids 8-14 in 1987:', seropositive_kids_in_1987
        print 'Fraction seropositive:', seropositive_kids_in_1987/kids_in_1987
        print 'Empirical fraction seropositive: ~0.6'
        print 'Expansion factor:', 0.6/(seropositive_kids_in_1987/kids_in_1987)


print "\t\t\tDone."
print "Writing to disk . . ."
fo = open('immunity-yucatan.txt', 'w')
#fo.write('pid imm1 imm2 imm3 imm4\n')
fo.write('pid age imm1 imm2 imm3 imm4\n')
for person in pop:
    #fo.write(' '.join([person['pid']] + person['state']) + '\n')
    states = [str(has_immunity(i)) for i in person['state']]
    fo.write(' '.join([person['pid'], str(person['age'])] + states) + '\n')

fo.close()

print "\t\t\tDone."
