#!/usr/bin/python
from random import random

# Reported DF + DHF cases, 1997-2011 (inclusive)
# Values are estimated from Fig. 1 of Hectors baseline studies manuscript
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

#print cases
   
EXPANSION_FACTOR = 12 

# read in population data
pop = []
header = True

print "Reading in population . . . "
for line in file('population-yucatan_final.txt'):
    if header:
        header = False
        continue
    p = line.split()
    pop.append({'pid':p[0], 'age':int(p[2]), 'state':['0','0','0','0']})

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
    for idx in range(len(pop)):
        person = pop[idx]
        if person['age'] >= ya:
            surviving_N[i].append(idx)
    N_then = len(surviving_N[i])
    expected_num_of_cases = round(cases[i]*float(N_then)/N)
    print YEARS[i], ya, N_then, "cases (# among current pop):", cases[i], expected_num_of_cases
    already_infected = {} # already infected THIS YEAR
    for s in range(4):
        print "\t" + str(s+1) + ":", serotype_wt[i][s] * expected_num_of_cases
        if serotype_wt[i][s] == 0.0:
            continue
        vulnerable_people = []
        for idx in surviving_N[i]:
            if pop[idx]['state'][s] == '0' and idx not in already_infected:
                vulnerable_people.append(idx)
        print "\t\tvulnerable ct:", len(vulnerable_people)
        hazard = EXPANSION_FACTOR*expected_num_of_cases/len(vulnerable_people)
        for idx in vulnerable_people:
            if random() < hazard:
                pop[idx]['state'][s] = '1'
    
    if YEARS[i] == 1987:
        for idx in surviving_N[i]:
            if pop[idx]['age']-ya >= 8 and pop[idx]['age']-ya <= 14:
                kids_in_1987 += 1.0
                seropos = sum(map(int,pop[idx]['state']))
                if seropos > 0:
                    seropositive_kids_in_1987 += 1.0
        print 'Number of kids 8-14 in 1987:', kids_in_1987
        print 'Number of seropositive kids 8-14 in 1987:', seropositive_kids_in_1987
        print 'Fraction seropositive:', seropositive_kids_in_1987/kids_in_1987
        print 'Empirical fraction seropositive: ~0.6'
        print 'Expansion factor:', 0.6/(seropositive_kids_in_1987/kids_in_1987)

print "\t\t\tDone."
print "Writing to disk . . ."
fo = open('immunity-yucatan_prelim.txt', 'w')
#fo.write('pid imm1 imm2 imm3 imm4\n')
fo.write('pid age imm1 imm2 imm3 imm4\n')
for person in pop:
    #fo.write(' '.join([person['pid']] + person['state']) + '\n')
    fo.write(' '.join([person['pid'], str(person['age'])] + person['state']) + '\n')

fo.close()

print "\t\t\tDone."
