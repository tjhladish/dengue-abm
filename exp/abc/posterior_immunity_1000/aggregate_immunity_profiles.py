#!/usr/bin/python
from sys import argv, exit

# we're looking to compare the immunity profiles output from running the simulator
# to what has been recently found in Yucatan.  We record a positive or negative infection
# history, and then bin in 5-year groups, with everyone over 60 binned together.

'''
pid age imm1 imm2 imm3 imm4
1 31 0 0 -9354 0 
2 29 0 0 0 -5239 
'''

if len(argv) < 3:
    print "\n\tUsage: ./aggregate_immunity_profiles.py <output_filename> <input_filename[s]>\n\n"
    exit()


# grab IDs for people in merida
merida_ids = set()
#for line in file('../../pop-merida/using_yucatan_location_ids/population-merida.txt'): # if yucatan was simulated
for line in file('../../pop-merida/population-merida.txt'): # if merida was simulated
    if line[0:3] == 'pid':
        continue
    p = line.strip().split()
    merida_ids.add(int(p[0]))

output_filename = argv[1]
fo = open(output_filename, 'w')
fo.write('rep age_class pos neg total\n')
input_filenames = argv[2:]

# bins: 0-4, 5-9, ... 50-59, 60+
bin_upper_bounds = [4,9,14,19,29,39,49,59,()] # () is a hackish idiom for infinity
nclass = 9 # number of age classes
rep = -1   # file counter 

for filename in input_filenames:
    print "processing file ", filename
    rep += 1
    tally = [{} for i in range(nclass)]
    for age_class in range(len(tally)):
        tally[age_class] = {'pos':0, 'neg':0, 'total':0}

    header = True
    for line in file(filename):
        if header:
            header = False
            continue
        p = map(int, line.strip().split())
        ID = p[0]
        if ID not in merida_ids:                # KEEP ONLY MERIDA RESIDENTS
            continue
        age = p[1]
        age_class = 0
        while age > bin_upper_bounds[age_class]:
            age_class += 1
        tally[age_class]['total'] += 1
        if p[2]!=0 or p[3]!=0 or p[4]!=0 or p[5]!=0: # any infections history
            tally[age_class]['pos'] += 1
        else:
            tally[age_class]['neg'] += 1

    
    for age_class in range(len(tally)):
        fo.write(' '.join(map(str, [rep, age_class, tally[age_class]['pos'], tally[age_class]['neg'], tally[age_class]['total']])) + '\n')

fo.close()


