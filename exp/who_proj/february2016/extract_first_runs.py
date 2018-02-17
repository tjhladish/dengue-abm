#!/usr/bin/python
from sys import argv
from glob import glob

'''
0 begin 651dd824 1 56b2fd5f 10.601800 4.436440 0.600000 0.130000 0.257537 -1.684510 86.548500 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 1.000000 
SCENARIO 41155428 10.6018 4.43644 0.6 0.13 0.257537 -1.68451 86.5485 0 0 0 1 0 0 1
839659 people
651dd824 year: 0 9 22 31 18 0 
651dd824 year: 1 9 18 27 12 0 
651dd824 year: 2 14 131 145 75 2 
651dd824 year: 3 6 441 447 212 12 
651dd824 year: 4 9 87 96 37 2 
651dd824 year: 5 3 298 301 127 3 
651dd824 year: 6 11 128 139 69 2 
[thladish@dev1 february2016]$ tail raw_logs 
COHORT,842433542,29,1,0,2,7,9
COHORT,842433542,29,1,0,3,0,9
COHORT,842433542,29,1,1,0,9363,9
COHORT,842433542,29,1,1,1,433,9
COHORT,842433542,29,1,1,2,49,9
COHORT,842433542,29,1,1,3,3,9
8ed8872f year: 129 24 91936 91960 32991 5352 
WARNING: Daily output file already exists: /scratch/lfs/thladish/who-feb-2016/daily.-1898412241.842433542
WARNING: Aborting write.
0 end 8ed8872f 8814 20.7441 4.94111 0.6 0.13 0.271104 -0.521552 50.9905 4 1 0 1 0 0 1 102228 496 2273 71438 114005 44782 69714 91041 12987 68886 119395 29769 95245 29160 14094 126816 35592 16993 68103 70456 25244 60953 94604 42463 27875 88822 49488 32174 40232 91960 0.904605 
'''

for filename in glob('auto_output/burn80.err*'):
    #print filename
    buffer = ''
    good_record = True

    for line in file(filename):
        buffer = buffer + line
        p = line.strip().split()
        if len(p) > 1 and p[1] == 'end':
            if good_record:
                print buffer,
            else:
                good_record = True 
            buffer = ''
        elif len(p) > 0 and p[0] == 'WARNING:':
            good_record = False
        elif len(p) > 0 and p[0] == 'ERROR:':
            buffer = ''
            good_record = True 
