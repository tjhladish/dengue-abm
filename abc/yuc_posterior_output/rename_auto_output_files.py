#!/usr/bin/python
from glob import glob
from os import rename

'''
[thladish@dev1 auto_output]$ head junk.err-54
0 begin f1441f0f 0 52.364300 2.704760 0.095186 -1.688130 73.382400 0.224527 
f1441f0f day: 0 0 0 0 0 0 18 10.9318
f1441f0f day: 1 0 0 0 0 0 18 10.8576
f1441f0f day: 2 0 0 0 0 0 18 10.7917
f1441f0f day: 3 0 0 0 0 0 18 10.7348
f1441f0f day: 4 0 0 0 0 0 18 10.6877
f1441f0f day: 5 0 0 0 0 0 17 10.651
f1441f0f day: 6 0 0 0 0 0 17 10.6257
f1441f0f day: 7 1 0 1 0 0 17 10.6123
f1441f0f day: 8 1 0 1 0 0 16 10.6118

[thladish@dev1 auto_output]$ head junk.out-54
f1441f0f T: 364 annual: 7 8925 8932 1749 196 
f1441f0f T: 729 annual: 8 168028 168036 36857 4272 
f1441f0f T: 1094 annual: 6 257380 257386 55101 6344 
f1441f0f T: 1459 annual: 2 301455 301457 63978 7346 
'''

seen = set()
day_min = 116*365 - 100 # Jan 1, 1995
day_max = 133*365 - 100 # Jan 1, 2012

for filename in glob('../auto_output/abc_yucatan-w-daily.err-*'):
    tag = ''
    ef_mild = 1
    ef_severe = 1
    fo = None
    day_buffer = []
    year_buffer = []
    for line in file(filename):
        p = line.strip().split()
        if (p[1] == 'begin'):
            tag = p[2] + '.' + p[3]
            ef_mild = p[5]   # these used to be fields 4 and 5; now 5 and 6
            ef_severe = p[6]
            day_buffer.append(' '.join(['tag', 'day', 'intros', 'local', 'infec', 'cases', 'severe', 'eip', 'nmos', 'ef_mild', 'ef_severe']) + '\n')
            year_buffer.append(' '.join(['tag', 'year', 'intros', 'local', 'infec', 'cases', 'severe', 'ef_mild', 'ef_severe']) + '\n')
        elif (p[2] == 'day:' and int(p[3]) >= day_min and int(p[3]) < day_max):
            #fo.write( ' '.join([p[0]] + p[2:] + [ef_mild, ef_severe]) + '\n' )
            day_buffer.append(' '.join([p[0]] + p[2:] + [ef_mild, ef_severe]) + '\n')
        elif p[2] == 'year:':
            year_buffer.append(' '.join([p[0]] + p[2:] + [ef_mild, ef_severe]) + '\n')
        elif (p[1] == 'end'):
            if (tag not in seen):
                print "writing", tag
                seen.add(tag)
                fo = open(tag + '.daily', 'w')
                fo.writelines(day_buffer)
                fo.close()
                fo = open(tag + '.annual', 'w')
                fo.writelines(year_buffer)
                fo.close()
            day_buffer = []
            year_buffer = []

#for filename in glob('../auto_output/abc_yucatan-w-daily.out-*'):
#    fo = open(filename)
#    p = fo.readline().strip().split()
#    fo.close()
#    tag = p[0]
#    rename(filename, tag + '.out')
