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

0 begin 1db5ee93 617476 1 11.737700 4.572130 0.537741 0.117897 0.292393 -1.146230 38.385300 0.163296 0.000000 0.000000 9.000000 30.000000 0.000000 0.000000 0.000000 0.000000 
1db5ee93 617476 day: 0 0 0 0 0 0 7.60738 5.69055
1db5ee93 617476 day: 1 0 0 0 0 0 7.62665 5.6519
1db5ee93 617476 day: 2 0 0 0 0 0 7.64217 5.61759
1db5ee93 617476 day: 3 0 0 0 0 0 7.6182 5.58799
1db5ee93 617476 day: 4 0 0 0 0 0 7.53833 5.56345
1db5ee93 617476 day: 5 0 0 0 0 0 7.51761 5.54437
1db5ee93 617476 day: 6 0 0 0 0 0 7.44435 5.53116
1db5ee93 617476 day: 7 0 0 0 0 0 7.31339 5.52423
1db5ee93 617476 day: 8 0 0 0 0 0 7.1134 5.52398
'''

seen = set()
day_min = 116*365 - 100 # Jan 1, 1995
day_max = 133*365 - 100 # Jan 1, 2012

ctr = 0
for filename in glob('./vac_merida*'):
    ctr += 1
    if ctr % 10 == 0:
        print ctr
    tag = ''
    ef_mild = 1
    ef_severe = 1
    fo = None
    for line in file(filename):
        p = line.strip().split()
        # is this the beginning of a no vaccine, no vector control log?
        if (p[1] == 'begin' and line.strip()[-72:] == '0.000000 0.000000 9.000000 30.000000 0.000000 0.000000 0.000000 0.000000'):
            copy_log = True
            tag = p[2]
            if (tag not in seen):
                seen.add(tag)
                fo = open('./asdf/' + tag + '.err', 'w')
                # NB: the cases column was previously, and incorrectly, called mild
                fo.write(' '.join(['tag', 'serial', 'day', 'intros', 'local', 'infec', 'cases', 'severe', 'eip', 'nmos', 'ef_mild', 'ef_severe']) + '\n')

            else:
                break
            ef_mild = p[5] 
            ef_severe = p[6] 
        elif (len(p) > 2 and p[2] == 'day:' and int(p[3]) >= day_min and int(p[3]) < day_max):
            fo.write( ' '.join(p[0:2] + p[3:] + [ef_mild, ef_severe]) + '\n' )
        elif (len(p) > 1 and p[1] == 'end' and fo != None):
            fo.close()

#for filename in glob('junk.out*'):
#    fo = open(filename)
#    p = fo.readline().strip().split()
#    fo.close()
#    tag = p[0]
#    rename(filename, tag + '.out')
