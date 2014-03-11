#!/usr/bin/python
from sys import argv

'''
time,type,id,location,serotype,symptomatic,withdrawn
0,p,1337222,284866,2,0,0
0,p,1276049,272183,3,0,0
0,p,1076962,231003,4,1,0
1,p,1670965,348400,3,1,0
1,p,1495506,315605,2,1,0
'''

header = True
seen = set()

for line in file(argv[1]):
    if header:
        print line,
        header = False
    else:
        p = line.strip().split(',')
        if p[1] == 'p' and (p[2], p[4]) not in seen:
            print line,
            seen.add((p[2], p[4]))
