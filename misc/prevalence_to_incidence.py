#!/usr/bin/python
from sys import argv

'''
time,type,id,location,serotype,symptomatic,withdrawn,new_infection
3,p,246402,48932,3,1,0,1
4,p,246402,48932,3,1,0,0
5,p,1424384,301949,4,0,0,1
5,p,246402,48932,3,1,0,0
6,p,1424384,301949,4,0,0,0
6,p,246402,48932,3,1,0,0
7,p,1424384,301949,4,0,0,0
7,p,652516,136580,3,1,0,1
7,p,326446,65402,4,1,0,1
'''

header = True
#seen = set()

for line in file(argv[1]):
    if header:
        print line,
        header = False
    else:
        p = line.strip().split(',')
        if p[1] == 'p' and p[7]=='1':
        #if p[1] == 'p' and (p[2], p[4]) not in seen:
            print line,
#            seen.add((p[2], p[4]))
