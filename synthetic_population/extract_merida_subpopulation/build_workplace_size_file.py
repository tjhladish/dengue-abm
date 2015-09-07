#!/usr/bin/python

size_lookup = dict()
# build lookup of workplace sizes based on lat, long coordinates
for line in file('schools_and_workplaces-merida.out'):
    p = line.strip().split()
    if p[0] == 'w':
        # map lat, long to workplace size
        size_lookup[(p[2], p[3])] = p[1]

# try to map all workplaces in locations file to the size data we have
fo = open('./using_yucatan_location_ids/workplace_size_lookup.txt', 'w')
for line in file('./using_yucatan_location_ids/locations-merida.txt'):
    p = line.strip().split()
    if p[1] == 'work':
        if (p[2], p[3]) in size_lookup:
            fo.write(p[0] + ' ' + size_lookup[(p[2], p[3])] + '\n')
        else:
            print "Failed match:", line
fo.close()
