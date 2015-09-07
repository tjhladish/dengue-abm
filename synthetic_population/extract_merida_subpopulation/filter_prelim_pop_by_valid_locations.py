#!/usr/bin/python

# build list of valid locations (i.e. those in Merida)
merida_locs = set()
for line in file('./using_yucatan_location_ids/locations-merida.txt'):
    p = line.strip().split()
    if p[0] == 'id':
        continue
    else:
        merida_locs.add(p[0])

# filter prelim pop file using valid households
fo = open('./using_yucatan_location_ids/population-merida_prelim.txt', 'w')
for line in file('./using_yucatan_location_ids/population-yucatan_prelim.txt'):
    p = line.strip().split()
    if p[0] == 'pid':
        fo.write(line)
        continue
    else:
        if p[1] in merida_locs:
            fo.write(line)
fo.close()

# filter network file using all valid locations
fo = open('./using_yucatan_location_ids/network-merida.txt', 'w')
for line in file('../../pop-yucatan/network-yucatan.txt'):
    p = line.strip().split()
    if p[0] == 'locationID1':
        fo.write(line)
        continue
    else:
        if p[0] in merida_locs and p[1] in merida_locs:
            fo.write(line)
fo.close()
