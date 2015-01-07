#!/usr/bin/python

locids = dict()
pids   = dict()

for line in file('./using_yucatan_location_ids/yuc-mer_location_id_map.txt'):
    p = line.strip().split()
    locids[p[0]] = p[1]

for line in file('./using_yucatan_location_ids/yuc-mer_person_id_map.txt'):
    p = line.strip().split()
    pids[p[0]] = p[1]


# process population file
fo = open('population-merida.txt', 'w')
for line in file('./using_yucatan_location_ids/population-merida.txt'):
    p = line.strip().split()
    p[0] = pids[p[0]]
    p[1] = locids[p[1]]
    p[6] = locids[p[6]]
    fo.write(' '.join(p) + '\n')

fo.close()

# process location file
# we also remove locations that no longer have people, because they all either
# worked outside of merida, or exclusively had employees/students who lived
# outside of merida
fo = open('locations-merida.txt', 'w')
for line in file('./using_yucatan_location_ids/locations-merida.txt'):
    p = line.strip().split()
    if p[0] in locids:
        p[0] = locids[p[0]]
        fo.write(' '.join(p) + '\n')

fo.close()
