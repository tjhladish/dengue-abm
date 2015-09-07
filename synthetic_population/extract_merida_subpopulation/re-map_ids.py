#!/usr/bin/python

loc_id_map = dict()
person_ctr = 0

first_loc_id = 1
first_person_id = 1

# re-map location IDs starting at 1
fo  = open('locations-merida.txt', 'w')
fo2 = open('yucatan-to-merida_loc_id_map.txt','w')
for line in file('./using_yucatan_location_ids/locations-merida.txt'):
    p = line.strip().split()
    if p[0] == 'id': # header
        fo.write(line)
        continue
    else:
        new_id = str(len(loc_id_map) + first_loc_id)
        fo2.write(p[0] + ' ' + new_id + '\n')
        loc_id_map[p[0]] = new_id
        p[0] = new_id # overwriting makes joining on the next line simpler
        fo.write(' '.join(p) + '\n')
fo.close()
fo2.close()

# update network file using new IDs
fo = open('network-merida.txt', 'w')
for line in file('./using_yucatan_location_ids/network-merida.txt'):
    p = line.strip().split()
    if p[0] == 'locationID1':
        fo.write(line)
        continue
    else:
        fo.write(loc_id_map[p[0]] + ' ' + loc_id_map[p[1]] + '\n')
fo.close()

# re-map people and update using new location IDs
fo  = open('population-merida.txt', 'w')
fo2 = open('yucatan-to-merida_population_id_map.txt','w')
for line in file('./using_yucatan_location_ids/population-merida.txt'):
    '''
    pid hid age sex hh_serial pernum workid
    5919 1199 55 1 2748540000 1 394452
    5920 1199 45 2 2748540000 2 1199
    5921 1199 74 2 2748540000 3 1199
    5922 1199 26 1 2748540000 4 393957
    '''
    p = line.strip().split()
    if p[0] == 'pid': # header
        fo.write(line)
        continue
    else:
        new_id = str(person_ctr + first_person_id)
        fo2.write(p[0] + ' ' + new_id + '\n')
        person_ctr += 1
        p[0] = new_id
        p[1] = loc_id_map[p[1]]
        p[6] = loc_id_map[p[6]]
        fo.write(' '.join(p) + '\n')
fo.close()
fo2.close()
