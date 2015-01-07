#!/usr/bin/python
from sys import exit

locations  = set()
houses     = dict()
workplaces = dict()
schools    = dict()

# process locations file that has already been pared down with something like
# awk '$3 >= -89.752178 && $3 <= -89.504299 && $4 >= 20.847583 && $4 <= 21.076652 { print $0 }' locations-yucatan.txt > locations-merida.txt

for line in file("./using_yucatan_location_ids/locations-merida.txt"):
    p = line.strip().split()
    locations.add(p[0])
    if p[1] == 'house':
        houses[p[0]] = -1
    elif p[1] == 'work':
        workplaces[p[0]] = -1
    elif p[1] == 'school':
        schools[p[0]] = -1
    else:
        print "Unknown location type:", line
        exit()


# filter network file
fo = open('./using_yucatan_location_ids/network-merida.txt','w')
for line in file("../pop-yucatan/network-yucatan.txt"):
    p = line.strip().split()
    if p[0] == 'locationID1':
        fo.write(line) # looking at header
    elif p[0] in locations and p[1] in locations:
        fo.write(line) # looking at an edge within merida
fo.close()

# counters
num_home    = 0
num_day_loc = 0
num_both    = 0
occupied_home   = set()
occupied_work   = set()
occupied_school = set()

# assign new sequential ids to people
# eventually these should all start at 0, but need to fix c++ code first
person_serial = 1

fo = open("./using_yucatan_location_ids/population-merida.txt", 'w')
fo3 = open("./using_yucatan_location_ids/yuc-mer_person_id_map.txt", 'w')

for line in file("../pop-yucatan/population-yucatan.txt"):
    p = line.strip().split()
    if p[0] == "pid":
        continue
    pid = p[0]
    hid = p[1]
    wid = p[6]
    home = False
    day_loc = False
    if hid in locations:
        home = True
        num_home += 1
    if wid in locations:
        day_loc = True
        num_day_loc += 1

    if home and day_loc: # this is a person who lives and works/schools/stays home in Merida -- keeping this
        occupied_home.add(hid)
        fo3.write(pid + ' ' + str(person_serial) + '\n')
        person_serial += 1
        num_both += 1
        fo.write(line)
        if wid in workplaces: # day location is a workplace
            occupied_work.add(wid)
        elif wid in schools:       # day location is a school
            occupied_school.add(wid)

fo.close()
fo3.close()

# assign new sequential ids to locations
# eventually these should all start at 0, but need to fix c++ code first
house_serial  = 1
work_serial   = len(occupied_home) + 1
school_serial = len(occupied_home) + len(occupied_work) + 1

fo2 = open("./using_yucatan_location_ids/yuc-mer_location_id_map.txt", 'w')
for hid in occupied_home:
    if houses[hid] == -1:
        houses[hid] = house_serial
        house_serial += 1
        fo2.write(hid + ' ' + str(houses[hid]) + '\n')
for wid in occupied_work:
    if workplaces[wid] == -1:
        workplaces[wid] = work_serial
        work_serial += 1
        fo2.write(wid + ' ' + str(workplaces[wid]) + '\n')
for sid in occupied_school:
    if schools[sid] == -1:
        schools[sid] = school_serial
        school_serial += 1
        fo2.write(sid + ' ' + str(schools[sid]) + '\n')

fo2.close()

print "home:", num_home, ",", float(num_both)/num_home
print "work/school/day loc:", num_day_loc, ",", float(num_both)/num_day_loc
print "both:", num_both





