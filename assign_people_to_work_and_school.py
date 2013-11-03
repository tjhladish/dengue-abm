#!/usr/bin/python
from math import radians, cos, sin, asin, sqrt
from random import shuffle
from collections import defaultdict
from sys import exit

pixel_size = 0.00416667
min_x_center = -90.40409499
min_y_center = 19.72078911

# number of workplaces to look at when selecting a workplace
workplace_neighborhood = 1000 

# ratio for Mexico, according to World Bank
student_teacher_ratio = 28

# IPUMS values for EMPSTATD (detailed employment status) variable
work_codes   = [110, 112, 113, 116, 120]
home_codes   = [114, 200, 310, 320, 340, 390, 999]
school_codes = [111, 330]
child_code   = 0 # Used, inexplicably, for children 0-11
school_age   = 5 # children under this stay home


def haversine(lon1, lat1, lon2, lat2):
    '''
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    '''
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    km = 6371 * c
    return km 


def lookup_location_code(age, code):
    code = int(code)
    loc = ''
    if code in school_codes or (code == 0 and age >= school_age):
        loc = 's'
    elif code in work_codes:
        loc = 'w'
    elif code in home_codes or (code == 0 and age < school_age):
        loc = 'h'
    else:
        print "Error: encountered bad employment (EMPSTATD) code:", code
        exit()
    return loc

def binary_search(val_list, val, _lt, _gt, _eq, label, bound): 
#    print "val_list length:", len(val_list)
#    print "val, val_list[0][label], val_list[-1][label]:", val, val_list[0][label], val_list[-1][label]
    if _gt(val_list[0], val, label):
        return 0
    if _lt(val_list[-1], val, label):
        return len(val_list)

    imax = len(val_list)
    imin = 0 

    while (imin < imax):
        imid = imin + (imax-imin)/2
        if _lt(val_list[imid], val, label):
            imin = imid+1
        else:
            imax = imid

    if (imax == imin) and _eq(val_list[imin], val, label):
        return imin
    else:
        if bound == 'lower':
            return imin
        if bound == 'upper':
            return imin + 1

def get_nearby_places(pxi, pyi, loc_type, workplaces_and_schools, num_loc_needed):
    #print "searching for schools:", loc_type
    commute_range = -1
    positions_found = 0
    nearby_places = [] # by index in workplaces_and_schools
    if loc_type == 'w':
        pos_type = 'workers'
    else:
        pos_type = 'students'
    while (len(nearby_places) < num_loc_needed) or (positions_found <= 0):
        nearby_places = []
        commute_range += 1
        #print "\n\nEnvelope size:", 2*commute_range + 1, 'x', 2*commute_range + 1
        for x_val in range(pxi-commute_range, pxi+commute_range+1):
            pos_xmin = binary_search(workplaces_and_schools, x_val, dict_val_lt, dict_val_gt, dict_val_eq, 'xi', 'lower')
            pos_xmax = binary_search(workplaces_and_schools, x_val+1, dict_val_lt, dict_val_gt, dict_val_eq, 'xi', 'upper')
            if pos_xmin == pos_xmax:
                continue

            pos_ymin = binary_search(workplaces_and_schools[pos_xmin:pos_xmax], pyi-commute_range, dict_val_lt, dict_val_gt, dict_val_eq, 'yi', 'lower')
            pos_ymax = binary_search(workplaces_and_schools[pos_xmin:pos_xmax], pyi+commute_range+1, dict_val_lt, dict_val_gt, dict_val_eq, 'yi', 'upper')
            
            for i,w in enumerate(workplaces_and_schools[pos_xmin:pos_xmax][pos_ymin:pos_ymax]):
                raw_idx = pos_xmin + pos_ymin + i
                nearby_places.append(raw_idx)
                if pos_type == 'students':
                    positions_found += 1
                else:
                    positions_found += w[pos_type]

        #print "local positions:", positions_found            
        #print "workplaces found:", len(nearby_places) 
    # print "Envelope size:", 2*commute_range + 1, 'x', 2*commute_range + 1
    return nearby_places, positions_found

def select_nearest_school(px, py, nearby_places, workplaces_and_schools):
    closest_school, min_dist = nearby_places[0], () # () evaluates to greater than any number
    for idx in nearby_places:
        s = workplaces_and_schools[idx]
        d = haversine(px, py, s['x'], s['y'])
        if d < min_dist:
            min_dist = d
            closest_school = idx
    return closest_school, min_dist

def x_to_col_num(x):
    return int(round((x - min_x_center)/pixel_size))

def y_to_row_num(y):
    return int(round((y - min_y_center)/pixel_size))

def xy_cmp(a,b):
    if a['xi'] < b['xi']:
        return -1
    elif a['xi'] > b['xi']:
        return 1
    else:
        if a['yi'] < b['yi']:
            return -1
        elif a['yi'] > b['yi']:
            return 1
        else:
            return 0

def dict_val_lt(A,b,label):
    if A[label] < b:
        return True
    else:
        return False

def dict_val_gt(A,b,label):
    if A[label] > b:
        return True
    else:
        return False

def dict_val_eq(A,b,label):
    if A[label] == b:
        return True
    else:
        return False


hh_loc = dict()
header = True
print "reading locations"

# At this point, only includes houses
for line in file('locations-yucatan.txt'):
    '''
    id type x y x_ctr y_ctr
    1 house -89.6910220003 20.7015302009 -89.69159442 20.69995656
    2 house -89.700145483 20.6641877526 -89.69992776 20.66245653
    3 house -89.758249035 20.6954360352 -89.75826114 20.69578989
    4 house -89.6974011142 20.6551249287 -89.69576109 20.65412319
    '''

    if header:
        header = False
        continue

    p = line.split()
    hh_loc[p[0]] = {'x':float(p[2]), 'y':float(p[3])}

fo = file('population-yucatan_final.txt','w')

total_raw_size = {'w':0, 's':0}

header = True
workplaces_and_schools = []
print "reading workplaces & schools"
for line in file('schools_and_workplaces.out'):
    '''
    W 1 -89.6264173747 20.9599660422
    W 2 -89.6116996054 20.964832419
    W 1 -89.6428676041 20.9771890368
    W 2 -89.6405575542 20.9678584161
    W 1 -89.6255278043 20.9746128005
    '''

    if header:
        header = False
        continue

    p = line.split()
    w = {'type':p[0].lower(), 'raw_size':int(p[1]), 'x':float(p[2]), 'y':float(p[3]), 'workers':0, 'students':0}
    total_raw_size[w['type']] += w['raw_size']
    w['xi'] = x_to_col_num(w['x'])
    w['yi'] = y_to_row_num(w['y'])
    workplaces_and_schools.append(w)
    
workplaces_and_schools.sort(cmp=xy_cmp)

day_loc_ctr = {'h':0, 'w':0, 's':0} # home/work/school
pop_ids = []
pop = dict()
i = -1 

print "reading population"
header = True
for line in file('population-yucatan.txt'):
    i += 1
    #if i % 10000 == 0:
    #    print i
    '''
    pid hid age sex gridx gridy workid hh_serial pernum empstatd
    1 1 31 1 0 0 -1 2748179000 1 110
    2 1 29 2 0 0 -1 2748179000 2 110
    3 1 10 2 0 0 -1 2748179000 3 0
    4 2 32 1 0 0 -1 2748114000 1 110
    '''
    if header:
        fo.write(line)
        header = False
        continue

    p = line.split()
    hid = p[1]
    age = int(p[2])
    day_loc = lookup_location_code(age, p[9])
    day_loc_ctr[day_loc] += 1
    #pop[p[0]] = {'hid':hid, 'age':age, 'sex':p[3], 'x':hh_loc[hid]['x'], 'y':hh_loc[hid]['y'], 'day_loc':day_loc }
    # it doesn't seem like we're actually using all of these, and this program uses a lot of RAM
    pop[p[0]] = {'hid':hid, 'x':hh_loc[hid]['x'], 'y':hh_loc[hid]['y'], 'day_loc':day_loc, 'workid':-1}
    pop_ids.append(p[0])

print "Population size:", len(pop)
print
print "Total number of workers (IPUMS):", day_loc_ctr['w']
print "Total number of students (IPUMS):", day_loc_ctr['s']
print "Total number of homebodies (IPUMS):", day_loc_ctr['h']
print
print "Total number of non-teacher jobs (DENUE):", total_raw_size['w']
print "Total number of student + teacher positions (Min. of Ed.):", total_raw_size['s']
print "Student:Teacher ratio (WHO):", student_teacher_ratio

shuffle(pop_ids)    # UNCOMMENT AFTER DEBUGGING

# send kids to nearest school
# needs to happen first so we know how many teachers are needed
students_allocated = 0
for pid in pop_ids:
    loc_type = pop[pid]['day_loc']
    if loc_type != 's':
        continue
    num_loc_needed = 1

    px, py = pop[pid]['x'], pop[pid]['y']
    pxi, pyi = x_to_col_num(px), y_to_row_num(py)

    nearby_places, positions_found = get_nearby_places(pxi, pyi, loc_type, workplaces_and_schools, num_loc_needed)
    school, distance = select_nearest_school(px, py, nearby_places, workplaces_and_schools)
    # 'school' is an index in the workplaces_and_schools list; needs to be increased
    # by the number of households before output
    pop[pid]['workid'] = school 
    students_allocated += 1
    if students_allocated % 1000 == 0:
        print "students sent to school:", students_allocated
    #print px, py, school, distance

# normalize school sizes
'''
student_fraction = float(student_teacher_ratio) / (student_teacher_ratio + 1)
total_raw_students = total_raw_size['s'] * student_fraction
enrollment_rescaling_factor = day_loc_ctr['s'] / total_raw_students
total_teachers = 0

for place in workplaces_and_schools:
    if place['type'] == 's':
        currently_supported_students = place['raw_size'] * student_fraction
        required_enrollment = currently_supported_students * enrollment_rescaling_factor
        place['students'] = required_enrollment
        place['workers']  = required_enrollment / student_teacher_ratio
        total_teachers += place['workers']
'''

exit()

# normalize workplace sizes
total_jobs_still_needed = day_loc_ctr['w'] - total_teachers
employment_rescaling_factor = float(total_jobs_still_needed) / total_raw_size['w']
for place in workplaces_and_schools:
    if place['type'] == 'w':
        place['workers'] = place['raw_size'] * employment_rescaling_factor 
        place['students'] = 0



for pid in pop_ids[:10000]:
    loc_type = pop[pid]['day_loc']
    if loc_type == 'h':
        # I think we ultimately need to output their daytime location -- T0DO
        continue
    elif loc_type == 'w':
        num_loc_needed = 1000
    elif loc_type == 's':
        num_loc_needed = 1

    px, py = pop[pid]['x'], pop[pid]['y']
    pxi, pyi = x_to_col_num(px), y_to_row_num(py)

    nearby_places, positions_found = get_nearby_places(pxi, pyi, loc_type, workplaces_and_schools, num_loc_needed)
    print len(nearby_places), positions_found 
    
