#!/usr/bin/python
from math import radians, cos, sin, asin, sqrt
from random import random, shuffle
from collections import defaultdict, OrderedDict
from sys import exit
from copy import deepcopy

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

# This is an optimization to cut down RAM by >50% and maintain
# speed by using lists with an index lookup rather than a dictionary
# to represent each person
field_idx = dict(zip('hid age sex x y workid hh_serial pernum day_loc'.split(), range(9)))

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
    while (len(nearby_places) < num_loc_needed and len(nearby_places) < len(workplaces_and_schools)) or (positions_found <= 0):
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
                #raw_idx = pos_xmin + pos_ymin + i
                if pos_type == 'students':
                    positions_found += 1
                    nearby_places.append(w['workid'])
                elif w[pos_type] > 0:
                        positions_found += w[pos_type]
                        nearby_places.append(w['workid'])

        #print "local positions:", positions_found            
        #print "workplaces found:", len(nearby_places) 
    # print "Envelope size:", 2*commute_range + 1, 'x', 2*commute_range + 1
    return nearby_places, positions_found

def select_nearest_school(px, py, nearby_places, workplaces_and_schools, workplace_lookup):
    # nearby_places is a list of indeces for workplaces_and_schools
    closest_school, min_dist = workplaces_and_schools[workplace_lookup[nearby_places[0]]], () # () evaluates to greater than any number
    for workid in nearby_places:
        s = workplaces_and_schools[workplace_lookup[workid]]
        d = haversine(px, py, s['x'], s['y'])
        if d < min_dist:
            min_dist = d
            closest_school = s
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

def import_workplace_sizes(filename):
    size_lookup = dict()
    total_raw_workers = 0
    for line in file(filename):
        p = line.strip().split()
        loc_id, size = int(p[0]), int(p[1])
        size_lookup[loc_id] = size
        total_raw_workers += size
    return size_lookup, total_raw_workers

def import_locations(locations_filename, workplace_size_filename, hh_loc, workplaces_and_schools):
    size_lookup, total_raw_workers = import_workplace_sizes(workplace_size_filename)
    workplace_lookup = dict()

    header = True
    print "reading merida locations from preliminary locations file"

    # At this point, only includes houses
    for line in file(locations_filename):
        '''
        id type x y x_ctr y_ctr
        1199 house -89.5207175949 20.8953195532 -89.52076095 20.89579005
        1203 house -89.5127969723 20.8952596464 -89.51242761 20.89579005
        1204 house -89.5312553802 20.9054693985 -89.53326096 20.90412339
        1216 house -89.5369291537 20.8734076091 -89.53742763 20.8749567
        1217 house -89.5047143269 20.8843110157 -89.50409427 20.88329004
        '''

        if header:
            header = False
            continue

        p = line.split()
        loc_type = p[1].lower()
        loc_id = int(p[0])
        if loc_type == 'house':
            hh_loc[loc_id] = {'x':float(p[2]), 'y':float(p[3])}
        else: # a workplace or school
            workers = 0
            w = {'workid':loc_id, 'type':loc_type, 'x':float(p[2]), 'y':float(p[3]), 'workers':0, 'students':0}
            if loc_type == 'work':
                w['raw_workers'] = size_lookup[loc_id]
            else:
                w['raw_workers'] = 0
            w['xi'] = x_to_col_num(w['x'])
            w['yi'] = y_to_row_num(w['y'])
            workplaces_and_schools.append(w)
            workplace_lookup[w['workid']] = len(workplaces_and_schools) - 1 # look up workplace idx by id

    workplaces_and_schools.sort(cmp=xy_cmp)
    return total_raw_workers, workplace_lookup


#def import_workplaces_and_schools(filename, workplaces_and_schools, current_max_loc_id):
#    header = True
#    print "reading workplaces & schools"
#    loc_id = current_max_loc_id + 1
#    total_raw_workers = 0
#    workplace_lookup = dict()
#    for line in file(filename):
#        '''
#        W 1 -89.6264173747 20.9599660422
#        W 2 -89.6116996054 20.964832419
#        W 1 -89.6428676041 20.9771890368
#        W 2 -89.6405575542 20.9678584161
#        W 1 -89.6255278043 20.9746128005
#        '''
#
#        if header:
#            header = False
#            continue
#
#        p = line.split()
#        w = {'workid':loc_id, 'type':p[0].lower(), 'x':float(p[2]), 'y':float(p[3]), 'workers':0, 'students':0}
#        if w['type'] == 'w':
#            w['raw_workers'] = int(p[1])
#            total_raw_workers += w['raw_workers']
#        else:
#            w['raw_workers'] = 0
#            
#        w['xi'] = x_to_col_num(w['x'])
#        w['yi'] = y_to_row_num(w['y'])
#        workplaces_and_schools.append(w)
#        workplace_lookup[w['workid']] = len(workplaces_and_schools) - 1
#        loc_id += 1
#        
#    workplaces_and_schools.sort(cmp=xy_cmp)
#    return total_raw_workers, workplace_lookup

def import_population(filename, pop, pop_ids, day_loc_ctr):
    ''' This function imports a preliminary version of the population file
    that doesn't have people assigned to workplaces yet (perhaps among
    other things.'''

    print "reading population"
    header = True
    i = -1 
    for line in file(filename):
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
#            fo.write(line)
            header = False
            continue

        p = map(int, line.split())

        hid = p[1]
        age = p[2]
        day_loc = lookup_location_code(age, p[9])
        day_loc_ctr[day_loc] += 1
        #           hid, age, sex, x, y, workid, hh_serial, pernum, day_loc
        pop[p[0]] = p[1:4] + [hh_loc[hid]['x'], hh_loc[hid]['y'], -1, p[7], p[8], day_loc]
        #pop[p[0]] = {'hid':hid, 'age':age, 'sex':p[3], 'x':hh_loc[hid]['x'], 'y':hh_loc[hid]['y'], 'day_loc':day_loc, 'workid':-1}
        pop_ids.append(p[0])

    print "Population size:", len(pop)
    print
    print "Total number of workers (IPUMS):", day_loc_ctr['w']
    print "Total number of students (IPUMS):", day_loc_ctr['s']
    print "Total number of homebodies (IPUMS):", day_loc_ctr['h']
    print
    return

def send_kids_to_school(pop, pop_ids, workplaces_and_schools, total_raw_workers, workplace_lookup):
    students_allocated = 0
    for pid in pop_ids:
        loc_type = pop[pid][field_idx['day_loc']]
        if loc_type != 's':
            continue
        num_loc_needed = 1

        px, py = pop[pid][field_idx['x']], pop[pid][field_idx['y']]
        pxi, pyi = x_to_col_num(px), y_to_row_num(py)

        nearby_places, positions_found = get_nearby_places(pxi, pyi, loc_type, workplaces_and_schools, num_loc_needed)
        school, distance = select_nearest_school(px, py, nearby_places, workplaces_and_schools, workplace_lookup)
        # 'school' is an element in the workplaces_and_schools list
        pop[pid][field_idx['workid']] = school['workid'] 
        school['students'] += 1
        school['raw_workers'] += 1.0/student_teacher_ratio
        total_raw_workers     += 1.0/student_teacher_ratio
        students_allocated += 1
        if students_allocated % 10000 == 0:
            print "students sent to school:", students_allocated
        #print px, py, school, distance
    return total_raw_workers

def choose_workplace(px, py, nearby_places, workplaces_and_schools, workplace_lookup):
    raw_weights = [0.0 for i in range(len(nearby_places))]
    for i, workid in enumerate(nearby_places):
        w = workplaces_and_schools[workplace_lookup[workid]]
        dist = haversine(px, py, w['x'], w['y'])
        size = w['workers']
        raw_weights[i] = size / dist**2

    # normalize weights
    probs = []
    total = sum(raw_weights)
    for wt in raw_weights:
        probs.append(wt/total)
    
    r = random()
    for i,p in enumerate(probs):
        if r < p:
            return workplaces_and_schools[workplace_lookup[nearby_places[i]]]
        r -= p

    return workplaces_and_schools[workplace_lookup[nearby_places[-1]]]


# Import location data
hh_loc = dict()
workplaces_and_schools = []
locations_filename = 'using_yucatan_location_ids/locations-merida.txt'
workplace_sizes_filename = 'using_yucatan_location_ids/workplace_size_lookup.txt'
total_raw_workers, workplace_lookup = import_locations(locations_filename, workplace_sizes_filename, hh_loc, workplaces_and_schools)

# Import workplace and school location data
#
# We are using workplace size (# employees) as a weight, but ignoring
# school size data currently, as we feel it is more realistic to send
# students to the nearest school.
# max_loc_id = max(hh_loc.keys())
# total_raw_workers, workplace_lookup = import_workplaces_and_schools('schools_and_workplaces.out', workplaces_and_schools, max_loc_id)

# Import population data
pop = OrderedDict()
pop_ids = []
day_loc_ctr = {'h':0, 'w':0, 's':0} # home / work / school
import_population('using_yucatan_location_ids/population-merida_prelim.txt', pop, pop_ids, day_loc_ctr)

print "Total number of non-teacher jobs (DENUE):", total_raw_workers
print "Student:Teacher ratio (World Bank):", student_teacher_ratio
print "Total number of needed teachers:", day_loc_ctr['s']/student_teacher_ratio
print "Total raw number of jobs (DENUE + needed teachers):", total_raw_workers + day_loc_ctr['s']/student_teacher_ratio

shuffle(pop_ids)

# Send kids to nearest school
# Needs to happen first so we know how many teachers are needed
total_raw_workers = send_kids_to_school(pop, pop_ids, workplaces_and_schools, total_raw_workers, workplace_lookup)

# normalize workplace sizes
employment_rescaling_factor = day_loc_ctr['w'] / float(total_raw_workers)
for place in workplaces_and_schools:
    place['workers'] = place['raw_workers'] * employment_rescaling_factor 

# Filehandle for file we're going to write
fo = file('using_yucatan_location_ids/population-merida.txt','w')
fo.write('pid hid age sex hh_serial pernum workid\n')

# Make a copy so we can delete places from the original data structure as they fill up
W_AND_S_COPY = deepcopy(workplaces_and_schools)
ctr = 0 

# For each person
for pid in pop.keys():
    person = pop[pid]
    loc_type = person[field_idx['day_loc']]
    # Have them stay at home if they don't work or go to school
    if loc_type == 'h':
        person[field_idx['workid']] = person[field_idx['hid']] # person stays home
    # If they work, probabilistically choose a workplace based on where they live
    # and how many positions are available at each workplace
    elif loc_type == 'w':
        px, py = person[field_idx['x']], person[field_idx['y']]
        pxi, pyi = x_to_col_num(px), y_to_row_num(py)

        nearby_places, positions_found = get_nearby_places(pxi, pyi, loc_type, workplaces_and_schools, workplace_neighborhood)
        #print len(nearby_places), positions_found 
        workplace = choose_workplace(px, py, nearby_places, W_AND_S_COPY, workplace_lookup)
        workplace['workers'] -= 1 # remove one available job
        # If the selected workplace no longer has openings,
        # remove it from the list, so we don't have to consider it again
        if workplace['workers'] <= 0:
            for i,v in enumerate(workplaces_and_schools):
                if v['workid'] == workplace['workid']:
                    del workplaces_and_schools[i]
                    break
        person[field_idx['workid']] = workplace['workid'] # assign worker

    # Students already have the "workid" (prob should be called day_loc_id to avoid
    # confusion), so we don't have to do much for them, just output their info
    fo.write(' '.join(map(str,[
                       pid, 
                       person[field_idx['hid']],  
                       person[field_idx['age']],  
                       person[field_idx['sex']], 
                       person[field_idx['hh_serial']],
                       person[field_idx['pernum']],
                       person[field_idx['workid']],
                      ])) + '\n') 
    ctr += 1
    if ctr % 1000 == 0:
        print "placed", ctr, "people"
        print "workplace list size:", len(workplaces_and_schools)
    '''
    pid hid age sex hh_serial pernum workid 
    1 1 31 1 0 0 -1 2748179000 1 110
    2 1 29 2 0 0 -1 2748179000 2 110
    3 1 10 2 0 0 -1 2748179000 3 0
    4 2 32 1 0 0 -1 2748114000 1 110
    '''
fo.close()
