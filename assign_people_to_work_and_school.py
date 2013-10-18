#!/usr/bin/python
from random import shuffle
from sys import exit

pixel_size = 0.00416667
min_x_center = -90.40409499
min_y_center = 19.72078911

# ratio for Mexico, according to World Bank
student_teacher_ratio = 28

# IPUMS values for EMPSTATD (detailed employment status) variable
work_codes   = [110, 112, 113, 116, 120]
home_codes   = [114, 200, 310, 320, 340, 390, 999]
school_codes = [111, 330]
child_code   = 0 # Used, inexplicably, for children 0-11
school_age   = 5 # children under this stay home

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

def binary_search(val_list, val, _lt, _gt, _eq, label): 
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
        return -1

def x_to_col_num(x):
    return round(x - min_x_center/pixel_size)

def y_to_row_num(y):
    return round((y - min_y_center)/pixel_size)

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
header = True
workplaces = []
#print "reading workplaces"

total_raw_size = {'w':0, 's':0}
for line in file('schools_and_workplaces.out'):
    '''
    W 1 -89.6264173747 20.9599660422
    W 2 -89.6116996054 20.964832419
    W 1 -89.6428676041 20.9771890368
    W 2 -89.6405575542 20.9678584161
    W 1 -89.6255278043 20.9746128005
    '''

    p = line.split()
    w = {'type':p[0].lower(), 'raw_size':int(p[1]), 'x':float(p[2]), 'y':float(p[3]), 'size':0}
    total_raw_size[w['type']] += w['raw_size']
    w['xi'] = x_to_col_num(w['x'])
    w['yi'] = y_to_row_num(w['y'])
    workplaces.append(w)
    
workplaces.sort(cmp=xy_cmp)

day_loc_ctr = {'h':0, 'w':0, 's':0} # home/work/school
pop_ids = []
pop = dict()
i = -1 

print "reading population"
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
    pop[p[0]] = {'hid':hid, 'age':age, 'sex':p[3], 'x':hh_loc[hid]['x'], 'y':hh_loc[hid]['y'], 'day_loc':day_loc }
    pop_ids.append(p[0])

# normalize workplace sizes
total_raw_jobs = 0
for w in workplaces:
    total_raw_jobs += w['raw_size']

print "Population size:", len(pop)
print
print "Total number of workers:", day_loc_ctr['w']
print "Total number of students:", day_loc_ctr['s']
print "Total number of homebodies:", day_loc_ctr['h']
print
print "Total number of non-teacher jobs:", total_raw_size['w']
print "Total number of student + teacher positions:", total_raw_size['s']
print "Student:Teacher ratio:", student_teacher_ratio

# number of classrooms plus extra staff
num_teachers = (total_raw_size['s'] / student_teacher_ratio) + total_raw_size['s'] % student_teacher_ratio
num_students = total_raw_size['s'] - num_teachers
print "Total number of teachers:", num_teachers
print "Total number of students:", num_students

exit()
for w in workplaces:
    w['size'] = float(len(pop)) * w['raw_size'] / total_raw_jobs

shuffle(pop_ids)

for pid in pop_ids:
    px, py = pop[pid]['x'], pop[pid]['y']
    pxi, pyi = int(px/pixel_size), int(py/pixel_size)
    commute_range = -1
    nearby_workplaces = []
    
    while len(nearby_workplaces) < 1000:
        commute_range += 1
        xmin = pxi-commute_range
        xmax = pxi+commute_range+1 # range is [xmin, xmax), thus +1
        pos_xmin = binary_search(workplaces, xmin, dict_val_lt, dict_val_gt, dict_val_eq, 'xi')
        pos_xmax = binary_search(workplaces, xmax, dict_val_lt, dict_val_gt, dict_val_eq, 'xi')
        #print "Person loc:", pxi, pyi, workplaces[pos_xmin]['xi'], workplaces[pos_xmax]['xi']
        #for w in workplaces[pos_xmin:pos_xmax]:
        #    print w

