#define _USE_MATH_DEFINES

#include <math.h>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctype.h>
#include <assert.h>
#include <algorithm>

using namespace std;

double pixel_size   = 0.00416667;
double min_x_center = -90.40409499;
double min_y_center = 19.72078911;

// number of workplaces to look at when selecting a workplace
const int workplace_neighborhood = 1000;

// ratio for Mexico, according to World Bank
const int student_teacher_ratio = 28;

// IPUMS values for EMPSTATD (detailed employment status) variable
const set<int> work_codes   = {110, 112, 113, 116, 120};
const set<int> home_codes   = {114, 200, 310, 320, 340, 390, 999};
const set<int> school_codes = {111, 330};
const int child_code   = 0; // Used, inexplicably, for children 0-11
const int school_age   = 5; // children under this stay home

// translated from python; probably still speeds things up
/*map<string, int> field_idx = {{"hid",       0},
                              {"age",       1},
                              {"sex",       2},
                              {"x",         3},
                              {"y",         4},
                              {"workid",    5},
                              {"hh_serial", 6},
                              {"pernum",    7},
                              {"day_loc",   8}};*/ 

struct HouseType {
    int hid;
    double x;
    double y;
};

struct LocationType {
    int workid;
    char type;
    double x; // longitude
    double y; // latitude
    map<string,int> pixel = {{"xi",0}, {"yi",0}};
    //int xi;
    //int yi;
    int raw_workers;
    int workers;
    int students;
};

struct PersonType {
    //    pid hid age sex gridx gridy workid hh_serial pernum empstatd
    //    1 1 31 1 0 0 -1 2748179000 1 110
    int pid;
    int hid;
    int age;
    int sex;
    double x;
    double y;
    int workid;
    long int hh_serial;
    int pernum;
    char day_loc;
    int empstat;
};

double deg_to_rad(double degree) { return M_PIl*degree/180; }
double rad_to_deg(double radians) { return 180*radians/M_PIl; }
int x_to_col_num(double x) { return (int) round((x - min_x_center)/pixel_size); }
int y_to_row_num(double y) { return (int) round((y - min_y_center)/pixel_size); }


bool xy_cmp(LocationType* a, LocationType* b) {
    if (a->pixel["xi"] < b->pixel["xi"]) {
        return true;
    } else if (a->pixel["xi"] > b->pixel["xi"]) {
        return false;
    } else {
        if (a->pixel["yi"] < b->pixel["yi"]) {
            return true;
        } else if (a->pixel["yi"] > b->pixel["yi"]) {
            return false;
        } else {
            return true;
        }
    }
}


char lookup_location_code(int age, int code) {
    char loc = (char) 0;
    if (school_codes.count(code) > 0 or (code == child_code and age >= school_age)) {
        loc = 's';
    } else if (work_codes.count(code) > 0) {
        loc = 'w';
    } else if (home_codes.count(code) or (code == child_code and age < school_age)) {
        loc = 'h';
    } else {
        cerr << "Error: encountered bad employment (EMPSTATD) code: " << code << endl;
        exit(-1);
    }
    return loc;
}


vector<HouseType*> import_households(string filename) {
    vector<HouseType*> hh_loc;
    cerr << "reading household locations\n";
    istringstream line;
    ifstream fh_in(filename);
	string buffer;
    int hh_ctr = 0;
    if (fh_in) {
        while (getline(fh_in, buffer)) {
            line.clear();
            line.str(buffer);
          //At this point, only includes houses
          /*id type x y x_ctr y_ctr
            0 house -89.6910220003 20.7015302009 -89.69159442 20.69995656
            1 house -89.700145483 20.6641877526 -89.69992776 20.66245653
            2 house -89.758249035 20.6954360352 -89.75826114 20.69578989
            3 house -89.6974011142 20.6551249287 -89.69576109 20.65412319*/

            int hid;
            string type;
            double x;
            double y;

            if (line >> hid >> type >> x >> y) {
                if (type != "house") {
                    cerr << "Expecting only houses in the location file at this point.  Found: " << buffer << endl;
                    exit(-1);
                }
                assert(hid == hh_ctr); // make sure hids are sequential integers starting with 0, because anything else would be stupid
                HouseType* h = new HouseType();
                h->hid = hid;
                h->x   = x;
                h->y   = y;

                hh_loc.push_back(h);
                hh_ctr++;
            }
        }
    }
    fh_in.close();
    return hh_loc;
}

/*struct LocationType {
    int workid;
    char type;
    double x; // longitude
    double y; // latitude
    map<string,int> pixel = {{"xi",0}, {"yi",0}};
    //int xi;
    //int yi;
    int workers;
    int students;
};*/


vector<LocationType*> import_workplaces_and_schools(const string filename, const int current_max_loc_id, map<int, int>& workplace_lookup, int& total_raw_workers) {
    vector<LocationType*> workplaces_and_schools;
    cerr << "reading workplaces & schools\n";
    int loc_id = current_max_loc_id + 1;
    total_raw_workers = 0;

    istringstream line;
    ifstream fh_in(filename);
	string buffer;
	ofstream fh_out("w_and_s_locs_w_id.txt");
    if (fh_in) {
        while (getline(fh_in, buffer)) {
            line.clear();
            line.str(buffer);
          /*W 1 -89.6264173747 20.9599660422
            W 2 -89.6116996054 20.964832419
            W 1 -89.6428676041 20.9771890368
            W 2 -89.6405575542 20.9678584161
            W 1 -89.6255278043 20.9746128005*/

            char type;
            int size;
            double x;
            double y;

            if (line >> type >> size >> x >> y) {
                LocationType* w = new LocationType();
                w->workid   = loc_id;
                w->type     = tolower(type);
                w->x        = x;
                w->y        = y;
                w->workers  = 0;
                w->students = 0;

                fh_out << w->workid << " " << w->type << " " << w->x << " " << w->y << endl;
                if (w->type == 'w') {
                    w->raw_workers = size;
                    total_raw_workers += w->raw_workers;
                } else {
                    w->raw_workers = 0;
                }

                w->pixel["xi"] = x_to_col_num(w->x);
                w->pixel["yi"] = y_to_row_num(w->y);
                workplaces_and_schools.push_back(w);
                loc_id++;
            }
        }
    }

    fh_in.close();
    fh_out.close();

    sort(workplaces_and_schools.begin(), workplaces_and_schools.end(), xy_cmp);

    for (unsigned int i = 0; i < workplaces_and_schools.size(); ++i) {
        workplace_lookup[workplaces_and_schools[i]->workid] = i;
    }

    return workplaces_and_schools;
}

vector<PersonType*> import_population(string filename, vector<HouseType*>& hh_loc, map<char, int>& day_loc_ctr) {
    vector<PersonType*> pop;
    /* This function imports a preliminary version of the population file
    that doesn't have people assigned to workplaces yet (perhaps among
    other things.)*/

    istringstream line;
    ifstream fh_in(filename);
    string buffer;
	//stringstream ss;
	//ofstream fh_out("w_and_s_locs_w_id.txt");
    int per_ctr = 0;
    if (fh_in) {
        cerr << "reading population\n";
        while (getline(fh_in, buffer)) {
            line.clear();
            line.str(buffer);
            int pid;
            int hid;
            int age;
            int sex;

            int _dummy_x;
            int _dummy_y;

            int workid;
            long int hh_serial;
            int pernum;
            int empstat;

       /* 
        pid hid age sex gridx gridy workid hh_serial pernum empstatd
        1 1 31 1 0 0 -1 2748179000 1 110
        2 1 29 2 0 0 -1 2748179000 2 110
        3 1 10 2 0 0 -1 2748179000 3 0
        4 2 32 1 0 0 -1 2748114000 1 110
       */ 

            if (line >> pid >> hid >> age >> sex >> _dummy_x >> _dummy_y >> workid >> hh_serial >> pernum >> empstat) {
                assert(pid == per_ctr); // make sure pids are sequential integers starting with 0, because anything else would be stupid
                PersonType* p = new PersonType();
                p->pid       = pid;
                p->hid       = hid;
                p->age       = age;
                p->sex       = sex;
                p->x         = hh_loc[hid]->x;
                p->y         = hh_loc[hid]->y;
                p->workid    = workid;
                p->hh_serial = hh_serial;
                p->pernum    = pernum;
                p->day_loc   = lookup_location_code(age, empstat);
                p->empstat   = empstat;

                day_loc_ctr[p->day_loc] += 1;
                pop.push_back(p);
                ++per_ctr;
            }
        }
    }
    cerr << "Population size: " << pop.size() << endl;
    cerr << "Total number of workers (IPUMS): " << day_loc_ctr['w'] << endl;
    cerr << "Total number of students (IPUMS): " << day_loc_ctr['s'] << endl;
    cerr << "Total number of homebodies (IPUMS): " << day_loc_ctr['h'] << endl;
    return pop;
}


unsigned int binary_search(vector<LocationType*> loc_vec, int search_coord, string label, unsigned int start, unsigned int end_plus_one) {
    if (start >= end_plus_one) return end_plus_one;
    if (loc_vec[start]->pixel[label] > search_coord)  return start;
    if (loc_vec[end_plus_one-1]->pixel[label] < search_coord) return end_plus_one;

    unsigned int imin = start; 
    unsigned int imax = end_plus_one;

    while (imin < imax) {
        unsigned int imid = imin + (imax-imin)/2; // TODO - verify that integer math is safe here
        if (loc_vec[imid]->pixel[label] < search_coord) {
            imin = imid+1;
        } else {
            imax = imid;
        }
    }

    return imin;
}

vector<int> get_nearby_places(int pxi, int pyi, char loc_type, vector<LocationType*> workplaces_and_schools, int num_loc_needed, int& positions_found) {
    int commute_range = -1;
    positions_found = 0;
    vector<int> nearby_places; // by index in workplaces_and_schools
    string pos_type = loc_type == 'w' ? "workers" : "students";

    while ((nearby_places.size() < num_loc_needed and nearby_places.size() < workplaces_and_schools.size()) or (positions_found <= 0)) {
        nearby_places.clear();
        commute_range += 1;
        //#print "\n\nEnvelope size:", 2*commute_range + 1, 'x', 2*commute_range + 1
        for (int x_val = pxi-commute_range; x_val <= pxi+commute_range; ++x_val) {
            unsigned int start_pos = 0;
            unsigned int end_pos_plus_one = workplaces_and_schools.size();
            start_pos = binary_search(workplaces_and_schools, x_val, "xi", start_pos, end_pos_plus_one);
            end_pos_plus_one = binary_search(workplaces_and_schools, x_val+1, "xi", start_pos, end_pos_plus_one);

            if (start_pos == end_pos_plus_one) { continue; }

            start_pos = binary_search(workplaces_and_schools, pyi-commute_range, "yi", start_pos, end_pos_plus_one);
            end_pos_plus_one = binary_search(workplaces_and_schools, pyi+commute_range+1, "yi", start_pos, end_pos_plus_one);

            for (unsigned int i = start_pos; i < end_pos_plus_one; ++i) {
                LocationType* w = workplaces_and_schools[i]; 
                if (pos_type == "students") {
                    positions_found++;
                    nearby_places.push_back(w->workid);
                } else if (w->workers > 0) {
                    positions_found += w->workers;
                    nearby_places.push_back(w->workid);
                }
            }
        }
    }
    return nearby_places;
}

/*
def select_nearest_school(px, py, nearby_places, workplaces_and_schools, workplace_lookup):
    //# nearby_places is a list of indeces for workplaces_and_schools
    closest_school, min_dist = workplaces_and_schools[workplace_lookup[nearby_places[0]]], () //# () evaluates to greater than any number
    for workid in nearby_places:
        s = workplaces_and_schools[workplace_lookup[workid]]
        //#print 'candidate:', workid, s
        d = haversine(px, py, s['x'], s['y'])
        //#print "dist:", d
        if d < min_dist:
            min_dist = d
            closest_school = s
    //#print "shortest distance:", min_dist
    return closest_school, min_dist


def send_kids_to_school(pop, pop_ids, workplaces_and_schools, total_raw_workers, workplace_lookup):
    schools = [location for location in workplaces_and_schools if location['type'] == 's']

    fo = open('student_placement.log', 'w')
    students_allocated = 0
    for pid in pop_ids:
        loc_type = pop[pid][field_idx['day_loc']]
        if loc_type != 's':
            continue
        //#print "looking at pid: ", pid
        num_loc_needed = 1

        px, py = pop[pid][field_idx['x']], pop[pid][field_idx['y']]
        pxi, pyi = x_to_col_num(px), y_to_row_num(py)

        nearby_places, positions_found = get_nearby_places(pxi, pyi, loc_type, schools, num_loc_needed)
        //#for s in nearby_places:
        //#    print s

        school, distance = select_nearest_school(px, py, nearby_places, workplaces_and_schools, workplace_lookup)
        //# 'school' is an element in the workplaces_and_schools list
        pop[pid][field_idx['workid']] = school['workid']
        school['students'] += 1
        school['raw_workers'] += 1.0/student_teacher_ratio
        total_raw_workers     += 1.0/student_teacher_ratio
        students_allocated += 1
        if students_allocated % 10000 == 0:
            print "students sent to school:", students_allocated
        //#print px, py, school, distance
        fo.write(' '.join(map(str,[pid, px, py, school['workid'], school['x'], school['y'], distance])) + '\n')
    fo.close()
    return total_raw_workers

def choose_workplace(px, py, nearby_places, workplaces_and_schools, workplace_lookup):
    raw_weights = [0.0 for i in range(len(nearby_places))]
    for i, workid in enumerate(nearby_places):
        w = workplaces_and_schools[workplace_lookup[workid]]
        dist = haversine(px, py, w['x'], w['y'])
        size = w['workers']
        raw_weights[i] = size / dist**2

    //# normalize weights
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
*/


double haversine(double lon1, double lat1, double lon2, double lat2) {
    // Calculate the great circle distance between two points 
    // on the earth (specified in decimal degrees)

    // convert decimal degrees to radians 
    lon1 = deg_to_rad(lon1);
    lon2 = deg_to_rad(lon2);
    lat1 = deg_to_rad(lat1);
    lat2 = deg_to_rad(lat2);
    // haversine formula 
    double dlon = lon2 - lon1;
    double dlat = lat2 - lat1;
    double a = pow(sin(dlat/2),2) + cos(lat1) * cos(lat2) * pow(sin(dlon/2),2);
    double c = 2 * asin(sqrt(a));
    double km = 6371 * c;
    return km;
}


int main() {
    // Import household location data
    vector<HouseType*> hh_loc = import_households("locations-yucatan_prelim.txt");
cerr << hh_loc.size() << endl;
    // Import workplace and school location data
    //
    // We are using workplace size (# employees) as a weight, but ignoring
    // school size data currently, as we feel it is more realistic to send
    // students to the nearest school.
    /*max_loc_id = hh_loc.size() - 1;
    map<int, int> workplace_lookup;
    int total_raw_workers = 0;
    vector<LocationType*> workplaces_and_schools = import_workplaces_and_schools("schools_and_workplaces.out", max_loc_id, workplace_lookup, total_raw_workers);

    // Import population data
    //pop = OrderedDict()
    vector<int> pop_ids;
    map<char, int> day_loc_ctr = {{'h', 0}, {'w', 0}, {'s', 0}}; // home / work / school
    map<> pop = import_population("population-yucatan_prelim.txt", pop_ids, hh_loc, day_loc_ctr);

    print "Total number of non-teacher jobs (DENUE):", total_raw_workers
    print "Student:Teacher ratio (World Bank):", student_teacher_ratio
    print "Total number of needed teachers:", day_loc_ctr['s']/student_teacher_ratio
    print "Total raw number of jobs (DENUE + needed teachers):", total_raw_workers + day_loc_ctr['s']/student_teacher_ratio

    shuffle(pop_ids)

    // Send kids to nearest school
    // Needs to happen first so we know how many teachers are needed
    total_raw_workers = send_kids_to_school(pop, pop_ids, workplaces_and_schools, total_raw_workers, workplace_lookup)

    // normalize workplace sizes
    employment_rescaling_factor = day_loc_ctr['w'] / float(total_raw_workers)
    for place in workplaces_and_schools:
        place['workers'] = place['raw_workers'] * employment_rescaling_factor

    // Filehandle for file we're going to write
    //fo = file('population-yucatan_no_copy.txt','w')
    fo = file('population-yucatan-silvio.txt','w')
    fo.write('pid hid age sex hh_serial pernum workid\n')

    // Make a copy so we can delete places from the original data structure as they fill up
    W_AND_S_COPY = deepcopy(workplaces_and_schools)
    ctr = 0

    // For each person
    for(const auto& [pid, person]: pop) {
    //for pid in pop.keys():
        //person = pop[pid]
        loc_type = person[field_idx['day_loc']]
        // Have them stay at home if they don't work or go to school
        if loc_type == 'h':
            person[field_idx['workid']] = person[field_idx['hid']] # person stays home
        // If they work, probabilistically choose a workplace based on where they live
        // and how many positions are available at each workplace
        elif loc_type == 'w':
            px, py = person[field_idx['x']], person[field_idx['y']]
            pxi, pyi = x_to_col_num(px), y_to_row_num(py)

            nearby_places, positions_found = get_nearby_places(pxi, pyi, loc_type, workplaces_and_schools, workplace_neighborhood)
            //print len(nearby_places), positions_found
            workplace = choose_workplace(px, py, nearby_places, W_AND_S_COPY, workplace_lookup)
            workplace['workers'] -= 1 # remove one available job
            // If the selected workplace no longer has openings,
            // remove it from the list, so we don't have to consider it again
            if workplace['workers'] <= 0:
                for i,v in enumerate(workplaces_and_schools):
                    if v['workid'] == workplace['workid']:
                        del workplaces_and_schools[i]
                        break
            person[field_idx['workid']] = workplace['workid'] # assign worker

        // Students already have the "workid" (prob should be called day_loc_id to avoid
        // confusion), so we don't have to do much for them, just output their info
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
    }
    fo.close()*/
}
