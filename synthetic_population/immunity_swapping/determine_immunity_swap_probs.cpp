#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <cmath>
#include <assert.h>

/******************************************************************************

    This program reads in a population file and a locations file, and calc-
    ulates the distances to NUM_NEIGHBORS nearest people who are 1 year
    younger in the synthetic population.  PIDs for those neighbors are also
    reported.
    
    It is intended to be used to calculate immunity swap probabilities, so
    that spatial structure can be maintained between years in the simulator.

******************************************************************************/


using namespace std;

struct Location { double x; double y; };
struct Person { int id; Location* loc; int age; string sex; };
struct Distance { int id; double dist; };

class comp {
    public:
        bool operator() (const Distance& lhs, const Distance& rhs) const {
            return (lhs.dist>rhs.dist);
        }
};

const int MAX_AGE = 100;

vector<Location*> locations;
vector<Person*> people;
vector< vector<Person*> > people_by_age(MAX_AGE + 1, vector<Person*>(0));

/*'''
tjhladish@capybara:~/work/dengue$ head population-bangphae.txt 
pid hid age sex gridx gridy workid
1 1 35 F 1 1 55560
2 1 10 M 1 1 55060
3 2 42 M 1 1 55559
4 2 37 F 1 1 55625
5 2 15 F 1 1 55793
6 2 12 F 1 1 54505
7 2 12 M 1 1 54505
8 3 39 M 1 1 55557
9 3 35 F 1 1 56035

tjhladish@roc:~/work/dengue$ head population-yucatan.txt 
pid hid age sex hh_serial pernum workid
1 1 31 1 2748179000 1 442670
2 1 29 2 2748179000 2 395324
3 1 10 2 2748179000 3 468423
4 2 32 1 2748114000 1 397104
5 2 30 2 2748114000 2 396166

'''*/

double radians(double deg) {
    // M_PIl is the long double PI from cmath
    return M_PIl * deg / 180.0;
}


double haversine(double lon1, double lat1, double lon2, double lat2) {
    /* 
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    */
    // convert decimal degrees to radians 
    lon1 = radians(lon1);
    lat1 = radians(lat1);
    lon2 = radians(lon2);
    lat2 = radians(lat2);
    // haversine formula 
    double dlon = lon2 - lon1;
    double dlat = lat2 - lat1;
    double a = pow(sin(dlat/2),2) + cos(lat1) * cos(lat2) * pow(sin(dlon/2),2);
    double c = 2 * asin(sqrt(a));
    double km = 6371 * c;
    return km;
}


double haversine(Person* p1, Person* p2) {
    return haversine(p1->loc->x, p1->loc->y, p2->loc->x, p2->loc->y);
}


bool loadPopulation(string popFilename) {
    ifstream iss(popFilename.c_str());
    if (!iss) {
        cerr << "ERROR: " << popFilename << " not found." << endl;
        return false;
    }
    while (iss) {
        char buffer[500];
        iss.getline(buffer,500);
        istringstream line(buffer);
        int id, age, hid, pernum_ignored;
        long int workid;
        string sex, hh_serial_ignored;
        if (line >> id >> hid >> age >> sex >> hh_serial_ignored >> pernum_ignored >> workid) {
            Person* peep = new Person();
            peep->id = id;
            peep->loc = locations[hid-1];
            peep->age = age;
            peep->sex = sex;
            people.push_back(peep);
            people_by_age[age].push_back(people[id-1]);
        }
    }
    iss.close();
    cerr << "loadPopulation: " << people.size() << endl;
    return true;
}


bool loadLocations(string locFilename) {
    /*
    Expected structure of locations file:
    tjhladish@capybara:~/work/dengue$ head locations-bangphae.txt 
    id type x y
    1 house 0.849513080669567 0.304251517867669
    2 house 0.00842146878130734 0.657222841167822
    3 house 0.369299583602697 0.349579052301124
    4 house 0.460470566991717 0.277130988193676
    5 house 0.0256796989124268 0.712766533251852
    6 house 0.431369476253167 0.20579392160289
    7 house 0.355606281897053 0.540817143395543
    8 house 0.996641094330698 0.151982805691659
    9 house 0.327530618291348 0.324621923966333
    */
    ifstream iss(locFilename.c_str());
    if (!iss) {
        cerr << "ERROR: " << locFilename << " not found." << endl;
        return false;
    }
    while (iss) {
        char buffer[500];
        iss.getline(buffer,500);
        istringstream line(buffer);
        int locID;
        string type;
        double xcoord, ycoord;
        if (line >> locID >> type >> xcoord >> ycoord) {
            Location* loc = new Location();
            loc->x = xcoord;
            loc->y = ycoord;
            if (locations.size() != locID - 1) cerr << "something is screwed up with location ID's and indexing\n";
            locations.push_back(loc);
        }
    }
    iss.close();
    cerr << "loadLocations: " << locations.size() << endl;
    return true;
}


int main(int argc, char** argv) { 

    if (argc != 3) {
        cerr << "\n\tUsage: ./generate_swap_file <locations_file> <population_file> > swap_filename.txt\n\n";
        return -5;
    }

    string locations_filename  = argv[1];
    string population_filename = argv[2];

  //const int NUM_PEOPLE = 1819497; // yucatan
  //const int NUM_PEOPLE =  207591; // bangphae 
    const int NUM_PEOPLE =  839659; // merida
    cerr << "Assuming a population size of " << NUM_PEOPLE << endl;
    const int NUM_NEIGHBORS = 1;           // Number of nearest neighbors to report for swap file
    const double MIN_DISPLACEMENT = 0.005;  // Distance (in km) to use between people who live in the same house

    loadLocations(locations_filename);
    loadPopulation(population_filename);

    assert (people.size() == NUM_PEOPLE);

    // Initialize matrix to hold pids of nearest neighbors
    int** pid_mat = new int*[NUM_PEOPLE];
    for(int i = 0; i < NUM_PEOPLE; ++i) pid_mat[i] = new int[NUM_NEIGHBORS];
 
    // Initialize matrix to hold distances of nearest neighbors
    float** dist_mat = new float*[NUM_PEOPLE];
    for(int i = 0; i < NUM_PEOPLE; ++i) dist_mat[i] = new float[NUM_NEIGHBORS];
    
    // Initialize matrix to hold probabilities of selecting each nearest neighbor
    float** prob_mat = new float*[NUM_PEOPLE];
    for(int i = 0; i < NUM_PEOPLE; ++i) prob_mat[i] = new float[NUM_NEIGHBORS];

    #pragma omp parallel for num_threads(10) // change num_threads based on available cores
    // For each person in the population
    for(int i = 0; i < people.size(); i++) {
        if (i % 1000 == 0) cerr << i << endl;
        Person* p1 = people[i];
        // If they're at least 1 year old
        if(!p1 || p1->age<1) continue;

        // Calculate the euclidean distances from this person to all individuals
        // who are one year younger.  Store than in a priority queue sorted
        // by distance, closest first.
        // TODO: Should use great circle distance, not simple euclidean
        priority_queue<Distance, vector<Distance>, comp > distances;
        for(int j = 0; j < people_by_age[p1->age-1].size(); j++) {
            Person* p2 = people_by_age[p1->age-1][j];
            
            Distance dist;
            dist.dist = haversine(p1, p2 );
            dist.id   = p2->id;
            distances.push(dist);
        }

        // Take the first NUM_NEIGHBORS people in the distances queue, and store
        // their information in the pid and distance matricies
        int j = 0;
        double weight_total = 0.0;
        int p1_idx = p1->id - 1;

        for (j = 0; j<NUM_NEIGHBORS && !distances.empty(); ++j) {
            Distance dist = distances.top(); distances.pop();
            if (dist.dist <= MIN_DISPLACEMENT) {
                // fprintf(stderr, "Setting displacement of %f to %f km\n", dist.dist, MIN_DISPLACEMENT);
                dist.dist = 0.005;
            }
            double weight = 1.0 / pow(dist.dist, 2);
            weight_total += weight;
            prob_mat[p1_idx][j] = weight;
            dist_mat[p1_idx][j] = dist.dist;
            pid_mat[p1_idx][j] = dist.id;
        }
    
        // normalize probabilities
        for (int k = 0; k<j; ++k) {
            prob_mat[p1_idx][k] /= weight_total; 
        }
        
        // In case fewer than NUM_NEIGHBORS people exist in the population who are
        // one year younger, fill with dummy values.
        while (j++<NUM_NEIGHBORS) {
            prob_mat[p1_idx][j] = -1;
            dist_mat[p1_idx][j] = -1;
            pid_mat[p1_idx][j] = -1;
        }
    }

    // Output each person's pid, the neighbor's pid, and the distance to them
    for(int i=0; i < NUM_PEOPLE;i++ ) {
        for (int j = 0; j<NUM_NEIGHBORS; j++) {
            printf("%d %d %f %f\n", i+1, pid_mat[i][j], prob_mat[i][j], dist_mat[i][j]);
        }
    }
    return 0;
}
