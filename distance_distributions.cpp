#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <math.h>
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
        int id, age, hid, sex, pernum;
        long int workid;
        string hh_serial;
        if (line >> id >> hid >> age >> sex >> hh_serial >> pernum >> workid) {
            Person* peep = new Person();
            peep->id = id;
            peep->loc = locations[hid-1];
            peep->age = age;
            peep->sex = sex;
            //if (people.size() != id - 1) cerr << "something is screwed up with person ID's and indexing\n";
            people.push_back(peep);
            people_by_age[age].push_back(people[id-1]);
        }
    }
    iss.close();
    cerr << "loadPopulation: " << people.size() << endl;
    return true;
}

/*
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

bool loadLocations(string locFilename) {
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

double euclidean(Person* p1, Person* p2) {
    return sqrt(pow(p1->loc->x - p2->loc->x, 2) + pow(p1->loc->y - p2->loc->y,2));
}


int main() { 
    loadLocations("locations-yucatan.txt");
    loadPopulation("population-yucatan.txt");

    const int NUM_PEOPLE = 1819497;
    const int NUM_NEIGHBORS = 10;

    int** pid_mat = new int*[NUM_PEOPLE];
    for(int i = 0; i < NUM_PEOPLE; ++i) pid_mat[i] = new int[NUM_NEIGHBORS];
 
    float** dist_mat = new float*[NUM_PEOPLE];
    for(int i = 0; i < NUM_PEOPLE; ++i)dist_mat[i] = new float[NUM_NEIGHBORS];
    
    #pragma omp parallel for num_threads(10)
    for(int i=0; i < people.size();i++ ) {
        Person* p1 = people[i];
        if(!p1 || p1->age<1) continue;

        priority_queue<Distance, vector<Distance>, comp > distances;
        for(int j=0; j <people_by_age[p1->age-1].size();j++ ) {
            Person* p2 = people_by_age[p1->age-1][j];
            
            Distance dist;
            dist.dist = euclidean(p1, p2 );
            dist.id   = p2->id;
            distances.push(dist);
        }
        int j = 0;
        for (j = 0; j<NUM_NEIGHBORS && !distances.empty(); j++) {
            Distance dist = distances.top(); distances.pop();
            dist_mat[p1->id - 1][j] = dist.dist;
            pid_mat[p1->id - 1][j] = dist.id;
        //    printf("%d %d %f\n", p1->id, dist.id, dist.dist);
        }
        while (j++<NUM_NEIGHBORS) {
            dist_mat[p1->id - 1][j] = -1;
            pid_mat[p1->id - 1][j] = -1;
        }
    }

    for(int i=0; i < NUM_PEOPLE;i++ ) {
        for (int j = 0; j<NUM_NEIGHBORS; j++) {
            printf("%d %d %f\n", i+1, pid_mat[i][j], dist_mat[i][j]);
        }
    }
    return 0;
}

        
/*        
for pid1 in people.keys():
    dists = list()
    age = people[pid1]['age']
    if age >= 1:
        for pid2 in people_by_age[age - 1]: # people 1 year younger
            pair = (euclidean(pid1, pid2), pid2) 
            dists.append(pair)
        dists = sorted(dists, key=lambda pair: pair[0])
        for close_person in dists[:100]:
            print pid1, close_person[1], close_person[0]
*/            
