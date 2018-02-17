#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>


using namespace std;


vector<int> read_pop_ids (string filename) {
    ifstream is(filename);
    istream_iterator<double> start(is), end;
    vector<int> ids(start, end);
    cout << "Read " << ids.size() << " numbers" << std::endl;

    return ids;
}

int main() {
    vector<int> ids = read_pop_ids("8-14_merida_ids.txt");

    for (unsigned int i = 0; i<5; i++) cout << ids[i] << endl;
    cout << "--------------\n";
    for (unsigned int i = ids.size() - 5; i<ids.size(); i++) cout << ids[i] << endl;

}
