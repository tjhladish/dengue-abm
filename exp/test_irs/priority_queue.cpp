#include <iostream>
#include <queue>
#include <vector>

using namespace std;

struct ITE {
    //ITE() {};
    ITE(int d, int s, float e) : duration(d), start(s), efficacy(e) {};
    int duration;
    int start;
    float efficacy;
};

inline bool operator<(const ITE lhs, const ITE rhs) {
    return lhs.start > rhs.start;
}

int main() {
    priority_queue<ITE> q;

    q.emplace( 10,3,0.1 );
    q.emplace( 7,4,0.2  );
    q.emplace( 4,2,0.3  );
    q.emplace( 5,7,0.4  );
    q.emplace( 14,1,0.5 );

    while (!q.empty()) {
        cout << q.top().duration << ' ' << q.top().start << ' ' << q.top().efficacy << endl;
        q.pop();
    }
}
