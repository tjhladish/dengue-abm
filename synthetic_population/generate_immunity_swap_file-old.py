#!/usr/bin/python
from sys import exit
 
'''
This python code has been rewritten as C++.  That version should be used,
although this may provide useful snippets for new scripts.

Functionality:
Determine, for everyone in the synthetic population over 0 years old, which 
of the individuals who are one year younger are closest.

Also, based on simple euclidean distance, calculate the relative probability
of inheriting the immune status from each of those people based on a simple
gravity model.
'''

'''
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
'''

locations = dict()
people = dict()
people_by_age = [[] for i in range(100)]

def euclidean(p1, p2):
    return ((people[p1]['x'] - people[p2]['x'])**2 + (people[p1]['y'] - people[p2]['y'])**2)**.5

def convert_to_probabilities(dists):
    probs = [(1.0/(pair[0]**2), pair[1]) for pair in dists]
    return probs

def normalize_probabilities(probs):
    cumsum = sum([pair[0] for pair in probs])
    normed = [(pair[0]/cumsum, pair[1]) for pair in probs]
    return normed

# Read in house locations to get coordinates
for line in file("locations-bangphae.txt"):
    if line[:2] == 'id': continue
    parts = line.split()
    locations[parts[0]] = {'x':float(parts[2]), 'y':float(parts[3])}

# Read in population, figure out where each person is and note the age
for line in file("population-bangphae.txt"):
    if line[:3] == 'pid': continue
    parts = line.split()
    pid = parts[0]
    people[pid] = {'x':locations[parts[1]]['x'], 'y':locations[parts[1]]['y'], 'age':int(parts[2]), 'sex':parts[3]}
    people_by_age[people[pid]['age']].append(pid) 

# For each person, get those that are 1 year younger
for pid1 in people.keys():
    dists = list()
    age = people[pid1]['age']
    if age >= 1:
        for pid2 in people_by_age[age - 1]:
            # And calculate the distance to that person
            pair = (euclidean(pid1, pid2), pid2)
            # For people in same location, make distance == 5 meters
            if pair[0] == 0:
                pair = (0.005, pair[1]) # units == km
            dists.append(pair)
        # Now sort them based on proximity (closest first)
        dists = sorted(dists, key=lambda pair: pair[0])

        # Forget about all but the closest 100 people
        dists = dists[:100]

        # Convert distances to gravity-like probabilities
        probs = convert_to_probabilities(dists)
        probs = normalize_probabilities(probs)
        
        for prob, pid2 in probs[:100]:
            print pid1, pid2, prob 
