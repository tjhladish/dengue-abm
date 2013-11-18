#!/usr/bin/python
from sys import argv

'''
tjhladish@roc:~/work/dengue$ head swap_probs-yucatan.txt 
1 3348 0.002366
1 2319 0.002478
1 4131 0.003218
1 3289 0.003819
1 5858 0.004130
1 2883 0.004452
1 1945 0.004505
1 4235 0.004649
1 3007 0.005890
1 2982 0.005942
'''
def output_probs(id1, weights, neighbors):
    tot = sum([w for w in weights if w > 0.0])
    probs = []
    if tot > 0:
        probs = [w/tot for w in weights]
    else:
        probs = [0.0] * len(weights)
    
    for i, p in enumerate(probs):
        if probs[i] > 0: 
            fo.write(id1 + ' ' + neighbors[i] + ' ' + str(p) + '\n')
        else:
            fo.write(id1 + ' -1 0\n')
    weights = []
    neighbors = []
    return weights, neighbors

if len(argv) != 3:
    print "\n\tUsage:\t./normalize_swap_probs.py <unnormalized_swap_probs> <normalized_swap_probs>\n\n"

fo = open(argv[2], 'w')

weights = []
neighbors = []
last_id = -1
for line in file(argv[1]):
    id1, id2, wt = line.split()
    wt = float(wt)
    if id1 != last_id and last_id != -1:
        weights, neighbors = output_probs(last_id, weights, neighbors)
    last_id = id1
    weights.append(wt)
    neighbors.append(id2)

output_probs(id1, weights, neighbors)

fo.close()
