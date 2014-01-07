#!/usr/bin/python

locfile  = 'locations-yucatan_final.txt'
edgefile = 'tmp_yucatan.1.edge'
distfile = 'neighbor_distances.out'
new_edgefile = 'yucatan.edgelist'

locs = dict()
distances = []

header = True
for line in file(locfile):
    if header:
        header = False
        continue
    p = line.split()
    locs[p[0]] = [float(p[2]), float(p[3])] # map loc ID to long, lat

fo = file(distfile, 'w')
fo2 = file(new_edgefile, 'w')
fo2.write('x1 y1 x2 y2\n')

header = True
for line in file(edgefile):
    if line[0] == '#':
        continue
    if header:
        header = False
        continue
 
    p = line.split()
    L1, L2 = p[1], p[2]
    d = ((locs[L1][0] - locs[L2][0])**2 + (locs[L1][1] - locs[L2][1])**2)**0.5
    fo.write(str(d) + '\n')
    fo2.write(' '.join(map(str, [locs[L1][0], locs[L1][1], locs[L2][0], locs[L2][1], '\n'])))

    distances.append(d)

print "Average distance between neighbors in Delaunay triangulation:", float(sum(distances))/len(distances)
fo.close()
fo2.close()
