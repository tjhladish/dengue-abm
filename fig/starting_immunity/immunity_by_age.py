#!/usr/bin/python

#expecting input files names immunity.1, immunity.2 ...
for k in range(0,51):
    print "running", k
    fo = open('imm_vs_age.year'+str(k),'w')
    pop_size = 0
    imm = [[0,0,0,0,0] for i in range(101)]

    header = True
    for line in file('immunity.year'+str(k)):
        '''
        pid age imm1 imm2 imm3 imm4
        1 31 1 0 1 0
        2 29 1 1 0 1
        3 10 1 0 0 0
        4 32 1 1 0 0
        5 30 1 0 0 1
        '''
        if header:
            header = False
            continue

        p = map(int,line.split())

        # haven't resolved what to do with people aged '999'
        if p[1] > 100:
            continue
        imm[p[1]][0] += 1
        for i, state in enumerate(p[2:6]):
            if state == 0:
                imm[p[1]][i+1] += 0
            else:
                imm[p[1]][i+1] += 1

    for age, vals in enumerate(imm):
        fo.write(' '.join([str(age)] + map(str,vals)) + '\n')

    fo.close()
