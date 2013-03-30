#!/usr/bin/python
from glob import glob
from Levenshtein import distance, jaro
from string import lowercase

def match_locality(string, localities):
    ''' Try to figure out which locality 'string' is by
    finding the known localities that have the minimum 
    (levenshtein) edit distances '''
    #best = () # for comparisons, () is infinity
    best = 0
    matches = []
    for loc in localities:
        similarity = jaro(string.lower(), loc['loc_name'].lower())
        #d = distance(string.lower(), loc['loc_name'].lower())
        if similarity > best:
        #if d < best:
            matches = []
            best = similarity
            #best = d
        if similarity == best: 
        #if d == best: 
            hit = dict(loc)
            hit['score'] = similarity
            matches.append(hit)
    return matches


def remove_funny_chars(string):
    newstr = []
    for c in string.lower():
        if c in lowercase:
            newstr.append(c)
    return ''.join(newstr)


def find_idxs_of_min(L):
    idx = [0]
    val = L[0]
    for i in range(1,len(L)):
        if val < L[i]:
            idx = [i]
            val = L[i]
        elif val == L[i]:
            idx.append(i)
    return idx


def guess_again(orig_string, matches):
    ''' Sometimes the mismatch is due to inconsistent use of
    accent marks.  See if any of the matches are particularly
    close when accent marks are disregarded '''
    alt_scores = []
    string = remove_funny_chars(orig_string) 
    for m in matches:
        alt_scores.append(distance(string, remove_funny_chars(m['loc_name'])))
    return alt_scores
    #return find_idxs_of_min(alt_scores)


# import locality --> municipality map (many-to-one mapping)
localities = [] 
for line in file("/home/tjhladish/work/dengue/raw_data/yucatan_locality_names_and_numbers"):
    parts = line.strip().split(',')
    localities.append(dict(zip(['state_num', 'state_name', 'muni_num', 'muni_name', 'loc_num', 'loc_name'], parts)))

total_cases_filenames = sorted(glob('*t.csv'))
fo = open("total_yucatan_cases.csv", 'w')

for filename in total_cases_filenames:
    year = filename[12:16]
    #print year
    start = False
    for line in file(filename):
        parts = line.strip().split(';')
        if parts[0] == 'YUCATAN':
            #print "  Found Yucatan"
            start = True
        elif parts[0] == 'ZACATECAS' or parts[0] == 'Total general':
            #print "  Ending Yucatan"
            start = False
        elif start:
            matches = match_locality(parts[0], localities)
            print parts[0], year
            best_guesses = [0] 
            alt_scores = ['NA']
            if len(matches) > 1:
                alt_scores = guess_again(parts[0], matches)
                best_guesses = find_idxs_of_min(alt_scores)
            for i in range(len(matches)):
                match = matches[i]
                print '\t', match['score'], alt_scores[i], match['loc_num'], match['loc_name'], match['muni_num'], match['muni_name'],
                if len(matches) == 1:
                    print "***"
                elif len(best_guesses) == 1 and best_guesses[0] == i:
                    print "**"
                elif i in best_guesses:
                    print "*"
                else:
                    print
            #fo.write(year + ',' + ','.join(parts[:53]) + '\n')
            print
        else:
            continue
    print "#################################################"
fo.close()
