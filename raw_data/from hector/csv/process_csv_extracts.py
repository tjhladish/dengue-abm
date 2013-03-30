#!/usr/bin/python
import unicodedata
from glob import glob
from Levenshtein import distance, jaro
from string import lowercase


def strip_accents(s):
   return ''.join(c for c in unicodedata.normalize('NFD', s)
                     if unicodedata.category(c) != 'Mn')


def match_locality(string, localities):
    ''' Try to figure out which locality 'string' is by
    finding the known localities that have the maximum 
    (jaro) similarity score '''
    if ',' in string:
        parts = string.split(',')
        string = parts[1].strip() + ' ' + parts[0]
    best = 0
    matches = []
    jaro_hits = []
    substring_hits = []
    for loc in localities:
        ulower_string = strip_accents(string.lower().decode("utf-8"))
        ulower_loc_name = strip_accents(loc['loc_name'].lower().decode("utf-8"))
        
        smaller, bigger = ulower_string, ulower_loc_name
        if len(bigger) < len(smaller):
            smaller, bigger = bigger, smaller
        if len(smaller) < len(bigger): # they might actually be the same size
            if smaller in bigger:
                hit = dict(loc)
                hit['score'] = 'sub' 
                substring_hits.append(hit)

        similarity = jaro(ulower_string, ulower_loc_name)
        if similarity > best:
            matches = []
            best = similarity
        if similarity == best: 
            hit = dict(loc)
            hit['score'] = similarity
            matches.append(hit)
    
    jaro_hits = [(m['loc_name'],m['muni_name']) for m in matches]
    for s in substring_hits:
        if (s['loc_name'],s['muni_name']) not in jaro_hits:
            matches.append(s)
    #print matches
    return matches


#def remove_funny_chars(string):
#    newstr = []
#    for c in string.lower():
#        if c in lowercase:
#            newstr.append(c)
#    return ''.join(newstr)


def find_idxs_of_max(scores):
    L = [s if s is not 'sub' else -1 for s in scores]
    idx = [0]
    val = L[0]
    for i in range(1,len(L)):
        if L[i] > val:
            idx = [i]
            val = L[i]
        elif L[i] == val:
            idx.append(i)
    return idx


#def guess_again(orig_string, matches):
#    ''' Sometimes the mismatch is due to inconsistent use of
#    accent marks.  See if any of the matches are particularly
#    close when accent marks are disregarded '''
#    alt_scores = []
#    string = remove_funny_chars(orig_string) 
#    for m in matches:
#        alt_scores.append(distance(string, remove_funny_chars(m['loc_name'])))
#    return alt_scores


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
#            print parts[0], year
            best_guesses = [0] 
            scores = [m['score'] for m in matches]
            alt_scores = ['NA']
            if len(matches) > 1:
                #alt_scores = guess_again(parts[0], matches)
                best_guesses = find_idxs_of_max(scores)
            for i in range(len(matches)):
                if i == 0:
                    print year,';', parts[0],';',
                else:
                    print ' ; ;',
                match = matches[i]
                print match['loc_name'],';', match['loc_num'],';', match['muni_name'],';', match['muni_num'],';', match['score'], 
                #print '\t', match['score'], alt_scores[i], match['loc_num'], match['loc_name'], match['muni_num'], match['muni_name'],
                if len(matches) == 1:
                    print "; ***"
                elif len(best_guesses) == 1 and best_guesses[0] == i:
                    print "; **"
                elif i in best_guesses:
                    print "; *"
                else:
                    print "; "
            #fo.write(year + ',' + ','.join(parts[:53]) + '\n')
            #print
        else:
            continue
fo.close()
