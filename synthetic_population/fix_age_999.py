#!/usr/bin/python
from sys import argv, exit
import random

if len(argv) < 4:
    print "\n\tUsage: ./fix_age_999.py age_column input_population_filename output_population_filename\n"
    print "\tIf input and output filenames are the same, file will be overwritten."
    print "\tage_column should be the index (starting with 0) of the column where age data is located in the input file\n"
    exit()

peeps = []
header = ''

for line in file(argv[2]):
    if not header:
        header = line
        continue
        
    peeps.append(line.split())

fo = file(argv[3], 'w')
fo.write(header)
age = int(argv[1])

for p in peeps:
    while p[age] == '999':
        p[age] = random.choice(peeps)[age]
    fo.write(' '.join(p) + '\n')
