#!/usr/bin/python

tag = ''
printing = False
tags = set()
for line in file('baseline_time_series.monthly'):
    p = line.strip().split()
    if len(p) == 1:
        tag = p[0]
        if tag in tags:
            printing = False
        else:
            printing = True
        tags.add(tag)
    elif printing:
        print tag, ' '.join(p[1:])
