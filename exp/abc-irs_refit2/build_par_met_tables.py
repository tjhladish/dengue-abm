#!/usr/bin/python
from sys import argv
import numpy as np
import json, sqlite3
import to_precision as prec

sig_figs = 3
iq_range = 95.0
q_min    = (100.0-iq_range)/2
q_max    = 100.0 - (100.0-iq_range)/2
pad      = 8
sep      = '&'
eol      = '\\\\'

if (len(argv) != 3):
    print "\n\tUsage:\tbuild_par_met_tables.py <json_filename> <sqlite_db_filename>\n"
    exit(-1)

def fmt(val):
    return prec.std_notation(val, sig_figs).ljust(pad)

# parse json file
fi = file('abc-irs_refit2.json', 'r')
jdat = json.load(fi)

# slurp parameters from db
con = sqlite3.connect(argv[2])
con.row_factory = sqlite3.Row
c = con.cursor()
c.execute('select * from upar;')
sql_pars = c.fetchall()
pars = np.array(sql_pars)
pars = np.delete(pars, [0,1], 1).astype(float) # delete serials from data (0 specifies which row/col, 1 specifies col and not row)

# build stats matrix
stats = {} 
stats['median'] = np.percentile(pars, 50, axis=0)
stats['q_min'] = np.percentile(pars, q_min, axis=0)
stats['q_max'] = np.percentile(pars, q_max, axis=0)

for par_idx in range(len(jdat['parameters'])):
    j_par = jdat['parameters'][par_idx]
    print j_par['name'].ljust(35), sep, #fmt(j_par['value']), sep,
    print fmt(stats['median'][par_idx]), sep,
    print '(' + fmt(stats['q_min'][par_idx]) + ', ' + fmt(stats['q_max'][par_idx]) + ')' + eol





# slurp metrics from db
c.execute('select * from met;')
sql_mets = c.fetchall()
#labels = sql_mets[0].keys()
mets = np.array(sql_mets)
mets = np.delete(mets, 0, 1) # delete serials from data (0 specifies which row/col, 1 specifies col and not row)

# build stats matrix
stats = {} 
stats['mean'] = np.mean(mets, axis=0)
stats['median'] = np.percentile(mets, 50, axis=0)
stats['q_min'] = np.percentile(mets, q_min, axis=0)
stats['q_max'] = np.percentile(mets, q_max, axis=0)

for met_idx in range(len(jdat['metrics'])):
    j_met = jdat['metrics'][met_idx]
    print j_met['name'].ljust(20), sep, fmt(j_met['value']), sep,
    print fmt(stats['mean'][met_idx]), sep,
    print fmt(stats['median'][met_idx]), sep,
    print '(' + fmt(stats['q_min'][met_idx]) + ', ' + fmt(stats['q_max'][met_idx]) + ')' + eol
