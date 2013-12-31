#!/usr/bin/python
import codecs
from ipums_xml_parser import *
from xml.dom import minidom
from sys import exit
from collections import defaultdict
from random import random
from sys import exit

class Census:
    def __init__(self):
        self.data = dict()
        self.muni_lookup = dict()

    def insert(self, name, num, pop):
        self.data[num] = pop
        self.muni_lookup[name] = num

    def get(self, key):
        pop = 0
        try:
            int(key)
            pop = self.data[key]
        except ValueError:
            pop = self.data[self.muni_lookup[key]]
        return pop

    def display(self):
        for k, v in sorted(self.muni_lookup.iteritems()):
            print k, v, self.data[v]

class Ipums_Person:
    fields = list()
    def __init__(self, vals):
        self.data = list(vals)

class Ipums_Household:
    fields = list()
    def __init__(self, hh_vals):
        self.data   = list(hh_vals)
        self.people = list()
        self.serial = hh_vals[Ipums_Household.fields.index('SERIAL')]
        self.normed_wt = 0

    def insert(self, person_vals):
        self.people.append(Ipums_Person(person_vals))

class Ipums:
    fields = list()
    hh_field_pos = list()
    person_field_pos = list()
    households_by_muni = defaultdict(list)

    def __init__(self):
        self.households = dict()
    
    def insert(self, row, pos):
        serial = row[pos['SERIAL']]
        muni = row[pos['MUNIMX']]
        if serial not in self.households:
            self.households[serial] = Ipums_Household([row[i] for i in Ipums.hh_field_pos])
        self.households[serial].insert([row[i] for i in Ipums.person_field_pos])
        Ipums.households_by_muni[muni].append(self.households[serial])

    def get_households(self, muni):
        return Ipums.households_by_muni[muni]

    def normalize_household_weights(self):
        # location of household weight in data structure
        wt_idx = Ipums_Household.fields.index('WTHH') 
        for m in Ipums.households_by_muni.keys():
            total_wt = 0
            for hh in Ipums.households_by_muni[m]:
                total_wt += float(hh.data[wt_idx])
            for hh in Ipums.households_by_muni[m]:
                hh.normed_wt = float(hh.data[wt_idx]) / total_wt

    def sample_households(self, m):
        rand = random()
        chosen_hh = Ipums.households_by_muni[m][0]
        for hh in Ipums.households_by_muni[m]:
            if rand < hh.normed_wt:
                chosen_hh = hh
                break
            else:
                rand -= hh.normed_wt
        return chosen_hh


def import_census(filename):
    census = Census()
    for line in codecs.open(filename, encoding='ISO-8859-1'):
        p = line.strip().split('\t')
        if p[5] == 'TOTAL MUNICIPAL':
            census.insert(p[3], p[2], int(p[9]))
    return census

def import_ipums_data(filename):
    known_household_fields = ['CNTRY','YEAR','SAMPLE','SERIAL','PERSONS','WTHH','STATEMX','MUNIMX','SIZEMX'] 
    known_person_fields    = ['PERNUM','WTPER','AGE','SEX','SCHOOL','EDATTAN','EDATTAND','YRSCHL','EDUCMX','EMPSTAT',
                              'EMPSTATD','OCCISCO','OCC','INDGEN','IND','CLASSWK','CLASSWKD','HRSWRK1','HLTHFAC'] 
    required_fields        = ['SERIAL','PERSONS','WTHH','MUNIMX','AGE'] 
    ipums = Ipums()
    header = list()
    pos    = dict() # holds position (index) of each field
    ctr = -1
    for line in file(filename):
        ctr += 1
        if ctr % 50000 == 0:
            print "parsed", ctr, "IPUMS records"
        if len(header) == 0: # we're looking at the first line/header
            header = [field.strip('"') for field in line.strip().split(',')]
            error = False
            for i,f in enumerate(header):
                if f not in known_household_fields and f not in known_person_fields:
                    print "Error: unknown field found in IPUMS csv file:", f
                    error = True
                else:
                    if f in known_household_fields:
                        Ipums_Household.fields.append(f)
                        Ipums.hh_field_pos.append(i)
                    if f in known_person_fields:
                        Ipums_Person.fields.append(f)
                        Ipums.person_field_pos.append(i)

            for r in required_fields:
                if r not in header:
                    print "Error: required field not found in IPUMS csv file:", r
                    error = True
            if error == True:
                print "Error: Failed to parse IPUMS csv file"
                sys.exit(-1)
            else:
                Ipums.fields = list(header)
                pos = dict(zip(header, range(len(header))))
        else: # we're looking at a data row
            ipums.insert(line.strip().split(','), pos)

    ipums.normalize_household_weights()
    return ipums


class Pixel:
    width = 0 # pixels are square, so height == width
    def __init__(self, _long, _lat, _lum, _muni_num, _muni_name):
        self.long      = _long
        self.lat       = _lat
        self.lum       = _lum
        self.muni_num  = _muni_num
        self.muni_name = _muni_name
        self.normed_wt = 0

    def locate_hh(self):
#        x = self.long
#        y = self.lat
        x = self.long - Pixel.width/2.0 + Pixel.width*random()
        y = self.lat  - Pixel.width/2.0 + Pixel.width*random()
        return x, y

class Municipality_Weights:
    # lookup municipality number used in gis data using other keys
    lookup_by_gis_name   = dict()
    lookup_by_census_num = dict()
    lookup_by_ipums_name = dict()
    lookup_by_ipums_num  = dict()

    def __init__(self):
        self.pixels_by_muni = defaultdict(list) # key is gis muni num, value is list of Pixels
        #self.pixels_by_muni = dict() # key is gis muni num, value is list of Pixels

    def insert_pixel(self, px):
       # if px.muni_num not in self.pixels_by_muni:
       #     self.pixels_by_muni[px.muni_num] = list()

        self.pixels_by_muni[px.muni_num].append(px)

    def define_muni_lookups(self, filename):
        header = True
        for line in file(filename):
            if header:
                header = False
                continue
            gis_num, gis_name, ipums_num, ipums_name = line.strip().split('\t')
            Municipality_Weights.lookup_by_gis_name[gis_name]        = gis_num
            # ipums form is state+municipality, e.g. 31001 (31 is yucatan).  Census is just municipality.
            Municipality_Weights.lookup_by_census_num[ipums_num[2:]] = gis_num
            Municipality_Weights.lookup_by_ipums_name[ipums_name]    = gis_num
            Municipality_Weights.lookup_by_ipums_num[ipums_num]      = gis_num


    def normalize_pixel_weights(self):
#        print "normalizing px weights"
        for gis_num in self.pixels_by_muni.keys():
            muni_total_luminosity = 0
            for px in self.pixels_by_muni[gis_num]:
                muni_total_luminosity += px.lum
            # normalize pixels
            for px in self.pixels_by_muni[gis_num]:
                px.normed_wt = px.lum / muni_total_luminosity


    def sample_pixels(self, m):
#        print "sampling"
        rand = random()
        chosen_px = self.pixels_by_muni[m][0]
        for px in self.pixels_by_muni[m]:
#            print rand, px.normed_wt
            if rand < px.normed_wt:
                chosen_px = px 
                break
            else:
                rand -= px.normed_wt
        return chosen_px



def import_nighttime_light_data(filename, luminosity_threshold):
    muni_wt = Municipality_Weights()
    for line in file(filename):
        _long, _lat, _lum, _muni_num, _muni_name = line.strip().split('\t')
        _long, _lat, _lum = map(float, [_long, _lat, _lum])  
        
        # Subtract off the natural albedo of the surface
        if _lum < luminosity_threshold:
            _lum = 0
        else:
            _lum -= luminosity_threshold

        px = Pixel(_long, _lat, _lum, _muni_num, _muni_name)
        muni_wt.insert_pixel(px)
    
    muni_wt.normalize_pixel_weights()
    return muni_wt


census = import_census("./raw_data/mexico_census/iter_31TXT.txt")
#census.display()
#print census.get(95)

# Import categorical information about the IPUMS data
IPUMS_XML_Parser("./raw_data/ipums/ipumsi_00005.xml")
ipums_lookup = IPUMS_XML_Parser.data
#print sorted(ipums_lookup['MUNIMX'])

# ipums_data is a 2D dictionary with the column name as the first key and the categorical values as the second key
# The returned value is what that maps to, e.g. ipums_data['SEX']['2'] maps to 'Female'
ipums_data = import_ipums_data("./raw_data/ipums/ipumsi_00005.csv")
age_field = Ipums_Person.fields.index('AGE')
sex_field = Ipums_Person.fields.index('SEX')
pernum_field = Ipums_Person.fields.index('PERNUM')
empstatd_field = Ipums_Person.fields.index('EMPSTATD')

# Import nighttime light data
# Luminosity values are given for pixel centered on the given coordinates
# Pixel height and width is expected to be 0.00416667 degrees
Pixel.width = 0.00416667
print "Warning: assuming square pixels with height and width ==", Pixel.width, "degrees"

luminosity_threshold = 0.16 # Anything dimmer than this is likely not a "light"

muni_weights = import_nighttime_light_data('./raw_data/all_yucatan_pixels.out', luminosity_threshold)
muni_weights.define_muni_lookups('./raw_data/muni_lookup')

size_idx = Ipums_Household.fields.index('PERSONS') # PERSONS is household size

hh_output_file = file('locations-yucatan.txt','w')
hh_output_file.write('id type x y x_ctr y_ctr\n')
pop_output_file = file('population-yucatan.txt','w')
pop_output_file.write('pid hid age sex gridx gridy workid hh_serial pernum empstatd\n')

# person ID counter.  Previous model uses IDs that start at 1; not sure whether that's necessary
pid = 1
hid = 1

# For each municipality
for m in sorted(Municipality_Weights.lookup_by_ipums_num.keys()):
    # Get expected number of households
    expected_hh_size = 0
    for hh in ipums_data.get_households(m):
        # Some households are more "common" according to IPUMS, and therefore
        # should be weighted more heavily.
        expected_hh_size += int(hh.data[size_idx]) * hh.normed_wt
    muni_pop = census.get(m[2:]) # municipality numbers don't include the state prefix in the census
    expected_num_hh = int(round(muni_pop / expected_hh_size))
    print m, ipums_lookup['MUNIMX'][unicode(m)], expected_hh_size, muni_pop, expected_num_hh

    # Locate each household
    for i in range(expected_num_hh):
        if i % 10000 == 0:
            print '\t' + str(i)
        # Sample IPUMS households (weighted by IPUMS authors)
        hh = ipums_data.sample_households(m)

        # Sample pixels (weighted by luminosity)
        px = muni_weights.sample_pixels(Municipality_Weights.lookup_by_ipums_num[m])

        # Sample lat-long uniformly within pixel
        x, y = px.locate_hh()
        #print x, y
        
        # Report household & people data
        hh_output_file.write(' '.join(map(str,[hid,'house',x,y,px.long,px.lat])) + '\n')
        # for each person in the household . . .
        for p in hh.people:
            field_vals = [pid, hid, p.data[age_field], p.data[sex_field], 0, 0, -1, hh.serial, p.data[pernum_field], p.data[empstatd_field]]
            field_vals = map(str,field_vals)
            pop_output_file.write(' '.join(field_vals) + '\n')
            pid += 1
        
        hid += 1

hh_output_file.close()
pop_output_file.close()
