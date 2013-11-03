#!/usr/bin/python
import xml.parsers.expat

class IPUMS_XML_Parser:
    vartags = list()
    taglist = list()
    data = dict()

    @staticmethod
    def start_element(name, attrs):
        IPUMS_XML_Parser.taglist.append(name)
        if name == 'var':
            #print 'starting var:', attrs['ID']
            IPUMS_XML_Parser.vartags.append({'id':attrs['ID'],'categories':[]})
            open_var = True
            
    @staticmethod
    def end_element(name):
        if (len(IPUMS_XML_Parser.taglist) > 0):
            IPUMS_XML_Parser.taglist.pop()

    @staticmethod
    def char_data(xml_data):
        if len(IPUMS_XML_Parser.taglist) > 2 and IPUMS_XML_Parser.taglist[-2] == "catgry":
            if IPUMS_XML_Parser.taglist[-1] == "catValu":
                #IPUMS_XML_Parser.vartags[-1]['categories'].append({'value':repr(xml_data)})
                IPUMS_XML_Parser.vartags[-1]['categories'].append({'value':xml_data})
            elif IPUMS_XML_Parser.taglist[-1] == "labl":
                #IPUMS_XML_Parser.vartags[-1]['categories'][-1]['label'] = repr(xml_data)
                IPUMS_XML_Parser.vartags[-1]['categories'][-1]['label'] = xml_data

    def __init__(self, filename):

        p = xml.parsers.expat.ParserCreate()

        p.StartElementHandler  = IPUMS_XML_Parser.start_element
        p.EndElementHandler    = IPUMS_XML_Parser.end_element
        p.CharacterDataHandler = IPUMS_XML_Parser.char_data

        xml_data = file(filename).read()
        p.Parse(xml_data)

        for tag in IPUMS_XML_Parser.vartags:
            IPUMS_XML_Parser.data[tag['id']] = dict()
            # IPUMS_XML_Parser.data['CNTRY']['484']
            if len(tag['categories']) > 0:
                for cat in tag['categories']:
                    IPUMS_XML_Parser.data[tag['id']][cat['value']] = cat['label']


#  Usage:
# IPUMS_XML_Parser('./raw_data/ipums/ipumsi_00005.xml')
# ipums = IPUMS_XML_Parser.data

# for c in ipums['CNTRY']:
#     print c, ipums['CNTRY'][c]
#print IPUMS_XML_Parser.vartags       
