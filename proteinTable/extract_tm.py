import csv
import requests
from lxml import etree
import re
import pandas as pd

proteins = {}
data = etree.parse('nextprot_chromosome_MT.xml')
data = etree.tostring(data)
doc = etree.XML(data)
for att in doc.xpath('//entry'):
    counter = 0
    match = re.match(r'NX_(.*)', att.attrib['accession'])
    identifier = match.group(1)
    list = att.xpath('.//cv-term[@terminology="nextprot-topology-cv"]/text()')
    for item in list:
        if item == "Transmembrane region":
            counter += 1
            
    proteins[identifier] = {}
    proteins[identifier]['Identifier'] = identifier
    proteins[identifier]['n_TMs'] = counter

df = pd.DataFrame(proteins)
df_t = df.transpose()
df_t = df_t[['Identifier', 'n_TMs']]
df_t.columns = ['Identifier', 'n_TMRs']
df_t.to_excel('tmr_table.xlsx', index = False)
        
    
    