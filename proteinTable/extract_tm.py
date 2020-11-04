import sys
import csv
import requests
from lxml import etree
import re
import pandas as pd
import gzip

proteins = {}

# Test with just one chromosome, then do the final run with the whole proteome
# Input files can be either compressed on uncompressed
input_file = 'nextprot_chromosome_MT.xml'
#input_file = 'G:/Data/nextprot/nextprot_release_2020-01-17/xml/nextprot_chromosome_21.xml.gz'
#input_file = 'G:/Data/nextprot/nextprot_release_2020-01-17/xml/nextprot_all.xml.gz'

# If the file ends with gz, then open it as a gzipped stream
is_compressed = False
if input_file.endswith('.gz'):
    infile = gzip.open(input_file)
    is_compressed = True
else:
    infile = open(input_file, 'r', encoding="utf-8", errors="replace")

buffer = ''
state = 'outside_entry'

for line in infile:

    # If we're reading a compressed stream, each line is actually a byte array, so convert to a string
    if is_compressed:
        line = line.decode('utf-8')

    # If we're not inside an entry and we see the start of a entry, begin buffering
    if state == 'outside_entry' and '<entry accession' in line:
        state = 'in_entry'

    # If we're inside an entry and find the end tag, end the buffering
    elif state == 'in_entry' and '</entry>' in line:
        state = 'parse_entry'
        buffer += line

    # If inside an entry, add this line to the buffer
    if state == 'in_entry':
        buffer += line

    # If we just finished an entry, then parse the XML in the buffer and extract what we need
    if state == 'parse_entry':
        doc = etree.fromstring(buffer)
        for entry in doc.xpath('//entry'):
            counter = 0
            match = re.match(r'NX_(.*)', entry.attrib['accession'])
            identifier = match.group(1)

            # Loop over the CV terms
            cvterms = entry.xpath('.//cv-term[@terminology="nextprot-topology-cv"]/text()')
            for item in cvterms:
                if item == "Transmembrane region":
                    counter += 1
                    
            proteins[identifier] = {}
            proteins[identifier]['Identifier'] = identifier
            proteins[identifier]['n_TMs'] = counter
            print(f"{identifier}\t{counter}")

        # Done parsing this entry, reset the buffer
        state = 'outside_entry'
        buffer = ''

# Convert the dict to a dataframe and write as xlsx
df = pd.DataFrame(proteins)
df_t = df.transpose()
df_t = df_t[['Identifier', 'n_TMs']]
df_t.columns = ['Identifier', 'n_TMRs']
df_t.to_excel('tmr_table.xlsx', index = False)


    
    