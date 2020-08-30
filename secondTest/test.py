import requests
import pandas as pd
import os.path
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from tabulate import tabulate
import textwrap

# if os.path.isfile('test.xls')==False:
#     link = 'https://api.nextprot.org/export/entries.xls?query=*'
#     resp = requests.get(link)
#     with open('test.xls', 'wb') as p_output:
#         p_output.write(resp.content)

# reading = pd.read_excel('test.xls')
# pe_list = reading['PE']

# u_list=[]
# for a in pe_list.index:
#     if pe_list[a] == 'Predicted':
#         u_list.append(reading['acc. code'][a])

if os.path.isfile('nP20k.fasta')==False:
    link = 'http://www.peptideatlas.org/tmp/nP20k.fasta.gz'
    resp = requests.get(link)
    with open('nP20k.fasta', 'wb') as f_output:
        f_output.write(resp.content)

info = []
with open('nP20k.fasta', 'rU') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        if record.description.__contains__('PE=2'):
            analyzed_seq = ProteinAnalysis(str(record.seq))
            tup = (record.id, '2', analyzed_seq.gravy(), textwrap.fill(record.description, 20))
            info.append(tup)
        if record.description.__contains__('PE=3'):
            analyzed_seq = ProteinAnalysis(str(record.seq))
            tup = (record.id, '3', analyzed_seq.gravy(), textwrap.fill(record.description, 20))
            info.append(tup)
        if record.description.__contains__('PE=4'):
            analyzed_seq = ProteinAnalysis(str(record.seq))
            tup = (record.id, '4', analyzed_seq.gravy(), textwrap.fill(record.description, 20))
            info.append(tup)


print(tabulate(info, headers = ['Identifier', 'PE', 'GRAVY', 'Description']))
        