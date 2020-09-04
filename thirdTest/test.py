import requests
import pandas as pd
import os.path
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import textwrap
import matplotlib.pyplot as plt

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
        recordstr = str(record.seq)
        if recordstr.__contains__('U'):
            recordstr = recordstr.replace('U', 'L')
        if record.description.__contains__('PE=2'):
            analyzed_seq = ProteinAnalysis(recordstr)
            tup = (record.id, '1', analyzed_seq.gravy(), textwrap.fill(record.description, 20))
            info.append(tup)
        if record.description.__contains__('PE=3'):
            analyzed_seq = ProteinAnalysis(recordstr)
            tup = (record.id, '3', analyzed_seq.gravy(), textwrap.fill(record.description, 20))
            info.append(tup)
        if record.description.__contains__('PE=4'):
            analyzed_seq = ProteinAnalysis(recordstr)
            tup = (record.id, '4', analyzed_seq.gravy(), textwrap.fill(record.description, 20))
            info.append(tup)

column_name_list = [ 'Identifier', 'PE', 'GRAVY', 'Description' ]
proteins = pd.DataFrame(info, columns = column_name_list)
proteins.to_csv('PE1_proteins.tsv',sep='\t')


min = -2
max = 2
binsize = 0.1
n_bins = int((max-min) / binsize)
count, x_floor, patches = plt.hist(proteins['GRAVY'], n_bins, [min,max], density = False, facecolor = 'r', alpha = 0.5)
plt.title('Distribution of hydrophobicity of PE 2,3,4 proteins')
plt.xlabel('Hydrophobicity (GRAVY score)')
plt.ylabel('Proteins')
#plt.xlim(min,2)
#plt.ylim(0,200)
plt.grid(True)
plt.show()


        