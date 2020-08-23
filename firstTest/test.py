import requests
import pandas as pd
import os.path
import Bio.SeqIO

if os.path.isfile('test.xls')==False:
    link = 'https://api.nextprot.org/export/entries.xls?query=*'
    resp = requests.get(link)
    with open('test.xls', 'wb') as p_output:
        p_output.write(resp.content)


reading = pd.read_excel('test.xls')
pe_list = reading['PE']

u_list=[]
for a in pe_list.index:
    if pe_list[a] == 'Predicted':
        u_list.append(reading['acc. code'][a])
print(u_list)

if os.path.isfile('test.fasta')==False:
    link = 'https://api.nextprot.org/export/entries.fasta?query=*'
    resp = requests.get(link)
    with open('test.fasta', 'wb') as f_output:
        f_output.write(resp.content)

with open('test.fasta') as file:
    for values in Bio.SeqIO.FastaIO.SimpleFastaParser(file):
        for u in u_list:
            if values[0].__contains__(u):
                print(values)


# list = []
# for a in pe_list.index:
#     info = {}
#     if pe_list[a] == 'Predicted':
#         for column in reading:
#             info[column]=reading[column][a]
#         list.append(info)
# print(list)
