import requests
import pandas as pd
import os.path
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import re 
from pyteomics import parser

proteins={}
with open ('nextprot_all.peff') as file:
    print('INFO: Reading nextprot_all.peff')
    for record in SeqIO.parse(file, 'fasta'):
        trypnum = 0
        recordstr = str(record.seq)
        fixed_sequence = re.sub(r'U','G',recordstr)           
        analyzed_seq = ProteinAnalysis(fixed_sequence)
        match = re.match(r'nxp:NX_(.+?)-(\d+)', record.description)
        identifier = match.group(1)
        match2 = re.search(r'\\PName=(.+?) isoform Iso (\d+) \\', record.description)
        if match2:
            name = match2.group(1)
        match3 = re.search(r'\\GName=(.+?) \\N', record.description)
        if match3:
            gene = match3.group(1)
        match4 = re.search(r'\\PE=(.+?) \\', record.description)
        if match4:
            pe = match4.group(1)
        match5 = re.search(r'\((\d+)\|(\d+)\|PEFF:\d+\|mature protein\)', record.description)
        if match5:
            start = int(match5.group(1))
            end = int(match5.group(2))
            mature = fixed_sequence[start:(end+1)]
            peptides = parser.cleave(mature, 'trypsin')
            for peptide in peptides:
                if 9 <= len(peptide) <= 30:
                    trypnum += 1
        if identifier not in proteins:
            proteins[identifier] = {}
            proteins[identifier]['Identifier'] = identifier
            proteins[identifier]['Name'] = name
            proteins[identifier]['Symbol'] = gene
            proteins[identifier]['PE'] = pe
            proteins[identifier]['Length'] = len(record.seq)
            proteins[identifier]['Gravy'] = analyzed_seq.gravy()
            proteins[identifier]['Chr'] = ''
            proteins[identifier]['Trypsin'] = trypnum
dirs = os.listdir()
chromosome = ''
print('INFO: Reading chromosome*.txt files')
for file in dirs:
    chrom = re.search(r'chromosome_(.+)\.txt', file)
    if chrom:
        with open (file) as infile:
            for line in infile:
                if line.__contains__('NX_'):
                    string = line[28:34] + ' ' + line[88:130]
                    match1 = re.search(r'NX_(.+?) ', line)
                    match2 = re.search(r'(\d+)', string)
                    match3 = re.match(r'([XYMu])', string)
                    info = string.split()
                    identifier = match1.group(1)
                    if match2:
                        match_n = match2.group(1)
                        chromosome = match_n
                    if match3:
                        match_l = match3.group(1)
                        if match_l[0] == 'u':
                            match_l= match_l.replace('u', '?')
                        chromosome = match_l
                    proteomics = info[1]
                    antibody = info[2]
                    three_d = info[3]
                    disease = info[4]
                    isoforms = info[5]
                    variants = info[6]
                    ptms = info[7] 
                    if proteins[identifier]['Chr'] != '' and proteins[identifier]['Chr'].split(',').__contains__(chromosome)==False:
                        chromosome = proteins[identifier]['Chr'] + ',' + chromosome
                    elif len(proteins[identifier]['Chr'].split(','))>=2:
                        chromosome = proteins[identifier]['Chr']
                    proteins[identifier]['Chr'] = chromosome
                    proteins[identifier]['Proteomics'] = proteomics
                    proteins[identifier]['Ab'] = antibody
                    proteins[identifier]['3D'] = three_d
                    proteins[identifier]['Disease'] = disease
                    proteins[identifier]['n_Isos'] = isoforms
                    proteins[identifier]['n_Var'] = variants
                    proteins[identifier]['n_PTMs'] = ptms

print('INFO: Writing final result: protein_table.xlsx')
df = pd.DataFrame.from_dict(proteins)
df_t = df.T
df_t = df_t[['Identifier', 'Chr', 'Symbol', 'PE', 'Name', 'n_Isos', 'Length', 'Gravy', 'n_PTMs', 'n_Var', 'Proteomics', 'Ab', '3D', 'Disease', 'Trypsin']]
df_t.to_excel('protein_table.xlsx')

                        
                        
        

