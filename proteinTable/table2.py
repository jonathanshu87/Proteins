import requests
import pandas as pd
import os.path
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import re 

proteins={}
with open ('nextprot_all.peff') as file:
    for record in SeqIO.parse(file, 'fasta'):
        recordstr = str(record.seq)
        fixed_sequence = re.sub(r'U','G',recordstr)
        analyzed_seq = ProteinAnalysis(fixed_sequence)
        match = re.match(r'nxp:NX_(.+?)-(\d+)', record.description)
        identifier = match.group(1)
        match2 = re.search(r'\\PName=(.+?) isoform Iso (\d+) \\', record.description)
        if match2:
            name = match2.group(1)
            isoform = match2.group(2)
        match3 = re.search(r'\\GName=(.+?) \\N', record.description)
        if match3:
            gene = match3.group(1)
        match4 = re.search(r'\\PE=(.+?) \\', record.description)
        if match4:
            pe = match4.group(1)
        if identifier not in proteins:
            proteins[identifier] = {}
            proteins[identifier]['Isoform_number'] = isoform
            proteins[identifier]['Name'] = name
            proteins[identifier]['Gene'] = gene
            proteins[identifier]['PE_number'] = pe
            proteins[identifier]['Length'] = len(record.seq)
            proteins[identifier]['Gravy'] = analyzed_seq.gravy()

dirs = os.listdir()
for file in dirs:
    chrom = re.search(r'chromosome_(.+)\.txt', file)
    if chrom:
        with open (file) as infile:
            for line in infile:
                if line.__contains__('NX_'):
                    string = line[28:34] + ' ' + line[88:130]
                    match = re.search(r'NX_(.+?) ', line)
                    identifier = match.group(1)
                    info = string.split()
                    chromosome = info[0][0:1]
                    proteomics = info[1]
                    antibody = info[2]
                    three_d = info[3]
                    disease = info[4]
                    isoforms = info[5]
                    variants = info[6]
                    ptms = info[7] 
                    proteins[identifier]['Chromosome'] = chromosome
                    proteins[identifier]['Proteomics'] = proteomics
                    proteins[identifier]['Antibody'] = antibody
                    proteins[identifier]['3D'] = three_d
                    proteins[identifier]['Disease'] = disease
                    proteins[identifier]['Isoforms'] = isoforms
                    proteins[identifier]['Variants'] = variants
                    proteins[identifier]['PTMs'] = ptms
df = pd.DataFrame.from_dict(proteins)
df_t = df.T
df_t.to_excel('protein_table.xlsx')

                        
                        
        

