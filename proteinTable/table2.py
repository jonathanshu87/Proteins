import requests
import pandas as pd
import os.path
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import re 
from pyteomics import parser
import csv 
import xlrd
import matplotlib.pyplot as plt
import statistics


proteins={}
geneError = 0
matureError = 0
with open ('nextprot_all.peff') as file:
    print('INFO: Reading nextprot_all.peff')
    for record in SeqIO.parse(file, 'fasta'):
        trypnum = 0
        recordstr = str(record.seq)
        fixed_sequence = re.sub(r'U','G',recordstr)           
        analyzed_seq = ProteinAnalysis(fixed_sequence)
        match = re.match(r'nxp:NX_(.+?)-(\d+)', record.description)
        identifier = match.group(1)
        match2 = re.search(r'\\PName=(.+?) isoform (.*?) (.*?) \\', record.description)
        if match2:
            name = match2.group(1)
        else:
            name = ''
            print(f'WARNING: Unable to find isoform information in {record.description}')
        match3 = re.search(r'\\GName=(.+?) \\N', record.description)
        if match3:
            gene = match3.group(1)
        else:
            gene = ''
            geneError += 1
        match4 = re.search(r'\\PE=(.+?) \\', record.description)
        if match4:
            pe = match4.group(1)
        else:
            pe = ''
            print(f'WARNING: Unable to find PE value in {record.description}')
        match5 = re.search(r'\((\d+)\|(\d+)\|PEFF:\d+\|mature protein\)', record.description)
        if match5:
            start = int(match5.group(1))
            end = int(match5.group(2))
            mature = fixed_sequence[start:(end+1)]
            peptides = parser.cleave(mature, 'trypsin')
            for peptide in peptides:
                if 9 <= len(peptide) <= 30:
                    trypnum += 1
        else:
            matureError += 1
        if identifier not in proteins:
            proteins[identifier] = {}
            proteins[identifier]['Identifier'] = identifier
            proteins[identifier]['Name'] = name
            proteins[identifier]['Symbol'] = gene
            proteins[identifier]['PE'] = int(pe)
            proteins[identifier]['Length'] = len(record.seq)
            proteins[identifier]['Gravy'] = analyzed_seq.gravy()
            proteins[identifier]['Chr'] = ''
            proteins[identifier]['n_Tryptic'] = trypnum

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
                    if match2 == False and match3 == False:
                        chromosome = ''
                        print(f'WARNING: Unable to find chromosome information in {file}')
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
                    proteins[identifier]['n_Isos'] = int(isoforms)
                    proteins[identifier]['n_Var'] = int(variants)
                    proteins[identifier]['n_PTMs'] = int(ptms)

for ident in proteins:
    if proteins[ident]['Chr'].isdigit():
        proteins[ident]['Chr'] = int(proteins[ident]['Chr'])

if os.path.isfile('peptideatlas.tsv')==False:
    link = 'https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetNextProtChromMapping?atlas_build_id=491&nextprot_mapping_id=55&apply_action=QUERY&output_mode=tsv'
    resp = requests.get(link)
    with open('peptideatlas.tsv', 'wb') as f_output:
        f_output.write(resp.content)
tsv_file = open('peptideatlas.tsv')
print('INFO: Reading peptideatlas.tsv')
read_tsv = csv.reader(tsv_file, delimiter = '\t')
next(read_tsv)
for row in read_tsv:
    identifier = row[1]
    if identifier in proteins:
        proteins[identifier]['PA_category'] = row[9]
        proteins[identifier]['PA_n_peptides'] = int(row[10])
tsv_file.close()

for identifier in proteins:
    proteins[identifier]['Edman'] = ''
    proteins[identifier]['textbook knowledge'] = ''
    proteins[identifier]['SP curated PPI'] = ''
    proteins[identifier]['IntAct PPI GOLD'] = ''
    proteins[identifier]['Exp function'] = ''
    proteins[identifier]['mutagenesis'] = ''
    proteins[identifier]['protein_category'] = ''

xls = pd.read_excel('Master_table_for_HPP_to_send[1].xls')
print('INFO: Reading Master_table_for_HPP_to_send[1].xls')
for index, row in xls.iterrows():
    match1 = re.search(r'NX_(.+)', row['acc. code'])
    if match1:
        identifier = match1.group(1)
    proteins[identifier]['Edman'] = row['Edman']
    proteins[identifier]['textbook knowledge'] = row['textbook knowledge']
    proteins[identifier]['SP curated PPI'] = row['SP curated PPI']
    proteins[identifier]['IntAct PPI GOLD'] = row['IntAct PPI GOLD']
    proteins[identifier]['Exp function'] = row['Exp function']
    proteins[identifier]['mutagenesis'] = row['mutagenesis']

xlsx = pd.read_excel('PE1_forGil.xlsx', sheet_name = [0,1,7], header = None)
print('INFO: Reading PE1_forGil.xlsx')
for sheet in xlsx:
    for line in xlsx[sheet]:
        for item in xlsx[sheet][line]:
            match1 = re.search(r'NX_(.+)', item)
            if match1:
                identifier = match1.group(1)
            if sheet == 0 and identifier in proteins:
                proteins[identifier]['MS_PA'] = 'Y'
            if sheet == 1 and identifier in proteins:
                proteins[identifier]['MS_MSV'] = 'Y'
            if sheet == 7 and identifier in proteins:
                 proteins[identifier]['MS_nP'] = 'Y'

for identifier in proteins:
    if 'Olfactory receptor' in proteins[identifier]['Name'] or 'olfactory receptor' in proteins[identifier]['Name']:
        proteins[identifier]['protein_category'] = 'OR'
    if 'GPCR' in proteins[identifier]['Name'] or 'G-protein coupled receptor' in proteins[identifier]['Name']:
        proteins[identifier]['protein_category'] = 'GPCR'
    if 'defensin' in proteins[identifier]['Name']:
        proteins[identifier]['protein_category'] = 'defensin'
    if 'zinc finger' in proteins[identifier]['Name'] or 'Zinc finger' in proteins[identifier]['Name']:
        proteins[identifier]['protein_category'] = 'zinc finger'
    if 'LINC' in proteins[identifier]['Symbol']:
        proteins[identifier]['protein_category'] = 'LINC RNA'
    if proteins[identifier]['Gravy'] >= 0.5 and proteins[identifier]['protein_category'] == '':
        proteins[identifier]['protein_category'] = 'very hydrophobic'
    if proteins[identifier]['Gravy'] >= 0.0 and proteins[identifier]['protein_category'] == '':
        proteins[identifier]['protein_category'] = 'hydrophobic'

print('INFO: Reading tmr_table.xlsx')
loc = ("tmr_table.xlsx")
wb = xlrd.open_workbook(loc)
sheet = wb.sheet_by_index(0)
for i in range(sheet.nrows):
    if i != 0:
        identifier = (sheet.cell_value(i,0))
        proteins[identifier]['n_TMRs'] = (sheet.cell_value(i,1))

print('INFO: Reading normal_tissue.tsv')
tissue = {}
tsv_file2 = open('normal_tissue.tsv')
read_tsv2 = csv.reader(tsv_file2, delimiter = '\t')
next(read_tsv2)
gene1 = 'TSPAN6'
counter = 0
for row in read_tsv2:
    gene2 = row[1]
    tissue[gene2] = {}
    tissue[gene2]['reliability'] = row[5]
    if gene1 != gene2:
        gene1 = gene2
        counter = 0
    if row[4] == 'High':
        counter += 1
    tissue[gene2]['high'] = counter
tsv_file2.close()

print('INFO: Reading rna_tissue_consensus.tsv')
tsv_file3 = open('rna_tissue_consensus.tsv')
read_tsv3 = csv.reader(tsv_file3, delimiter = '\t')
next(read_tsv3)
gene1 = 'TSPAN6'
data = []
for row in read_tsv3:
    gene2 = row[1]
    data.append(float(row[3]))
    if gene1 != gene2:
        gene1 = gene2
        data = []
    if len(data) > 0:
        if gene2 in tissue:
            tissue[gene2]['median'] = statistics.median(data)
            tissue[gene2]['max'] = max(data)
        else:
            tissue[gene2] = {}
            tissue[gene2]['reliability'] = ''
            tissue[gene2]['high'] = ''
            tissue[gene2]['median'] = statistics.median(data)
            tissue[gene2]['max'] = max(data)
tsv_file3.close()

for identifier in proteins:
    if proteins[identifier]['Symbol'] in tissue:
        proteins[identifier]['HPA_Ab_Reliability'] = tissue[proteins[identifier]['Symbol']]['reliability']
        proteins[identifier]['HPA_nHigh'] = tissue[proteins[identifier]['Symbol']]['high']
        proteins[identifier]['HPAcRNA_Median'] = tissue[proteins[identifier]['Symbol']]['median']
        proteins[identifier]['HPAcRNA_Max'] = tissue[proteins[identifier]['Symbol']]['max']

print('No gene symbol:', geneError)
print('No mature proteins:', matureError)

print('INFO: Writing final result: protein_table.xlsx')
df = pd.DataFrame(proteins)
df_t = df.transpose()
df_t = df_t[['Identifier', 'Chr', 'Symbol', 'PE', 'Name', 'n_Isos', 'Length', 'Gravy', 'n_PTMs', 'n_Var', 'Proteomics', 'Ab', '3D', 'Disease', 'n_Tryptic' , 'PA_category', 'PA_n_peptides', 'Edman', 'textbook knowledge', 'SP curated PPI', 'IntAct PPI GOLD', 'Exp function', 'mutagenesis', 'MS_PA', 'MS_MSV', 'MS_nP', 'protein_category', 'n_TMRs', 'HPA_Ab_Reliability', 'HPA_nHigh', 'HPAcRNA_Median', 'HPAcRNA_Max']]
df_t.columns = ['Identifier', 'Chr', 'Symbol', 'PE', 'Name', 'n_Isos', 'Length', 'Gravy', 'n_PTMs', 'n_Var', 'Proteomics', 'Ab', '3D', 'Disease', 'n_Tryptic' , 'PA_category', 'PA_n_peptides', 'Edman', 'textbook knowledge', 'SP curated PPI', 'IntAct PPI GOLD', 'Exp function', 'mutagenesis', 'MS_PA', 'MS_MSV', 'MS_nP', 'protein_category', 'n_TMRs', 'HPA_Ab_Reliability', 'HPA_nHigh','HPAcRNA_Median', 'HPAcRNA_Max']
df_t.to_excel('protein_table.xlsx', index = False)

gravy = []
for identifier in proteins:
    gravy.append(proteins[identifier]['Gravy'])

gravy2 = []
for identifier in proteins:
    if proteins[identifier]['PE'] == 2 or proteins[identifier]['PE'] ==3 or proteins[identifier]['PE'] ==4:
        gravy2.append(proteins[identifier]['Gravy'])

tryptic = []
for identifier in proteins:
    tryptic.append(proteins[identifier]['n_Tryptic'])

tryptic2 = []
for identifier in proteins:
    if proteins[identifier]['PE'] == 2 or proteins[identifier]['PE'] ==3 or proteins[identifier]['PE'] ==4:
        tryptic2.append(proteins[identifier]['n_Tryptic'])

tmr = {}
for identifier in proteins:
    tmr_n = proteins[identifier]['n_TMRs']
    if tmr_n not in tmr:
        tmr[tmr_n]= []
    tmr[tmr_n].append(proteins[identifier]['Gravy'])

if os.path.isdir('Histograms')==False:
    os.mkdir('Histograms')

print('INFO: Creating histograms')
min = -2
max = 2
binsize = 0.1
n_bins = int((max-min) / binsize)
plt.hist(gravy, n_bins, [min,max], density = False, facecolor = 'r', alpha = 0.5)
plt.title('Distribution of hydrophobicity of proteins')
plt.xlabel('Hydrophobicity (GRAVY score)')
plt.ylabel('Proteins')
plt.grid(True)
plt.savefig('Histograms/gravy.png')
plt.show()

plt.hist(gravy2, n_bins, [min,max], density = False, facecolor = 'b', alpha = 0.5)
plt.title('Distribution of hydrophobicity of PE 2, 3, 4 proteins')
plt.xlabel('Hydrophobicity (GRAVY score)')
plt.ylabel('Proteins')
plt.grid(True)
plt.savefig('Histograms/gravy2.png')
plt.show()

plt.hist(tryptic, 50, [0,50], density = False, facecolor = 'r', alpha = 0.5)
plt.grid(True)
plt.title('Distribution of tryptic peptides')
plt.xlabel('Number of peptides')
plt.ylabel('Proteins')
plt.savefig('Histograms/tryptic.png')
plt.show()
                        
plt.hist(tryptic2, 50, [0,50], density = False, facecolor = 'b', alpha = 0.5)
plt.grid(True)
plt.title('Distribution of tryptic peptides for PE 2, 3, 4 proteins')
plt.xlabel('Number of peptides')
plt.ylabel('Proteins')
plt.savefig('Histograms/tryptic2.png')
plt.show()

index = []
data = []
for i_tmr in range(16):
    if i_tmr in tmr:
        data.append(tmr[i_tmr])
        index.append(i_tmr)

plt.violinplot(data, index, showmedians = True, showextrema = False)
plt.xlabel('n_TMRs')
plt.ylabel('GRAVY')
plt.ylim(-2,2)
plt.savefig('Histograms/tmr.png')
plt.show()

        

