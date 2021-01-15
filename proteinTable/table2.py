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
from zipfile import ZipFile
import numpy as np
import mygene


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
    if 'T cell receptor' in proteins[identifier]['Name']:
        if 'variable' in proteins[identifier]['Name'] or 'diversity' in proteins[identifier]['Name'] or 'joining' in proteins[identifier]['Name']:
            proteins[identifier]['protein_category'] = 'TCR_VDJ'
    if 'Taste receptor' in proteins[identifier]['Name']:
        proteins[identifier]['protein_category'] = 'TasteRecep'
    if 'Endogenous retrovirus' in  proteins[identifier]['Name'] or 'endogenous retrovirus' in proteins[identifier]['Name']:
        proteins[identifier]['protein_category'] = 'ERVM'
    if 'PRAME' in proteins[identifier]['Name']:
        proteins[identifier]['protein_category'] = 'PRAME'
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

mg = mygene.MyGeneInfo()
print('INFO: Reading normal_tissue.tsv')
if os.path.isfile('normal_tissue.tsv.zip')==False:
    link = 'https://www.proteinatlas.org/download/normal_tissue.tsv.zip'
    resp = requests.get(link)
    with open('normal_tissue.tsv.zip', 'wb') as f_output:
        f_output.write(resp.content)
    zip = ZipFile('normal_tissue.tsv.zip')
    zip.extractall()
tsv_file2 = open('normal_tissue.tsv')
tissue = {}
read_tsv2 = csv.reader(tsv_file2, delimiter = '\t')
next(read_tsv2)
gene1 = 'TSPAN6'
counter = 0
ensg1 = []
ensg1.append('ENSG00000000003')
for row in read_tsv2:
    gene2 = row[1]
    tissue[gene2] = {}
    tissue[gene2]['reliability'] = row[5]
    if gene1 != gene2:
        ensg1.append(row[0])
        gene1 = gene2
        counter = 0
    if row[4] == 'High':
        counter += 1
    tissue[gene2]['high'] = counter

missing1 = []
uniprot1 = []

ensg_convert1 = mg.querymany(ensg1, scopes = 'ensembl.gene', fields = 'uniprot', species = 'human')
for item in ensg_convert1:
    if 'uniprot' in item:
        if 'Swiss-Prot' in item['uniprot']:
            uniprot1.append((item['uniprot']['Swiss-Prot']))
    else: 
        missing1.append(item['query'])

symbol1 = []
tsv_file2.seek(0)
for row in read_tsv2:
    for missing in missing1:
        if missing == row[0]:
            head, sep, tail = row[1].partition('.')
            x = head
            if x not in symbol1:
                symbol1.append(x)

tsv_file2.close()

counter = 0
ident = []
gene = []
for identifier in proteins:
    ident.append(proteins[identifier]['Identifier'])
    gene.append(proteins[identifier]['Symbol'])
for x in uniprot1:
    if x not in ident:
        counter += 1
for y in symbol1:
    if y not in gene:
        counter += 1
print('Unmatched:', counter)   
print('INFO: Reading rna_tissue_consensus.tsv')
if os.path.isfile('rna_tissue_consensus.tsv.zip')==False:
    link = 'https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip'
    resp = requests.get(link)
    with open('rna_tissue_consensus.tsv.zip', 'wb') as f_output:
        f_output.write(resp.content)
    zip = ZipFile('rna_tissue_consensus.tsv.zip')
    zip.extractall()
tsv_file3 = open('rna_consensus.tsv')
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
            if statistics.median(data) != 0:
                tissue[gene2]['median'] = statistics.median(data)
            else: 
                tissue[gene2]['median'] = ''
            tissue[gene2]['max'] = max(data)
            positive = 0
            for number in data:
                if number > 0:
                    positive += 1
            tissue[gene2]['positive'] = positive
        else:
            tissue[gene2] = {}
            tissue[gene2]['reliability'] = ''
            tissue[gene2]['high'] = ''
            if statistics.median(data) != 0:
                tissue[gene2]['median'] = statistics.median(data)
            else: 
                tissue[gene2]['median'] = ''
            tissue[gene2]['max'] = max(data)
tsv_file3.close()

for identifier in proteins:
    if proteins[identifier]['Symbol'] in tissue:
        proteins[identifier]['HPA_Ab_Reliability'] = tissue[proteins[identifier]['Symbol']]['reliability']
        proteins[identifier]['HPA_nHigh'] = tissue[proteins[identifier]['Symbol']]['high']
        proteins[identifier]['HPAcRNA_gt0Median'] = tissue[proteins[identifier]['Symbol']]['median']
        proteins[identifier]['HPAcRNA_Ngt0'] = tissue[proteins[identifier]['Symbol']]['positive']
        proteins[identifier]['HPAcRNA_Max'] = tissue[proteins[identifier]['Symbol']]['max']

print('No gene symbol:', geneError)
print('No mature proteins:', matureError)

print('INFO: Writing final result: protein_table.xlsx')
df = pd.DataFrame(proteins)
df_t = df.transpose()
df_t = df_t[['Identifier', 'Chr', 'Symbol', 'PE', 'Name', 'n_Isos', 'Length', 'Gravy', 'n_PTMs', 'n_Var', 'Proteomics', 'Ab', '3D', 'Disease', 'n_Tryptic' , 'PA_category', 'PA_n_peptides', 'Edman', 'textbook knowledge', 'SP curated PPI', 'IntAct PPI GOLD', 'Exp function', 'mutagenesis', 'MS_PA', 'MS_MSV', 'MS_nP', 'protein_category', 'n_TMRs', 'HPA_Ab_Reliability', 'HPA_nHigh', 'HPAcRNA_gt0Median', 'HPAcRNA_Ngt0', 'HPAcRNA_Max']]
df_t.columns = ['Identifier', 'Chr', 'Symbol', 'PE', 'Name', 'n_Isos', 'Length', 'Gravy', 'n_PTMs', 'n_Var', 'Proteomics', 'Ab', '3D', 'Disease', 'n_Tryptic' , 'PA_category', 'PA_n_peptides', 'Edman', 'textbook knowledge', 'SP curated PPI', 'IntAct PPI GOLD', 'Exp function', 'mutagenesis', 'MS_PA', 'MS_MSV', 'MS_nP', 'protein_category', 'n_TMRs', 'HPA_Ab_Reliability', 'HPA_nHigh','HPAcRNA_gt0Median', 'HPAcRNA_Ngt0', 'HPAcRNA_Max']
df_t.to_excel('protein_table.xlsx', index = False)

gravy = []
for identifier in proteins:
    if proteins[identifier]['PE'] != 5:
        gravy.append(proteins[identifier]['Gravy'])

gravy2 = []
for identifier in proteins:
    if proteins[identifier]['PE'] == 2 or proteins[identifier]['PE'] ==3 or proteins[identifier]['PE'] ==4:
        gravy2.append(proteins[identifier]['Gravy'])

tryptic = []
for identifier in proteins:
    if proteins[identifier]['PE'] != 5:
        tryptic.append(proteins[identifier]['n_Tryptic'])

tryptic2 = []
for identifier in proteins:
    if proteins[identifier]['PE'] == 2 or proteins[identifier]['PE'] ==3 or proteins[identifier]['PE'] ==4:
        tryptic2.append(proteins[identifier]['n_Tryptic'])

tmr = {}
for identifier in proteins:
    if proteins[identifier]['PE'] != 5:
        tmr_n = proteins[identifier]['n_TMRs']
        if tmr_n not in tmr:
            tmr[tmr_n]= []
        tmr[tmr_n].append(proteins[identifier]['Gravy'])

reliability = {}
reliability['blank'] = 0
for identifier in proteins:
    if proteins[identifier]['PE'] != 5:
        if 'HPA_Ab_Reliability' in proteins[identifier]:
            level = proteins[identifier]['HPA_Ab_Reliability']
            if level == '':
                reliability['blank'] += 1
            else:
                if level not in reliability:
                    reliability[level] = 0
                reliability[level] += 1

reliability_pe = {}
reliability_pe['blank'] = {}
for identifier in proteins:
    if proteins[identifier]['PE'] != 5:
        if 'HPA_Ab_Reliability' in proteins[identifier]:
            level = proteins[identifier]['HPA_Ab_Reliability']
            pe = str(proteins[identifier]['PE'])
            if level == '':
                if pe not in reliability_pe['blank']:
                    reliability_pe['blank'][pe] = 0
                reliability_pe['blank'][pe] += 1
            else:
                if level not in reliability_pe:
                    reliability_pe[level] = {}
                if pe not in reliability_pe[level]:
                    reliability_pe[level][pe] = 0
                reliability_pe[level][pe] += 1
reliability_pe_list = []
reliability_pe_list.append(reliability_pe['Enhanced'])
reliability_pe_list.append(reliability_pe['Supported'])
reliability_pe_list.append(reliability_pe['Approved'])
reliability_pe_list.append(reliability_pe['Uncertain'])
reliability_pe_list.append(reliability_pe['blank'])

if os.path.isdir('Histograms')==False:
    os.mkdir('Histograms')

print('INFO: Creating plots')
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

names = ['Enhanced', 'Supported', 'Approved', 'Uncertain']
values = [reliability['Enhanced'],reliability['Supported'],reliability['Approved'],reliability['Uncertain']]
plt.bar(names, values)
plt.savefig('Histograms/reliability')
plt.show()

category = ['Enhanced', 'Supported', 'Approved', 'Uncertain', 'blank']
pe1 = []
pe2 = []
pe3 = []
pe4 = []

for status in reliability_pe_list:      
    if '1' in status:
        pe1.append(status['1'])
    else: 
        pe1.append(0)
    if '2' in status:
        pe2.append(status['2'])
    else: 
        pe2.append(0)
    if '3' in status:
        pe3.append(status['3'])
    else: 
        pe3.append(0)
    if '4' in status:
        pe4.append(status['4'])
    else: 
        pe4.append(0)

r = range(len(category))
plt.bar(r, pe1, color = '#C1FFC1', label = 'PE1')
plt.bar(r, pe2, bottom = pe1, color = '#FFDEAD', label = 'PE2')
plt.bar(r, pe3, bottom = np.add(pe1,pe2), color = '#FF9999', label = 'PE3')
plt.bar(r, pe4, bottom = np.add(np.add(pe1,pe2), pe3), color = '#CAE1FF', label = 'PE4')
plt.xticks(r, category)
plt.legend()
plt.savefig('Histograms/reliability_pe')
plt.show()

category2 = ['TCR_VDJ', 'TasteRecep', 'ERVM', 'PRAME']
category_pe = {}
for identifier in proteins:
    if proteins[identifier]['protein_category'] in category2:
        cat = proteins[identifier]['protein_category']
        pe = str(proteins[identifier]['PE'])
        if cat not in category_pe:
            category_pe[cat] = {}
        if pe not in category_pe[cat]:
            category_pe[cat][pe] = 0
        category_pe[cat][pe] += 1

pe1 = []
pe2 = []
pe3 = []
pe4 = []
for c in category_pe:      
    if '1' in category_pe[c]:
        pe1.append(category_pe[c]['1'])
    else: 
        pe1.append(0)
    if '2' in category_pe[c]:
        pe2.append(category_pe[c]['2'])
    else: 
        pe2.append(0)
    if '3' in category_pe[c]:
        pe3.append(category_pe[c]['3'])
    else: 
        pe3.append(0)
    if '4' in category_pe[c]:
        pe4.append(category_pe[c]['4'])
    else: 
        pe4.append(0)

r = range(len(category2))
plt.bar(r, pe1, color = '#C1FFC1', label = 'PE1')
plt.bar(r, pe2, bottom = pe1, color = '#FFDEAD', label = 'PE2')
plt.bar(r, pe3, bottom = np.add(pe1,pe2), color = '#FF9999', label = 'PE3')
plt.bar(r, pe4, bottom = np.add(np.add(pe1,pe2), pe3), color = '#CAE1FF', label = 'PE4')
plt.xticks(r, category2)
plt.legend()
plt.savefig('Histograms/category1_pe')
plt.show()

category3 = ('PE2,3,4', 'protein_category', '')
pe_234 = {}
for identifier in proteins:
    if proteins[identifier]['PE']==2 or proteins[identifier]['PE']==3 or proteins[identifier]['PE']==4:
        num = str(proteins[identifier]['PE'])
        if num not in pe_234:
            pe_234[num] = 0
        pe_234[num] += 1
pc = {}
pc['other'] = 0
for identifier in proteins:
    if proteins[identifier]['PE']==2 or proteins[identifier]['PE']==3 or proteins[identifier]['PE']==4:
        p = proteins[identifier]['protein_category']
        if p == '':
            pc['other'] += 1
        else:
            if p not in pc:
                pc[p] = 0
            pc[p] += 1
 
plt.bar(0, pe_234['2'], color = '#FFDEAD', label = 'PE2')
plt.bar(0, pe_234['3'], bottom = pe_234['2'], color = '#FF9999', label = 'PE3')
plt.bar(0, pe_234['4'], bottom = np.add(pe_234['2'],pe_234['3']), color = '#CAE1FF', label = 'PE4')
plt.bar(1, pc['OR'], color = '#B00923', label = 'OR')
plt.bar(1, pc['GPCR'], bottom = pc['OR'], color = '#F09729', label = 'GPCR')
plt.bar(1, pc['defensin'], bottom = np.add(pc['OR'], pc['GPCR']), color = '#3B7D5C', label = 'defensin')
plt.bar(1, pc['zinc finger'], bottom = np.add(np.add(pc['OR'], pc['GPCR']), pc['defensin']), color = '#038DB2', label = 'zinc finger')
plt.bar(1, pc['LINC RNA'], bottom = np.add(np.add(np.add(pc['OR'], pc['GPCR']), pc['defensin']), pc['zinc finger']), color = '#6E6B95', label = 'LINC RNA')
plt.bar(1, pc['TCR_VDJ'], bottom = np.add(np.add(np.add(np.add(pc['OR'], pc['GPCR']), pc['defensin']), pc['zinc finger']), pc['LINC RNA']), color = '#1C2440', label = 'TCR_VDJ')
plt.bar(1, pc['TasteRecep'], bottom = np.add(np.add(np.add(np.add(np.add(pc['OR'], pc['GPCR']), pc['defensin']), pc['zinc finger']), pc['LINC RNA']), pc['TCR_VDJ']), color = '#3C1249', label = 'TasteRecep')
plt.bar(1, pc['ERVM'], bottom = np.add(np.add(np.add(np.add(np.add(np.add(pc['OR'], pc['GPCR']), pc['defensin']), pc['zinc finger']), pc['LINC RNA']), pc['TCR_VDJ']), pc['TasteRecep']), color = '#CBA884', label = 'ERVM')
plt.bar(1, pc['PRAME'], bottom = np.add(np.add(np.add(np.add(np.add(np.add(np.add(pc['OR'], pc['GPCR']), pc['defensin']), pc['zinc finger']), pc['LINC RNA']), pc['TCR_VDJ']), pc['TasteRecep']), pc['ERVM']), color = '#ACA296', label = 'PRAME')
plt.bar(1, pc['hydrophobic'], bottom = np.add(np.add(np.add(np.add(np.add(np.add(np.add(np.add(pc['OR'], pc['GPCR']), pc['defensin']), pc['zinc finger']), pc['LINC RNA']), pc['TCR_VDJ']), pc['TasteRecep']), pc['ERVM']), pc['PRAME']), color = '#678b8b', label = 'hydrophobic')
plt.bar(1, pc['very hydrophobic'], bottom = np.add(np.add(np.add(np.add(np.add(np.add(np.add(np.add(np.add(pc['OR'], pc['GPCR']), pc['defensin']), pc['zinc finger']), pc['LINC RNA']), pc['TCR_VDJ']), pc['TasteRecep']), pc['ERVM']), pc['PRAME']), pc['hydrophobic']), color = '#09311C', label = 'very hydrophobic')
plt.bar(1, pc['other'], bottom = np.add(np.add(np.add(np.add(np.add(np.add(np.add(np.add(np.add(np.add(pc['OR'], pc['GPCR']), pc['defensin']), pc['zinc finger']), pc['LINC RNA']), pc['TCR_VDJ']), pc['TasteRecep']), pc['ERVM']), pc['PRAME']), pc['hydrophobic']), pc['very hydrophobic']), color = '#BE33FF', label = 'other')
plt.bar(2, 0, width = 1.5)

plt.xticks(range(len(category3)), category3)
handles, labels = plt.gca().get_legend_handles_labels()
order = [2,1,0,14,13,12,11,10,9,8,7,6,5,4,3]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc= 'upper right')
plt.savefig('Histograms/category2_pe')
plt.show()