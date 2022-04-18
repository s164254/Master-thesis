from curses.ascii import isdigit
from os import path
import pandas as pd
import csv
import math

def unique_peptides_valid_value(value):
    return (value != '1') and (isdigit(value[0]) or value[0]=='.')

def unique_peptides_valid_row(row_data):
    return any( (1 for v in row_data if unique_peptides_valid_value(v)))

def label_free_quant_transform(value, unique_peptide_value):
    return unique_peptide_value != '1' and math.log2(float(value)) or '0'

LABEL_FREE_QUANT = 'Label-Free Quant'
UNIQUE_PEPTIDES_COL_PREFIX = 'UniquePeptides'
SKIP_COLUMNS = ('PG.ProteinDescriptions',)

GENE_LOOKUP_FILE='gene_lookup.txt'
MISSING_GENE_CODES_FILE='missing_gene_codes.txt'
GENE_CODE_COLUMN_NAME = 'PG.Genes'

script_dir = path.dirname(path.realpath(__file__))
fname = 'Batch2_data.csv'

gene_dict = {}
if path.exists(path.join(script_dir, GENE_LOOKUP_FILE)):
    with open(path.join(script_dir, GENE_LOOKUP_FILE), 'r') as f:
        gene_dict = dict([(c[0],1) for c in f.readlines().split()])


with open(path.join(script_dir, fname)) as csvfile:
    rdr = csv.reader(csvfile, delimiter=',')
    rows = [tuple(row) for row in rdr]

header = rows[0]
data = rows[1:]

unique_peptides_idxs = [i for i, c in enumerate(header) if c.find(UNIQUE_PEPTIDES_COL_PREFIX) > 0]
idx1 = unique_peptides_idxs[0]
idx2 = unique_peptides_idxs[-1]
valid_rows = [row for row in data if unique_peptides_valid_row(row[idx1:idx2])]

gene_code_column_idx = [i for i, c in enumerate(header) if c.find(GENE_CODE_COLUMN_NAME) >= 0][0]
missing_gene_codes = [row[gene_code_column_idx] for row in valid_rows if not row[gene_code_column_idx] in gene_dict]
if missing_gene_codes:
    with open(path.join(script_dir, MISSING_GENE_CODES_FILE), 'w') as f:
        f.write('\n'.join(missing_gene_codes))

label_free_quant_idxs = [i for i, c in enumerate(header) if c.find(LABEL_FREE_QUANT) >= 0]
unique_peptides_adder = unique_peptides_idxs[0] - label_free_quant_idxs[0]
idx1 = label_free_quant_idxs[0]
idx2 = label_free_quant_idxs[-1] + 1 

label_free_quant_transformed = [[label_free_quant_transform(row[lfq_idx], row[lfq_idx+unique_peptides_adder]) for lfq_idx in range(idx1,idx2)] for row in valid_rows]

l = len(valid_rows)