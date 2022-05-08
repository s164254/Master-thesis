from math import isnan
from os import path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import fileutils as ft

def get_column_range(item_attrs,attr_name,df):
    attr_range = [x for x in item_attrs if x[2]==attr_name]
    return df.columns[attr_range[0][0]:attr_range[-1][0]]

def get_columns(item_attrs,attr_name):
    #return [df.columns[x[0]] for x in item_attrs if x[2]==attr_name]
    return [x[3] for x in item_attrs if x[2]==attr_name]

def is_unique_peptides_nan(value):
    return value == 0 or value == 1 or isnan(value)

# select if any col is 1
ATTR_RE = '_([a-zA-Z0-9]+)\.raw\.PG\.(.+)'
UNIQUE_PEPTIDES = 'UniquePeptides'
LABEL_FREE_QUANT = 'Label-Free Quant'

script_dir = ft.get_script_dir(__file__)
fname = 'Batch2_data.csv'
csv=pd.read_csv(path.join(script_dir, fname))

# create (column index, sample name, attr name) list
item_attrs = [(i, m.groups()[0], m.groups()[1], c) for i, m, c in [(i, re.search(ATTR_RE, c), c) for i, c in enumerate(csv.columns)] if m and len(m.groups()) == 2]

# get list of unique peptides column names to be used in filter column range selection  
unique_peptides_col_names = get_columns(item_attrs,UNIQUE_PEPTIDES)

# remove rows where all unique peptides values are either 0, 1 or NAN
filtered = csv[(csv[unique_peptides_col_names].applymap(lambda value: not is_unique_peptides_nan(value)).any(1))]

# for the remaining rows set label-free quant to 0 if the corresponding unique peptides value is either 0, 1 or NAN
abundance_col_names = get_columns(item_attrs,LABEL_FREE_QUANT)
for abundance_col_name,  unique_peptides_col_name in zip( abundance_col_names, unique_peptides_col_names):
    filtered[abundance_col_name] = filtered.apply(lambda x: 0 if is_unique_peptides_nan(x[unique_peptides_col_name]) else x[abundance_col_name], axis=1)
l=len(filtered)
