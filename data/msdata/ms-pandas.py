from math import isnan
from os import path
from tkinter import font
import numpy as np
import pandas as pd
import re
import fileutils as ft
import matplotlib.pyplot as plt

def get_column_names(item_attrs,attr_name):
    return [x[3] for x in item_attrs if x[2]==attr_name]

#def filter_nlargest(df,column_name,n):
#    return df[df[column_name].isin(df.nlargest(,coln).index)]

def get_sample_names(item_attrs):
    return list(set([x[1] for x in item_attrs]))

def get_column_name(item_attrs,sample_name,attr_name):
    return [x[3] for x in item_attrs if x[2]==attr_name and x[1]==sample_name][0]

def is_unique_peptides_nan(value):
    return value == 0 or value == 1 or isnan(value)

def plot_hist(df, sample, plot_title, x_label, y_label):
    sort_indices = m.sort_indices(sample, LABEL_FREE_QUANT, True)[:10]
    values = m.get_values(sample, LABEL_FREE_QUANT, sort_indices)
    values = values / np.max(values)

    # plot
    fig, ax = plt.subplots()
    ax.bar(m.get_prot(sort_indices),
           values,
           width=1,
           edgecolor="white",
           linewidth=0.7)
    plt.title(plot_title % (sample,))
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()

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
unique_peptides_col_names = get_column_names(item_attrs,UNIQUE_PEPTIDES)

# remove rows where all unique peptides values are either 0, 1 or NAN
filtered = csv[(csv[unique_peptides_col_names].applymap(lambda value: not is_unique_peptides_nan(value)).any(1))]

# put newlines in protein description
filtered[filtered.columns[2]] = filtered.apply(lambda x: re.sub('\s+', '\n', x[2]), axis=1)

# for the remaining rows set label-free quant to 0 if the corresponding unique peptides value is either 0, 1 or NAN
abundance_col_names = get_column_names(item_attrs,LABEL_FREE_QUANT)
for abundance_col_name,  unique_peptides_col_name in zip( abundance_col_names, unique_peptides_col_names):
    filtered[abundance_col_name] = filtered.apply(lambda x: 0 if is_unique_peptides_nan(x[unique_peptides_col_name]) else x[abundance_col_name], axis=1)

for sample_abundance in [get_column_name(item_attrs,sample_name,LABEL_FREE_QUANT) for sample_name in get_sample_names(item_attrs)]:
    #filtered[sample_abundance].sort_values(ascending=0)[:10].plot(x=filtered.columns[2], y=sample_abundance, kind='bar', title=sample_abundance)
    ax = filtered.nlargest(10, sample_abundance).set_index(filtered.columns[2]).plot(y=sample_abundance, kind='bar', rot=0)
    #plot_df = filter_nlargest(filtered, sample_abundance, 10)
    #plot_df.plot(x=filtered.columns[2], y=sample_abundance, kind='bar', title=sample_abundance)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    #plt.title(label=sample_abundance,fontsize=8)
    plt.show()
    l=len(filtered)
