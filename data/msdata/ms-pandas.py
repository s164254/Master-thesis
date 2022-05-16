from math import isnan
from os import path
from tkinter import font
from matplotlib import legend
import numpy as np
import pandas as pd
import re
import fileutils as ft
import matplotlib.pyplot as plt


def get_columns_by_attr_name(item_attrs, attr_name):
    return [(x[1], x[3]) for x in item_attrs if x[2] == attr_name]


def get_column_names(item_attrs, attr_name, sample_names=[]):
    return [
        x[3] for x in item_attrs if x[2] == attr_name and (
            (len(sample_names) == 0) or (x[1] in sample_names))
    ]


def is_unique_peptides_nan(value):
    return value == 0 or value == 1 or isnan(value)


# select if any col is 1
ATTR_RE = '_([a-zA-Z0-9]+)\.raw\.PG\.(.+)'
UNIQUE_PEPTIDES = 'UniquePeptides'
LABEL_FREE_QUANT = 'Label-Free Quant'
PG_PROTEINDESCRIPTIONS = 'PG.ProteinDescriptions'

script_dir = ft.get_script_dir(__file__)
fname = 'Batch2_data.csv'
csv = pd.read_csv(path.join(script_dir, fname))

# create (column index, sample name, attr name) list
item_attrs = [(i, m.groups()[0], m.groups()[1], c)
              for i, m, c in [(i, re.search(ATTR_RE, c), c)
                              for i, c in enumerate(csv.columns)]
              if m and len(m.groups()) == 2]

# get list of unique peptides column names to be used in filter column range selection
unique_peptides_col_names = get_column_names(item_attrs, UNIQUE_PEPTIDES)

# remove rows where all unique peptides values are either 0, 1 or NAN
filtered = csv[(csv[unique_peptides_col_names].applymap(
    lambda value: not is_unique_peptides_nan(value)).any(1))]

# put newlines in protein description
filtered[filtered.columns[2]] = filtered.apply(
    lambda x: re.sub('\s+', '\n', x[2]), axis=1)

# for the remaining rows set label-free quant to 0 if the corresponding unique peptides value is either 0, 1 or NAN
abundance_col_names = get_column_names(item_attrs, LABEL_FREE_QUANT)
for abundance_col_name, unique_peptides_col_name in zip(
        abundance_col_names, unique_peptides_col_names):
    filtered[abundance_col_name] = filtered.apply(
        lambda x: 0 if is_unique_peptides_nan(x[unique_peptides_col_name]
                                              ) else x[abundance_col_name],
        axis=1)

if False:
    # plot 10 largest abundance value with x axis label = protein normalized to 1
    for sample_name, abundance_sample_colname in get_columns_by_attr_name(
            item_attrs, LABEL_FREE_QUANT):  # loop over each sample
        # normalize abundance sample column
        col_max = filtered[abundance_sample_colname].max()
        filtered[abundance_sample_colname] = filtered[
            abundance_sample_colname].div(col_max)

        # take the 10 largest protein description as x axis label and make the bar plot
        chrt = filtered.nlargest(10, abundance_sample_colname).set_index(
            filtered.columns[2]).plot(y=abundance_sample_colname,
                                      kind='bar',
                                      rot=0,
                                      legend=False)

        # use a smaller font size than the default
        plt.xticks(fontsize=6)
        plt.yticks(fontsize=6)
        plt.title(label=sample_name, fontsize=6)

        # show plot and wait for user to close the plot
        plt.show(block=True)

# fold change
fold_frame = filtered[(filtered[abundance_col_names].applymap(
    lambda value: value > 0 and not isnan(value)).all(1))]

# calculate mean of abundance for each of the sample groups
all_mean_cols = (('adult_mean', ('M1', 'M2')), ('piglet_mean', ('P1', 'P2')))
for mean_col, sample_names in all_mean_cols:
    col_names = get_column_names(item_attrs, LABEL_FREE_QUANT, sample_names)
    fold_frame[mean_col] = fold_frame[col_names].mean(
        axis=1)  # remember axis=1

RATIO = 'ratio'
# calculate fraction of the abundance mean columns
fold_frame[RATIO] = fold_frame[all_mean_cols[0][0]] / fold_frame[
    all_mean_cols[1][0]]

#fold_frame['index'] = range(len(fold_frame))
#fold_frame[fold_frame[RATIO]<2](plot.scatter(x='index',y='ratio')
#plt.show(block=True)

# keep rows with ratio > 2 or ratio < 0.5
above = fold_frame[RATIO] > 2
below = fold_frame[RATIO] < 0.5
fold_frame_filtered = fold_frame[above | below]

# only collagen
fold_frame_filtered = fold_frame_filtered[
    fold_frame_filtered[PG_PROTEINDESCRIPTIONS].str.contains('Collagen',
                                                             case=False)]
fold_frame_filtered.set_index(filtered.columns[2]).plot(y=RATIO,
                                                        kind='bar',
                                                        rot=0,
                                                        legend=False)

# use a smaller font size than the default
plt.xticks(fontsize=4)
plt.yticks(fontsize=4)
plt.title(label='fold', fontsize=6)
plt.show(block=True)

l, l1 = len(fold_frame), len(filtered)
