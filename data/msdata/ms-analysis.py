from statistics import mean
from os import path
import re
import fileutils as ft
import utils
import numpy as np
import matplotlib.pyplot as plt

ATTR_RE = '_([a-zA-Z0-9]+)\.raw\.PG\.(.+)'


def common_proteins_criteria(gene_values):
    return all(v > 0 for v in gene_values)

def value_list(dict,samples,attr_name,i):
    return [dict[sample][attr_name][i] for sample in samples]


def is_fold(v1, v2, n):
    return max(v1,v2) / min(v1,v2) >= n

def is_two_fold(v1, v2):
    return is_fold(v1, v2, 2)


UNIQUE_PEPTIDES = 'UniquePeptides'
LABEL_FREE_QUANT = 'Label-Free Quant'


def plot_hist(m, sample, plot_title, x_label, y_label):
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


class MsAnalysis:

    def __init__(self, csv_file_name):
        self.header, rows = ft.load_csv(csv_file_name, ',', True)

        item_attrs = [(i, m.groups()[0], m.groups()[1])
                      for i, m in [(i, re.search(ATTR_RE, c))
                                   for i, c in enumerate(self.header)]
                      if m and len(m.groups()) == 2]
        self.attrs = tuple(
            set([attr_name for i, sample, attr_name in item_attrs]))
        self.attr_first_column_index = dict([
            (attr,
             min([
                 i for i, sample, attr_name in item_attrs
                 if attr_name == attr
             ])) for attr in self.attrs
        ])
        self.samples = tuple(
            set([sample for i, sample, attr_name in item_attrs]))
        self.item_attrs = item_attrs

        unique_peptides_first_column_index = self.attr_first_column_index[
            UNIQUE_PEPTIDES]
        sample_len = len(self.samples)
        # remove 'empty' rows
        rows = [
            row for row in rows if utils.unique_peptides_valid_row(
                row[unique_peptides_first_column_index:
                    unique_peptides_first_column_index + sample_len])
        ]

        sample_data = {}
        for col_index, sample, attr_name in self.item_attrs:
            if attr_name == UNIQUE_PEPTIDES:
                continue

            if not sample in sample_data:
                sample_data[sample] = {}

            NOT_VALID_VALUE = 0
            unique_peptides_offset = self.get_offset(col_index, attr_name,
                                                     sample,
                                                     UNIQUE_PEPTIDES)

            for row_index, row in enumerate(rows):
                if not utils.unique_peptides_valid_value(
                        row[col_index + unique_peptides_offset]):
                    row[col_index] = '0'

            sample_data[sample][attr_name] = tuple([
                utils.parse_csv_cell(row[col_index], NOT_VALID_VALUE)
                for row in rows
            ])

        self.UniprotIds = tuple([row[0] for row in rows])
        self.Genes = tuple([row[1] for row in rows])
        self.ProteinDescriptions = tuple([row[2] for row in rows])

        self.rows = rows
        self.sample_data = sample_data

    def get_offset(self, from_col_index, from_attr_name, from_sample,
                   to_attr_name_):
        return [
            to_col_index - from_col_index
            for to_col_index, to_sample, to_attr_name in self.item_attrs
            if to_attr_name == to_attr_name_ and to_sample == from_sample
        ][0]

    def unique_peptides_valid_row_criteria(self, row):
        return utils.unique_peptides_valid_row(
            row[self.unique_peptides_first_column_index:self.
                unique_peptides_first_column_index + len(self.samples)])
        #column_values = row[]

    def filter_helper(self, samples, attr_name, criteria):
        sd = self.sample_data
        
        flattened = [value_list(sd,samples,attr_name,i) for i in range(len(self.Genes))]
        return [(self.UniprotIds[i], self.Genes[i], self.ProteinDescriptions[i],values) for i,values in enumerate(flattened) if criteria(values)]

    def common_proteins(self, samples, attr_name):
        return self.filter_helper(samples, attr_name,
                                  common_proteins_criteria)

    def grouped_common_proteins(self, samples, attr_name, groups):
        rows = self.common_proteins(self.samples,LABEL_FREE_QUANT)
        avgs = [[mean((row[-1][i1],row[-1][i2])) for i1,i2 in groups] for row in rows]
        indices = [i for i,avg in enumerate(avgs) if is_fold(avg[0],avg[1],2)]
        return self.filter_helper(samples, attr_name,common_proteins_criteria)

    def fold_change(self, samples, attr_name):
        common = self.common_proteins(samples,attr_name)
        return self.filter_helper(samples, attr_name, is_two_fold)

    def save_csv(self, csv_file_name, attr_name):
        """"""
        first_col_index = self.attr_first_column_index[attr_name]
        sample_len = len(self.samples)

        ft.write_csv(
            csv_file_name, self.header[:3] +
            self.header[first_col_index:first_col_index + sample_len], [
                row[:3] + row[first_col_index:first_col_index + sample_len]
                for row in self.rows
            ])

    def sort_indices(self, sample, attr_name, reverse):
        indices = [
            i[0]
            for i in sorted(enumerate(self.sample_data[sample][attr_name]),
                            key=lambda x: x[1])
        ]
        return reversed and indices.reverse() or indices

    def get_values(self, sample, attr_name, indices):
        data = self.sample_data[sample][attr_name]
        return [data[i] for i in indices]

    def get_values_by_genes(self, sample, attr_name, genes):
        return self.get_values(sample,attr_name,[self.get_index(g) for g in genes])

    def get_prot(self, indices):
        return [
            re.sub('\s+', '\n', self.ProteinDescriptions[i])[:30]
            for i in indices
        ]
    
    def get_index(self,gene):
        return [i for i,x in enumerate(self.Genes) if x==gene][0]

script_dir = ft.get_script_dir(__file__)
fname = 'Batch2_data.csv'

m = MsAnalysis(path.join(script_dir, fname))
common_proteins = m.grouped_common_proteins(m.samples,LABEL_FREE_QUANT,((0,1),(2,3)))

l = m.sample_data['M1'][LABEL_FREE_QUANT]
l = m.sample_data['M2'][LABEL_FREE_QUANT]

m.save_csv(path.join(script_dir, 'csv.csv'), LABEL_FREE_QUANT)
#common_proteins = m.common_proteins(samples, LABEL_FREE_QUANT)
#fold_change = m.fold_change(samples, LABEL_FREE_QUANT)

sample = 'M1'
genes = ('COL6A1','COL1A2')
print(m.get_values_by_genes(sample,LABEL_FREE_QUANT,genes))
plot_hist(m, 'M1', '10 proteins with largest abundance values for %s',
        'Protein description', 'Normalized abundance')
for sample in m.samples:
    plot_hist(m, sample, '10 proteins with largest abundance values for %s',
          'Protein description', 'Normalized abundance')
i = 1
