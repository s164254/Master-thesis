from math import isnan
import pandas as pd
import re
import plotutils


def is_unique_peptides_nan(value):
    return value == 0 or value == 1 or isnan(value)


ATTR_RE = '([_a-zA-Z0-9]+)\.raw\.PG\.(.+)'
UNIQUE_PEPTIDES = 'UniquePeptides'
LABEL_FREE_QUANT = 'Label-Free Quant'
PG_PROTEINDESCRIPTIONS = 'PG.ProteinDescriptions'
PG_GENES = 'PG.Genes'
RATIO = 'ratio'


class CsvToPandas:
    def __init__(self, csv_file) -> None:
        try:
            csv = pd.read_csv(csv_file)
        except:
            csv = pd.read_csv(csv_file, delimiter='\t')

        col_info = [(i, m.groups()[0], m.groups()[1], c)
                    for i, m, c in [(i, re.search(ATTR_RE, c), c)
                                    for i, c in enumerate(csv.columns)]
                    if m and len(m.groups()) == 2]

        self.sample_names = tuple(set([c[1] for c in col_info]))
        self.attr_names = tuple(set([c[2] for c in col_info]))

        self.col_info = col_info

        unique_peptides_col_names = self.get_column_names(UNIQUE_PEPTIDES)

        # remove rows where all unique peptides values are either 0, 1 or NAN
        filtered = csv[(csv[unique_peptides_col_names].applymap(
            lambda value: not is_unique_peptides_nan(value)).any(1))]

        # put newlines in protein description
        protein_desc_column = csv.columns.tolist().index(
            PG_PROTEINDESCRIPTIONS)
        filtered[filtered.columns[protein_desc_column]] = filtered.apply(
            lambda x: re.sub('\s+', '\n', x[protein_desc_column]), axis=1)

        # for the remaining rows set label-free quant to 0 if the corresponding unique peptides value is either 0, 1 or NAN
        abundance_col_names = self.get_column_names(LABEL_FREE_QUANT)
        for abundance_col_name, unique_peptides_col_name in zip(
                abundance_col_names, unique_peptides_col_names):
            filtered[abundance_col_name] = filtered.apply(
                lambda x: 0 if is_unique_peptides_nan(x[
                    unique_peptides_col_name]) else x[abundance_col_name],
                axis=1)

        self.abundance_col_names = abundance_col_names
        self.unique_peptides_col_names = unique_peptides_col_names
        self.filtered = filtered

    def get_column_names(self, attr_name, sample_names=[]):
        return [
            x[3] for x in self.col_info if x[2] == attr_name and (
                (len(sample_names) == 0) or (x[1] in sample_names))
        ]

    def fold_analysis(self, groups, protein_description_filters):
        fold_frame = self.filtered[(
            self.filtered[self.abundance_col_names].applymap(
                lambda value: value > 0 and not isnan(value)).all(1))]

        # calculate mean of abundance for each of the sample groups
        for mean_col, sample_names in groups:
            col_names = self.get_column_names(LABEL_FREE_QUANT, sample_names)
            fold_frame[mean_col] = fold_frame[col_names].mean(
                axis=1)  # remember axis=1

        # calculate fraction of the abundance mean columns
        fold_frame[RATIO] = fold_frame[groups[0][0]] / fold_frame[groups[1][0]]

        # keep rows with ratio > 2 or ratio < 0.5
        above = fold_frame[RATIO] > 2
        below = fold_frame[RATIO] < 0.5
        fold_frame_ratio_filtered = fold_frame[above | below]

        # make a plot for each protein_description_filter
        for protein_description_filter in protein_description_filters:
            fold_frame_filtered = fold_frame_ratio_filtered[
                protein_description_filter(
                    fold_frame_ratio_filtered[PG_PROTEINDESCRIPTIONS])]

            plotutils.dataframe_plot(
                fold_frame_filtered,
                lambda df: df.set_index(fold_frame_filtered.columns[2]).plot(
                    y=RATIO, kind='bar', rot=0, legend=False),
                '',
                block=True)

    def to_gene_list(self, gen_fname):
        gene_list = []
        for column_name in self.get_column_names(LABEL_FREE_QUANT):
            # sort by values in label-free quant column for sample_name
            df = self.filtered.sort_values(by=[column_name], ascending=False)
            df = df[df[column_name] > 0]
            genes = df[PG_GENES]
            genes = [x for x in genes.values if isinstance(x, str)]
            gene_list.append(genes)
            with open(gen_fname('%s.csv' % (column_name,)), 'w') as f:
                f.write('\n'.join(genes))
        
        with open(gen_fname('common.csv'), 'w') as f:
            f.write('\n'.join(sorted(set.intersection(*map(set,gene_list)))))

    def to_csv(self, csv_file) -> None:
        self.filtered.to_csv(csv_file)

    def group_analysis(self, sample_names, protein_description_filters, N):
        """
        A bar plot is made for each protein_description_filter of the N most frequently occuring proteins
        """
        for protein_description_filter in protein_description_filters:
            # select the rows that match the protein_description_filter
            filtered = self.filtered[protein_description_filter(
                self.filtered[PG_PROTEINDESCRIPTIONS])]

            column_names = self.get_column_names(LABEL_FREE_QUANT,
                                                 sample_names)
            all_column_values = [
                sorted(list(
                    zip(filtered[PG_PROTEINDESCRIPTIONS],
                        filtered[column_name])),
                       key=lambda x: x[1],
                       reverse=True) for column_name in column_names
            ]

            # continue N common proteins are found across all samples
            n = N
            n_max = min([len(x) for x in all_column_values])
            num_samples = len(all_column_values)
            common = set.intersection(
                *map(set, [[x[0] for x in column_values[:n]]
                           for column_values in all_column_values]))
            while n < n_max and len(common) < N:
                n += 1
                common = set.intersection(
                    *map(set, [[x[0] for x in column_values[:n]]
                               for column_values in all_column_values]))

            # create a new dataframe based on the list of sample values for each of the N found proteins
            d = {PG_PROTEINDESCRIPTIONS: list(common)}
            for sample_name, column_values in zip(sample_names,
                                                  all_column_values):
                d.update({
                    sample_name:
                    [[x[1] for x in column_values if x[0] == protein][0]
                     for protein in common]
                })

            nmost = pd.DataFrame(data=d)

            plotutils.dataframe_plot(
                nmost,
                lambda df: df.plot(x=PG_PROTEINDESCRIPTIONS,
                                   y=list(sample_names),
                                   kind='bar',
                                   rot=0,
                                   legend=False),
                'N most common proteins',
                block=True)
