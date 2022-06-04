from locale import normalize
from math import isnan
from operator import index
import pandas as pd
import re
import plotutils
import experiment_args
import fileutils as ft


def is_unique_peptides_nan(value):
    return value == 0 or value == 1 or isnan(value)


ATTR_RE = '([_a-zA-Z0-9]+)\.raw\.PG\.(.+)'
UNIQUE_PEPTIDES = 'UniquePeptides'
LABEL_FREE_QUANT = 'Label-Free Quant'
PG_PROTEINDESCRIPTIONS = 'PG.ProteinDescriptions'
PG_PROTEINDESCRIPTIONS_NEWLINE = 'ProteinDescriptions'
PG_GENES = 'PG.Genes'
RATIO = 'ratio'


def nmost_common(lists, N, common_column_idx, df_column_names, normalize=False):
    n = N
    n_max = min([len(l) for l in lists])
    common = set.intersection(
        *map(set, [[row[common_column_idx] for row in lst[:n]]
                   for lst in lists]))
    while n < n_max and len(common) < N:
        n += 1
        common = set.intersection(
            *map(set, [[row[common_column_idx] for row in lst[:n]]
                       for lst in lists]))

    #return [[row for row in lst[:n] if row[common_column_idx] in common]
    #        for lst in lists]
    # sorteret

    # create dict for resuting dataframe
    rows = []
    for gene in sorted(common):
        row_values = []
        for lst in lists:
            row_values.append([row[(common_column_idx+1) % 2] for row in lst[:n] if row[common_column_idx] == gene][0])
        if normalize:
            mx = max(row_values)
            row_values = [v / mx for v in row_values] # brug numpy
        rows.append(row_values)

    common = list(sorted(common))
    d = {df_column_names[0]: common}
    idx = 0
    for lst, column_name in zip(lists, df_column_names[1:]):
        d[column_name] = [row[idx] for row in rows]
        idx += 1
    return pd.DataFrame(d)


class CsvToPandas:
    def __init__(self, args) -> None:
        args = isinstance(
            args, str) and experiment_args.to_experiment_args(args) or args
        self.args = args
        try:
            csv = pd.read_csv(args.filename)
        except:
            csv = pd.read_csv(args.filename, delimiter='\t')

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
        filtered[PG_PROTEINDESCRIPTIONS_NEWLINE] = filtered.apply(
            lambda x: re.sub('\s+', '\n', x[protein_desc_column]), axis=1)

        # for the remaining rows set label-free quant to 0 if the corresponding unique peptides value is either 0, 1 or NAN
        abundance_col_names = self.get_column_names(LABEL_FREE_QUANT)
        for abundance_col_name, unique_peptides_col_name in zip(
                abundance_col_names, unique_peptides_col_names):
            filtered[abundance_col_name] = filtered.apply(
                lambda x: 0 if is_unique_peptides_nan(x[
                    unique_peptides_col_name]) else x[abundance_col_name],
                axis=1)

        # set abundance value to 0 if it is NAN
        for abundance_col_name in abundance_col_names:
            filtered[abundance_col_name] = filtered.apply(lambda x: isnan(x[
                abundance_col_name]) and 0 or x[abundance_col_name],
                                                          axis=1)

        self.abundance_col_names = abundance_col_names
        self.unique_peptides_col_names = unique_peptides_col_names
        self.filtered = filtered

    def gene_analysis(self, N, title, xlabel, ylabel, sample_names=None):
        gene_abundance_list = []
        if not sample_names:
            sample_names = self.sample_names
        sample_names = sorted(sample_names)

        sample_names_key = '_'.join([sn.lower() for sn in sample_names])
        common_proteins = self.args.common_proteins[sample_names_key]

        column_names = self.get_column_names(LABEL_FREE_QUANT, sample_names)
        for column_name in column_names:
            # create dataframe where all rows have a value > 0 in all abundance columns
            df = self.filtered
            df = df[df[column_name] > 0].copy()

            # sort with largest abundance value first
            df = df.sort_values(by=[column_name], ascending=False)

            # extract gene and abundance value to a list of tuples
            l = list(zip(df[PG_GENES], df[column_name]))
            gene_abundance_list.append(
                [row for row in l if row[0] in common_proteins])

        res = nmost_common(gene_abundance_list, N, 0,
                           [PG_GENES] + column_names, normalize=True)
        fig_filename = self.args.fig_filename(
                'batch_to_batch-common-cellular-analysis.%s.png' %
                (sample_names_key, ))
        plotutils.dataframe_plot(
            res,
            lambda df: df.plot(
                x=PG_GENES, y=column_names, kind='bar', rot=0, legend=True, ylim=(0,1.2)),
            title,
            axis_setup_func=None,  #lambda ax: ax.get_xaxis().set_ticklabels([]),
            plot_setup_func=None,
            xlabel=xlabel,
            ylabel=ylabel,
            block=True,
            fig_filename=None #fig_filename
        ) 

    def get_column_names(self, attr_name, sample_names=None):
        return [
            x[3] for x in self.col_info if x[2] == attr_name
            and self.sample_name_match(x[1], sample_names)
        ]

    def sample_name_match(self, sample_name, sample_names):
        if not sample_names:
            return True

        matches = [
            sn for sn in sample_names if sample_name.find('_%s' % (sn, )) > 0
        ]
        if len(matches) > 1:
            raise Exception(
                'sample_name_match: More than one matched sample: sample_name:%s, sample_names:%s'
                % (sample_name, sample_names))

        return matches and matches[0] or None

    def columnname_to_samplename(self, column_name):
        csv_samplename = [c for c in self.col_info if c[3] == column_name][0]
        pass

    def fold_analysis(self, sample_names, group_name, use_filter,
                      protein_description_filter):
        # create dataframe where all rows have a value > 0 in all abundance columns
        column_names = self.get_column_names(LABEL_FREE_QUANT, sample_names)

        df = self.filtered[(self.filtered[column_names].applymap(
            lambda value: value > 0).all(1))].copy()

        # calculate fraction of the abundance mean columns
        df[RATIO] = df[column_names[0]] / df[column_names[1]]

        # keep rows with ratio > 2 or ratio < 0.5
        if use_filter:
            above = df[RATIO] > 2
            below = df[RATIO] < 0.5
            df = df[above | below]
        else:
            df = df[df[RATIO] > 0]

        # make a plot for each protein_description_filter
        fig_samplenames = '_'.join([sn.lower() for sn in sample_names])
        try:
            plot_df = df[protein_description_filter(
                df[PG_PROTEINDESCRIPTIONS])].copy()
        except Exception as ex:
            df.to_csv(
                self.args.fig_filename('ecm_foldchange.%s.%s.err.csv' %
                                       (fig_samplenames, group_name)))
            return

        if len(plot_df) == 0:
            return

        plotutils.dataframe_plot(
            plot_df,
            lambda x: x.set_index(plot_df[PG_PROTEINDESCRIPTIONS_NEWLINE]).
            plot(y=RATIO, kind='bar', rot=0, legend=False),
            '',
            axis_setup_func=lambda ax: ax.bar_label(ax.containers[0]),
            fig_filename=self.args.fig_filename('ecm_foldchange.%s.%s.png' %
                                                (fig_samplenames, group_name)),
            block=True)

    def to_gene_list(self):
        for column_name in self.get_column_names(LABEL_FREE_QUANT):
            # sort by values in label-free quant column for sample_name
            df = self.filtered.sort_values(by=[column_name],
                                           ascending=False).copy()

            # copy rows with abundance value > 0 to a new dataframe
            df = df[df[column_name] > 0].copy()

            # get list of corresponding gene id's but remove rows with invalid gene id's (NAN)
            genes = [(g, d) for g, d in zip(df[PG_GENES].values,
                                            df[PG_PROTEINDESCRIPTIONS].values)
                     if isinstance(g, str)]

            # write list of genes to CSV file
            # todo: write to file in experiment output dir.
            ft.to_file(
                self.args.gene_filename('%s.genelistdesc.csv' %
                                        (column_name, )),
                '\n'.join(['%s,%s' % x for x in genes]))
            ft.to_file(
                self.args.gene_filename('%s.genelist.csv' % (column_name, )),
                '\n'.join([x[0] for x in genes]))

    def to_csv(self, csv_file) -> None:
        self.filtered.to_csv(csv_file)

    def group_analysis(self, sample_names, protein_description_filters, N):
        """
        A bar plot is made for each protein_description_filter of the N most frequently occuring proteins
        """
        for protein_description_filter in protein_description_filters:
            # select the rows that match the protein_description_filter
            filtered = self.filtered[protein_description_filter(
                self.filtered[PG_PROTEINDESCRIPTIONS_NEWLINE])]

            column_names = self.get_column_names(LABEL_FREE_QUANT,
                                                 sample_names)
            all_column_values = [
                sorted(list(
                    zip(filtered[PG_PROTEINDESCRIPTIONS_NEWLINE],
                        filtered[column_name])),
                       key=lambda x: x[1],
                       reverse=True) for column_name in column_names
            ]

            # continue N common proteins are found across all samples
            n = N
            n_max = min([len(x) for x in all_column_values])
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
