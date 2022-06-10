from os import path
from math import isnan
import pandas as pd
import re
import plotutils
import experiment_args
import fileutils as ft
import numpy as np
import warnings


def is_unique_peptides_nan(value):
    return value == 0 or value == 1 or isnan(value)


ATTR_RE = '([_a-zA-Z0-9]+)\.raw\.PG\.(.+)'
UNIQUE_PEPTIDES = 'UniquePeptides'
LABEL_FREE_QUANT = 'Label-Free Quant'
PG_PROTEINDESCRIPTIONS = 'PG.ProteinDescriptions'
PG_PROTEINDESCRIPTIONS_NEWLINE = 'ProteinDescriptions'
UNIPROTID_PROTEINDESCRIPTIONS_NEWLINE = 'UniprotidProteinDescriptions'
PG_GENES = 'PG.Genes'
PG_UNIPROTID = 'PG.UniprotIds'
RATIO = 'ratio'

CELLULAR_SEARCH_COLUMNS = [
    'PG.CellularComponent', 'PG.BiologicalProcess', 'PG.MolecularFunction',
    PG_PROTEINDESCRIPTIONS
]


def legend_outside_chart(ax):
    pos = ax.get_position()
    ax.set_position([pos.x0, pos.y0, pos.width * 1.0, pos.height])
    ax.legend(loc='center right', bbox_to_anchor=(1.25, 0.5))


def write_gene_2_genedesc(genes, gene_2_genedesc, fname):
    ft.to_file(
        fname,
        '\n'.join('%s=%s' % (gene, gene_2_genedesc[gene]) for gene in genes))


def parse_uniprotid(uniprotid):
    return uniprotid.split(';')[0]


def dataframe_applymap_on_rows(df, column_names, column_func, all_column_func):
    return df[all_column_func(df[column_names].applymap(column_func))]

    # csv[(csv[unique_peptides_col_names].applymap(
    #         lambda value: not is_unique_peptides_nan(value)).any(1))]


def possible_nan_2_str(v):
    return isinstance(v, float) and 'NAN' or v


def set_bar_labels(ax, fmt='%.2f'):
    for c in ax.containers:
        ax.bar_label(c, fmt=fmt)


def nmost_common(lists, N, common_column_idx, df_column_names):
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
            row_values.append([
                row[(common_column_idx + 1) % 2] for row in lst[:n]
                if row[common_column_idx] == gene
            ][0])
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
        warnings.filterwarnings("ignore")
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
        test = dataframe_applymap_on_rows(
            csv, unique_peptides_col_names,
            lambda value: not is_unique_peptides_nan(value),
            lambda x: x.any(1)).copy()

        filtered = csv[(csv[unique_peptides_col_names].applymap(
            lambda value: not is_unique_peptides_nan(value)).any(1))].copy()

        # put newlines in protein description
        protein_desc_column = csv.columns.tolist().index(
            PG_PROTEINDESCRIPTIONS)
        filtered[PG_PROTEINDESCRIPTIONS_NEWLINE] = filtered.apply(
            lambda x: re.sub('\s+', '\n', x[protein_desc_column]), axis=1)
        filtered[UNIPROTID_PROTEINDESCRIPTIONS_NEWLINE] = filtered.apply(
            lambda x: '%s\n%s' %
            (x[PG_GENES], x[PG_PROTEINDESCRIPTIONS_NEWLINE]),
            axis=1)

        # for the remaining rows set label-free quant to 0 if the corresponding unique peptides value is either 0, 1 or NAN
        abundance_col_names = self.get_column_names(LABEL_FREE_QUANT)
        for abundance_col_name, unique_peptides_col_name in zip(
                abundance_col_names, unique_peptides_col_names):
            filtered[abundance_col_name] = filtered.apply(
                lambda x: 0
                if (isnan(x[abundance_col_name]) or is_unique_peptides_nan(x[
                    unique_peptides_col_name])) else x[abundance_col_name],
                axis=1)

        self.abundance_col_names = abundance_col_names
        self.unique_peptides_col_names = unique_peptides_col_names
        self.filtered = filtered

        self.gene_2_genedesc = dict(
            self.filtered[[PG_GENES, PG_PROTEINDESCRIPTIONS]].values)

        self.experiment_lookup = {}

    def cellular_analysis(self, N, title, xlabel, ylabel, sample_names=None):
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

        fname_base = self.args.fig_filename(
            'batch_to_batch-common-cellular-analysis.%s' %
            (sample_names_key, ))
        display_column_sample_names = [
            self.get_output_name(cn, fname_base) for cn in column_names
        ]
        res = nmost_common(gene_abundance_list, N, 0,
                           [PG_GENES] + display_column_sample_names)
        
        write_gene_2_genedesc(res[PG_GENES].values, self.gene_2_genedesc,
                              fname_base + '.txt')
        
        fig_filename = fname_base + '.png'
        fig_filename = ''
        mx = max(res[display_column_sample_names].max())
        plotutils.dataframe_plot(
            res,
            lambda df: df.plot(x=PG_GENES,
                               y=display_column_sample_names,
                               kind='bar',
                               rot=90,
                               legend=True,
                               logy=True),
            title,
            axis_setup_func=None,  #lambda ax: ax.get_xaxis().set_ticklabels([]),
            plot_setup_func=None,
            xlabel=xlabel,
            ylabel=ylabel,
            block=True,
            fig_filename=fig_filename,
            ylim=(0, mx * 1.15))

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

    def get_output_name(self, column_name, lookup_basename=None):
        csv_samplename = [c for c in self.col_info if c[3] == column_name]
        if len(csv_samplename) != 1:
            raise Exception(
                'get_output_name: expected 1 sample name matching %s, got %d' %
                (column_name, len(csv_samplename)))

        lookup_dict = self.args.sample_name_lookup
        if lookup_basename:
            if not lookup_basename in self.experiment_lookup:
                lookup_filename = self.args.lookup_filename(lookup_basename +
                                                            '.txt')
                if path.exists(lookup_filename):
                    content =ft.to_dict(
                        lookup_filename, '=')
                    if content:
                        self.experiment_lookup[lookup_basename] = content
                    else:
                        print('empty lookup file:' + lookup_filename)
 
                else:
                    ft.to_file(lookup_filename, '')
                    print('MISSING lookup file:' + lookup_filename)
            if lookup_basename in self.experiment_lookup:
                lookup_dict = self.experiment_lookup[lookup_basename]

        matches = [
            v for k, v in lookup_dict.items()
            if csv_samplename[0][3].find('%s.' % (k, )) > 0
        ]
        if len(matches) != 1:
            raise Exception(
                'get_output_name: expected 1 matching lookup sample name key, got %d'
                % (csv_samplename[0], len(matches)))

        return matches[0]

    def fold_analysis(self, sample_names, group_name, use_ratio_filter,
                      protein_description_filter, remove_non_existing):
        column_names = self.get_column_names(LABEL_FREE_QUANT, sample_names)

        # apply protein_description_filter and copy to new dataframe
        df = self.filtered[protein_description_filter(
            self.filtered[PG_PROTEINDESCRIPTIONS])].copy()

        for column_name in column_names:
            if True or remove_non_existing:  # remove True if transformation below is fixed
                df = df[df[column_name] > 0].copy()
            else:
                pass
                # need to fix this df[column_name] = df[df[column_name] <= 0] = 0.00001

        # calculate fraction of the abundance mean columns
        df[RATIO] = df[column_names[0]] / df[column_names[1]]

        # keep rows with ratio > 2 or ratio < 0.5
        if use_ratio_filter:
            above = df[RATIO] > 2
            below = df[RATIO] < 0.5
            df = df[above | below]
        else:
            df = df[df[RATIO] > 0]

        if len(df) == 0:
            return

        # try:
        #     plot_df = df[protein_description_filter(
        #         df[PG_PROTEINDESCRIPTIONS])].copy()
        # except Exception as ex:
        #     df.to_csv(
        #         self.args.fig_filename('ecm_foldchange.%s.%s.err.csv' %
        #                                (fig_samplenames, group_name)))
        #     return

        fig_samplenames = '_'.join([sn.lower() for sn in sample_names])
        fname_base = self.args.fig_filename('ecm_foldchange.%s.%s' %
                                            (fig_samplenames, group_name))

        write_gene_2_genedesc(df[PG_GENES].values, self.gene_2_genedesc,
                              fname_base + '.txt')
        fig_filename = fname_base + '.png'
        #fig_filename = ''

        mx = max(df[column_names].max())
        self.to_output_and_plot(
            df = df,
            columns = [PG_GENES],
            sample_names = column_names,
            fname_base=fname_base,
            plot_func=lambda inp: plotutils.dataframe_plot(
                inp[0],
                lambda x: x.set_index(inp[0][PG_GENES]).plot(
                    y=inp[1],
                    kind='bar',
                    rot=90,
                    legend=True,
                    logy=True
                    #ylim=normalize and (0, 1.2) or None
                ),
                '',
                axis_setup_func=None,
                fig_filename=fig_filename,
                block=True,
                ylim=(0, mx * 1.2)))

    def ecm_common(self, sample_names, group_name, use_filter,
                   protein_description_filter):
        column_names = self.get_column_names(LABEL_FREE_QUANT, sample_names)
        for column_name in column_names:
            other_column_names = [
                cn for cn in column_names if cn != column_name
            ]
            column_display_names = map(self.get_output_name,
                                       [column_name] + other_column_names)

            # start by removing rows where sample does not have a LABEL_FREE_QUANT value
            df = self.filtered[self.filtered[column_name] > 0].copy()

            df = df[protein_description_filter(
                df[PG_PROTEINDESCRIPTIONS])].copy()

            uniprotid = df[PG_UNIPROTID].values
            genenames = df[PG_GENES].values
            protein_desc = df[PG_PROTEINDESCRIPTIONS].values

            other_values = df[other_column_names].values
            other_missing = [(i, ' '.join([
                self.get_output_name(other_column_name)
                for value, other_column_name in zip(row, other_column_names)
                if value <= 0
            ])) for i, row in enumerate(other_values)]
            other_missing = [
                ','.join((uniprotid[x[0]], possible_nan_2_str(genenames[x[0]]),
                          possible_nan_2_str(protein_desc[x[0]]), x[1]))
                for x in other_missing if x[1]
            ]
            ft.to_file(
                self.args.common_filename(
                    '%s.%s.csv' %
                    (group_name, '_'.join(column_display_names))),
                '\n'.join(other_missing))

        # one file with all sample names
        column_names = self.get_column_names(LABEL_FREE_QUANT)
        df = self.filtered[protein_description_filter(
            self.filtered[PG_PROTEINDESCRIPTIONS])].copy()
        column_display_names = map(self.get_output_name, column_names)
        for column_name, column_display_name in zip(column_names,
                                                    column_display_names):
            df[column_name] = df[column_name].apply(
                lambda value: value > 0 and column_display_name or '')
        output_columns = [PG_UNIPROTID, PG_GENES, PG_PROTEINDESCRIPTIONS
                          ] + column_names
        output = df[output_columns].values
        ft.to_file(
            self.args.common_filename('%s.csv' % (group_name, )), '\n'.join([
                ','.join([possible_nan_2_str(v) for v in row])
                for row in output
            ]))

    def to_gene_list(self):
        for column_name in self.get_column_names(LABEL_FREE_QUANT):
            # sort by values in label-free quant column for sample_name
            df = self.filtered.sort_values(by=[column_name],
                                           ascending=False).copy()

            # copy rows with abundance value > 0 to a new dataframe
            df = df[df[column_name] > 0].copy()

            # get list of corresponding gene id's but remove rows with invalid gene id's (NAN)
            genes = [
                (parse_uniprotid(u), possible_nan_2_str(g), d)
                for u, g, d in zip(df[PG_UNIPROTID].values, df[PG_GENES].
                                   values, df[PG_PROTEINDESCRIPTIONS].values)
                if True or isinstance(g, str)
            ]  # add 'or True' to get NAN as well

            # write list of genes to CSV file
            # todo: write to file in experiment output dir.
            ft.to_file(
                self.args.gene_filename('%s.genelistdesc.csv' %
                                        (column_name, )),
                '\n'.join(['%s,%s,%s' % x for x in genes]))
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
                sorted(list(zip(filtered[PG_GENES], filtered[column_name])),
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
            d = {PG_GENES: list(common)}
            for sample_name, column_values in zip(sample_names,
                                                  all_column_values):
                d.update({
                    sample_name:
                    [[x[1] for x in column_values if x[0] == protein][0]
                     for protein in common]
                })

            nmost = pd.DataFrame(data=d)

            plotutils.dataframe_plot(nmost,
                                     lambda df: df.plot(x=PG_GENES,
                                                        y=list(sample_names),
                                                        kind='barh',
                                                        rot=90,
                                                        legend=True),
                                     'N most common proteins',
                                     block=True)

    def get_cellular_dataframe(self,gene_column_name):
        if not any([
                c for c in self.filtered.columns
                if c.upper() == 'PG.CellularComponent'.upper()
        ]):
            raise Exception(
                'get_cellular_dataframe: missing search columns in dataframe '
                + self.args.experiment_name)

        search_for_text = [
            s.lower()
            for s in ('extracellular', 'basement', 'collagen', 'fibronectin',
                      'laminin', 'Elastin', 'Proteoglycan',
                      'Uncharacterized protein', 'fibrinogen', 'growth factor',
                      'vimentin', 'Inter-alpha-trypsin inhibitor heavy chain',
                      'hemoglobin', 'myoglobin', 'transglutaminase',
                      'Cartilage intermediate layer protein',
                      'glycogen phosphorylase','nidogen')
        ]
        has_match = lambda txt: any(
            (1 for search_for in search_for_text
             if isinstance(txt, str) and txt.lower().find(search_for) >= 0))

        self.filtered[CELLULAR_SEARCH_COLUMNS].to_csv(
            'CELLULAR_SEARCH_COLUMNS')

        df = self.filtered[(self.filtered[CELLULAR_SEARCH_COLUMNS].applymap(
            lambda value: not has_match(value)).all(1))].copy()

        # remove rows with invalid uniprotid
        df = df[df[gene_column_name].apply(lambda x: isinstance(x, str))]

        return df.copy()

    def generate_cellular_files(self):
        column_names = self.get_column_names(LABEL_FREE_QUANT)
        display_column_sample_names = [
            self.get_output_name(cn) for cn in column_names
        ]
        f=self.filtered.copy()

        for column_name, display_column_name in zip( column_names, display_column_sample_names):
                    # get dataframe with cellular rows only
            df = self.get_cellular_dataframe(PG_UNIPROTID)

            df = df[df[column_name]>0]
            # df = dataframe_applymap_on_rows(df, [column_name],
            #                                 lambda value: value > 0,
            #                                 lambda x: x.any(1))

            df = df[df[PG_UNIPROTID].apply(lambda x: isinstance(x, str))]

            cellular_uniprotids = df[PG_UNIPROTID].values
            total = len(f[f[column_name]>0])
            print('Cellular:%s, Extracellular:%d' %
                (len(cellular_uniprotids),
                total - len(cellular_uniprotids)))
            ft.to_file(
                ft.msdata_filename('%s.%s.cellular.txt' %
                                (self.args.experiment_name, display_column_name )),
                '|'.join(cellular_uniprotids))

    def cellular_analysis_2(self, sample_names=None):
        '''cellular_analysis on csv file with additional cellular description columns'''
        # get dataframe with cellular rows only
        df = self.get_cellular_dataframe()

        all_rows = self.filtered[self.filtered[PG_GENES].apply(
            lambda x: isinstance(x, str))].copy()

        for column_sample_name in self.get_column_names(
                LABEL_FREE_QUANT, sample_names):

            # row count for sample with value > 0
            total = len(all_rows[all_rows[column_sample_name] > 0])

            # remove rows with sample value <= 0
            df_sample = df[df[column_sample_name] > 0]

            # print total for sample
            print('%s: sample cellular rows:%d, all cellular rows:%d' %
                  (self.get_output_name(column_sample_name), len(df_sample),
                   total))

            # order by abundance value desc
            #df_sample = df_sample.sort_values(by=[column_sample_name], ascending=False)

    def cellular_analysis_3(self, sample_names, title_func, xlabel, ylabel):
        '''cellular_analysis on csv file with additional cellular description columns'''
        # get dataframe with cellular rows only
        df = self.get_cellular_dataframe(PG_GENES)

        uniprotid_abundance_lists = []
        column_sample_names = self.get_column_names(LABEL_FREE_QUANT,
                                                    sample_names)
        sample_names_key = '_'.join(sample_names)
        for column_sample_name in column_sample_names:
            # remove rows with sample value <= 0
            df_sample = df[df[column_sample_name] > 0]

            # order by abundance value desc
            df_sample = df_sample.sort_values(by=[column_sample_name],
                                              ascending=False)

            uniprotid_abundance_lists.append(
                df_sample[[PG_GENES, column_sample_name]].values)

        fname_base = self.args.fig_filename(
            'batch_to_batch-common-cellular-analysis.%s' %
            (sample_names_key, ))

        display_column_sample_names = [
            self.get_output_name(cn, fname_base) for cn in column_sample_names
        ]
        res = nmost_common(uniprotid_abundance_lists, 20, 0,
                           [PG_GENES] + display_column_sample_names)

        write_gene_2_genedesc(res[PG_GENES].values, self.gene_2_genedesc,
                              fname_base + '.txt')

        fig_filename = fname_base + '.png'
        #fig_filename = None

        mx = max(res[display_column_sample_names].max())
        plotutils.dataframe_plot(
            res,
            lambda df: df.plot(x=PG_GENES,
                               y=display_column_sample_names,
                               kind='bar',
                               rot=90,
                               legend=True,
                               logy=True),
            title_func(self.args.experiment_name, display_column_sample_names),
            axis_setup_func=None,  #legend_outside_chart,
            plot_setup_func=None,
            xlabel=xlabel,
            ylabel=ylabel,
            block=True,
            fig_filename=fig_filename,
            ylim=(0, mx * 1.15))

    def to_output_dataframe(self, df, columns, sample_names, fname_base):
        res = pd.DataFrame(df[columns])
        output_sample_names = []
        for sample_name in sample_names:
            output_sample_name = self.get_output_name(sample_name, fname_base)
            output_sample_names.append(output_sample_name)
            res[output_sample_name] = df[sample_name]
        return res, output_sample_names

    def to_output_and_plot(self, df, columns, sample_names, plot_func,
                           fname_base):
        ret = self.to_output_dataframe(df, columns, sample_names, fname_base)
        plot_func(ret)
