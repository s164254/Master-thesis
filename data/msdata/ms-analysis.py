from os import path
import re
import fileutils as ft
import utils

ATTR_RE = '_([a-zA-Z0-9]+)\.raw\.PG\.(.+)'


def common_proteins_criteria(v1, v2):
    return v1 > 0 and v2 > 0


def is_fold(v1, v2, n):
    ma = max(v1, v2)
    if ma == 0:
        return False

    mi = min(v1, v2)
    if mi == 0:
        return False  #True

    return ma / mi >= 2


def is_two_fold(v1, v2):
    return is_fold(v1, v2, 2)


UNIQUE_PEPTIDES = 'UniquePeptides'
LABEL_FREE_QUANT = 'Label-Free Quant'


class MsAnalysis:

    def __init__(self, csv_file_name):
        self.header, rows = ft.load_csv(csv_file_name, ',', True)

        item_attrs = [(i, m.groups()[0], m.groups()[1])
                      for i, m in [(i, re.search(ATTR_RE, c))
                                   for i, c in enumerate(self.header)]
                      if m and len(m.groups()) == 2]
        self.attrs = tuple(
            set([attr_name for i, item_name, attr_name in item_attrs]))
        self.attr_first_column_index = dict([
            (attr,
             min([
                 i for i, item_name, attr_name in item_attrs
                 if attr_name == attr
             ])) for attr in self.attrs
        ])
        unique_peptides_first_column_index = self.attr_first_column_index[
            UNIQUE_PEPTIDES]
        self.item_names = tuple(
            set([item_name for i, item_name, attr_name in item_attrs]))
        self.item_attrs = item_attrs

        item_name_len = len(self.item_names)
        l1 = len(rows)
        # remove 'empty' rows
        rows = [
            row for row in rows if utils.unique_peptides_valid_row(
                row[unique_peptides_first_column_index:
                    unique_peptides_first_column_index + item_name_len])
        ]
        l2 = len(rows)

        item_data = {}
        for col_index, item_name, attr_name in self.item_attrs:
            if attr_name == UNIQUE_PEPTIDES:
                continue

            if not item_name in item_data:
                item_data[item_name] = {}

            NOT_VALID_VALUE = 0
            unique_peptides_offset = self.get_offset(col_index, attr_name,
                                                     item_name,
                                                     UNIQUE_PEPTIDES)

            for row_index, row in enumerate(rows):
                if not utils.unique_peptides_valid_value(
                        row[col_index + unique_peptides_offset]):
                    row[col_index] = '0'

            item_data[item_name][attr_name] = tuple([
                utils.parse_csv_cell(row[col_index], NOT_VALID_VALUE)
                for row in rows
            ])

        self.UniprotIds = tuple([row[0] for row in rows])
        self.Genes = tuple([row[1] for row in rows])
        self.ProteinDescriptions = tuple([row[1] for row in rows])

        item_data = {}
        for col_index, item_name, attr_name in self.item_attrs:
            if attr_name == UNIQUE_PEPTIDES:
                continue

            if not item_name in item_data:
                item_data[item_name] = {}

            NOT_VALID_VALUE = 0
            unique_peptides_offset = self.get_offset(col_index, attr_name,
                                                     item_name,
                                                     UNIQUE_PEPTIDES)

        self.rows = rows
        self.item_data = item_data

    def get_offset(self, from_col_index, from_attr_name, from_item_name,
                   to_attr_name_):
        return [
            to_col_index - from_col_index
            for to_col_index, to_item_name, to_attr_name in self.item_attrs
            if to_attr_name == to_attr_name_ and to_item_name == from_item_name
        ][0]

    def unique_peptides_valid_row_criteria(self, row):
        return utils.unique_peptides_valid_row(
            row[self.unique_peptides_first_column_index:self.
                unique_peptides_first_column_index + len(self.item_names)])
        #column_values = row[]

    def filter_helper(self, item_names, attr_name, criteria):
        item1, item2 = item_names
        v1 = self.item_data[item1][attr_name]
        v2 = self.item_data[item2][attr_name]
        temp = [(i, v1[i], v2[i]) for i in range(len(self.Genes))
                if criteria(v1[i], v2[i])]
        return [(self.UniprotIds[i], self.Genes[i],
                 self.ProteinDescriptions[i], v1, v2) for i, v1, v2 in temp]

    def common_proteins(self, item_names, attr_name):
        return self.filter_helper(item_names, attr_name,
                                  common_proteins_criteria)

    def fold_change(self, item_names, attr_name):
        return self.filter_helper(item_names, attr_name, is_two_fold)

    def save_csv(self, csv_file_name, attr_name):
        first_col_index = self.attr_first_column_index[attr_name]
        item_name_len = len(self.item_names)

        ft.write_csv(
            csv_file_name, self.header[:3] +
            self.header[first_col_index:first_col_index + item_name_len], [
                row[:3] + row[first_col_index:first_col_index + item_name_len]
                for row in self.rows
            ])


script_dir = ft.get_script_dir(__file__)
fname = 'Batch2_data.csv'

m = MsAnalysis(path.join(script_dir, fname))
m.save_csv(path.join(script_dir, 'csv.csv'), LABEL_FREE_QUANT)
item_names = ['M2', 'P2']
common_proteins = m.common_proteins(item_names, LABEL_FREE_QUANT)
fold_change = m.fold_change(item_names, LABEL_FREE_QUANT)
ft.write_csv(path.join(script_dir, 'csv.csv'), m.header[:3] + item_names,
             common_proteins)
