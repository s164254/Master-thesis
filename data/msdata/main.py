from csvtopandas import CsvToPandas
import fileutils as ft
from os import path

def col_starts_with(startswith):
    return lambda col: col.str.startswith(startswith)

fname = 'Batch2_data.csv'
p = CsvToPandas(ft.relative_to_script_dir(__file__, fname))

fold_groups = (('adult_mean', ('M1', 'M2')), ('piglet_mean', ('P1', 'P2')))
protein_groups = ('Collagen','Laminin')
protein_description_filters = map(lambda s: col_starts_with(s), protein_groups)
#p.fold_analysis(fold_groups, protein_description_filters)
N = 10
p.group_analysis(('M1', 'M2', 'P1', 'P2'), protein_description_filters, N)
