import csvtopandas
import fileutils as ft
from os import path
import plotutils as pl

def col_starts_with(startswith):
    return lambda col: col.str.startswith(startswith)

prots = ft.read_file(ft.msdata_filename('proteomics_experiment_1_common_manual_cellular_proteins.txt'))

fname = 'Proteomics_experiment_2.tsv'
fname = 'Batch2_data.csv'
p = csvtopandas.CsvToPandas(ft.msdata_filename(fname))
df = p.filtered[p.filtered['PG.Genes'].isin(prots)]
pl.dataframe_plot(df, lambda df: df.plot(x=csvtopandas.PG_PROTEINDESCRIPTIONS,
                                   y=p.abundance_col_names,
                                   kind='bar',
                                   rot=0,
                                   legend=False), 'A title')
#p.to_csv(ft.msdata_csv_filename(fname))
p.to_gene_list(ft.msdata_gene_filename)

#fold_groups = (('adult_mean', ('M1', 'M2')), ('piglet_mean', ('P1', 'P2')))
#protein_groups = ('Collagen','Laminin')
protein_description_filters = map(lambda s: col_starts_with(s), protein_groups)
#p.fold_analysis(fold_groups, protein_description_filters)
N = 10
p.group_analysis(('M1', 'M2', 'P1', 'P2'), protein_description_filters, N)
