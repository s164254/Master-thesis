import csvtopandas
import fileutils as ft
from os import path
import plotutils as pl

def col_starts_with(startswith):
    return lambda col: col.str.startswith(startswith)

prots = ft.read_file(ft.msdata_filename('proteomics_experiment_1_common_manual_cellular_proteins.txt'))

fname = 'Proteomics_experiment_2.tsv'
fname = 'Batch2_data.csv'

# csv with diuplicates
#df = p.filtered[p.filtered.duplicated(['PG.Genes'])].sort_values(by=['PG.Genes'])
#df.to_csv(ft.msdata_gene_filename('dups'))

p = csvtopandas.CsvToPandas(ft.msdata_filename(fname))
df = p.filtered[p.filtered['PG.Genes'].isin(prots)]
m1_p1 = ['M1','P1']
sns = [sn for sn in p.abundance_col_names if any([x for x in m1_p1 if sn.find(x)>0])]
pl.dataframe_plot(df, lambda df: df.plot(x=csvtopandas.PG_PROTEINDESCRIPTIONS,
                                   y=sns,
                                   kind='bar',
                                   rot=0,
                                   legend=False), 'A title', lambda ax: ax.legend(m1_p1))
#p.to_csv(ft.msdata_csv_filename(fname))
p.to_gene_list()

#fold_groups = (('adult_mean', ('M1', 'M2')), ('piglet_mean', ('P1', 'P2')))
#protein_groups = ('Collagen','Laminin')
protein_description_filters = map(lambda s: col_starts_with(s), protein_groups)
#p.fold_analysis(fold_groups, protein_description_filters)
N = 10
p.group_analysis(('M1', 'M2', 'P1', 'P2'), protein_description_filters, N)
