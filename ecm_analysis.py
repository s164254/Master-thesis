import csvtopandas
import warnings

warnings.filterwarnings("ignore")

def col_starts_with(startswith):
    return lambda col: col.str.startswith(startswith)

exp = csvtopandas.CsvToPandas('proteomics_experiment_1')

N = 10
protein_groups = ('Collagen','Laminin','Fibronectin','Elastin','Proteoglycan')
for protein_group in protein_groups:
    protein_description_filter = col_starts_with(protein_group)
    exp.fold_analysis('M1', 'M2')
