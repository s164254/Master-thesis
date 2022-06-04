import csvtopandas
import warnings

warnings.filterwarnings("ignore")

def col_starts_with(startswith):
    return lambda col: col.str.startswith(startswith)

def col_contains(match_regex):
    return lambda col: col.str.contains(match_regex, case=False)

exp = csvtopandas.CsvToPandas('proteomics_experiment_1')

N = 10
# I,III,IV,V,VI
collagen_regex = 'Collagen.+ (I|V)[I|V]*[^\S|$]'
#collagen_regex = 'Collagen.+ [(I)|(III)|(IV)|(V)|(VI)][^\S|$]'
protein_groups = ('Collagen','Laminin','Fibronectin','Elastin','Proteoglycan')

protein_description_filters = [col_contains(collagen_regex)] + list(map( col_starts_with, protein_groups[1:]))
sample_names = ('M1', 'P1')
exp.fold_analysis(sample_names, protein_description_filters)
