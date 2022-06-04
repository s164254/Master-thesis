import csvtopandas
import warnings

warnings.filterwarnings("ignore")


def col_starts_with(startswith):
    return lambda col: col.str.startswith(startswith, na=False)


def col_contains(match_regex):
    return lambda col: col.str.contains(match_regex, case=False, na=False)


exp = csvtopandas.CsvToPandas('proteomics_experiment_1')

N = 10
# I,III,IV,V,VI
collagen_regex = 'Collagen.+ (I|V){1,3}[^\S|$]'
protein_groups = (
    ('Collagen', True),
    ('Laminin', False),
    ('Fibronectin', False),
    ('Elastin', False),
    ('Proteoglycan', False),
)

protein_description_filters = [col_contains(collagen_regex)] + list(
    map(col_starts_with, protein_groups[1:]))

all_samples = (
    ('M1', 'P1'),
    ('M1', 'M2'),
    #('M1', 'M2'),
)

for samples in all_samples:
    for g, f in zip(protein_groups, protein_description_filters):
        g, use_filter = g
        exp.fold_analysis(samples, g, use_filter, f)
