import csvtopandas
import warnings

warnings.filterwarnings("ignore")


def col_starts_with(startswith):
    return lambda col: col.str.startswith(startswith, na=False)


def col_contains(match_regex):
    return lambda col: col.str.contains(match_regex, case=False, na=False)

experiments = (
    ('proteomics_experiment_1', (('M1', 'P1'), ('M1', 'M2'))),
    ('proteomics_experiment_2', (('Fa1a', 'Fa2a'), ('Fa1a', 'T1', 'T2',
                                                    'T3'))),
)

N = 10

# I,III,IV,V,VI
collagen_regex = 'Collagen.+ (I|V){1,3}[^\S|$]'

def run_analysis(protein_groups, analysis_func):
    protein_description_filters = [col_contains(collagen_regex)] + list(
        map(col_starts_with, protein_groups[1:]))

    for experiment_name, all_samples in experiments:
        exp = csvtopandas.CsvToPandas(experiment_name)
        for samples in all_samples:
            for protein_group, filter_func in zip(protein_groups, protein_description_filters):
                protein_group, use_filter = protein_group
                analysis_func(exp, samples, protein_group, use_filter, filter_func)
