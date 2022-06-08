import csvtopandas


def col_starts_with(startswith):
    return lambda col: col.str.startswith(startswith, na=False)

def col_contains(match_regex):
    return lambda col: col.str.contains(match_regex, case=False, na=False)

experiments = (
    ('proteomics_experiment_1', (('M1', 'P1'), ('M1', 'M2'))),
    ('proteomics_experiment_2', (('Fa1a', 'Fa2a'), ('Fa1a', 'T1', 'T2',
                                                    'T3'))),
)

# I,III,IV,V,VI
collagen_regex = 'Collagen.+ (I|V){1,3}[^\S|$]'
fibronectin_regex = '^Fibronectin$'
protein_description_filter_dict = {
    'Collagen': col_contains(collagen_regex),
    'Fibronectin': col_contains(fibronectin_regex),
    'Laminin': col_starts_with('Laminin'),
    'Elastin': col_starts_with('Elastin'),
    'Proteoglycan': col_starts_with('Proteoglycan')
}

def run_analysis(protein_groups, analysis_func):
    protein_description_filters = [protein_description_filter_dict[protein_group] for protein_group, use_ratio_filter in protein_groups]

    for experiment_name, all_samples in experiments:
        exp = csvtopandas.CsvToPandas(experiment_name)
        for samples in all_samples:
            for protein_group, filter_func in zip(protein_groups, protein_description_filters):
                protein_group, use_ratio_filter = protein_group
                analysis_func(exp, samples, protein_group, use_ratio_filter, filter_func)
