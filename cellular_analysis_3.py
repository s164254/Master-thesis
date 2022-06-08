import csvtopandas

exp = csvtopandas.CsvToPandas('proteomics_experiment_2')
all_sample_groups = (
    ('Fa1a_500ng', 'Fa2a_500ng'),
    ('Fa1a_500ng', 'T1_500ng', 'T2_500ng', 'T3_500ng'),
)

def title_func(experiment_name, display_columns): 
    return 'Cellular analysis plot. Experiment:%s, sample names:%s' % (experiment_name,', '.join(display_columns))

for sample_names in all_sample_groups:
    exp.cellular_analysis_3(sample_names, title_func, 'X label', 'Y label')
