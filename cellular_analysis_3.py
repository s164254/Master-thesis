import csvtopandas

exp = csvtopandas.CsvToPandas('proteomics_experiment_2')
all_sample_groups = (
    ('Fa1a_500ng', 'Fa2a_500ng'),
    ('Fa1a_500ng', 'T1_500ng', 'T2_500ng', 'T3_500ng'),
)

for sample_names in all_sample_groups:
    exp.cellular_analysis_3(sample_names=sample_names, N=20)
