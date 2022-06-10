import csvtopandas

exp = csvtopandas.CsvToPandas('proteomics_experiment_2')
all_sample_groups = (
    (('Fa1a_500ng', 'Fa2a_500ng'), 'Abundance of most abundant cellular proteins for the samples'),
    (('T1_500ng', 'T2_500ng', 'T3_500ng'), 'Abundance of most abundant cellular proteins for the samples'),
     (('Fa1a_500ng', 'T3_500ng'), 'Abundance of most abundant cellular proteins for the samples'),
)

def title_func(experiment_name, display_columns): 
    return 'Cellular analysis plot. Experiment:%s, sample names:%s' % (experiment_name,', '.join(display_columns))

def title_closure(title): 
    def f(experiment_name, display_columns):
        return title
    return f

for sample_names, chart_title in all_sample_groups:
    exp.cellular_analysis_3(sample_names, title_closure(chart_title), 'Gene names', 'Abundance')


#Have option to remove gene names from plot.