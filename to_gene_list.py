import csvtopandas

for e in ('proteomics_experiment_1','proteomics_experiment_2'):
    exp = csvtopandas.CsvToPandas(e)
    exp.to_gene_list()
