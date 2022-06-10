import csvtopandas

N = 20
exp = csvtopandas.CsvToPandas('proteomics_experiment_1')
exp.cellular_analysis(N, 'Abundance of selected cellular proteins for the samples', 'Gene names', 'Abundance', ('M1', 'P1'))
exp.cellular_analysis(N, 'Abundance of selected cellular proteins for the samples', 'Gene names', 'Abundance', ('M1', 'M2'))

