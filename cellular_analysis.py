import csvtopandas

N = 20
exp = csvtopandas.CsvToPandas('proteomics_experiment_1')
exp.cellular_analysis(N, 'a title', 'bla bla', 'bli bli', ('M1', 'P1'))
exp.cellular_analysis(N, 'a title', 'bla bla', 'bli bli', ('M1', 'M2'))
