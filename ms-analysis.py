import csvtopandas

N=10
exp = csvtopandas.CsvToPandas('proteomics_experiment1')
exp.gene_analysis(N, ('M1', 'P1'))

exp = csvtopandas.CsvToPandas('proteomics_experiment2')
exp.gene_analysis(N, ('ABC', 'DEF'))
