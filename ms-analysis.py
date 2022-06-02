import csvtopandas

exp = csvtopandas.CsvToPandas('proteomics_experiment1')
#exp.fold_analysis(None,None)
#exp.to_gene_list()
N = 10
exp.gene_analysis(N, ('M1', 'P1'))
