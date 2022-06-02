import csvtopandas
import warnings
warnings.filterwarnings("ignore")

for e in ('proteomics_experiment1','proteomics_experiment2'):
    exp = csvtopandas.CsvToPandas(e)
    exp.to_gene_list()
