import csvtopandas
import warnings
warnings.filterwarnings("ignore")

N=10
exp = csvtopandas.CsvToPandas('proteomics_experiment1')
exp.gene_analysis(N, ('M1', 'P1'))
