import csvtopandas
import warnings

warnings.filterwarnings("ignore")

exp = csvtopandas.CsvToPandas('proteomics_experiment_2')
exp.cellular_analysis_2(N=20, sample_names=None)
