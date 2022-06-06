import csvtopandas
import warnings

warnings.filterwarnings("ignore")

exp = csvtopandas.CsvToPandas('proteomics_experiment_2')
exp.cellular_analysis_2()
