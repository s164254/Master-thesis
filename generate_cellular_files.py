import csvtopandas
import warnings

warnings.filterwarnings("ignore")

exp = csvtopandas.CsvToPandas('proteomics_experiment_2')
exp.generate_cellular_file()