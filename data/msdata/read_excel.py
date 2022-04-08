from os import path
import pandas as pd

script_dir = path.dirname(path.realpath(__file__))
fname = 'Batch2_data.csv'
data = pd.read_csv(path.join(script_dir, fname),header='infer')

idxs = [i for i,c in enumerate(data.columns) if c.find('UniquePeptides')>0]
header = data.values[0]
l = len(data)

