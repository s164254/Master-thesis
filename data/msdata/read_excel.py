from os import path
import pandas as pd

script_dir = path.dirname(path.realpath(__file__))
fname = 'Batch2_data.csv'
data = pd.read_csv(path.join(script_dir, fname), header='infer')

idxs = [i for i, c in enumerate(data.columns) if c.find('UniquePeptides') > 0]
rng = range(len(data.values))
valid_rows = [
    r for r in rng
    if [idx for idx in idxs if data.values[r][idx] > 1]
]
