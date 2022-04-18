from os import path
import pandas as pd
import csv

UNIQUE_PEPTIDES_COL_PREFIX = 'UniquePeptides'
SKIP_COLUMNS = ('PG.ProteinDescriptions',)

script_dir = path.dirname(path.realpath(__file__))
fname = 'Batch2_data.csv'

with open(path.join(script_dir, fname)) as csvfile:
    rdr = csv.reader(csvfile, delimiter=';', quotechar='')
    rows = [row for row in rdr]

idxs = [i for i, c in enumerate(df.columns) if c.find('UNIQUE_PEPTIDES_COL_PREFIX') > 0]
rng = range(len(df.values))
valid_rows = [
    r for r in rng
    if [idx for idx in idxs if df.values[r][idx] > 1]
]

out_columns = set(df.columns) - set(SKIP_COLUMNS)
skip_rows = set(rng) - set(valid_rows)
df.drop(skip_rows)

csv_out_name = path.join(script_dir, '%s.output.csv' % (path.splitext(fname)[0],))
df.to_csv(csv_out_name, columns=out_columns)
