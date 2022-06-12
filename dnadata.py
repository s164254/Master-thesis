import pandas as pd
from os import path
import fileutils as ft
import plotutils



files = [
    ('PB220307cd_samples_batch_3','PDE3'),
    ('PB220307ad_samples_B1_A1_A2_modified_for_report','PDE1'),
    ('PB220307bd_samples_P1_M1_M2_equal_replicates','PDE2'),
    ('PB220419dd_batch_4','PDE4'),
    ('PB220419dd_batch_5','PDE5'),
    ('PB220522fd_batch_6_real','PDE6'),
#    ('PB220522fd_batch_5_and_6','PDE5 and PDE6'),
    ('PB220522fd_batch_final_corrected','PDE7')
]

plot_files = files + [('PB220522fd_batch_5_and_6','PDE5 and PDE6')]

SAMPLE_ID = 'Sample ID'
NUCLEIC_ACID_CONC = 'Nucleic Acid Conc.'
DNA_CONC_NG_MG = 'DNA concentration'
WEIGHT = 'Weight (mg)'
MULT_FACTOR = 'MultFactor'

mult_factor_dict = dict( (
    ('PB220522fd_batch_final_corrected', 'Fa1a'),
    ('PB220522fd_batch_final_corrected', 'Fa2a'),
))

def get_mult_factor(key):
    return mult_factor_dict.get(key, 400)

def to_csv_1():
    for fname, title in files:
        full_fname = ft.dnadata_filename(fname)
        if not path.exists(full_fname):
            continue

        data = pd.read_excel(ft.dnadata_filename(fname))
        m = data.groupby([SAMPLE_ID])[[NUCLEIC_ACID_CONC, WEIGHT]].mean()
        #key = (fname,)
        m[MULT_FACTOR] = 400
        m.to_csv(ft.dnadata_filename('%s_1' % (fname,) , 'csv'))

def to_csv_2():
    for fname, title in files:
        df = pd.read_csv(ft.dnadata_filename('%s_1' % (fname,), 'csv'))
        df[DNA_CONC_NG_MG] = df[MULT_FACTOR] * df[NUCLEIC_ACID_CONC] / df[WEIGHT]
        df.to_csv(ft.dnadata_filename('%s_2' % (fname,), 'csv'))


def plot_all():
    for fname, title in plot_files:
        df = pd.read_csv(ft.dnadata_filename('%s_2' % (fname,), 'csv'))
        plotutils.dataframe_plot(
            df,
            lambda df: df.plot(
                x=SAMPLE_ID, y=DNA_CONC_NG_MG, kind='bar', rot=0, legend=False
            ),
            '%s %s' % (DNA_CONC_NG_MG, title),
            lambda ax: ax.bar_label(ax.containers[0]),
            block=True)


# create csv with calculated mean and mult factor
#to_csv_1()

# create csv with DNA conc. after correcting mult. factor
#to_csv_2()

plot_all()