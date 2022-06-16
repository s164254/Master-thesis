import pandas as pd
from os import path
import fileutils as ft
import plotutils
import re



files = [
   ('PB220522fd_batch_5_and_6','in subset of samples from experiments FilterDe and SDSDe'),
    ('PB220307cd_samples_batch_3','in samples from experiment BlendDe'),
    ('PB220307ad_samples_B1_A1_A2_modified_for_report','in samples from experiment IniDe'),
    ('PB220307bd_samples_P1_M1_M2_equal_replicates','in samples from experiment PatiDe'),
    ('PB220419dd_batch_4','in samples from experiment FineDe'),
    ('PB220419dd_batch_5','in samples from experiment FilterDe'),
    ('PB220522fd_batch_6_real','in samples from experiment SDSDe'),
    ('PB220522fd_batch_final_corrected','in samples from experiment HarshDe')
]

plot_files = files #+ [('PB220522fd_batch_5_and_6','PDE5 and PDE6')]

SAMPLE_ID = 'Sample ID'
DISPLAY_SAMPLE_ID = 'Samples'
NUCLEIC_ACID_CONC = 'Nucleic Acid Conc.'
DNA_CONC_NG_MG = 'DNA concentration'
DISPLAY_DNA_CONC_NG_MG = 'DNA concentration'
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
        if path.exists(full_fname):
            continue

        data = pd.read_excel(ft.dnadata_filename(fname))
        m = data.groupby([SAMPLE_ID])[[NUCLEIC_ACID_CONC, WEIGHT]].mean()
        #key = (fname,)
        m[MULT_FACTOR] = 400
        m.to_csv(ft.dnadata_filename('%s_1' % (fname,) , 'csv'))

def to_csv_2():
    for fname, title in files:
        out_fname = ft.dnadata_filename('%s_2' % (fname,), 'csv')
        if path.exists(out_fname):
            continue

        df = pd.read_csv(ft.dnadata_filename('%s_1' % (fname,), 'csv'))
        df[DNA_CONC_NG_MG] = df[MULT_FACTOR] * df[NUCLEIC_ACID_CONC] / df[WEIGHT]
        df[DNA_CONC_NG_MG] = df[DNA_CONC_NG_MG].apply(int)
        df.to_csv(out_fname)


def plot_all():
    for fname, title in plot_files:
        df = pd.read_csv(ft.dnadata_filename('%s_2' % (fname,), 'csv'))
        df.rename( columns={SAMPLE_ID:DISPLAY_SAMPLE_ID, DNA_CONC_NG_MG:DISPLAY_DNA_CONC_NG_MG}, inplace=True)
        df[DISPLAY_SAMPLE_ID] = df.apply(
            lambda x: re.sub('\s+', '\n', x[DISPLAY_SAMPLE_ID]), axis=1)

        plotutils.dataframe_plot(
            df,
            lambda df: df.plot(
                x=DISPLAY_SAMPLE_ID, y=DISPLAY_DNA_CONC_NG_MG, kind='bar', rot=0, legend=False
            ),
            '%s %s' % (DISPLAY_DNA_CONC_NG_MG, title),
            ylabel='DNA concentration [ng/mg tissue weight]',
            axis_setup_func=lambda ax: ax.bar_label(ax.containers[0]),
            block=True,
            fig_filename=ft.dnadata_chart_filename('%s.png' % (fname,)),
            minor_ticks=False)


# create csv with calculated mean and mult factor
#to_csv_1()

# create csv with DNA conc. after correcting mult. factor
to_csv_2()

plot_all()