import pandas as pd
from os import path
import fileutils as ft
import plotutils



files = (
    ('PB220307cd_samples_batch_3','PDE3'),
    ('PB220307ad_samples_B1_A1_A2_modified_for_report','PDE1'),
    ('PB220307bd_samples_P1_M1_M2_equal_replicates','PDE2'),
    ('PB220419dd_batch_4','PDE4'),
    ('PB220419dd_batch_5','PDE5'),
    ('PB220522fd_batch_6_real','PDE6')
)

SAMPLE_ID = 'Sample ID'
NUCLEIC_ACID_CONC = 'Nucleic Acid Conc.'
DNA_CONC_NG_MG = 'DNA Conc ng/mg'
WEIGHT = 'Weight (mg)'


def to_csv():
    for fname, title in files[-1:]:
        data = pd.read_excel(resulting_filename(fname))
        m = data.groupby([SAMPLE_ID])[[NUCLEIC_ACID_CONC, WEIGHT]].mean()
        m[DNA_CONC_NG_MG] = 400 * m[NUCLEIC_ACID_CONC] / m[WEIGHT]
        m.to_csv(resulting_filename(fname, 'csv'))


def plot_all():
    for fname, title in files:
        df = pd.read_csv(resulting_filename(fname, 'csv'))
        plotutils.dataframe_plot(
            df,
            lambda df: df.plot(
                x=SAMPLE_ID, y=DNA_CONC_NG_MG, kind='bar', rot=0, legend=False
            ),
            '%s %s' % (DNA_CONC_NG_MG, title),
            lambda ax: ax.bar_label(ax.containers[0]),
            block=True)


#to_csv()
plot_all()