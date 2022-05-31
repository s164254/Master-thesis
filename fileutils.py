from ntpath import join
from os import path, makedirs
import csv

DNA_DATA_PATH = path.join('data', 'DNA_data')
MS_DATA_PATH = path.join('data', 'msdata')
MS_DATA_CSV_PATH = path.join(MS_DATA_PATH, 'csv')
MS_DATA_GENE_PATH = path.join(MS_DATA_PATH, 'gene')


def get_script_dir(_file_):
    return path.dirname(path.realpath(_file_))

def read_file(fname):
    with open(fname,'r') as f:
        return [l.strip() for l in f.readlines() if l.strip()]

def relative_to_script_dir(_file_, fname, ext, data_path=''):
    fdir = path.join(get_script_dir(_file_), data_path)
    if not path.exists(fdir):
        makedirs(fdir)
    return path.join(fdir,
                     fname.find('.') > 0 and fname or '%s.%s' % (fname, ext))


def dnadata_filename(fname, ext='xlsx'):
    return relative_to_script_dir(__file__, fname, ext, DNA_DATA_PATH)


def msdata_filename(fname, ext='xlsx'):
    return relative_to_script_dir(__file__, fname, ext, MS_DATA_PATH)


def msdata_csv_filename(fname):
    return relative_to_script_dir(__file__, fname, 'csv', MS_DATA_CSV_PATH)

def msdata_gene_filename(fname):
    return relative_to_script_dir(__file__, fname, 'csv', MS_DATA_GENE_PATH)
