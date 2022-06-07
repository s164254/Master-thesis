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
    if not path.exists(fname):
        return None
    with open(fname, 'r') as f:
        return [l.strip() for l in f.readlines() if l.strip()]


def relative_to_script_dir(_file_, fname, ext, data_path=''):
    fdir = path.join(get_script_dir(_file_), data_path)
    if not path.exists(fdir):
        makedirs(fdir, exist_ok=True)
    return path.join(fdir,
                     fname.find('.') >= 0 and fname or '%s.%s' % (fname, ext))


def dnadata_filename(fname, ext='xlsx'):
    return relative_to_script_dir(__file__, fname, ext, DNA_DATA_PATH)


def msdata_filename(fname, ext='xlsx'):
    return relative_to_script_dir(__file__, fname, ext, MS_DATA_PATH)


def msdata_csv_filename(fname):
    return relative_to_script_dir(__file__, fname, 'csv', MS_DATA_CSV_PATH)


def msdata_gene_filename(fname):
    return relative_to_script_dir(__file__, fname, 'csv', MS_DATA_GENE_PATH)


def to_dict(fname, sep):
    content = read_file(fname)
    return content and dict([x.split(sep) for x in content]) or None


def ensure_path_exists(fname):
    fname = fname.replace(' ', '_')
    dirname = path.dirname(fname)
    if dirname and not path.exists(dirname):
        makedirs(dirname, exist_ok=True)
    return fname


def to_list(fname, sep=None):
    content = read_file(fname)
    if not content:
        return None
    return content[0].find(sep) >= 0 and content[0].split(sep) or content


def to_file(fname, content):
    fn = ensure_path_exists(fname)
    print(fn)
    with open(fn, 'w') as f:
        f.write(content)
