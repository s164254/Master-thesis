from os import path
import csv

DNA_DATA_PATH = path.join('data', 'DNA_data')
MS_DATA_PATH = path.join('data', 'msdata')


def get_script_dir(_file_):
    return path.dirname(path.realpath(_file_))


def relative_to_script_dir(_file_, fname, ext, data_path=''):
    return path.join(get_script_dir(_file_), data_path,
                     fname.find('.') > 0 and fname or '%s.%s' % (fname, ext))


def dnadata_filename(fname, ext='xlsx'):
    return relative_to_script_dir(__file__, fname, ext, DNA_DATA_PATH)


def msdata_filename(fname, ext='xlsx'):
    return relative_to_script_dir(__file__, fname, ext, MS_DATA_PATH)
