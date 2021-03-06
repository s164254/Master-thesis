import fileutils as ft
from dotdict import dotdict
from os import path, makedirs, listdir


def load_common_proteins(experiment_name):
    res = {}
    msdata_dir = ft.msdata_filename('.')
    cellular_files = [
        path.join(msdata_dir, f) for f in listdir(msdata_dir)
        if f.startswith('%s.cellular.' % (experiment_name, ))
    ]
    for f in cellular_files:
        content = ft.to_list(f, '|')
        key = f.split('.')[-3]
        res[key] = content
    return res


attr_actions = (
    ('sample_name_lookup',
     lambda n: ft.to_dict(ft.msdata_filename('%s.samplenames.txt' % (n, )), '=')),
    ('common_proteins', load_common_proteins),
)


def to_experiment_args(experiment_name):
    res = dotdict({})
    res['experiment_name'] = experiment_name
    res['filename'] = ft.msdata_filename(experiment_name, 'csv')
    for attr_name, action in attr_actions:
        res[attr_name] = action(experiment_name)
    output_dir = path.join(ft.msdata_filename('.'), 'output')
    input_dir = path.join(ft.msdata_filename('.'), 'input')

    if not path.exists(output_dir):
        makedirs(output_dir)

    res['experiment_output_dir'] = ft.ensure_path_exists(path.join(output_dir, experiment_name))
    res['experiment_input_dir'] = ft.ensure_path_exists(path.join(input_dir, experiment_name))
    res['gene_filename'] = lambda fn: ft.ensure_path_exists(path.join(res.experiment_output_dir,
                                                'gene', fn))
    res['fig_filename'] = lambda fn: ft.ensure_path_exists(path.join(res.experiment_output_dir, fn))
    res['common_filename'] = lambda fn: ft.ensure_path_exists(path.join(res.experiment_output_dir,
                                                'ecm_common', fn))
    res['lookup_filename'] = lambda fn: ft.ensure_path_exists(path.join(res.experiment_input_dir, path.basename(fn)))

    return res
