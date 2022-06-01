import fileutils as ft
from dotdict import dotdict
from os import path, makedirs

attr_actions = (
    ('sample_name_lookup',
     lambda n: ft.to_dict(ft.msdata_filename('%s.samplesnames' % (n, )), '=')),
    ('common_proteins',
     lambda n: ft.to_list(ft.msdata_filename('%s.cellular.common' %
                                             (n, )), '|')),
)


def to_experiment_args(experiment_name):
    res = dotdict({})
    res['filename'] = ft.msdata_filename(experiment_name, 'csv')
    for attr_name, action in attr_actions:
        res[attr_name] = action(experiment_name)
    output_dir = path.join(ft.msdata_filename('.'), 'output')

    if not path.exists(output_dir):
        makedirs(output_dir)

    res['experiment_output_dir'] = path.join(output_dir, experiment_name)
    if not path.exists(res.experiment_output_dir):
        makedirs(res.experiment_output_dir)
    
    return res
