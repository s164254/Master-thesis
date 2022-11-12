import fileutils as ft
from dotdict import dotdict

attr_actions = (
    ('sample_name_lookup',lambda n: ft.to_dict('%s.samplenames' % (n,),'=')),
    ('common_proteins',lambda n: ft.to_list('%s.cellular.common' % (n,),'|')),
)

def to_experiment_args(experiment_name):
    res = {}
    res['filename'] = ft.msdata_filename(experiment_name,'csv')
    for attr_name, action in attr_actions:
        res[attr_name] = action(experiment_name)
    return dotdict(res)
